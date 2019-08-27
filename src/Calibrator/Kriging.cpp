#include "Kriging.h"
#include "../Util.h"
#include "../Parameters.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <math.h>

CalibratorKriging::CalibratorKriging(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mEfoldDist(30000),
      mMaxElevDiff(Util::MV),
      mAuxVariable(""),
      mLowerThreshold(Util::MV),
      mUpperThreshold(Util::MV),
      mCrossValidate(false),
      mRadius(30000),
      mOperator(Util::OperatorAdd),
      mUseApproxDistance(true),
      mKrigingType(TypeCressman),
      mWindow(0) {
   iOptions.getValue("efoldDist", mEfoldDist);
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("maxElevDiff", mMaxElevDiff);
   std::string parfilename;
   if(mEfoldDist< 0) {
      Util::error("CalibratorKriging: 'efoldDist' must be >= 0");
   }
   if(Util::isValid(mMaxElevDiff) && mMaxElevDiff < 0) {
      Util::error("CalibratorKriging: 'maxElevDiff' must be >= 0");
   }
   if(mRadius < 0) {
      Util::error("CalibratorKriging: 'radius' must be >= 0");
   }

   std::string type;
   if(iOptions.getValue("type", type)) {
      if(type == "cressman") {
         mKrigingType = TypeCressman;
      }
      else if(type == "barnes") {
         mKrigingType = TypeBarnes;
      }
      else {
         Util::error("CalibratorKriging: 'type' not recognized");
      }
   }

   iOptions.getValue("approxDist", mUseApproxDistance);

   std::string op;
   if(iOptions.getValue("operator", op)) {
      if(op == "add") {
         mOperator = Util::OperatorAdd;
      }
      else if(op == "subtract") {
         mOperator = Util::OperatorSubtract;
      }
      else if(op == "multiply") {
         mOperator = Util::OperatorMultiply;
      }
      else if(op == "divide") {
         mOperator = Util::OperatorDivide;
      }
      else {
         Util::error("CalibratorKriging: 'operator' not recognized");
      }
   }

   std::string auxVariable;
   if(iOptions.getValue("auxVariable", mAuxVariable)) {
      std::vector<float> range;
      if(iOptions.getValues("range", range)) {
         if(range.size() != 2) {
            Util::error("CalibratorKriging: 'range' must be of the form lower,upper");
         }
         mLowerThreshold = range[0];
         mUpperThreshold = range[1];
      }
      else {
         Util::error("CalibratorKriging: 'range' required if using 'auxVariable'.");
      }

      iOptions.getValue("window", mWindow);
      if(mLowerThreshold > mUpperThreshold) {
         Util::error("CalibratorKriging: the lower value must be less than upper value in 'range'");
      }
   }
   iOptions.getValue("crossValidate", mCrossValidate);
   iOptions.check();
}

bool CalibratorKriging::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   // Check if this method can be applied
   bool hasValidGridpoint = false;
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         if(Util::isValid(lats[i][j]) && Util::isValid(lons[i][j]) && Util::isValid(elevs[i][j])) {
            hasValidGridpoint = true;
         }
      }
   }
   if(!hasValidGridpoint) {
      Util::warning("There are no gridpoints with valid lat/lon/elev values. Skipping kriging...");
      return false;
   }

   // Precompute weights from auxillary variable
   std::vector<std::vector<std::vector<std::vector<float> > > > auxWeights;
   if(mAuxVariable != "") {
      // Initialize
      auxWeights.resize(nLat);
      for(int i = 0; i < nLat; i++) {
         auxWeights[i].resize(nLon);
         for(int j = 0; j < nLon; j++) {
            auxWeights[i][j].resize(nEns);
            for(int e = 0; e < nEns; e++) {
               auxWeights[i][j][e].resize(nTime, 0);
            }
         }
      }
      // Load auxillarcy variable
      std::vector<FieldPtr> auxFields;
      auxFields.resize(nTime);
      for(int t = 0; t < nTime; t++) {
         auxFields[t] = iFile.getField(mAuxVariable, t);
      }

      // Compute auxillary weights
      for(int t = 0; t < nTime; t++) {
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               for(int e = 0; e < nEns; e++) {
                  float total = 0;
                  int start = std::max(t-mWindow,0);
                  int end   = std::min(nTime-1,t+mWindow);
                  int numValid = 0;
                  for(int tt = start; tt <= end; tt++) {
                     float aux = (*auxFields[tt])(i,j,e);
                     if(Util::isValid(aux)) {
                        if(aux >= mLowerThreshold && aux <= mUpperThreshold) {
                           total++;
                        }
                        numValid++;
                     }
                  }
                  int windowSize = end - start + 1;
                  if(numValid == 0)
                     auxWeights[i][j][e][t] = 1;
                  else
                     auxWeights[i][j][e][t] += total / numValid;
               }
            }
         }
      }
   }

   if(!iParameterFile->isLocationDependent()) {
      std::stringstream ss;
      ss << "Kriging requires a parameter file with spatial information";
      Util::error(ss.str());
   }
   std::vector<Location> obsLocations = iParameterFile->getLocations();

   // General proceedure for a given gridpoint:
   // S              = matrix * weights
   // weights        = (matrix)^-1 * S
   // gridpoint_bias = weights' * bias (scalar)
   // where
   // matrix:  The obs-to-obs covariance matrix (NxN)
   // S:       The obs-to-current_grid_point covariance (Nx1)
   // bias:    The bias at each obs location (Nx1)
   //
   // Note that when computing the weights, we can take a shortcut since most values in
   // S are zero. However, the weights will still have the length of all stations (not the
   // number of nearby stations), since when computing the bias we still need to 
   // weight in far away biases (because they can covary with the nearby stations)

   // Compute obs-obs covariance-matrix once
   vec2 matrix;
   int N = obsLocations.size();
   std::cout << "      Point locations: " << N << std::endl;
   matrix.resize(N);
   for(int ii = 0; ii < N; ii++) {
      matrix[ii].resize(N,0);
   }
   for(int ii = 0; ii < N; ii++) {
      Location iloc = obsLocations[ii];
      // The diagonal is 1, since the distance from a point to itself
      // is 0, therefore its weight is 1.
      matrix[ii][ii] = 1;
      // The matrix is symmetric, so only compute one of the halves
      for(int jj = ii+1; jj < N; jj++) {
         Location jloc = obsLocations[jj];
         // Improve conditioning of matrix when you have two or more stations
         // that are very close
         float factor = 0.414 / 0.5;
         float R = calcCovar(iloc, jloc)*factor;
         // Store the number in both halves
         matrix[ii][jj] = R;
         matrix[jj][ii] = R;
      }
   }

   // Compute (matrix)^-1
   std::cout << "      Precomputing inverse of obs-to-obs covariance matrix: ";
   std::cout.flush();
   double s1 = Util::clock();
   vec2 inverse = Util::inverse(matrix);
   double e1 = Util::clock();
   std::cout << e1 - s1 << " seconds" << std::endl;

   // Compute grid-point to obs-point covariances
   std::cout << "      Precomputing gridpoint-to-obs covariances: ";
   std::cout.flush();
   double s2 = Util::clock();
   // Store the covariances of each gridpoint to every obs-point. To save memory, 
   // only store values that are above 0. Store the index of the obs-point.
   // This means that Sindex does not have the same size for every gridpoint.
   std::vector<std::vector<std::vector<float> > > S; // lat, lon, obspoint
   std::vector<std::vector<std::vector<int> > > Sindex;
   S.resize(nLat);
   Sindex.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      S[i].resize(nLon);
      Sindex[i].resize(nLon);
   }

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         float lat = lats[i][j];
         float lon = lons[i][j];
         float elev = elevs[i][j];
         const Location gridPoint(lat, lon, elev);
         for(int ii = 0; ii < N; ii++) {
            Location obsPoint = obsLocations[ii];
            float covar = calcCovar(obsPoint, gridPoint);
            if(covar > 0) {
               S[i][j].push_back(covar);
               Sindex[i][j].push_back(ii);
            }
         }
      }
   }
   double e2 = Util::clock();
   std::cout << e2 - s2 << " seconds" << std::endl;

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      FieldPtr accum = iFile.getEmptyField(0);
      FieldPtr weights = iFile.getEmptyField(0);

      // Arrange all the biases for all stations into one vector
      std::vector<float> bias(N,0);
      for(int k = 0; k < obsLocations.size(); k++) {
         Location loc = obsLocations[k];
         Parameters parameters = iParameterFile->getParameters(t, loc, false);
         if(parameters.size() > 0) {
            float currBias = parameters[0];
            if(Util::isValid(currBias)) {
               // For * and /, operate on the flucuations areound a mean of 1
               if(mOperator == Util::OperatorMultiply) {
                  currBias = currBias - 1;
               }
               else if(mOperator == Util::OperatorDivide) {
                  currBias = currBias - 1;
               }
            }
            bias[k] = currBias;
         }
      }

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {

            std::vector<float> currS = S[i][j];
            std::vector<int> currI = Sindex[i][j];
            int currN = currS.size();

            // No correction if there are no nearby stations
            if(currN == 0)
               continue;

            // Don't use the nearest station when cross validating
            float maxCovar = Util::MV;
            int ImaxCovar = Util::MV;
            if(mCrossValidate) {
               for(int ii = 0; ii < currS.size(); ii++) {
                  if(!Util::isValid(ImaxCovar) || currS[ii] > maxCovar) {
                     ImaxCovar = ii;
                     maxCovar = currS[ii];
                  }
               }
               currS[ImaxCovar] = 0;
               vec2 cvMatrix = matrix;
               for(int ii = 0; ii < currN; ii++) {
                  cvMatrix[ImaxCovar][ii] = 0;
                  cvMatrix[ii][ImaxCovar] = 0;
               }
               cvMatrix[ImaxCovar][ImaxCovar] = 1;
               inverse = Util::inverse(cvMatrix);
            }

            // Compute weights (matrix-vector product)
            std::vector<float> weights;
            weights.resize(N, 0);
            for(int ii = 0; ii < N; ii++) {
               // Only loop over non-zero values in the vector
               for(int jj = 0; jj < currN; jj++) {
                  int JJ = currI[jj];
                  weights[ii] += inverse[ii][JJ] * currS[jj];
               }
            }
            // Set the weight of the nearest location to 0 when cross-validating
            if(mCrossValidate) {
               weights[ImaxCovar] = 0;
            }

            // Compute final bias (dot product of bias and weights)
            float finalBias = 0;
            for(int ii = 0; ii < N; ii++) {
               float currBias = bias[ii];
               if(!Util::isValid(currBias)) {
                  finalBias = Util::MV;
                  break;
               }
               finalBias += bias[ii]*weights[ii];
            }

            if(Util::isValid(finalBias)) {
               // Reconstruct the factor/divisor by adding the flucuations
               // onto the mean of 1
               if(mOperator == Util::OperatorMultiply)
                  finalBias = finalBias + 1;
               else if(mOperator == Util::OperatorDivide)
                  finalBias = finalBias - 1;

               // Apply bias to each ensemble member
               for(int e = 0; e < nEns; e++) {
                  float rawValue = (*field)(i,j,e);

                  // Adjust bias based on auxillary weight
                  if(mAuxVariable != "") {
                     float weight = auxWeights[i][j][e][t];
                     if(mOperator == Util::OperatorAdd || mOperator == Util::OperatorSubtract) {
                        finalBias = finalBias * weight;
                     }
                     else {
                        finalBias = pow(finalBias, weight);
                     }
                  }

                  if(mOperator == Util::OperatorAdd) {
                     (*field)(i,j,e) += finalBias;
                  }
                  else if(mOperator == Util::OperatorSubtract) {
                     (*field)(i,j,e) -= finalBias;
                  }
                  else if(mOperator == Util::OperatorMultiply) {
                     // TODO: How do we ensure that the matrix is positive definite in this
                     // case?
                     (*field)(i,j,e) *= finalBias;
                  }
                  else if(mOperator == Util::OperatorDivide) {
                     // TODO: How do we ensure that the matrix is positive definite in this
                     // case?
                     (*field)(i,j,e) /= finalBias;
                  }
                  else {
                     Util::error("Unrecognized operator in CalibratorKriging");
                  }
               }
            }
         }
      }
   }
   return true;
}

float CalibratorKriging::calcCovar(const Location& loc1, const Location& loc2) const {
   float weight = Util::MV;
   bool isValidLoc1 = Util::isValid(loc1.lat()) && Util::isValid(loc1.lon()) && Util::isValid(loc1.elev());
   bool isValidLoc2 = Util::isValid(loc2.lat()) && Util::isValid(loc2.lon()) && Util::isValid(loc2.elev());
   if(isValidLoc1 && isValidLoc2) {
      float horizDist = Util::getDistance(loc1.lat(), loc1.lon(), loc2.lat(), loc2.lon(), mUseApproxDistance);
      float vertDist  = fabs(loc1.elev() - loc2.elev());
      if(horizDist >= mRadius || (Util::isValid(mMaxElevDiff) && vertDist >= mMaxElevDiff)) {
         // Outside radius of influence
         return 0;
      }
      if(mKrigingType == TypeCressman) {
         if(horizDist > mEfoldDist || (Util::isValid(mMaxElevDiff) && vertDist > mMaxElevDiff))
            return 0;
         float horizWeight = (mEfoldDist*mEfoldDist - horizDist*horizDist) /
                             (mEfoldDist*mEfoldDist + horizDist*horizDist);
         float vertWeight = 1;
         if(Util::isValid(mMaxElevDiff)) {
            vertWeight  = (mMaxElevDiff*mMaxElevDiff - vertDist*vertDist) /
                                (mMaxElevDiff*mMaxElevDiff + vertDist*vertDist);
         }
         weight = horizWeight * vertWeight;
      }
      else if(mKrigingType == TypeBarnes) {
         float horizWeight = exp(-horizDist*horizDist/(2*mEfoldDist*mEfoldDist));
         float vertWeight  = 1;
         if(Util::isValid(mMaxElevDiff)) {
            vertWeight  = exp(-vertDist*vertDist/(2*mMaxElevDiff*mMaxElevDiff));
         }
         weight = horizWeight * vertWeight;
         return weight;
      }
   }
   return weight;
}

Parameters CalibratorKriging::train(const std::vector<ObsEns>& iData) const {
   double timeStart = Util::clock();
   float totalObs = 0;
   float totalFcst = 0;
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      float obs = iData[i].first;
      std::vector<float> ens = iData[i].second;
      float mean = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(mean)) {
         totalObs += obs;
         totalFcst += mean;
         counter++;
      }
   }
   float bias = Util::MV;
   if(counter <= 0) {
      std::stringstream ss;
      ss << "CalibratorKriging: No valid data, no correction will be made.";
      Util::warning(ss.str());
      bias = 0;
   }
   else {
      if(mOperator == Util::OperatorAdd) {
         bias = (totalObs - totalFcst)/counter;
      }
      else if(mOperator == Util::OperatorSubtract) {
         bias = (totalFcst - totalObs)/counter;
      }
      else if(mOperator == Util::OperatorMultiply) {
         bias = totalObs / totalFcst;
      }
      else if(mOperator == Util::OperatorDivide) {
         bias = totalFcst / totalObs;
      }
      else {
         Util::error("CalibratorKriging: 'operator' not recognized");
      }
   }

   std::vector<float> values(1, bias);
   Parameters par(values);

   double timeEnd = Util::clock();
   std::cout << "Time: " << timeEnd - timeStart << std::endl;
   return par;
}

std::string CalibratorKriging::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c kriging","Spreads bias in space by using kriging. A parameter file is required, which must have one column with the bias.")<< std::endl;
      ss << Util::formatDescription("   radius=30000","Only use values from locations within this radius (in meters). Must be >= 0.") << std::endl;
      ss << Util::formatDescription("   efoldDist=30000","How fast should the weight of a station reduce with distance? For cressman: linearly decrease to this distance (in meters); For barnes: reduce to 1/e after this distance (in meters). Must be >= 0.") << std::endl;
      ss << Util::formatDescription("   maxElevDiff=undef","What is the maximum elevation difference (in meters) that bias can be spread to? Must be >= 0. Leave undefined if no reduction of bias in the vertical is desired.") << std::endl;
      ss << Util::formatDescription("   auxVariable=undef","Should an auxilary variable be used to turn off kriging? For example turn off kriging where there is precipitation.") << std::endl;
      ss << Util::formatDescription("   range=undef","What range of the auxillary variable should kriging be turned on for? For example use 0,0.3 to turn kriging off for precip > 0.3.") << std::endl;
      ss << Util::formatDescription("   window=0","Use a time window to allow weighting of the kriging. Use the fraction of timesteps within +- window where the auxillary variable is within the range. Use 0 for no window.") << std::endl;
      ss << Util::formatDescription("   type=cressman","Weighting function used in kriging. One of 'cressman', or 'barnes'.") << std::endl;
      ss << Util::formatDescription("   operator=add","How should the bias be applied to the raw forecast? One of 'add', 'subtract', 'multiply', 'divide'. For add/subtract, the mean of the field is assumed to be 0, and for multiply/divide, 1.") << std::endl;
      ss << Util::formatDescription("   approxDist=true","When computing the distance between two points, should the equirectangular approximation be used to save time? Should be good enough for most kriging purposes.") << std::endl;
      ss << Util::formatDescription("   crossValidate=false","If true, then don't use the nearest point in the kriging. The end result is a field that can be verified against observations at the kriging points.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c kriging","Spreads bias in space by using kriging")<< std::endl;
   return ss.str();
}
