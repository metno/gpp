#include "Oi.h"
#include "../Util.h"
#include "../Parameters.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <math.h>
#include <armadillo>
#include "Neighbourhood.h"

CalibratorOi::CalibratorOi(Variable iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mVLength(100),
      mHLength(30000),
      mHLengthC(10000),
      mMu(0.9),
      mMinRho(0.0013),
      mEpsilon(0.5),
      mEpsilonC(0.5),
      mElevGradient(-999),
      mBiasVariable(""),
      mSigma(1),
      mSigmaC(1),
      mDelta(1),
      mC(1.03),
      mSaveDiff(false),
      mDeltaVariable(""),
      mUseEns(true),
      // Add mDeltaVariable
      mX(Util::MV),
      mY(Util::MV),
      mMaxLocations(20),
      mNumVariable(""),
      mMaxElevDiff(200),
      // Default model error variance
      mMinValidEns(5),
      mNewDeltaVar(1),
      mExtrapolate(false),
      mDiagnose(false),
      mWMin(0.5),
      mLambda(0.5),
      mCrossValidate(false),
      mMaxBytes(6.0 * 1024 * 1024 * 1024),
      mLandOnly(false),
      mTransformType(TransformTypeNone),
      mDiaFile(""),
      mGamma(0.25),
      mRhoType(RhoTypeGaussian),
      mBoxCoxThreshold(Util::MV) {
   iOptions.getValue("biasVariable", mBiasVariable);
   iOptions.getValue("d", mHLength);
   iOptions.getValue("h", mVLength);
   iOptions.getValue("dc", mHLengthC);
   iOptions.getValue("maxLocations", mMaxLocations);
   iOptions.getValue("sigma", mSigma);
   if(!iOptions.getValue("sigmaC", mSigmaC))
      mSigmaC = mSigma;
   iOptions.getValue("delta", mDelta);
   iOptions.getValue("deltaVariable", mDeltaVariable);
   iOptions.getValue("gamma", mGamma);
   iOptions.getValue("mu", mMu);
   iOptions.getValue("x", mX);
   iOptions.getValue("y", mY);
   iOptions.getValue("extrapolate", mExtrapolate);
   iOptions.getValue("minRho", mMinRho);
   iOptions.getValue("saveDiff", mSaveDiff);
   iOptions.getValue("maxBytes", mMaxBytes);
   iOptions.getValue("minEns", mMinValidEns);
   iOptions.getValue("numVariable", mNumVariable);
   iOptions.getValue("elevGradient", mElevGradient);
   iOptions.getValue("useEns", mUseEns);
   iOptions.getValue("wmin", mWMin);
   iOptions.getValue("epsilon", mEpsilon);
   if(iOptions.getValue("epsilonC", mEpsilonC))
      mEpsilonC = mEpsilon;
   iOptions.getValue("c", mC);
   iOptions.getValue("lambda", mLambda);
   iOptions.getValue("crossValidate", mCrossValidate);
   iOptions.getValue("landOnly", mLandOnly);
   iOptions.getValue("diaFile", mDiaFile);
   iOptions.getValue("diagnose", mDiagnose);
   iOptions.getValue("newDeltaVar", mNewDeltaVar);
   iOptions.getValue("boxCoxThreshold", mBoxCoxThreshold);
   iOptions.getValue("maxElevDiff", mMaxElevDiff);  // Don't use obs that are further than this from their nearest neighbour
   std::string transformType;
   if(iOptions.getValue("transform", transformType)) {
      if(transformType == "boxcox")
         mTransformType = TransformTypeBoxCox;
      else {
         std::stringstream ss;
         ss << "Could not recognize transform=" << transformType << std::endl;
         Util::error(ss.str());
      }
   }
   std::string rhoType;
   if(iOptions.getValue("rhoType", rhoType)) {
      if(rhoType == "gaussian")
         mRhoType = RhoTypeGaussian;
      else if(rhoType == "soar")
         mRhoType = RhoTypeSoar;
      else {
         std::stringstream ss;
         ss << "Could not recognize rhoType=" << rhoType << std::endl;
         Util::error(ss.str());
      }
   }

   iOptions.check();

   // Gamma: The error covariance of the bias is this fraction of the background error
}

// Set up convenient functions for debugging in gdb
template<class Matrix>
void print_matrix(Matrix matrix) {
       matrix.print(std::cout);
}

template void print_matrix<CalibratorOi::mattype>(CalibratorOi::mattype matrix);
template void print_matrix<CalibratorOi::cxtype>(CalibratorOi::cxtype matrix);

bool CalibratorOi::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Single-member SOAR variables
   float sigmaThreshold = 0.001;
   std::vector<CalibratorNeighbourhood> smoothers;
   smoothers.push_back(CalibratorNeighbourhood(Variable(), Options("radius=25 stat=mean fast=0")));
   smoothers.push_back(CalibratorNeighbourhood(Variable(), Options("radius=5 stat=mean fast=0")));
   smoothers.push_back(CalibratorNeighbourhood(Variable(), Options("radius=3 stat=mean fast=0")));
   CalibratorNeighbourhood outputSmoother = CalibratorNeighbourhood(Variable(), Options("radius=3 stat=mean fast=0"));

   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 lafs = iFile.getLandFractions();

   // Check if this method can be applied
   bool hasValidGridpoint = false;
   for(int y = 0; y < nY; y++) {
      for(int x = 0; x < nX; x++) {
         if(Util::isValid(lats[y][x]) && Util::isValid(lons[y][x]) && Util::isValid(elevs[y][x])) {
            hasValidGridpoint = true;
         }
      }
   }
   if(!hasValidGridpoint) {
      Util::warning("There are no gridpoints with valid lat/lon/elev values. Skipping oi...");
      return false;
   }

   // Check parameters
   int maxNumParameters = 3;
   if(iParameterFile->getNumParameters() > maxNumParameters) {
      std::stringstream ss;
      ss << "Parameter file has " << iParameterFile->getNumParameters() << " parameters, which is greater than " << maxNumParameters;
      Util::error(ss.str());
   }

   // Find the grid configuration: regular or irregular
   // For regular grids, we can do certain optimizations
   std::vector<Location> gLocations = iParameterFile->getLocations();
   int gS = gLocations.size();
   float gridSize = Util::MV; // meters between each gridbox
   if(lats.size() > 1 && lats[0].size() > 1) {
      gridSize = Util::getDistance(lats[0][0], lons[0][0], lats[1][0], lons[1][0]);
      std::stringstream ss;
      ss << "Grid size: " << gridSize << " m";
      Util::info(ss.str());
   }
   else {
      std::stringstream ss;
      ss << "Could not determine grid size. Treating as irregular grid.";
      Util::info(ss.str());
   }
   bool isRegularGrid = Util::isValid(gridSize);

   // Log file for diagnostics of rejected observations
   std::ofstream diaFile;
   if(mDiaFile != "") {
     diaFile.open(mDiaFile.c_str());
   }

   // Loop over each observation, find the nearest gridpoint and place the obs into all gridpoints
   // in the vicinity of the nearest neighbour. This is only meant to be an approximation, but saves
   // considerable time instead of doing a loop over each grid point and each observation.

   // Store the indicies (into the gLocations array) that a gridpoint has available
   std::vector<std::vector<std::vector<int> > > gLocIndices; // Y, X, obs indices
   std::vector<float> gYi(gS, Util::MV);
   std::vector<float> gXi(gS, Util::MV);
   std::vector<float> gLafs(gS, Util::MV);
   std::vector<float> gElevs(gS, Util::MV);
   std::vector<float> gCi(gS, 1);
   std::vector<float> gRadarL(gS, 0);
   std::vector<float> gObs(gS, Util::MV);
   gLocIndices.resize(nY);
   for(int y = 0; y < nY; y++) {
      gLocIndices[y].resize(nX);
   }

   // Calculate the factor that the horizontal decorrelation scale should be multiplied by
   // to get the localization radius. For minRho=0.0013, the factor is 3.64
   float radiusFactor = sqrt(-2*log(mMinRho));

   // Spread each observation out to this many gridpoints from the nearest neighbour

   // Check that we do not run out of memory
   int gridpointRadius = Util::MV;
   if(isRegularGrid) {
      gridpointRadius = radiusFactor * mHLength / gridSize;

      // When large radiuses are used, the process becomes memory-intensive. Try to fail here
      // if we expect to use more memory than desired. The true memory is roughly
      // 1 GB + expectedBytes * F
      int bytesPerValue = 4;
      float expectedBytes = float(gridpointRadius * gridpointRadius) * 4 * bytesPerValue * gS;
      std::cout << "Expected MB: " << 1000 + expectedBytes / 1024 / 1024 << std::endl;
      if(Util::isValid(mMaxBytes) && expectedBytes > mMaxBytes) {
         std::stringstream ss;
         ss << "Number of required bytes (" << expectedBytes / 1024 / 1024 << " MB) exceeds threshold "
            << float(mMaxBytes) / 1024 / 1024 << " MB. Use maxBytes= to set threshold.";
         Util::error(ss.str());
      }
   }

   int numParameters = iParameterFile->getNumParameters();

   // For each gridpoint, find which observations are relevant. Parse the observations and only keep
   // those that pass certain checks
   double time_s = Util::clock();
   KDTree searchTree(iFile.getLats(), iFile.getLons());
   for(int i = 0; i < gS; i++) {
      if(i % 1000 == 0) {
         std::stringstream ss;
         ss << i;
         Util::info(ss.str());
      }
      Parameters parameters = iParameterFile->getParameters(0, gLocations[i]);
      assert(parameters.size() > 0);
      gObs[i] = transform(parameters[0]);
      if(parameters.size() == 2)
         gCi[i] = parameters[1];
      if(parameters.size() == 3)
         gRadarL[i] = parameters[2];

      gElevs[i] = gLocations[i].elev();
      int Y, X;
      searchTree.getNearestNeighbour(gLocations[i].lat(), gLocations[i].lon(), Y, X);
      gYi[i] = Y;
      gXi[i] = X;
      gLafs[i] = lafs[Y][X];

      if(!Util::isValid(gridSize) || (X > 0 && X < nX - 1 && Y > 0 && Y < nY - 1)) {
         if(Util::isValid(gObs[i])) {
            // Check if the elevation of the station roughly matches the reference grid elevation
            bool hasValidElev = !Util::isValid(mMaxElevDiff) || Util::isValid(gLocations[i].elev());
            if(Util::isValid(mMaxElevDiff) && hasValidElev) {
               float elevDiff = abs(gLocations[i].elev() - elevs[Y][X]);
               hasValidElev = elevDiff < mMaxElevDiff;
            }
            if(hasValidElev) {

               // Don't include an observation if it is in the ocean and landOnly=1
               bool wrongLaf = Util::isValid(gLafs[i]) && mLandOnly && gLafs[i] == 0;
               if(!wrongLaf) {
                  if(isRegularGrid) {
                     for(int y = std::max(0, Y - gridpointRadius); y < std::min(nY, Y + gridpointRadius); y++) {
                        for(int x = std::max(0, X - gridpointRadius); x < std::min(nX, X + gridpointRadius); x++) {
                           gLocIndices[y][x].push_back(i);
                        }
                     }
                  }
                  else {
                     // Do a rough localization based on distance
                     for(int y = 0; y < nY; y++) {
                        for(int x = 0; x < nX; x++) {
                           float dist = Util::getDistance(gLocations[i].lat(), gLocations[i].lon(), lats[y][x], lons[y][x], true);
                           if(dist < radiusFactor * mHLength)
                              gLocIndices[y][x].push_back(i);
                        }
                     }
                  }
               }
               else {
                  std::stringstream ss;
                  ss << "Removing station (" << gLocations[i].lat() << " " << gLocations[i].lon() << ") because laf=0";
                  Util::warning(ss.str());
                  if(mDiaFile != "") {
                     diaFile << gLocations[i].lon() << ";" << gLocations[i].lat() << ";2;" << std::endl;  
                  }
               }
            }
            else {
               std::stringstream ss;
               ss << "Removing station (" << gLocations[i].lat() << " " << gLocations[i].lon() << ") because elevation (" << gLocations[i].elev() << " m) is too far from grid (" << elevs[Y][X] << " m)";
               Util::warning(ss.str());
               if(mDiaFile != "") {
                  diaFile << gLocations[i].lon() << ";" << gLocations[i].lat() << ";1;" << std::endl;  
               }
            }
         }
      }
      else {
         std::stringstream ss;
         ss << "Removing station (" << gLocations[i].lat() << " " << gLocations[i].lon() << ") because outside domain";
         Util::warning(ss.str());
      }
   }
   double time_e = Util::clock();
   std::cout << "Assigning locations " << time_e - time_s << std::endl;


   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      FieldPtr output = iFile.getEmptyField();
      FieldPtr bias;
      FieldPtr newbias;
      FieldPtr delta;
      FieldPtr newdelta;

      // Transform the background
      #pragma omp parallel for
      for(int x = 0; x < nX; x++) {
         for(int y = 0; y < nY; y++) {
            for(int e = 0; e < nEns; e++) {
               float value = (*field)(y, x, e);
               if(Util::isValid(value))
                  (*field)(y, x, e) = transform(value);
            }
         }
      }

      bool useBias = mBiasVariable != "";
      if(useBias) {
         bias = iFile.getField(mBiasVariable, t);

         // Scale the bias
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               if(Util::isValid((*bias)(y, x, 0)))
                  (*bias)(y, x, 0) *= mMu;
               else
                  (*bias)(y, x, 0) = 0;
            }
         }
         newbias = iFile.getEmptyField(0);
      }
      bool useDelta = mDeltaVariable != "";
      if(useDelta) {
         delta = iFile.getField(mDeltaVariable, t);
         // Initialize if missing
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               if(!Util::isValid((*delta)(y, x, 0)))
                  (*delta)(y, x, 0) = 1;
            }
         }
         newdelta = iFile.getEmptyField(0);
      }

      FieldPtr num;
      if(mNumVariable != "") {
         num = iFile.getField(mNumVariable, t);
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               for(int e = 0; e < nEns; e++) {
                  (*num)(y, x, e) = 0;
               }
            }
         }
      }

      // Compute Y
      vec2 gY(gS);
      std::vector<float> gYhat(gS);
      for(int i = 0; i < gS; i++) {
         gY[i].resize(nEns, 0);
         float elevCorr = 0;
         if(Util::isValid(mElevGradient) && mElevGradient != 0) {
            float nnElev = elevs[gYi[i]][gXi[i]];
            assert(Util::isValid(nnElev));
            assert(Util::isValid(gLocations[i].elev()));
            float elevDiff = gLocations[i].elev() - nnElev;
            elevCorr = mElevGradient * elevDiff;
         }
         float total = 0;
         int count = 0;
         for(int e = 0; e < nEns; e++) {
            float value = (*field)(gYi[i], gXi[i], e);
            if(Util::isValid(value)) {
               value += elevCorr;
               gY[i][e] = value;
               total += value;
               count++;
            }
         }
         float mean = total / count;
         for(int e = 0; e < nEns; e++) {
            float value = gY[i][e];
            if(Util::isValid(value)) {
               gY[i][e] -= mean;
            }
         }
         gYhat[i] = mean;
         if(useBias) {
            float currBias = (*bias)(gYi[i], gXi[i], 0);
            gYhat[i] -= currBias;
         }
      }

      // Calculate number of valid members
      int nValidEns = 0;
      std::vector<int> validEns;
      for(int e = 0; e < nEns; e++) {
         int numInvalid = 0;
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               float value = (*field)(y, x, e);
               if(!Util::isValid(value))
                  numInvalid++;
            }
         }
         if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
         }
      }
      std::cout << "Number of valid ensemble members: " << nValidEns << std::endl;
      bool singleMemberMode = !mUseEns || nValidEns < mMinValidEns;

      // Temporary field for single-member mode when using a transform
      FieldPtr sigmaTransformed = iFile.getEmptyField(0);

      #pragma omp parallel for
      for(int x = 0; x < nX; x++) {
         for(int y = 0; y < nY; y++) {
            float lat = lats[y][x];
            float lon = lons[y][x];
            float elev = elevs[y][x];
            float laf = lafs[y][x];

            //
            // Create list of locations for this gridpoint
            //
            std::vector<int> lLocIndices0 = gLocIndices[y][x];
            std::vector<int> lLocIndices;
            lLocIndices.reserve(lLocIndices0.size());
            std::vector<std::pair<float,int> > lRhos0;
            lRhos0.reserve(lLocIndices0.size());
            for(int i = 0; i < lLocIndices0.size(); i++) {
               int index = lLocIndices0[i];
               float hdist = Util::getDistance(gLocations[index].lat(), gLocations[index].lon(), lat, lon, true);
               float vdist = Util::MV;
               if(Util::isValid(gLocations[index].elev() && Util::isValid(elev)))
                  vdist = gLocations[index].elev() - elev;
               float lafdist = 0;
               if(Util::isValid(gLafs[index]) && Util::isValid(laf))
                  lafdist = gLafs[index] - laf;
               float rho = calcRho(hdist, vdist, lafdist, mRhoType);
               int X = gXi[index];
               int Y = gYi[index];
               // Only include observations that are within the domain
               if(!isRegularGrid || (X > 0 && X < lats[0].size()-1 && Y > 0 && Y < lats.size()-1)) {
                  if(rho > mMinRho) {
                     lRhos0.push_back(std::pair<float,int>(rho, i));
                  }
               }
            }

            // Don't include the best value if cross-validating
            if(mCrossValidate) {
               float maxRho = 0;
               float maxRhoIndex = Util::MV;
               for(int i = 0; i < lRhos0.size(); i++) {
                  if(lRhos0[i].first > maxRho) {
                     maxRho = lRhos0[i].first;
                     maxRhoIndex = i;
                  }
               }
               if(Util::isValid(maxRhoIndex)) {
                  int ii = lRhos0[maxRhoIndex].second;
                  int index = lLocIndices0[ii];
                  float obsLat = gLocations[index].lat();
                  float obsLon = gLocations[index].lon();
                  lRhos0.erase(lRhos0.begin() + maxRhoIndex);
                  std::stringstream ss;
                  ss << "Omitting " << obsLat << "," << obsLon
                     << " for forecast point " << lat << "," << lon << " due to cross-validation";
                  Util::info(ss.str());
               }
            }
            arma::vec lRhos;
            if(lRhos0.size() > mMaxLocations) {
               // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
               // Otherwise, just use the last locations added
               lRhos = arma::vec(mMaxLocations);
               std::sort(lRhos0.begin(), lRhos0.end(), Util::sort_pair_first<float,int>());
               for(int i = 0; i < mMaxLocations; i++) {
                  // The best values start at the end of the array
                  int index = lRhos0[lRhos0.size() - 1 - i].second;
                  lLocIndices.push_back(lLocIndices0[index]);
                  lRhos(i) = lRhos0[lRhos0.size() - 1 - i].first;
               }
            }
            else {
               lRhos = arma::vec(lRhos0.size());
               for(int i = 0; i < lRhos0.size(); i++) {
                  int index = lRhos0[i].second;
                  lLocIndices.push_back(lLocIndices0[index]);
                  lRhos(i) = lRhos0[i].first;
               }
            }

            int lS = lLocIndices.size();
            if(x == mX && y == mY) {
               std::cout << "Number of local stations: " << lS << std::endl;
            }

            if(lS == 0) {
               // If we have too few observations though, then use the background
               for(int e = 0; e < nEns; e++) {
                  if(mSaveDiff)
                     (*output)(y, x, e) = Util::MV;
                  else
                     (*output)(y, x, e) = (*field)(y, x, e);
               }
               continue;
            }

            vectype lObs(lS);
            vectype lElevs(lS);
            vectype lLafs(lS);
            for(int i = 0; i < lLocIndices.size(); i++) {
               int index = lLocIndices[i];
               lObs[i] = gObs[index];
               lElevs[i] = gElevs[index];
               lLafs[i] = gLafs[index];
            }

            // Compute Y (model at obs-locations)
            mattype lY(lS, nValidEns);
            vectype lYhat(lS);

            for(int i = 0; i < lS; i++) {
               // Use the nearest neighbour for this location
               int index = lLocIndices[i];
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  lY(i, e) = gY[index][ei];
               }
               lYhat(i) = gYhat[index];

            }

            ////////////////////////////////////////////////////////////////////////////////////////
            // Single-member mode:                                                                //
            // Revert to static structure function when there is not enough ensemble information  //
            ////////////////////////////////////////////////////////////////////////////////////////
            if(singleMemberMode) {
               // Current grid-point to station error covariance matrix
               mattype lG(1, lS, arma::fill::zeros);
               // Station to station error covariance matrix
               mattype lP(lS, lS, arma::fill::zeros);
               // Station variance
               mattype lR(lS, lS, arma::fill::zeros);
               for(int i = 0; i < lS; i++) {
                  int index = lLocIndices[i];
                  lR(i, i) = gCi[index];
                  float hdist = Util::getDistance(gLocations[index].lat(), gLocations[index].lon(), lat, lon, true);
                  float vdist = Util::MV;
                  if(Util::isValid(gLocations[index].elev() && Util::isValid(elev)))
                     vdist = gLocations[index].elev() - elev;
                  float lafdist = 0;
                  if(Util::isValid(gLafs[index]) && Util::isValid(laf))
                     lafdist = gLafs[index] - laf;
                  float rho = calcRho(hdist, vdist, lafdist, mRhoType);
                  lG(0, i) = rho;
                  for(int j = 0; j < lS; j++) {
                     int index_j = lLocIndices[j];
                     float hdist = Util::getDistance(gLocations[index].lat(), gLocations[index].lon(), gLocations[index_j].lat(), gLocations[index_j].lon(), true);
                     float vdist = Util::MV;
                     if(Util::isValid(gLocations[index].elev() && Util::isValid(gLocations[index_j].elev())))
                        vdist = gLocations[index].elev() - gLocations[index_j].elev();
                           float lafdist = 0;
                     if(Util::isValid(gLafs[index]) && Util::isValid(laf))
                        lafdist = gLafs[index] - gLafs[index_j];

                     lP(i, j) = calcRho(hdist, vdist, lafdist, mRhoType);
                  }
               }
               mattype lGSR;  // Kalman gain
               mattype lSRinv;
               // TODO: This will be different for precipitation
               if(useBias) {
                  lSRinv = arma::inv(lP + 1 / (1 + mGamma) * mEpsilon * mEpsilon * lR);
                  lGSR = lG * lSRinv;
               }
               else {
                  lSRinv = arma::inv(lP + mEpsilon * mEpsilon * lR);
                  lGSR = lG * lSRinv;
               }

               // This should loop over nValidEns. And use ei.
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  if (Util::isValid((*field)(y, x, ei))) {
                     vectype currFcst = lY.col(e) + lYhat;
                     vectype dx = lGSR * (lObs - currFcst);

                     // Store sigma in transformed space
                     (*output)(y, x, ei) = (*field)(y, x, ei) + dx[0];
                     if(mTransformType != TransformTypeNone) {
                        if( (*output)(y, x, ei)< -1.0 / mLambda) {
                           (*output)(y, x, ei) = -1.0 / mLambda;
                        }
                        if( (*output)(y, x, ei) >= transform(mBoxCoxThreshold)) {
                           vectype incrementAtObsPoints = lP * (lSRinv * (lObs - currFcst));
                           float total = 0;
                           float totalDiagR = 0;
                           float lGSRG = 0;
                           for(int s = 0; s < lS; s++) {
                              total += (lObs[s] - currFcst[s]) * (lObs[s] - currFcst[s] - incrementAtObsPoints[s]);
                              // CL old: totalDiagR += mEpsilon * mEpsilon * lR[s];
                              totalDiagR += mEpsilon * mEpsilon * lR(s,s);
                              lGSRG += lGSR[s] * lG[s];
                           }
                           float sigmaObs = total / lS;
                           float meanDiagR = totalDiagR / lS;
                           float sigmaB = sigmaObs / meanDiagR;
                           (*sigmaTransformed)(y, x, ei) = std::max(sigmaThreshold,sigmaB * (1 - lGSRG));
                           if(x == mX && y == mY) {
                              std::cout << "sigmaObs: " << sigmaObs << std::endl;
                              std::cout << "meanDiagR: " << meanDiagR << std::endl;
                              std::cout << "sigmaB: " << sigmaB << std::endl;
                              std::cout << "sigmaTransformed: " << (*sigmaTransformed)(y, x, ei) << std::endl;
                           }
                        }
                     }
                     if(x == mX && y == mY) {
                        std::cout << "Lat: " << lat << std::endl;
                        std::cout << "Lon: " << lon << " " << lat << " " << std::endl;
                        std::cout << "Elev: " << elev << std::endl;
                        std::cout << "P:" << std::endl;
                        print_matrix<mattype>(lP);
                        std::cout << "R:" << std::endl;
                        print_matrix<mattype>(lR);
                        std::cout << "GSR:" << std::endl;
                        print_matrix<mattype>(lGSR);
                        std::cout << "Obs:" << std::endl;
                        print_matrix<mattype>(lObs);
                        std::cout << "Current forecast: " << std::endl;
                        print_matrix<mattype>(currFcst);
                        std::cout << "Increment" << std::endl;
                        print_matrix<mattype>(lObs - currFcst);
                        std::cout << "Yhat" << std::endl;
                        print_matrix<mattype>(lYhat);
                        std::cout << "dx: " << dx[0] << std::endl;
                     }
                  }
               }

               // Update bias
               if(useBias) {
                  float biasTotal = 0;
                  (*newbias)(y, x, 0) = (*bias)(y, x, 0) - mGamma / (1 + mGamma) * biasTotal;
               }
            }
            ////////////////////////////////////////////////////////////////////////////////////////
            // Ensemble-member mode:                                                              //
            // Use ensemble covariance structure                                                  //
            ////////////////////////////////////////////////////////////////////////////////////////
            else {
               // Compute Rinv
               mattype Rinv(lS, lS, arma::fill::zeros);
               if(numParameters == 2) {
                  for(int i = 0; i < lS; i++) {
                     int index = lLocIndices[i];
                     Rinv(i, i) = lRhos[i] / (mSigma * mSigma * gCi[index]);
                     if(x == mX && y == mY) {
                        std::cout << "R(" << i << ") " << Rinv(i, i) << std::endl;
                     }
                  }
               }
               else if(numParameters == 3) {
                  // Inverting the matrix is more complicated, since the radar observations
                  // have covariances. Therefore invert the covariance matrix for the radar part and
                  // insert the values into the big inverse matrix.
                  // std::cout << "Computing R matrix" << std::endl;
                  // R = get_precipitation_r(gRadarL, gCi, lLocIndices, lRhos);
                  // Compute little R
                  std::vector<int> gRadarIndices;
                  gRadarIndices.reserve(lS);
                  std::vector<int> lRadarIndices;
                  lRadarIndices.reserve(lS);
                  for(int i = 0; i < lS; i++) {
                     int index = lLocIndices[i];
                     if(gRadarL[index] > 0) {
                        gRadarIndices.push_back(index);
                        lRadarIndices.push_back(i);
                     }
                  }
                  int lNumRadar = gRadarIndices.size();

                  // Compute R tilde r
                  mattype radarR(lNumRadar, lNumRadar, arma::fill::zeros);
                  for(int i = 0; i < lNumRadar; i++) {
                     for(int j = 0; j < lNumRadar; j++) {
                        int gIndex_i = gRadarIndices[i];
                        int gIndex_j = gRadarIndices[j];
                        int lIndex_i = lRadarIndices[i];
                        int lIndex_j = lRadarIndices[j];
                        if(i == j) {
                           radarR(i, i) = 1;
                        }
                        else {
                           // Equation 5
                           float dist = Util::getDistance(gLocations[gIndex_i].lat(), gLocations[gIndex_i].lon(), gLocations[gIndex_j].lat(), gLocations[gIndex_j].lon(), true);
                           float h = dist / mHLengthC;
                           float rho = (1 + h) * exp(-h);
                           radarR(i, j) = rho;
                        }
                     }
                  }
                  if(x == mX && y == mY) {
                     std::cout << "Number of radar points: " << " " << lNumRadar << std::endl;
                     if(lNumRadar > 0) {
                        print_matrix<mattype>(radarR);
                     }
                  }

                  float cond = arma::rcond(radarR);
                  if(cond <= 0) {
                     std::stringstream ss;
                     ss << "Condition number of " << cond << " for radar values. Using raw values";
                     Util::warning(ss.str());
                     for(int e = 0; e < nEns; e++) {
                        (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
                     }
                     continue;
                  }

                  mattype radarRinv(lNumRadar, lNumRadar, arma::fill::zeros); 
                  radarRinv = arma::inv(radarR);

                  for(int i = 0; i < lS; i++) {
                     int index = lLocIndices[i];
                     Rinv(i, i) = lRhos[i] / (mSigma * mSigma * gCi[index]);
                  }
                  // Overwrite where we have radar pixels
                  for(int i = 0; i < lNumRadar; i++) {
                     int ii = lRadarIndices[i];
                     for(int j = 0; j < lNumRadar; j++) {
                        int jj = lRadarIndices[j];
                        Rinv(ii, jj) = sqrt(lRhos[ii] * lRhos[jj]) / (mSigmaC * mSigmaC) * radarRinv(i, j);
                     }
                  }
               }
               else {
                  abort();
               }

               // Compute C matrix
               // k x gS * gS x gS
               mattype C(nValidEns, lS);
               C = lY.t() * Rinv;

               mattype Pinv(nValidEns, nValidEns);
               float currDelta = 1;
               if(useDelta)
                  currDelta = (*delta)(y, x, 0);
               else
                  currDelta = mDelta;
               float diag = 1 / currDelta * (nValidEns - 1);
               if(useBias)
                  diag = 1 / currDelta / (1 + mGamma) * (nValidEns - 1);

               Pinv = C * lY + diag * arma::eye<mattype>(nValidEns, nValidEns);
               float cond = arma::rcond(Pinv);
               if(cond <= 0) {
                  std::stringstream ss;
                  ss << "Condition number of " << cond << ". Using raw values";
                  Util::warning(ss.str());
                  for(int e = 0; e < nEns; e++) {
                     (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
                  }
                  continue;
               }

               // Compute sqrt of matrix. Armadillo 6.6 has this function, but on many systems this
               // is not available. Therefore, compute sqrt using the method found in 6.6
               // cxtype Wcx(nValidEns, nValidEns);
               // status = arma::sqrtmat(Wcx, (nValidEns - 1) * P);
               // mattype W = arma::real(Wcx);

               mattype P = arma::inv(Pinv);
               vectype eigval;
               mattype eigvec;
               bool status = arma::eig_sym(eigval, eigvec, (nValidEns - 1) * P);
               if(!status) {
                  std::cout << "Cannot find eigenvector:" << std::endl;
                  std::cout << "Lat: " << lat << std::endl;
                  std::cout << "Lon: " << lon << std::endl;
                  std::cout << "Elev: " << elev << std::endl;
                  std::cout << "Laf: " << laf << std::endl;
                  std::cout << "Pinv" << std::endl;
                  print_matrix<mattype>(Pinv);
                  std::cout << "P" << std::endl;
                  print_matrix<mattype>(P);
                  std::cout << "Y:" << std::endl;
                  print_matrix<mattype>(lY);
                  std::cout << "lObs:" << std::endl;
                  print_matrix<mattype>(lObs);
                  std::cout << "Yhat" << std::endl;
                  print_matrix<mattype>(lYhat);
               }
               eigval = sqrt(eigval);
               mattype Wcx = eigvec * arma::diagmat(eigval) * eigvec.t();
               mattype W = arma::real(Wcx);

               if(W.n_rows == 0) {
                  std::stringstream ss;
                  ss << "Could not find the real part of W. Using raw values.";
                  Util::warning(ss.str());
                  for(int e = 0; e < nEns; e++) {
                     (*output)(y, x, e) = (*field)(y, x, e);
                  }
                  continue;
               }

               // Compute PC
               mattype PC(nValidEns, lS);
               PC = P * C;

               // Compute w
               vectype w(nValidEns);
               if(mDiagnose)
                  w = PC * (arma::ones<vectype>(lS));
               else
                  w = PC * (lObs - lYhat);

               // Add w to W
               for(int e = 0; e < nValidEns; e++) {
                  for(int e2 = 0; e2 < nValidEns; e2 ++) {
                     W(e, e2) = W(e, e2) + w(e) ;
                  }
               }

               // Compute X (perturbations about model mean)
               vectype X(nValidEns);
               float total = 0;
               int count = 0;
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  float value = (*field)(y, x, ei);
                  if(Util::isValid(value)) {
                     X(e) = value;
                     total += value;
                     count++;
                  }
                  else {
                     std::cout << "Invalid value " << y << " " << x << " " << e << std::endl;
                  }
               }
               float ensMean = total / count;
               for(int e = 0; e < nValidEns; e++) {
                  X(e) -= ensMean;
               }

               // Write debugging information
               if(x == mX && y == mY) {
                  std::cout << "Lat: " << lat << std::endl;
                  std::cout << "Lon: " << lon << " " << lat << " " << std::endl;
                  std::cout << "Elev: " << elev << std::endl;
                  std::cout << "Laf: " << laf << std::endl;
                  std::cout << "Num obs: " << lS << std::endl;
                  std::cout << "Num ens: " << nValidEns << std::endl;
                  std::cout << "rhos" << std::endl;
                  print_matrix<mattype>(lRhos);
                  std::cout << "P" << std::endl;
                  print_matrix<mattype>(P);
                  std::cout << "C" << std::endl;
                  print_matrix<mattype>(C);
                  std::cout << "C * lY" << std::endl;
                  print_matrix<mattype>(C * lY);
                  std::cout << "PC" << std::endl;
                  print_matrix<mattype>(PC);
                  std::cout << "W" << std::endl;
                  print_matrix<mattype>(W);
                  std::cout << "w" << std::endl;
                  print_matrix<mattype>(w);
                  std::cout << "Y:" << std::endl;
                  print_matrix<mattype>(lY);
                  std::cout << "Yhat" << std::endl;
                  print_matrix<mattype>(lYhat);
                  std::cout << "lObs" << std::endl;
                  print_matrix<mattype>(lObs);
                  std::cout << "lObs - Yhat" << std::endl;
                  print_matrix<mattype>(lObs - lYhat);
                  std::cout << "X" << std::endl;
                  print_matrix<mattype>(X);
                  std::cout << "elevs" << std::endl;
                  print_matrix<mattype>(lElevs);
                  std::cout << "lafs" << std::endl;
                  print_matrix<mattype>(lLafs);
                  std::cout << "Analysis increment:" << std::endl;
                  print_matrix<mattype>(X.t() * W);
                  std::cout << "My: " << arma::mean(arma::dot(lObs - lYhat, lRhos) / lS) << std::endl;
               }

               // Compute analysis
               for(int e = 0; e < nValidEns; e++) {
                  int ei = validEns[e];
                  float total = 0;
                  for(int k = 0; k < nValidEns; k++) {
                     total += X(k) * W(k, e);
                  }

                  float currIncrement = total;

                  if(mSaveDiff)
                     (*output)(y, x, ei) = currIncrement;
                  else {
                     float raw = ensMean;
                     if(useBias) {
                        raw -= (*bias)(y, x, 0);
                     }

                     ///////////////////////////////
                     // Anti-extrapolation filter //
                     ///////////////////////////////
                     if(!mExtrapolate) {
                        // Don't allow a final increment that is larger than any increment
                        // at station points
                        float maxInc = arma::max(lObs - (lY[e] + lYhat));
                        float minInc = arma::min(lObs - (lY[e] + lYhat));
                        if(x == mX && y == mY) {
                           std::cout << "Increments: " << maxInc << " " << minInc << " " << currIncrement << std::endl;
                        }

                        // The increment for this member. currIncrement is the increment relative to
                        // ensemble mean
                        float memberIncrement = currIncrement - X(e);
                        // Adjust increment if it gives a member increment that is outside the range
                        // of the observation increments
                        if(x == mX && y == mY) {
                           std::cout << "Analysis increment: " << memberIncrement << " " << ensMean << " " << currIncrement << " " << X(e) << std::endl;
                        }
                        if(maxInc > 0 && memberIncrement > maxInc) {
                           currIncrement = maxInc + X(e);
                        }
                        else if(maxInc < 0 && memberIncrement > 0) {
                           currIncrement = 0 + X(e);
                        }
                        else if(minInc < 0 && memberIncrement < minInc) {
                           currIncrement = minInc + X(e);
                        }
                        else if(minInc > 0 && memberIncrement < 0) {
                           currIncrement = 0 + X(e);
                        }
                        if(x == mX && y == mY) {
                           std::cout << "Final increment: " << currIncrement << " " << currIncrement - X(e) << std::endl;
                        }
                     }
                     (*output)(y, x, ei) = ensMean + currIncrement;
                  }

                  if(mNumVariable != "") {
                     (*num)(y, x, ei) = lS;
                  }
               }

               // Update bias
               if(useBias) {
                  float biasTotal = 0;
                  for(int e = 0; e < nValidEns; e++) {
                     int ei = validEns[e];
                     biasTotal += (*field)(y, x, ei) * w(e);
                  }
                  (*newbias)(y, x, 0) = (*bias)(y, x, 0) - mGamma / (1 + mGamma) * biasTotal;
               }

               // Update delta
               /*
               float deltaVar = mC - 1;
               float trace = arma::trace(lY * lY.t());
               float numerator = mSigma * mSigma / mEpsilon / mEpsilon;
               float denomenator = 1.0 / lS / (nValidEns - 1) * trace;
               float currDeltaEvidence = numerator / denomenator;
               float weightOld = deltaVar;
               float weightNew = mNewDeltaVar;
               (*newdelta)(y, x, 0) = ((*delta)(y, x, 0) * weightNew + currDeltaEvidence * weightOld) / (weightOld + weightNew);
               */
            }
         }
      }

      // Back-transform
      Util::info("Back transform");
      if(singleMemberMode) {
         // Single member mode needs to deal with variance
         if(mTransformType != TransformTypeNone) {
            // Smothen sigmaTranformed several times
            for(int r = 0; r < smoothers.size(); r++) {
               #pragma omp parallel for
               for(int x = 0; x < nX; x++) {
                  for(int y = 0; y < nY; y++) {
                     for(int e = 0; e < nEns; e++) {
                       if( (*output)(y, x, e) < transform(mBoxCoxThreshold)) {
                          (*sigmaTransformed)(y, x, e)=0;
                       }
                     }
                  }
               }
               smoothers[r].calibrateField(*sigmaTransformed, *sigmaTransformed);
            }

            #pragma omp parallel for
            for(int x = 0; x < nX; x++) {
               for(int y = 0; y < nY; y++) {
                  for(int e = 0; e < nEns; e++) {
                     float value = (*output)(y, x, e);
                     float sigmaValue = (*sigmaTransformed)(y, x, e);
                     // (*output)(y, x, e) = (*sigmaTransformed)(y, x, e);
                     if(Util::isValid(value)) {
                        if(sigmaValue==0 || value<transform(mBoxCoxThreshold)) {
                          (*output)(y, x, e) = invTransform(value);
                        }
                        else {
                           float f = pow(mLambda * value + 1, 1 / mLambda);
                           float f2 = pow((1 - mLambda) * (mLambda * value + 1), 1 / mLambda - 2);
                           if(x == mX && y == mY) {
                              // TODO: Is it * f2?
                              std::cout <<  (*sigmaTransformed)(y, x, e) << " " <<  (*output)(y, x, e)
                                 << " " << f << " " << f2 <<  " " << f + 0.5 * (*sigmaTransformed)(y, x, e) *
                                 f2 << std::endl;
                           }
                           (*output)(y, x, e) = f + 0.5 * sigmaValue * f2;
                        }
                     }
                  }
               }
            } // end for Backtransf
            // Smooth the output field
            outputSmoother.calibrateField(*output, *output);
         }
      }
      else {
         #pragma omp parallel for
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               for(int e = 0; e < nEns; e++) {
                  float value = (*output)(y, x, e);
                  if(Util::isValid(value))
                     (*output)(y, x, e) = invTransform(value);
               }
            }
         }
      }


      iFile.addField(output, mVariable, t);
      if(mNumVariable != "")
         iFile.addField(num, Variable(mNumVariable), t);
      if(useBias) {
         iFile.addField(newbias, Variable(mBiasVariable), t);
      }
      if(useDelta) {
         float oldDelta = (*delta)(0, 0, 0);
         float value = calcDelta(oldDelta, gY);
         for(int x = 0; x < nX; x++) {
            for(int y = 0; y < nY; y++) {
               (*newdelta)(y, x, 0) = value;
            }
         }

         iFile.addField(newdelta, Variable(mDeltaVariable), t);
      }
   }
   if(mDiaFile != "") {
     diaFile.close();
   }
   return true;
}

float CalibratorOi::calcDelta(float iOldDelta, const vec2& iY) const {
   float deltaVar = mC - 1;
   float trace = 0;
   float numValidS = 0;
   int S = iY.size();
   int nEns = iY[0].size();
   for(int s = 0; s < S; s++) {
      // Compute value in coord s,s
      float value = 0;
      int count = 0;
      for(int e = 0; e < nEns; e++) {
         if(Util::isValid(iY[s][e])) {
            value += iY[s][e] * iY[s][e];
            count++;
         }
      }
      if(count > 1) {
         value = value / (count-1);
         trace += value;
         numValidS++;
      }
   }
   float numerator = mSigma * mSigma / mEpsilon / mEpsilon;
   float denomenator = 1.0 / numValidS * trace;
   float currDeltaEvidence = numerator / denomenator;
   float weightOld = deltaVar;
   float weightNew = mNewDeltaVar;
   return (iOldDelta * weightNew + currDeltaEvidence * weightOld) / (weightOld + weightNew);
}

float CalibratorOi::calcRho(float iHDist, float iVDist, float iLDist, RhoType iType) const {
   float h = (iHDist/mHLength);
   float rho = 1;
   if(iType == RhoTypeGaussian)
      rho = exp(-0.5 * h * h);
   else
      rho = (1 + h) * exp(-h);
   if(Util::isValid(mVLength)) {
      if(!Util::isValid(iVDist)) {
         rho = 0;
      }
      else {
         float v = (iVDist/mVLength);
         rho *= exp(-0.5 * v * v);
      }
   }
   if(Util::isValid(mWMin)) {
      rho *= 1 - (1 - mWMin) * std::abs(iLDist);
   }
   return rho;
}

float CalibratorOi::transform(float iValue) const {
   if(mTransformType == TransformTypeNone)
      return iValue;

   if(iValue <= 0)
      iValue = 0;
   if(mLambda == 0)
      return log(iValue);
   else
      return (pow(iValue, mLambda) - 1) / mLambda;
}

float CalibratorOi::invTransform(float iValue) const {
   if(mTransformType == TransformTypeNone)
      return iValue;

   float rValue = 0;
   if(mLambda == 0)
      rValue = exp(iValue);
   else {
      if(iValue < -1.0 / mLambda) {
         iValue = -1.0 / mLambda;
      }
      rValue = pow(1 + mLambda * iValue, 1 / mLambda);
   }
   if(rValue <= 0)
      rValue = 0;
   return rValue;
}

std::string CalibratorOi::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c oi","Merge observations from parameter file with background field using Optimal interpolation. A parameter file is required where the first parameter is an observation.")<< std::endl;
      ss << Util::formatDescription("   transform=none","One of 'none', 'boxcox'. Use 'boxcox' for 'precipitation'.") << std::endl;
      ss << Util::formatDescription("   d=30000","Horizontal decorrelation distance (in meters). Must be >= 0.") << std::endl;
      ss << Util::formatDescription("   h=100","Vertical decorrelation distance (in meters). Use -999 to disable.") << std::endl;
      ss << Util::formatDescription("   dc=10000","Decorrelation distance (in meters) for observations with spatial correlation.") << std::endl;
      ss << Util::formatDescription("   maxLocations=20","Don't use more than this many locations within the localization region. Sort the stations by rho and use the best ones.") << std::endl;
      ss << Util::formatDescription("   sigma=1","Average standard error of observations") << std::endl;
      ss << Util::formatDescription("   sigmaC=1","Average standard error of observations with spatial correlation") << std::endl;
      ss << Util::formatDescription("   delta=1","Background variance inflation factor") << std::endl;
      // ss << Util::formatDescription("   gamma=0.25","") << std::endl;
      // ss << Util::formatDescription("   mu=0.9","") << std::endl;
      ss << Util::formatDescription("   x=undef","Turn on debug info for this x-coordinate") << std::endl;
      ss << Util::formatDescription("   y=undef","Turn on debug info for this y-coordinate") << std::endl;
      ss << Util::formatDescription("   extrapolate=0","Allow OI to extrapolate increments. If 0, then increments are bounded by the increments at the observation sites.") << std::endl;
      ss << Util::formatDescription("   minRho=0.0013","Perform localization by requiring this minimum rho value") << std::endl;
      ss << Util::formatDescription("   maxBytes=6442450944","Don't allocate more than this many bytes when creating the localization information") << std::endl;
      ss << Util::formatDescription("   minEns=5","Switch to single-member mode if fewer than this number of members available") << std::endl;
      ss << Util::formatDescription("   elevGradient=0","Elevation gradient when downscaling background to obs. Use -0.0065 for temperature.") << std::endl;
      ss << Util::formatDescription("   useEns=1","Enable ensemble-mode. If 0, use single-member mode.") << std::endl;
      ss << Util::formatDescription("   wmin=0.5","") << std::endl;
      ss << Util::formatDescription("   epsilon=0.5","") << std::endl;
      ss << Util::formatDescription("   epsilonC=0.2916","") << std::endl;
      ss << Util::formatDescription("   lambda=0.5","") << std::endl;
      ss << Util::formatDescription("   boxCoxThreshold=undef","") << std::endl;
      ss << Util::formatDescription("   rhoType=gaussian","One of 'gaussian', 'soar'") << std::endl;
      ss << Util::formatDescription("   diagnose=0","") << std::endl;
      ss << Util::formatDescription("   maxElevDiff=200","Remove stations that are further away from the background elevation than this (in meters)") << std::endl;
      ss << Util::formatDescription("   landOnly=0","Remove stations that are not on land (laf > 0)") << std::endl;
      ss << Util::formatDescription("   diaFile=undef","If defined, write information about removed stations to this filename") << std::endl;
      ss << Util::formatDescription("   crossValidate=0","If 1, then don't use the nearest point in the kriging. The end result is a field that can be verified against observations at the kriging points.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c oi","Merging of observations and background using Optimal interpolation")<< std::endl;
   return ss.str();
}
