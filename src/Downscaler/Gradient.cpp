#include "Gradient.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerGradient::DownscalerGradient(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions),
      mElevRadius(3),
      mMinElevGradient(Util::MV),
      mMaxElevGradient(Util::MV),
      mLogTransform(false),
      mConstantElevGradient(Util::MV),
      mDefaultElevGradient(0),
      mMinElevDiff(30),
      mAverageNeighbourhood(false),
      mSaveGradient(""),
      mMinLafForElevGradient(0),
      mMinLafDiff(0.1),
      mLafRadius(0),
      mDownscalerName("nearestNeighbour"),
      mMinNumPoints(2),
      mLimit(false),
      mMinFracSeaPoints(0),
      mElevGradientVariableName("") {

   iOptions.getValue("constantElevGradient", mConstantElevGradient);
   iOptions.getValue("elevRadius", mElevRadius);
   iOptions.getValue("minElevDiff", mMinElevDiff);
   iOptions.getValue("defaultElevGradient", mDefaultElevGradient);
   iOptions.getValue("minElevGradient", mMinElevGradient);
   iOptions.getValue("maxElevGradient", mMaxElevGradient);
   iOptions.getValue("minLafForElevGradient", mMinLafForElevGradient);
   iOptions.getValue("minNumPoints", mMinNumPoints);
   if(!iOptions.getValue("elevGradientVariable", mElevGradientVariableName)) {
      mElevGradientVariableName = mInputVariable.name();
   }
   iOptions.getValue("lafRadius", mLafRadius);
   iOptions.getValue("minLafDiff", mMinLafDiff);
   iOptions.getValue("minFracSeaPoints", mMinFracSeaPoints);
   iOptions.getValue("minLafGradient", mMinLafGradient);
   iOptions.getValue("maxLafGradient", mMaxLafGradient);
   iOptions.getValues("lafSearchRadii", mLafSearchRadii);
   iOptions.getValues("lafWeights", mLafWeights);

   iOptions.getValue("logTransform", mLogTransform);
   iOptions.getValue("averageNeighbourhood", mAverageNeighbourhood);
   iOptions.getValue("saveGradient", mSaveGradient);
   iOptions.getValue("downscaler", mDownscalerName);
   iOptions.getValue("limit", mLimit);
   // mDownscaler = Downscaler::getScheme(downscalerName, iVariable, Options());

   iOptions.getValues("lafSearchRadii", mLafSearchRadii);

   if(!iOptions.getValues("lafWeights", mLafWeights)) {
      for(int i = 0; i < mLafSearchRadii.size(); i++) {
         mLafWeights.push_back(1.0/mLafSearchRadii.size());
      }
   }
   assert(mLafSearchRadii.size() == mLafWeights.size());
   iOptions.check();
}

DownscalerGradient::~DownscalerGradient() {
   // delete mDownscaler;
}

void DownscalerGradient::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumY();
   int nLon = iOutput.getNumX();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   vec2 ilats  = iInput.getLats();
   vec2 ilons  = iInput.getLons();
   vec2 ielevs = iInput.getElevs();
   vec2 ilafs = iInput.getLandFractions();
   vec2 olats  = iOutput.getLats();
   vec2 olons  = iOutput.getLons();
   vec2 oelevs = iOutput.getElevs();
   vec2 olafs = iOutput.getLandFractions();

   float minAllowed = mOutputVariable.min();
   float maxAllowed = mOutputVariable.max();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);
      Field& gfield = *iInput.getField(mElevGradientVariableName, t);

      vec2 elevsInterp;
      vec2 lafsInterp;
      if(mDownscalerName == "nearestNeighbour") {
         elevsInterp = DownscalerNearestNeighbour::downscaleVec(ielevs, ilats, ilons, olats, olons, nearestI, nearestJ);
         lafsInterp = DownscalerNearestNeighbour::downscaleVec(ilafs, ilats, ilons, olats, olons, nearestI, nearestJ);
         DownscalerNearestNeighbour::downscaleField(ifield, ofield, ilats, ilons, olats, olons, nearestI, nearestJ);
      }
      else if (mDownscalerName == "bilinear") {
         elevsInterp = DownscalerBilinear::downscaleVec(ielevs, ilats, ilons, olats, olons, nearestI, nearestJ);
         lafsInterp = DownscalerBilinear::downscaleVec(ilafs, ilats, ilons, olats, olons, nearestI, nearestJ);
         DownscalerBilinear::downscaleField(ifield, ofield, ilats, ilons, olats, olons, nearestI, nearestJ);
      }

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int Icenter = nearestI[i][j];
            int Jcenter = nearestJ[i][j];
            assert(Icenter < ielevs.size());
            assert(Jcenter < ielevs[Icenter].size());
            for(int e = 0; e < nEns; e++) {
               float currElev = oelevs[i][j];
               float currLaf = Util::MV;
               if(mLafRadius == 0) {
                  currLaf = olafs[i][j];
               }
               else {
                  // Find the average LAF in a neighbourhood
                  int count = 0;
                  float total = 0;
                  for(int ii = std::max(0, i-mLafRadius); ii <= std::min(nLat-1, i+mLafRadius); ii++) {
                     for(int jj = std::max(0, j-mLafRadius); jj <= std::min(nLon-1, j+mLafRadius); jj++) {
                        if(Util::isValid(olafs[ii][jj])) {
                           total += olafs[ii][jj];
                           count++;
                        }
                     }
                  }
                  if(count > 0)
                     currLaf = total / count;
               }

               // TODO: What base value do we use? For elevation gradient, we want to use the nearest
               // neighbour, but for LAF gradient, we want to use the lowest LAF point?
               // float nearestElev = ielevs[Icenter][Jcenter];
               // float nearestLaf = ilafs[Icenter][Jcenter];
               // float nearestValue = ifield(Icenter, Jcenter, e);
               float baseElev = elevsInterp[i][j];
               float baseLaf = lafsInterp[i][j];
               float baseValue = ofield(i, j, e);

               // Missing ensemble member
               if(!Util::isValid(baseValue)) {
                  continue;
               }

               // Alternatively, use a neighbourhood mean
               // int averagingRadius = mAverageNeighbourhood * mElevRadius;
               // calcNeighbourhoodMean(ifield, ielevs, i, j, e, Icenter, Jcenter, averagingRadius, baseValue, baseElev, baseLaf);

               float elevGradient = calcElevGradient(i, j, e, Icenter, Jcenter, ifield, gfield, ielevs, ilafs);
               float lafGradient = Util::MV;
               if(mLafSearchRadii.size() > 0)
                  lafGradient = calcLafGradient(i, j, e, Icenter, Jcenter, ifield, ilafs, ielevs, elevGradient);

               if(mSaveGradient == "elev") {
                  ofield(i, j, e) = elevGradient;
               }
               else if(mSaveGradient == "laf") {
                  ofield(i, j, e) = lafGradient;
               }
               else {
                  /* General approach
                     1) Compute the elevation gradient, by ignoring any ocean points
                     2) Compute the LAF gradient, by adjusting neighbourhood using elevation gradient
                     3) Adjust the point with the lower LAF by adding on the LAF gradient and elevation gradient
                  */
                  float dElev = currElev - baseElev;
                  float dLaf = currLaf - baseLaf;

                  float value = baseValue + dElev * elevGradient;
                  if(Util::isValid(lafGradient))
                     value += dLaf * lafGradient;

                  if((Util::isValid(minAllowed) && value < minAllowed) || (Util::isValid(maxAllowed) && value > maxAllowed)) {
                     value = baseValue;
                  }

                  // Ensure that the value is inside the range in the neighbourhood
                  if(mLimit) {
                     float min = Util::MV;
                     float max = Util::MV;
                     for(int ii = std::max(0, Icenter-mElevRadius); ii <= std::min(ifield.getNumY()-1, Icenter+mElevRadius); ii++) {
                        for(int jj = std::max(0, Jcenter-mElevRadius); jj <= std::min(ifield.getNumX()-1, Jcenter+mElevRadius); jj++) {
                           float curr = ifield(ii, jj, e);
                           if(Util::isValid(curr)) {
                              if(!Util::isValid(min) || curr < min)
                                 min = curr;
                              else if(!Util::isValid(max) || curr > max)
                                 max = curr;
                           }
                        }
                     }
                     if(value > max)
                        value = max;
                     if(value < min)
                        value = min;
                  }
                  ofield(i,j,e)  = value;
               }
            }
         }
      }
   }
}

bool DownscalerGradient::calcBaseValues(const Field& iField, const vec2& iElevs, const vec2& iLafs, int i, int j, int e, int Icenter, int Jcenter, int iRadius, float& iValueMean, float& iElevMean, float& iLafMean) const {
   float totalValue = 0;
   float totalElev = 0;
   float totalLaf = 0;
   int counter = 0;
   int counterLaf = 0;
   for(int ii = std::max(0, Icenter-iRadius); ii <= std::min(iField.getNumY()-1, Icenter+iRadius); ii++) {
      for(int jj = std::max(0, Jcenter-iRadius); jj <= std::min(iField.getNumX()-1, Jcenter+iRadius); jj++) {
         float currValue = iField(ii,jj,e);
         float currElev  = iElevs[ii][jj];
         float currLaf  = iLafs[ii][jj];
         if(Util::isValid(currValue) && Util::isValid(currElev)) {
            totalValue += currValue;
            totalElev  += currElev;
            counter++;
         }
         if(Util::isValid(currLaf)) {
            totalLaf += currLaf;
            counterLaf++;
         }
      }
   }
   if(counter == 0) {
      iValueMean = Util::MV;
      iElevMean = Util::MV;
   }
   else {
      iValueMean = totalValue / counter;
      iElevMean = totalElev / counter;
   }
   if(counterLaf == 0) {
      iLafMean = Util::MV;
   }
   else {
      iLafMean = totalLaf / counter;
   }

   return counter > 0;
}
bool DownscalerGradient::calcNeighbourhoodMean(const Field& iField, const vec2& iElevs, int i, int j, int e, int Icenter, int Jcenter, int iRadius, float& iValueMean, float& iElevMean) const {
   float totalValue = 0;
   float totalElev = 0;
   int counter = 0;
   for(int ii = std::max(0, Icenter-iRadius); ii <= std::min(iField.getNumY()-1, Icenter+iRadius); ii++) {
      for(int jj = std::max(0, Jcenter-iRadius); jj <= std::min(iField.getNumX()-1, Jcenter+iRadius); jj++) {
         float currValue = iField(ii,jj,e);
         float currElev  = iElevs[ii][jj];
         if(Util::isValid(currValue) && Util::isValid(currElev)) {
            totalValue += currValue;
            totalElev  += currElev;
            counter++;
         }
      }
   }
   if(counter == 0) {
      iValueMean = Util::MV;
      iElevMean = Util::MV;
   }
   else {
      iValueMean = totalValue / counter;
      iElevMean = totalElev / counter;
   }

   return counter > 0;
}
float DownscalerGradient::calcLafGradient(int i, int j, int e, int Icenter, int Jcenter, const Field& iField, const vec2& iLafs, const vec2& iElevs, float iElevGradient) const {
   /* Here is the formula:
      value = minT + (maxT  - minT) / (maxLaf - minLaf) * (laf - minLaf)
      value = minT + gradient * dLaf
   */
   // Compute the average temperature and LAF in a small radius
   float minLaf = Util::MV;
   float maxLaf = Util::MV;
   float minElev = Util::MV;
   float maxElev = Util::MV;
   float minT = Util::MV;
   float maxT = Util::MV;
   // std::vector<float> values(mSearchRadii.size(), Util::MV);
   float value = 0;
   int counter = 0;
   float totalWeight = 0;
   for(int r = 0; r < mLafSearchRadii.size(); r++) {
      int numSeaPoints = 0;
      int numPoints = 0;
      int searchRadius = mLafSearchRadii[r];
      for(int ii = std::max(0, Icenter-searchRadius); ii <= std::min(iField.getNumY()-1, Icenter+searchRadius); ii++) {
         for(int jj = std::max(0, Jcenter-searchRadius); jj <= std::min(iField.getNumX()-1, Jcenter+searchRadius); jj++) {
            float currLaf  = iLafs[ii][jj];
            float currElev = iElevs[ii][jj];
            float currValue = iField(ii,jj,e);
            if(Util::isValid(currValue)) {
               numSeaPoints += currLaf < 0.5;
               numPoints++;
            }
            if(Util::isValid(currValue) && Util::isValid(currLaf)) {
               if(!Util::isValid(minLaf) || currLaf < minLaf) {
                  minLaf = currLaf;
                  minElev = currElev;
                  minT = currValue;
               }
               if(!Util::isValid(maxLaf) || currLaf > maxLaf) {
                  maxLaf = currLaf;
                  maxElev = currElev;
                  maxT = currValue;
               }
            }
         }
      }
      // std::cout << minLaf << " " << maxLaf << " " << minT << " " << maxT << std::endl;
      if(Util::isValid(maxLaf) && Util::isValid(minLaf) && (maxLaf - minLaf) >= mMinLafDiff && (float) numSeaPoints / numPoints >= mMinFracSeaPoints) {
         float lafGradient = (maxT - minT) / (maxLaf - minLaf);
         // Altitude-adjust maxT so that it is on the same elevation as minT
         // Before computing the gradient
         if(Util::isValid(maxElev) && Util::isValid(minElev)) {
            float maxTcorr = maxT - iElevGradient * (maxElev - minElev);
            lafGradient = (maxTcorr - minT) / (maxLaf - minLaf);
            value += lafGradient * mLafWeights[r];
            totalWeight += mLafWeights[r];
            counter++;
         }
         /*
         // Safety check
         if(!Util::isValid(lafGradient))
            lafGradient = 0;
         // Check against minimum and maximum gradients
         if(Util::isValid(mMinElevGradient) && lafGradient < mMinElevGradient)
            lafGradient = mMinElevGradient;
         if(Util::isValid(mMaxElevGradient) && lafGradient > mMaxElevGradient)
            lafGradient = mMaxElevGradient;
         // std::cout << lafGradient << std::endl;
         // float currValue = nearestValue + dLaf * lafGradient;
         float currValue = minT + dLaf * lafGradient;

         // Altitude-adjust the current point
         if(Util::isValid(currElev) && Util::isValid(minElev)) {
            float dElev = currElev - minElev;
            currValue = currValue + dElev * mElevGradient;
         }
         // std::cout << currValue << std::endl;
         value += currValue * mWeights[r];
         totalWeight += mWeights[r];
         counter++;
         if(j == 43 && i == 56 && t == 0) {
            float dElev = currElev - minElev;
            float elevCorrection = dElev * mElevGradient;
            std::cout << minLaf << " " << maxLaf << " " << currLaf << "\n"
                      << minT << " " << maxT << " " << nearestValue << "\n"
                      << minElev << " " << maxElev << " " << currElev << "\n"
                      << currValue << " " << nearestValue << " " << value << "\n"
                      << lafGradient << " " << elevCorrection << std::endl;
         }
         */
      }
   }
   // std::cout << value << std::endl;
   float gradient;
   if(counter == 0 || !Util::isValid(value)) {
      gradient = Util::MV;
   }
   else {
      gradient = value / totalWeight;
   }
   return gradient;
}
float DownscalerGradient::calcElevGradient(int i, int j, int e, int Icenter, int Jcenter, const Field& iField, const Field& iGfield, const vec2& iElevs, const vec2& iLafs) const {
   float elevGradient = Util::MV;
   if(Util::isValid(mConstantElevGradient)) {
      elevGradient = mConstantElevGradient;
   }
   else {
      /* Compute the model's gradient:
         The gradient is computed by using linear regression on forecast ~ elevation
         using all forecasts within a neighbourhood. To produce stable results, there
         is a requirement that the elevation within the neighbourhood has a large
         range (see mMinElevDiff).

         For bounded variables (e.g. wind speed), the gradient approach could cause
         forecasts to go outside its domain (e.g. negative winds). If this occurs,
         the nearest neighbour is used.
      */
      float meanXY  = 0; // elev*T
      float meanX   = 0; // elev
      float meanY   = 0; // T
      float meanXX  = 0; // elev*elev
      int   counter = 0;
      float min = Util::MV;
      float max = Util::MV;
      for(int ii = std::max(0, Icenter-mElevRadius); ii <= std::min(iField.getNumY()-1, Icenter+mElevRadius); ii++) {
         for(int jj = std::max(0, Jcenter-mElevRadius); jj <= std::min(iField.getNumX()-1, Jcenter+mElevRadius); jj++) {
            assert(ii < iElevs.size());
            assert(jj < iElevs[ii].size());
            float x = iElevs[ii][jj];
            float y = iGfield(ii,jj,e);
            float currLaf = iLafs[ii][jj];
            if(mLogTransform) {
               y = log(y);
            }
            if(Util::isValid(x) && Util::isValid(y) &&(!Util::isValid(currLaf) || currLaf >= mMinLafForElevGradient)) {
               meanXY += x*y;
               meanX  += x;
               meanY  += y;
               meanXX += x*x;
               counter++;
               // Found a new min
               if(!Util::isValid(min) || x < min)
                  min = x;
               // Found a new max
               if(!Util::isValid(max) || x > max)
                  max = x;
            }
         }
      }
      // Compute elevation difference within neighbourhood
      float elevDiff = Util::MV;
      if(Util::isValid(min) && Util::isValid(max)) {
         assert(max >= min);
         elevDiff = max - min;
      }

      // Use model gradient if:
      // 1) sufficient elevation difference in neighbourhood
      // 2) regression parameters are stable enough
      if(counter >= mMinNumPoints && Util::isValid(elevDiff) && elevDiff >= mMinElevDiff && meanXX != meanX*meanX) {
         // Estimate lapse rate
         meanXY /= counter;
         meanX  /= counter;
         meanY  /= counter;
         meanXX /= counter;
         if(meanXX - meanX*meanX != 0) {
            elevGradient = (meanXY - meanX*meanY)/(meanXX - meanX*meanX);
         }
      }
   }

   // Safety check
   if(!Util::isValid(elevGradient))
      elevGradient = mDefaultElevGradient;
   // Check against minimum and maximum gradients
   if(Util::isValid(mMinElevGradient) && elevGradient < mMinElevGradient)
      elevGradient = mMinElevGradient;
   if(Util::isValid(mMaxElevGradient) && elevGradient > mMaxElevGradient)
      elevGradient = mMaxElevGradient;
   return elevGradient;
}
std::string DownscalerGradient::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-d gradient", "Adjusts the nearest neighbour based on elevation or coastal gradients in a neighbourhood. The gradient is applied using the difference of the elevation and /or land-area fraction of the output gridpoint to the nearest neighbour: T = T(nn) + elev_gradient * dElev + laf_gradient * dLaf. If the gradient puts the forecast outside the domain of the variable (e.g. negative precipitation) then the nearest neighbour is used.") << std::endl;
      ss << Util::formatDescription("   constantElevGradient=undef", "Fix elevation gradient to this value. Units are units_of_variable / meter. Positive values mean an increase with height. If unspecified, computes the gradient by regression of points in a neighbourhood. The remaining options are ignored.") << std::endl;
      ss << Util::formatDescription("   elevRadius=3", "Compute elevation gradient in a neighbourhood box of points within within +- radius in both east-west and north-south direction. Also use this neighbourhood to compute the vertical gradient.") << std::endl;
      ss << Util::formatDescription("   minElevDiff=30", "Minimum elevation range (in meters) required within the neighbourhood to compute gradient. Use nearest neighbour otherwise.") << std::endl;
      ss << Util::formatDescription("   defaultElevGradient=0", "If the elevation gradient is not computable (too unstable), use this default gradient.") << std::endl;
      ss << Util::formatDescription("   minElevGradient=undef", "Do not allow elevation gradient to be smaller than this value. If undefined, do not alter gradient.") << std::endl;
      ss << Util::formatDescription("   maxElevGradient=undef", "Do not allow elevation gradient to be larger than this value. If undefined, do not alter gradient.") << std::endl;
      ss << Util::formatDescription("   minLafForElevGradient=0", "When computing the elevation gradient, only include points with a LAF above this value.") << std::endl;
      ss << Util::formatDescription("   minNumPoints=2", "Minimum number of points needed to compute elevation gradient.") << std::endl;
      ss << Util::formatDescription("   elevGradientVariable=undef", "Which variable to use for the gradient? If undefined, use the same variables as the forecast variable.") << std::endl;

      ss << Util::formatDescription("   lafRadius=1", "Smoothen the LAF output with this neighbourhood radius.") << std::endl;
      ss << Util::formatDescription("   minLafDiff=0.1", "Minimum LAF range (between 0 and 1) required within the neighbourhood to compute gradient. Use nearest neighbour otherwise.") << std::endl;
      ss << Util::formatDescription("   minFracSeaPoints=0", "Minimum fraction of points in neighbourhood with a LAF < 0.5 to compute a LAF gradient.") << std::endl;
      ss << Util::formatDescription("   minLafGradient=undef", "Do not allow LAF gradient to be smaller than this value. If undefined, do not alter gradient.") << std::endl;
      ss << Util::formatDescription("   maxLafGradient=undef", "Do not allow LAF gradient to be larger than this value. If undefined, do not alter gradient.") << std::endl;
      ss << Util::formatDescription("   lafSearchRadii=1,2,3", "Compute LAF gradients using neighbourhoods with these sizes.") << std::endl;
      ss << Util::formatDescription("   lafWeights=undef", "Weights when computing the average LAF gradient in different neighbourhoods.") << std::endl;

      ss << Util::formatDescription("   logTransform=0", "Should the variable be log-transformed first? I.e should a linear gradient be applied to log(variable)? T = T(nn) * exp(gradient * dElev) * exp(gradient * dLaf). Useful for pressure variables. Can be used together with constantElevGradient.") << std::endl;
      ss << Util::formatDescription("   averageNeighbourhood=0", "Should the average forecast, elevation, and LAF within the search radius be used when determining what value to apply the gradient to?") << std::endl;
      ss << Util::formatDescription("   saveGradient=""", "Store the gradient instead of the value for debugging purposes. Use elev to store the elevation gradient; Use laf to store the laf gradient.") << std::endl;
      ss << Util::formatDescription("   downscaler=nearestNeighbour", "Use this downscaler on the field and on the altitudes before applying gradients.") << std::endl;
      ss << Util::formatDescription("   limit=0", "Should the downscaled value be forced to be within the min/max range of the elevation neighbourhood?") << std::endl;
   }
   else
      ss << Util::formatDescription("-d gradient", "Adjust nearest neighbour based on elevation and coastal gradients") << std::endl;
   return ss.str();
}
