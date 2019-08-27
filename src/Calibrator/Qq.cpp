#include "Qq.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorQq::CalibratorQq(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mLowerQuantile(0),
      mUpperQuantile(1),
      mPolicy(ExtrapolationPolicy::OneToOne) {

   std::string extrapolation;
   if(iOptions.getValue("extrapolation", extrapolation)) {
      if(extrapolation == "1to1")
         mPolicy = ExtrapolationPolicy::OneToOne;
      else if(extrapolation == "meanSlope")
         mPolicy = ExtrapolationPolicy::MeanSlope;
      else if(extrapolation == "nearestSlope")
         mPolicy = ExtrapolationPolicy::NearestSlope;
      else if(extrapolation == "zero")
         mPolicy = ExtrapolationPolicy::Zero;
      else {
         Util::error("CalibratorQq: value for 'extrapolation' not recognized");
      }
   }
   iOptions.getValues("quantiles", mQuantiles);
   iOptions.getValues("extraObs", mExtraObs);
   iOptions.getValues("extraFcst", mExtraFcst);
   if(mExtraObs.size() != mExtraFcst.size()) {
      Util::error("extraObs must have the same length as extraFcst");
   }
   iOptions.check();
}
bool CalibratorQq::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(iParameterFile->getNumParameters() % 2 != 0) {
      std::stringstream ss;
      ss << "Parameter file '" + iParameterFile->getFilename() + "' must have an even number of datacolumns (currently " << iParameterFile->getNumParameters() << " parameters)";
      Util::error(ss.str());
   }

   const int nLat = iFile.getNumY();
   const int nLon = iFile.getNumX();
   const int nEns = iFile.getNumEns();
   const int nTime = iFile.getNumTime();
   const vec2 lats = iFile.getLats();
   const vec2 lons = iFile.getLons();
   const vec2 elevs = iFile.getElevs();

   for(int t = 0; t < nTime; t++) {
      const FieldPtr field = iFile.getField(mVariable, t);

      // Retrieve the calibration parameters for this time
      // Overwrite them later if they are location dependent
      std::vector<float> obsVecGlobal, fcstVecGlobal;
      if(!iParameterFile->isLocationDependent()) {
         Parameters parameters = iParameterFile->getParameters(t);
         separate(parameters, obsVecGlobal, fcstVecGlobal);
      }
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            std::vector<float> obsVec, fcstVec;
            if(iParameterFile->isLocationDependent()) {
               Parameters parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
               separate(parameters, obsVec, fcstVec);
            }
            else {
               obsVec = obsVecGlobal;
               fcstVec = fcstVecGlobal;
            }
            int N = obsVec.size();
            if(obsVec.size() < 1) {
               Util::error("CalibratorQq cannot use parameters with size less than 2");
            }
            // Only process if all parameters are valid
            bool hasValidParameters = true;
            for(int p = 0; p < obsVec.size(); p++) {
               hasValidParameters = hasValidParameters && Util::isValid(obsVec[p]);
            }
            for(int p = 0; p < fcstVec.size(); p++) {
               hasValidParameters = hasValidParameters && Util::isValid(fcstVec[p]);
            }
            if(hasValidParameters) {
               for(int e = 0; e < nEns; e++) {
                  float raw = (*field)(i,j,e);
                  float value = Util::MV;
                  if(Util::isValid(raw)) {
                     float smallestObs  = obsVec[0];
                     float smallestFcst = fcstVec[0];
                     float largestObs   = obsVec[obsVec.size()-1];
                     float largestFcst  = fcstVec[fcstVec.size()-1];

                     // Linear interpolation within curve
                     if(raw > smallestFcst && raw < largestFcst) {
                        value = Util::interpolate(raw, fcstVec, obsVec);
                     }
                     // Extrapolate outside curve
                     else {
                        float nearestObs;
                        float nearestFcst;
                        if(raw <= smallestFcst) {
                           nearestObs  = smallestObs;
                           nearestFcst = smallestFcst;
                        }
                        else {
                           nearestObs  = largestObs;
                           nearestFcst = largestFcst;
                        }
                        float slope = 1;
                        if(mPolicy == ExtrapolationPolicy::Zero) {
                           slope = 0;
                        }
                        if(mPolicy == ExtrapolationPolicy::OneToOne || N <= 1) {
                           slope = 1;
                        }
                        else if(mPolicy == ExtrapolationPolicy::MeanSlope) {
                           float dObs  = largestObs - smallestObs;
                           float dFcst = largestFcst - smallestFcst;
                           slope = dObs / dFcst;
                        }
                        else if(mPolicy == ExtrapolationPolicy::NearestSlope) {
                           float dObs;
                           float dFcst;
                           if(raw <= smallestFcst) {
                              dObs  = obsVec[1] - obsVec[0];
                              dFcst = fcstVec[1] - fcstVec[0];
                           }
                           else {
                              dObs  = obsVec[N-1] - obsVec[N-2];
                              dFcst = fcstVec[N-1] - fcstVec[N-2];
                           }
                           slope = dObs / dFcst;
                        }
                        value = nearestObs + slope * (raw - nearestFcst);
                     }
                     (*field)(i,j,e) = value;
                  }
                  else {
                     (*field)(i,j,e)  = Util::MV;
                  }
               }
            }
         }
      }
   }
   return true;
}

Parameters CalibratorQq::train(const std::vector<ObsEns>& iData) const {
   double timeStart = Util::clock();
   std::vector<float> obs, fcst;
   std::vector<float> values;
   values.reserve(2*iData.size());
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      float currObs = iData[i].first;
      std::vector<float> ens = iData[i].second;
      float currFcst = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(currObs) && Util::isValid(currFcst)) {
         obs.push_back(currObs);
         fcst.push_back(currFcst);
      }
   }
   assert(obs.size() == fcst.size());

   if(obs.size() == 0 || fcst.size() == 0) {
      // std::stringstream ss;
      // ss << "CalibratorQq: No valid data, no correction will be made.";
      // Util::warning(ss.str());
   }

   // Sort
   std::sort(obs.begin(), obs.end());
   std::sort(fcst.begin(), fcst.end());
   std::vector<float> obsQ;
   std::vector<float> fcstQ;
   if(obs.size() == 0) {

   }
   else if(mQuantiles.size() > 0) {
      int N = obs.size();
      std::vector<float> x(N);
      for(int i = 0; i < N; i++) {
         x[i] = (float) i / (N-1);
      }
      for(int i = 0; i < mQuantiles.size(); i++) {
         float currObs  = Util::interpolate(mQuantiles[i], x, obs);
         float currFcst = Util::interpolate(mQuantiles[i], x, fcst);
         obsQ.push_back(currObs);
         fcstQ.push_back(currFcst);
      }
   }
   else {
      for(int i = 0; i < obs.size(); i++) {
         obsQ.push_back(obs[i]);
         fcstQ.push_back(fcst[i]);
      }
   }

   // Add extra points
   if(mExtraObs.size() > 0) {
      for(int i = 0; i < mExtraObs.size(); i++) {
         obsQ.push_back(mExtraObs[i]);
         fcstQ.push_back(mExtraFcst[i]);
      }
      std::sort(obsQ.begin(), obsQ.end());
      std::sort(fcstQ.begin(), fcstQ.end());
   }

   // Put into parameters
   for(int i = 0; i < obsQ.size(); i++) {
      values.push_back(obsQ[i]);
      values.push_back(fcstQ[i]);
   }

   Parameters par(values);

   double timeEnd = Util::clock();
   // std::cout << "Time: " << timeEnd - timeStart << std::endl;
   return par;
}

std::string CalibratorQq::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c qq", "Quantile-quantile mapping. Calibrates forecasts based on a map of sorted observations and sorted forecasts. For a given raw forecast, the quantile within the historical forecasts is found. Then the observation at the same quantile is used as the calibrated forecast. A parameter file is required with an even number of columns as follows") << std::endl;
      ss << Util::formatDescription("", "[obs0 fcs0 obs1 fcst1 .. obsN fcstN") << std::endl;
      ss << Util::formatDescription("", "Note that observations and forecasts must be sorted. I.e obs0 does not necessarily correspond to the time when fcst0 was issued. If a parameter set has one or more missing values, then the ensemble using this parameter set is not processed.") << std::endl;
      ss << Util::formatDescription("   extrapolation=1to1", "If a forecast is outside the curve, how should extrapolation be done? '1to1': Use a slope of 1, i.e. preserving the bias at the nearest end point; 'meanSlope': Use the average slope from the lower to upper endpoints; 'nearestSlope': Use the slope through the two nearest points; 'zero': Use a slope of 0, meaning that the forecast will equal the max/min observation.") << std::endl;
      ss << Util::formatDescription("   quantiles=undef", "If creating training parameters, should specific quantiles be stored, instead of all observations and forecasts? Can be a vector of values between 0 and 1.") << std::endl;
      ss << Util::formatDescription("   extraObs=undef", "Add these extra observations points to the curve. Must be the same length as extraFcst.") << std::endl;
      ss << Util::formatDescription("   extraFcst=undef", "Add these extra forecasts points to the curve at the end. Must be the same length as extraObs.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c qq", "Quantile-quantile mapping") << std::endl;
   return ss.str();
}

void CalibratorQq::separate(const Parameters& iParameters, std::vector<float>& iObs, std::vector<float>& iFcst) const {
   iObs.clear();
   iFcst.clear();
   int N = iParameters.size() / 2;
   iObs.resize(N, Util::MV);
   iFcst.resize(N, Util::MV);
   for(int i = 0; i < N; i++) {
      iObs[i] = iParameters[2*i];
      iFcst[i] = iParameters[2*i+1];
   }
   if(mExtraObs.size() > 0) {
      for(int i = 0; i < mExtraObs.size(); i++) {
         iObs.push_back(mExtraObs[i]);
         iFcst.push_back(mExtraFcst[i]);
      }
      std::sort(iObs.begin(), iObs.end());
      std::sort(iFcst.begin(), iFcst.end());
   }
}
