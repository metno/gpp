#include "Regression.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorRegression::CalibratorRegression(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mOrder(1),
      mIntercept(true) {
   iOptions.getValue("order", mOrder);
   iOptions.getValue("intercept", mIntercept);
   iOptions.getValues("variables", mVariables);
   iOptions.check();
}
bool CalibratorRegression::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   if(iParameterFile->getNumParameters() == 0) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have at least one dataacolumns");
   }
   if(mVariables.size() > 0 && iParameterFile->getNumParameters() != mVariables.size()) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have at the same number of parameters as number of variables in regression");
   }

   bool multiVariate = mVariables.size() > 0;

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parametersGlobal;
      if(!iParameterFile->isLocationDependent())
         parametersGlobal = iParameterFile->getParameters(t);
      const FieldPtr field = iFile.getField(mVariable, t);
      std::vector<FieldPtr> fields;
      if(multiVariate) {
         for(int i = 0; i < mVariables.size(); i++) {
            if(mVariables[i] != "1") {
               fields.push_back(iFile.getField(mVariables[i], t));
            }
            else {
               fields.push_back(iFile.getEmptyField(1));
            }
         }
      }
      else {
         fields.push_back(field);
      }

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {

            Parameters parameters;
            if(iParameterFile->isLocationDependent())
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
            else
               parameters = parametersGlobal;

            for(int e = 0; e < nEns; e++) {
               if(Util::isValid((*field)(i,j,e))) {
                  if(multiVariate) {
                     float total = 0;
                     // Accumulate a + b * var1 + c * var2 ...
                     for(int p = 0; p < parameters.size(); p++) {
                        float coeff = parameters[p];
                        if(!Util::isValid(coeff)) {
                           total = Util::MV;
                           break;
                        }
                        total += coeff*(*fields[p])(i,j,e);
                     }
                     (*field)(i,j,e)  = total;
                  }
                  else {
                     float total = 0;
                     // Accumulate a + b * fcst + c * fcst^2 ...
                     for(int p = 0; p < parameters.size(); p++) {
                        float coeff = parameters[p];
                        if(!Util::isValid(coeff)) {
                           total = Util::MV;
                           break;
                        }
                        total += coeff*pow((*fields[0])(i,j,e), p);
                     }
                     (*field)(i,j,e)  = total;
                  }
               }
               else {
                  (*field)(i,j,e)  = Util::MV;
               }
            }
         }
      }
   }
   return true;
}

Parameters CalibratorRegression::train(const std::vector<ObsEns>& iData) const {
   if(mOrder > 1) {
      std::stringstream ss;
      ss << "CalibratorRegression: Cannot train regression of order greater than 1 (i.e a + bx)";
      Util::error(ss.str());
   }
   if(mVariables.size() > 0) {
      std::stringstream ss;
      ss << "CalibratorRegression: Cannot train multivariate regression models";
      Util::error(ss.str());
   }

   float totalObs = 0;
   float totalForecast = 0;
   float totalForecast2 = 0;
   float totalObsForecast = 0;
   int counter = 0;
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      float obs = iData[i].first;
      std::vector<float> ens = iData[i].second;
      float ensMean = Util::calculateStat(ens, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(ensMean)) {
         totalObs += obs;
         totalForecast += ensMean;
         totalForecast2 += ensMean*ensMean;
         totalObsForecast += obs*ensMean;
         counter++;
      }
   }

   if(counter <= 0) {
      std::stringstream ss;
      ss << "CalibratorRegression: Cannot train regression, no valid data";
      Util::error(ss.str());
   }

   std::vector<float> values;
   float meanObs = totalObs / counter;
   float meanForecast = totalForecast / counter;
   float meanForecast2 = totalForecast2 / counter;
   float meanObsForecast = totalObsForecast / counter;
   if(mOrder == 0) {
      values.push_back(meanObs);
   }
   else if(mOrder == 1) {
      float intercept;
      float slope;
      if(mIntercept) {
         slope = (meanObsForecast - meanForecast*meanObs)/(meanForecast2 - meanForecast*meanForecast);
         intercept = meanObs - slope*meanForecast;
      }
      else {
         intercept = 0;
         slope = meanObsForecast / meanForecast2;
      }
      values.push_back(intercept);
      values.push_back(slope);
   }
   else {
      abort();
   }

   Parameters par(values);

   return par;
}

std::string CalibratorRegression::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c regression", "Polynomial or multivariate regression equation. Set 'order' for polynomial regression and 'variables' for multivariate regression. Applies the following models: newForecast = a + b * forecast + c * forecast^2 ... or newForecast = a * var1 + b * var2.... A parameter file is required with the values [a b c ... ]") << std::endl;
      ss << Util::formatDescription("   order=1", "What order is the polynomial regression? 0th order: a; 1st order: a + bx; 2nd order a + bc + cx^2; etc.") << std::endl;
      ss << Util::formatDescription("   variables=undef", "Use these variables as predictors. Use '1' to represent a constant term.") << std::endl;
      ss << Util::formatDescription("   intercept=1", "Applies only to training of coefficients for polynomial regression. Should the regression include an intercept term? If 1, then yes.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c regression", "Polynomial or multivariate regression") << std::endl;
   return ss.str();
}
