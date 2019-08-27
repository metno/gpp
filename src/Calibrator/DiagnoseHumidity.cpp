#include "DiagnoseHumidity.h"
#include "../Downscaler/Pressure.h"
#include "../Util.h"
#include "../File/File.h"
#include <cmath>
CalibratorDiagnoseHumidity::CalibratorDiagnoseHumidity(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mTemperature(""),
      mRh(""),
      mDewpoint(""),
      mPressure(""),
      mCompute("") {
   iOptions.getRequiredValue("temperature", mTemperature);
   iOptions.getValue("rh", mRh);
   iOptions.getValue("dewpoint", mDewpoint);
   iOptions.getValue("pressure", mPressure);
   iOptions.getRequiredValue("compute", mCompute);
   if(mCompute == "dewpoint" && mRh == "")
      Util::error("Dewpoint requires rh=");
   if(mCompute == "rh" && mDewpoint == "")
      Util::error("Rh requires dewpoint=");
   if(mCompute == "wetbulb" && mRh == "")
      Util::error("Wetbulb requires rh=");
}
bool CalibratorDiagnoseHumidity::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 elevs = iFile.getElevs();

   // Get all fields
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      if(mCompute == "dewpoint") {
         if(!iFile.hasVariable(mTemperature))
            Util::error("Cannot diagnose dewpoint, since temperature field is missing");
         if(!iFile.hasVariable(mRh))
            Util::error("Cannot diagnose dewpoint, since rh field is missing");
         const Field& temperature = *iFile.getField(mTemperature, t);
         const Field& rh = *iFile.getField(mRh, t);

         #pragma omp parallel for
         for(int y = 0; y < nY; y++) {
            for(int x = 0; x < nX; x++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(temperature(y, x, e)) && Util::isValid(rh(y, x, e)))
                     output(y, x, e) = computeDewpoint(temperature(y, x, e), rh(y, x, e));
               }
            }
         }
      }

      else if(mCompute == "rh") {
         if(!iFile.hasVariable(mTemperature))
            Util::error("Cannot diagnose rh, since temperature field is missing");
         if(!iFile.hasVariable(mDewpoint))
            Util::error("Cannot diagnose rh, since dewpoint field is missing");
         const Field& temperature = *iFile.getField(mTemperature, t);
         const Field& dewpoint = *iFile.getField(mDewpoint, t);

         #pragma omp parallel for
         for(int y = 0; y < nY; y++) {
            for(int x = 0; x < nX; x++) {
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(temperature(y, x, e)) && Util::isValid(dewpoint(y, x, e)))
                     output(y, x, e) = computeRh(temperature(y, x, e), dewpoint(y, x, e));
               }
            }
         }
      }

      else if(mCompute == "wetbulb") {
         if(!iFile.hasVariable(mTemperature))
            Util::error("Cannot diagnose wetbulb, since temperature field is missing");
         if(!iFile.hasVariable(mRh))
            Util::error("Cannot diagnose wetbulb, since rh field is missing");
         const Field& temperature = *iFile.getField(mTemperature, t);
         const Field& rh = *iFile.getField(mRh, t);
         FieldPtr pressure;
         if(mPressure != "") {
            if(!iFile.hasVariable(mPressure))
               Util::error("Cannot diagnose wetbulb, since pressure field is missing");
            pressure = iFile.getField(mPressure, t);
         }
         else {
            if(!iFile.hasElevs()) {
               Util::error("Cannot diagnose wetbulb, since pressure can be diagnosed since altitude is missing");
            }
         }

         #pragma omp parallel for
         for(int y = 0; y < nY; y++) {
            for(int x = 0; x < nX; x++) {
               float currElev = elevs[y][x];
               for(int e = 0; e < nEns; e++) {
                  if(Util::isValid(temperature(y, x, e)) && Util::isValid(rh(y, x, e))) {
                     float currPressure = 0;
                     if(mPressure != "")
                        currPressure = (*pressure)(y, x, e);
                     else {
                        currPressure = DownscalerPressure::calcPressure(0, 101325, currElev);
                     }
                     output(y, x, e) = computeWetbulb(temperature(y, x, e), currPressure, rh(y, x, e));
                  }
               }
            }
         }
      }
   }
   return true;
}

float CalibratorDiagnoseHumidity::computeDewpoint(float iTemperature, float iRelativeHumidity) {
   if(Util::isValid(iTemperature) && Util::isValid(iRelativeHumidity)) {
      // Taken from https://github.com/metno/wdb2ts
      float tempC = iTemperature - 273.15;
      float e = (iRelativeHumidity)*0.611*exp( (17.63 * tempC) / (tempC + 243.04) );
      float tdC = (116.9 + 243.04 * log( e ))/(16.78-log( e ));
      float td = tdC + 273.15;
      return (td<=iTemperature ? td : iTemperature);
      // Taken from https://github.com/WFRT/comps
      // float es = 0.611*exp(5423.0*(1/273.15 - 1/(iTemperature)));
      // float e  = iRelativeHumidity*(es);
      // float td = 1/(1/273.15 - 1.844e-4*log(e/0.611));
      // return td;
   }
   else
      return Util::MV;
}

float CalibratorDiagnoseHumidity::mEwt[41] = { .000034,.000089,.000220,.000517,.001155,.002472,
                            .005080,.01005, .01921, .03553, .06356, .1111,
                            .1891,  .3139,  .5088,  .8070,  1.2540, 1.9118,
                            2.8627, 4.2148, 6.1078, 8.7192, 12.272, 17.044,
                            23.373, 31.671, 42.430, 56.236, 73.777, 95.855,
                            123.40, 157.46, 199.26, 250.16, 311.69, 385.56,
                            473.67, 578.09, 701.13, 845.28, 1013.25 };

float CalibratorDiagnoseHumidity::computeRh(float iTemperature, float iDewPointTemperature) {
   // Taken from https://github.com/metno/wdb2ts
   if(Util::isValid(iTemperature) && Util::isValid(iDewPointTemperature)) {
      float x, et, etd, rh;
      int l;

      x = (iTemperature - 173.16) * 0.2;
      l = int(x);
      et = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));

      x = (iDewPointTemperature - 173.16) * 0.2;
      l = int(x);
      etd = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));
      rh = etd / et;

      if(rh < 0)
         rh = 0;
      if(rh > 1)
         rh = 1;

      return rh;
   }
   else
      return Util::MV;
}

float CalibratorDiagnoseHumidity::computeWetbulb(float iTemperature, float iPressure, float iRelativeHumidity) {
   float temperatureC = iTemperature - 273.15;
   if(temperatureC <= -243.04 || iRelativeHumidity <= 0)
      return Util::MV;
   if(Util::isValid(temperatureC) && Util::isValid(iPressure) && Util::isValid(iRelativeHumidity)) {
      float e  = (iRelativeHumidity)*0.611*exp((17.63*temperatureC)/(temperatureC+243.04));
      float Td = (116.9 + 243.04*log(e))/(16.78-log(e));
      float gamma = 0.00066 * iPressure/1000;
      float delta = (4098*e)/pow(Td+243.04,2);
      if(gamma + delta == 0)
         return Util::MV;
      float wetbulbTemperature   = (gamma * temperatureC + delta * Td)/(gamma + delta);
      float wetbulbTemperatureK  = wetbulbTemperature + 273.15;
      return wetbulbTemperatureK;
   }
   else {
      return Util::MV;
   }
}

std::string CalibratorDiagnoseHumidity::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c diagnoseHumidity","Compute dewpoint (from temperature, rh), wetbulb (from temperature, rh, pressure), or relative humidity (from temperature, dewpoint).") << std::endl;
      ss << Util::formatDescription("   temperature=required","Temperature variable name") << std::endl;
      ss << Util::formatDescription("   rh=undef","Relative humidity variable name") << std::endl;
      ss << Util::formatDescription("   dewpoint=undef","Dewpoint temperature variable name") << std::endl;
      ss << Util::formatDescription("   pressure=undef","Pressure variable name") << std::endl;
      ss << Util::formatDescription("   compute=required","Which variable type should be diagnosed? One of dewpoint, wetbulb, rh.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c diagnoseHumidity","Diagnose dewpoint, wetbulb, or relative humidity") << std::endl;
   return ss.str();
}

