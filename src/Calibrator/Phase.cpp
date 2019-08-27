#include "Phase.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorPhase::CalibratorPhase(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mSnowThreshold(273.15),
      mRainThreshold(274.15),
      mMinPrecip(0.2) {
   iOptions.getRequiredValue("temperature", mTemperatureVariable);
   iOptions.getRequiredValue("precipitation", mPrecipitationVariable);
   iOptions.getValue("snowThreshold", mSnowThreshold);
   iOptions.getValue("rainThreshold", mRainThreshold);
   iOptions.check();
}
bool CalibratorPhase::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   int nTime = iFile.getNumTime();
   vec2 elevs = iFile.getElevs();


   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const FieldPtr temp = iFile.getField(mTemperatureVariable, t);
      const FieldPtr precip = iFile.getField(mPrecipitationVariable, t);
      FieldPtr phase = iFile.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               float currTemp     = (*temp)(i,j,e);
               float currPrecip   = (*precip)(i,j,e);
               if(Util::isValid(currTemp) && Util::isValid(currPrecip)) {
                  if(currPrecip <= mMinPrecip)
                     (*phase)(i,j,e)  = PhaseNone;
                  else if(!Util::isValid(currTemp))
                     (*phase)(i,j,e)  = Util::MV;
                  else if(currTemp <= mSnowThreshold)
                     (*phase)(i,j,e)  = PhaseSnow;
                  else if(currTemp <= mRainThreshold)
                     (*phase)(i,j,e)  = PhaseSleet;
                  else
                     (*phase)(i,j,e)  = PhaseRain;
               }
               else {
                  (*phase)(i,j,e) = Util::MV;
               }
            }
         }
      }
      Variable variable;
      iFile.addField(phase, mVariable, t);
   }
   return true;
}

std::string CalibratorPhase::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c phase", "Compute precipitation phase based on temperature, with values encoded by:") << std::endl;
      ss << Util::formatDescription("", "* 0 = no precipitation (Precip < minPrecip)") << std::endl;
      ss << Util::formatDescription("", "* 1 = rain (b < T)") << std::endl;
      ss << Util::formatDescription("", "* 2 = sleet (a < T <= b)") << std::endl;
      ss << Util::formatDescription("", "* 3 = snow (T < a)") << std::endl;
      ss << Util::formatDescription("", "T can be either a regular temperature or a wetbulb temperature") << std::endl;
      ss << Util::formatDescription("   temperature=required", "Name of temperature variable to use.") << std::endl;
      ss << Util::formatDescription("   precipitation=required", "Name of precipitation variable to use.") << std::endl;
      ss << Util::formatDescription("   minPrecip=0.2", "Minimum precip (in mm) needed to be considered as precipitation.") << std::endl;
      ss << Util::formatDescription("   snowThreshold=273.15", "Temperature threshold between snow and sleet.") << std::endl;
      ss << Util::formatDescription("   rainThreshold=274.15", "Temperature threshold between sleet and rain.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c phase", "Compute precipitation phase based on temperature") << std::endl;
   return ss.str();
}
