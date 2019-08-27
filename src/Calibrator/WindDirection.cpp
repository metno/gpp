#include "WindDirection.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Parameters.h"
CalibratorWindDirection::CalibratorWindDirection(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions) {
   iOptions.getRequiredValue("directionVariable", mDirectionVariable);
   iOptions.check();
}

bool CalibratorWindDirection::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(iParameterFile->getNumParameters() != 9) {
      Util::error("CalibratorWindDirection: ParameterFile must have 9 parameters");
   }

   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& wind      = *iFile.getField(mVariable, t);
      Field& direction = *iFile.getField(mDirectionVariable, t);

      Parameters parametersGlobal;
      if(!iParameterFile->isLocationDependent())
         parametersGlobal = iParameterFile->getParameters(t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            Parameters parameters;
            if(iParameterFile->isLocationDependent())
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
            else
               parameters = parametersGlobal;
            for(int e = 0; e < nEns; e++) {
               float currDirection = direction(i,j,e);
               if(Util::isValid(currDirection)) {
                   float factor = getFactor(currDirection, parameters);
                   wind(i,j,e) = factor*wind(i,j,e);
               }
            }
         }
      }
   }
   return true;
}

std::string CalibratorWindDirection::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c windDirection", "Multiply a variable by a factor based on the wind-direction:") << std::endl;
      ss << "                                factor = a + b*sin(dir)   + c*cos(dir)   + d*sin(2*dir) + e*cos(2*dir)" << std::endl;
      ss << "                                           + f*sin(3*dir) + g*cos(3*dir) + h*sin(4*dir) + i*cos(4*dir)" << std::endl;
      ss << Util::formatDescription("", "A parameter file is required, with the values [a b c d e f g h i].") << std::endl;
      ss << Util::formatDescription("   directionVariable=required", "Variable name to use as wind direction.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c windDirection", "Wind-direction-based bias-correction") << std::endl;
   return ss.str();
}

float CalibratorWindDirection::getFactor(float iWindDirection, const Parameters& iPar) {
   if(!Util::isValid(iWindDirection) || !iPar.isValid())
      return Util::MV;
   float dir = Util::deg2rad(iWindDirection);
   float factor = 1;
   factor = iPar[0] +
            iPar[1] * sin(dir) +
            iPar[2] * cos(dir) +
            iPar[3] * sin(2*dir) +
            iPar[4] * cos(2*dir) +
            iPar[5] * sin(3*dir) +
            iPar[6] * cos(3*dir) +
            iPar[7] * sin(4*dir) +
            iPar[8] * cos(4*dir);
   if(factor < 0)
      factor = 0;
   return factor;
}
