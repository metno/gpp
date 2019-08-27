#include "Qnh.h"
#include "../Util.h"
#include "../File/File.h"
#include <math.h>
CalibratorQnh::CalibratorQnh(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.getRequiredValue("pressureVariable", mPressureVariable);
   iOptions.check();
}
bool CalibratorQnh::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   vec2 lats  = iFile.getLats();
   vec2 lons  = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const Field& input = *iFile.getField(mPressureVariable, t);
      Field& output      = *iFile.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            float currElev = elevs[i][j];
            for(int e = 0; e < nEns; e++) {
               float currPressure = input(i,j,e);
               if(Util::isValid(currPressure) && Util::isValid(currElev)) {
                  output(i,j,e)  = calcQnh(currElev, currPressure);
               }
            }
         }
      }
   }
   return true;
}
std::string CalibratorQnh::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-c qnh", "Adjusts the surface pressure down to sea-level based on a standard atmosphere (ICAO) producing the QNH variable.") << std::endl;
   else
      ss << Util::formatDescription("-c qnh", "Compute QNH") << std::endl;
   return ss.str();
}
float CalibratorQnh::calcQnh(float iElev, float iPressure) {
   if(iPressure == 0)
         return 0;
   else if(Util::isValid(iElev) && Util::isValid(iPressure)) {
      float g  = 9.80665;   // m/s2
      float T0 = 288.15;    // K
      float L  = 0.0065;    // K/m
      // Method 1:
      // float dElev = 0 - iElev;
      // float M  = 0.0289644; // kg/mol
      // float R  = 8.31447;   // J/(mol•K)
      // float cp = 1007;      // J/(kg•K)
      // float constant = -g*M/R/T0;
      // float qnh = iPressure * pow(1 - L*dElev/T0, g*M/R/L);

      // Method 2: http://www.hochwarth.com/misc/AviationCalculator.html
      float CRGas = 287.053; // [m^2/(s^2*K)] = [J/(kg*K)]
      float p0    = 101325;  // pa
      float qnh   = p0*pow(pow((iPressure/p0), (CRGas*L)/g) + (iElev*L)/T0, g/(CRGas*L));
      return qnh;
   }
   else {
      return Util::MV;
   }
}
