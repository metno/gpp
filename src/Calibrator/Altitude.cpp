#include "Altitude.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
CalibratorAltitude::CalibratorAltitude(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.check();
}
bool CalibratorAltitude::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Cannot use a location independent parameter file to update the altitudes");
   }
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         Location loc(Util::MV, Util::MV, Util::MV);
         iParameterFile->getNearestLocation(0, Location(lats[i][j], lons[i][j], elevs[i][j]), loc);
         elevs[i][j] = loc.elev();
      }
   }
   iFile.setElevs(elevs);
   return true;
}

std::string CalibratorAltitude::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-c altitude", "Changes the altitudes to the altitudes in the parameter file. If the file does not have an altitude field, then a new one is not created.") << std::endl;
   else
      ss << Util::formatDescription("-c altitude", "Changes the altitudes to the altitudes in the parameter file") << std::endl;
   return ss.str();
}
