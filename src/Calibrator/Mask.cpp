#include "Mask.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
CalibratorMask::CalibratorMask(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mUseNearestOnly(false),
      mKeep(true) {
   iOptions.getValue("keep", mKeep);
   iOptions.check();
}
bool CalibratorMask::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Cannot use a location independent parameter file in mask calibrator");
   }
   if(!iParameterFile->isFixedSize()) {
      Util::error("Cannot use a parameter file without a constant number of parameters");
   }
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 isWithinRadius;

   // Stores if a certain gridpoint is within some radius of influence
   isWithinRadius.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      isWithinRadius[i].resize(nLon, 0);
   }

   for(int t = 0; t < nTime; t++) {
      FieldPtr field = iFile.getField(mVariable, t);
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            if(t == 0 || iParameterFile->isTimeDependent()) {
               Location loc(Util::MV, Util::MV, Util::MV);
               Location currLocation(lats[i][j], lons[i][j], elevs[i][j]);
               if(mUseNearestOnly) {
                  iParameterFile->getNearestLocation(t, currLocation, loc);
                  Parameters parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
                  float dist = currLocation.getDistance(loc);
                  if(Util::isValid(dist) && dist < parameters[0]) {
                     isWithinRadius[i][j] = 1;
                  }
               }
               else {
                  std::vector<Location> locations = iParameterFile->getLocations();
                  for(int k = 0; k < locations.size(); k++) {
                     Parameters parameters = iParameterFile->getParameters(t, locations[k]);
                     float dist = currLocation.getDistance(locations[k]);
                     if(Util::isValid(dist) && dist < parameters[0]) {
                        isWithinRadius[i][j] = 1;
                     }
                  }
               }
            }
            // Remove the point if keep and we are not within the radius or we don't keep and are
            // inside
            bool remove = mKeep != isWithinRadius[i][j];
            if(remove) {
               for(int e = 0; e < nEns; e++) {
                  (*field)(i,j,e) = Util::MV;
               }
            }
         }
      }
   }
   iFile.setElevs(elevs);
   return true;
}

std::string CalibratorMask::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c mask", "Allows removal of gridpoints that are not within (or alternatively within) a radius of a set of points provided in a parameter file. The removed gridpoints are set to missing. The parameter file must be spatial and conain one parameter representing the radius in meters. Each gridpoint is checked against each parameter point to see if the gridpoint is within its radius.") << std::endl;
      ss << Util::formatDescription("   keep=1", "If 1, then all gridpoints within the radius of any parameter point will be kept. If 0, then all gridpoints within the radius of any parameter point will be removed.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c mask", "Remove or keep gridpoints near parameter locations") << std::endl;
   return ss.str();
}
