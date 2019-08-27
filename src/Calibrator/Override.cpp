#include "Override.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorOverride::CalibratorOverride(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mRadius(0),
      mMaxElevDiff(Util::MV) {
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("maxElevDiff", mMaxElevDiff);
   iOptions.check();
}
bool CalibratorOverride::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   int nX = iFile.getNumX();
   int nY = iFile.getNumY();
   vec2 elevs = iFile.getElevs();

   if(iParameterFile->getNumParameters() != 1) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have one dataacolumn");
   }
   if(!iParameterFile->isLocationDependent()) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must be location dependent");
   }

   std::vector<Location> pointLocations = iParameterFile->getLocations();
   KDTree searchTree(iFile.getLats(), iFile.getLons());
   std::vector<int> Ys(pointLocations.size(), 0);
   std::vector<int> Xs(pointLocations.size(), 0);
   for(int k = 0; k < pointLocations.size(); k++) {
      searchTree.getNearestNeighbour(pointLocations[k].lat(), pointLocations[k].lon(), Ys[k], Xs[k]);
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const FieldPtr field = iFile.getField(mVariable, t);
      for(int k = 0; k < pointLocations.size(); k++) {
         float currElev = pointLocations[k].elev();
         Parameters parameters = iParameterFile->getParameters(t, pointLocations[k]);
         float value = parameters[0];
         int Ynn = Ys[k];
         int Xnn = Xs[k];
         if(Util::isValid(value)) {
            if(Ynn > 0 && Ynn < nY - 1 && Xnn > 0 && Xnn < nX - 1) {
               // Don't do this for points that are on the boundary, because then it is possible that
               // they are outside the domain. This means that some points that are ligitimately inside
               // but near the boundary may be removed.
               for(int e = 0; e < nEns; e++) {
                  for(int y = std::max(0, Ynn - mRadius); y <= std::min(nY-1, Ynn + mRadius); y++) {
                     for(int x = std::max(0, Xnn - mRadius); x <= std::min(nX-1, Xnn + mRadius); x++) {
                        if(!Util::isValid(mMaxElevDiff)) {
                           // Ignore elevation information
                           (*field)(y, x, e) = value;
                        }
                        else if(Util::isValid(elevs[y][x]) && Util::isValid(currElev)) {
                           // Check if within elevation bound
                           float elevDiff = abs(elevs[y][x] - currElev);
                           if(elevDiff < mMaxElevDiff) {
                              (*field)(y, x, e) = value;
                           }
                        }
                        else {
                           // Do not update
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return true;
}

std::string CalibratorOverride::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-c override", "Overrides points in the file based on values at point locations.") << std::endl;
   if(full) {
      ss << Util::formatDescription("   radius=0", "Write values into all gridpoints in a square box with this radius (in number of gridpoints).") << std::endl;
      ss << Util::formatDescription("   maxElevDiff=undef", "Only write in gridpoints that are within this elevation difference to the point. If undefined, then do not do an elevation check.") << std::endl;
   }
   return ss.str();
}
