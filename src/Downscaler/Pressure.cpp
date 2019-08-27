#include "Pressure.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>
const float DownscalerPressure::mConstant = -1.21e-4;

DownscalerPressure::DownscalerPressure(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerPressure::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumY();
   int nLon = iOutput.getNumX();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   vec2 ilats  = iInput.getLats();
   vec2 ilons  = iInput.getLons();
   vec2 ielevs = iInput.getElevs();
   vec2 olats  = iOutput.getLats();
   vec2 olons  = iOutput.getLons();
   vec2 oelevs = iOutput.getElevs();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int Icenter = nearestI[i][j];
            int Jcenter = nearestJ[i][j];
            assert(Icenter < ielevs.size());
            assert(Jcenter < ielevs[Icenter].size());
            for(int e = 0; e < nEns; e++) {
               float currElev = oelevs[i][j];
               float nearestElev = ielevs[Icenter][Jcenter];
               if(!Util::isValid(currElev) || !Util::isValid(nearestElev)) {
                  // Can't adjust if we don't have an elevation, use nearest neighbour
                  ofield(i,j,e) = ifield(Icenter,Jcenter,e);
               }
               else {
                  float nearestPressure = ifield(Icenter,Jcenter,e);
                  float currPressure = calcPressure(nearestElev, nearestPressure, currElev);
                  ofield(i,j,e) = currPressure;
               }
            }
         }
      }
   }
}
std::string DownscalerPressure::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-d pressure", "Adjusts the pressure of the nearest neighbour based on the elevation difference and a standard atmosphere.") << std::endl;
   else
      ss << Util::formatDescription("-d pressure", "Adjusts the pressure based on the elevation differences") << std::endl;
   return ss.str();
}

float DownscalerPressure::calcPressure(float iElev0, float iPressure0, float iElev1) {
   if(Util::isValid(iElev0) && Util::isValid(iPressure0) && Util::isValid(iElev1)) {
      float dElev = iElev1 - iElev0;
      return iPressure0 * exp(mConstant * (dElev));
   }
   else {
      return Util::MV;
   }
}
