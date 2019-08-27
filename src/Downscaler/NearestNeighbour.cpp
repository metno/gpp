#include "NearestNeighbour.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

// std::map<const File*, std::map<const File*, std::pair<vec2Int, vec2Int> > > DownscalerNearestNeighbour::mNeighbourCache;

DownscalerNearestNeighbour::DownscalerNearestNeighbour(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerNearestNeighbour::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumY();
   int nLon = iOutput.getNumX();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   double s = Util::clock();
   getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            int I = nearestI[i][j];
            int J = nearestJ[i][j];
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(I) && Util::isValid(J))
                  ofield(i,j,e) = ifield(I,J,e);
               else
                  ofield(i,j,e) = Util::MV;
            }
         }
      }
   }
}
void DownscalerNearestNeighbour::downscaleField(const Field& iInput, Field& iOutput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ) {

   int nLat = iOutput.getNumY();
   int nLon = iOutput.getNumX();
   int nEns = iOutput.getNumEns();

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         int I = nearestI[i][j];
         int J = nearestJ[i][j];
         for(int e = 0; e < nEns; e++) {
            if(Util::isValid(I) && Util::isValid(J)) {
               iOutput(i,j,e) = iInput(I, J, e);
            }
            else
               iOutput(i,j,e) = Util::MV;
         }
      }
   }
}

vec2 DownscalerNearestNeighbour::downscaleVec(const vec2& iInput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ) {

   int nLat = iOutputLats.size();
   int nLon = iOutputLats[0].size();

   assert(nearestI.size() == iOutputLats.size());
   assert(nearestJ.size() == iOutputLats.size());
   assert(nearestI[0].size() == iOutputLats[0].size());
   assert(nearestJ[0].size() == iOutputLats[0].size());

   vec2 output;
   output.resize(nLat);
   for(int i = 0; i < nLat; i++)
      output[i].resize(nLon);

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         assert(i < nearestI.size());
         assert(j < nearestI[i].size());
         assert(i < nearestJ.size());
         assert(j < nearestJ[i].size());
         int I = nearestI[i][j];
         int J = nearestJ[i][j];
         assert(I < iInput.size());
         assert(J < iInput[I].size());
         assert(I >= 0);
         assert(J >= 0);
         assert(i < output.size());
         assert(j < output[i].size());
         assert(i >= 0);
         assert(j >= 0);
         if(Util::isValid(I) && Util::isValid(J)) {
            output[i][j] = iInput[I][J];
         }
         else
            output[i][j] = Util::MV;
      }
   }
   return output;
}


std::string DownscalerNearestNeighbour::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-d nearestNeighbour", "Uses the nearest gridpoint in curved distance") << std::endl;
   return ss.str();
}
