#include "Upscale.h"
#include "../File/File.h"
#include "../Util.h"
#include "../KDTree.h"
#include <math.h>

// std::map<const File*, std::map<const File*, std::pair<vec2Int, vec2Int> > > DownscalerUpscale::mNeighbourCache;

DownscalerUpscale::DownscalerUpscale(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      mStatType(Util::StatTypeMean),
      mQuantile(0),
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   std::string statType;
   if(iOptions.getValue("stat", statType)) {
      bool status = Util::getStatType(statType, mStatType);
      if(!status) {
         std::stringstream ss;
         ss << "Could not recognize stat=" << statType;
         Util::error(ss.str());
      }
   }
   if(mStatType == Util::StatTypeQuantile) {
      iOptions.getRequiredValue("quantile", mQuantile);
      if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
         Util::error("'quantile' must be on the interval [0,1]");
      }
   }
   iOptions.check();
}

void DownscalerUpscale::downscaleCore(const File& iInput, File& iOutput) const {
   int nLat = iOutput.getNumY();
   int nLon = iOutput.getNumX();
   int nEns = iOutput.getNumEns();
   int nTime = iInput.getNumTime();

   // Create a map from each input point to the output grid
   vec2Int I, J;
   KDTree searchTree(iOutput.getLats(), iOutput.getLons());
   searchTree.getNearestNeighbour(iInput, I, J);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            vec2 values;
            values.resize(nEns);
            for(int ii = 0; ii < I.size(); ii++) {
               for(int jj = 0; jj < I[0].size(); jj++) {
                  if(I[ii][jj] == i && J[ii][jj] == j) {
                     for(int e = 0; e < nEns; e++) {
                        float curr = ifield(ii, jj, e);
                        if(Util::isValid(curr)) {
                           values[e].push_back(curr);
                        }
                     }
                  }
               }
            }
            for(int e = 0; e < nEns; e++) {
               ofield(i, j, e) = Util::calculateStat(values[e], mStatType, mQuantile);
            }
         }
      }
   }
}

std::string DownscalerUpscale::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-d upscale", "Aggregate all neighbours who are closest to the lookup point.") << std::endl;
      ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'quantile', 'std', or 'sum'. 'std' is the population standard deviation.") << std::endl;
      ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
   return ss.str();
}
