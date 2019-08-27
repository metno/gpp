#include "Sort.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorSort::CalibratorSort(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.check();
}
bool CalibratorSort::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& field = *iFile.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            std::vector<float> values = field(i,j);
            std::sort(values.begin(), values.end());

            // Create a new array with all the missing values
            // at the end
            std::vector<float> missingLast(nEns, Util::MV);
            int counter = 0;
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(values[e])) {
                  missingLast[counter] = values[e];
                  counter++;
               }
            }
            for(int e = 0; e < nEns; e++) {
               field(i,j,e) = missingLast[e];
            }
         }
      }
   }
   return true;
}

std::string CalibratorSort::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-c sort", "Sorts ensemble members from lowest to highest. Any missing members are placed last.") << std::endl;
   else
      ss << Util::formatDescription("-c sort", "Sorts ensemble members from lowest to highest") << std::endl;
   return ss.str();
}
