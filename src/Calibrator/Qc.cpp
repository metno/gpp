#include "Qc.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
CalibratorQc::CalibratorQc(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mMin(Util::MV),
      mMax(Util::MV) {
   iOptions.getValue("min", mMin);
   iOptions.getValue("max", mMax);
   if(!Util::isValid(mMin) && !Util::isValid(mMax))
      Util::warning("CalibratorQc: both 'min' and 'max' are missing, therefore no correction is applied.");
   iOptions.check();
}

bool CalibratorQc::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& field = *iFile.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {

            for(int e = 0; e < nEns; e++) {
               float value = field(i,j,e);
               if(Util::isValid(value)) {
                  if(Util::isValid(mMin) && value < mMin)
                     value = mMin;
                  else if(Util::isValid(mMax) && value > mMax)
                     value = mMax;
                  field(i,j,e) = value;
               }
            }
         }
      }
   }
   return true;
}

std::string CalibratorQc::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c qc", "Apply quality control, ensuring the values are within the bounds of the variable. If the original value is missing, no correction is applied.") << std::endl;
      ss << Util::formatDescription("   min=undef", "Force the minimum allowed value. If not provided, don't force a minimum.") << std::endl;
      ss << Util::formatDescription("   max=undef", "Force the maximum allowed value. If not provided, don't force a maximum.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c qc", "Ensures values are within bounds") << std::endl;
   return ss.str();
}
