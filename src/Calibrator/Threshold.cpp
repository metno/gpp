#include "Threshold.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
CalibratorThreshold::CalibratorThreshold(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.getRequiredValues("thresholds", mThresholds);
   iOptions.getRequiredValues("values", mValues);
   if(!iOptions.getValues("equals", mEquals)) {
      mEquals.resize(mThresholds.size(), 0);
   }
   if(mValues.size() != mThresholds.size() + 1) {
      std::stringstream ss;
      ss << "Length of 'values' must be one longer than the length of 'thresholds'";
      Util::error(ss.str());
   }
   if(mEquals.size() != mThresholds.size()) {
      std::stringstream ss;
      ss << "Length of 'equals' must be the same as the length of 'thresholds'";
      Util::error(ss.str());
   }
   iOptions.check();
}
bool CalibratorThreshold::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   int nX = iFile.getNumX();
   int nY = iFile.getNumY();
   vec2 elevs = iFile.getElevs();
   int nThresholds = mThresholds.size();

   for(int t = 0; t < nTime; t++) {
      const FieldPtr field = iFile.getField(mVariable, t);
      for(int e = 0; e < nEns; e++) {
         for(int y = 0; y < nY; y++) {
            for(int x = 0; x < nX; x++) {
               float value = (*field)(y, x, e);
               if(Util::isValid(value)) {
                  (*field)(y, x, e) = mValues[nThresholds];
                  for(int p = 0; p < nThresholds; p++) {
                     if(value < mThresholds[p]) {
                        (*field)(y, x, e) = mValues[p];
                        break;
                     }
                     else if(value == mThresholds[p] && mEquals[p] == 1) {
                        (*field)(y, x, e) = mValues[p];
                        break;
                     }
                  }
               }
            }
         }
      }
   }
   return true;
}

std::string CalibratorThreshold::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-c threshold", "Apply thresholding to the field, converting ranges into specified values.") << std::endl;
   if(full) {
      ss << Util::formatDescription("   thresholds=required", "List of thresholds (x1, x2, x3, ..., xn). Must be sorted.") << std::endl;
      ss << Util::formatDescription("   values=required", "List of values to set the ranges to. The first value is used for field values below the first threshold. Must have length one longer than threhsolds") << std::endl;
      ss << Util::formatDescription("   equals=undef", "Should the intervals include the upper threshold? One number for each threshold. 1 means yes, 0 means no.") << std::endl;
   }
   return ss.str();
}
