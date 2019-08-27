#ifndef CALIBRATOR_THRESHOLD_H
#define CALIBRATOR_THRESHOLD_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Thresholds values from parameter file into gridpoints
class CalibratorThreshold : public Calibrator {
   public:
      CalibratorThreshold(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "threshold";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::vector<float> mThresholds;
      std::vector<float> mValues;
      std::vector<int> mEquals;
};
#endif
