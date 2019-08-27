#ifndef CALIBRATOR_OVERRIDE_H
#define CALIBRATOR_OVERRIDE_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Overrides values from parameter file into gridpoints
class CalibratorOverride : public Calibrator {
   public:
      CalibratorOverride(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "override";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      int mRadius;
      float mMaxElevDiff;
};
#endif
