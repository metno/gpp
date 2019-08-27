#ifndef CALIBRATOR_MASK_H
#define CALIBRATOR_MASK_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Sets gridpoints that are far away from parameter points to missing
class CalibratorMask : public Calibrator {
   public:
      CalibratorMask(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "mask";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      bool mUseNearestOnly;
      bool mKeep;
};
#endif
