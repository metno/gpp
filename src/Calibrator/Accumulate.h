#ifndef CALIBRATOR_ACCUMULATE_H
#define CALIBRATOR_ACCUMULATE_H
#include "Calibrator.h"
#include "../Variable.h"

// Accumlates a certain variable
class CalibratorAccumulate : public Calibrator {
   public:
      CalibratorAccumulate(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "accumulate";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
};
#endif
