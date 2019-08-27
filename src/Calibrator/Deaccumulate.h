#ifndef CALIBRATOR_DEACCUMULATE_H
#define CALIBRATOR_DEACCUMULATE_H
#include "Calibrator.h"
#include "../Variable.h"

// Deaccumlates a certain variable
class CalibratorDeaccumulate : public Calibrator {
   public:
      CalibratorDeaccumulate(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "deaccumulate";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
};
#endif
