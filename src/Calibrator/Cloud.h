#ifndef CALIBRATOR_CLOUD_H
#define CALIBRATOR_CLOUD_H
#include "Calibrator.h"
#include "../Variable.h"

//! Ensures that if a member has precip it also has full cloud cover
class CalibratorCloud : public Calibrator {
   public:
      CalibratorCloud(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "cloud";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::string mPrecipVariable;
      float mValue;
};
#endif
