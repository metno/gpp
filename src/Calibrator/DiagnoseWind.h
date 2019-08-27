#ifndef CALIBRATOR_DIAGNOSE_WIND_H
#define CALIBRATOR_DIAGNOSE_WIND_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnoseWind : public Calibrator {
   public:
      CalibratorDiagnoseWind(const Variable& iVariable, const Options& iOptions);
      std::string name() const {return "diagnoseWind";};
      bool requiresParameterFile() const { return false;};
      static std::string description(bool full=true);
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::string mX;
      std::string mY;
      std::string mSpeed;
      std::string mDirection;
      std::string mCompute;
};
#endif
