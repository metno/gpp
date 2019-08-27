#ifndef CALIBRATOR_QNH_H
#define CALIBRATOR_QNH_H
#include "Calibrator.h"
#include "../Variable.h"

//! Creates the QNH field by using surface pressure
class CalibratorQnh : public Calibrator {
   public:
      CalibratorQnh(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "qnh";};
      static float calcQnh(float iElev, float iPressure);
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::string mPressureVariable;
};
#endif
