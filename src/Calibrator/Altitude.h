#ifndef CALIBRATOR_ALTITUDE_H
#define CALIBRATOR_ALTITUDE_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Changes the altitudes in iFile to the altitudes in the parameter file
class CalibratorAltitude : public Calibrator {
   public:
      CalibratorAltitude(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "altitude";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
};
#endif
