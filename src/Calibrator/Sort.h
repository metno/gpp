#ifndef CALIBRATOR_SORT_H
#define CALIBRATOR_SORT_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Applies polynomial regression to forecasts
class CalibratorSort : public Calibrator {
   public:
      CalibratorSort(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "sort";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
};
#endif
