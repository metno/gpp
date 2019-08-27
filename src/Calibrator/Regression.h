#ifndef CALIBRATOR_REGRESSION_H
#define CALIBRATOR_REGRESSION_H
#include "Calibrator.h"
class ParameterFile;
class Parameters;

//! Applies polynomial regression to forecasts
class CalibratorRegression : public Calibrator {
   public:
      CalibratorRegression(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "regression";};
      Parameters train(const std::vector<ObsEns>& iData) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      int mOrder;
      bool mIntercept;
      std::vector<std::string> mVariables;
};
#endif
