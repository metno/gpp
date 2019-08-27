#ifndef CALIBRATOR_BCT_H
#define CALIBRATOR_BCT_H
#include "Calibrator.h"
#include "../Variable.h"

class ParameterFile;
class Parameters;

//! Ensemble calibration using a Box-Cox t-distribution. Its predictors are:
//! - ensemble mean
//! - ensemble standard deviation
//! Designed for windspeed
class CalibratorBct : public Calibrator {
   public:
      CalibratorBct(const Variable& iMainPredictor, const Options& iOptions);
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsStd, const Parameters& iParameters);

      static std::string description(bool full=true);
      std::string name() const {return "bct";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      float mMaxEnsMean;
};
#endif
