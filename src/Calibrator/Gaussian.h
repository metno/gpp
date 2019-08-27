#ifndef CALIBRATOR_GAUSSIAN_H
#define CALIBRATOR_GAUSSIAN_H
#include "Calibrator.h"
#include "../Variable.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

class ParameterFile;
class Parameters;

//! Ensemble calibration using zero-adjusted gamma distribution. Its predictors are:
//! - ensemble mean
//! - ensemble fraction
//! Designed for precip
class CalibratorGaussian : public Calibrator {
   public:
      CalibratorGaussian(const Variable& iVariable, const Options& iOptions);
      static float getInvCdf(float iQuantile, float iEnsMean, float iEnsSpread, const Parameters& iParameters);
      static float getCdf(float iThreshold, float iEnsMean, float iEnsSpread, const Parameters& iParameters);
      static float getPdf(float iThreshold, float iEnsMean, float iEnsSpread, const Parameters& iParameters);

      static std::string description(bool full=true);
      std::string name() const {return "gaussian";};
      Parameters train(const std::vector<ObsEns>& iData) const;
   private:
      static double my_f(const gsl_vector *v, void *params);
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      int  mNeighbourhoodSize;
      float mLogLikelihoodTolerance;
      static const int mNumParameters = 2;
};
#endif
