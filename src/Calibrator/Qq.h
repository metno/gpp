#ifndef CALIBRATOR_QQ_H
#define CALIBRATOR_QQ_H
#include "Calibrator.h"
class ParameterFile;

//! Applies polynomial regression to forecasts
class CalibratorQq : public Calibrator {
   public:
      CalibratorQq(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "qq";};
      Parameters train(const std::vector<ObsEns>& iData) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      float mLowerQuantile;
      float mUpperQuantile;
      std::vector<float> mQuantiles;
      struct ExtrapolationPolicy {
         enum Policy {
            OneToOne = 0,
            MeanSlope = 10,
            NearestSlope = 20,
            Zero = 30,
         };
      };
      ExtrapolationPolicy::Policy mPolicy;
      // Separate the vector of obs,fcst,obs,fcst,... into separate vectors
      void separate(const Parameters& iParameters, std::vector<float>& iObs, std::vector<float>& iFcst) const;
      std::vector<float> mExtraObs;
      std::vector<float> mExtraFcst;
};
#endif
