#ifndef CALIBRATOR_PHASE_H
#define CALIBRATOR_PHASE_H
#include "Calibrator.h"
class ParameterFile;

//! Creates a precipitation-phase field
class CalibratorPhase : public Calibrator {
   public:
      CalibratorPhase(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "phase";};

      //! Precipitation phase
      enum Phase {
         PhaseNone  = 0,
         PhaseRain  = 1,
         PhaseSleet = 2,
         PhaseSnow  = 3
      };
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      float mMinPrecip;
      float mSnowThreshold;
      float mRainThreshold;
      std::string mTemperatureVariable;
      std::string mPrecipitationVariable;
};
#endif
