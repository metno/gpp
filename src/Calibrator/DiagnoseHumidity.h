#ifndef CALIBRATOR_DIAGNOSE_HUMIDITY_H
#define CALIBRATOR_DIAGNOSE_HUMIDITY_H
#include "Calibrator.h"
#include "../Variable.h"

class CalibratorDiagnoseHumidity : public Calibrator {
   public:
      CalibratorDiagnoseHumidity(const Variable& iVariable, const Options& iOptions);
      std::string name() const {return "diagnoseHumidity";};
      bool requiresParameterFile() const { return false;};
      static std::string description(bool full=true);
      // Temperature in K, RH in [0,1], Returns TD in K
      static float computeDewpoint(float iTemperature, float iRelativeHumidity);
      // Temperatures in K, Returns RH in [0,1]
      static float computeRh(float iTemperature, float iDewPointTemperature);
      //! Compute wetbulb temperature
      //! @param iTemperature Temperature in K
      //! @param iPressure Pressure in pa
      //! @param iRelativeHumidity Relative humidity (out of 1)
      //! @return Wetbulb temperature in K
      static float computeWetbulb(float iTemperature, float iPressure, float iRelativeHumidity);
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      static float mEwt[41];
      std::string mTemperature;
      std::string mRh;
      std::string mDewpoint;
      std::string mPressure;
      std::string mCompute;
};
#endif
