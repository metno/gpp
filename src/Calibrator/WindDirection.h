#ifndef CALIBRATOR_WIND_DIRECTION_H
#define CALIBRATOR_WIND_DIRECTION_H
#include "Calibrator.h"
#include "../Variable.h"

class ParameterFile;
class Parameters;

//! Multiply a variable by a factor based on the wind-direction
//! factor = a + b*sin(dir)   + c*cos(dir)   + d*sin(2*dir) + e*cos(2*dir)
//!            + f*sin(3*dir) + g*cos(3*dir) + h*sin(4*dir) + i*cos(4*dir)
class CalibratorWindDirection : public Calibrator {
   public:
      CalibratorWindDirection(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "windDirection";};
      //! Get multiplication factor for given wind direction
      //! @param iWindDirection in degrees, meteorological wind direction (0 degrees is from North)
      static float getFactor(float iWindDirection, const Parameters& iPar);
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      std::string mDirectionVariable;
};
#endif
