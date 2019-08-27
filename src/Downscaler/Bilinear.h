#ifndef DOWNSCALER_BILINEAR_H
#define DOWNSCALER_BILINEAR_H
#include <map>
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
#include "../Field.h"
class File;
class DownscalerBilinear : public Downscaler {
   public:
      DownscalerBilinear(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);

      //! Bilinearly interpolate four irregularly placed points:
      //! v1 v3
      //! v0 v2
      //! with corresponding x,y coordinates.
      //! @param x Interpolate to this x-coordinate
      //! @param y Interpolate to this y-coordinate
      //! Two or more points cannot be colocated
      static float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3);

      // Interpolate a whole field
      static void downscaleField(const Field& iInput, Field& iOutput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ);

      // Interpolate a whole vec2
      static vec2 downscaleVec(const vec2& iInput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ);

      //! Find which I/J coordinates surround a lookup point
      //! Returns false if the lookup point is outside the grid. In this case, I1, I2, J1, J2 values
      //! cannot be used.
      static bool findCoords(float iLat, float iLon, const vec2& iLats, const vec2& iLons, int I, int J, int& I1, int& J1, int& I2, int& J2);

      static std::string description(bool full=true);
      std::string name() const {return "bilinear";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      static bool calcParallelogram(float x, float y, float X1, float X2, float X3, float X4, float Y1, float Y2, float Y3, float Y4, float &t, float &s);
      static bool calcGeneral(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float &t, float &s);
};
#endif
