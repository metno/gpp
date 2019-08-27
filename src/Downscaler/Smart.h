#ifndef DOWNSCALER_SMART_H
#define DOWNSCALER_SMART_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
typedef std::vector<std::vector<std::vector<int> > > vec3Int; // lat, lon, neighbour index

//! Downscale using neighbours at similar elevation
class DownscalerSmart : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerSmart(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "smart";};

      //! Method may return fewer than num smart neighbours
      void getSmartNeighbours(const File& iFrom, const File& iTo, vec3Int& iI, vec3Int& iJ) const;
      static int getNumSearchPoints(int iSearchRadius) ;
             int getNumSearchPoints() const;
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      int mRadius;
      int mNum;
      float mMinElevDiff;
};
#endif
