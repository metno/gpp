#ifndef GRID_H
#define GRID_H
#include <vector>
typedef std::vector<std::vector<float> > vec2;

class Grid {
   public:
      Grid();
      vec2 lats() const;
      vec2 lons() const;
      vec2 elevs() const;
      vec2 landFractions() const;
      void lats(vec2 iLats);
      void lons(vec2 iLons);
      void elevs(vec2 iElevs);
      void landFractions(vec2 iLandFractions);
   private:
      vec2 mLats;
      vec2 mLons;
      vec2 mElevs;
      vec2 mLandFractions;
};
#endif
