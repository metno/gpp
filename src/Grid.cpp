#include "Grid.h"

Grid::Grid() {

}

vec2 Grid::lats() const {
   return mLats;
}
vec2 Grid::lons() const {
   return mLons;
}
vec2 Grid::elevs() const {
   return mElevs;
}
vec2 Grid::landFractions() const {
   return mLandFractions;
}
void Grid::lats(vec2 iLats) {
   mLats = iLats;
}
void Grid::lons(vec2 iLons) {
   mLons = iLons;
}
void Grid::elevs(vec2 iElevs) {
   mElevs = iElevs;
}
void Grid::landFractions(vec2 iLandFractions) {
   mLandFractions = iLandFractions;
}
