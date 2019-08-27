#include <boost/scoped_ptr.hpp>
#include <cmath>

#include "Downscaler.h"
#include "../File/File.h"
#include "../KDTree.h"

std::map<Uuid, std::map<Uuid, std::pair<vec2Int, vec2Int> > > Downscaler::mNeighbourCache;

Downscaler::Downscaler(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) : Scheme(iOptions),
      mInputVariable(iInputVariable),
      mOutputVariable(iOutputVariable) {
}

bool Downscaler::downscale(const File& iInput, File& iOutput) const {
   if(iInput.getNumTime() != iOutput.getNumTime())
      return false;

   downscaleCore(iInput, iOutput);
   return true;
}

Downscaler* Downscaler::getScheme(std::string iName, const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) {
   Downscaler* d;
   if(iName == "nearestNeighbour") {
      d = new DownscalerNearestNeighbour(iInputVariable, iOutputVariable, iOptions);
   }
   else if(iName == "smart") {
      d = new DownscalerSmart(iInputVariable, iOutputVariable, iOptions);
   }
   else if(iName == "bypass") {
      d = new DownscalerBypass(iInputVariable, iOutputVariable, iOptions);
   }
   else if(iName == "pressure") {
      d = new DownscalerPressure(iInputVariable, iOutputVariable, iOptions);
   }
   else if(iName == "gradient") {
      d = new DownscalerGradient(iInputVariable, iOutputVariable, iOptions);
   }
   else if(iName == "bilinear") {
      d = new DownscalerBilinear(iInputVariable, iOutputVariable, iOptions);
   }
   else if(iName == "upscale") {
      d = new DownscalerUpscale(iInputVariable, iOutputVariable, iOptions);
   }
   else {
      Util::error("Could not instantiate downscaler of type '" + iName + "'");
      return NULL;
   }
   return d;
}

void Downscaler::getNearestNeighbourBruteForce(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(isCached(iFrom, iTo)) {
       getFromCache(iFrom, iTo, iI, iJ);
       return;
   }

   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumX();
   int nLat = iTo.getNumY();

   iI.resize(nLat);
   iJ.resize(nLat);

   // Check if the grid is the same
   if(iFrom.getNumY() == iTo.getNumY() && iFrom.getNumX() == iTo.getNumX()) {
      if(ilats == olats && ilons == olons) {
         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
         Util::info("Grids are identical, short cut in finding nearest neighbours");
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
      for(int j = 0; j < nLon; j++) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            getNearestNeighbourBruteForce(iFrom, olons[i][j], olats[i][j], iI[i][j], iJ[i][j]);
         }
      }
   }
   addToCache(iFrom, iTo, iI, iJ);
}

void Downscaler::getNearestNeighbour(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(isCached(iFrom, iTo)) {
       getFromCache(iFrom, iTo, iI, iJ);
       return;
   }

   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();
   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   int nLon = iTo.getNumX();
   int nLat = iTo.getNumY();

   // Check if the grid is the same
   if(iFrom.getNumY() == iTo.getNumY() && iFrom.getNumX() == iTo.getNumX()) {
      if(ilats == olats && ilons == olons) {

         iI.resize(nLat);
         iJ.resize(nLat);

         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, 0);
            iJ[i].resize(nLon, 0);
            for(int j = 0; j < nLon; j++) {
               iI[i][j] = i;
               iJ[i][j] = j;
            }
         }
         Util::info("Grids are identical, short cut in finding nearest neighbours");
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }
   // Check if the grid is regular lat/lon
   bool checkRegular = false;
   if(checkRegular) {
      bool isRegular = true;
      for(int y = 0; y < iFrom.getNumY(); y++) {
         for(int x = 0; x < iFrom.getNumX(); x++) {
            if(ilats[y][x] != ilats[y][0])
               isRegular = false;
            if(ilons[y][x] != ilons[0][x])
               isRegular = false;
         }
         if(!isRegular)
            break;
      }
      if(isRegular) {
         Util::info("Input grid is regular lat/lon, short cut in finding nearest neighbours");
         iI.resize(nLat);
         iJ.resize(nLat);

         std::vector<float> lats;
         std::vector<float> lons;
         lats.resize(iFrom.getNumY());
         lons.resize(iFrom.getNumX());
         for(int i = 0; i < iFrom.getNumY(); i++) {
            lats[i] = ilats[i][0];
         }
         for(int j = 0; j < iFrom.getNumX(); j++) {
            lons[j] = ilons[0][j];
         }

         for(int i = 0; i < nLat; i++) {
            iI[i].resize(nLon, Util::MV);
            iJ[i].resize(nLon, Util::MV);
            for(int j = 0; j < nLon; j++) {
               if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
                  int iNearest = Util::MV;
                  int jNearest = Util::MV;
                  float iMinDist = 1e9;
                  for(int ii = 0; ii < lats.size(); ii++) {
                     if(Util::isValid(lats[ii])) {
                        float dist = fabs(lats[ii] - olats[i][j]);
                        if(dist < iMinDist) {
                           iMinDist = dist;
                           iNearest = ii;
                        }
                     }
                  }
                  float jMinDist = 1e9;
                  for(int jj = 0; jj < lons.size(); jj++) {
                     if(Util::isValid(lons[jj])) {
                        float dist = fabs(lons[jj] - olons[i][j]);
                        if(dist < jMinDist) {
                           jMinDist = dist;
                           jNearest = jj;
                        }
                     }
                  }
                  if(Util::isValid(iNearest) && Util::isValid(jNearest)) {
                     iI[i][j] = iNearest;
                     iJ[i][j] = jNearest;
                  }
               }
            }
         }
         addToCache(iFrom, iTo, iI, iJ);
         return;
      }
   }

   KDTree searchTree(iFrom.getLats(), iFrom.getLons());
   searchTree.getNearestNeighbour(iTo, iI, iJ);

   addToCache(iFrom, iTo, iI, iJ);
}

bool Downscaler::isCached(const File& iFrom, const File& iTo) {
   std::map<Uuid, std::map<Uuid, std::pair<vec2Int, vec2Int> > >::const_iterator it = mNeighbourCache.find(iFrom.getUniqueTag());
   if(it == mNeighbourCache.end()) {
      return false;
   }
   std::map<Uuid, std::pair<vec2Int, vec2Int> >::const_iterator it2 = it->second.find(iTo.getUniqueTag());
   if(it2 == it->second.end()) {
      return false;
   }
   return true;
}

void Downscaler::addToCache(const File& iFrom, const File& iTo, vec2Int iI, vec2Int iJ) {
   std::pair<vec2Int, vec2Int> pair(iI, iJ);
   mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()] = pair;
}
bool Downscaler::getFromCache(const File& iFrom, const File& iTo, vec2Int& iI, vec2Int& iJ) {
   if(!isCached(iFrom, iTo))
      return false;
   iI = mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()].first;
   iJ = mNeighbourCache[iFrom.getUniqueTag()][iTo.getUniqueTag()].second;
   return true;
}

void Downscaler::getNearestNeighbourBruteForce(const File& iFrom, float iLon, float iLat, int& iI, int &iJ) {
   vec2 ilats = iFrom.getLats();
   vec2 ilons = iFrom.getLons();

   iI = Util::MV;
   iJ = Util::MV;

   if(Util::isValid(iLat) && Util::isValid(iLon)) {
      float minDist = Util::MV;
      for(int i = 0; i < ilats.size(); i++) {
         for(int j = 0; j < ilats[0].size(); j++) {
            if(Util::isValid(ilats[i][j]) && Util::isValid(ilons[i][j])) {
               float currDist = Util::getDistance(ilats[i][j], ilons[i][j], iLat, iLon);
               if(!Util::isValid(minDist) || currDist < minDist) {
                  iI = i;
                  iJ = j;
                  minDist = currDist;
               }
            }
         }
      }
   }
}

std::string Downscaler::getDescriptions(bool full) {
   std::stringstream ss;
   ss << DownscalerNearestNeighbour::description(full);
   ss << DownscalerGradient::description(full);
   ss << DownscalerBilinear::description(full);
   ss << DownscalerPressure::description(full);
   ss << DownscalerSmart::description(full);
   ss << DownscalerUpscale::description(full);
   ss << DownscalerBypass::description(full);
   return ss.str();
}

void Downscaler::clearCache() {
   mNeighbourCache.clear();
}
