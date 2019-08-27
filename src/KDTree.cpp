#include "KDTree.h"
KDTree::KDTree() : mRoot(NULL) {
}

KDTree::KDTree(const vec2& iLats, const vec2& iLons) {
   build(iLats, iLons);
}

KDTree& KDTree::operator=(KDTree other) {
   std::swap(mLats, other.mLats);
   std::swap(mLons, other.mLons);

   mRoot.swap(other.mRoot);

   return *this;
}

KDTree::KDTree(const KDTree& other) {
   if(other.mRoot != NULL) {
      mLats = other.mLats;
      mLons = other.mLons;
      mRoot.reset(new KDTree::TreeNode(*other.mRoot));
   }
}

void KDTree::build(const vec2& iLats, const vec2& iLons) {
   mLats = iLats;
   mLons = iLons;
   if(iLats.size() != iLons.size())
      Util::error("Cannot initialize KDTree, lats and lons not the same size");

   size_t nLat = iLats.size();
   if(nLat == 0)
      Util::error("Cannot initialize KDTree, no valid locations");

   size_t nLon = iLats[0].size();
   if(nLon == 0)
      Util::error("Cannot initialize KDTree, no valid locations");

   indexdVec lons(nLon*nLat);

   indexdVec::iterator currLon = lons.begin();
   size_t to = -1;
   for(size_t i = 0; i < iLats.size(); ++i) {
     for(size_t j = 0; j < iLats[0].size(); ++j) {
        if(Util::isValid(iLons[i][j]) && Util::isValid(iLats[i][j])) {
           *(currLon++) = Indexed(iLons[i][j], iLats[i][j], i, j);
           ++to;
        }
     }
   }

   if(to == -1) {
      Util::error("Cannot initialize KDTree, no valid locations");
   }
   if(to >= 0) subTree(lons, 0, to, true, NULL, mRoot);
}

bool compareLons (const KDTree::Indexed& l, const KDTree::Indexed& r){ return (l.lon < r.lon); }
bool compareLats (const KDTree::Indexed& l, const KDTree::Indexed& r){ return (l.lat < r.lat); }

void KDTree::subTree(indexdVec& iLonLat,
                     const size_t from,
                     const size_t to,
                     const bool xsection,
                     const TreeNode* parent,
                     unode& root) {

   root.reset(new KDTree::TreeNode());

   root->parent = parent;
   root->xsection = xsection;

   size_t len = to - from + 1;

   if(len == 1) {
      root->lon = iLonLat[from].lon;
      root->lat = iLonLat[from].lat;
      root->ipos = iLonLat[from].index1;
      root->jpos = iLonLat[from].index2;
      return;
   }

   size_t med = from + len/2;
   if(xsection) {
      std::sort(iLonLat.begin() + from, iLonLat.begin() + to + 1, compareLons);
      while(med < to && iLonLat[med].lon == iLonLat[med+1].lon) {
         assert(iLonLat.size() > med+1);
         ++med;
      }
   }
   else {
      std::sort(iLonLat.begin() + from, iLonLat.begin() + to + 1, compareLats);
      while(med < to && iLonLat[med].lat == iLonLat[med+1].lat) {
         assert(iLonLat.size() > med+1);
         ++med;
      }
   }

   root->lon = iLonLat[med].lon;
   root->lat = iLonLat[med].lat;
   root->ipos = iLonLat[med].index1;
   root->jpos = iLonLat[med].index2;

   if(med > from) subTree(iLonLat, from, med - 1, !xsection, root.get(), root->left);
   if(med < to) subTree(iLonLat, med + 1, to, !xsection, root.get(), root->right);

}

const KDTree::TreeNode* KDTree::firstGuess(const unode& root, const float lon, const float lat) {
   if(root->xsection){
      if(lon <= root->lon) {
         if(root->left) return firstGuess(root->left, lon, lat);
      } else {
         if(root->right) return firstGuess(root->right, lon, lat);
      }
   } else {
      if(lat <= root->lat) {
         if(root->left) return firstGuess(root->left, lon, lat);
      } else {
         if(root->right) return firstGuess(root->right, lon, lat);
      }
   }
   return root.get();
}

const KDTree::TreeNode* KDTree::nearestNeighbour(const unode& root, const float lon, const float lat) {

   const TreeNode * nearLeaf = firstGuess(root, lon, lat);

   float currBest = Util::getDistance(lat, lon, nearLeaf->lat, nearLeaf->lon);

   const TreeNode * currentLeaf = nearLeaf;
   const TreeNode * nearestLeaf = nearLeaf;

   while (currentLeaf) {
      float dist = Util::getDistance(lat, lon, currentLeaf->lat, currentLeaf->lon);

      if(dist < currBest) {
         currBest = dist;
         nearestLeaf = currentLeaf;
      }

      bool cross;
      if(currentLeaf->xsection) {
         cross = Util::getDistance(lat, lon, lat, currentLeaf->lon) < currBest;
      } else {
         cross = Util::getDistance(lat, lon, currentLeaf->lat, lon) < currBest;
      }

      if(cross) {
         if(currentLeaf->left.get() == nearLeaf) {
            nearLeaf = NULL;
            if(currentLeaf->right) nearLeaf = nearestNeighbour( currentLeaf->right, lon, lat);
         } else {
            nearLeaf = NULL;
            if(currentLeaf->left) nearLeaf = nearestNeighbour(currentLeaf->left, lon, lat);
         }

         if(nearLeaf) {
            float dist = Util::getDistance(lat, lon, nearLeaf->lat, nearLeaf->lon);
            if(dist < currBest) {
               currBest = dist;
               nearestLeaf = nearLeaf;
            }
         }
      }

      nearLeaf = currentLeaf;
      if(currentLeaf == root.get()) {
         currentLeaf = NULL;
      } else {
         currentLeaf = currentLeaf->parent;
      }
   }

   return nearestLeaf;
}

void KDTree::getNearestNeighbour(const File& iTo, vec2Int& iI, vec2Int& iJ) const {

   vec2 olats = iTo.getLats();
   vec2 olons = iTo.getLons();
   size_t nLon = iTo.getNumX();
   size_t nLat = iTo.getNumY();

   iI.resize(nLat);
   iJ.resize(nLat);

   if(!mRoot) {
      for(size_t i = 0; i < nLat; ++i) {
         iI[i].resize(nLon, Util::MV);
         iJ[i].resize(nLon, Util::MV);
      }
      return;
   }

   const TreeNode * nearest;

   for(size_t i = 0; i < nLat; ++i) {
      iI[i].resize(nLon, Util::MV);
      iJ[i].resize(nLon, Util::MV);
   }
   #pragma omp parallel for private(nearest)
   for(size_t i = 0; i < nLat; ++i) {
      for(size_t j = 0; j < nLon; ++j) {
         if(Util::isValid(olats[i][j]) && Util::isValid(olons[i][j])) {
            // Find the nearest neighbour from input grid (ii, jj)
            nearest = nearestNeighbour(mRoot, olons[i][j], olats[i][j]);
            iI[i][j] = nearest->ipos;
            iJ[i][j] = nearest->jpos;
         }
      }
   }
}

void KDTree::getNearestNeighbour(float iLat, float iLon, int& iI, int& iJ) const {
   const TreeNode * nearest = nearestNeighbour(mRoot, iLon, iLat);
   iI = nearest->ipos;
   iJ = nearest->jpos;
   return;
}
