#include "Fake.h"
FileFake::FileFake(const Options& iOptions) :
      File("", iOptions) {
   int nLat = 10;
   int nLon = 10;
   mNEns = 2;
   int nTime = 10;

   iOptions.getValue("nLat", nLat);
   iOptions.getValue("nLon", nLon);
   iOptions.getValue("nEns", mNEns);
   iOptions.getValue("nTime", nTime);
   if(!Util::isValid(nLat) || nLat <= 0) {
      Util::error("FileFake: Invalid number of latitudes");
   }
   if(!Util::isValid(nLat) || nLon <= 0) {
      Util::error("FileFake: Invalid number of longitudes");
   }
   if(!Util::isValid(mNEns) || mNEns <= 0) {
      Util::error("FileFake: Invalid number of ensemble members");
   }
   if(!Util::isValid(nTime) || nTime <= 0) {
      Util::error("FileFake: Invalid number of times");
   }
   vec2 lats;
   lats.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      lats[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         lats[i][j] = 50 + 10.0 * i / nLat;
      }
   }
   vec2 lons;
   lons.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      lons[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         lons[i][j] = 0 + 10.0 * j / nLon;
      }
   }
   bool successLats = setLats(lats);
   if(!successLats) {
      std::stringstream ss;
      ss << "Could not set latitudes in " << getFilename();
      Util::error(ss.str());
   }
   bool successLons = setLons(lons);
   if(!successLons) {
      std::stringstream ss;
      ss << "Could not set longitudes in " << getFilename();
      Util::error(ss.str());
   }

   vec2 elevs;
   elevs.resize(getNumY());
   for(int i = 0; i < nLat; i++) {
      elevs[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         elevs[i][j] = 0;
      }
   }
   setElevs(elevs);
   mLandFractions.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      mLandFractions[i].resize(nLon);
      for(int j = 0; j < nLon; j++) {
         mLandFractions[i][j] = (float) i / nLat / 2 + (float) j / nLon;
      }
   }
   std::vector<double> times(nTime, 0);
   for(int i = 0; i < nTime; i++)
      times[i] = i;
   setTimes(times);
}

FieldPtr FileFake::getFieldCore(const Variable& iVariable, int iTime) const {
   FieldPtr field = getEmptyField();

   for(int i = 0; i < getNumY(); i++) {
      for(int j = 0; j < getNumX(); j++) {
         for(int e = 0; e < getNumEns(); e++) {
            (*field)(i,j,e) = i + j + e;
         }
      }
   }
   return field;
}
void FileFake::writeCore(std::vector<Variable> iVariables, std::string iMessage) {
   Util::warning("Cannot write file using the 'Fake' format");
}

bool FileFake::hasVariableCore(const Variable& iVariable) const {
   return true;
}
std::string FileFake::description() {
   std::stringstream ss;
   ss << "   type=fake                    Fake file" << std::endl;
   return ss.str();
}
