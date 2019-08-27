#include "Point.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"

FilePoint::FilePoint(std::string iFilename, const Options& iOptions) :
      File(iFilename, iOptions) {
   float lat,lon,elev;
   if(!iOptions.getValue("lat", lat)) {
      Util::error("Missing 'lat' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("lon", lon)) {
      Util::error("Missing 'lon' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("elev", elev)) {
      Util::error("Missing 'elev' option for '" + iFilename + "'");
   }
   if(lat < -90 || lat > 90) {
      std::stringstream ss;
      ss << "Invalid latitude: " << lat;
      Util::error(ss.str());
   }
   std::vector<float> lat0(1, lat);
   std::vector<float> lon0(1, lon);
   std::vector<float> elev0(1, elev);
   std::vector<float> landFraction0(1, Util::MV);
   vec2 lats, lons;
   lats.push_back(lat0);
   lons.push_back(lon0);
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
   elevs.push_back(elev0);
   setElevs(elevs);
   mLandFractions.push_back(landFraction0);
   mNEns = 1;
   std::vector<double> times;
   mNEns = Util::MV;

   // Determine time and ensemble dimension if possible
   std::ifstream ifs(getFilename().c_str());
   if(ifs.good()) {
      while(ifs.good()) {
         char line[10000];
         ifs.getline(line, 10000, '\n');
         if(ifs.good() && line[0] != '#') {
            std::stringstream ss(line);
            double time;
            ss >> time;
            times.push_back(time);
            double value;
            mNEns = 0;
            while(ss >> value) {
               mNEns++;
            }
         }
      }
   }

   // Otherwise get the time or ensemble dimension from options
   iOptions.getValue("ens", mNEns);
   int numTimes;
   if(iOptions.getValue("time", numTimes)) {
      times.clear();
      // Empty file, probably used as output only
      for(int i = 0; i < numTimes; i++)
         times.push_back(i);
   }

   // Check that we got time and ensemble dimension
   if(!Util::isValid(numTimes)) {
      Util::error("Missing 'time' option for empty file '" + iFilename + "'");
   }
   if(!Util::isValid(mNEns)) {
      Util::error("Missing 'ens' option for empty file '" + iFilename + "'");
   }
   setTimes(times);
}

FilePoint::~FilePoint() {
}

FieldPtr FilePoint::getFieldCore(const Variable& iVariable, int iTime) const {
   std::ifstream ifs(getFilename().c_str());
   FieldPtr field = getEmptyField();

   int row = 0;
   while(ifs.good()) {
      char line[10000];
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         double time;
         ss >> time;

         if(time == iTime) {
            double value;
            int e = 0;
            while(ss >> value) {
               if(e >= mNEns) {
                  std::stringstream ss;
                  ss << "Row " << row << " in file '" << getFilename() << "' has too many members (expecting " << mNEns << ")" << std::endl;
                  Util::error(ss.str());
               }
               (*field)(0,0,e) = value;
               e++;
            }
            if(e != mNEns) {
               std::stringstream ss;
               ss << "Row " << row << " in file '" << getFilename() << "' has too many members (expecting " << mNEns << ")" << std::endl;
               Util::error(ss.str());
            }
         }
      }
      row++;
   }

   ifs.close();
   return field;
}

void FilePoint::writeCore(std::vector<Variable> iVariables, std::string iMessage) {
   std::ofstream ofs(getFilename().c_str());
   if(iVariables.size() == 0) {
      Util::warning("No variables to write");
      return;
   }
   for(int i = 0; i < getNumTime(); i++) {
      ofs.precision(0);
      FieldPtr field = getField(iVariables[0], i);
      if(field != NULL) {
         ofs << (long) getTimes()[i];
         ofs.precision(2);
         for(int e = 0; e < getNumEns(); e++) {
            ofs << std::fixed << " " << (*field)(0,0,e);
         }
         ofs << std::endl;
      }
   }
   ofs.close();
}

std::string FilePoint::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=point", "Point file for one location. Each row contains columns where the first column is the UNIX time and the second and onward columns are an ensemble of forecast values (each member has one column).") << std::endl;
   ss << Util::formatDescription("   lat=required", "Latitude (in degrees, north is positive)") << std::endl;
   ss << Util::formatDescription("   lon=required", "Longitude (in degrees, east is positive)") << std::endl;
   ss << Util::formatDescription("   elev=required", "Elevation (in meters)") << std::endl;
   ss << Util::formatDescription("   time=undef", "Number of times. Required if the file does not exist.") << std::endl;
   ss << Util::formatDescription("   ens=1", "Number of ensemble members.") << std::endl;
   return ss.str();
}
