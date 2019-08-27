#include "Text.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"
#include "../Location.h"

FileText::FileText(std::string iFilename, const Options& iOptions) :
      File(iFilename, iOptions) {

   std::ifstream ifs(getFilename().c_str(), std::ifstream::in);
   if(!ifs.good()) {
      return;
   }
   std::set<int> timesSet;
   std::set<Location> locationsSet;
   std::map<int, std::map<Location, std::vector<float> > > values;
   mNEns = Util::MV;
   while(ifs.good()) {
      char line[10000];
      ifs.getline(line, 10000, '\n');
      if(ifs.good() && line[0] != '#') {
         std::stringstream ss(line);
         // Loop over each value
         std::vector<float> currValues;
         int time;
         bool status = ss >> time;
         if(!status) {
            Util::error("Could not read time from file '" + iFilename + "'");
         }
         timesSet.insert(time);

         Location location(Util::MV, Util::MV, Util::MV);
         float lat;
         status = ss >> lat;
         if(!status) {
            Util::error("Could not read lat from file '" + iFilename + "'");
         }

         float lon;
         status = ss >> lon;
         if(!status) {
            Util::error("Could not read lon from file '" + iFilename + "'");
         }

         float elev;
         status = ss >> elev;
         if(!status) {
            Util::error("Could not read elev from file '" + iFilename + "'");
         }
         location = Location(lat, lon, elev);

         while(ss.good()) {
            float value;
            bool status  = ss >> value;
            if(!status) {
               Util::error("Could not read value from file '" + iFilename + "'");
            }
            currValues.push_back(value);
         }
         if(mNEns == Util::MV)
            mNEns = currValues.size();
         else if(currValues.size() != mNEns) {
            std::stringstream ss;
            ss << "File '" + getFilename() + "' is corrupt, because it does not have the same"
               << " number of columns on each line" << std::endl;
            Util::error(ss.str());
         }
         values[time][location] = currValues;
         locationsSet.insert(location);
      }
   }
   ifs.close();

   std::vector<double> times(timesSet.begin(), timesSet.end());
   std::vector<Location> locations(locationsSet.begin(), locationsSet.end());
   std::sort(times.begin(), times.end());

   vec2 lats, lons;
   lats.resize(locations.size());
   lons.resize(locations.size());
   vec2 elevs;
   elevs.resize(locations.size());
   for(int i = 0; i < locations.size(); i++) {
      lats[i].resize(1);
      lons[i].resize(1);
      elevs[i].resize(1);
      lats[i][0] = locations[i].lat();
      lons[i][0] = locations[i].lon();
      elevs[i][0] = locations[i].elev();
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
   setElevs(elevs);

   setTimes(times);

   // Put data into mLocalFields
   mLocalFields.resize(times.size());
   for(int t = 0; t < times.size(); t++) {
      std::map<int, std::map<Location, std::vector<float> > >::const_iterator it = values.find(times[t]);
      FieldPtr field = getEmptyField();
      for(int i = 0; i < locations.size(); i++) {
         std::map<Location, std::vector<float> >::const_iterator it2 = it->second.find(locations[i]);
         if(it2 != it->second.end()) {
            for(int e = 0; e < mNEns; e++) {
               (*field)(i,0,e) = it2->second[e];
            }
         }

      }
      mLocalFields[t] = field;
   }
}

FieldPtr FileText::getFieldCore(const Variable& iVariable, int iTime) const {
   return mLocalFields[iTime];
}

void FileText::writeCore(std::vector<Variable> iVariables, std::string iMessage) {
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

std::string FileText::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=text", "Text file format where each row has the columns:") << std::endl;
   ss << Util::formatDescription("", "time lat lon elev ens0 ens1 ... ensN.") << std::endl;
   return ss.str();
}
