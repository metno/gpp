#include "NorcomQnh.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include "../Util.h"
#include <iomanip>

FileNorcomQnh::FileNorcomQnh(std::string iFilename, const Options& iOptions) :
      File(iFilename, iOptions) {
   vec2 lats, lons;
   lats.resize(1);
   lons.resize(1);
   vec2 elevs;
   elevs.resize(1);
   mLandFractions.resize(1);
   mNEns = 1;
   if(!iOptions.getValues("lats", lats[0])) {
      Util::error("Missing 'lats' option for '" + iFilename + "'");
   }
   if(!iOptions.getValues("lons", lons[0])) {
      Util::error("Missing 'lons' option for '" + iFilename + "'");
   }
   if(!iOptions.getValues("elevs", elevs[0])) {
      Util::error("Missing 'elevs' option for '" + iFilename + "'");
   }
   mLandFractions[0].resize(elevs[0].size(), Util::MV);
   if(!iOptions.getValues("names", mNames)) {
      Util::error("Missing 'names' option for '" + iFilename + "'");
   }
   int numTimes = Util::MV;
   if(!iOptions.getValue("numTimes", numTimes)) {
      Util::error("Missing 'numTimes' option for '" + iFilename + "'");
   }
   if(lats[0].size() != lons[0].size() || lats[0].size() != elevs[0].size() || lats[0].size() != mNames.size()) {
      Util::error("FileNorcomQnh: 'lats', 'lons', 'elevs', 'names' must be the same size");
   }
   for(int i = 0; i < lats[0].size(); i++) {
      float lat = lats[0][i];
      if(lat < -90 || lat > 90) {
         std::stringstream ss;
         ss << "Invalid latitude: " << lat;
         Util::error(ss.str());
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
   setElevs(elevs);

   std::vector<double> times;
   for(int i = 0; i < numTimes; i++)
      times.push_back(i);

   // Determine the times for this filetype.
   if(!iOptions.getValue("startTime", mStartTime)) {
      Util::error("Missing 'startTime' option for '" + iFilename + "'");
   }
   if(!iOptions.getValue("endTime", mEndTime)) {
      Util::error("Missing 'endTime' option for '" + iFilename + "'");
   }
   if(mStartTime > mEndTime) {
      Util::error("FileNorcomQnh: 'startTime' must be <= 'endTime'");
   }
   setTimes(times);
}

FileNorcomQnh::~FileNorcomQnh() {
}

FieldPtr FileNorcomQnh::getFieldCore(const Variable& iVariable, int iTime) const {
   FieldPtr field = getEmptyField();
   return field;
}

void FileNorcomQnh::writeCore(std::vector<Variable> iVariables, std::string iMessage) {
   std::ofstream ofs(getFilename().c_str());
   if(iVariables.size() == 0) {
      Util::warning("No variables to write");
      return;
   }
   Variable variable = iVariables[0];
   if(iVariables.size() > 1) {
      std::stringstream ss;
      ss <<"Output NorcomQnh can only write one variables, several given. Will write variable ";
      ss << variable.name();
      Util::warning(ss.str());
   }

   // Find the length of the longest station name
   int maxNameSize = 0;
   for(int j = 0; j < mNames.size(); j++) {
      if(mNames[j].size() > maxNameSize)
         maxNameSize = mNames[j].size();
   }

   // Write header
   time_t currUnixTime = Util::getCurrentUnixTime();
   std::string currTimeStamp = getNorcomTimeStamp(currUnixTime);
   std::vector<double> times = getTimes();
   std::string startTime = getNorcomTimeStamp(times[mStartTime]);
   std::string endTime   = getNorcomTimeStamp(times[mEndTime]);
   ofs << "FBNO52 ENNC " << currTimeStamp << "\r\r" << std::endl;
   ofs << "VALID " << startTime << " - " << endTime << " UTC." << "\r" << std::endl;

   ofs.precision(0);
   // Write one line for each station
   vec2 lats = getLats();
   for(int j = 0; j < lats[0].size(); j++) {
      std::string locationName = mNames[j];
      ofs << "EST MIN QNH ";
      ofs << std::setfill(' ') << std::setw(maxNameSize) << std::left << locationName << ": ";

      // Find minimum
      int valuePa = Util::MV;
      for(int t = mStartTime; t <= mEndTime; t++) {
         FieldPtr field = getField(variable, t);
         int currValue = (int) (*field)(0,j,0);
         if(!Util::isValid(valuePa)|| (Util::isValid(currValue) && currValue < valuePa))
            valuePa = currValue;
      }
      if(!Util::isValid(valuePa)) {
         Util::error("Invalid value when writing QNH to NorcomQnh");
      }
      int valueHpa = valuePa / 100;
      ofs << std::setfill('0') << std::setw(4) << std::right << valueHpa << " HPA" << "\r" << std::endl;
   }
   ofs.close();
}

std::string FileNorcomQnh::getNorcomTimeStamp(time_t iUnixTime) const {
   int date = Util::getDate(iUnixTime);
   int time = Util::getTime(iUnixTime);
   std::stringstream ss;
   int day = date % 100;
   int HHMM = time / 100;
   ss << std::setfill('0') << std::right << std::setw(2) << day;
   ss << std::setfill('0') << std::right << std::setw(4) << HHMM;

   return ss.str();
}

std::string FileNorcomQnh::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=norcomQnh", "Output format for sending minimum QNH values to Norcom. Finds the minimum QNH values on the interval [startTime,endTime]. QNH must either exist in the input file, or created using a calibrator (such as -c qnh)") << std::endl;
   ss << Util::formatDescription("   lats=required", "Comma-separated list of latitudes (lat1,lat2,lat3,...). Values in degrees, north is positive.") << std::endl;
   ss << Util::formatDescription("   lons=required", "Longitudes (in degrees, east is positive)") << std::endl;
   ss << Util::formatDescription("   elevs=required", "Elevations (in meters)") << std::endl;
   ss << Util::formatDescription("   names=required", "Comma-separated list of station names.") << std::endl;
   ss << Util::formatDescription("   numTimes=undef", "Number of times. Set this equal to the number of times in the input file.") << std::endl;
   ss << Util::formatDescription("   startTime=undef", "First time index to find minimum over.") << std::endl;
   ss << Util::formatDescription("   endTime=undef", "Last time index to find minimum over.") << std::endl;
   return ss.str();
}
