#include "ParameterFile.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include <assert.h>
#include <set>
#include <fstream>

ParameterFile::ParameterFile(const Options& iOptions, bool iIsNew) :
      Scheme(iOptions),
      mFilename(""),
      mOptions(iOptions),
      mIsTimeDependent(false),
      mMaxTime(0),
      mIsNew(iIsNew) {
   iOptions.getValue("file", mFilename);
}

void ParameterFile::recomputeTree() const {
   if(isLocationDependent()) {
      vec2 lats, lons;
      LocationParameters::const_iterator it = mParameters.begin();
      for(it = mParameters.begin(); it != mParameters.end(); it++) {
         const Location loc = it->first;
         std::vector<float> lat(1, loc.lat());
         std::vector<float> lon(1, loc.lon());
         lats.push_back(lat);
         lons.push_back(lon);
         mLocations.push_back(loc);
      }
      mNearestNeighbourTree.build(lats, lons);
   }
}

ParameterFile* ParameterFile::getScheme(std::string iName, const Options& iOptions, bool iIsNew) {
   ParameterFile* p;
   if(iName == "text") {
      p = new ParameterFileText(iOptions, iIsNew);
   }
   else if(iName == "netcdf") {
      p = new ParameterFileNetcdf(iOptions, iIsNew);
   }
   else {
      Util::error("Parameter file type '" + iName + "' not recognized");
   }
   return p;
}

Parameters ParameterFile::getParameters(int iTime) const {
   if(iTime < 0) {
      std::stringstream ss;
      ss << "Could not load parameters for time " << iTime;
      Util::error(ss.str());
   }

   int time = iTime;
   if(!isTimeDependent())
      time = 0;

   if(time > mMaxTime) {
      std::stringstream ss;
      ss << "Could not load parameters for time " << time << " (max " << mMaxTime << ")";
      Util::error(ss.str());
   }

   if(isLocationDependent()) {
      Util::error("Cannot retrieve location-independent parameters for a location-dependent file");
   }
   // Find the right location to use
   LocationParameters::const_iterator it;
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
       it = mParameters.begin();
   }

   // Find the right time to use
   const std::vector<Parameters>& timeParameters = it->second;
   if(timeParameters.size() > time)
      return timeParameters[time];
   else {
      return Parameters();
   }
}
Parameters ParameterFile::getParameters(int iTime, const Location& iLocation, bool iAllowNearestNeighbour) const {
   if(iTime < 0) {
      std::stringstream ss;
      ss << "Could not load parameters for time " << iTime;
      Util::error(ss.str());
   }

   int time = iTime;
   if(!isTimeDependent())
      time = 0;

   if(time > mMaxTime) {
      std::stringstream ss;
      ss << "Could not load parameters for time " << time << " (max " << mMaxTime << ")";
      Util::error(ss.str());
   }

   if(mParameters.size() == 0)
      return Parameters();
   // Find the right location to use
   Location loc = iLocation;
   if(iAllowNearestNeighbour) {
      bool found = getNearestLocation(time, iLocation, loc);
      if(!found)
         return Parameters();
   }

   // Find the right time to use
   LocationParameters::const_iterator it = mParameters.find(loc);
   if(it == mParameters.end())
      return Parameters();

   const std::vector<Parameters>& timeParameters = it->second;
   if(timeParameters.size() > time)
      return timeParameters[time];
   else {
      return Parameters();
   }
}

bool ParameterFile::getNearestLocation(int iTime, const Location& iLocation, Location& iNearestLocation) const {
   if(mParameters.size() == 1) {
      // One set of parameters for all locations
      LocationParameters::const_iterator it = mParameters.begin();
      iNearestLocation = it->first;
   }
   else {
      // Try to see if we have an exact location
      LocationParameters::const_iterator it2 = mParameters.find(iLocation);
      if(it2 == mParameters.end()) {
         // If not, use the nearest neighbour
         int I, J;
         mNearestNeighbourTree.getNearestNeighbour(iLocation.lat(), iLocation.lon(), I, J);
         const Location& loc = mLocations[I];
         it2 = mParameters.find(loc);
      }
      bool hasAtThisTime = it2->second.size() > iTime && it2->second[iTime].size() != 0;
      if(hasAtThisTime) {
         iNearestLocation = it2->first;
         return true;
      }
      else {
         LocationParameters::const_iterator it = mParameters.find(iLocation);
         if(it == mParameters.end() || it->second.size() <= iTime || it->second[iTime].size() == 0) {
            LocationParameters::const_iterator it;
            float minDist = Util::MV;
            for(it = mParameters.begin(); it != mParameters.end(); it++) {
               const Location& currLoc = it->first;
               float dist = iLocation.getDistance(currLoc);
               if(!Util::isValid(minDist) || (Util::isValid(dist) && dist < minDist)) {
                  // Check that this location actually has parameters available for this time
                  bool hasAtThisTime = it->second.size() > iTime && it->second[iTime].size() != 0;
                  if(hasAtThisTime) {
                     minDist = dist;
                     iNearestLocation = it->first;
                  }
               }
            }
            if(!Util::isValid(minDist)) {
               return false;
            }
            else {
            }
         }
         else {
            iNearestLocation = it->first;
         }
      }
   }
   return true;
}

void ParameterFile::setParameters(Parameters iParameters, int iTime, const Location& iLocation) {
   mMaxTime = std::max(mMaxTime, iTime);
   LocationParameters::const_iterator it = mParameters.find(iLocation);
   if(mParameters[iLocation].size() <= iTime) {
      mParameters[iLocation].resize(mMaxTime+1);
   }
   mParameters[iLocation][iTime] = iParameters;
   mIsTimeDependent = mIsTimeDependent || iTime > 0;
}
void ParameterFile::setParameters(Parameters iParameters, int iTime) {
   setParameters(iParameters, iTime, Location(Util::MV,Util::MV,Util::MV));
}

std::string ParameterFile::getFilename() const {
   return mFilename;
}

std::vector<Location> ParameterFile::getLocations() const {
   return mLocations;
}

std::vector<int> ParameterFile::getTimes() const {
   std::set<int> times;
   LocationParameters::const_iterator it;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      for(int i = 0; i < it->second.size(); i++) {
         times.insert(i);
      }
   }
   return std::vector<int>(times.begin(), times.end());
}

bool ParameterFile::isLocationDependent() const {
   bool locationDependent = mParameters.size() > 1;
   if(!locationDependent) {
      Location location = mParameters.begin()->first;
      if(Util::isValid(location.lat()) && Util::isValid(location.lon()))
         locationDependent = true;
   }
   return locationDependent;
}
// TODO
bool ParameterFile::isTimeDependent() const {
   return mIsTimeDependent;
}

int ParameterFile::getNumParameters() const {
   int size = Util::MV;

   LocationParameters::const_iterator itLoc;
   for(itLoc = mParameters.begin(); itLoc != mParameters.end(); itLoc++) {
      const std::vector<Parameters>& parvec = itLoc->second;
      for(int i = 0; i < parvec.size(); i++) {
         int currSize = parvec[i].size();
         if(currSize != 0) {
            if(Util::isValid(size) && currSize != size) {
               std::cout << "Missmatch in parameter size" << size << " " << currSize << std::endl;
               return Util::MV;
            }
            size = currSize;
         }
         else
            return 0;
      }
   }
   return size;
}

void ParameterFile::setFilename(std::string iFilename) {
   mFilename = iFilename;
}

std::string ParameterFile::getDescriptions(bool full) {
   std::stringstream ss;
   ss << ParameterFileText::description(full) << std::endl;
   ss << ParameterFileNetcdf::description(full) << std::endl;
   return ss.str();
}

long ParameterFile::getCacheSize() const {
   std::map<Location, std::vector<Parameters>, Location::CmpIgnoreElevation >::const_iterator it;
   long total = 0;
   for(it = mParameters.begin(); it != mParameters.end(); it++) {
      for(int i = 0; i < it->second.size(); i++) {
         total++;
      }
   }
   return total;
}

void ParameterFile::initializeEmpty(const std::vector<Location>& iLocations, int iNumTimes, int iNumParameters) {
   std::vector<float> params(iNumParameters, Util::MV);
   Parameters parameters(params);
   std::vector<Parameters> pars(iNumTimes, params);
   for(int i = 0; i < iLocations.size(); i++) {
      mParameters[iLocations[i]] = pars;
   }
}

void ParameterFile::setIsTimeDependent(bool iFlag) {
   mIsTimeDependent = iFlag;
}

void ParameterFile::setMaxTime(int iMaxTime) {
   mMaxTime = iMaxTime;
}

int ParameterFile::getMaxTime() const {
   return mMaxTime;
}

Options ParameterFile::getOptions() const {
   return mOptions;
}
