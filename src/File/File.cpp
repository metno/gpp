#include "File.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include "../Util.h"
#include "../Options.h"
Uuid File::mNextTag = 0;

File::File(std::string iFilename, const Options& iOptions) :
      mFilename(iFilename),
      mHasElevs(false),
      mReferenceTime(Util::MV) {
   createNewTag();

}

File* File::getScheme(std::string iFilename, const Options& iOptions, bool iReadOnly) {
   File* file;
   // Determine the filetype, either through user-specified option type=...
   // or by autodetecting.
   std::string type = "";
   if(!iOptions.getValue("type", type)) {
      // Autodetect type based on content
      if(FileNetcdf::isValid(iFilename, iOptions)) {
         type = "netcdf";
      }
   }

   // Instantiate the file
   if(type == "") {
      if(!Util::exists(iFilename)) {
         Util::warning("File '" + iFilename + " does not exist");
         return NULL;
      }
      else {
         Util::warning("Could not find suitable parser for '" + iFilename + "'");
         return NULL;
      }
   }
   else if(type == "netcdf") {
      return new FileNetcdf(iFilename, iOptions, iReadOnly);
   }
   else if(type == "point") {
      return new FilePoint(iFilename, iOptions);
   }
   else if(type == "text") {
      return new FileText(iFilename, iOptions);
   }
   else if(type == "norcomQnh") {
      return new FileNorcomQnh(iFilename, iOptions);
   }
   else {
      Util::error("Could not understand file type " + type);
   }
}


FieldPtr File::getField(std::string iVariable, int iTime) const {
   // Check internal variables first
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].name() == iVariable)
         return getField(mVariables[i], iTime);
   }

   // Check aliases
   std::map<std::string, Variable>::const_iterator it = mVariableAliases.find(iVariable);
   if(it != mVariableAliases.end()) {
      std::stringstream ss;
      ss << "Retrieving variable alias " << iVariable << " for timestep " << iTime;
      Util::info(ss.str());
      return getField(it->second, iTime);
   }

   // Create new variable
   return getField(Variable(iVariable), iTime);
}
FieldPtr File::getField(const Variable& iVariable, int iTime, bool iSkipRead) const {
   // Determine if values have been cached
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   bool needsReading = it == mFields.end();
   if(!needsReading) {
      // The variable has at least been partly read
      if(mFields[iVariable].size() <= iTime) {
         // This is an internal error. The variable has been partly read, but space has not
         // been allocated for all timesteps.
         std::stringstream ss;
         ss << "Attempted to access variable '" << iVariable.name() << "' for time " << iTime
            << " in file '" << getFilename() << "'";
         Util::error(ss.str());
      }

      needsReading = mFields[iVariable][iTime] == NULL;
   }
   else {
      // The variable has never been read or diagnosed. Allocate space for it.
      if(getNumTime() <= iTime) {
         std::stringstream ss;
         ss << "Attempted to access variable '" << iVariable.name() << "' for time " << iTime
            << " in file '" << getFilename() << "'";
         Util::error(ss.str());
      }
      mFields[iVariable].resize(getNumTime());
   }

   if(needsReading) {
      // Load non-derived variable from file
      if(!iSkipRead && hasVariableCore(iVariable)) {
         addField(getFieldCore(iVariable, iTime), iVariable, iTime);
      }
      else {
         std::string variableType = iVariable.name();
         Util::warning(variableType + " not available in '" + getFilename() + "'");
         for(int t = 0; t < getNumTime(); t++) {
            FieldPtr field = getEmptyField();
            addField(field, iVariable, t);
         }
      }
   }
   if(mFields[iVariable].size() <= iTime) {
      std::stringstream ss;
      ss << "Attempted to access variable '" << iVariable.name() << "' for time " << iTime
         << " in file '" << getFilename() << "'";
      Util::error(ss.str());
   }
   FieldPtr field = mFields[iVariable][iTime];
   if(!hasDefinedVariable(iVariable))
      mVariables.push_back(iVariable);
   return field;
}

File::~File() {
}

void File::write(std::vector<Variable> iVariables, std::string iMessage) {
   writeCore(iVariables, iMessage);
   // mCache.clear();
}


FieldPtr File::getEmptyField(float iFillValue) const {
   return getEmptyField(getNumY(), getNumX(), getNumEns(), iFillValue);
}
FieldPtr File::getEmptyField(int nY, int nX, int nEns, float iFillValue) const {
   FieldPtr field = FieldPtr(new Field(nY, nX, nEns, iFillValue));
   return field;
}

void File::addField(FieldPtr iField, const Variable& iVariable, int iTime) const {
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   if(it == mFields.end()) {
      mFields[iVariable].resize(getNumTime());
   }
   if(!hasDefinedVariable(iVariable))
      mVariables.push_back(iVariable);

   mFields[iVariable][iTime] = iField;
}

bool File::hasSameDimensions(const File& iOther) const {
   if(getNumY() == iOther.getNumY()
         && getNumX() == iOther.getNumX()
         && getNumEns() == iOther.getNumEns()
         && getNumTime() == iOther.getNumTime())
      return true;
   return false;
}

std::string File::getDimenionString() const {
   std::stringstream ss;
   ss << "[" << getNumTime() << " " << getNumEns() << " " << getNumY() << " " << getNumX()<< "]";
   return ss.str();
}

std::string File::getFilename() const {
   return mFilename;
}

void File::initNewVariable(const Variable& iVariable) {
   if(!hasVariable(iVariable)) {
      for(int t = 0; t < getNumTime(); t++) {
         addField(getEmptyField(), iVariable, t);
      }
   }
}
void File::setVariables(std::vector<Variable> iVariables) {
   mVariables = iVariables;
}
bool File::hasVariable(const Variable& iVariable) const {
   bool status = hasVariableCore(iVariable);
   if(status)
      return true;

   // Check if field has been initialized
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it = mFields.find(iVariable);
   return it != mFields.end();
}

void File::clear() {
   mFields.clear();
}

long File::getCacheSize() const {
   long size = 0;
   std::map<Variable, std::vector<FieldPtr> >::const_iterator it;
   for(it = mFields.begin(); it != mFields.end(); it++) {
      size += it->second.size() * getNumY()*getNumX()*getNumEns()*sizeof(float);
   }
   return size;
}

Uuid File::getUniqueTag() const {
   return mTag;
}
bool File::setLats(vec2 iLats) {
   bool uninitialized = mLats.size() == 0;
   if(!uninitialized && (iLats.size() != getNumY() || iLats[0].size() != getNumX()))
      return false;
   if(mLats != iLats)
      createNewTag();
   mLats = iLats;

   // Check that latitudes are valid
   for(int i = 0; i < mLats.size(); i++) {
      for(int j = 0; j < mLats[i].size(); j++) {
         if(mLats[i][j] < -90 || mLats[i][j] > 90) {
            std::stringstream ss;
            ss << "Invalid latitude in file '" << getFilename() << "': " << mLats[i][j] << std::endl;
            Util::error(ss.str());
         }
      }
   }
   return true;
}
bool File::setLons(vec2 iLons) {
   bool uninitialized = mLons.size() == 0;
   if(!uninitialized && (iLons.size() != getNumY() || iLons[0].size() != getNumX()))
      return false;
   if(mLons != iLons)
      createNewTag();
   mLons = iLons;

   // Check that longitudes are valid
   for(int i = 0; i < mLons.size(); i++) {
      for(int j = 0; j < mLons[i].size(); j++) {
         if(mLons[i][j] < -360 || mLons[i][j] > 360) {
            std::stringstream ss;
            ss << "Invalid longitude in file '" << getFilename() << "': " << mLons[i][j] << std::endl;
            Util::error(ss.str());
         }
      }
   }

   for(int i = 0; i < mLons.size(); i++) {
      for(int j = 0; j < mLons[i].size(); j++) {
         float lon = mLons[i][j];
         if(Util::isValid(lon)) {
            // Ensure lon is between -180 and 180
            int sign = lon / fabs(lon);
            lon = fabs(lon);
            lon = fmod(lon,360); // lon is between 0 and 360
            lon = sign * lon; // lon is between -360 and 360
            if(lon > 180)
               lon = lon - 360;
            else if(lon < -180)
               lon = lon + 360;
            mLons[i][j] = lon;
            assert(mLons[i][j] >= -180.0001 && mLons[i][j] <= 180.0001);
         }
      }
   }

   return true;
}
bool File::setElevs(vec2 iElevs) {
   if(iElevs.size() != getNumY() || iElevs[0].size() != getNumX())
      return false;
   mElevs = iElevs;
   mHasElevs = true;
   return true;
}
bool File::setLandFractions(vec2 iLandFractions) {
   if(iLandFractions.size() != getNumY() || iLandFractions[0].size() != getNumX())
      return false;
   mLandFractions = iLandFractions;
   return true;
}
vec2 File::getLats() const {
   return mLats;
}
vec2 File::getLons() const {
   return mLons;
}
vec2 File::getElevs() const {
   // Elevations not set, return a grid of missing values
   if(mElevs.size() == 0) {
      vec2 elevs;
      elevs.resize(getNumY());
      for(int i = 0; i < elevs.size(); i++)
         elevs[i].resize(getNumX(), Util::MV);
      return elevs;
   }
   else
      return mElevs;
}
vec2 File::getLandFractions() const {
   return mLandFractions;
}
int File::getNumY() const {
   return mLats.size();
}
int File::getNumX() const {
   return mLats[0].size();
}
int File::getNumEns() const {
   return mNEns;
}
int File::getNumTime() const {
   return mTimes.size();
}
void File::createNewTag() const {
   mTag = mNextTag; //boost::uuids::random_generator()();
   mNextTag++;
}
void File::setReferenceTime(double iTime) {
   mReferenceTime = iTime;
}
double File::getReferenceTime() const {
   return mReferenceTime;
}
void File::setTimes(std::vector<double> iTimes) {
   mTimes = iTimes;
}

bool File::setNumEns(int iNum) {
   std::stringstream ss;
   ss << "Setting number of ensemble members to " << iNum;
   // if(Util::isValid(mNEns)) {
   //    ss << ". Already set.";
   //    Util::info(ss.str());
   //    return false;
   // }
   mNEns = iNum;
   Util::info(ss.str());
   return true;
}
std::vector<double> File::getTimes() const {
   return mTimes;
}

bool File::getVariable(std::string iVariableName, Variable& iVariable) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i].name() == iVariableName) {
         iVariable = mVariables[i];
         return true;
      }
   }
   return false;
}

bool File::hasDefinedVariable(Variable iVariable) const {
   for(int i = 0; i < mVariables.size(); i++) {
      if(mVariables[i] == iVariable) {
         return true;
      }
   }
   return false;
}

bool File::hasElevs() const {
   return mHasElevs;
}

void File::addVariableAlias(std::string iAlias, Variable iVariable) {
   mVariableAliases[iAlias] =  iVariable;
}

std::string File::getDescriptions() {
   std::stringstream ss;
   ss << FileNetcdf::description();
   ss << FilePoint::description();
   ss << FileNorcomQnh::description();
   ss << FileText::description();
   return ss.str();
}
