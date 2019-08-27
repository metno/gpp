#include "Variable.h"
#include <sstream>
#include "Util.h"

Variable::Variable() {
}
Variable::Variable(std::string iName, float iMin, float iMax, std::string iUnits, std::string iStandardName, int iLevel):
      mName(iName),
      mMin(iMin),
      mMax(iMax),
      mUnits(iUnits),
      mStandardName(iStandardName),
      mLevel(iLevel) {
}
Variable::Variable(std::string iName, std::string iUnits, std::string iStandardName, int iLevel):
      mName(iName),
      mMin(Util::MV),
      mMax(Util::MV),
      mUnits(iUnits),
      mStandardName(iStandardName),
      mLevel(iLevel) {
}
Variable::Variable(const Options& iOptions):
      mName(""),
      mMin(Util::MV),
      mMax(Util::MV),
      mLevel(Util::MV),
      mStandardName(""),
      mUnits("") {
   iOptions.getRequiredValue("name", mName);
   iOptions.getValue("standardName", mStandardName);
   iOptions.getValue("units", mUnits);
   iOptions.getValue("min", mMin);
   iOptions.getValue("max", mMax);
   iOptions.getValue("level", mLevel);
   iOptions.check();
}

float Variable::min() const {
   return mMin;
}
void Variable::min(float iValue) {
   mMin = iValue;
}

float Variable::max() const {
   return mMax;
}
void Variable::max(float iValue) {
   mMax = iValue;
}

float Variable::level() const {
   return mLevel;
}
void Variable::level(int iValue) {
   mLevel = iValue;
}

std::string Variable::name() const {
   return mName;
}
void Variable::name(std::string iValue) {
   mName = iValue;
}

std::string Variable::units() const {
   return mUnits;
}
void Variable::units(std::string iValue) {
   mUnits = iValue;
}

std::string Variable::standardName() const {
   return mStandardName;
}
void Variable::standardName(std::string iValue) {
   mStandardName = iValue;
}

void Variable::add(const Options& iOptions) {
   std::string value;
   if(iOptions.getValue("units", value)) {
      mUnits = value;
   }
   if(iOptions.getValue("standardName", value)) {
      mStandardName = value;
   }
   float min;
   if(iOptions.getValue("min", min)) {
      mMin = min;
   }
   float max;
   if(iOptions.getValue("min", max)) {
      mMax = max;
   }
}

bool Variable::operator<(const Variable &right) const {
   return mName < right.name();
}

bool Variable::operator==(const Variable &right) const {
   return mName == right.name();
}
