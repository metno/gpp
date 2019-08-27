#ifndef VARIABLE_H
#define VARIABLE_H
#include <iostream>
#include <vector>
#include "Options.h"

//! Represents a meteorological variable and its metadata
class Variable {
   public:

      Variable();
      Variable(std::string iName, float iMin, float iMax, std::string iUnits, std::string iStandardName, int iLevel=Util::MV);
      Variable(std::string iName, std::string iUnits="", std::string iStandardName="", int iLevel=Util::MV);
      Variable(const Options& iOptions);

      // Name of variable
      std::string name() const;
      void name(std::string);

      //! Get the minimum possible attainable value for this variable
      float min() const;
      void min(float iValue);

      //! Get the maximum possible attainable value for this variable
      float max() const;
      void max(float iValue);

      //! Get the level in the heigh dimension
      float level() const;
      void level(int iValue);

      //! Returns the units of the variable
      std::string units() const;
      void units(std::string iValue);

      //! Returns the NetcdfCF standard name
      std::string standardName() const;
      void standardName(std::string iValue);

      void add(const Options& iOptions);

      bool operator<(const Variable &right) const;
      bool operator==(const Variable &right) const;

   private:
      float mMin;
      float mMax;
      std::string mUnits;
      std::string mStandardName;
      std::string mName;
      int mLevel;
};
#endif
