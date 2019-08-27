#ifndef PARAMETER_FILE_TEXT_H
#define PARAMETER_FILE_TEXT_H
#include <iostream>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"

//! Represents a collection of parameters, one set for each forecast time
//! Parameters are read from a specified text file with the following format:
//! time lat lon elev param1 param2 param3
//! 0 60 10 3.1 3.4 2.1 5.2 12 41
//! 1 60 10 3.1 3.4 2.1 5.2 12 41
//! 2 60 10 3.1 3.4 2.1 5.2 12 41
//! Each line represents one forecast time. The first column is the forecast timestep (an index
//! not the number of hours), starting at 0. The remaning columns are parameters that can be used
//! in post-processing methods. The number of columns most be constant. If a file has only one
//! line, then the parameters are used for all forecast hours.
class ParameterFileText : public ParameterFile {
   public:
      //! Read parameters from this file
      ParameterFileText(const Options& iOptions, bool iIsNew=false);

      bool isFixedSize() const;

      std::vector<int> getTimes() const;
      static bool isValid(std::string iFilename);
      bool isReadable() const;
      bool isLocationDependent() const;

      //! Write parameter file to disk
      //! @param iFilename Write to this filename. If empty, write to the file that was read from.
      void write(const std::string& iFilename) const;
      void write() const;
      static std::string description(bool full=true);
      std::string name() const {return "text";};
   private:
      std::vector<int> mTimes;
};
#endif
