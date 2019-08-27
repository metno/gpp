#ifndef FILE_NORCOM_QNH_H
#define FILE_NORCOM_QNH_H
#include <vector>
#include <map>
#include "File.h"
#include "../Variable.h"
#include "../Options.h"

//! Represents a point-based text file suitable for sending minimum QNH values to Norcom. 
//! The output file looks as follows:
//! FBNO52 ENNC 290545                  // day,hour,minute
//! VALID 290600 - 290900 UTC.
//!  ST MIN QNH FANARAAKEN  : 0987 HPA
//!  ST MIN QNH GAUSTATOPPEN: 0989 HPA
class FileNorcomQnh : public File {
   public:
      FileNorcomQnh(std::string iFilename, const Options& iOptions);
      ~FileNorcomQnh();
      static std::string description();
      std::string name() const {return "norcom";};
   protected:
      FieldPtr getFieldCore(const Variable& iVariable, int iTime) const;
      void writeCore(std::vector<Variable> iVariables, std::string iMessage="");
      bool hasVariableCore(const Variable& iVariable) const {return true;};
      int mStartTime;
      int mEndTime;
      std::vector<std::string> mNames;
      // Create a formated string with format DDHHMM
      std::string getNorcomTimeStamp(time_t iUnixTime) const;
};
#endif
