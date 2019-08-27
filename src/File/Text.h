#ifndef FILE_TEXT_H
#define FILE_TEXT_H
#include <vector>
#include <map>
#include "File.h"
#include "../Variable.h"
#include "../Options.h"

class Location;

//! Represents a point-based text file. First column is the time (in seconds since 1970)
//! and the second column is the forecast value.
class FileText : public File {
   public:
      FileText(std::string iFilename, const Options& iOptions);
      static std::string description();
      std::string name() const {return "text";};
   protected:
      FieldPtr getFieldCore(const Variable& iVariable, int iTime) const;
      void writeCore(std::vector<Variable> iVariables, std::string iMessage="");
      bool hasVariableCore(const Variable& iVariable) const {return true;};
      std::vector<FieldPtr> mLocalFields;
};
#endif
