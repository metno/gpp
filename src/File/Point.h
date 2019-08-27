#ifndef FILE_POINT_H
#define FILE_POINT_H
#include <vector>
#include <map>
#include "File.h"
#include "../Variable.h"
#include "../Options.h"

//! Represents a point-based text file. First column is the time (in seconds since 1970)
//! and the second column is the forecast value.
class FilePoint : public File {
   public:
      FilePoint(std::string iFilename, const Options& iOptions);
      ~FilePoint();
      static std::string description();
      std::string name() const {return "point";};
   protected:
      FieldPtr getFieldCore(const Variable& iVariable, int iTime) const;
      void writeCore(std::vector<Variable> iVariables, std::string iMessage="");
      bool hasVariableCore(const Variable& iVariable) const {return true;};
};
#endif
