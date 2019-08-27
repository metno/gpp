#ifndef FILE_FAKE_H
#define FILE_FAKE_H
#include "File.h"

//! A file type giving fake data, useful for testing
class FileFake : public File {
   public:
      FileFake(const Options& iOptions);
      static std::string description();
      std::string name() const {return "fake";};
   protected:
      void writeCore(std::vector<Variable> iVariables, std::string iMessage="");
      FieldPtr getFieldCore(const Variable& iVariable, int iTime) const;
      bool hasVariableCore(const Variable& iVariable) const;
};
#endif
