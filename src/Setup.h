#ifndef METCAL_SETUP_H
#define METCAL_SETUP_H
#include <map>
#include <string>
#include <vector>
#include "Variable.h"
#include "Options.h"
class File;
class Calibrator;
class Downscaler;
class ParameterFile;
class Variable;

//! Represents the post-processing of one variable
struct VariableConfiguration {
   //! Which variable should be post-processed?
   Variable inputVariable;
   Variable outputVariable;
   //! Which downscaler should be used
   Downscaler* downscaler;
   //! Which calibrators should be use
   std::vector<Calibrator*> calibrators;
   Options inputVariableOptions;
   Options outputVariableOptions;
   std::vector<ParameterFile*> parameterFileCalibrators;
   ParameterFile* parameterFileDownscaler;
};

//! Represents what and how the post-processing should be done. Includes which input file to
//! post-process, which output file to place the results in, which variables to post-process,
//! and what post-processing methods to invoke on each variable.
//! Aborts if one of the input files does not exits/cannot be parsed.
class Setup {
   public:
      std::vector<File*> inputFiles;
      std::vector<File*> outputFiles;
      Options inputOptions;
      Options outputOptions;
      std::vector<VariableConfiguration> variableConfigurations;
      Setup(const std::vector<std::string>& argv);
      ~Setup();
      static std::string defaultDownscaler();
      std::map<std::string, Variable> variableAliases;
   private:
      // In some cases, it is not possible to open the same file first as readonly and then writeable
      // (for NetCDF). Therefore, use the same filehandle for both if the files are the same. Remember
      // to not free the memory of both files.
      std::map<std::string, File*> mFileMap; // filename, file handle
      bool hasFile(std::string iFilename) const;
};
#endif
