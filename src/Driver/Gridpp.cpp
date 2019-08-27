#include <iostream>
#include <string>
#include <string.h>
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Calibrator.h"
#include "../Downscaler/Downscaler.h"
#include "../Util.h"
#include "../Options.h"
#include "../Setup.h"

void writeUsage(bool full) {
   std::cout << "Post-processes gridded forecasts. For more information see https://github.com/metno/gridpp." << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp inputs [options] outputs [options] [-v var [options] [-d downscaler [options] [-p parameters [options]]] [-c calibrator [options] [-p parameters [options]]]*]+ [--debug <level>]" << std::endl;
   std::cout << "        gridpp [--version]" << std::endl;
   std::cout << "        gridpp [--help]" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << "   inputs        One or more input files. Use comma to separate names. When" << std::endl;
   std::cout << "                 string is surrounded by double quotes, wildcards can be used (e.g. *)." << std::endl;
   std::cout << "   outputs       One or more output files. File must already exist. Must" << std::endl;
   std::cout << "                 contain lat/lon information. If multiple input files are" << std::endl;
   std::cout << "                 used, then the same number of outputs must be used." << std::endl;
   std::cout << "   -v var        One of the variables below." << std::endl;
   std::cout << "   -d downscaler One of the downscalers below." << std::endl;
   std::cout << "   -c calibrator One of the calibrators below." << std::endl;
   std::cout << "   -p parameters One of the parameter formats below." << std::endl;
   std::cout << "   options       Options of the form key=value" << std::endl;
   std::cout << "   --version     Print the program's version" << std::endl;
   std::cout << "   --debug lvl   Set debug level: quiet, error, warn (default), info" << std::endl;
   std::cout << "   --help        Print usage information including all options" << std::endl;
   std::cout << std::endl;
   std::cout << "Inputs/Outputs:" << std::endl;
   std::cout << "   I/O types are autodetected, but can be specified using:" << std::endl;
   std::cout << File::getDescriptions();
   std::cout << std::endl;
   std::cout << "Variable options (and default values):" << std::endl;
   std::cout << Util::formatDescription("write=1", "Set to 0 to prevent the variable to be written to output") << std::endl;
   std::cout << Util::formatDescription("units=undef", "Write this to the units attribute in the output.") << std::endl;
   std::cout << Util::formatDescription("standardName=1", "Write this to the standard_name attribute in the output.") << std::endl;
   std::cout << std::endl;
   if(full)
      std::cout << "Downscalers with options (and default values):" << std::endl;
   else
      std::cout << "Downscalers:" << std::endl;
   std::cout << Downscaler::getDescriptions(full);
   std::cout << std::endl;
   if(full)
      std::cout << "Calibrators with options (and default values):" << std::endl;
   else
      std::cout << "Calibrators:" << std::endl;
   std::cout << Calibrator::getDescriptions(full);
   std::cout << std::endl;
   if(full)
      std::cout << "Parameter formats with Options (and default values):" << std::endl;
   else
      std::cout << "Parameter formats:" << std::endl;
   std::cout << ParameterFile::getDescriptions(full);
}

int main(int argc, const char *argv[]) {
   double start = Util::clock();

   // Parse command line attributes
   for(int i = 1; i < argc; i++) {
      if(strcmp(argv[i], "--version") == 0) {
         std::cout << "gridpp version " << Util::gridppVersion() << std::endl;
         return 0;
      }
      else if(strcmp(argv[i], "--help") == 0) {
         writeUsage(true);
         return 0;
      }

   }
   if(argc < 3) {
      writeUsage(false);
      return 0;
   }
   // Retrieve setup
   std::vector<std::string> args;
   std::string debugMode = "warn";
   Util::setShowError(true);
   for(int i = 1; i < argc; i++) {
      if(std::string(argv[i]) == "--debug") {
         i++;
         if(argc <= i) {
            Util::error("Missing debug level");
         }
         debugMode = std::string(argv[i]);
      }
      else {
         args.push_back(std::string(argv[i]));
      }
   }

   // Set logging levels
   Util::setShowStatus(true);
   Util::setShowWarning(true);
   Util::setShowInfo(true);
   if(debugMode == "quiet") {
      Util::setShowError(false);
      Util::setShowStatus(false);
      Util::setShowWarning(false);
      Util::setShowInfo(false);
   }
   else if(debugMode == "error") {
      Util::setShowStatus(false);
      Util::setShowWarning(false);
      Util::setShowInfo(false);
   }
   else if(debugMode == "warn") {
      Util::setShowInfo(false);
   }
   else if(debugMode == "info") {
   }

   Setup setup(args);
   for(int f = 0; f < setup.inputFiles.size(); f++) {
      Util::info("Input type:  " + setup.inputFiles[f]->name());
      Util::info("Output type: " + setup.outputFiles[f]->name());
      Util::info( "Input file '" + setup.inputFiles[f]->getFilename() + "' has dimensions " + setup.inputFiles[f]->getDimenionString());
      Util::info( "Output file '" + setup.outputFiles[f]->getFilename() + "' has dimensions " + setup.outputFiles[f]->getDimenionString());

      setup.outputFiles[f]->setTimes(setup.inputFiles[f]->getTimes());
      setup.outputFiles[f]->setNumEns(setup.inputFiles[f]->getNumEns());
      setup.outputFiles[f]->setReferenceTime(setup.inputFiles[f]->getReferenceTime());

      // Post-process file
      std::vector<Variable> writeVariables;
      for(int v = 0; v < setup.variableConfigurations.size(); v++) {
         double s = Util::clock();
         VariableConfiguration varconf = setup.variableConfigurations[v];
         Variable inputVariable = varconf.inputVariable;
         Variable outputVariable = varconf.outputVariable;

         bool write = 1;
         varconf.outputVariableOptions.getValue("write", write);
         if(write) {
            writeVariables.push_back(outputVariable);
         }
         setup.outputFiles[f]->initNewVariable(outputVariable);

         Util::status("Processing " + outputVariable.name());

         // Downscale
         Util::status("   Downscaler " +  varconf.downscaler->name() + ": ", false);
         double ss = Util::clock();
         varconf.downscaler->downscale(*setup.inputFiles[f], *setup.outputFiles[f]);
         double ee = Util::clock();
         std::stringstream ss0;
         ss0 << ee-ss << " seconds";
         Util::status(ss0.str());

         // Calibrate
         for(int c = 0; c < varconf.calibrators.size(); c++) {
            double s = Util::clock();
            Util::status("   Calibrator " + varconf.calibrators[c]->name() + ": ", false);
            varconf.calibrators[c]->calibrate(*setup.outputFiles[f], varconf.parameterFileCalibrators[c]);
            double e = Util::clock();
            std::stringstream ss;
            ss << e-s << " seconds";
            Util::status(ss.str());
         }
         double e = Util::clock();
         std::stringstream ss1;
         ss1 << "   Total: " << e-s << " seconds";
         Util::status(ss1.str());

         std::stringstream ss2;
         ss2 << "Mem usage input: " << setup.inputFiles[f]->getCacheSize() / 1e6;
         Util::status(ss2.str());

         std::stringstream ss3;
         ss3 << "Mem usage output: " << setup.outputFiles[f]->getCacheSize() / 1e6;
         Util::status(ss3.str());
         // setup.inputFile->clear();
      }

      // Write to output
      double s = Util::clock();
      std::stringstream ssMessage;
      for(int i = 1; i < argc; i++) {
         if(i > 1)
            ssMessage << " ";
         ssMessage << argv[i];
      }

      setup.outputFiles[f]->write(writeVariables, ssMessage.str());
      double e = Util::clock();
      std::stringstream ss1;
      ss1 << "Writing file: " << e-s << " seconds";
      Util::status(ss1.str());

      std::stringstream ss2;
      ss2 << "Total time:   " << e-start << " seconds";
      Util::status(ss2.str());
      setup.inputFiles[f]->clear();
      setup.outputFiles[f]->clear();
   }
   return 0;
}
