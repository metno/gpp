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

void writeUsage() {
   std::cout << "Post-processes gridded forecasts" << std::endl;
   std::cout << std::endl;
   std::cout << "usage:  gridpp input [options] output [options] [-v var [options] [-d downscaler [options] [-p parameters [options]]] [-c calibrator [options] [-p parameters [options]]]*]+" << std::endl;
   std::cout << "        gridpp [--version]" << std::endl;
   std::cout << "        gridpp [--help]" << std::endl;
   std::cout << std::endl;
   std::cout << "Arguments:" << std::endl;
   std::cout << "   input         Input file with NetCDF data." << std::endl;
   std::cout << "   output        Output file with NetCDF data. File must already exist. Must" << std::endl;
   std::cout << "                 contain lat/lon information." << std::endl;
   std::cout << "   -v var        One of the variables below." << std::endl;
   std::cout << "   -d downscaler One of the downscalers below." << std::endl;
   std::cout << "   -c calibrator One of the calibrators below." << std::endl;
   std::cout << "   -p parameters One of the parameter formats below." << std::endl;
   std::cout << "   options       Options of the form key=value" << std::endl;
   std::cout << "   --version     Print the program's version" << std::endl;
   std::cout << "   --help        Print usage information" << std::endl;
   std::cout << std::endl;
   std::cout << "Notes:" << std::endl;
   std::cout << "   - At least one variable must be specified." << std::endl;
   std::cout << "   - If multiple downscalers are specified for one variable, the last is used." << std::endl;
   std::cout << "   - If the same variable is specified multiple times, the first definition is used." << std::endl;
   std::cout << "   - Multiple identical calibrators are allowed for a single variable." << std::endl;
   std::cout << "   - Only one parameter format can be specified for given a downscaler or calibrator." << std::endl;
   std::cout << std::endl;
   std::cout << "Example:" << std::endl;
   std::cout << "   ./gridpp testing/files/10x10.nc testing/files/10x10_copy.nc -v T -d nearestNeighbour" << std::endl;
   std::cout << std::endl;
   std::cout << "Inputs/Outputs:" << std::endl;
   std::cout << "   I/O types are autodetected, but can be specified using:" << std::endl;
   std::cout << File::getDescriptions();
   std::cout << std::endl;
   std::cout << "Variables:" << std::endl;
   std::cout << Variable::getDescriptions();
   std::cout << std::endl;
   std::cout << "Variable options (and default values):" << std::endl;
   std::cout << Util::formatDescription("write=1", "Set to 0 to prevent the variable to be written to output") << std::endl;
   std::cout << std::endl;
   std::cout << "Downscalers with options (and default values):" << std::endl;
   std::cout << Downscaler::getDescriptions();
   std::cout << std::endl;
   std::cout << "Calibrators with options (and default values):" << std::endl;
   std::cout << Calibrator::getDescriptions();
   std::cout << std::endl;
   std::cout << "Parameter formats with Options (and default values)" << std::endl;
   std::cout << ParameterFile::getDescriptions();
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
         writeUsage();
         return 0;
      }

   }
   if(argc < 3) {
      writeUsage();
      return 0;
   }
   Util::setShowError(true);
   Util::setShowWarning(true);
   Util::setShowStatus(false);

   // Retrieve setup
   std::vector<std::string> args;
   for(int i = 5; i < argc; i++) {
      args.push_back(std::string(argv[i]));
   }
   int startDate, endDate;
   std::stringstream ss0(argv[1]);
   std::stringstream ss1(argv[2]);
   std::string inprefix(argv[3]);
   std::string outprefix(argv[4]);
   ss0 >> startDate;
   ss1 >> endDate;
   Setup setup(args);
   std::cout << "Input type:  " << setup.inputFile->name() << std::endl;
   std::cout << "Output type: " << setup.outputFile->name() << std::endl;

   for(int date = startDate; date <= endDate; date = Util::calcDate(date, 24)) {
      std::cout << "Date: " << date << std::endl;
      delete setup.inputFile;
      delete setup.outputFile;
      std::stringstream ssin;
      ssin << inprefix << date << ".nc";
      std::stringstream ssout;
      ssout << outprefix << date << ".nc";
      setup.inputFile = File::getScheme(ssin.str(), setup.inputOptions, true);
      setup.outputFile = File::getScheme(ssout.str(), setup.outputOptions, false);

      setup.outputFile->setTimes(setup.inputFile->getTimes());
      setup.outputFile->setReferenceTime(setup.inputFile->getReferenceTime());

      // Post-process file
      std::vector<Variable::Type> writeVariables;
      for(int v = 0; v < setup.variableConfigurations.size(); v++) {
         double s = Util::clock();
         VariableConfiguration varconf = setup.variableConfigurations[v];
         Variable::Type variable = varconf.variable;

         bool write = 1;
         varconf.variableOptions.getValue("write", write);
         if(write) {
            writeVariables.push_back(variable);
         }
         setup.outputFile->initNewVariable(variable);

         std::cout << "Processing " << Variable::getTypeName(variable) << std::endl;

         // Downscale
         std::cout << "   Downscaler " << varconf.downscaler->name() << std::endl;
         varconf.downscaler->downscale(*setup.inputFile, *setup.outputFile);
         for(int c = 0; c < varconf.calibrators.size(); c++) {
            // Calibrate
            std::cout << "   Calibrator " << varconf.calibrators[c]->name() << std::endl;
            varconf.calibrators[c]->calibrate(*setup.outputFile, varconf.parameterFileCalibrators[c]);
         }
         double e = Util::clock();
         std::cout << "   " << e-s << " seconds" << std::endl;
         std::cout << "Current mem usage input: " << setup.inputFile->getCacheSize() / 1e6<< std::endl;
         std::cout << "Current mem usage output: " << setup.outputFile->getCacheSize() / 1e6<< std::endl;
         // setup.inputFile->clear();
      }

      // Write to output
      double s = Util::clock();
      setup.outputFile->write(writeVariables);
      double e = Util::clock();
      std::cout << "Writing file: " << e-s << " seconds" << std::endl;
      std::cout << "Total time:   " << e-start << " seconds" << std::endl;
   }

   return 0;
}
