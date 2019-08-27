#include "Setup.h"
#include "File/File.h"
#include "Calibrator/Calibrator.h"
#include "Downscaler/Downscaler.h"

Setup::Setup(const std::vector<std::string>& argv) {

   std::string inputFilename = "";
   std::string outputFilename = "";

   // Process input/output filenames and options
   int index = 0;
   while(index < argv.size()) {
      std::string arg = argv[index];
      if(inputFilename == "") {
         inputFilename = arg;
      }
      else if(outputFilename == "") {
         if(Util::hasChar(arg, '=')) {
            inputOptions.addOptions(arg);
         }
         else if(arg[0] == '-') {
            // No output file specified, use the same input and output file
            outputFilename = inputFilename;
            outputOptions = inputOptions;
            break;
         }
         else {
            outputFilename = arg;
         }
      }
      else {
         if(Util::hasChar(arg, '=')) {
            outputOptions.addOptions(arg);
         }
         else {
            break;
         }
      }
      index++;
   }
   std::vector<std::string> inputFilenames = Util::glob(inputFilename);
   std::vector<std::string> outputFilenames = Util::glob(outputFilename);
   if(inputFilenames.size() != outputFilenames.size()) {
      std::stringstream ss;
      ss << "Unequal number of input (" << inputFilenames.size() << ") and output ("
         << outputFilenames.size() << ")";
      Util::error(ss.str());
   }
   if(inputFilenames.size() == 0) {
      Util::error("No valid input files");
   }
   if(outputFilenames.size() == 0) {
      Util::error("No valid output files");
   }

   for(int f = 0; f < outputFilenames.size(); f++) {

      File* outputFile;
      if(!hasFile(outputFilenames[f])) {
         outputFile  = File::getScheme(outputFilenames[f], outputOptions, false);
         if(outputFile == NULL) {
            Util::error("File '" + outputFilenames[f] + " is invalid");
         }
         mFileMap[outputFilenames[f]] = outputFile;
      }
      else {
         outputFile = mFileMap[outputFilenames[f]];
      }
      outputFiles.push_back(outputFile);

      File* inputFile;
      if(!hasFile(inputFilenames[f])) {
         inputFile  = File::getScheme(inputFilenames[f], inputOptions, true);
         if(inputFile == NULL) {
            Util::error("File '" + inputFilenames[f] + " is invalid");
         }
         mFileMap[inputFilenames[f]] = inputFile;
      }
      else {
         if(inputOptions != outputOptions) {
            Util::error("Input and output options cannot be different if they apply to the same filename");
         }
         inputFile = mFileMap[inputFilenames[f]];
      }
      inputFiles.push_back(inputFile);
   }

   // Implement a finite state machine
   enum State {START = 0, VAR = 1, VAROPT = 2, NEWVAR = 3, VARI = 4, VARIOPT = 5, VARALIAS = 6, VARALIASOPT = 7, NEWVARALIAS = 8, DOWN = 10, DOWNOPT = 15, CAL = 20, NEWCAL = 22, CALOPT = 25, PARDOWN = 30, PAROPTDOWN = 35, PARCAL = 40, PAROPTCAL = 45, END = 90, ERROR = 100};
   State state = START;
   State prevState = START;

   std::string variableName = "";
   std::string variableInputName = "";
   Options vOptions;
   Options viOptions;
   Options dOptions;
   std::vector<Options> cOptions;
   Options pOptions;
   std::string errorMessage = "";

   std::string downscaler = defaultDownscaler();
   std::string calibrator = "";
   Options cOpt;
   std::string parameterFile = "";
   std::vector<std::string> calibratorNames;
   std::vector<ParameterFile*> parameterFileCalibrators;
   ParameterFile* parameterFileDownscaler = NULL;
   std::string variableAlias = "";
   Options aOptions;
   while(true) {
      std::stringstream ss;
      ss << "State: " << state;
      Util::info(ss.str());

      if(state != ERROR)
         prevState = state;
      if(state == START) {
         if(index < argv.size() && argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else if(index < argv.size() && argv[index] == "-vi") {
            state = VARI;
            index++;
         }
         else if(index < argv.size() && argv[index] == "-va") {
            state = VARALIAS;
            index++;
         }
         else {
            errorMessage = "No variables defined";
            state = ERROR;
         }
      }
      else if(state == VAR) {
         if(argv.size() <= index) {
            // -v but nothing after it
            errorMessage = "No variable after '-v'";
            state = ERROR;
         }
         else if(argv[index][0] == '-') {
            errorMessage = "No variable name after '-v'";
            state = ERROR;
         }
         else {
            variableName = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(variableInputName == "" && argv[index] == "-vi") {
               state = VARI;
               index++;
            }
            else if(argv[index] == "-vi") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-p") {
               errorMessage = "-p must be after a -d or -c";
               state = ERROR;
            }
            else {
               state = VAROPT;
            }
         }
      }
      else if(state == VAROPT) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-d") {
            state = DOWN;
            index++;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-vi") {
            index++;
            state = VARI;
         }
         else if(argv[index] == "-p") {
            errorMessage = "-p must be after a -d or -c";
            state = ERROR;
         }
         else {
            // Process variable options
            vOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == VARI) {
         if(argv.size() <= index) {
            // -vi but nothing after it
            errorMessage = "No variable after '-vi'";
            state = ERROR;
         }
         else if(argv[index][0] == '-') {
            errorMessage = "No variable name after '-v'";
            state = ERROR;
         }
         else {
            variableInputName = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-vi") {
               state = NEWVAR;
            }
            else if(variableName == "" && argv[index] == "-v") {
               state = VAR;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-p") {
               errorMessage = "-p must be after a -d or -c";
               state = ERROR;
            }
            else {
               state = VARIOPT;
            }
         }
      }
      else if(state == VARIOPT) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-d") {
            state = DOWN;
            index++;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-vi") {
            state = NEWVAR;
         }
         else if(argv[index] == "-v") {
            state = VAR;
            index++;
         }
         else if(argv[index] == "-p") {
            errorMessage = "-p must be after a -d or -c";
            state = ERROR;
         }
         else {
            // Process variable options
            viOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWVAR) {
         if(variableName == "") {
            errorMessage = "-vi but no -v";
            state = ERROR;
         }
         else {
            // Check that we haven't added the variable before
            bool alreadyExists = false;
            for(int i = 0; i < variableConfigurations.size(); i++) {
               if(variableConfigurations[i].outputVariable.name() == variableName)
                  alreadyExists = true;
            }
            std::string useVariableInputName = variableInputName;
            Options useVariableInputOptions = viOptions;
            if(useVariableInputName == "") {
               useVariableInputName = variableName;
               useVariableInputOptions = vOptions;
            }

            if(!alreadyExists) {
               VariableConfiguration varconf;
               varconf.parameterFileDownscaler = parameterFileDownscaler;
               // Find input variable
               bool foundVarInput = inputFiles[0]->getVariable(useVariableInputName, varconf.inputVariable);
               varconf.inputVariable.add(useVariableInputOptions);
               if(!foundVarInput && downscaler != "bypass") {
                  std::stringstream ss;
                  ss << "Input file does not contain variable with name '" << useVariableInputName << "'. Use -d bypass if the variable will be diagnosed by a calibrator.";
                  Util::error(ss.str());
               }

               // Find output variable
               bool foundVarOutput = outputFiles[0]->getVariable(variableName, varconf.outputVariable);

               if(!foundVarOutput) {
                  foundVarOutput = inputFiles[0]->getVariable(variableName, varconf.outputVariable);
                  if(!foundVarOutput)
                     varconf.outputVariable = Variable(variableName);
               }
               varconf.outputVariable.add(vOptions);

               std::vector<Calibrator*> calibrators;
               for(int c = 0; c < calibratorNames.size(); c++) {
                  Calibrator* calibrator = Calibrator::getScheme(calibratorNames[c], varconf.outputVariable, cOptions[c]);
                  calibrators.push_back(calibrator);
               }
               varconf.calibrators = calibrators;
               varconf.parameterFileCalibrators = parameterFileCalibrators;
               varconf.inputVariableOptions = useVariableInputOptions;
               varconf.outputVariableOptions = vOptions;
               Downscaler* d = Downscaler::getScheme(downscaler, varconf.inputVariable, varconf.outputVariable, dOptions);
               varconf.downscaler = d;

               variableConfigurations.push_back(varconf);
            }
            else {
               std:: stringstream ss;
               ss << "Variable '" << variableName << "' already read. Using first instance.";
               Util::warning(ss.str());
            }

            // Reset to defaults
            vOptions.clear();
            viOptions.clear();
            variableName = "";
            variableInputName = "";
            downscaler = defaultDownscaler();
            parameterFileDownscaler = NULL;
            dOptions.clear();
            cOptions.clear();
            calibratorNames.clear();
            parameterFileCalibrators.clear();

            if(argv.size() <= index) {
               state = END;
            }
            else if(argv[index] == "-v"){
               state = VAR;
               index++;
            }
            else if(argv[index] == "-vi"){
               state = VARI;
               index++;
            }
         }
      }
      else if(state == VARALIAS) {
         if(argv.size() <= index) {
            // -va but nothing after it
            errorMessage = "No variable after '-va'";
            state = ERROR;
         }
         else if(argv[index][0] == '-') {
            errorMessage = "No variable name after '-va'";
            state = ERROR;
         }
         else {
            variableAlias = argv[index];
            index++;
            state = VARALIASOPT;
         }
      }
      else if(state == VARALIASOPT) {
         if(argv.size() <= index) {
            state = NEWVARALIAS;
         }
         else if(argv[index] == "-v") {
            state = NEWVARALIAS;
         }
         else if(argv[index] == "-vi") {
            state = NEWVARALIAS;
         }
         else {
            aOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWVARALIAS) {
         // Check that we haven't added the variable before
         std::map<std::string, Variable>::const_iterator it = variableAliases.find(variableAlias);
         if(it == variableAliases.end()) {
            std::string name = "";
            if(!aOptions.getValue("name", name)) {
               Util::error("Variable alias must have a name= option");
            }
            variableAliases[variableAlias] = Variable(aOptions);
            variableAlias = "";
            aOptions.clear();
         }
         else {
            std:: stringstream ss;
            ss << "Variable alias '" << variableAlias << "' already read. Using first instance.";
            Util::warning(ss.str());
         }
         if(argv.size() <= index) {
            state = END;
         }
         else if(argv[index] == "-v"){
            state = VAR;
            index++;
         }
         else if(argv[index] == "-vi"){
            state = VARI;
            index++;
         }
      }
      else if(state == DOWN) {
         if(argv.size() <= index) {
            // -d but nothing after it
            errorMessage = "No downscaler after '-d'";
            state = ERROR;
         }
         else {
            downscaler = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-vi") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               // Two downscalers defined for one variable
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-p") {
               state = PARDOWN;
            }
            else {
               state = DOWNOPT;
            }
         }
      }
      else if(state == DOWNOPT) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-vi") {
            state = NEWVAR;
         }
         else if(argv[index] == "-p") {
            state = PARDOWN;
         }
         else {
            // Process downscaler options
            dOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == PARDOWN) {
         if(argv.size() <= index) {
            // -p but nothing after it
            errorMessage = "No parameter file after '-p'";
            state = ERROR;
         }
         else {
            parameterFile = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-vi") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               // Two downscalers defined for one variable
               state = DOWN;
               index++;
            }
            else if(argv[index] == "-p") {
               // Two parameter files defined for one variable
               errorMessage = "Two or more -p used for one downscaler";
               state = ERROR;
               index++;
            }
            else {
               state = PAROPTDOWN;
            }
         }
      }
      else if(state == PAROPTDOWN) {
         if(argv.size() <= index) {
            state = NEWVAR;
         }
         else if(argv[index] == "-c") {
            state = CAL;
            index++;
         }
         else if(argv[index] == "-v") {
            state = NEWVAR;
         }
         else if(argv[index] == "-vi") {
            state = NEWVAR;
         }
         else if(argv[index] == "-p") {
            errorMessage = "Two or more -p used for one calibrator";
            state = ERROR;
         }
         else {
            // Process parameter file options
            pOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == CAL) {
         if(argv.size() <= index) {
            // -c but nothing after it
            state = ERROR;
            errorMessage = "No calibrator after '-c'";
         }
         else {
            calibrator = argv[index];
            index++;
            if(argv.size() <= index) {
               state = NEWCAL;
            }
            else if(argv[index] == "-v") {
               state = NEWCAL;
            }
            else if(argv[index] == "-vi") {
               state = NEWCAL;
            }
            else if(argv[index] == "-c") {
               state = NEWCAL;
            }
            else if(argv[index] == "-d") {
               state = NEWCAL;
            }
            else if(argv[index] == "-p") {
               state = PARCAL;
            }
            else {
               state = CALOPT;
            }
         }
      }
      else if(state == CALOPT) {
         if(argv.size() <= index) {
            state = NEWCAL;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else if(argv[index] == "-v") {
            state = NEWCAL;
         }
         else if(argv[index] == "-vi") {
            state = NEWCAL;
         }
         else if(argv[index] == "-d") {
            state = NEWCAL;
         }
         else if(argv[index] == "-p") {
            state = PARCAL;
         }
         else {
            // Process calibrator options
            cOpt.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == PARCAL) {
         if(argv.size() <= index) {
            // -p but nothing after it
            state = ERROR;
            errorMessage = "No parameter file after '-p'";
         }
         else {
            parameterFile = argv[index];
            if(argv.size() <= index) {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-vi") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-c") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-d") {
               state = NEWCAL;
               index++;
            }
            else if(argv[index] == "-p") {
               state = PARCAL;
               index++;
            }
            else {
               state = PAROPTCAL;
            }
         }
      }
      else if(state == PAROPTCAL) {
         if(argv.size() <= index) {
            state = NEWCAL;
         }
         else if(argv[index] ==  "-c") {
            state = NEWCAL;
         }
         else if(argv[index] == "-v") {
            state = NEWCAL;
         }
         else if(argv[index] == "-vi") {
            state = NEWCAL;
         }
         else if(argv[index] == "-d") {
            state = NEWCAL;
         }
         else if(argv[index] == "-p") {
            state = PARCAL;
         }
         else {
            // Process calibrator options
            pOptions.addOptions(argv[index]);
            index++;
         }
      }
      else if(state == NEWCAL) {
         // We do not need to check that the same calibrator has been added for this variable
         // since this is perfectly fine (e.g. smoothing twice).

         ParameterFile* p = NULL;
         if(parameterFile != "") {
            std::string schemeName;
            if(!pOptions.getValue("type", schemeName)) {
               state = ERROR;
               errorMessage = "Parameter file missing 'type': " + pOptions.toString();
            }
            else {
               pOptions.addOption("file", parameterFile);
               p = ParameterFile::getScheme(schemeName, pOptions);
               if(!p->isReadable()) {
                  state = ERROR;
                  errorMessage = "Could not open parameter file: " + pOptions.toString();
               }
            }
         }
         if(state != ERROR) {
            calibratorNames.push_back(calibrator);
            cOptions.push_back(cOpt);
            parameterFileCalibrators.push_back(p);

            // Reset
            calibrator = "";
            cOpt.clear();
            parameterFile = "";
            if(argv.size() <= index) {
               state = NEWVAR;
            }
            else if(argv[index] == "-c") {
               state = CAL;
               index++;
            }
            else if(argv[index] == "-v") {
               state = NEWVAR;
            }
            else if(argv[index] == "-vi") {
               state = NEWVAR;
            }
            else if(argv[index] == "-d") {
               state = DOWN;
               index++;
            }
            else {
               state = ERROR;
               errorMessage = "No recognized option after '-c calibrator'";
            }
         }
      }
      else if(state == END) {
         break;
      }
      else if(state == ERROR) {
         std::stringstream ss;
         ss << "Invalid command line arguments: " << errorMessage << ".";
         Util::error(ss.str());
      }
   }

   // Add aliases to input files
   for(int i = 0; i < inputFiles.size(); i++) {
      std::map<std::string, Variable>::const_iterator it = variableAliases.find(variableAlias);
      for(it = variableAliases.begin(); it != variableAliases.end(); it++) {
         inputFiles[i]->addVariableAlias(it->first, it->second);
      }
   }
}
Setup::~Setup() {
   std::map<std::string, File*>::const_iterator it = mFileMap.begin();
   for(it = mFileMap.begin(); it != mFileMap.end(); it++)
      delete it->second;

   for(int i = 0; i < variableConfigurations.size(); i++) {
      delete variableConfigurations[i].downscaler;
      for(int c = 0; c < variableConfigurations[i].calibrators.size(); c++) {
         delete variableConfigurations[i].calibrators[c];
      }
   }
}
std::string Setup::defaultDownscaler() {
   return "nearestNeighbour";
}

bool Setup::hasFile(std::string iFilename) const {
   std::map<std::string, File*>::const_iterator it = mFileMap.find(iFilename);
   return it != mFileMap.end();
}
