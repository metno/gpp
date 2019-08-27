#ifndef METCAL_SETUP_TRAIN_H
#define METCAL_SETUP_TRAIN_H
#include <string>
#include <vector>
#include "Variable.h"
#include "Options.h"
class Calibrator;
class Downscaler;
class ParameterFile;
class File;

//! Represents what post-processing method should be trained. Includes which data file to
//! read data from, which output file to place the parameters in, which variable to post-process,
//! and what post-processing method to train.
class SetupTrain {
   public:
      ParameterFile* output;
      std::vector<File*> forecasts;
      std::vector<File*> observations;
      Calibrator* method;
      Downscaler* downscaler;
      Variable variable;
      SetupTrain(const std::vector<std::string>& argv);
      ~SetupTrain();
};
#endif
