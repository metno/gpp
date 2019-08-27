#ifndef CALIBRATOR_H
#define CALIBRATOR_H
#include <string>
#include <vector>
#include "../Scheme.h"
#include "../Parameters.h"
#include "../Field.h"
#include "../Variable.h"

typedef std::vector<float> Ens;
typedef std::pair<float,Ens> ObsEns;
typedef std::pair<FieldPtr,FieldPtr> ObsEnsField;
class File;
class Options;
class ParameterFile;
class Grid;

//! Abstract calibration class
class Calibrator : public Scheme {
   public:
      //! The calibrator does not free the memory of iParameterFile
      Calibrator(const Variable& iVariable, const Options& iOptions);
      virtual ~Calibrator() {};
      //! \brief Calibrate one or more fields in iFile
      //! @return true if calibration was successful, false otherwise
      bool calibrate(File& iFile, const ParameterFile* iParameterFile=NULL) const;

      //! Instantiates a calibrator with name iName
      static Calibrator* getScheme(std::string iName, Variable iVariable, const Options& iOptions);

      //! \brief Ensure that values in iAfter are in the same order as in iBefore.
      //! If missing values are encountered in iBefore or iAfter, then iAfter is left unchanged.
      //! If the sizes are different, then iAfter is left unchanged.
      static void  shuffle(const std::vector<float>& iBefore, std::vector<float>& iAfter);

      //! Returns the name of this calibrator
      virtual std::string name() const = 0;

      static std::string getDescriptions(bool full=true);

      virtual Parameters train(const std::vector<ObsEns>& iData) const;
      // Defaults to using the regular train method
      virtual Parameters train(const std::vector<ObsEnsField>& iData, const Grid& iObsGrid, const Grid& iEnsGrid, int iIobs, int iJobs, int iIEns, int iJEns) const;

      // Does this calibrator require a parameter file?
      virtual bool requiresParameterFile() const { return true;};
      Options getOptions() const;
   protected:
      virtual bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const = 0;
      Variable mVariable;
   private:
      Options mOptions;
};
// #include "Wind.h"
#include "Accumulate.h"
#include "Altitude.h"
#include "Bct.h"
#include "Cloud.h"
#include "Coastal.h"
#include "Deaccumulate.h"
#include "DiagnoseHumidity.h"
#include "DiagnoseWind.h"
#include "Gaussian.h"
#include "Kriging.h"
#include "Mask.h"
#include "Neighbourhood.h"
#include "Oi.h"
#include "Override.h"
#include "Phase.h"
#include "Qc.h"
#include "Qnh.h"
#include "Qq.h"
#include "Regression.h"
#include "Sort.h"
#include "Temperature.h"
#include "Threshold.h"
#include "WindDirection.h"
#include "Window.h"
#include "Zaga.h"
#endif
