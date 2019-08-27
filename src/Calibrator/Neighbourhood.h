#ifndef CALIBRATOR_NEIGHBOURHOOD_H
#define CALIBRATOR_NEIGHBOURHOOD_H
#include "Calibrator.h"
#include "../Variable.h"
#include "../Util.h"

class ParameterFile;
class Parameters;

//! Applies a statistical operator to a neighbourhood
class CalibratorNeighbourhood : public Calibrator {
   public:
      CalibratorNeighbourhood(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "neighbourhood";};
      bool requiresParameterFile() const { return false;};
      void calibrateField(const Field& iInput, Field& iOutput, const ParameterFile* iParameterFile=NULL, int t=0) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      int mRadius;
      Util::StatType mStatType;
      float mQuantile;
      bool mFast;
      bool mApprox;
      int numMissingValues(const Field& iField, int iEnsIndex) const;
};
#endif
