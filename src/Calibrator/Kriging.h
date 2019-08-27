#ifndef CALIBRATOR_KRIGING_H
#define CALIBRATOR_KRIGING_H
#include "Calibrator.h"
#include "../ParameterFile/ParameterFile.h"
class Obs;
class Forecast;
class Parameters;

class CalibratorKriging : public Calibrator {
   public:
      CalibratorKriging(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      float calcCovar(const Location& loc1, const Location& loc2) const;
      enum Type {
         TypeCressman = 10,
         TypeBarnes   = 20
      };
      //! Compute the bias at the training point
      Parameters train(const std::vector<ObsEns>& iData) const;
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;

      float mRadius;
      float mMaxElevDiff;
      float mEfoldDist;
      std::string name() const {return "kriging";};
      File* mPrevious;
      std::string mAuxVariable;
      float mLowerThreshold;
      float mUpperThreshold;
      Type mKrigingType;
      bool mUseApproxDistance;
      Util::Operator mOperator;
      bool mCrossValidate;
      int mWindow;
};
#endif
