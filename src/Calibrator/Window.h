#ifndef CALIBRATOR_WINDOW_H
#define CALIBRATOR_WINDOW_H
#include "Calibrator.h"
#include "../Variable.h"

// Applies a statistical operator to values within a temporal window
class CalibratorWindow : public Calibrator {
   public:
      CalibratorWindow(const Variable& iVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "window";};
      bool requiresParameterFile() const { return false;};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      int mLength;
      Util::StatType mStatType;
      float mQuantile;
      bool mBefore;
      bool mKeepMissing;
      std::vector<float> mWeights;
      enum EdgePolicy {
         EdgePolicyCompute = 0,
         EdgePolicyMissing = 10
      };
      EdgePolicy mEdgePolicy;
};
#endif
