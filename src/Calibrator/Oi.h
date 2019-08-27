#ifndef CALIBRATOR_OI_H
#define CALIBRATOR_OI_H
#include <armadillo>
#include "Calibrator.h"
#include "../ParameterFile/ParameterFile.h"
class Obs;
class Forecast;
class Parameters;

class CalibratorOi : public Calibrator {
   public:
      CalibratorOi(Variable iVariable, const Options& iOptions);
      // Compute rho. A rho of 0 is returned if a vertical distance is missing when the
      // vertical scale is defined.
      enum RhoType {RhoTypeGaussian, RhoTypeSoar};
      float calcRho(float iHdist, float iVdist, float iLdist, RhoType iType=RhoTypeGaussian) const;
      static std::string description(bool full=true);
      std::string name() const {return "oi";};
   private:
      bool calibrateCore(File& iFile, const ParameterFile* iParameterFile) const;
      enum Type {TypeTemperature, TypePrecipitation};
      enum TransformType {TransformTypeNone, TransformTypeBoxCox};
      float mVLength;
      float mHLength;
      float mHLengthC;
      float mMu;
      float mGamma;
      std::string mBiasVariable;
      int mMaxLocations;
      float mSigma;
      float mSigmaC;
      float mDelta;
      std::string mDeltaVariable;
      std::string mNumVariable;
      float mC;
      float mEpsilonC;
      float mEpsilon;
      int mX;
      int mY;
      float mMinRho;
      float mMaxBytes;
      int mMinValidEns;
      bool mSaveDiff;
      float mElevGradient;
      bool mExtrapolate;
      float mWMin;
      float mNewDeltaVar;
      float mMaxElevDiff;
      bool mDiagnose;
      bool mLandOnly;
      std::string mDiaFile;
      bool mUseEns;
      typedef arma::mat mattype;
      typedef arma::vec vectype;
      typedef arma::cx_mat cxtype;
      float mLambda;
      bool mCrossValidate;
      Type mType;
      TransformType mTransformType;
      float calcDelta(float iOldDelta, const vec2& iY) const;
      float transform(float iValue) const;
      float invTransform(float iValue) const;
      RhoType mRhoType;
      float mBoxCoxThreshold;
};
#endif
