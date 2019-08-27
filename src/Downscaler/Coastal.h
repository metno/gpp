#ifndef DOWNSCALER_COASTAL_H
#define DOWNSCALER_COASTAL_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;
class DownscalerCoastal : public Downscaler {
   public:
      //! Downscale the specified variable
      DownscalerCoastal(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);
      std::vector<int> getSearchRadii() const;
      float getMinLafDiff() const;
      float getMinGradient() const;
      float getMaxGradient() const;
      static std::string description(bool full=true);
      std::string name() const {return "coastal";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      std::vector<int>  mSearchRadii;
      std::vector<float>  mWeights;
      int mLafRadius;
      float mMinLafDiff; // Minimum elevation difference within neighbourhood to use gradient
      float mMinGradient;
      float mMaxGradient;
      float mElevGradient;
};
#endif
