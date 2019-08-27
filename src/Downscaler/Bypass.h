#ifndef DOWNSCALER_BYPASS_H
#define DOWNSCALER_BYPASS_H
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
typedef std::vector<std::vector<int> > vec2Int;

//! Skip the downscaling stage
class DownscalerBypass : public Downscaler {
   public:
      DownscalerBypass(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);
      static std::string description(bool full=true);
      std::string name() const {return "bypass";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
};
#endif
