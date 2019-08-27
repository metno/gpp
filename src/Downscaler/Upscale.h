#ifndef DOWNSCALER_UPSCALE_H
#define DOWNSCALER_UPSCALE_H
#include <map>
#include "Downscaler.h"
#include "../Variable.h"
#include "../Util.h"
#include "../Field.h"
class File;
class DownscalerUpscale : public Downscaler {
   public:
      DownscalerUpscale(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions);

      static std::string description(bool full=true);
      std::string name() const {return "upscale";};
   private:
      void downscaleCore(const File& iInput, File& iOutput) const;
      Util::StatType mStatType;
      float mQuantile;
};
#endif
