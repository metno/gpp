#include "Bypass.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerBypass::DownscalerBypass(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerBypass::downscaleCore(const File& iInput, File& iOutput) const {
   if(!iOutput.hasVariable(mOutputVariable)) {
      // int nTime = iInput.getNumTime();
      // for(int t = 0; t < nTime; t++) {
      iOutput.initNewVariable(mOutputVariable);
      // }
   }
}

std::string DownscalerBypass::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-d bypass", "Skip downscaling, used when diagnosed fields") << std::endl;
   return ss.str();
}
