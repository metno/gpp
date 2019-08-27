#include "Cloud.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorCloud::CalibratorCloud(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mValue(1),
      mPrecipVariable("") {
   iOptions.getRequiredValue("precipVariable", mPrecipVariable);
   iOptions.getValue("value", mValue);
   iOptions.check();
}
bool CalibratorCloud::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      const Field& precip = *iFile.getField(mPrecipVariable, t);
      Field& cloud        = *iFile.getField(mVariable, t);

      // TODO: Figure out which cloudless members to use. Ideally, if more members
      // need precip, we should pick members that already have clouds, so that we minimize
      // our effect on the cloud cover field.

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // Turn on clouds if needed, i.e don't allow a member to
            // have precip without cloud cover.
            for(int e = 0; e < nEns; e++) {
               float currPrecip = precip(i,j,e);
               float currCloud  = cloud(i,j,e);
               if(Util::isValid(currPrecip) && Util::isValid(currCloud)) {
                  cloud(i,j,e)  = currCloud;
                  if(currPrecip > 0 && currCloud < mValue) {
                     cloud(i,j,e) = mValue;
                  }
               }
            }
         }
      }
   }
   return true;
}
std::string CalibratorCloud::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c cloud", "Ensure a minimum cloud cover value precipitation is present") << std::endl;
      ss << Util::formatDescription("   precipVariable=undef", "Name of precipitation variable") << std::endl;
      ss << Util::formatDescription("   value=1", "Minimum cloud cover value allowed") << std::endl;
   }
   else
      ss << Util::formatDescription("-c cloud", "Ensure clouds when precip is present") << std::endl;
   return ss.str();
}
