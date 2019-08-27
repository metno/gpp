#include "Accumulate.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorAccumulate::CalibratorAccumulate(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions) {
   iOptions.check();
}
bool CalibratorAccumulate::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nY = iFile.getNumY();
   int nX = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   std::vector<FieldPtr> fields(nTime);
   std::vector<FieldPtr> fieldsAcc(nTime);
   for(int t = 0; t < nTime; t++) {
      fields[t]    = iFile.getField(mVariable, t);
      fieldsAcc[t] = iFile.getEmptyField();
   }

   for(int t = 0; t < nTime; t++) {
      #pragma omp parallel for
      for(int y = 0; y < nY; y++) {
         for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nEns; e++) {
               if(t == 0) {
                  (*fieldsAcc[t])(y, x, e) = 0;
               }
               else {
                  float previous = (*fieldsAcc[t - 1])(y, x, e);
                  float current  = (*fields[t])(y, x, e);
                  if(Util::isValid(current) && Util::isValid(previous)) {
                     (*fieldsAcc[t])(y, x, e) = current + previous;
                  }
                  else {
                     (*fieldsAcc[t])(y, x, e) = Util::MV;
                  }
               }
            }
         }
      }

      iFile.addField(fieldsAcc[t], mVariable, t);
   }
   return true;
}
std::string CalibratorAccumulate::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-c accumulate","Accumlates a variable over time") << std::endl;
   return ss.str();
}
