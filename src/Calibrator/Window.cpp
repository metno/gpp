#include "Window.h"
#include "../Util.h"
#include "../File/File.h"
CalibratorWindow::CalibratorWindow(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mLength(7),
      mBefore(false),
      mStatType(Util::StatTypeMean),
      mEdgePolicy(EdgePolicyCompute),
      mKeepMissing(false),
      mQuantile(Util::MV) {
   iOptions.getValue("length", mLength);
   iOptions.getValue("before", mBefore);
   iOptions.getValue("keepMissing", mKeepMissing);
   std::string edgePolicy = "compute";
   iOptions.getValue("edgePolicy", edgePolicy);
   if(edgePolicy == "compute")
      mEdgePolicy = EdgePolicyCompute;
   else if(edgePolicy == "missing")
      mEdgePolicy = EdgePolicyMissing;
   else {
      std::stringstream ss;
      ss << "Cannot understand edgePolicy=" << edgePolicy;
      Util::error(ss.str());
   }

   if(mLength < 1) {
      std::stringstream ss;
      ss << "CalibratorWindow: 'length' (" << mLength << ") must be > 0";
      Util::error(ss.str());
   }

   std::string op;
   if(iOptions.getValue("stat", op)) {
      bool status = Util::getStatType(op, mStatType);
      if(!status) {
         std::stringstream ss;
         ss << "Could not recognize stat=" << op;
         Util::error(ss.str());
      }
   }
   if(mStatType == Util::StatTypeQuantile) {
      iOptions.getRequiredValue("quantile", mQuantile);
      if(!Util::isValid(mQuantile) || mQuantile < 0 || mQuantile > 1) {
         Util::error("'quantile' must be on the interval [0,1]");
      }
   }
   iOptions.check();
}
bool CalibratorWindow::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();

   // Get all fields
   std::vector<FieldPtr> fields(nTime);
   std::vector<FieldPtr> fieldsOrig(nTime);
   for(int t = 0; t < nTime; t++) {
      fieldsOrig[t]     = iFile.getField(mVariable, t);
      fields[t] = iFile.getEmptyField();
   }

   for(int t = 0; t < nTime; t++) {
      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            for(int e = 0; e < nEns; e++) {
               std::vector<float> window;
               window.resize(mLength, Util::MV);
               int index = 0;
               int before = (mLength-1) / 2;
               int after = mLength - before - 1;
               int start = std::max(0, t - before);
               int end = std::min(nTime-1, t + after);
               if(mBefore) {
                  start = std::max(0, t - (mLength - 1));
                  end = std::min(nTime-1, t);
               }
               assert(end - start + 1 <= mLength);
               if(mEdgePolicy == EdgePolicyMissing && (end - start + 1 != mLength)) {
                  (*fields[t])(i,j,e) = Util::MV;
               }
               else {
                  bool numMissing = 0;
                  for(int tt = start; tt <= end; tt++) {
                     float curr = (*fieldsOrig[tt])(i,j,e);
                     numMissing += !Util::isValid(curr);
                     window[index] = curr;
                     index++;
                  }
                  if(mKeepMissing && numMissing > 0)
                     (*fields[t])(i,j,e) = Util::MV;
                  else
                     (*fields[t])(i,j,e) = Util::calculateStat(window, mStatType, mQuantile);
               }
            }
         }
      }
   }
   for(int t = 0; t < nTime; t++) {
      iFile.addField(fields[t], mVariable, t);
   }
   return true;
}

std::string CalibratorWindow::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-c window","Applies a statistical operator to values within a temporal window. Any missing values are ignored when computing the statistic.") << std::endl;
   if(full) {
      ss << Util::formatDescription("   length=required", "Length of the window (in number of timesteps) to apply operator on (must be 0 or greater).") << std::endl;
      ss << Util::formatDescription("   before=0", "If 0, the window is centered on each leadtime (if length is an even number, then it is shifted such that it includes one extra future leadtime). If 1, then the window ends and includes the leadtime.") << std::endl;
      ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'quantile', 'std', or 'sum'. 'std' is the population standard deviation.") << std::endl;
      ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
      ss << Util::formatDescription("   keepMissing=0", "Should missing values cause the window to be missing? If not, then ignore missing values in calculation.") << std::endl;
      ss << Util::formatDescription("   edgePolicy=compute", "What policy should be used on edges? Either 'compute' to compute as usual, or 'missing' to set missing value if the window is not full.") << std::endl;
   }
   return ss.str();
}
