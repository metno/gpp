#include "Coastal.h"
#include <cmath>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Downscaler/Pressure.h"
#include "../Grid.h"
CalibratorCoastal::CalibratorCoastal(const Variable& iVariable, const Options& iOptions) :
      Calibrator(iVariable, iOptions),
      mSearchRadius(3),
      mMinLafDiff(0.1),
      mUseNN(false) {
   iOptions.getValue("searchRadius", mSearchRadius);
   iOptions.getValue("minLafDiff", mMinLafDiff);
   iOptions.getValue("useNN", mUseNN);
   iOptions.check();
}
bool CalibratorCoastal::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();
   vec2 laf = iFile.getLandFractions();

   if(iParameterFile->getNumParameters() == 0) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must have at least one dataacolumns");
   }
   if(iParameterFile->isLocationDependent() == 0) {
      Util::error("Parameter file '" + iParameterFile->getFilename() + "' must be spatial");
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Parameters parametersGlobal;
      if(!iParameterFile->isLocationDependent())
          parametersGlobal = iParameterFile->getParameters(t);
      const FieldPtr field = iFile.getField(mVariable, t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            Parameters parameters;
            if(iParameterFile->isLocationDependent())
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
            else
               parameters = parametersGlobal;
            for(int e = 0; e < nEns; e++) {
               float lowerValue = Util::MV;
               float upperValue = Util::MV;
               float minLaf = 2;
               float maxLaf = -1;
               float nearestNeighbour = (*field)(i,j,e);
               // Find the range within a neighbourhood
               for(int ii = std::max(0, i-mSearchRadius); ii <= std::min(nLat-1, i+mSearchRadius); ii++) {
                  for(int jj = std::max(0, j-mSearchRadius); jj <= std::min(nLon-1, j+mSearchRadius); jj++) {
                     float currLaf = laf[ii][jj];
                     if(currLaf < minLaf) {
                        lowerValue = (*field)(ii,jj,e);
                        minLaf = currLaf;
                     }
                     if(currLaf > maxLaf) {
                        upperValue = (*field)(ii,jj,e);
                        maxLaf = currLaf;
                     }
                  }
               }

               float a = parameters[0];
               float b = parameters[1];
               float c = parameters[2];
               // If there is not sufficient variability in LSM in the neighbourhood, then
               // force don't use the range in the regression
               float gradient = 0;
               if(Util::isValid(lowerValue) && Util::isValid(upperValue) && (maxLaf - minLaf > mMinLafDiff))
                  gradient = (upperValue - lowerValue) / (maxLaf - minLaf);
               assert(Util::isValid(gradient));
               assert(Util::isValid(nearestNeighbour));

               float value = Util::MV;
               if(mUseNN)
                  value = a + b * nearestNeighbour + c * gradient;
               else
                  value = a + b * lowerValue + c * gradient;
               assert(Util::isValid(value));
               (*field)(i,j,e) = value;
            }
         }
      }
   }
   return true;
}

Parameters CalibratorCoastal::train(const std::vector<ObsEnsField>& iData, const Grid& iObsGrid, const Grid& iEnsGrid, int iIobs, int iJobs, int iIens, int iJens) const {
   std::vector<float> y;
   int nLat = iData[0].second->getNumY();
   int nLon = iData[0].second->getNumX();
   vec2 laf = iEnsGrid.landFractions();

   // Find nearest land and sea point
   float lowerValue = Util::MV;
   float upperValue = Util::MV;
   int Imin = Util::MV;
   int Imax = Util::MV;
   int Jmin = Util::MV;
   int Jmax = Util::MV;
   float minLaf = 2;
   float maxLaf = -1;
   // Find the range within a neighbourhood
   for(int ii = std::max(0, iIens-mSearchRadius); ii <= std::min(nLat-1, iIens+mSearchRadius); ii++) {
      for(int jj = std::max(0, iJens-mSearchRadius); jj <= std::min(nLon-1, iJens+mSearchRadius); jj++) {
         float currLaf = laf[ii][jj];
         if(currLaf < minLaf) {
            Imin = ii;
            Jmin = jj;
            minLaf = currLaf;
         }
         if(currLaf > maxLaf) {
            Imax = ii;
            Jmax = jj;
            maxLaf = currLaf;
         }
      }
   }

   bool useRange = (maxLaf - minLaf) > mMinLafDiff;
   if(!useRange && 0)
      std::cout << "not using range for " << iIens << "," << iJens << " range = " << (maxLaf - minLaf) << std::endl;
   // useRange = 0;

   assert(Util::isValid(Imax));
   assert(Util::isValid(Jmax));

   // Compute predictors in model
   std::vector<float> values;
   std::vector<std::vector<float> > x(2);
   for(int i = 0; i < iData.size(); i++) {
      float obs = (*iData[i].first)(iIobs, iJobs, 0);
      std::vector<float> ensMin = (*iData[i].second)(Imin, Jmin);
      float valueMin = Util::calculateStat(ensMin, Util::StatTypeMean);
      std::vector<float> ensMax = (*iData[i].second)(Imax, Jmax);
      float valueMax = Util::calculateStat(ensMax, Util::StatTypeMean);
      std::vector<float> ensNn = (*iData[i].second)(iIens, iJens);
      float valueNN = Util::calculateStat(ensNn, Util::StatTypeMean);
      if(Util::isValid(obs) && Util::isValid(valueNN)) {
         float gradient = 0;
         if(Util::isValid(valueMin) && Util::isValid(valueMax))
            gradient = (valueMax - valueMin)/(maxLaf - minLaf);
         if(mUseNN)
            y.push_back(obs-valueNN);
         else
            y.push_back(obs-valueMin);
         x[0].push_back(1);
         x[1].push_back(gradient);
         if(iIobs == 54 && iJobs == 16) {
            std::cout << obs-valueNN << " " << gradient << std::endl;
         }
      }
   }

   if(y.size() <= 0) {
      std::stringstream ss;
      ss << "CalibratorCoastal: Cannot train coastal, no valid data";
      Util::error(ss.str());
   }

   // a + 1*Tnn + c * range
   // a + 1*Tsea + c * range
   if(useRange) {
      std::vector<float> values2 = Util::regression(y, x, false);
      values.push_back(values2[0]);
      values.push_back(1);
      values.push_back(values2[1]);
      // values.push_back(0);
      // values.push_back(1);
      // values.push_back(0);
   }
   // a + 1*Tnn + 0 * range
   else {
      float total = 0;
      int counter = 0;
      for(int i = 0; i < y.size(); i++) {
         if(Util::isValid(y[i])) {
            total += y[i];
            counter++;
         }
      }
      float a = 0;
      if(counter > 0)
         a = total / counter;
      values.push_back(a);
      values.push_back(1);
      values.push_back(0);
      assert(Util::isValid(a));
   }

   Parameters par(values);

   return par;
}

std::string CalibratorCoastal::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-c coastal", "Weights land and sea forecasts") << std::endl;
   if(full) {
      ss << Util::formatDescription("   searchRadius=3", "Radius (in number of gridpoints) to search for land and sea points.") << std::endl;
      ss << Util::formatDescription("   minLafDiff=0.1", "Minimum difference in land area fraction in order to use the  land and sea points") << std::endl;
      ss << Util::formatDescription("   useNN=False", "Use the nearest neighbour as the basevalue, instead of the point with the lowest LAF.") << std::endl;
   }
   return ss.str();
}
