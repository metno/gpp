#include "Neighbourhood.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
CalibratorNeighbourhood::CalibratorNeighbourhood(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mRadius(3),
      mStatType(Util::StatTypeMean),
      mFast(true),
      mApprox(false),
      mQuantile(Util::MV) {
   iOptions.getValue("radius", mRadius);
   iOptions.getValue("fast", mFast);
   iOptions.getValue("approx", mApprox);
   if(mRadius < 0) {
      std::stringstream ss;
      ss << "CalibratorNeighbourhood: 'radius' (" << mRadius << ") must be >= 0";
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

bool CalibratorNeighbourhood::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nTime = iFile.getNumTime();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& output = *iFile.getField(mVariable, t);
      Field raw = output;

      calibrateField(raw, output, iParameterFile);
   }
   return true;
}
int CalibratorNeighbourhood::numMissingValues(const Field& iField, int iEnsIndex) const {
   int count = 0;
   for(int x = 0; x < iField.getNumX(); x++) {
      for(int y = 0; y < iField.getNumY(); y++) {
         count += !Util::isValid(iField(y, x, iEnsIndex));
      }
   }
   return count;
}
void CalibratorNeighbourhood::calibrateField(const Field& iInput, Field& iOutput, const ParameterFile* iParameterFile, int t) const {
   double start_time = Util::clock();
   int radius = mRadius;
   int nEns = iInput.getNumEns();
   int nLat = iInput.getNumY();
   int nLon = iInput.getNumX();
   if(iParameterFile != NULL) {
      if(iParameterFile->isLocationDependent()) {
         Util::error("Cannot use a location dependent parameter file for CalibratorNeighbourhood");
      }
      radius = iParameterFile->getParameters(t)[0];
   }

   int count_stat = 0;
   for(int e = 0; e < nEns; e++) {
      if(mFast && (mStatType == Util::StatTypeMean || mStatType == Util::StatTypeSum)) {
         vec2 values;
         vec2 counts;
         values.resize(nLat);
         counts.resize(nLat);
         for(int i = 0; i < nLat; i++) {
            values[i].resize(nLon, 0);
            counts[i].resize(nLon, 0);
         }
         // Compute accumulated values
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               float value = iInput(i, j, e);
               if(j == 0 && i == 0) {
                  // Lower corner
                  if(Util::isValid(value)) {
                     values[i][j] = iInput(i, j, e);
                     counts[i][j] = 1;
                  }
               }
               else if(j == 0) {
                  // Lower row
                  if(Util::isValid(value)) {
                     values[i][j] = values[i-1][j] + iInput(i,j,e);
                     counts[i][j] = counts[i-1][j] + 1;
                  }
                  else {
                     values[i][j] = values[i-1][j];
                     counts[i][j] = counts[i-1][j];
                  }
               }
               else if(i == 0) {
                  // Left column
                  if(Util::isValid(value)) {
                     values[i][j] = values[i][j-1] + iInput(i,j,e);
                     counts[i][j] = counts[i][j-1] + 1;
                  }
                  else {
                     values[i][j] = values[i][j-1];
                     counts[i][j] = counts[i][j-1];
                  }

               }
               else {
                  if(Util::isValid(value)) {
                     values[i][j] = values[i][j-1] + values[i-1][j] - values[i-1][j-1] + iInput(i,j,e);
                     counts[i][j] = counts[i][j-1] + counts[i-1][j] - counts[i-1][j-1] + 1;
                  }
                  else {
                     values[i][j] = values[i][j-1] + values[i-1][j] - values[i-1][j-1];
                     counts[i][j] = counts[i][j-1] + counts[i-1][j] - counts[i-1][j-1];
                  }
               }
            }
         }
         // Put neighbourhood into vector
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               int i1 = std::min(nLat-1, i + mRadius);
               int j1 = std::min(nLon-1, j + mRadius);
               int i0 = i - mRadius - 1;
               int j0 = j - mRadius - 1;
               float value11 = values[i1][j1];
               float value00 = 0;
               float value10 = 0;
               float value01 = 0;
               int count11 = counts[i1][j1];
               int count00 = 0;
               int count10 = 0;
               int count01 = 0;
               if(i0 >= 0 && j0 >= 0) {
                  value00 = values[i0][j0];
                  value10 = values[i1][j0];
                  value01 = values[i0][j1];
                  count00 = counts[i0][j0];
                  count10 = counts[i1][j0];
                  count01 = counts[i0][j1];
               }
               else if(j0 >= 0) {
                  value10 = values[i1][j0];
                  count10 = counts[i1][j0];
               }
               else if(i0 >= 0) {
                  value01 = values[i0][j1];
                  count01 = counts[i0][j1];
               }
               float value = value11 + value00 - value10 - value01;
               int count = count11 + count00 - count10 - count01;
               if(count > 0) {
                  if(mStatType == Util::StatTypeMean) {
                     value /= count;
                  }
                  iOutput(i,j,e) = value;
               }
            }
         }
      }
      else if(numMissingValues(iInput, e) == 0 && (
               (mFast && (mStatType == Util::StatTypeMin || mStatType == Util::StatTypeMax)) ||
               (mApprox && (mStatType == Util::StatTypeMedian || mStatType == Util::StatTypeQuantile)))) {
         // Compute min/max quickly or any other quantile in a faster, but approximate way
         vec2 values;
         values.resize(nLat);
         for(int i = 0; i < nLat; i++) {
            values[i].resize(nLon, 0);
         }
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            if(i < mRadius || i >= nLat - mRadius) {
               // Regular way
               for(int j = 0; j < nLon; j++) {
                  // Put neighbourhood into vector
                  std::vector<float> neighbourhood;
                  int Ni = std::min(nLat-1, i+radius) - std::max(0, i-radius) + 1;
                  int Nj = std::min(nLon-1, j+radius) - std::max(0, j-radius) + 1;
                  assert(Ni > 0);
                  assert(Nj > 0);
                  neighbourhood.resize(Ni*Nj, Util::MV);
                  int index = 0;
                  for(int ii = std::max(0, i-radius); ii <= std::min(nLat-1, i+radius); ii++) {
                     for(int jj = std::max(0, j-radius); jj <= std::min(nLon-1, j+radius); jj++) {
                        float value = iInput(ii,jj,e);
                        assert(index < Ni*Nj);
                        neighbourhood[index] = value;
                        index++;
                     }
                  }
                  assert(index == Ni*Nj);
                  values[i][j] = Util::calculateStat(neighbourhood, mStatType, mQuantile);
                  count_stat += neighbourhood.size();
               }
            }
            else {
               // Fast way: Compute stats on each sliver
               std::vector<float> slivers(nLon, 0);
               for(int j = 0; j < nLon; j++) {
                  std::vector<float> sliver(2*mRadius+1, 0);
                  int count = 0;
                  for(int ii = i - mRadius; ii <= i + mRadius; ii++) {
                     sliver[count] = iInput(ii, j, e);
                     count++;
                  }
                  slivers[j] = Util::calculateStat(sliver, mStatType, mQuantile);
                  count_stat += sliver.size();
               }
               for(int j = 0; j < nLon; j++) {
                  std::vector<float> curr;
                  for(int jj = std::max(0, j - mRadius); jj <= std::min(nLon-1, j + mRadius); jj++) {
                     curr.push_back(slivers[jj]);
                  }
                  values[i][j] = Util::calculateStat(curr, mStatType, mQuantile);
                  count_stat += curr.size();
               }
            }
         }
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               iOutput(i,j,e) = values[i][j];
            }
         }
      }
      else {
         // Compute by brute force
         vec2 values;
         values.resize(nLat);
         for(int i = 0; i < nLat; i++) {
            values[i].resize(nLon, 0);
         }
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               // Put neighbourhood into vector
               std::vector<float> neighbourhood;
               int Ni = std::min(nLat-1, i+radius) - std::max(0, i-radius) + 1;
               int Nj = std::min(nLon-1, j+radius) - std::max(0, j-radius) + 1;
               assert(Ni > 0);
               assert(Nj > 0);
               neighbourhood.resize(Ni*Nj, Util::MV);
               int index = 0;
               for(int ii = std::max(0, i-radius); ii <= std::min(nLat-1, i+radius); ii++) {
                  for(int jj = std::max(0, j-radius); jj <= std::min(nLon-1, j+radius); jj++) {
                     float value = iInput(ii,jj,e);
                     assert(index < Ni*Nj);
                     neighbourhood[index] = value;
                     index++;
                  }
               }
               assert(index == Ni*Nj);
               values[i][j] = Util::calculateStat(neighbourhood, mStatType, mQuantile);
               count_stat += neighbourhood.size();
            }
         }
         #pragma omp parallel for
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               iOutput(i,j,e) = values[i][j];
            }
         }
      }
   }
   std::stringstream ss;
   ss << "Number of neighbourhood stat calculations: " << count_stat << " " << Util::clock() - start_time;
   Util::info(ss.str());
}

std::string CalibratorNeighbourhood::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood (example by averaging across a neighbourhood thereby smoothing the field).") << std::endl;
      ss << Util::formatDescription("   radius=3", "Use gridpoints within this number of points within in both east-west and north-south direction. The radius can alternatively be specified using a location-independent parameter file, with one parameter.") << std::endl;
      ss << Util::formatDescription("   stat=mean", "What statistical operator should be applied to the neighbourhood? One of 'mean', 'median', 'min', 'max', 'quantile', 'std', or 'sum'. 'std' is the population standard deviation.") << std::endl;
      ss << Util::formatDescription("   quantile=undef", "If stat=quantile is selected, what quantile (number on the interval [0,1]) should be used?") << std::endl;
      ss << Util::formatDescription("   fast=1", "Use shortcuts to compute 'max', 'min', 'mean' or 'sum' faster.") << std::endl;
      ss << Util::formatDescription("   approx=0", "Use approximations to compute 'median' or 'quantile' faster.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c neighbourhood", "Applies a statistical operator on a neighbourhood") << std::endl;
   return ss.str();
}
