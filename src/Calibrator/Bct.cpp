#include "Bct.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/students_t.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Parameters.h"
CalibratorBct::CalibratorBct(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mMaxEnsMean(100) {
   iOptions.check();
}

bool CalibratorBct::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      Field& field = *iFile.getField(mVariable, t);

      Parameters parametersGlobal;
      if(!iParameterFile->isLocationDependent())
         parametersGlobal = iParameterFile->getParameters(t);

      #pragma omp parallel for
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            Parameters parameters;
            if(iParameterFile->isLocationDependent())
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
            else
               parameters = parametersGlobal;

            // Compute ensemble mean
            float ensMean = 0;
            int counter = 0;
            bool isValid = true;
            for(int e = 0; e < nEns; e++) {
               float value = field(i,j,e);
               if(Util::isValid(value) && Util::isValid(ensMean)) {
                  ensMean += value;
                  counter++;
               }
               else {
                  ensMean = Util::MV;
               }
            }
            if(Util::isValid(ensMean))
               ensMean = ensMean / counter;

            // Calculate ensemble standard deviation
            float ensStd = Util::MV;
            if(Util::isValid(ensMean)) {
               float deviation = 0;
               int counter = 0;
               for(int e = 0; e < nEns; e++) {
                  float value = field(i,j,e);
                  if(Util::isValid(value) && Util::isValid(deviation)) {
                     deviation += pow(value - ensMean, 2);
                     counter++;
                  }
                  else {
                     deviation = Util::MV;
                  }
               }
               if(Util::isValid(deviation))
                  ensStd = sqrt(deviation / counter);
            }

            // Make a copy of the raw values for later sorting
            const std::vector<float>& raw = field(i,j);

            // Only calibrate the ensemble if all members are available. Otherwise
            // use the raw members.
            if(Util::isValid(ensMean) && Util::isValid(ensStd) && counter > 0) {
               // Limit the input to the calibration model to prevent it
               // from creating very extreme values.
               if(Util::isValid(mMaxEnsMean) && ensMean > mMaxEnsMean) {
                  ensMean = mMaxEnsMean;
               }

               // Compute ensemble
               // Calibrate
               std::vector<std::pair<float,int> > pairs(nEns);
               std::vector<float> valuesCal(nEns);
               for(int e = 0; e < nEns; e++) {
                  float quantile = ((float) e+0.5)/nEns;
                  float valueCal = getInvCdf(quantile, ensMean, ensStd, parameters);
                  field(i,j,e) = valueCal;
                  if(!Util::isValid(valueCal)) {
                     isValid = false;
                  }
               }
               if(isValid) {
                  std::vector<float> cal = field(i,j);
                  Calibrator::shuffle(raw, cal);
                  for(int e = 0; e < nEns; e++) {
                     field(i,j,e) = cal[e];
                  }
               }
               else {
                  // Calibrator produced some invalid members. Revert to the raw values.
                  for(int e = 0; e < nEns; e++) {
                     field(i,j,e) = raw[e];
                  }
               }
            }
            else {
               // One or more members are missing, don't calibrate
               for(int e = 0; e < nEns; e++) {
                  field(i,j,e) = raw[e];
               }
            }
         }
      }
   }
   return true;
}

float CalibratorBct::getInvCdf(float iQuantile, float iEnsMean, float iEnsStd, const Parameters& iParameters) {
   if(iQuantile < 0 || iQuantile >= 1) {
      Util::warning("Quantile must be in the interval [0,1)");
      return Util::MV;
   }
   // Check that predictors are valid
   if(!Util::isValid(iEnsMean) || iEnsMean < 0 || !Util::isValid(iEnsStd) || iEnsStd < 0) {
      return Util::MV;
   }

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i])) {
         return Util::MV;
      }
   }

   if(iQuantile == 0)
      return 0;

   // Compute mu
   float a = iParameters[0];
   float b = iParameters[1];
   double mu = a + b * iEnsMean;

   // Compute sigma
   float c = iParameters[2];
   float d = iParameters[3];
   double sigma = exp(c + d * pow(iEnsStd, 1.0/3));

   // Compute nu
   float e = iParameters[4];
   float f = iParameters[5];
   double nu = e + f * iEnsMean;

   // Compute tau
   float g = iParameters[6];
   // boost throws a rounding exception when tau is too large (over e^22).
   // Therefore ensure it is not too large. The t-distribution approaches
   // a normal distribution for large values of tau so this approximation 
   // shouldn't have a big effect.
   if(g > 10)
      g = 10;
   double tau = exp(g);

   boost::math::students_t_distribution<> dist(tau);

   // Compute quantile in truncated-t distribution (Z)
   float qz = Util::MV;
   if(nu <= 0) {
      qz = iQuantile * boost::math::cdf(dist, 1.0 / (sigma * fabs(nu)));
   }
   else {
      qz = 1 - (1 - iQuantile) * boost::math::cdf(dist, 1.0 / (sigma * fabs(nu)));
   }

   // Compute threshold corresponding to this quantile
   float z = boost::math::quantile(dist, qz);

   // Compute threshold in BCT distribution
   float value = Util::MV;
   if(nu != 0) {
      value =  mu * pow(1 + sigma * nu * z, 1.0/nu);
   }
   else {
      value = mu * exp(sigma * z);
   }

   if(!Util::isValid(value))
      return Util::MV;
   return value;
}

std::string CalibratorBct::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c bct", "Calibrates an ensemble using a Box-Cox t-distribution, suitable for parameters like windspeed. The distribution has four parameters:") << std::endl;
      ss << Util::formatDescription("", "* mean  = a + b * ensmean") << std::endl;
      ss << Util::formatDescription("", "* sigma = exp(c + d * ensstd^(1/3)") << std::endl;
      ss << Util::formatDescription("", "* nu = e + f * ensmean") << std::endl;
      ss << Util::formatDescription("", "* tau = exp(g)") << std::endl;
      ss << Util::formatDescription("", "where ensmean is the ensemble mean, and ensstd is the ensemble standard deviation. The parameter set must contain 7 columns with the values [a b c d e f g].") << std::endl;
      ss << Util::formatDescription("   maxEnsMean=100", "Upper limit of what the ensemble mean is allowed to be when passed into the distribution. This effectively prevents the distribution from giving very high values.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c bct", "Calibrates an ensemble using a Box-Cox t-distribution") << std::endl;
   return ss.str();
}
