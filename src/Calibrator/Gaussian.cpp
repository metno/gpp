#include "Gaussian.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/normal.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Parameters.h"
CalibratorGaussian::CalibratorGaussian(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mNeighbourhoodSize(0),
      mLogLikelihoodTolerance(1e-5) {

   iOptions.getValue("neighbourhoodSize", mNeighbourhoodSize);
   if(mNeighbourhoodSize < 0) {
      std::stringstream ss;
      ss << "CalibratorGaussian: neighbourhoodSize (" << mNeighbourhoodSize << ") must be >= 0";
      Util::error(ss.str());
   }
   iOptions.check();
}

bool CalibratorGaussian::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
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

            // Compute model variables
            float total2 = 0;
            float total = 0;
            int counter = 0;
            bool isValid = true;
            for(int e = 0; e < nEns; e++) {
               // Create a neighbourhood ensemble
               for(int ii = std::max(0, i-mNeighbourhoodSize); ii <= std::min(nLat-1, i+mNeighbourhoodSize); ii++) {
                  for(int jj = std::max(0, j-mNeighbourhoodSize); jj <= std::min(nLon-1, j+mNeighbourhoodSize); jj++) {
                     float value = field(ii,jj,e);
                     if(Util::isValid(value)) {
                        total += value;
                        total2 += value*value;
                        counter++;
                     }
                  }
               }
            }
            const std::vector<float>& raw = field(i,j);

            // Only calibrate the ensemble if all members are available. Otherwise
            // use the raw members.
            if(counter > 0) {
               float ensMean = Util::MV;
               float ensSpread = Util::MV;
               if(counter > 0) {
                  ensMean = total / counter;
                  ensSpread = sqrt(total2/counter - (ensMean*ensMean));
               }

               // Calibrate
               std::vector<std::pair<float,int> > pairs(nEns);
               std::vector<float> valuesCal(nEns);
               for(int e = 0; e < nEns; e++) {
                  float quantile = ((float) e+0.5)/nEns;
                  float valueCal   = getInvCdf(quantile, ensMean, ensSpread, parameters);
                  field(i,j,e) = valueCal;
                  if(!Util::isValid(valueCal))
                     isValid = false;
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

float CalibratorGaussian::getInvCdf(float iQuantile, float iEnsMean, float iEnsSpread, const Parameters& iParameters) {
   if(iQuantile < 0 || iQuantile >= 1) {
      Util::warning("Quantile must be in the interval [0,1)");
      return Util::MV;
   }
   if(!Util::isValid(iEnsMean) || !Util::isValid(iEnsSpread))
      return Util::MV;

   if(iEnsMean < 0 || iEnsSpread < 0)
      return Util::MV;

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }

   float sa  = iParameters[0];
   float sb  = iParameters[1];

   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = iEnsMean;
   float sigma = exp(sa  + sb * iEnsSpread);

   if(mu <= 0 || sigma <= 0)
      return Util::MV;
   if(!Util::isValid(mu) || !Util::isValid(sigma))
      return Util::MV;

   boost::math::normal dist(mu, sigma);
   float value = boost::math::quantile(dist, iQuantile);
   if(!Util::isValid(value))
      return Util::MV;
   return value;
}
float CalibratorGaussian::getCdf(float iThreshold, float iEnsMean, float iEnsSpread, const Parameters& iParameters) {
   if(!Util::isValid(iThreshold) || !Util::isValid(iEnsMean) || !Util::isValid(iEnsSpread))
      return Util::MV;

   if(iEnsMean < 0 || iEnsSpread < 0)
      return Util::MV;

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }

   float sa  = iParameters[0];
   float sb  = iParameters[1];

   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = iEnsMean;
   float sigma = exp(sa + sb * iEnsSpread);

   if(mu <= 0 || sigma <= 0)
      return Util::MV;
   if(!Util::isValid(mu) || !Util::isValid(sigma))
      return Util::MV;

   // Parameters in boost and wikipedia
   float shape = 1/(sigma*sigma); // k
   float scale = sigma*sigma*mu;  // theta
   if(!Util::isValid(scale) || !Util::isValid(shape))
      return Util::MV;

   // std::cout << mu << " " << sigma << " " << P0 << " " << shape << " " << scale << std::endl;
   boost::math::normal dist(mu, sigma);
   float cdf = boost::math::cdf(dist, iThreshold) ;
   // std::cout << cdf << " " << P0 << " " << iThreshold << " " << contCdf << " " << shape << " " << scale << std::endl;
   if(!Util::isValid(cdf))
      return Util::MV;
   assert(cdf <= 1);
   assert(cdf >= 0);
   return cdf;
}

float CalibratorGaussian::getPdf(float iThreshold, float iEnsMean, float iEnsSpread, const Parameters& iParameters) {
   float sa  = iParameters[0];
   float sb  = iParameters[1];

   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = iEnsMean;
   float sigma = exp(sa + sb * iEnsSpread);

   if(mu <= 0 || sigma <= 0)
      return Util::MV;
   if(!Util::isValid(mu) || !Util::isValid(sigma))
      return Util::MV;

   boost::math::normal dist(mu, sigma);
   float pdf = boost::math::pdf(dist, iThreshold) ;
   return pdf;
}

Parameters CalibratorGaussian::train(const std::vector<ObsEns>& iData) const {
   if(iData.size() == 0) {
      std::cout << "No data to train on...";
      return Parameters();
   }
   std::vector<float> obs, mean, spread;
   obs.resize(iData.size(), Util::MV);
   mean.resize(iData.size(), Util::MV);
   spread.resize(iData.size(), Util::MV);
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      obs[i] = iData[i].first;
      std::vector<float> ens = iData[i].second;
      mean[i] = Util::calculateStat(ens, Util::StatTypeMean);
      spread[i] = Util::calculateStat(ens, Util::StatTypeStd);
   }

   int N = mean.size();
   double* p = new double[1+3*N]; 
   p[0] = N;
   for(int n = 0; n < N; n++) {
      p[1+n] = obs[n];
      p[1+n+N] = mean[n];
      p[1+n+2*N] = spread[n];
   }

   /*
   gsl_multimin_function_fdf my_func;
   double p[8] = { 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 }; 

   my_func.n = 8;
   my_func.f = &CalibratorGaussian::my_f;
   my_func.df = &CalibratorGaussian::my_df;
   my_func.fdf = &CalibratorGaussian::my_fdf;
   my_func.params = (void *)p;

   */

   gsl_multimin_function my_func;
   my_func.n = mNumParameters;
   my_func.f = &CalibratorGaussian::my_f;
   my_func.params = (void *)p;

   // Initialize parameters
   gsl_vector* x = gsl_vector_alloc (mNumParameters);
   gsl_vector_set (x, 0, 0);
   gsl_vector_set (x, 1, 1);

   const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
   gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, mNumParameters);
   gsl_vector *ss = gsl_vector_alloc (mNumParameters);
   gsl_vector_set_all (ss, 0.01);
   gsl_multimin_fminimizer_set (s, &my_func, x, ss);

   int iter = 0;
   int status = GSL_CONTINUE;
   do
   {
      iter++;
      gsl_multimin_fminimizer_iterate (s);

      double size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, mLogLikelihoodTolerance);
      /*
      for(int i = 0; i < 4; i++) {
         std::cout << gsl_vector_get (s->x, i) << " ";
      }
      std::cout << std::endl;
      */

   }
   while (status == GSL_CONTINUE && iter < 5000);

   std::vector<float> values(mNumParameters,0);
   for(int i = 0; i < mNumParameters; i++) {
      values[i] = gsl_vector_get (s->x, i);
   }

   // gsl_multimin_fdfminimizer_free (s);
   gsl_multimin_fminimizer_free (s);
   gsl_vector_free (x);
   gsl_vector_free (ss);

   Parameters par(values);

   std::cout << "Num of iterations: " << iter << std::endl;
   return par;
}

double CalibratorGaussian::my_f(const gsl_vector *v, void *params) {
   double x, y;
   double *p = (double *)params;

   int N = p[0];
   double* obs = p + 1;
   double* mean = p + 1 + N;
   double* spread = p + 1 + 2*N;

   const double* arr = gsl_vector_const_ptr(v, 0);
   std::vector<float> vec(arr, arr+mNumParameters);
   Parameters par(vec);
   float sa  = gsl_vector_get(v, 0);
   float sb  = gsl_vector_get(v, 1);

   float total = 0;
   for(int n = 0; n < N; n++) {
      float pdf = getPdf(obs[n], mean[n], spread[n], par);
      // std::cout << "   " << obs[n] << " " << mean[n] << std::endl;
      if(pdf == 0 || !Util::isValid(pdf))
         pdf = 0.00001;
      assert(pdf > 0);
      total += log(pdf);
   }
   // std::cout << "Log likelihood: " << total << std::endl;
   return -total; 
}
/*
void CalibratorGaussian::my_df(const gsl_vector *v, void *params, gsl_vector *df) {
   double x, y;
   double *p = (double *)params;

   int N = p[0];
   double* obs = p + 1;
   double* mean = p + 1 + N;
   double* frac = p + 1 + 2*N;

   float mua = gsl_vector_get(v, 0);
   float mub = gsl_vector_get(v, 1);
   float sa  = gsl_vector_get(v, 2);
   float sb  = gsl_vector_get(v, 3);
   float a   = gsl_vector_get(v, 4);
   float b   = gsl_vector_get(v, 5);
   float c   = gsl_vector_get(v, 6);
   float d   = gsl_vector_get(v, 7);

   float total = 0;
   double deriv[8] = {0,0,0,0,0,0,0,0};
   int counter = 0;
   for(int n = 0; n < N; n++) {
      float logit = a + b * mean[n] + c * frac[n] + d * pow(mean[n], 1.0/3);
      float P0 = Util::invLogit(logit);
      if(obs[n] == 0) {
         deriv[4] += (1-P0);
         deriv[5] += mean[n]*(1-P0);
         deriv[6] += frac[n]*(1-P0);
         deriv[7] += pow(mean[n], 1.0/3)*(1-P0);
      }
      else {
         // Compute parameters of distribution (in same way as done in gamlss in R)
         float mu    = exp(mua + mub * pow(mean[n], 1.0/3));
         float sigma = exp(sa + sb * mean[n]);

         if(mu <= 0 || sigma <= 0)
            abort();
         if(!Util::isValid(mu) || !Util::isValid(sigma))
            abort();

         // Parameters in boost and wikipedia
         float shape = 1/(sigma*sigma); // k
         float scale = sigma*sigma*mu;  // theta
         if(!Util::isValid(scale) || !Util::isValid(shape))
            abort();

         // std::cout << mu << " " << sigma << " " << P0 << " " << shape << " " << scale << std::endl;
         boost::math::gamma_distribution<> dist(shape, scale);
         float contPdf = boost::math::pdf(dist, obs[n]) ;
         float pdf = (1 - P0)*contPdf;
         // TODO:
         deriv[0] += 1;
      }
      counter++;
   }
   gsl_vector_set(df, 0, 1);
   gsl_vector_set(df, 1, 1);
}

void CalibratorGaussian::my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f = my_f(x, params); 
  my_df(x, params, df);
}
*/

std::string CalibratorGaussian::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c gaussian", "Ensemble spread calibration based on a Gaussian distribution. The distribution has two parameters:") << std::endl;
      ss << Util::formatDescription("", "* mean  = ensmean") << std::endl;
      ss << Util::formatDescription("", "* sigma = exp(a + b * ensspread)") << std::endl;
      ss << Util::formatDescription("", "where ensmean is the ensemble mean, and ensspread is the standard deviation of members. The parameter set must contain 8 columns with the values [a b].") << std::endl;
      ss << Util::formatDescription("   neighbourhoodSize=0", "Increase the ensemble by taking all gridpoints within a neighbourhood. A value of 0 means no neighbourhood is used.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c gaussian", "Ensemble spread calibration based on a Gaussian distribution") << std::endl;
   return ss.str();
}
