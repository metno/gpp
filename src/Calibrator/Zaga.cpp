#include "Zaga.h"
#include <algorithm>
#include <math.h>
#include <boost/math/distributions/gamma.hpp>
#include "../Util.h"
#include "../File/File.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Parameters.h"
CalibratorZaga::CalibratorZaga(const Variable& iVariable, const Options& iOptions):
      Calibrator(iVariable, iOptions),
      mFracThreshold(0.5),
      mNeighbourhoodSize(0),
      mPopThreshold(0.5),
      mPrecipLowQuantile(Util::MV),
      mPrecipMiddleQuantile(Util::MV),
      mPrecipHighQuantile(Util::MV),
      mLowVariable(""),
      mMiddleVariable(""),
      mHighVariable(""),
      mMaxEnsMean(100),
      mLogLikelihoodTolerance(1e-4),
      m6h(false) {

   iOptions.getValue("fracThreshold", mFracThreshold);
   iOptions.getValue("precipLowQuantile", mPrecipLowQuantile);
   iOptions.getValue("precipMiddleQuantile", mPrecipMiddleQuantile);
   iOptions.getValue("precipHighQuantile", mPrecipHighQuantile);
   iOptions.getValue("lowVariable", mLowVariable);
   iOptions.getValue("middleVariable", mMiddleVariable);
   iOptions.getValue("highVariable", mHighVariable);
   iOptions.getValue("popVariable", mPopVariable);
   iOptions.getValue("neighbourhoodSize", mNeighbourhoodSize);
   iOptions.getValue("popThreshold", mPopThreshold);
   iOptions.getValue("maxEnsMean", mMaxEnsMean);
   iOptions.getValue("6h", m6h);
   if(!Util::isValid(mMaxEnsMean) || mMaxEnsMean <= 0) {
      Util::error("CalibratorZaga maxEnsMean must be > 0");
   }
   if(mNeighbourhoodSize < 0) {
      std::stringstream ss;
      ss << "CalibratorZaga neighbourhoodSize (" << mNeighbourhoodSize << ") must be >= 0";
      Util::error(ss.str());
   }
   iOptions.check();
}

bool CalibratorZaga::calibrateCore(File& iFile, const ParameterFile* iParameterFile) const {
   int nLat = iFile.getNumY();
   int nLon = iFile.getNumX();
   int nEns = iFile.getNumEns();
   int nTime = iFile.getNumTime();
   vec2 lats = iFile.getLats();
   vec2 lons = iFile.getLons();
   vec2 elevs = iFile.getElevs();

   int startTime = 0;
   int timeWindow = 1;
   if(m6h) {
      startTime = 5;
      timeWindow = 6;
   }

   std::vector<FieldPtr> precips;
   for(int t = 0; t < nTime; t++) {
      precips.push_back(iFile.getField(mVariable, t));
   }

   // Loop over offsets
   for(int t = 0; t < nTime; t++) {
      int numInvalidRaw = 0;
      int numInvalidCal = 0;

      Parameters parametersGlobal;
      if(!iParameterFile->isLocationDependent())
         parametersGlobal = iParameterFile->getParameters(t);

      // Load the POP output field, if needed
      FieldPtr pop;
      if(mPopVariable != "") {
         pop = iFile.getField(mPopVariable, t);
      }
      FieldPtr precipLow;
      if(Util::isValid(mPrecipLowQuantile)) {
         precipLow = iFile.getField(mLowVariable, t);
      }
      FieldPtr precipMiddle;
      if(Util::isValid(mPrecipMiddleQuantile)) {
         precipMiddle = iFile.getField(mMiddleVariable, t);
      }
      FieldPtr precipHigh;
      if(Util::isValid(mPrecipHighQuantile)) {
         precipHigh = iFile.getField(mHighVariable, t);
      }

      #pragma omp parallel for reduction(+:numInvalidRaw, numInvalidCal)
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            Parameters parameters;
            if(iParameterFile->isLocationDependent())
               parameters = iParameterFile->getParameters(t, Location(lats[i][j], lons[i][j], elevs[i][j]));
            else
               parameters = parametersGlobal;

            // for Pop6h, the first few hours are undefined, since we cannot do a 6h accumulation
            if(mPopVariable != "" && t < startTime) {
               for(int e = 0; e < nEns; e++) {
                  (*pop)(i,j,e) = Util::MV;
               }
            }
            else {
               // Compute model variables
               float ensMean = 0;
               float ensFrac = 0;
               int counter = 0;
               bool isValid = true;
               for(int e = 0; e < nEns; e++) {
                  // Create a neighbourhood ensemble
                  for(int ii = std::max(0, i-mNeighbourhoodSize); ii <= std::min(nLat-1, i+mNeighbourhoodSize); ii++) {
                     for(int jj = std::max(0, j-mNeighbourhoodSize); jj <= std::min(nLon-1, j+mNeighbourhoodSize); jj++) {

                        float total = 0;
                        // Accumulate over the time window, if necessary
                        for(int tt = t-timeWindow + 1; tt <= t; tt++) {
                           Field& precip = *precips[tt];

                           float value = precip(ii,jj,e);
                           if(Util::isValid(total) && Util::isValid(value)) {
                              total += value;
                           }
                           else {
                              total = Util::MV;
                           }
                        }
                        if(Util::isValid(total) && Util::isValid(ensMean)) {
                           ensMean += total;
                           ensFrac += (total <= mFracThreshold);
                           counter++;
                        }
                        else {
                           ensMean = Util::MV;
                           ensFrac = Util::MV;
                        }
                     }
                  }
               }
               Field& precip = *precips[t];
               const std::vector<float>& precipRaw = precip(i,j);

               // Only calibrate the ensemble if all members are available. Otherwise
               // use the raw members.
               if(Util::isValid(ensMean) && counter > 0) {
                  ensMean = ensMean / counter;
                  ensFrac = ensFrac / counter;

                  if(i == 0 && j == 0) {
                     // std::cout << ensMean << " " << ensFrac << " " << counter << std::endl;
                  }

                  // Limit the input to the calibration model to prevent it
                  // from creating very extreme values.
                  if(Util::isValid(mMaxEnsMean) && ensMean > mMaxEnsMean) {
                     ensMean = mMaxEnsMean;
                  }

                  // Compute POP
                  if(mPopVariable != "") {
                     for(int e = 0; e < nEns; e++) {
                        float cdf = getCdf(mPopThreshold, ensMean, ensFrac, parameters);
                        // std::cout << ensMean << " " << ensFrac << " " << cdf << std::endl;
                        if(Util::isValid(cdf))
                           (*pop)(i,j,e) = 1 - cdf;
                        else
                           (*pop)(i,j,e) = Util::MV;
                     }
                     if(Util::isValid(mPrecipLowQuantile)) {
                        for(int e = 0; e < nEns; e++) {
                           (*precipLow)(i,j,e) = getInvCdf(mPrecipLowQuantile, ensMean, ensFrac, parameters);
                        }
                     }
                     if(Util::isValid(mPrecipMiddleQuantile)) {
                        for(int e = 0; e < nEns; e++) {
                           (*precipMiddle)(i,j,e) = getInvCdf(mPrecipMiddleQuantile, ensMean, ensFrac, parameters);
                        }
                     }
                     if(Util::isValid(mPrecipHighQuantile)) {
                        for(int e = 0; e < nEns; e++) {
                           (*precipHigh)(i,j,e) = getInvCdf(mPrecipHighQuantile, ensMean, ensFrac, parameters);
                        }
                     }
                  }
                  // Compute precip ensemble
                  else {
                     // Calibrate
                     std::vector<std::pair<float,int> > pairs(nEns);
                     std::vector<float> valuesCal(nEns);
                     for(int e = 0; e < nEns; e++) {
                        float quantile = ((float) e+0.5)/nEns;
                        float valueCal   = getInvCdf(quantile, ensMean, ensFrac, parameters);
                        precip(i,j,e) = valueCal;
                        if(!Util::isValid(valueCal))
                           isValid = false;
                     }
                     if(isValid) {
                        std::vector<float> precipCal = precip(i,j);
                        Calibrator::shuffle(precipRaw, precipCal);
                        for(int e = 0; e < nEns; e++) {
                           precip(i,j,e) = precipCal[e];
                        }
                     }
                     else {
                        // Calibrator produced some invalid members. Revert to the raw values.
                        for(int e = 0; e < nEns; e++) {
                           precip(i,j,e) = precipRaw[e];
                        }
                     }
                  }
               }
               else {
                  // One or more members are missing, don't calibrate
                  for(int e = 0; e < nEns; e++) {
                     precip(i,j,e) = precipRaw[e];
                  }
               }
            }
         }
      }
      if(numInvalidRaw > 0) {
         std::stringstream ss;
         ss << "File '" << iFile.getFilename() << "' missing " << numInvalidRaw
            << "/" << nLat * nLon << " ensembles for timestep " << t << ".";
         Util::warning(ss.str());
      }
      if(numInvalidCal > 0) {
         std::stringstream ss;
         ss << "Calibrator produces '" << numInvalidRaw
            << "/" << nLat * nLon << " invalid ensembles for timestep " << t << ".";
         Util::warning(ss.str());
      }
   }
   return true;
}

float CalibratorZaga::getInvCdf(float iQuantile, float iEnsMean, float iEnsFrac, const Parameters& iParameters) {
   if(iQuantile < 0 || iQuantile >= 1) {
      Util::warning("Quantile must be in the interval [0,1)");
      return Util::MV;
   }
   if(!Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac))
      return Util::MV;

   if(iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1)
      return Util::MV;

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }

   if(iQuantile == 0)
      return 0;

   // Check if we are in the discrete mass
   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);
   if(!Util::isValid(P0))
      return Util::MV;
   if(iQuantile < P0)
      return 0;

   float mua = iParameters[0];
   float mub = iParameters[1];
   float sa  = iParameters[2];
   float sb  = iParameters[3];

   float quantileCont = (iQuantile-P0)/(1-P0);
   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3));
   float sigma = exp(sa + sb * iEnsMean);

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
   boost::math::gamma_distribution<> dist(shape, scale);
   float value = boost::math::quantile(dist, quantileCont);
   if(!Util::isValid(value))
      return Util::MV;
   return value;
}
float CalibratorZaga::getCdf(float iThreshold, float iEnsMean, float iEnsFrac, const Parameters& iParameters) {
   if(!Util::isValid(iThreshold) || !Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac))
      return Util::MV;

   if(iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1)
      return Util::MV;

   if(iThreshold < 0)
      return 0;

   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }

   // Check if we are in the discrete mass
   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);
   if(!Util::isValid(P0))
      return Util::MV;
   if(iThreshold == 0)
      return P0;

   float mua = iParameters[0];
   float mub = iParameters[1];
   float sa  = iParameters[2];
   float sb  = iParameters[3];

   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3));
   float sigma = exp(sa + sb * iEnsMean);

   // std::cout << mua << " " << mub << " " << sa << " " << sb << " " << mu << " " << sigma << std::endl;

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
   boost::math::gamma_distribution<> dist(shape, scale);
   float contCdf = boost::math::cdf(dist, iThreshold) ;
   float cdf = P0 + (1 - P0)*contCdf;
   // std::cout << cdf << " " << P0 << " " << iThreshold << " " << contCdf << " " << shape << " " << scale << std::endl;
   if(!Util::isValid(cdf))
      return Util::MV;
   assert(cdf <= 1);
   assert(cdf >= 0);
   return cdf;
}

float CalibratorZaga::getPdf(float iThreshold, float iEnsMean, float iEnsFrac, const Parameters& iParameters) {
   float mua = iParameters[0];
   float mub = iParameters[1];
   float sa  = iParameters[2];
   float sb  = iParameters[3];

   float P0 = getP0(iEnsMean, iEnsFrac, iParameters);

   // Compute parameters of distribution (in same way as done in gamlss in R)
   float mu    = exp(mua + mub * pow(iEnsMean, 1.0/3));
   float sigma = exp(sa + sb * iEnsMean);

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
   boost::math::gamma_distribution<> dist(shape, scale);
   float contPdf = boost::math::pdf(dist, iThreshold) ;
   float pdf = (1 - P0)*contPdf;
   return pdf;
}
float CalibratorZaga::getP0(float iEnsMean, float iEnsFrac, const Parameters& iParameters) {
   if(!Util::isValid(iEnsMean) || !Util::isValid(iEnsFrac) || iEnsMean < 0 || iEnsFrac < 0 || iEnsFrac > 1)
      return Util::MV;
   // Check that parameters are valid
   for(int i =0; i < iParameters.size(); i++) {
      if(!Util::isValid(iParameters[i]))
         return Util::MV;
   }
   float a = iParameters[4];
   float b = iParameters[5];
   float c = iParameters[6];
   float d = iParameters[7];
   float logit = a + b * iEnsMean + c * iEnsFrac + d * pow(iEnsMean, 1.0/3);
   float P0 = Util::invLogit(logit);
   return P0;
}

Parameters CalibratorZaga::train(const std::vector<ObsEns>& iData) const {
   if(iData.size() == 0) {
      std::cout << "No data to train on...";
      return Parameters();
   }
   std::vector<float> obs, mean, frac;
   obs.resize(iData.size(), Util::MV);
   mean.resize(iData.size(), Util::MV);
   frac.resize(iData.size(), Util::MV);
   // Compute predictors in model
   for(int i = 0; i < iData.size(); i++) {
      obs[i] = iData[i].first;
      std::vector<float> ens = iData[i].second;
      mean[i] = Util::calculateStat(ens, Util::StatTypeMean);
      int total = 0;
      int valid = 0;
      for(int j = 0; j < ens.size(); j++) {
         if(Util::isValid(ens[j])) {
            total += ens[j] <= mFracThreshold;
            valid++;
         }
      }
      if(valid > 0)
         frac[i] = total / valid;
      else
         abort();
      assert(frac[i] >= 0 && frac[i] <= 1);
   }

   int N = mean.size();
   double* p = new double[1+3*N]; 
   p[0] = N;
   for(int n = 0; n < N; n++) {
      p[1+n] = obs[n];
      p[1+n+N] = mean[n];
      p[1+n+2*N] = frac[n];
   }

   /*
   gsl_multimin_function_fdf my_func;
   double p[8] = { 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 }; 

   my_func.n = 8;
   my_func.f = &CalibratorZaga::my_f;
   my_func.df = &CalibratorZaga::my_df;
   my_func.fdf = &CalibratorZaga::my_fdf;
   my_func.params = (void *)p;

   */

   gsl_multimin_function my_func;
   my_func.n = 8;
   my_func.f = &CalibratorZaga::my_f;
   my_func.params = (void *)p;

   // Initialize parameters
   gsl_vector* x = gsl_vector_alloc (8);
   gsl_vector_set (x, 0, 0.1);
   gsl_vector_set (x, 1, 0.1);
   gsl_vector_set (x, 2, 0.1);
   gsl_vector_set (x, 3, 0.1);
   gsl_vector_set (x, 4, 0.1);
   gsl_vector_set (x, 5, 0.1);
   gsl_vector_set (x, 6, 0.1);
   gsl_vector_set (x, 7, 0.1);

   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   // gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 2);
   // gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

   const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
   gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, 8);
   gsl_vector *ss = gsl_vector_alloc (8);
   gsl_vector_set_all (ss, 0.01);
   gsl_multimin_fminimizer_set (s, &my_func, x, ss);


   int iter = 0;
   int status = GSL_CONTINUE;
   do
   {
      iter++;
      // status = gsl_multimin_fdfminimizer_iterate (s);
      // status = gsl_multimin_fminimizer_iterate (s);
      gsl_multimin_fminimizer_iterate (s);

      // if (status)
      //    break;

      double size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, mLogLikelihoodTolerance);

      // status = gsl_multimin_test_gradient (s->gradient, 1e-3);
      // if (status == GSL_SUCCESS)
      //    printf ("Minimum found at:\n");

      if(0) {
         std::cout << iter;
         for(int i = 0; i < 8; i++) {
         std::cout << " " << gsl_vector_get(s->x, i);
         assert(Util::isValid(gsl_vector_get(s->x, i)));
         }
         std::cout << std::endl;
      }
   }
   while (status == GSL_CONTINUE && iter < 5000);

   std::vector<float> values(8,0);
   for(int i = 0; i < 8; i++) {
      values[i] = gsl_vector_get (s->x, i);
   }

   // gsl_multimin_fdfminimizer_free (s);
   gsl_multimin_fminimizer_free (s);
   gsl_vector_free (x);
   gsl_vector_free (ss);

   Parameters par(values);

   std::cout << "Iterations: " << iter << std::endl;
   return par;
}

double CalibratorZaga::my_f(const gsl_vector *v, void *params) {
   double x, y;
   double *p = (double *)params;

   int N = p[0];
   double* obs = p + 1;
   double* mean = p + 1 + N;
   double* frac = p + 1 + 2*N;

   const double* arr = gsl_vector_const_ptr(v, 0);
   std::vector<float> vec(arr, arr+8);
   Parameters par(vec);
   float mua = gsl_vector_get(v, 0);
   float mub = gsl_vector_get(v, 1);
   float sa  = gsl_vector_get(v, 2);
   float sb  = gsl_vector_get(v, 3);
   float a   = gsl_vector_get(v, 4);
   float b   = gsl_vector_get(v, 5);
   float c   = gsl_vector_get(v, 6);
   float d   = gsl_vector_get(v, 7);

   float total = 0;
   for(int n = 0; n < N; n++) {
      // float mean03 = pow(mean[n], 1.0/3);
      // float logit = a + b * mean[n] + c * frac[n] + d * mean03;
      // float P0 = Util::invLogit(logit);
      float P0 = getP0(mean[n], frac[n], par);
      if(obs[n] == 0) {
         assert(P0 > 0);
         total += log(P0);
      }
      else {
         /*
         // Compute parameters of distribution (in same way as done in gamlss in R)
         float mu    = exp(mua + mub * mean03);
         float sigma = exp(sa + sb * mean[n]);

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
         boost::math::gamma_distribution<> dist(shape, scale);
         float contPdf = boost::math::pdf(dist, obs[n]) ;
         float pdf = (1 - P0)*contPdf;
         */
         float pdf = getPdf(obs[n], mean[n], frac[n], par);
         if(pdf == 0 || !Util::isValid(pdf))
            pdf = 0.00001;
         assert(pdf > 0);
         total += log(pdf);
      }
   }
   // std::cout << "Log likelihood: " << total << std::endl;
   return -total; 
}
/*
void CalibratorZaga::my_df(const gsl_vector *v, void *params, gsl_vector *df) {
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

void CalibratorZaga::my_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f = my_f(x, params); 
  my_df(x, params, df);
}
*/

float CalibratorZaga::logLikelihood(float obs, float iEnsMean, float iEnsFrac, const Parameters& iParameters) {
   assert(Util::isValid(obs));
   float pdf = getCdf(obs, iEnsMean, iEnsFrac, iParameters);
   return log(pdf);
}

std::string CalibratorZaga::description(bool full) {
   std::stringstream ss;
   if(full) {
      ss << Util::formatDescription("-c zaga", "Calibrates an ensemble using a zero-adjusted gamma distribution, suitable for parameters like precipitation. The distribution has three parameters:") << std::endl;
      ss << Util::formatDescription("", "* mean  = exp(a + b * ensmean^(1/3)") << std::endl;
      ss << Util::formatDescription("", "* sigma = exp(c + d * ensmean") << std::endl;
      ss << Util::formatDescription("", "* logit(p0) = e + f * ensmean + g * ensfrac + h * ensmean^(1/3)") << std::endl;
      ss << Util::formatDescription("", "where ensmean is the ensemble mean, and ensfrac is the fraction of members above a certain threshold (use fracThreshold option). The parameter set must contain 8 columns with the values [a b c d e f g h].") << std::endl;
      ss << Util::formatDescription("   fracThreshold=0.5", "Threshold defining precip/no-precip boundary when computing fraction of members with precip.") << std::endl;
      ss << Util::formatDescription("   neighbourhoodSize=0", "Increase the ensemble by taking all gridpoints within a neighbourhood. A value of 0 means no neighbourhood is used.") << std::endl;
      ss << Util::formatDescription("   precipLowQuantile=undef", "If set, write values to the PrecipLow variable using this quantile (number between 0 and 1)") << std::endl;
      ss << Util::formatDescription("   precipMiddleQuantile=undef", "If set, write values to the PrecipMiddle variable using this quantile (number between 0 and 1)") << std::endl;
      ss << Util::formatDescription("   precipHighQuantile=undef", "If set, write values to the PrecipHigh variable using this quantile (number between 0 and 1)") << std::endl;
      ss << Util::formatDescription("   lowVariable=undef", "If set, write precip low values to this variable") << std::endl;
      ss << Util::formatDescription("   middleVariable=undef", "If set, write precip middle values to this variable name.") << std::endl;
      ss << Util::formatDescription("   highVariable=undef", "If set, write precip high values to this variablei name.") << std::endl;
      ss << Util::formatDescription("   popVariable=undef", "iF set, write POP values to this variable name.") << std::endl;
      ss << Util::formatDescription("   popThreshold=0.5", "If POP is written, what threshold should be used?") << std::endl;
      ss << Util::formatDescription("   maxEnsMean=100", "Upper limit of what the ensemble mean is allowed to be when passed into the distribution. This effectively prevents the distribution to yield very high values.") << std::endl;
      ss << Util::formatDescription("   6h=0", "If POP is produced, should it be based on the precip in the last 6 hours? If so, the Pop6h variable is written.") << std::endl;
   }
   else
      ss << Util::formatDescription("-c zaga", "Calibrates an ensemble using a zero-adjusted gamma distribution") << std::endl;
   return ss.str();
}
