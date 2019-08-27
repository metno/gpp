#include "Util.h"
#include "Version.h"
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <assert.h>
namespace Cglob {
#include <glob.h>
}
#include <boost/date_time/gregorian/gregorian_types.hpp>
#include <boost/date_time/date_duration.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <execinfo.h>
#include <signal.h>
#include <fstream>
#include <istream>
#include <iomanip>
#include <cstdio>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#ifdef DEBUG
extern "C" void __gcov_flush();
#endif
bool Util::mShowError = true;
bool Util::mShowWarning = false;
bool Util::mShowStatus = false;
bool Util::mShowInfo = false;
float Util::MV = -999;
float Util::pi  = 3.14159265;
double Util::radiusEarth = 6.378137e6;

void Util::error(std::string iMessage) {
#ifdef DEBUG
   if(mShowError) {
      std::cout << "Error: " << iMessage << std::endl;
      void *array[10];
      size_t size = backtrace(array, 10);
      std::cout << "Stack trace:" << std::endl;
      backtrace_symbols_fd(array, size, 2);
   }
   __gcov_flush();
#else
   if(mShowError) {
      std::cout << "Error: " << iMessage << std::endl;
   }
#endif
   abort();
}
void Util::warning(std::string iMessage) {
   if(mShowWarning)
      std::cout << "Warning: " << iMessage << std::endl;
}
void Util::status(std::string iMessage, bool iNewLine) {
   if(mShowStatus) {
      std::cout << "" << iMessage;
      if(iNewLine) {
         std::cout << std::endl;
      }
      else {
         std::fflush(stdout);
      }
   }
}

void Util::info(std::string iMessage) {
   if(mShowInfo)
      std::cout << "Info: " << iMessage << std::endl;
}

double Util::clock() {
   timeval t;
   gettimeofday(&t, NULL);
   double sec = (t.tv_sec);
   double msec= (t.tv_usec);
   return sec + msec/1e6;
}

void Util::setShowError(bool flag) {
   mShowError = flag;
}

void Util::setShowWarning(bool flag) {
   mShowWarning = flag;
}

void Util::setShowStatus(bool flag) {
   mShowStatus = flag;
}

void Util::setShowInfo(bool flag) {
   mShowInfo = flag;
}

bool Util::exists(const std::string& iFilename) {
   std::ifstream infile(iFilename.c_str());
   return infile.good();
}

float Util::getDistance(float lat1, float lon1, float lat2, float lon2, bool approx) {
   if(!Util::isValid(lat1) || !Util::isValid(lat2) ||
      !Util::isValid(lon1) || !Util::isValid(lon2)) {
      return Util::MV;
   }
   if(!(fabs(lat1) <= 90 && fabs(lat2) <= 90 && fabs(lon1) <= 360 && fabs(lon2) <= 360)) {
      std::stringstream ss;
      ss  <<" Cannot calculate distance, invalid lat/lon: (" << lat1 << "," << lon1 << ") (" << lat2 << "," << lon2 << ")";
      Util::error(ss.str());
   }

   if(lat1 == lat2 && lon1 == lon2)
      return 0;

   double lat1r = deg2rad(lat1);
   double lat2r = deg2rad(lat2);
   double lon1r = deg2rad(lon1);
   double lon2r = deg2rad(lon2);
   // Calculate distance according to: http://www.movable-type.co.uk/scripts/latlong.html
   if(approx) {
      float dx2 = pow(cos((lat1r+lat2r)/2),2)*(lon1r-lon2r)*(lon1r-lon2r);
      float dy2 = (lat1r-lat2r)*(lat1r-lat2r);
      return radiusEarth*sqrt(dx2+dy2);
   }
   else {
      double ratio = cos(lat1r)*cos(lon1r)*cos(lat2r)*cos(lon2r)
                   + cos(lat1r)*sin(lon1r)*cos(lat2r)*sin(lon2r)
                   + sin(lat1r)*sin(lat2r);
      double dist = acos(ratio)*radiusEarth;
      return (float) dist;
   }
}

float Util::deg2rad(float deg) {
   return (deg * Util::pi / 180);
}
float Util::rad2deg(float rad) {
   return (rad * 180 / Util::pi);
}
std::vector<std::string> Util::glob(std::string iFilenames) {
   // Split on commas
   std::vector<std::string> returnFiles;
   std::vector<std::string> files = Util::split(iFilenames, ",");
   for(int k = 0; k < files.size(); k++) {
      int flags = GLOB_TILDE | GLOB_NOMAGIC;
      Cglob::glob_t results;
      Cglob::glob(files[k].c_str(), flags, NULL, &results);
      for(int i = 0; i < results.gl_pathc; i++) {
         returnFiles.push_back(results.gl_pathv[i]);
      }
   }
   return returnFiles;
}

std::string Util::gridppVersion() {
   return GRIDPP_VERSION;
}
int Util::getDate(time_t iUnixTime) {
   boost::gregorian::date epoch(1970,1,1);
   boost::gregorian::date_duration diff(iUnixTime/86400);
   boost::gregorian::date newDate = epoch + diff;

   return newDate.year() * 10000 + newDate.month() * 100 + newDate.day();
}
int Util::getTime(time_t iUnixTime) {
   int numSeconds = iUnixTime % 86400;
   int hh = numSeconds / 3600;
   int mm = (numSeconds / 60) % 60;
   int ss = numSeconds % 60;

   return hh * 10000 + mm*100 + ss;
}
time_t Util::getUnixTime(int iDate, int iTime) {
   int year  = iDate / 10000;
   int month = (iDate / 100) % 100;
   int day   = iDate % 100;
   int hh  = iTime / 10000;
   int mm = (iTime / 100) % 100;
   int ss   = iTime % 100;
   int numSeconds = hh*3600 + mm*60 + ss;

   boost::gregorian::date time(year, month, day);
   boost::gregorian::date epoch(1970,1,1);
   boost::gregorian::date_duration diff = time - epoch;
   time_t days = diff.days();
   time_t unixTime = days*86400 + ((time_t) numSeconds);

   return unixTime;
}
int Util::getCurrentDate() {
    boost::gregorian::date today = boost::gregorian::day_clock::universal_day();
    return today.year()*10000 + today.month()*100 + today.day();
}

int Util::getCurrentTime() {
    boost::posix_time::ptime today = boost::posix_time::second_clock::universal_time();
    return today.time_of_day().hours()*10000 + today.time_of_day().minutes()*100 + today.time_of_day().seconds();
}
time_t Util::getCurrentUnixTime() {
   return time(NULL);
}

std::string Util::getCurrentTimeStamp() {
    // boost::gregorian::date today = boost::gregorian::day_clock::universal_day();
    boost::posix_time::ptime today = boost::posix_time::second_clock::universal_time();
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(4) << today.date().year() << "-"
       << std::setfill('0') << std::setw(2) << today.date().month() << "-"
       << std::setfill('0') << std::setw(2) << today.date().day() << " "
       << std::setfill('0') << std::setw(2) << today.time_of_day().hours() << ":"
       << std::setfill('0') << std::setw(2) << today.time_of_day().minutes() << ":"
       << std::setfill('0') << std::setw(2) << today.time_of_day().seconds();
    return ss.str();
}

int Util::calcDate(int iDate, float iAddHours) {
   int year  = floor(iDate / 10000);
   int month = floor(iDate % 10000)/100;
   int day   = floor(iDate % 100);

   int offDay = floor((iAddHours)/24);

   boost::gregorian::date currDate(year, month, day);
   boost::gregorian::date_duration diff(offDay);
   boost::gregorian::date newDate = currDate + diff; 
   int newYear  = newDate.year();
   int newMonth = newDate.month();
   int newDay   = newDate.day();
   int returnDate = newYear * 10000 + newMonth * 100 + newDay;
   return returnDate;
}

std::vector<std::string> Util::split(std::string iString, std::string iDelims) {
   std::vector<std::string> strings;

   // Skip delimiters at beginning.
   std::string::size_type lastPos = iString.find_first_not_of(iDelims, 0);

   // Find first non-delimiter.
   std::string::size_type pos = iString.find_first_of(iDelims, lastPos);

   while(std::string::npos != pos || std::string::npos != lastPos) {
      // Found a string
      strings.push_back(iString.substr(lastPos, pos - lastPos));

      // Skip delimiters
      lastPos = iString.find_first_not_of(iDelims, pos);

      // Find next non-delimiter.
      pos = iString.find_first_of(iDelims, lastPos);
   }

   return strings;
}
float Util::logit(float p) {
   if(!Util::isValid(p) || p <= 0 || p >= 1)
      return Util::MV;
   return log(p/(1-p));
}
float Util::invLogit(float x) {
   if(!Util::isValid(x))
      return Util::MV;
   return exp(x)/(exp(x)+1);
}

bool Util::hasChar(std::string iString, char iChar) {
   return iString.find(iChar) != std::string::npos;
}

bool Util::copy(std::string iFrom, std::string iTo) {
   std::ifstream source(iFrom.c_str(), std::ios::binary);
   std::ofstream dest(iTo.c_str(), std::ios::binary);
   if(!source.good() || !dest.good())
      return false;

   dest << source.rdbuf();

   source.close();
   dest.close();
   return true;
}

bool Util::remove(std::string iFilename) {
   bool failure = std::remove(iFilename.c_str());
   return !failure;
}

std::string Util::formatDescription(std::string iTitle, std::string iMessage, int iTitleLength, int iMaxLength, int iTitleIndent) {
   // Invalid input
   if(iTitleLength >= iMaxLength ) {
      std::stringstream ss;
      ss << "Cannot format description: " << iTitle << ": " << iMessage << std::endl;
      Util::warning(ss.str());
      std::stringstream ss2;
      ss2 << iTitle << " " << iMessage;
      return ss2.str();
   }

   std::vector<std::string> words = Util::split(iMessage);
   std::stringstream ss;
   std::stringstream curr;
   // Indent title
   for(int i = 0; i < iTitleIndent; i++)
      curr << " ";
   curr  << iTitle;
   int N = iTitleLength-curr.str().length();
   // Pad the remainder of the title with spaces
   for(int k = 0; k < N; k++) {
      curr << " ";
   }

   for(int i = 0; i < words.size(); i++) {
      std::string word = words[i];
      int currLength = curr.str().length();
      if(currLength + word.length() > iMaxLength) {
         // Create a new line
         ss << curr.str();
         ss << std::endl;
         curr.str("");
         // Indent to the beginning of the message
         for(int k = 0; k < iTitleLength; k++)
            curr << " ";
      }
      else if(i != 0) {
         curr << " ";
      }
      curr << word;
   }
   ss << curr.str();
   return ss.str();
}

float Util::calculateStat(const std::vector<float>& iArray, Util::StatType iStatType, float iQuantile) {
   // Initialize to missing
   float value = Util::MV;
   if(iStatType == Util::StatTypeMean || iStatType == Util::StatTypeSum) {
      float total = 0;
      int count = 0;
      for(int n = 0; n < iArray.size(); n++) {
         if(Util::isValid(iArray[n])) {
            total += iArray[n];
            count++;
         }
      }
      if(count > 0) {
         if(iStatType == Util::StatTypeMean)
            value = total / count;
         else
            value = total;
      }
   }
   else if(iStatType == Util::StatTypeStd) {
      // STD = sqrt(E[X^2] - E[X]^2)
      // The above formula is unstable when the variance is small and the mean is large.
      // Use the property that VAR(X) = VAR(X-K). Provided K is any element in the array,
      // the resulting calculation of VAR(X-K) is stable. Set K to the first non-missing value.
      float total  = 0;
      float total2 = 0;
      float K = Util::MV;
      int count = 0;
      for(int n = 0; n < iArray.size(); n++) {
         if(Util::isValid(iArray[n])) {
            if(!Util::isValid(K))
               K = iArray[n];
            assert(Util::isValid(K));

            total  += iArray[n] - K;
            total2 += (iArray[n] - K)*(iArray[n] - K);
            count++;
         }
      }
      if(count > 0) {
         float mean  = total / count;
         float mean2 = total2 / count;
         float var   = mean2 - mean*mean;
         if(var < 0) {
            // This should never happen
            var = 0;
            Util::warning("CalibratorNeighbourhood: Problems computing std, unstable result. Setting value to 0");
         }
         float std = sqrt(var);
         value = std;
      }
   }
   else {
      if(iStatType == Util::StatTypeMin)
         iQuantile = 0;
      if(iStatType == Util::StatTypeMedian)
         iQuantile = 0.5;
      if(iStatType == Util::StatTypeMax)
         iQuantile = 1;
      // Remove missing
      std::vector<float> cleanHood;
      cleanHood.reserve(iArray.size());
      for(int i = 0; i < iArray.size(); i++) {
         if(Util::isValid(iArray[i]))
            cleanHood.push_back(iArray[i]);
      }
      int N = cleanHood.size();
      if(N > 0) {
         std::sort(cleanHood.begin(), cleanHood.end());
         int lowerIndex = floor(iQuantile * (N-1));
         int upperIndex = ceil(iQuantile * (N-1));
         float lowerQuantile = (float) lowerIndex / (N-1);
         float upperQuantile = (float) upperIndex / (N-1);
         float lowerValue = cleanHood[lowerIndex];
         float upperValue = cleanHood[upperIndex];
         if(lowerIndex == upperIndex) {
            value = lowerValue;
         }
         else {
            assert(upperQuantile > lowerQuantile);
            assert(iQuantile >= lowerQuantile);
            float f = (iQuantile - lowerQuantile)/(upperQuantile - lowerQuantile);
            assert(f >= 0);
            assert(f <= 1);
            value   = lowerValue + (upperValue - lowerValue) * f;
         }
      }
   }
   return value;
}
bool Util::getStatType(std::string iName, Util::StatType& iType) {
   if(iName == "mean") {
      iType = Util::StatTypeMean;
   }
   else if(iName == "min") {
      iType = Util::StatTypeMin;
   }
   else if(iName == "max") {
      iType = Util::StatTypeMax;
   }
   else if(iName == "median") {
      iType = Util::StatTypeMedian;
   }
   else if(iName == "quantile") {
      iType = Util::StatTypeQuantile;
   }
   else if(iName == "std") {
      iType = Util::StatTypeStd;
   }
   else if(iName == "sum") {
      iType = Util::StatTypeSum;
   }
   else {
      return false;
   }
   return true;
}

vec2 Util::inverse(const vec2 iMatrix) {
   int N = iMatrix.size();
   if(N == 0) {
      return vec2();
   }

   boost::numeric::ublas::matrix<float> inverseMatrix(N,N);

   // Convert to matrix
   boost::numeric::ublas::matrix<float> matrix(N,N);
   for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++) {
         float value = iMatrix[i][j];
         if(!Util::isValid(value)) {
            vec2 matrixVec;
            matrixVec.resize(N);
            // Create a missing matrix
            for(int ii = 0; ii < N; ii++) {
               matrixVec[ii].resize(N, Util::MV);
            }
            return matrixVec;
         }
         assert(iMatrix[i].size() > j);
         matrix(i,j) = iMatrix[i][j];
      }
   }

   // Taken from https://gist.github.com/2464434
   using namespace boost::numeric::ublas;
   typedef permutation_matrix<std::size_t> pmatrix; 

   // create a permutation matrix for the LU-factorization
   pmatrix pm(matrix.size1());

   // perform LU-factorization
   int res = lu_factorize(matrix,pm);

   if( res != 0 ) {
      Util::warning("Could not compute inverse, unstable.");
      vec2 matrixVec;
      matrixVec.resize(N);
      // Create a missing matrix
      for(int ii = 0; ii < N; ii++) {
         matrixVec[ii].resize(N, Util::MV);
      }
      return matrixVec;
   }

   // Initialize to identity matrix
   inverseMatrix.assign(identity_matrix<float>(matrix.size1()));

   // Compute the inverse
   lu_substitute(matrix, pm, inverseMatrix);

   // Convert to vec2
   vec2 inverseVec;
   inverseVec.resize(N);
   for(int i = 0; i < N; i++) {
      inverseVec[i].resize(N);
      for(int j = 0; j < N; j++) {
         inverseVec[i][j] = inverseMatrix(i,j);
      }
   }
   return inverseVec;
}

float Util::interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY) {
   float y = Util::MV;

   if(x > iX[iX.size()-1])
      return iY[iX.size()-1];
   if(x < iX[0])
      return iY[0];

   int i0   = Util::getLowerIndex(x, iX);
   int i1   = Util::getUpperIndex(x, iX);
   assert(Util::isValid(i0));
   assert(i0 >= 0);
   assert(Util::isValid(i1));
   assert(i1 < iX.size());
   float x0 = iX[i0];
   float x1 = iX[i1];
   float y0 = iY[i0];
   float y1 = iY[i1];

   if(x0 == x1)
      y = (y0+y1)/2;
   else {
      assert(x1 >= x0);
      y = y0 + (y1 - y0) * (x - x0)/(x1 - x0);
   }

   return y;

}

std::vector<float> Util::regression(const std::vector<float>& iPredictand, const std::vector<std::vector<float> >& iPredictors, bool iAddConst) {
   int N = iPredictors.size() + iAddConst; // Number of predictors
   int T = iPredictand.size();

   if(T < 2) {
      Util::error("Cannot do regression with only one data point");
   }
   // Check that we have the right number of rows
   for(int i = 0; i < iPredictors.size(); i++) {
      if(iPredictors[i].size() != T) {
         Util::error("One (or more) of the predictors does not have the same number of items as the predictand");
      }
   }

   boost::numeric::ublas::matrix<double> A(T,N);
   boost::numeric::ublas::matrix<double> y(T,1);
   for(int t = 0; t < T; t++) {
      for(int n = 0; n < N; n++) {
         if(iAddConst) {
            if(n == 0)
               A(t, n) = 1;
            else
               A(t, n) = iPredictors[n-1][t];
         }
         else {
            A(t,n) = iPredictors[n][t];
         }
         y(t,0) = iPredictand[t];
      }
   }

#if 0
   std::cout << "\nX" << std::endl;
   for(int i = 0; i < A.size1(); i++) {
      for(int j = 0; j < A.size2(); j++) {
         std::cout << A(i,j) << " ";
      }
      std::cout << std::endl;
   }
#endif

   using namespace boost::numeric::ublas;
   matrix<double> solution(T,1);
   matrix<double> XTX(T,T);
   matrix<double> XT = trans(A);
   XTX = prod(XT, A);

# if 0
   std::cout << "\nXTX" << std::endl;
   for(int i = 0; i < XTX.size1(); i++) {
      for(int j = 0; j < XTX.size2(); j++) {
         std::cout << XTX(i,j) << " ";
      }
      std::cout << std::endl;
   }
#endif
   matrix<double> XTXinv(N,N);
   XTXinv.assign(identity_matrix<double>(N)); 

   typedef permutation_matrix<std::size_t> pmatrix; 
   pmatrix pm(XTX.size1());
   int res = lu_factorize(XTX,pm);
   assert(res == 0);

#if 0
   std::cout << "\npm" << std::endl;
   for(int i = 0; i < XTX.size1(); i++) {
      std::cout << pm(i) << " ";
   }
   std::cout << std::endl;
#endif

   lu_substitute(XTX, pm, XTXinv);
   matrix<double> XTXinvXTy(N,N);
   matrix<double> XTXinvXT = prod(XTXinv,XT);
   solution = prod(XTXinvXT, y);

   //boost::numeric::ublas::identity_matrix<float> Ainv(A.size1());
   //boost::numeric::ublas::permutation_matrix<size_t> pm(A.size1());
   //boost::numeric::ublas::lu_factorize(A, pm);
   //boost::numeric::ublas::lu_substitute(A, pm, y);
   std::vector<float> values(N, 0);
   for(int n = 0; n < N; n++) {
      values[n] = solution(n,0);
   }
   return values;
}

int Util::getLowerIndex(float iX, const std::vector<float>& iValues) {
   int index = Util::MV;
   for(int i = 0; i < (int) iValues.size(); i++) {
      float currValue = iValues[i];
      if(Util::isValid(currValue)) {
         if(currValue < iX) {
            index = i;
         }
         else if(currValue == iX) {
            index = i;
            break;
         }
         else if(currValue > iX) {
            break;
         }
      }
   }
   return index;
}

int Util::getUpperIndex(float iX, const std::vector<float>& iValues) {
   int index = Util::MV;
   for(int i = iValues.size()-1; i >= 0; i--) {
      float currValue = iValues[i];
      if(Util::isValid(currValue)) {
         if(currValue > iX) {
            index = i;
         }
         else if(currValue == iX) {
            index = i;
            break;
         }
         else if(currValue < iX) {
            break;
         }
      }
   }
   return index;
}
