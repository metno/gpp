#ifndef UTIL_H
#define UTIL_H
#include <string>
#include <vector>
#include <set>
#include <cmath>

typedef std::vector<std::vector<float> > vec2; // Lat, Lon

//! Contains helper functions
class Util {
   public:
      //! Aborts the program with an error message to std::cout
      static void error(std::string iMessage);

      //! Issues a warning message to std::cout
      static void warning(std::string iMessage);

      //! Issues a status message to std::cout. Add new line to the end.
      static void status(std::string iMessage, bool iNewLine=true);

      //! Issues a info message to std::cout
      static void info(std::string iMessage);

      //! True will force calls to Util::error to show error messages.
      //! Calls to errors will cause abort regardless of this setting.
      static void setShowError(bool flag);

      //! True will force calls to Util::warning to show warning messages
      static void setShowWarning(bool flag);

      //! True will force calls to Util::status to show status messages
      static void setShowStatus(bool flag);

      //! True will force calls to Util::info to show status messages
      static void setShowInfo(bool flag);

      //! Missing value indicator
      static float MV;

      //! \brief Checks if a value is valid. Inlined for improved performance on some systems
      //! @return false if iValue equals Util::MV or is +-inf or nan. Returns true otherwise.
      inline static bool isValid(float iValue) {
         return !std::isnan(iValue) && !std::isinf(iValue) && iValue != Util::MV;
      }


      //! Checks if the file exists
      static bool exists(const std::string& iFilename);

      //! \brief Computes the geographical distance between two locations, assuming a spherical earth.
      //! @param lat1 Latitude in degrees of point 1
      //! @param lon1 Longitude in degrees of point 1
      //! @param lat2 Latitude in degrees of point 2
      //! @param lon2 Longitude in degrees of point 2
      //! @param approx Should an approximation be used? If true, use equirectangular approximation
      //!               Computes the result much faster.
      //! @return Distance in meters
      static float getDistance(float lat1, float lon1, float lat2, float lon2, bool approx=false);

      //! Returns the current unix time in seconds
      static double clock();
      
      //! Convert degrees to radians
      static float deg2rad(float deg);
      
      //! Convert radians to degrees
      static float rad2deg(float rad);

      static float pi;

      // Return a vector of filenames from a single string containing wildcards and commas
      static std::vector<std::string> glob(std::string iFilenames);

      static std::string gridppVersion();
      
      //! Radius of the earth in meters
      static double radiusEarth;
      
      //! \brief Convert unix timestamp to YYYYMMDD format
      static int getDate(time_t iUnixTime);

      //! \brief Convert unix timestamp to HHMMSS format
      static int getTime(time_t iUnixTime);
      
      //! \brief Get the current unix time (UTC)
      //! @param iDate in YYYYMMDD format
      //! @param iTime in HHMMSS format (can be negative to subtract time)
      static time_t getUnixTime(int iDate, int iTime);

      //! \brief Get the current date in YYYYMMDD format (UTC)
      static int getCurrentDate();

      //! \brief Get the current date in HHMMSS format (UTC)
      static int getCurrentTime();

      //! \brief Get the current date/time in unix time (UTC)
      static time_t getCurrentUnixTime();

      //! \brief Get the current date in YYYY-MM-DD HH:MM:SS format
      static std::string getCurrentTimeStamp();
      
      //! \brief Calculates new date by adding hours to an existing date
      //! @param iDate Date in format YYYYMMDD
      //! @param iAddHours Number of hours to add
      //! @return New date in format YYYYMMDD
      static int calcDate(int iDate, float iAddHours);
      
      //! \brief Creates a vector of strings by splitting string on each space
      //! Multiple consecutive spaces are treated as one space
      static std::vector<std::string> split(std::string iString, std::string iDelims=" ");

      //! \brief Computes the logit of p
      //! @param p Must lie in the interval (0,1)
      static float logit(float p);

      //! \brief Computes the inverse logit of x
      //! @return x A value on the interval (0,1)
      static float invLogit(float x);

      template <class T> static std::vector<T> combine(const std::vector<T>& i1, const std::vector<T>& i2) {
         std::set<T> allValues(i1.begin(), i1.end());
         for(int i = 0; i < i2.size(); i++) {
            allValues.insert(i2[i]);
         }
         return std::vector<T>(allValues.begin(), allValues.end());
      };
      //! Returns true if iString contains one or more of iChar
      static bool hasChar(std::string iString, char iChar);

      //! Create a formatted string useful for outputing descriptions of methods. Places the title
      //! on the left side and the message on the right, with automatic line wrapping.
      static std::string formatDescription(std::string iTitle, std::string iMessage, int iTitleLength=32, int iMaxLength=100, int iTitleIndent=3);

      enum StatType {
         StatTypeMean      = 0,
         StatTypeMin       = 10,
         StatTypeMedian    = 20,
         StatTypeMax       = 30,
         StatTypeQuantile  = 40,
         StatTypeStd       = 50,
         StatTypeSum       = 60
      };

      enum Operator {
         OperatorAdd       = 0,
         OperatorSubtract  = 10,
         OperatorMultiply  = 20,
         OperatorDivide    = 30,
      };

      //! Applies statistics operator to array. Missing values are ignored.
      static float calculateStat(const std::vector<float>& iArray, Util::StatType iStatType, float iQuantile=Util::MV);
      static bool getStatType(std::string iName, Util::StatType& iType);
      
      //! \brief Comparator class for sorting pairs using the first entry.
      //! Sorts from smallest to largest
      template<class T1, class T2> struct sort_pair_first {
         bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
         }
      };

      //! Computes the matrix inverse
      //! If any elements in iMatrix is missing, then all elements in the inverse will be missing.
      //! If the matrix is not invertable, then all elements in the inverse will be missing.
      static vec2 inverse(const vec2 iMatrix);

      static float interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY);

      //! @param addConst Add a constant coefficient to the predictors. The coefficient will appear
      //! first in the output vector.
      static std::vector<float> regression(const std::vector<float>& iPredictand, const std::vector<std::vector<float> >& iPredictors, bool iAddConst=false);

      // Array operations
      //! Note: iValues must be sorted
      //! Finds index into array pointing to the largest element smaller than or equal to iX.
      //! In case of a tie, return the index of the last tied elements. When the search value equals
      //! multiple values, the first element is returned. 
      static int getLowerIndex(float iX, const std::vector<float>&  iValues);
      //! Note: iValues must be sorted
      //! Finds index into array pointing to the smallest element greater than or equal to iX.
      //! In case of a tie, return the index of the first tied elements. When the search value equals
      //! multiple values, the last element is returned.
      static int getUpperIndex(float iX, const std::vector<float>& iValues);

      //! Copy the file with filename iFrom to filename iTo. Returns true if successful.
      static bool copy(std::string iFrom, std::string iTo);
     
      //! Remove the file with filename. Returns true if successful.
      static bool remove(std::string iFilename);
     
      //! \brief Comparator class for sorting pairs using the second entry.
      //! Sorts from smallest to largest
      template<class T1, class T2> struct sort_pair_second {
         bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.second < right.second;
         }
      };
   private:
      static bool mShowStatus;
      static bool mShowError;
      static bool mShowWarning;
      static bool mShowInfo;
};
#endif
