#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include <iostream>
#include <map>
#include "../Parameters.h"
#include "../Location.h"
#include "../Options.h"
#include "../Scheme.h"
#include "../KDTree.h"

//! Represents a collection of parameters, one set for each location and forecast time
//! Parameters can be missing for some locations/times
//! File can location and/or time independent, meaning that the same parameters are returned for all
//! locations and/or times. This occurs if setParameters is used without a location and/or if only
//! one time is available.
class ParameterFile : public Scheme {
   public:
      ParameterFile(const Options& iOptions, bool iIsNew=false);

      //! Get the parameter valid for specified forecast timestep. This is an index, not an hour.
      //! @param iAllowNearestNeighbour Use the nearest neighbour if the location isn't in the set
      Parameters getParameters(int iTime, const Location& iLocation, bool iAllowNearestNeighbour=true) const;
      //! Only use this if isLocationDependent() is false otherwise an error occurs
      Parameters getParameters(int iTime) const;

      static ParameterFile* getScheme(std::string iName, const Options& iOptions, bool iIsNew=false);
      //! Finds the nearest parameter location with valid data at time iTime. Returns true if a
      //! location is found and location is stored in iNearestLocation. If no locations are available,
      //! false is returned.
      bool getNearestLocation(int iTime, const Location& iLocation, Location& iNearestLocation) const;

      //! Does this file provide different parameters for different locations?
      virtual bool isLocationDependent() const;
      //! Does this file provide different parameters for different lead times?
      bool isTimeDependent() const;
      //! Does this file have the same number of parameters for all locations/lead times?
      virtual bool isFixedSize() const = 0;

      //! Set the parameter valid for specified time
      void setParameters(Parameters iParameters, int iTime, const Location& iLocation);
      void setParameters(Parameters iParameters, int iTime);
      //! After all parameters have been set, this function must be called
      void recomputeTree() const;

      std::vector<Location> getLocations() const;
      virtual std::vector<int> getTimes() const;
      virtual bool isReadable() const = 0;

      // Returns Util::MV if the number are not consistent
      int getNumParameters() const;

      std::string getFilename() const;
      virtual std::string name() const = 0;
      virtual void write() const {};

      static std::string getDescriptions(bool full=true);

      // Return the number of bytes in cache
      long getCacheSize() const;
      Options getOptions() const;
   protected:

      // Store all location-dependent parameters here
      typedef std::map<Location, std::vector<Parameters>, Location::CmpIgnoreElevation > LocationParameters;
      LocationParameters mParameters; // Location, Offset, Parameters
      std::string mFilename;
      void setFilename(std::string iFilename);
      bool mIsNew; // Should this file be created?
      void initializeEmpty(const std::vector<Location>& iLocations, int iNumTimes, int iNumParameters);
      void setIsTimeDependent(bool iFlag);
      void setMaxTime(int iMaxTime);
      int getMaxTime() const;
   private:
      bool mIsTimeDependent;
      int mMaxTime;

      // Storing nearest neighbour information. Create a tree with the locations so that lookup for
      // a location is fast. However, every time a new location is added to mParameters, the tree
      // must be recomputed.
      mutable KDTree mNearestNeighbourTree;
      // Locations in the tree
      mutable std::vector<Location> mLocations;
      Options mOptions;
};
#include "Text.h"
#include "Simple.h"
#include "Netcdf.h"
#endif
