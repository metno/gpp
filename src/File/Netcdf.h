#ifndef FILE_NETCDF_H
#define FILE_NETCDF_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "File.h"
#include "../Variable.h"
#include "../Options.h"

//! Represents an ensemble Netcdf data file from ECMWF
//! Must have:
//!    Time dim: time
//!    Latitude dim:  lat, latitude, y
//!    Longitude dim: lon, longitude, x
//!    Latitude var:  lat, latitude, y. With the latitude dim, and possibly the lon dim
//!    Longitude var: lon, longitude, x. With the longitude dim, and possibly the lat dim
//!    5 dimensional variables (time, *, ensemble_member, <lat dim>, <lon dim>)
class FileNetcdf : public File {
   public:
      FileNetcdf(std::string iFilename, const Options& iOptions=Options(), bool iReadOnly=false);
      ~FileNetcdf();

      static bool isValid(std::string iFilename, const Options& iOptions);

      //! Add attribute to a variable (overwrite if existing). Variable must exist
      void setAttribute(std::string iVariable, std::string iName, std::string iValue);

      //! Get string attribute belowing to a variable. Variable must exist.
      //! Returns "" if attribute does not exist.
      std::string getAttribute(std::string iVariable, std::string iName);
      std::string getAttribute(int iVar, std::string iName);

      //! Add global attribute to file (overwrite if existing)
      void setGlobalAttribute(std::string iName, std::string iValue);

      //! Add global attribute to file (append to attribute if existing)
      void appendGlobalAttribute(std::string iName, std::string iValue);

      //! Add global attribute to file (prepend to attribute if existing)
      void prependGlobalAttribute(std::string iName, std::string iValue);

      //! Get global string attribute. Returns "" if non-existant.
      std::string getGlobalAttribute(std::string iName);

      static std::string description();
      std::string name() const {return "netcdf";};

   protected:
      float getScale(int iVar) const;
      float getOffset(int iVar) const;
      void writeCore(std::vector<Variable> iVariables, std::string iMessage="");
      FieldPtr getFieldCore(const Variable& iVariable, int iTime) const;
      bool hasVariableCore(const Variable& iVariable) const;

      vec2 getGridValues(int iVariable) const;
      void writeAltitude() const;
      void defineAltitude();

      int mYDim;
      int mXDim;
      int mEnsDim;
      int mTimeDim;
      int mLatVar;
      int mLonVar;
      int mElevVar;
      int mLafVar;
      int mTimeVar;
      std::string mEnsDimName;

      int getDim(std::string iDim) const;
      std::string getDimName(int iDim) const;
      std::vector<int> getDims(int iVar) const;
      int getVar(std::string iVar) const;
      std::string getVarName(int iVar) const;
      int getDimSize(std::string iDim) const;
      int getDimSize(int iDim) const;
      int getNumDims(int iVar) const;
      std::string detectEnsDim() const;
      int detectYDim() const;
      int detectTimeDim() const;
      int detectTimeVar() const;
      int detectXDim() const;
      int detectLatVar() const;
      int detectLonVar() const;
      bool hasVar(std::string iVar) const;
      static bool hasVar(int iFile, std::string iVar);
      int mFile;
      void handleNetcdfError(int status, std::string message="") const;
      void startDefineMode() const;
      void startDataMode() const;
      mutable bool mInDataMode;
      float getMissingValue(int iVar) const;
      void  setMissingValue(int iVar, float iValue) const;
      //! Convert linear index 'i' to vector 'iInidices'. 'iCount' specifies the size of the data
      //! Using row-major ordering (last index varies fastest)
      int getIndex(const std::vector<int>& iCount, const std::vector<int>& iIndices) const;
      void getIndices(int i, const std::vector<int>& iCount, std::vector<int>& iIndices) const;
      void setAttribute(int iVar, std::string iName, std::string iValue);
      void defineTimes();
      void defineEns();
      void defineReferenceTime();
      void defineGlobalAttributes(std::string iMesssage);
      void writeTimes();
      void writeReferenceTime();
      bool hasDim(std::string iDim) const;
      vec2 getLatLonVariable(int iVariable) const;
      static bool hasDim(int iFile, std::string iDim);
      const static int mMaxAttributeLength = 100000000;

};
#endif
