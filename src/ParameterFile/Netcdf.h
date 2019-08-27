#ifndef PARAMETER_FILE_NETCDF_H
#define PARAMETER_FILE_NETCDF_H
#include <iostream>
#include <netcdf.h>
#include <map>
#include "ParameterFile.h"
#include "../Parameters.h"
#include "../Location.h"
#include "../File/File.h"
// TODO: if a location has missing parameters for some hours but not others, then
// the missing values are returned, not the nearest neighbours values.
class ParameterFileNetcdf : public ParameterFile {
   public:
      ParameterFileNetcdf(const Options& iOptions, bool iIsNew=false);
      ~ParameterFileNetcdf();

      bool isFixedSize() const {return true;};

      static bool isValid(std::string iFilename);
      bool isReadable() const;
      bool isLocationDependent() const { return true; };

      static std::string description(bool full=true);
      std::string name() const {return "netcdf";};

      void write() const;
   private:
      float mLocalMV;
      int    getLatDim(int iFile) const;
      int    getLonDim(int iFile) const;
      int    getTimeDim(int iFile) const;
      int    getCoefficientDim(int iFile) const;
      std::vector<int> getDims(int iFile, int iVar) const;

      int    getDim(int iFile, std::string iDim) const;
      int    createDim(int iFile, std::string iDim, int iLength) const;
      int    getVar(int iFile, std::string iVar) const;
      int    getDimSize(int iFile, int iDim) const;
      static bool   hasDim(int iFile, std::string iDim);
      static bool   hasVar(int iFile, std::string iVar);
      std::string mDimName;
      std::string mVarName;
      std::string mXDimName;
      std::string mYDimName;
      void handleNetcdfError(int status, std::string message="") const;
      //! Convert linear index 'i' to vector 'iInidices'. 'iCount' specifies the size of the data
      //! Using row-major ordering (last index varies fastest)
      std::vector<int> getIndices(int i, const std::vector<int>& iCount) const;
      int getIndex(const std::vector<int>& iIndices, const std::vector<int>& iCount) const;

      // Read variable from file, convert missing values
      // User must release memory
      float* getNcFloats(int iFile, int iVar);

      // Reads lat/lon/elev values from the variable, regardless of the order of the x, y dimensions
      // For lat/lon, the variable can have one or two dimensions. If the former, then the values
      // are assumed to be the same across the other dimension.
      vec2 getGridValues(int iFile, int iVariable) const;

      int mFile;

      void startDefineMode() const;
      void endDefineMode() const;
      mutable bool mInDefineMode;
};
#endif
