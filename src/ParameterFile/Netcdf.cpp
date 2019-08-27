#include "Netcdf.h"
#include <fstream>
#include <sstream>
#include "../Util.h"
#include "../NetcdfUtil.h"
#include <assert.h>
#include <set>
#include <fstream>
#include <math.h>

ParameterFileNetcdf::ParameterFileNetcdf(const Options& iOptions, bool iIsNew) : ParameterFile(iOptions, iIsNew),
      mDimName("coeff"),
      mVarName("coefficient"),
      mXDimName(""),
      mYDimName(""),
      mInDefineMode(false) {
   iOptions.getValue("dimName", mDimName);
   iOptions.getValue("varName", mVarName);
   iOptions.getValue("xDim", mXDimName);
   iOptions.getValue("yDim", mYDimName);

   if(iIsNew) {
      int status = nc_create(getFilename().c_str(), NC_NETCDF4, &mFile);
      handleNetcdfError(status, "Could not write file");
      mInDefineMode = true; // New files are in define mode when created
      endDefineMode();
      return;
   }

   int status = nc_open(getFilename().c_str(), NC_NOWRITE, &mFile);
   handleNetcdfError(status, "Could not read file");
   mInDefineMode = false;
   endDefineMode();

   int dLon  = getLonDim(mFile);
   int dLat  = getLatDim(mFile);
   int dCoeff = getCoefficientDim(mFile);
   int nLat   = getDimSize(mFile, dLat);
   int nLon   = getDimSize(mFile, dLon);
   int nCoeff = getDimSize(mFile, dCoeff);
   std::stringstream ss;
   ss << "Parameter file has " << nCoeff << " parameters";
   Util::info(ss.str());

   // TODO: Deal with case where lat and lon are one-dimensional variables
   // Get latitudes
   int vLat = Util::MV;
   if(hasVar(mFile, "lat"))
      vLat = getVar(mFile, "lat");
   else if(hasVar(mFile, "latitude"))
      vLat = getVar(mFile, "latitude");
   else
      Util::error("Could not determine latitude variable");
   vec2 lats = getGridValues(mFile, vLat);

   // Get longitudes
   int vLon = Util::MV;
   if(hasVar(mFile, "lon"))
      vLon = getVar(mFile, "lon");
   else if(hasVar(mFile, "longitude"))
      vLon = getVar(mFile, "longitude");
   else
      Util::error("Could not determine longitude variable");
   vec2 lons = getGridValues(mFile, vLon);

   // Get elevations
   vec2 elevs;
   if(hasVar(mFile, "altitude")) {
      int vElev = getVar(mFile, "altitude");
      elevs = getGridValues(mFile, vElev);
   }
   else {
      elevs.resize(nLat);
      for(int i = 0; i < nLat; i++) {
         elevs[i].resize(nLon);
         for(int j = 0; j < nLon; j++) {
            elevs[i][j]= Util::MV;
         }
      }
   }

   int nTime = 1;
   int dTime = Util::MV;
   if(hasVar(mFile, "time") && hasDim(mFile, "time")) {
      dTime = getTimeDim(mFile);
      nTime = getDimSize(mFile, dTime);
   }
   double* times = new double[nTime];
   if(hasVar(mFile, "time") && hasDim(mFile, "time")) {
      int vTime = getVar(mFile, "time");
      status = nc_get_var_double(mFile, vTime, times);
      handleNetcdfError(status, "could not get times");
   }
   else {
      times[0] = 0;
   }

   int var = getVar(mFile, mVarName);
   long totalNumParameters = nLat*nLon;
   totalNumParameters *= nTime*nCoeff;
   assert(totalNumParameters > 0);
   float* values = getNcFloats(mFile, var);

   // Initialize parameters to empty and then fill in later
   std::vector<Location> locations;
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         Location location(lats[i][j], lons[i][j], elevs[i][j]);
         locations.push_back(location);
      }
   }
   initializeEmpty(locations, nTime, nCoeff);
   Util::info("Done initializing empty parameters");

   /*
   Read parameters from file and arrange them in mParameters. This is a bit tricky because we do not
   force a certain ordering of dimensions in the coefficients variable. In addition, the variable
   may or may not have a time dimension. The algorithm is to figure out what order the various
   dimensions are and then loop over the retrieved parameters placing them into the right position
   in mParameters.
   */

   // Which index in list of dimension is each dimension?
   int latDimIndex = Util::MV;
   int lonDimIndex = Util::MV;
   int timeDimIndex = Util::MV;
   int coeffDimIndex = Util::MV;

   int ndims;
   status = nc_inq_varndims(mFile, var, &ndims);
   handleNetcdfError(status, "could not get number of dimensions of coefficients variable");

   // Fetch the dimensions of the coefficients variable
   int dims[ndims];
   status = nc_inq_vardimid(mFile, var, dims);
   handleNetcdfError(status, "could not get dimensions for coefficients variable");

   // Get dimensions sizes
   std::vector<int> sizes;
   for(int d = 0; d < ndims; d++) {
      size_t size;
      int status = nc_inq_dimlen(mFile, dims[d], &size);
      handleNetcdfError(status, "could not get size of a dimension for coefficients variable");
      sizes.push_back(size);
      if(dims[d] == dLat)
         latDimIndex = d;
      else if(dims[d] == dLon)
         lonDimIndex = d;
      else if(dims[d] == dTime)
         timeDimIndex = d;
      else if(dims[d] == dCoeff)
         coeffDimIndex = d;
   }

   // Check that coefficients has all required dimensions (does not need time)
   if(!Util::isValid(latDimIndex))
      Util::error("Coefficients in " + getFilename() + " is missing latitude dimension");
   if(!Util::isValid(lonDimIndex))
      Util::error("Coefficients in " + getFilename() + " is missing longitude dimension");
   if(!Util::isValid(coeffDimIndex))
      Util::error("Coefficients in " + getFilename() + " is missing coefficient dimension");

   // Loop over the parameters placing them into the right position
   Parameters par;
   std::stringstream ss0;
   ss0 << "Parameter sizes (lat, lon, time, coeff): " << nLat << " " << nLon << " " << nTime << " " << nCoeff;
   Util::info(ss0.str());
   std::vector<int> indices(4, 0);
   for(int i = 0; i < lats.size(); i++) {
      indices[latDimIndex] = i;
      for(int j = 0; j < lats[i].size(); j++) {
         indices[lonDimIndex] = j;
         Location location(lats[i][j], lons[i][j], elevs[i][j]);
         for(int t = 0; t < nTime; t++) {
            if(Util::isValid(timeDimIndex))
               indices[timeDimIndex] = t;
            std::vector<float> par(nCoeff, 0);
            for(int c = 0; c < nCoeff; c++) {
               indices[coeffDimIndex] = c;
               int index = getIndex(indices, sizes);
               par[c] = values[index];
            }
            mParameters[location][t] = Parameters(par);
         }
      }
   }
   if(nTime > 1)
      setIsTimeDependent(true);
   setMaxTime(nTime);
   assert(getNumParameters() > 0);

   delete[] values;
   delete[] times;

   endDefineMode();

   recomputeTree();
}

ParameterFileNetcdf::~ParameterFileNetcdf() {
   nc_close(mFile);
}

int ParameterFileNetcdf::getDim(int iFile, std::string iDim) const {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}

int ParameterFileNetcdf::createDim(int iFile, std::string iDim, int iLength) const {
   int id;
   int status = nc_def_dim(iFile, iDim.c_str(), iLength, &id);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "Could not create dimension '" << iDim << "'";
      handleNetcdfError(status, ss.str());
   }
   return id;
}

int ParameterFileNetcdf::getVar(int iFile, std::string iVar) const {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

bool ParameterFileNetcdf::hasDim(int iFile, std::string iDim) {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   return status == NC_NOERR;
}

bool ParameterFileNetcdf::hasVar(int iFile, std::string iVar) {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   return status == NC_NOERR;
}

int ParameterFileNetcdf::getDimSize(int iFile, int iDim) const {
   size_t len;
   int status = nc_inq_dimlen(iFile, iDim, &len);
   return len;
}

std::vector<int> ParameterFileNetcdf::getDims(int iFile, int iVar) const {
   int ndims;
   int status = nc_inq_varndims(iFile, iVar, &ndims);
   handleNetcdfError(status, "could not get number of dimensions of coefficients variable");

   // Fetch the dimensions of the coefficients variable
   int dims[ndims];
   status = nc_inq_vardimid(iFile, iVar, dims);
   handleNetcdfError(status, "could not get dimensions for coefficients variable");
   std::vector<int> dimsVec(dims, dims+ndims);
   return dimsVec;
}

bool ParameterFileNetcdf::isValid(std::string iFilename) {
   if(!Util::exists(iFilename))
      return false;

   int file;
   int status = nc_open(iFilename.c_str(), NC_NOWRITE, &file);
   nc_close(file);
   return status == NC_NOERR;
}

bool ParameterFileNetcdf::isReadable() const {
   return ParameterFileNetcdf::isValid(getFilename());
}

std::string ParameterFileNetcdf::description(bool full) {
   std::stringstream ss;
   ss << Util::formatDescription("-p netcdf", "Parameters stored in a Netcdf file. File must have contain: dimensions time, lat (or latitude or y), lon (or longitude or x), coefficient; variables with dims: time[time], latitude[lat,lon], longitude[lat,lon], coefficients[*]. the coefficient variable must have lat, lon, coefficient dimensions, and optimally time. These can be in any order. If there is no time dimension, then the parameters are used for all times. The number of parameters in a set must be constant and equals the size of the 'coefficient' dimension.") << std::endl;
   if(full) {
      ss << Util::formatDescription("   dimName=coefficient", "What is the name of the dimension representing different coefficients?") << std::endl;
      ss << Util::formatDescription("   varName=coefficients", "What is the name of the variable containing the coefficients?") << std::endl;
      ss << Util::formatDescription("   file=required", "Filename of file.") << std::endl;
   }
   return ss.str();
}
void ParameterFileNetcdf::write() const {
   startDefineMode();
   std::vector<Location> locations = getLocations();
   std::vector<int> times = getTimes();
   if(times.size() == 0 || locations.size() == 0) {
      Util::error("Cannot write parameter file '" + getFilename() + "'. No data to write.");
   }
   int nCoeff = 0;
   for(int i = 0; i < times.size(); i++) {
      Parameters par0 = getParameters(times[i], locations[0], false);
      nCoeff = std::max(nCoeff, par0.size());
   }

   if(nCoeff == 0) {
      Util::error("Cannot write a parameter file where all parameters are missing");
   }

   int nTime = times.size();
   int nLon = 1;
   int nLat = locations.size();

   // Create dimensions
   int dTime = createDim(mFile, "time", nTime);
   int dLon  = createDim(mFile, "lon", nLon);
   int dLat  = createDim(mFile, "lat", nLat);
   int dCoeff = createDim(mFile, "coeff", nCoeff);
   int dims[4]  = {dTime, dLat, dLon, dCoeff};
   int var = Util::MV;
   int status = nc_def_var(mFile, mVarName.c_str(), NC_FLOAT, 4, dims, &var);
   handleNetcdfError(status, "could not define variable");

   // Create lat/lon/elev
   int vLat = Util::MV;
   int vLon = Util::MV;
   int vElev = Util::MV;
   int vTime = Util::MV;
   int dimsGrid[2] = {dLat, dLon};
   status = nc_def_var(mFile, "latitude", NC_FLOAT, 2, dimsGrid, &vLat);
   handleNetcdfError(status, "could not define latitude");
   status = nc_def_var(mFile, "longitude", NC_FLOAT, 2, dimsGrid, &vLon);
   handleNetcdfError(status, "could not define longitude");
   status = nc_def_var(mFile, "altitude", NC_FLOAT, 2, dimsGrid, &vElev);
   handleNetcdfError(status, "could not define altitude");
   status = nc_def_var(mFile, "time", NC_FLOAT, 1, &dTime, &vTime);
   handleNetcdfError(status, "could not define time");

   endDefineMode();

   double s = Util::clock();
   float* lats  = new float[locations.size()];
   float* lons  = new float[locations.size()];
   float* elevs = new float[locations.size()];
   // Write locations
   for(int i = 0; i < locations.size(); i++) {
      lats[i] = locations[i].lat();
      lons[i] = locations[i].lon();
      elevs[i] = locations[i].elev();
   }
   status = nc_put_var_float(mFile, vLat, lats);
   handleNetcdfError(status, "could not write latitude");
   delete[] lats;
   status = nc_put_var_float(mFile, vLon, lons);
   handleNetcdfError(status, "could not write longitude");
   delete[] lons;
   status = nc_put_var_float(mFile, vElev, elevs);
   handleNetcdfError(status, "could not write altitude");
   delete[] elevs;

   // Write times
   int* timesAr = new int[nTime];
   for(int t = 0; t < times.size(); t++) {
      timesAr[t] = times[t];
   }
   size_t countTime = nTime;
   size_t startTime = 0;
   status = nc_put_vara_int(mFile, vTime, &startTime, &countTime, timesAr);
   handleNetcdfError(status, "could not write times");
   delete[] timesAr;

   double e = Util::clock();
   std::cout << "Writing times/locations: " << e - s << std::endl;

   // Write parameters
   float* values = new float[nTime*nCoeff*nLat];
   for(int t = 0; t < times.size(); t++) {
      int time = times[t];
      for(int i = 0; i < locations.size(); i++) {
         Parameters par = getParameters(times[t], locations[i], false);
         if(par.size() == nCoeff) {
            for(int k = 0; k < nCoeff; k++) {
               int index = t * nLat * nCoeff + i * nCoeff + k;
               values[index] = par[k];
            }
         }
         else if(par.size() == 0) {
            for(int k = 0; k < nCoeff; k++) {
               int index = t * nLat * nCoeff + i * nCoeff + k;
               values[index] = Util::MV;
            }
         }
         else {
            Util::error("Wrong parameter size");
         }
      }
   }
   double e2 = Util::clock();
   std::cout << "Rearranging parameters: " << e2 - e << std::endl;
   size_t count[4] = {nTime, nLat, 1, nCoeff};
   size_t start[4] = {0, 0, 0, 0};
   status = nc_put_vara_float(mFile, var, start, count, values);
   handleNetcdfError(status, "could not write data");
   double e3 = Util::clock();
   std::cout << "Writing parameters: " << e3 - e2 << std::endl;
   delete[] values;
}

void ParameterFileNetcdf::handleNetcdfError(int status, std::string message) const {
   if(status != NC_NOERR) {
      std::stringstream ss;
      if(message == "") {
         ss << "Netcdf error when reading/writing " << getFilename() << ". Netcdf error code: " << status << ".";
      }
      else {
         ss << "Netcdf error for file " << getFilename() << ": " << message << ". "
            << "Netcdf error code: " << status << ".";
      }
      Util::error(ss.str());
   }
}

int ParameterFileNetcdf::getLatDim(int iFile) const {
   int dLat;
   if(mYDimName != "")
      dLat = getDim(iFile, mYDimName);
   else if(hasDim(iFile, "lat"))
      dLat = getDim(iFile, "lat");
   else if(hasDim(iFile, "latitude"))
      dLat = getDim(iFile, "latitude");
   else if(hasDim(iFile, "y"))
      dLat = getDim(iFile, "y");
   else
      Util::error("Could not determine latitude dimension in " + getFilename());
   return dLat;
}

int ParameterFileNetcdf::getLonDim(int iFile) const {
   int dLon;
   if(mXDimName != "")
      dLon = getDim(iFile, mXDimName);
   else if(hasDim(iFile, "lon"))
      dLon = getDim(iFile, "lon");
   else if(hasDim(iFile, "longitude"))
      dLon = getDim(iFile, "longitude");
   else if(hasDim(iFile, "x"))
      dLon = getDim(iFile, "x");
   else
      Util::error("Could not determine longitude dimension in " + getFilename());
   return dLon;
}

int ParameterFileNetcdf::getTimeDim(int iFile) const {
   int dTime;
   if(hasDim(iFile, "time"))
      dTime = getDim(iFile, "time");
   else
      Util::error("Could not determine time dimension in " + getFilename());
   return dTime;
}

int ParameterFileNetcdf::getCoefficientDim(int iFile) const {
   int dCoeff;
   if(hasDim(iFile, mDimName))
      dCoeff = getDim(iFile,mDimName);
   else
      Util::error("Could not find the coefficient dimension in " + getFilename());
   return dCoeff;
}

std::vector<int> ParameterFileNetcdf::getIndices(int i, const std::vector<int>& iCount) const {
   // The last index changes fastest, the first slowest
   int numDims = iCount.size();
   std::vector<int> indices;
   indices.resize(numDims);
   int sizeSoFar = 1;
   for(int k = numDims-1; k >= 0; k--) {
      int index = floor(i / sizeSoFar);
      index     = index % iCount[k];
      indices[k] = index;
      sizeSoFar *= iCount[k];
   }
   return indices;
}

int ParameterFileNetcdf::getIndex(const std::vector<int>& iIndices, const std::vector<int>& iCount) const {
   // The last index changes fastest, the first slowest
   int index = 0;
   int numDims = iCount.size();
   int sizeSoFar = 1;
   for(int k = numDims-1; k >= 0; k--) {
      int currIndex = iIndices[k];
      index += currIndex * sizeSoFar;
      sizeSoFar *= iCount[k];
   }
   return index;
}

float* ParameterFileNetcdf::getNcFloats(int iFile, int iVar) {
   long size = NetcdfUtil::getTotalSize(iFile, iVar);
   assert(size > 0);
   float* values = new float[size];
   int status = nc_get_var_float(iFile, iVar, values);
   char name[512];
   status = nc_inq_varname(iFile, iVar, name);
   handleNetcdfError(status, "could not get variable name");

   std::stringstream ss;
   ss <<  "could not retrieve values from variable '" << name << "'";
   handleNetcdfError(status, ss.str());

   // Convert missing values
   float MV = NetcdfUtil::getMissingValue(iFile, iVar);
   for(int i = 0; i < size; i++) {
      if(values[i] == MV)
         values[i] = Util::MV;
   }
   return values;
}

vec2 ParameterFileNetcdf::getGridValues(int iFile, int iVar) const {
   // Initialize values
   vec2 grid;
   int dLat = getLatDim(iFile);
   int dLon = getLonDim(iFile);
   int nLat = getDimSize(iFile, dLat);
   int nLon = getDimSize(iFile, dLon);
   grid.resize(nLat);
   for(int i = 0; i < nLat; i++) {
      grid[i].resize(nLon, Util::MV);
   }

   // We have a lat/lon grid, where lat/lons are only provided along the pertinent dimension
   // Values are assumed to be constant across the other dimension.
   std::vector<int> dims = getDims(iFile, iVar);
   int numDims = dims.size();
   if(numDims == 1) {
      int dim = dims[0];
      long size = getDimSize(iFile, dim);
      float* values = new float[size];
      nc_get_var_float(iFile, iVar, values);
      // Latitude variable
      if(dim == getLatDim(iFile)) {
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               grid[i][j] = values[i];
            }
         }
      }
      // Longitude variable
      else if(dim == getLonDim(iFile)) {
         for(int i = 0; i < nLat; i++) {
            for(int j = 0; j < nLon; j++) {
               grid[i][j] = values[j];
            }
         }
      }
      else {
         std::stringstream ss;
         ss << "Missing lat or lon dimension";
         Util::error(ss.str());
      }
      delete[] values;
   }
   // We have a projected grid, where lat and lons are provided for each grid point
   else {
      long count[numDims];
      int size = 1;
      int indexLat = Util::MV;
      int indexLon = Util::MV;
      std::vector<int> dims = getDims(iFile, iVar);
      for(int i = 0; i < numDims; i++) {
         if(dims[i] == getLatDim(iFile)) {
            count[i] = nLat;
            size *= count[i];
            indexLat = i;
         }
         else if(dims[i] == getLonDim(iFile)) {
            count[i] = nLon;
            size *= count[i];
            indexLon = i;
         }
         else {
            count[i] = 1;
         }
      }
      if(!Util::isValid(indexLat) || !Util::isValid(indexLon)) {
         std::stringstream ss;
         ss << "Missing lat and/or lon dimensions";
         Util::error(ss.str());
      }
      float* values = new float[size];
      nc_get_var_float(iFile, iVar, values);
      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            // Latitude dimension is ordered first
            if(indexLat < indexLon) {
               grid[i][j] = values[i*nLon + j];
            }
            // Longitude dimension is ordered first
            else {
               grid[i][j] = values[j*nLat + i];
            }
         }
      }
      delete[] values;
   }
   return grid;
}
void ParameterFileNetcdf::startDefineMode() const {
   if(!mInDefineMode) {
      // Only call redefine if we went into data mode at some point
      int status = ncredef(mFile);
      handleNetcdfError(status, "could not put into define mode");
      mInDefineMode = true;
   }
}
void ParameterFileNetcdf::endDefineMode() const {
   if(mInDefineMode) {
      // Only end define mode if we at some point entered
      int status = ncendef(mFile);
      handleNetcdfError(status, "could not put into data mode");
      mInDefineMode = false;
   }
}
