#include "Netcdf.h"
#include <math.h>
#include <netcdf.h>
#include <assert.h>
#include <stdlib.h>
#include "../Util.h"

FileNetcdf::FileNetcdf(std::string iFilename, const Options& iOptions, bool iReadOnly) : File(iFilename, iOptions),
      mInDataMode(true)
{
   int status = nc_open(getFilename().c_str(), iReadOnly ? NC_NOWRITE: NC_WRITE, &mFile);
   if(status != NC_NOERR) {
      if(iReadOnly)
         Util::error("Could not open NetCDF file in read-only mode" + getFilename());
      else
         Util::error("Could not open NetCDF file in write mode" + getFilename());
   }
   // Get defaults
   std::string latVar, lonVar, timeVar, ensDim, yDim, xDim, timeDim;
   if(!iOptions.getValue("latVar", latVar))
      mLatVar = detectLatVar();
   else
      mLatVar = getVar(latVar);

   if(!iOptions.getValue("lonVar", lonVar))
      mLonVar = detectLonVar();
   else
      mLonVar = getVar(lonVar);

   if(!iOptions.getValue("timeVar", timeVar))
      mTimeVar = detectTimeVar();
   else
      mTimeVar = getVar(timeVar);

   // Dimensions
   if(!iOptions.getValue("ensDim", mEnsDimName))
      mEnsDimName = detectEnsDim();

   if(hasDim(mEnsDimName))
      mEnsDim = getDim(mEnsDimName);
   else
      mEnsDim = Util::MV;

   if(!iOptions.getValue("yDim", yDim))
      mYDim = detectYDim();
   else
      mYDim = getDim(yDim);

   if(!iOptions.getValue("xDim", xDim))
      mXDim = detectXDim();
   else
      mXDim = getDim(xDim);

   if(!iOptions.getValue("timeDim", timeDim))
      mTimeDim = detectTimeDim();
   else
      mTimeDim = getDim(timeDim);

   // Determine dimension sizes
   mNEns = 1;
   if(Util::isValid(mEnsDim))
      mNEns = getDimSize(mEnsDim);

   // Compute elevations
   std::string elevVar;
   mElevVar = Util::MV;
   float elevScale = 1;
   if(!iOptions.getValue("elevVar", elevVar)) {
      if(hasVar("altitude")) {
         mElevVar = getVar("altitude");
      }
      else if(hasVar("surface_geopotential")) {
         mElevVar = getVar("surface_geopotential");
         elevScale = 1.0 / 9.81;
      }
   }
   else {
      if(hasVar(elevVar))
         mElevVar = getVar(elevVar);
      else {
         std::stringstream ss;
         ss << "Cannot find altitude field '" << elevVar << "'. Altitudes not used.";
         Util::warning(ss.str());
      }

   }

   std::string lafVar;
   if(!iOptions.getValue("lafVar", lafVar)) {
      lafVar = "land_area_fraction";
   }
   iOptions.check();

   // Done reading options

   // Retrieve lat/lon grid
   bool successLats = setLats(getLatLonVariable(mLatVar));
   if(!successLats) {
      std::stringstream ss;
      ss << "Could not set latitudes in " << getFilename();
      Util::error(ss.str());
   }
   bool successLons = setLons(getLatLonVariable(mLonVar));
   if(!successLons) {
      std::stringstream ss;
      ss << "Could not set longitudes in " << getFilename();
      Util::error(ss.str());
   }

   if(Util::isValid(mElevVar)) {
      vec2 elevs = getLatLonVariable(mElevVar);
      for(int i = 0; i < getNumY(); i++) {
         for(int j = 0; j < getNumX(); j++) {
            elevs[i][j] *= elevScale;
         }
      }
      setElevs(elevs);
   }
   else {
      vec2 elevs;
      elevs.resize(getNumY());
      for(int i = 0; i < getNumY(); i++) {
         elevs[i].resize(getNumX());
         for(int j = 0; j < getNumX(); j++) {
            elevs[i][j] = Util::MV;
         }
      }
      // setElevs(elevs);
      Util::warning("No altitude field available in " + getFilename());
   }

   if(hasVar(lafVar)) {
      int land_area_fraction = getVar(lafVar);
      mLandFractions = getLatLonVariable(land_area_fraction);
   }
   else {
      mLandFractions.resize(getNumY());
      for(int i = 0; i < getNumY(); i++) {
         mLandFractions[i].resize(getNumX());
         for(int j = 0; j < getNumX(); j++) {
            mLandFractions[i][j] = Util::MV;
         }
      }
   }

   // Compute times
   if(hasVar("forecast_reference_time")) {
      double referenceTime = getReferenceTime();
      int vReferenceTime = getVar("forecast_reference_time");
      int numDims = getNumDims(vReferenceTime);
      if(numDims != 0) {
         std::stringstream ss;
         ss << "forecast_reference_time should not have any dimensions";
         Util::error(ss.str());
      }
      int status = nc_get_var_double(mFile, vReferenceTime, &referenceTime);
      handleNetcdfError(status, "could not get reference time");
      setReferenceTime(referenceTime);
   }

   if(Util::isValid(mTimeVar)) {
      int size = getDimSize(mTimeDim);
      if(Util::isValid(size)) {
         double* times = new double[size];
         int status = nc_get_var_double(mFile, mTimeVar, times);
         handleNetcdfError(status, "could not get times");
         setTimes(std::vector<double>(times, times+size));
         delete[] times;
      }
      else {
         // Time is a scalar
         double time = Util::MV;
         int status = nc_get_var_double(mFile, mTimeVar, &time);
         handleNetcdfError(status, "could not get time");
         std::vector<double> times(1, time);
         setTimes(times);
      }
   }
   else if(Util::isValid(getReferenceTime())) {
      std::stringstream ss;
      ss << "File does not contain times, using forecast reference time as the output time";
      Util::warning(ss.str());
      std::vector<double> times(1, getReferenceTime());
      setTimes(times);
   }
   else {
      std::stringstream ss;
      ss << "Cannot determine forecast times, since neither time variable nor reference time variable is available";
      Util::error(ss.str());
   }

   // Load variables
   int numVars = Util::MV;
   status = nc_inq_nvars(mFile, &numVars);
   int vars[numVars];
   status = nc_inq_varids(mFile, &numVars, vars);
   handleNetcdfError(status, "could not get list of variables");
   for(int v = 0; v < numVars; v++) {
      char nameChar[10000];
      status = nc_inq_varname (mFile, vars[v], nameChar);
      std::string name(nameChar);
      std::string units = getAttribute(vars[v], "units");
      std::string standardName = getAttribute(vars[v], "standard_name");

      Variable variable(name, units, standardName);
      mVariables.push_back(variable);
   }
}

FileNetcdf::~FileNetcdf() {
   nc_close(mFile);
}

FieldPtr FileNetcdf::getFieldCore(const Variable& iVariable, int iTime) const {
   startDataMode();
   std::string variableName = iVariable.name();
   int var = getVar(variableName);
   std::vector<int> dims = getDims(var);

   // Determine which slice to retrieve
   size_t count[dims.size()];
   size_t start[dims.size()];
   int ensPos = Util::MV;
   int yPos = Util::MV;
   int xPos = Util::MV;
   for(int d = 0; d < dims.size(); d++) {
      int dim = dims[d];
      count[d] = 1;
      start[d] = 0;
      if(dim == mTimeDim) {
         start[d] = iTime;
      }
      else if(dim == mEnsDim) {
         count[d] = getDimSize(dim);
         ensPos = d;
      }
      else if(dim == mYDim) {
         count[d] = getDimSize(dim);
         yPos = d;
      }
      else if(dim == mXDim) {
         count[d] = getDimSize(dim);
         xPos = d;
      }
      else {
         int level = iVariable.level();
         // Non-recognized dimension
         if(Util::isValid(level)) {
            int size = getDimSize(dim);
            if(size > 1) {
               if(level >= dim) {
                  std::stringstream ss;
                  char name[1000];
                  int status = nc_inq_varname(mFile, dim, name);
                  handleNetcdfError(status, "could not get dimension name of level dimension'");
                  ss << "Could not get level " << level
                     << " from dimension '" << name << "' since it only has a size of " << size;
                  Util::error(ss.str());
               }
               start[d] = level;
            }
         }
      }
   }
   if(!Util::isValid(ensPos) && Util::isValid(mEnsDim)) {
      std::stringstream ss;
      ss << "Ensemble dimension '" << getDimName(mEnsDim) << "' found in file, but variable '" << iVariable.name() << "' does not have this dimension";
      Util::warning(ss.str());
   }

   // Initialize vector
   size_t size = 1;
   for(int d = 0; d < dims.size(); d++) {
      size *= count[d];
   }
   float* values = new float[size];
   nc_get_vara_float(mFile, var, start, count, values);

   float MV = getMissingValue(var);
   float offset = getOffset(var);
   float scale = getScale(var);

   FieldPtr field = getEmptyField();
   std::vector<int> countVector(count, count + dims.size());
   int nEns = 1;
   if(Util::isValid(ensPos))
      nEns = count[ensPos];
   int nY = 1;
   if(Util::isValid(yPos))
      nY = count[yPos];
   int nX = 1;
   if(Util::isValid(xPos))
      nX = count[xPos];
   std::vector<int> indices(dims.size(), 0);
   for(int e = 0; e < nEns; e++) {
      if(Util::isValid(ensPos))
         indices[ensPos] = e;
      for(int y = 0; y < nY; y++) {
         if(Util::isValid(yPos))
            indices[yPos] = y;
         for(int x = 0; x < nX; x++) {
            if(Util::isValid(xPos))
               indices[xPos] = x;
            int i = getIndex(countVector, indices);
            float value = values[i];
            if(Util::isValid(MV) && value == MV) {
               // Field has missing value indicator and the value is missing
               // Save values using our own internal missing value indicator
               value = Util::MV;
            }
            else {
               value = scale*value + offset;
            }
            (*field)(y,x,e) = value;
            // std::cout << e << " " << y << " " << x << " " << i << " " << value << std::endl;
         }
      }
   }
   // Alternative, but slower way. This loops over the vector and fills values into
   // the 3D array, instead of the otherway around.
   /*
   std::vector<int> indices(dims.size(), 0);
   for(int i = 0; i < size; i++) {
      getIndices(i, countVector, indices);
      int e = 0;
      int y = 0;
      int x = 0;
      if(Util::isValid(ensPos))
         e = indices[ensPos];
      if(Util::isValid(yPos))
         y = indices[yPos];
      if(Util::isValid(xPos))
         x = indices[xPos];
      float value = values[i];
      if(Util::isValid(MV) && value == MV) {
         // Field has missing value indicator and the value is missing
         // Save values using our own internal missing value indicator
         value = Util::MV;
      }
      else {
         value = scale*value + offset;
      }
      (*field)(y,x,e) = value;
   }
  */
   delete[] values;
   return field;
}

void FileNetcdf::writeCore(std::vector<Variable> iVariables, std::string iMessage) {
   startDefineMode();

   // check if altitudes are valid
   bool isAltitudeValid = false;
   vec2 elevs = getElevs();
   for(int i = 0; i < elevs.size(); i++) {
      for(int j = 0; j < elevs[i].size(); j++) {
         isAltitudeValid = isAltitudeValid || Util::isValid(elevs[i][j]);
      }
   }
   if(isAltitudeValid && !hasVar("altitude")) {
      defineAltitude();
   }

   defineTimes();
   defineEns();
   defineReferenceTime();
   defineGlobalAttributes(iMessage);
   // Define variables
   for(int v = 0; v < iVariables.size(); v++) {
      Variable variable = iVariables[v];
      std::string variableName = variable.name();

      if(!hasVariableCore(variable)) {
         // Create variable
         int numDims = Util::isValid(mYDim) + Util::isValid(mXDim) + Util::isValid(mEnsDim) + Util::isValid(mTimeDim);
         int dims[numDims];
         int counter = 0;
         if(Util::isValid(mTimeDim)) {
            dims[counter] = mTimeDim;
            counter++;
         }
         if(Util::isValid(mEnsDim)) {
            dims[counter] = mEnsDim;
            counter++;
         }
         if(Util::isValid(mYDim)) {
            dims[counter] = mYDim;
            counter++;
         }
         if(Util::isValid(mXDim)) {
            dims[counter] = mXDim;
            counter++;
         }

         int var = Util::MV;
         int status = nc_def_var(mFile, variableName.c_str(), NC_FLOAT, numDims, dims, &var);
         handleNetcdfError(status, "could not define variable '" + variableName + "'");
      }
      int var = getVar(variableName);
      float MV = getMissingValue(var); // The output file's missing value indicator

      // Set the coordinate attribute based on the names of the lat lon variables
      std::stringstream ss;
      ss << getVarName(mLonVar) << " " << getVarName(mLatVar);
      setAttribute(var, "coordinates", ss.str());

      if(variable.units() != "")
         setAttribute(var, "units", variable.units());
      if(variable.standardName() != "")
         setAttribute(var, "standard_name", variable.standardName());
   }
   startDataMode();

   writeTimes();
   writeReferenceTime();
   if(isAltitudeValid) {
      writeAltitude();
   }
   for(int v = 0; v < iVariables.size(); v++) {
      Variable variable = iVariables[v];
      std::string variableName = variable.name();
      assert(hasVariableCore(variable));
      int var = getVar(variableName);
      float MV = getMissingValue(var); // The output file's missing value indicator
      size_t size = 1*1*mNEns*getNumY()*getNumX();
      float* values = new float[size];

      std::vector<int> dims = getDims(var);
      size_t count[dims.size()];
      int ensPos = Util::MV;
      int yPos = Util::MV;
      int xPos = Util::MV;
      int timePos = Util::MV;
      for(int d = 0; d < dims.size(); d++) {
         int dim = dims[d];
         count[d] = 1;
         if(dim == mTimeDim) {
            timePos = d;
         }
         else if(dim == mEnsDim) {
            count[d] = getDimSize(dim);
            ensPos = d;
         }
         else if(dim == mYDim) {
            count[d] = getDimSize(dim);
            yPos = d;
         }
         else if(dim == mXDim) {
            count[d] = getDimSize(dim);
            xPos = d;
         }
      }

      for(int t = 0; t < getNumTime(); t++) {
         size_t start[dims.size()];
         for(int d = 0; d < dims.size(); d++)
            start[d] = 0;
         if(Util::isValid(timePos))
            start[timePos] = t;

         float offset = getOffset(var);
         float scale = getScale(var);
         FieldPtr field = getField(variable, t);
         if(field != NULL) { // TODO: Can't be null if coming from reference
            std::vector<int> countVector(count, count + dims.size());

            int nEns = 1;
            if(Util::isValid(ensPos))
               nEns = count[ensPos];
            int nY = 1;
            if(Util::isValid(yPos))
               nY = count[yPos];
            int nX = 1;
            if(Util::isValid(xPos))
               nX = count[xPos];
            std::vector<int> indices(dims.size(), 0);
            for(int e = 0; e < nEns; e++) {
               if(Util::isValid(ensPos))
                  indices[ensPos] = e;
               for(int y = 0; y < nY; y++) {
                  if(Util::isValid(yPos))
                     indices[yPos] = y;
                  for(int x = 0; x < nX; x++) {
                     if(Util::isValid(xPos))
                        indices[xPos] = x;
                     int i = getIndex(countVector, indices);
                     float value = (*field)(y, x, e);
                     if(Util::isValid(MV) && !Util::isValid(value)) {
                        // Field has missing value indicator and the value is missing
                        // Save values using the file's missing indicator value
                        value = MV;
                     }
                     else {
                        value = (value - offset)/scale;
                     }
                     values[i] = value;
                  }
               }
            }
            int status = nc_put_vara_float(mFile, var, start, count, values);
            handleNetcdfError(status, "could not write variable " + variableName);
         }
      }
      delete[] values;
   }
}


bool FileNetcdf::isValid(std::string iFilename, const Options& iOptions) {
   bool isValid = true;
   int file;
   int status = nc_open(iFilename.c_str(), NC_NOWRITE, &file);
   if(status == NC_NOERR) {
      // Check dimensions
      // isValid = isValid && (hasDim(file, "time") || iOptions.hasValue("timeDim"));
      // isValid = isValid && (hasDim(file, "lat") || hasDim(file, "latitude") || hasDim(file, "y") || iOptions.hasValue("yDim"));
      // isValid = isValid && (hasDim(file, "lon") || hasDim(file, "longitude") || hasDim(file, "x") || iOptions.hasValue("xDim"));
      // isValid = isValid && (hasDim(file, "ensemble_member") || iOptions.hasValue("ensDim"));

      // Check variables
      // isValid = isValid && (hasVar(file, "time") || iOptions.hasValue("timeVar"));
      isValid = isValid && (hasVar(file, "lat") || hasVar(file, "latitude") || iOptions.hasValue("latVar"));
      isValid = isValid && (hasVar(file, "lon") || hasVar(file, "longitude") || iOptions.hasValue("lonVar"));

      nc_close(file);
   }
   else
      isValid = false;
   return isValid;
}
vec2 FileNetcdf::getGridValues(int iVar) const {
   // Initialize values
   vec2 grid;
   grid.resize(getNumY());
   for(int i = 0; i < getNumY(); i++) {
      grid[i].resize(getNumX(), Util::MV);
   }

   // We have a lat/lon grid, where lat/lons are only provided along the pertinent dimension
   // Values are assumed to be constant across the other dimension.
   int numDims = getNumDims(iVar);
   if(numDims == 1) {
      int dim;
      nc_inq_vardimid(mFile, iVar, &dim);
      long size = getDimSize(dim);
      float* values = new float[size];
      nc_get_var_float(mFile, iVar, values);
      // Latitude variable
      if(dim == mYDim) {
         for(int i = 0; i < getNumY(); i++) {
            for(int j = 0; j < getNumX(); j++) {
               grid[i][j] = values[i];
            }
         }
      }
      // Longitude variable
      else if(dim == mXDim) {
         for(int i = 0; i < getNumY(); i++) {
            for(int j = 0; j < getNumX(); j++) {
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
      int N = getNumDims(iVar);
      size_t count[N];
      size_t start[N];
      int size = 1;
      int indexY = Util::MV;
      int indexX = Util::MV;
      int dims[N];
      nc_inq_vardimid(mFile, iVar, dims);
      for(int i = 0; i < N; i++) {
         if(dims[i] == mYDim) {
            count[i] = getNumY();
            size *= count[i];
            indexY = i;
         }
         else if(dims[i] == mXDim) {
            count[i] = getNumX();
            size *= count[i];
            indexX = i;
         }
         else {
            size_t dimsize = 1;
            nc_inq_dimlen(mFile, dims[i], &dimsize);
            if(dimsize > 1){
               std::stringstream ss;
               ss << "Lat/lon/elev has an extra non-singleton dimension (dim " << i << ") of length " << dimsize << ". Using index 0 to extract lat/lon/elev.";
               Util::warning(ss.str());
            }
            count[i] = 1;
         }
         start[i] = 0;
      }
      if(!Util::isValid(indexY) || !Util::isValid(indexX)) {
         std::stringstream ss;
         ss << "Missing lat and/or lon dimensions";
         Util::error(ss.str());
      }
      float* values = new float[size];
      nc_get_vara_float(mFile, iVar, start, count, values);
      for(int i = 0; i < getNumY(); i++) {
         for(int j = 0; j < getNumX(); j++) {
            // Latitude dimension is ordered first
            if(indexY < indexX) {
               grid[i][j] = values[i*getNumX() + j];
            }
            // Longitude dimension is ordered first
            else {
               grid[i][j] = values[j*getNumY() + i];
            }
         }
      }
      delete[] values;
   }
   return grid;
}

int FileNetcdf::detectYDim() const {
   int dim = Util::MV;
   if(hasDim("y"))
      dim = getDim("y");
   else if(hasDim("latitude"))
      dim = getDim("latitude");
   else if(hasDim("lat"))
      dim = getDim("lat");
   return dim;
}
int FileNetcdf::detectXDim() const {
   int dim = Util::MV;
   if(hasDim("x"))
      dim = getDim("x");
   else if(hasDim("longitude"))
      dim = getDim("longitude");
   else if(hasDim("lon"))
      dim = getDim("lon");
   return dim;
}
int FileNetcdf::detectTimeDim() const {
   int dim = Util::MV;
   if(hasDim("time"))
      dim = getDim("time");
   return dim;
}
int FileNetcdf::detectLatVar() const {
   int var = Util::MV;
   if(hasVar("latitude"))
      var = getVar("latitude");
   else if(hasVar("lat"))
      var = getVar("lat");
   else {
      std::stringstream ss;
      ss << "File '" << getFilename() << " does not have a latitude varable";
      Util::error(ss.str());
   }
   return var;
}
int FileNetcdf::detectTimeVar() const {
   int var = Util::MV;
   if(hasVar("time"))
      var = getVar("time");
   return var;
}
int FileNetcdf::detectLonVar() const {
   int var = Util::MV;
   if(hasVar("longitude"))
      var = getVar("longitude");
   else if(hasVar("lon"))
      var = getVar("lon");
   else {
      std::stringstream ss;
      ss << "File '" << getFilename() << " does not have a longitude varable";
      Util::error(ss.str());
   }
   return var;
}

void FileNetcdf::defineAltitude() {
   // Determine the order of lat and lon dimensions in altitude field
   int indexY = Util::MV;
   int indexX = Util::MV;

   // First check if latitude variable has two dimensions
   int vUseDims = Util::MV;
   if(getNumDims(mLatVar) >= 2) {
      vUseDims = mLatVar;
   }
   else if(hasVar("surface_geopotential")) {
      vUseDims = getVar("surface_geopotential");
   }
   if(Util::isValid(vUseDims)) {
      int N = getNumDims(vUseDims);
      int dims[N];
      nc_inq_vardimid(mFile, vUseDims, dims);
      for(int i = 0; i < N; i++) {
         if(dims[i] == mYDim)
            indexY = i;
         if(dims[i] == mXDim)
            indexX = i;
      }
   }
   else {
      Util::warning("Could not determine lat/lon ordering when creating altitude variable. Using [lat, lon]");
      indexY = 0;
      indexX = 1;
   }

   int dims[2];
   int dLat = mYDim;
   int dLon = mXDim;
   if(indexY < indexX) {
      dims[0]  = dLat;
      dims[1]  = dLon;
   }
   else {
      dims[0]  = dLon;
      dims[1]  = dLat;
   }
   int var = Util::MV;
   int status = nc_def_var(mFile, "altitude", NC_FLOAT, 2, dims, &var);
   handleNetcdfError(status, "could not define altitude");
}

void FileNetcdf::writeAltitude() const {
   int vElev = getVar("altitude");
   int numDims = getNumDims(vElev);
   int numHorizontalDims = Util::isValid(mYDim) + Util::isValid(mXDim);

   if(numDims != numHorizontalDims) {
      std::stringstream ss;
      ss << "Altitude variable has " << numDims << " dimension. However, the file has "
         << numHorizontalDims << " horizontal dimensions. Cannot write altitude.";
      Util::error(ss.str());
   }
   vec2 elevs = getElevs();

   int N = getNumDims(vElev);
   size_t count[N];
   size_t start[N];
   int size = 1;
   int indexY = Util::MV;
   int indexX = Util::MV;
   int dims[N];
   nc_inq_vardimid(mFile, vElev, dims);
   for(int i = 0; i < N; i++) {
      if(dims[i] == mYDim) {
         count[i] = getNumY();
         size *= count[i];
         indexY = i;
      }
      else if(dims[i] == mXDim) {
         count[i] = getNumX();
         size *= count[i];
         indexX = i;
      }
      else {
         size_t dimsize = 1;
         nc_inq_dimlen(mFile, dims[i], &dimsize);
         if(dimsize > 1){
            std::stringstream ss;
            ss << "Lat/lon/elev has an extra non-singleton dimension (dim " << i << ") of length " << dimsize << ". Using index 0 to write lat/lon/elev.";
            Util::warning(ss.str());
         }
         count[i] = 1;
      }
      start[i] = 0;
   }
   float MV = getMissingValue(vElev);
   float* values = new float[size];
   for(int i = 0; i < getNumY(); i++) {
      for(int j = 0; j < getNumX(); j++) {
         float elev = elevs[i][j];
         if(Util::isValid(MV) && !Util::isValid(elev))
            elev = MV;
         // Latitude dimension is ordered first
         if(indexY < indexX) {
            values[i*getNumX() + j] = elev;
         }
         // Longitude dimension is ordered first
         else {
            values[j*getNumY() + i] = elev;
         }
      }
   }
   int status = nc_put_vara_float(mFile, vElev, start, count, values);
   handleNetcdfError(status, "could not write altitude");
   delete[] values;
}
int FileNetcdf::getVar(std::string iVar) const {
   int var;
   int status = nc_inq_varid(mFile, iVar.c_str(), &var);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return var;
}

std::string FileNetcdf::getVarName(int iVar) const {
   char var[1000];
   int status = nc_inq_varname(mFile, iVar, var);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have variable '" << iVar << "'";
      Util::error(ss.str());
   }
   return std::string(var);
}

std::string FileNetcdf::getDimName(int iDim) const {
   char name[1000];
   int status = nc_inq_dimname(mFile, iDim, name);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return std::string(name);
}

int FileNetcdf::getDimSize(std::string iDim) const {
   if(!hasDim(iDim))
      return 0;
   int dim = getDim(iDim);
   size_t len = Util::MV;
   int status = nc_inq_dimlen(mFile, dim, &len);
   return len;
}

int FileNetcdf::getDimSize(int iDim) const {
   size_t len = Util::MV;
   int status = nc_inq_dimlen(mFile, iDim, &len);
   return len;
}

int FileNetcdf::getNumDims(int iVar) const {
   int len = Util::MV;
   int status = nc_inq_varndims(mFile, iVar, &len);
   return len;
}
std::string FileNetcdf::detectEnsDim() const {
   std::string dim = "ensemble_member";  // default
   if(hasDim("ensemble_member"))
      dim = "ensemble_member";
   return dim;
}
bool FileNetcdf::hasVar(int iFile, std::string iVar) {
   int var;
   int status = nc_inq_varid(iFile, iVar.c_str(), &var);
   return status == NC_NOERR;
}
bool FileNetcdf::hasVar(std::string iVar) const {
   return hasVar(mFile, iVar);
}
void FileNetcdf::handleNetcdfError(int status, std::string message) const {
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
void FileNetcdf::startDefineMode() const {
   if(mInDataMode) {
      int status = ncredef(mFile);
      handleNetcdfError(status, "could not put into define mode");
   }
   mInDataMode = false;
}
void FileNetcdf::startDataMode() const {
   if(!mInDataMode) {
      int status = ncendef(mFile);
      handleNetcdfError(status, "could not put into data mode");
   }
   mInDataMode = true;
}
std::vector<int> FileNetcdf::getDims(int iVar) const {
   int ndims = getNumDims(iVar);
   assert(ndims >= 0);
   std::vector<int> dims;
   int ids[ndims];
   int status = nc_inq_vardimid(mFile, iVar, ids);
   handleNetcdfError(status, "Could not get dimension");
   for(int i = 0; i < ndims; i++) {
      dims.push_back(ids[i]);
   }
   return dims;
}
float FileNetcdf::getMissingValue(int iVar) const {
   float fillValue;
   int status = nc_get_att_float(mFile, iVar, "_FillValue", &fillValue);
   if(status != NC_NOERR)
      fillValue  = NC_FILL_FLOAT;
   return fillValue;
}
void FileNetcdf::setMissingValue(int iVar, float iValue) const {
   // TODO: Mysterious errors can occur if the existing file was written
   // using an older HDF5 implementation, and the variable has 8 or more
   // attributes already. For more information, see "Corruption Problem
   // In HDF5 1.8.0 through HDF5 1.8.4" on
   // http://www.hdfgroup.org/HDF5/release/known_problems/index.html
   // Does this still apply when using the C-interface to NetCDF?
   if(iValue != NC_FILL_FLOAT) {
      int status = nc_put_att(mFile, iVar, "_FillValue", NC_FLOAT, 1, &iValue);
      handleNetcdfError(status, "could not set missing value flag");
   }
}

float FileNetcdf::getScale(int iVar) const {
   float scale;
   int status = nc_get_att_float(mFile, iVar, "scale_factor", &scale);
   if(status != NC_NOERR)
      scale  = 1;
   return scale;
}
float FileNetcdf::getOffset(int iVar) const {
   float offset;
   int status = nc_get_att_float(mFile, iVar, "add_offset", &offset);
   if(status != NC_NOERR)
      offset  = 0;
   return offset;
}

int FileNetcdf::getIndex(const std::vector<int>& iCount, const std::vector<int>& iIndices) const {
   int numDims = iCount.size();
   int index = 0;
   int sizeSoFar = 1;
   for(int k = numDims-1; k >= 0; k--) {
      index += iIndices[k] * sizeSoFar;
      sizeSoFar *= iCount[k];
   }
   return index ;
}
void FileNetcdf::getIndices(int i, const std::vector<int>& iCount, std::vector<int>& iIndices) const {
   // This is a bottleneck
   // The last index changes fastest, the first slowest
   int numDims = iCount.size();
   // iIndices.resize(numDims);
   int sizeSoFar = 1;
   for(int k = numDims-1; k >= 0; k--) {
      iIndices[k] = int(i / sizeSoFar) % iCount[k];
      sizeSoFar *= iCount[k];
   }
}

void FileNetcdf::setAttribute(std::string iVariable, std::string iName, std::string iValue) {
   int var = getVar(iVariable);
   setAttribute(var, iName, iValue);
}

void FileNetcdf::setAttribute(int iVar, std::string iName, std::string iValue) {
   startDefineMode();
   int status = nc_put_att_text(mFile, iVar,iName.c_str(), iValue.size(), iValue.c_str());
   handleNetcdfError(status, "could not set attribute");
}

void FileNetcdf::setGlobalAttribute(std::string iName, std::string iValue) {
   startDefineMode();
   int status = nc_put_att_text(mFile, NC_GLOBAL, iName.c_str(), iValue.size(), iValue.c_str());
   handleNetcdfError(status, "could not set global attribute");
}

void FileNetcdf::appendGlobalAttribute(std::string iName, std::string iValue) {
   int id;
   int status = nc_inq_attid(mFile, NC_GLOBAL, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, NC_GLOBAL, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine global attribute length");
      std::stringstream ss;
      if(len > mMaxAttributeLength) {
         std::stringstream ss0;
         ss0 << "Cannot append to attribute with length greater than " << mMaxAttributeLength << ". Resetting to new value";
         Util::warning(ss0.str());
      }
      else {
         char value[len+1];
         status = nc_get_att_text(mFile, NC_GLOBAL, iName.c_str(), value);
         handleNetcdfError(status, "could not append global attribute");
         value[len] = '\0';
         ss << value << "\n";
      }
      ss << iValue;
      setGlobalAttribute(iName, ss.str());
   }
   else {
      setGlobalAttribute(iName, iValue);
   }
}

void FileNetcdf::prependGlobalAttribute(std::string iName, std::string iValue) {
   int id;
   int status = nc_inq_attid(mFile, NC_GLOBAL, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, NC_GLOBAL, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine global attribute length");

      std::stringstream ss;
      ss << iValue;
      if(len > mMaxAttributeLength) {
         std::stringstream ss0;
         ss0 << "Cannot prepend to attribute with length greater than " << mMaxAttributeLength << ". Resetting to new value";
         Util::warning(ss0.str());
      }
      else {
         char value[len+1];
         status = nc_get_att_text(mFile, NC_GLOBAL, iName.c_str(), value);
         handleNetcdfError(status, "could not get attribute when prepending a new value");
         value[len] = '\0';
         ss << "\n" << value;
      }
      setGlobalAttribute(iName, ss.str());
   }
   else {
      setGlobalAttribute(iName, iValue);
   }
}

std::string FileNetcdf::getAttribute(std::string iVariable, std::string iName) {
   int var = getVar(iVariable);
   return getAttribute(var, iName);
}

std::string FileNetcdf::getAttribute(int iVar, std::string iName) {
   std::string ret = "";
   int id;
   int status = nc_inq_attid(mFile, iVar, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, iVar, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine attribute length");

      if(len > mMaxAttributeLength) {
         std::stringstream ss;
         ss << "Cannot get attribute with length greater than " << mMaxAttributeLength;
         Util::warning(ss.str());
      }
      else {
         char* value = new char[len+1];
         status = nc_get_att_text(mFile, iVar, iName.c_str(), value);
         handleNetcdfError(status, "could not get attribute");
         value[len] = '\0';
         if(status == NC_NOERR)
            ret = std::string(value);
         delete[] value;
      }
   }
   return ret;
}

std::string FileNetcdf::getGlobalAttribute(std::string iName) {
   std::string ret = "";
   int id;
   int status = nc_inq_attid(mFile, NC_GLOBAL, iName.c_str(), &id);
   if(status == NC_NOERR) {
      size_t len;
      int status = nc_inq_attlen(mFile, NC_GLOBAL, iName.c_str(), &len);
      handleNetcdfError(status, "could not determine global attribute length");

      if(len > mMaxAttributeLength) {
         std::stringstream ss;
         ss << "Cannot get attribute with length greater than " << mMaxAttributeLength;
         Util::warning(ss.str());
      }
      else {
         char* value = new char[len+1];
         status = nc_get_att_text(mFile, NC_GLOBAL, iName.c_str(), value);
         handleNetcdfError(status, "could not get global attribute");
         value[len] = '\0';
         if(status == NC_NOERR)
            ret = std::string(value);
         delete[] value;
      }
   }
   return ret;
}

void FileNetcdf::defineTimes() {
   if(!hasVar("time")) {
      if(Util::isValid(mTimeDim)) {
         int id;
         int status = nc_def_var(mFile, "time", NC_DOUBLE, 1, &mTimeDim, &id);
         handleNetcdfError(status, "creating time variable");
      }
      else {
         int status = nc_def_dim(mFile, "time", NULL, &mTimeDim);
         handleNetcdfError(status, "creating time dimension");
         status = nc_def_var(mFile, "time", NC_DOUBLE, 1, &mTimeDim, &mTimeVar);
         handleNetcdfError(status, "creating time variable");
      }
   }
   int vTime = getVar("time");
   setAttribute(vTime, "long_name", "time");
   setAttribute(vTime, "standard_name", "time");
   setAttribute(vTime, "units", "seconds since 1970-01-01 00:00:00 +00:00");
}
void FileNetcdf::defineEns() {
   if(!Util::isValid(mEnsDim) && mNEns > 1) {
      int status = nc_def_dim(mFile, mEnsDimName.c_str(), mNEns, &mEnsDim);
      handleNetcdfError(status, "creating '" + mEnsDimName + "' dimension");
   }
}
void FileNetcdf::writeTimes() {
   std::vector<double> times = getTimes();
   if(times.size() != getNumTime()) {
      std::stringstream ss;
      ss << "The times specified for NetCDF file '" << getFilename() << "' has " << times.size()
         << " elements, but the time dimension is " << getNumTime() << ". Putting missing values.";
      Util::warning(ss.str());
      times = std::vector<double>(getNumTime(), Util::MV);
   }

   // Convert missing
   for(int i = 0; i < times.size(); i++) {
      if(!Util::isValid(times[i]))
         times[i] = NC_FILL_FLOAT;
   }
   int vTime = getVar("time");
   double timesArr[getNumTime()];
   for(int t = 0; t < getNumTime(); t++) {
      timesArr[t] = times[t];
   }

   size_t start = 0;
   size_t count = getNumTime();
   int status = nc_put_vara_double(mFile, vTime, &start, &count, timesArr);
   handleNetcdfError(status, "could not write times");
}
void FileNetcdf::defineReferenceTime() {
   if(!hasVar("forecast_reference_time")) {
      int id;
      int status = nc_def_var(mFile, "forecast_reference_time", NC_DOUBLE, 0, NULL, &id);
      handleNetcdfError(status, "writing reference time");
   }
   int vTime = getVar("forecast_reference_time");
   setAttribute(vTime, "standard_name", "forecast_reference_time");
   setAttribute(vTime, "units", "seconds since 1970-01-01 00:00:00 +00:00");
}
void FileNetcdf::writeReferenceTime() {
   int vTime = getVar("forecast_reference_time");
   double referenceTime = getReferenceTime();
   if(!Util::isValid(referenceTime))
      referenceTime = NC_FILL_DOUBLE;
   int status = nc_put_var_double(mFile, vTime, &referenceTime);
   handleNetcdfError(status, "could not write reference time");
}
void FileNetcdf::defineGlobalAttributes(std::string iMessage) {
   if(getGlobalAttribute("Conventions") != "")
      setGlobalAttribute("Conventions", "CF-1.0");
   std::stringstream ss;
   ss << Util::getCurrentTimeStamp() << ": gridpp " << iMessage;
   prependGlobalAttribute("history", ss.str());
}
bool FileNetcdf::hasDim(int iFile, std::string iDim) {
   int dim;
   int status = nc_inq_dimid(iFile, iDim.c_str(), &dim);
   return status == NC_NOERR;
}
bool FileNetcdf::hasDim(std::string iDim) const {
   return hasDim(mFile, iDim);
}
bool FileNetcdf::hasVariableCore(const Variable& iVariable) const {
   std::string name = iVariable.name();
   int var=Util::MV;
   int status = nc_inq_varid(mFile, name.c_str(), &var);
   return status == NC_NOERR;
}
vec2 FileNetcdf::getLatLonVariable(int iVar) const {
   float MV = getMissingValue(iVar);
   std::vector<int> dims = getDims(iVar);
   // Initialize grid
   vec2 grid;
   int numY = 1;
   if(Util::isValid(mYDim))
      numY = getDimSize(mYDim);
   int numX = 1;
   if(Util::isValid(mXDim))
      numX = getDimSize(mXDim);

   grid.resize(numY);
   for(int i = 0; i < numY; i++) {
      grid[i].resize(numX);
   }

   if(dims.size() == 1) {
      // 1D, try to expand vector onto a grid
      float* values;
      if(dims[0] == mXDim)
         values = new float[numX];
      else if(dims[0] == mYDim)
         values = new float[numY];
      else {
         std::stringstream ss;
         ss << "Could not load latitude/longitude/altitude in " << getFilename() << " . Variable does not have x or y dimensions.";
         Util::error(ss.str());
      }

      bool status = nc_get_var_float(mFile, iVar, values);
      handleNetcdfError(status, "could not get data from latitude/longitude/altitude variable");

      for(int i = 0; i < numY; i++) {
         for(int j = 0; j < numX; j++) {
            float value = Util::MV;
            if(dims[0] == mXDim)
               value = values[j];
            else
               value = values[i];
            grid[i][j] = value;
         }
      }
      delete[] values;
   }
   else {
      int yPos = Util::MV;
      int xPos = Util::MV;
      size_t count[dims.size()];
      size_t start[dims.size()];
      for(int d = 0; d < dims.size(); d++) {
         count[d] = 1;
         start[d] = 0;
         if(dims[d] == mYDim) {
            yPos = d;
            count[d] = numY;
         }
         else if(dims[d] == mXDim) {
            xPos = d;
            count[d] = numX;
         }
      }
      if(!Util::isValid(xPos)) {
         std::stringstream ss;
         ss << "Latitude/longitude/altitude field in " << getFilename() << " has multiple dimensions, but it does not have the x dimension";
         Util::error(ss.str());
      }
      if(!Util::isValid(yPos)) {
         std::stringstream ss;
         ss << "Latitude/longitude/altitude field in " << getFilename() << " has multiple dimensions, but it does not have the y dimension";
         Util::error(ss.str());
      }

      float* values = new float[numX*numY];
      bool status = nc_get_vara_float(mFile, iVar, start, count, values);
      handleNetcdfError(status, "could not get data from latitude/longitude/altitude variable");
      for(int i = 0; i < numY; i++) {
         for(int j = 0; j < numX; j++) {
            int index = 0;
            if(yPos < xPos) {
               index = i*numX + j;
            }
            else {
               index = j*numY + i;
            }
            float value = values[index];
            if(values[index] == MV)
               value = Util::MV;
            grid[i][j] = value;
            assert(index < numX*numY);
         }
      }
      delete[] values;
   }
   return grid;
}

int FileNetcdf::getDim(std::string iDim) const {
   int dim;
   int status = nc_inq_dimid(mFile, iDim.c_str(), &dim);
   if(status != NC_NOERR) {
      std::stringstream ss;
      ss << "File '" << getFilename() << "' does not have dimension '" << iDim << "'";
      Util::error(ss.str());
   }
   return dim;
}

std::string FileNetcdf::description() {
   std::stringstream ss;
   ss << Util::formatDescription("type=netcdf", "Netcdf file") << std::endl;
   ss << Util::formatDescription("   yDim=undef", "Name of y-axis dimension. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   xDim=undef", "Name of x-axis dimension. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   timeDim=undef", "Name of time dimension. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   ensDim=undef", "Name of ensemble dimension. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   latVar=undef", "Name of variable with latitudes. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   lonVar=undef", "Name of variable with longitudes. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   timeVar=undef", "Name of variable with time. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   elevVar=undef", "Name of altitude variable. If unspecified, 'altitude' or 'surface_geopotential' is used.") << std::endl;
   ss << Util::formatDescription("   lafVar=undef", "Name of land-area-fraction variable. Auto-detected if unspecified.") << std::endl;
   ss << Util::formatDescription("   variables=undef", "Variable definition file.") << std::endl;
   return ss.str();
}
