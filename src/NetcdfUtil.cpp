#include "NetcdfUtil.h"
#include <sstream>

float NetcdfUtil::getMissingValue(int iFile, int iVar) {
   float fillValue;
   int status = nc_get_att_float(iFile, iVar, "_FillValue", &fillValue);
   if(status != NC_NOERR)
      fillValue  = NC_FILL_FLOAT;
   return fillValue;
}

long NetcdfUtil::getTotalSize(int iFile, int iVar) {
   int ndims;
   int status = nc_inq_varndims(iFile, iVar, &ndims);
   handleNetcdfError(status, "Could not get number of dimensions");

   int dims[ndims];
   status = nc_inq_vardimid(iFile, iVar, dims);
   handleNetcdfError(status, "Could not get dimension ids");

   long total = 1;
   for(int i = 0; i < ndims; i++) {
      size_t size;
      int status = nc_inq_dimlen(iFile, dims[i], &size);
      handleNetcdfError(status,"Could not get size of dimension");
      total = total * size;
   }
   return total;

}

void NetcdfUtil::handleNetcdfError(int status, std::string message) {
   if(status != NC_NOERR) {
      std::stringstream ss;
      if(message == "") {
         ss << "Netcdf error code: " << status << ".";
      }
      else {
         ss << "Netcdf error: " << message << ". "
            << "Netcdf error code: " << status << ".";
      }
      Util::error(ss.str());
   }
}
