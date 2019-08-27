#ifndef NETCDF_UTIL_H
#define NETCDF_UTIL_H
#include "Util.h"
#include <netcdf.h>
#include <string>

class NetcdfUtil {
   public:
      static float getMissingValue(int iFile, int iVar);
      static long  getTotalSize(int iFile, int iVar);
      static void handleNetcdfError(int status, std::string message="");
};
#endif
