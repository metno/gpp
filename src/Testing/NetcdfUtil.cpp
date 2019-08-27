#include "../NetcdfUtil.h"
#include "../Variable.h"
#include "../Calibrator/Calibrator.h"
#include <algorithm>
#include <math.h>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <unistd.h>
#include <netcdf.h>

namespace {
   class NetcdfUtilTest : public ::testing::Test {
      protected:
   };

   TEST_F(NetcdfUtilTest, get) {
      int file;
      int status = nc_open("testing/files/10x10.nc", NC_NOWRITE, &file);
      EXPECT_EQ(status, NC_NOERR);
      int var;
      status = nc_inq_varid(file, "air_temperature_2m", &var);
      EXPECT_EQ(status, NC_NOERR);
      EXPECT_EQ(NC_FILL_FLOAT, NetcdfUtil::getMissingValue(file, var));
      EXPECT_EQ(200, NetcdfUtil::getTotalSize(file, var));
   }
   TEST_F(NetcdfUtilTest, error) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(NetcdfUtil::handleNetcdfError(NC_EBADID, "error"), ".*");
      EXPECT_DEATH(NetcdfUtil::handleNetcdfError(NC_EBADID, ""), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
