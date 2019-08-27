#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>


namespace {
   class ParameterFileNetcdfTest : public ::testing::Test {
      public:
         Parameters createParameters(float i1, float i2, float i3) {
            std::vector<float> values;
            values.push_back(i1);
            values.push_back(i2);
            values.push_back(i3);
            return Parameters(values);
         };
      protected:
   };
   TEST_F(ParameterFileNetcdfTest, singleTime) {
      ParameterFileNetcdf file(Options("file=testing/files/10x10_param.nc"));
      ASSERT_EQ(2, file.getTimes().size());
      ASSERT_EQ(100, file.getLocations().size());
      // Location loc(5,5,142.1456);
      Location loc(5,5,3);
      Parameters par = file.getParameters(0, loc);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(0.3, par[0]);
      EXPECT_FLOAT_EQ(2.3, par[1]);
      par = file.getParameters(1, loc);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(1, par[0]);
      EXPECT_FLOAT_EQ(2, par[1]);
   }
   TEST_F(ParameterFileNetcdfTest, singleTime_xy) {
      ParameterFileNetcdf file(Options("file=testing/files/10x10_param_xy.nc"));
      ASSERT_EQ(2, file.getTimes().size());
      ASSERT_EQ(100, file.getLocations().size());
      // Location loc(5,5,142.1456);
      Location loc(5,5,3);
      Parameters par = file.getParameters(0, loc);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(0.3, par[0]);
      EXPECT_FLOAT_EQ(2.3, par[1]);
   }

   // 10x10_param_xy.nc has longitude and coefficients with x,y ordering
   // Should still give you the same results
   TEST_F(ParameterFileNetcdfTest, xy_order) {
      ParameterFileNetcdf file_yx(Options("file=testing/files/10x10_param.nc"));
      ParameterFileNetcdf file_xy(Options("file=testing/files/10x10_param_xy.nc"));
      EXPECT_EQ(file_yx.getTimes(), file_xy.getTimes());
      std::vector<Location> loc_yx = file_yx.getLocations();
      std::vector<Location> loc_xy = file_yx.getLocations();
      for(int i = 0; i < loc_yx.size(); i++) {
         if(loc_yx[i].lat() == 5 && loc_yx[i].lon() == 5) {
            Parameters par_yx= file_yx.getParameters(0, loc_yx[i]);
            Parameters par_xy= file_xy.getParameters(0, loc_yx[i]);
            EXPECT_EQ(par_yx.getValues(), par_xy.getValues());
         }
      }
   }
   // Should be possible to write to file
   TEST_F(ParameterFileNetcdfTest, write) {
      // Write parameters (make sure it goes out of scope)
      {
         ParameterFileNetcdf file(Options("file=testing/files/test192837.nc"), true);
         file.setParameters(createParameters(1,2,3), 0, Location(3,2,5));
         file.setParameters(createParameters(4,5,6), 1, Location(3,2,5));
         file.setParameters(createParameters(7,8,9), 1, Location(1,9,2));
         file.recomputeTree();
         file.write();
      }

      // Read them again
      ParameterFileNetcdf file(Options("file=testing/files/test192837.nc"), false);
      std::vector<Location> locations = file.getLocations();
      ASSERT_EQ(2, locations.size());
      // Time 0
      Parameters par = file.getParameters(0, Location(3,2,5));
      ASSERT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(1, par[0]);
      EXPECT_FLOAT_EQ(2, par[1]);
      EXPECT_FLOAT_EQ(3, par[2]);
      /* Removed until netcdf parameter files can handle cases where parameters are available for
       * some times but not others
      par = file.getParameters(0, Location(1,9,2));
      ASSERT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(1, par[0]);
      EXPECT_FLOAT_EQ(2, par[1]);
      EXPECT_FLOAT_EQ(3, par[2]);
      */
      // Expect missing values for now, until above is fixed
      {
      Parameters par = file.getParameters(0, Location(1,9,2));
      ASSERT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(Util::MV, par[0]);
      EXPECT_FLOAT_EQ(Util::MV, par[1]);
      EXPECT_FLOAT_EQ(Util::MV, par[2]);
      }

      // Time 1
      {
      Parameters par = file.getParameters(1, Location(3,2,5));
      ASSERT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(4, par[0]);
      EXPECT_FLOAT_EQ(5, par[1]);
      EXPECT_FLOAT_EQ(6, par[2]);
      par = file.getParameters(1, Location(1,9,2));
      ASSERT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(7, par[0]);
      EXPECT_FLOAT_EQ(8, par[1]);
      EXPECT_FLOAT_EQ(9, par[2]);
      }

      // Util::remove("testing/files/test192837.nc");
   }
   TEST_F(ParameterFileNetcdfTest, invalidFiles) {
      EXPECT_FALSE(ParameterFileText::isValid("testing/files/parametersf98wey8y8y89rwe.nc"));
      ParameterFileText p(Options("testing/files/parametersf98wey8y8y89rwe.nc"));
      EXPECT_FALSE(p.isReadable());
   }
   TEST_F(ParameterFileNetcdfTest, description) {
      ParameterFileNetcdf::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
