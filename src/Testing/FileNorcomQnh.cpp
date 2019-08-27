#include "../File/NorcomQnh.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileNorcomQnhTest : public ::testing::Test {
      protected:
         virtual void SetUp() {
             mVariable = Variable("surface_air_pressure");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(FileNorcomQnhTest, options) {
      FileNorcomQnh file("testing/files/test.txt", Options("lats=1,2 lons=2,3 elevs=100,120 names=point1,point2 numTimes=2 startTime=0 endTime=2"));
      vec2 lats = file.getLats();
      vec2 lons = file.getLons();
      vec2 elevs = file.getElevs();
      ASSERT_EQ(1, lats.size());
      ASSERT_EQ(1, lons.size());
      ASSERT_EQ(1, elevs.size());
      ASSERT_EQ(2, lats[0].size());
      ASSERT_EQ(2, lons[0].size());
      ASSERT_EQ(2, elevs[0].size());
      EXPECT_FLOAT_EQ(1, lats[0][0]);
      EXPECT_FLOAT_EQ(2, lats[0][1]);
      EXPECT_FLOAT_EQ(2, lons[0][0]);
      EXPECT_FLOAT_EQ(3, lons[0][1]);
      EXPECT_FLOAT_EQ(100, elevs[0][0]);
      EXPECT_FLOAT_EQ(120, elevs[0][1]);
   }
   TEST_F(FileNorcomQnhTest, asOutput) {
      FileNetcdf from("testing/files/10x10.nc");
      FileNorcomQnh to("testing/files/test.txt", Options("lats=1,2 lons=2,3 elevs=100,120 names=point1,point2 numTimes=2 startTime=0 endTime=1"));
      DownscalerNearestNeighbour d = DownscalerNearestNeighbour(mVariable, mVariable, Options());
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      std::vector<Variable> variables(1, mVariable);
      to.write(variables);

      // FieldPtr field = to.getField(mVariable, 0);
      // EXPECT_FLOAT_EQ(1, (*field)(0,0,0));
   }
   TEST_F(FileNorcomQnhTest, description) {
      FileNorcomQnh::description();
   }
   TEST_F(FileNorcomQnhTest, valid) {
      FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=300 elevs=3 numTimes=2 startTime=0 endTime=1 names=test"));
   }
   TEST_F(FileNorcomQnhTest, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // lats, lons, elevs not the same
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1,2 lons=2 elevs=3 names=test numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=2 lons=2,3,2 elevs=3 names=test numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=2 lons=2 elevs=3,2 names=test numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 names=q,w numTimes=2 startTime=0 endTime=1")), ".*");

      // Missing attributes
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lons=2 elevs=3 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 elevs=3 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 numTimes=2 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 numTimes=2 startTime=0")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("")), ".*");

      // Invalid latitude
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=101 lons=2 elevs=3 numTimes=2 startTime=0 endTime=1")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=-91 lons=180 elevs=3 numTimes=2 startTime=0 endTime=1")), ".*");

      // Longitudes outside 360, -360 should not be allowed
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=370 elevs=3 numTimes=2 startTime=0 endTime=1 names=test")), ".*");
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=-370 elevs=3 numTimes=2 startTime=0 endTime=1 names=test")), ".*");

      // Start/end time not in ascending order
      EXPECT_DEATH(FileNorcomQnh("testing/files/test.txt", Options("lats=1 lons=2 elevs=3 names=q numTimes=2 startTime=1 endTime=0")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
