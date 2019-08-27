#include "../File/Point.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FilePointTest : public ::testing::Test {
      protected:
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(FilePointTest, asInput) {
      FilePoint file("testing/files/validPoint1.txt", Options("lat=1 lon=2 elev=3"));
      FieldPtr field0 = file.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(290, (*field0)(0,0,0));
      FieldPtr field1 = file.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(288, (*field1)(0,0,0));
   }
   TEST_F(FilePointTest, asOutput) {
      {
         FileNetcdf from("testing/files/10x10.nc");
         FilePoint to("testing/files/filePoint.txt", Options("lat=1 lon=2 elev=3 time=2 ens=1"));
         DownscalerNearestNeighbour d = DownscalerNearestNeighbour(mVariable, mVariable, Options());
         bool status = d.downscale(from, to);
         EXPECT_TRUE(status);
         std::vector<Variable> variables(1, mVariable);
         to.write(variables);
      }
      // Check that the right values are written
      FilePoint from("testing/files/filePoint.txt", Options("lat=1 lon=2 elev=3 time=2"));
      FieldPtr field = from.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(303, (*field)(0,0,0));
   }
   TEST_F(FilePointTest, asEnsemble) {
      FilePoint file("testing/files/validPoint2.txt", Options("lat=1 lon=2 elev=3"));
      FieldPtr field0 = file.getField(mVariable, 0);
      ASSERT_EQ(2, field0->getNumEns());
      EXPECT_FLOAT_EQ(290, (*field0)(0,0,0));
      EXPECT_FLOAT_EQ(291, (*field0)(0,0,1));
      FieldPtr field1 = file.getField(mVariable, 1);
      ASSERT_EQ(2, field1->getNumEns());
      EXPECT_FLOAT_EQ(288, (*field1)(0,0,0));
      EXPECT_FLOAT_EQ(300, (*field1)(0,0,1));
   }
   TEST_F(FilePointTest, validFiles) {
      FilePoint file1("testing/files/validPoint1.txt", Options("lat=1 lon=2 elev=3 time=67"));
      FilePoint file2("testing/files/validPoint2.txt", Options("lat=1 lon=2 elev=3 time=67"));
      FilePoint file3("testing/files/validPoint2.txt", Options("lat=89 lon=2 elev=3 time=67"));
      FilePoint file4("testing/files/validPoint2.txt", Options("lat=-89 lon=-180 elev=3 time=67"));
      FilePoint file5("testing/files/validPoint2.txt", Options("lat=-89 lon=180 elev=3 time=67"));
      FilePoint file6("testing/files/validPoint2.txt", Options("lat=-89 lon=180 elev=-32 time=67"));
      FilePoint file7("testing/files/validPoint1.txt", Options("lat=89 lon=200 elev=3 time=67"));
      FilePoint file8("testing/files/validPoint1.txt", Options("lat=89 lon=-200 elev=3 time=67"));
   }
   TEST_F(FilePointTest, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(FilePoint("testing/files/validPoint1.txt", Options("lon=2 elev=3 time=67")), ".*");
      EXPECT_DEATH(FilePoint("testing/files/validPoint1.txt", Options("lat=1 elev=3 time=67")), ".*");
      EXPECT_DEATH(FilePoint("testing/files/validPoint1.txt", Options("lat=1 lon=2 time=67")), ".*");
      // Invalid lat
      EXPECT_DEATH(FilePoint("testing/files/validPoint1.txt", Options("lat=91 lon=2 elev=3 time=67")), ".*");
      EXPECT_DEATH(FilePoint("testing/files/validPoint1.txt", Options("lat=-91 lon=2 elev=3 time=67")), ".*");
      // Missing time for non-existant file
      EXPECT_DEATH(FilePoint("testing/files/hd92h3d98h38.txt", Options("lat=1 lon=2 elev=3")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
