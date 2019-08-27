#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Mask.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorMask : public ::testing::Test {
      protected:
         TestCalibratorMask() {
         }
         virtual ~TestCalibratorMask() {
         }
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
         std::vector<float> getVector(float iArray[]) {
            return std::vector<float>(iArray, iArray + sizeof(iArray)/sizeof(float));
         }
   };
   // Mask out values
   TEST_F(TestCalibratorMask, 10x10_mask_out) {
      FileNetcdf from("testing/files/10x10.nc");
      // One point at 3,5 with 223 km radius and one point at 4,6 with 336 km radius
      ParameterFileText par(Options("file=testing/files/mask0.txt"));
      CalibratorMask cal = CalibratorMask(mVariable, Options("keep=0"));

      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(3,5,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(3,3,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(2,5,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(4,9,0));
      EXPECT_FLOAT_EQ(302, (*after)(2,3,0));
      EXPECT_FLOAT_EQ(310, (*after)(6,9,0));
   }
   // Mask in values
   TEST_F(TestCalibratorMask, 10x10_mask_in) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/mask0.txt"));
      CalibratorMask cal = CalibratorMask(mVariable, Options("mask=1"));

      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(Util::MV, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(302, (*after)(3,5,0));
      EXPECT_FLOAT_EQ(302, (*after)(3,3,0));
      EXPECT_FLOAT_EQ(302, (*after)(2,5,0));
      EXPECT_FLOAT_EQ(302, (*after)(4,9,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(2,3,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(6,9,0));
   }
   // Missing parameter file
   TEST_F(TestCalibratorMask, invalid2) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
	  FileNetcdf from("testing/files/10x10.nc");
      Util::setShowError(false);
      CalibratorMask calibrator(mVariable, Options());
      EXPECT_DEATH(calibrator.calibrate(from, NULL), ".*");
   }
   TEST_F(TestCalibratorMask, description) {
      CalibratorMask::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
