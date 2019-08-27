#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Accumulate.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorAccumulate : public ::testing::Test {
      protected:
         TestCalibratorAccumulate() {
         }
         virtual ~TestCalibratorAccumulate() {
         }
         virtual void SetUp() {
            mTemperature = Variable("air_temperature_2m");
            mPrecipitation = Variable("precipitation_amount");
            mMissing = Variable("qnh");
         }
         virtual void TearDown() {
         }
         Variable mTemperature;
         Variable mPrecipitation;
         Variable mMissing;
   };

   TEST_F(TestCalibratorAccumulate, 1x1) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorAccumulate cal = CalibratorAccumulate(mTemperature, Options());
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(0, (*from.getField(mTemperature, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20, (*from.getField(mTemperature, 1))(0,0,0));
      EXPECT_FLOAT_EQ(35, (*from.getField(mTemperature, 2))(0,0,0));
      EXPECT_FLOAT_EQ(56, (*from.getField(mTemperature, 3))(0,0,0));
      EXPECT_FLOAT_EQ(70, (*from.getField(mTemperature, 4))(0,0,0));
      EXPECT_FLOAT_EQ(100, (*from.getField(mTemperature, 5))(0,0,0));
      EXPECT_FLOAT_EQ(121, (*from.getField(mTemperature, 6))(0,0,0));
      EXPECT_FLOAT_EQ(140, (*from.getField(mTemperature, 7))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mTemperature, 8))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mTemperature, 9))(0,0,0));
   }
   TEST_F(TestCalibratorAccumulate, 10x10) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorAccumulate cal = CalibratorAccumulate(mPrecipitation, Options());
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(0, (*from.getField(mPrecipitation, 0))(5,2,0));
      EXPECT_FLOAT_EQ(0.539526, (*from.getField(mPrecipitation, 1))(5,2,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(mPrecipitation, 0))(5,9,0));
      EXPECT_FLOAT_EQ(6.929162, (*from.getField(mPrecipitation, 1))(5,9,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(mPrecipitation, 0))(0,9,0));
      EXPECT_FLOAT_EQ(5.442121, (*from.getField(mPrecipitation, 1))(0,9,0));
   }
   TEST_F(TestCalibratorAccumulate, description) {
      CalibratorAccumulate::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
