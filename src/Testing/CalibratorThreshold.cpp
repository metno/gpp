#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Threshold.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorThreshold : public ::testing::Test {
      protected:
         TestCalibratorThreshold() {
         }
         virtual ~TestCalibratorThreshold() {
         }
         virtual void SetUp() {
            mTemperature = Variable("air_temperature_2m");
            mPrecipitation = Variable("precipitation_amount_acc");
         }
         virtual void TearDown() {
         }
         Variable mTemperature;
         Variable mPrecipitation;
   };

   TEST_F(TestCalibratorThreshold, 1x1) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorThreshold cal = CalibratorThreshold(mTemperature, Options("thresholds=20 values=0,2"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(2, (*from.getField(mTemperature, 0))(0,0,0));
      EXPECT_FLOAT_EQ(2, (*from.getField(mTemperature, 1))(0,0,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(mTemperature, 2))(0,0,0));
      EXPECT_FLOAT_EQ(2, (*from.getField(mTemperature, 3))(0,0,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(mTemperature, 4))(0,0,0));
      EXPECT_FLOAT_EQ(2, (*from.getField(mTemperature, 5))(0,0,0));
      EXPECT_FLOAT_EQ(2, (*from.getField(mTemperature, 6))(0,0,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(mTemperature, 7))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mTemperature, 8))(0,0,0));
      EXPECT_FLOAT_EQ(2, (*from.getField(mTemperature, 9))(0,0,0));
   }
   TEST_F(TestCalibratorThreshold, 1x1_equals) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorThreshold cal = CalibratorThreshold(mPrecipitation, Options("thresholds=3,3.5,4 values=-5,11,0,2 equals=0,1,0"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(-5, (*from.getField(mPrecipitation, 0))(0,0,0)); // 0
      EXPECT_FLOAT_EQ(11, (*from.getField(mPrecipitation, 1))(0,0,0)); // 3
      EXPECT_FLOAT_EQ(2, (*from.getField(mPrecipitation, 2))(0,0,0)); // 4
      EXPECT_FLOAT_EQ(2, (*from.getField(mPrecipitation, 4))(0,0,0)); // 5.5
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mPrecipitation, 6))(0,0,0)); // MV
      EXPECT_FLOAT_EQ(2, (*from.getField(mPrecipitation, 7))(0,0,0)); // 12
   }
   TEST_F(TestCalibratorThreshold, 1x1_equals2) {
      // Check that equals=1 works on the upper limit
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorThreshold cal = CalibratorThreshold(mPrecipitation, Options("thresholds=3,3.5,10 values=-5,11,0,2 equals=1,0,1"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(-5, (*from.getField(mPrecipitation, 1))(0,0,0)); // 3
      EXPECT_FLOAT_EQ(0, (*from.getField(mPrecipitation, 5))(0,0,0)); // 10
      EXPECT_FLOAT_EQ(2, (*from.getField(mPrecipitation, 7))(0,0,0)); // 12
   }
   TEST_F(TestCalibratorThreshold, description) {
      CalibratorThreshold::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
