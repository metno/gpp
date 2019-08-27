#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Deaccumulate.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDeaccumulate : public ::testing::Test {
      protected:
         TestCalibratorDeaccumulate() {
         }
         virtual ~TestCalibratorDeaccumulate() {
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

   TEST_F(TestCalibratorDeaccumulate, 1x1) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorDeaccumulate cal = CalibratorDeaccumulate(mTemperature, Options());
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mTemperature, 0))(0,0,0));
      EXPECT_FLOAT_EQ(-3, (*from.getField(mTemperature, 1))(0,0,0));
      EXPECT_FLOAT_EQ(-5, (*from.getField(mTemperature, 2))(0,0,0));
      EXPECT_FLOAT_EQ(6, (*from.getField(mTemperature, 3))(0,0,0));
      EXPECT_FLOAT_EQ(-7, (*from.getField(mTemperature, 4))(0,0,0));
      EXPECT_FLOAT_EQ(16, (*from.getField(mTemperature, 5))(0,0,0));
      EXPECT_FLOAT_EQ(-9, (*from.getField(mTemperature, 6))(0,0,0));
      EXPECT_FLOAT_EQ(-2, (*from.getField(mTemperature, 7))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mTemperature, 8))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mTemperature, 9))(0,0,0));
   }
   TEST_F(TestCalibratorDeaccumulate, description) {
      CalibratorDeaccumulate::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
