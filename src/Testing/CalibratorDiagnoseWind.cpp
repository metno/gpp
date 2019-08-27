#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/DiagnoseWind.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorDiagnoseWind : public ::testing::Test {
      protected:
         virtual void SetUp() {
            mX = Variable("x_wind_10m");
            mY = Variable("y_wind_10m");
            mS = Variable("speed");
            mD = Variable("direction");
            calSpeed = new CalibratorDiagnoseWind(mS,Options("x=x_wind_10m y=y_wind_10m compute=speed"));
            calDir   = new CalibratorDiagnoseWind(mD,Options("x=x_wind_10m y=y_wind_10m compute=direction"));
            calX     = new CalibratorDiagnoseWind(mX,Options("speed=speed direction=direction compute=x"));
            calY     = new CalibratorDiagnoseWind(mY,Options("speed=speed direction=direction compute=y"));
         }
         virtual void TearDown() {
            delete calSpeed;
            delete calDir;
            delete calX;
            delete calY;
         }
         Variable mX;
         Variable mY;
         Variable mS;
         Variable mD;
         CalibratorDiagnoseWind* calSpeed;
         CalibratorDiagnoseWind* calDir;
         CalibratorDiagnoseWind* calX;
         CalibratorDiagnoseWind* calY;
   };
   TEST_F(TestCalibratorDiagnoseWind, speed) {
      FileNetcdf from("testing/files/10x10.nc");
      FieldPtr fieldX = from.getField(mX, 0);
      FieldPtr fieldY = from.getField(mY, 0);
      // Set some values to missing
      (*fieldX)(1, 2, 0) = Util::MV;
      (*fieldY)(2, 2, 0) = Util::MV;
      (*fieldY)(3, 2, 0) = Util::MV;

      calSpeed->calibrate(from, NULL);
      FieldPtr field = from.getField(mS, 0);
      EXPECT_FLOAT_EQ(1.8898070623831842, (*field)(5,2,0));  // x: -1.886 y: -0.108
      EXPECT_FLOAT_EQ(Util::MV, (*field)(1,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*field)(2,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*field)(3,2,0));
   }
   TEST_F(TestCalibratorDiagnoseWind, direction) {
      FileNetcdf from("testing/files/10x10.nc");
      FieldPtr fieldX = from.getField(mX, 0);
      FieldPtr fieldY = from.getField(mY, 0);
      // Set some values to missing
      (*fieldX)(1, 2, 0) = Util::MV;
      (*fieldY)(2, 2, 0) = Util::MV;
      (*fieldY)(3, 2, 0) = Util::MV;

      // Check that we don't get negative directions
      (*fieldX)(0, 0, 0) = 1;
      (*fieldY)(0, 0, 0) = 0;

      calDir->calibrate(from, NULL);
      FieldPtr field = from.getField(mD, 0);
      EXPECT_FLOAT_EQ(86.715614, (*field)(5,2,0));  // x: -1.886 y: -0.108
      EXPECT_FLOAT_EQ(Util::MV, (*field)(1,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*field)(2,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*field)(3,2,0));
      EXPECT_FLOAT_EQ(270, (*field)(0,0,0)); // x: 1 y: 0
   }
   TEST_F(TestCalibratorDiagnoseWind, xy) {
      // Diagnose windspeed and direction from u and v. Then change the windspeed to
      // 13 and check that rediagnosis of U/V yields correct winds
      FileNetcdf from("testing/files/10x10.nc");
      FieldPtr fieldX = from.getField(mX, 0);
      FieldPtr fieldY = from.getField(mY, 0);
      FieldPtr fieldS = from.getField(mS, 0);

      calSpeed->calibrate(from, NULL);
      calDir->calibrate(from, NULL);

      (*fieldS)(5,2,0) = 13;

      // Re-diagnose the new U and V
      calX->calibrate(from, NULL);
      calY->calibrate(from, NULL);

      FieldPtr field = from.getField(mD, 0);

      EXPECT_FLOAT_EQ(-12.978647, (*fieldX)(5,2,0));
      EXPECT_FLOAT_EQ(-0.74479485, (*fieldY)(5,2,0));
   }
   TEST_F(TestCalibratorDiagnoseWind, description) {
      CalibratorDiagnoseWind::description();
   }
   TEST_F(TestCalibratorDiagnoseWind, valid) {
      CalibratorDiagnoseWind(mX, Options("speed=s direction=d compute=x"));
      CalibratorDiagnoseWind(mY, Options("speed=s direction=d compute=y"));
      CalibratorDiagnoseWind(mS, Options("x=x y=y compute=speed"));
      CalibratorDiagnoseWind(mD, Options("x=x y=y compute=direction"));
   }
   TEST_F(TestCalibratorDiagnoseWind, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // Missing the correct fields to compute x
      EXPECT_DEATH(CalibratorDiagnoseWind(mX, Options("x=x y=y compute=x")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseWind(mX, Options("speed=speed compute=x")), ".*");
      // Missing the correct fields to compute y
      EXPECT_DEATH(CalibratorDiagnoseWind(mY, Options("x=x y=y compute=y")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseWind(mY, Options("speed=speed compute=y")), ".*");

      // Missing the correct fields to compute speed
      EXPECT_DEATH(CalibratorDiagnoseWind(mS, Options("x=x compute=speed")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseWind(mS, Options("x=x speed=speed compute=speed")), ".*");
      // Missing the correct fields to compute direction
      EXPECT_DEATH(CalibratorDiagnoseWind(mD, Options("x=x compute=direction")), ".*");
      EXPECT_DEATH(CalibratorDiagnoseWind(mD, Options("x=x speed=speed compute=direction")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
