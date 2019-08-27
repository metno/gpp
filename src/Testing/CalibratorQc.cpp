#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Qc.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorQc : public ::testing::Test {
      protected:
         TestCalibratorQc() {
         }
         virtual ~TestCalibratorQc() {
         }
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };
   TEST_F(TestCalibratorQc, 10x10) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorQc cal = CalibratorQc(mVariable, Options("min=304 max=305.8"));

      cal.calibrate(from);

      FieldPtr after = from.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(304, (*after)(5,2,0)); // 301 (row,col)
      EXPECT_FLOAT_EQ(304, (*after)(5,9,0)); // 304
      EXPECT_FLOAT_EQ(305.8, (*after)(0,9,0)); // 320
   }
   TEST_F(TestCalibratorQc, 10x10_nomax) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorQc cal = CalibratorQc(mVariable ,Options("max=307"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(304, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(307, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorQc, 10x10_nomin) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorQc cal = CalibratorQc(mVariable ,Options("min=303"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(303, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(304, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(320, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorQc, 10x10_no) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorQc cal = CalibratorQc(mVariable ,Options(""));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(304, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(320, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorQc, 10x10_missingValue) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorQc cal = CalibratorQc(mVariable ,Options("min=303 max=307"));

      FieldPtr after = from.getField(mVariable, 0);
      (*after)(5,2,0) = Util::MV;
      (*after)(5,9,0) = Util::MV;
      (*after)(0,9,0) = Util::MV;
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(Util::MV, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorQc, getRadius) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
   }
   TEST_F(TestCalibratorQc, description) {
      CalibratorQc::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
