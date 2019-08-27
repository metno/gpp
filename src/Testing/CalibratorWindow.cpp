#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Window.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorWindow : public ::testing::Test {
      protected:
         TestCalibratorWindow() {
         }
         virtual ~TestCalibratorWindow() {
         }
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };
   TEST_F(TestCalibratorWindow, length1) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=1 stat=mean"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20, (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(15, (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(21, (*from.getField(mVariable, 3))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(mVariable, 4))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(mVariable, 5))(0,0,0));
      EXPECT_FLOAT_EQ(21, (*from.getField(mVariable, 6))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(mVariable, 8))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, length5) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=5 stat=mean"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(19.333333, (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(19.75,     (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(18.6,      (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(20,        (*from.getField(mVariable, 3))(0,0,0));
      EXPECT_FLOAT_EQ(20.2,      (*from.getField(mVariable, 4))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(mVariable, 5))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(mVariable, 6))(0,0,0));
      EXPECT_FLOAT_EQ(23.25,     (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(mVariable, 8))(0,0,0));
      EXPECT_FLOAT_EQ(21,        (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, min) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=5 stat=min"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(15, (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(15, (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(mVariable, 3))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(mVariable, 4))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(mVariable, 5))(0,0,0));
      EXPECT_FLOAT_EQ(14, (*from.getField(mVariable, 6))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(mVariable, 8))(0,0,0));
      EXPECT_FLOAT_EQ(19, (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, max) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=5 stat=max"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(mVariable, 3))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(mVariable, 4))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(mVariable, 5))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(mVariable, 6))(0,0,0));
      EXPECT_FLOAT_EQ(30, (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 8))(0,0,0));
      EXPECT_FLOAT_EQ(23, (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, std) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=5 stat=std"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(3.299832,  (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(2.947457,  (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(3.4985712, (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(4.145781,  (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(2,         (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, quantile) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=5 stat=quantile quantile=0.5"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(20,   (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20.5, (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(20,   (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(22,   (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(21,   (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, median) {
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=5 stat=median"));
      cal.calibrate(from);

      EXPECT_FLOAT_EQ(20,   (*from.getField(mVariable, 0))(0,0,0));
      EXPECT_FLOAT_EQ(20.5, (*from.getField(mVariable, 1))(0,0,0));
      EXPECT_FLOAT_EQ(20,   (*from.getField(mVariable, 2))(0,0,0));
      EXPECT_FLOAT_EQ(22,   (*from.getField(mVariable, 7))(0,0,0));
      EXPECT_FLOAT_EQ(21,   (*from.getField(mVariable, 9))(0,0,0));
   }
   TEST_F(TestCalibratorWindow, length100) {
      // A large length forces all values to be the same
      FileNetcdf from("testing/files/1x1.nc");
      CalibratorWindow cal = CalibratorWindow(mVariable ,Options("length=201 stat=mean"));
      cal.calibrate(from);

      for(int t = 0; t < 10; t++) {
         EXPECT_FLOAT_EQ(20.666666, (*from.getField(mVariable, t))(0,0,0));
      }
   }
   TEST_F(TestCalibratorWindow, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      EXPECT_DEATH(CalibratorWindow(mVariable, Options("stat=invalidStatType")), ".*");

      // Negative length
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("length=-1")), ".*");
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("length=-1 stat=quantile")), ".*");
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("length=-1 stat=quantile quantile=0.5")), ".*");

      // Missing quantile
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("stat=quantile")), ".*");
      // Invalid quantile
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("stat=quantile quantile=-0.1")), ".*");
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("stat=quantile quantile=1.1")), ".*");
      EXPECT_DEATH(CalibratorWindow(mVariable, Options("stat=quantile quantile=-999")), ".*");
   }
   TEST_F(TestCalibratorWindow, description) {
      CalibratorWindow::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
