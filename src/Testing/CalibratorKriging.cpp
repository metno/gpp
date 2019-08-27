#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Kriging.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorKriging : public ::testing::Test {
      protected:
         TestCalibratorKriging() {
         }
         virtual ~TestCalibratorKriging() {
         }
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };
   // The parameter file has the following data:
   // 0 9 9 800  -1.1
   // 0 9 8 850  -0.3
   // 0 5 5 140  4.2
   // 1 5 5 140 -5.4
   // 1 0 0 150  4
   TEST_F(TestCalibratorKriging, calcCovar) {
      CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=300000 maxElevDiff=100 efoldDist=300000"));
      EXPECT_FLOAT_EQ(0.75794888, cal.calcCovar(Location(60,10,100),Location(61,10,100)));
      EXPECT_FLOAT_EQ(0.75794888, cal.calcCovar(Location(61,10,100),Location(60,10,100)));
      EXPECT_FLOAT_EQ(0.28969184, cal.calcCovar(Location(60,10,100),Location(62,10,100)));
      EXPECT_FLOAT_EQ(0.28969184, cal.calcCovar(Location(62,10,100),Location(60,10,100)));
      EXPECT_FLOAT_EQ(0.6,        cal.calcCovar(Location(60,10,100),Location(60,10,150)));
      EXPECT_FLOAT_EQ(0.6,        cal.calcCovar(Location(60,10,150),Location(60,10,100)));
      EXPECT_FLOAT_EQ(0.4547693,  cal.calcCovar(Location(60,10,100),Location(61,10,150)));
      EXPECT_FLOAT_EQ(0.4547693,  cal.calcCovar(Location(60,10,100),Location(61,10,150)));
      EXPECT_FLOAT_EQ(0.4547693,  cal.calcCovar(Location(61,10,150),Location(60,10,100)));
      EXPECT_FLOAT_EQ(0.4547693,  cal.calcCovar(Location(60,10,150),Location(61,10,100)));
      EXPECT_FLOAT_EQ(0.4547693,  cal.calcCovar(Location(61,10,100),Location(60,10,150)));

      // Identical locations
      EXPECT_FLOAT_EQ(1, cal.calcCovar(Location(60,10,100),Location(60,10,100)));
      EXPECT_FLOAT_EQ(1, cal.calcCovar(Location(20,0,300), Location(20,0,300)));

      // Missing locations
      EXPECT_FLOAT_EQ(Util::MV, cal.calcCovar(Location(60,10,Util::MV),Location(60,10,100)));

      // Outside the radius
      EXPECT_FLOAT_EQ(0, cal.calcCovar(Location(60,10,200),Location(60,10,400)));
      EXPECT_FLOAT_EQ(0, cal.calcCovar(Location(65,10,200),Location(60,10,200)));
      EXPECT_FLOAT_EQ(0, cal.calcCovar(Location(65,10,200),Location(60,10,400)));
   }
   TEST_F(TestCalibratorKriging, calcCovarElev) {
      // Same covariance, regardless of elevation difference
      {
         CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=300000 efoldDist=300000"));
         EXPECT_FLOAT_EQ(0.75794888, cal.calcCovar(Location(60,10,100),Location(61,10,100)));
         EXPECT_FLOAT_EQ(0.75794888, cal.calcCovar(Location(61,10,100),Location(60,10,200)));
         EXPECT_FLOAT_EQ(0.28969184, cal.calcCovar(Location(60,10,100),Location(62,10,100)));
         EXPECT_FLOAT_EQ(0.28969184, cal.calcCovar(Location(62,10,100),Location(60,10,200)));
      }
      {
         CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=300000 efoldDist=300000 maxElevDiff=-999"));
         EXPECT_FLOAT_EQ(0.75794888, cal.calcCovar(Location(60,10,100),Location(61,10,100)));
         EXPECT_FLOAT_EQ(0.75794888, cal.calcCovar(Location(61,10,100),Location(60,10,200)));
         EXPECT_FLOAT_EQ(0.28969184, cal.calcCovar(Location(60,10,100),Location(62,10,100)));
         EXPECT_FLOAT_EQ(0.28969184, cal.calcCovar(Location(62,10,100),Location(60,10,200)));
      }
   }
   // Test that we don't get negative weights with Cressman
   TEST_F(TestCalibratorKriging, radiusBiggerThanEfold) {
      CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=3000000 maxElevDiff=100 efoldDist=3000"));
      EXPECT_FLOAT_EQ(0, cal.calcCovar(Location(60,10,100),Location(61,10,100)));
   }
   TEST_F(TestCalibratorKriging, 10x10) {
      FileNetcdf from("testing/files/10x10.nc");
      // Kriging when each observation only affects one grid point at a time (radius 1m)
      ParameterFile* parFile = ParameterFile::getScheme("text", Options("file=testing/files/parametersKriging.txt"));
      CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=1 maxElevDiff=100 efoldDist=200000"));

      // Prior to kriging
      FieldPtr field0 = from.getField(mVariable, 0);
      FieldPtr field1 = from.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(303,      (*field0)(5,5,0));
      EXPECT_FLOAT_EQ(306,      (*field0)(4,5,0));
      EXPECT_FLOAT_EQ(284.1962, (*field1)(5,5,0));
      EXPECT_FLOAT_EQ(283.9548,  (*field1)(0,0,0));

      // After kriging
      cal.calibrate(from, parFile);
      FieldPtr after0 = from.getField(mVariable, 0);
      FieldPtr after1 = from.getField(mVariable, 1);
      EXPECT_NEAR(307.19666,     (*after0)(5,5,0), 1e-3); // Same lat/lon, diff elev. Bias = 1*0.9992003 * 4.2
      EXPECT_FLOAT_EQ(306,       (*after0)(4,5,0)); // No change
      EXPECT_FLOAT_EQ(300,       (*after0)(0,0,0)); // No change
      EXPECT_FLOAT_EQ(278.80118, (*after1)(5,5,0)); // 0.9992003 * -5.4
      EXPECT_FLOAT_EQ(287.86127, (*after1)(0,0,0)); // 0.9760893 * 4
   }
   /*
   TEST_F(TestCalibratorKriging, radius) {
      {
         FileNetcdf from("testing/files/10x10.nc");
         CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=120000 efoldDist=2e10 maxElevDiff=1000 parameters=testing/files/parametersKriging.txt fileType=textSpatial"));

         cal.calibrate(from, parFile);
         FieldPtr after0 = from.getField(mVariable, 0);
         FieldPtr after1 = from.getField(mVariable, 1);

         // Near 5,5
         EXPECT_FLOAT_EQ(302, (*after0)(3,5,0));
         EXPECT_NEAR(310.2,   (*after0)(4,5,0), 1e-3);
         EXPECT_NEAR(307.2,   (*after0)(5,5,0), 1e-3);
         EXPECT_NEAR(306.2,   (*after0)(6,5,0), 1e-3);
         EXPECT_FLOAT_EQ(302, (*after0)(7,5,0));
         EXPECT_FLOAT_EQ(302, (*after0)(5,3,0));
         EXPECT_NEAR(314.2,   (*after0)(5,4,0), 1e-3);
         EXPECT_NEAR(305.2,   (*after0)(5,6,0), 1e-3);
         EXPECT_FLOAT_EQ(300, (*after0)(5,7,0));

         // Near 9,9 (influenced by bias at 9.8 and 9,9
         EXPECT_FLOAT_EQ(306, (*after0)(7,9,0));       // No influence
         EXPECT_NEAR(303.3,   (*after0)(8,9,0), 1e-3); // Influenced by -1.1 and -0.3
         EXPECT_NEAR(303.3,   (*after0)(9,9,0), 1e-3);
         EXPECT_NEAR(301.3,   (*after0)(9,8,0), 1e-3);
         EXPECT_NEAR(300.7,   (*after0)(9,7,0), 1e-3); // Influenced by -0.3 only
      }
      {
         FileNetcdf from("testing/files/10x10.nc");
         CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=0 efoldDist=2e10 parameters=testing/files/parametersKriging.txt fileType=textSpatial"));

         cal.calibrate(from, parFile);
         FieldPtr after0 = from.getField(mVariable, 0);
         FieldPtr after1 = from.getField(mVariable, 1);

         EXPECT_FLOAT_EQ(302, (*after0)(3,5,0));
         EXPECT_FLOAT_EQ(306, (*after0)(4,5,0));
         EXPECT_NEAR(307.2,   (*after0)(5,5,0), 1e-3);
         EXPECT_FLOAT_EQ(302, (*after0)(6,5,0));
         EXPECT_FLOAT_EQ(302, (*after0)(7,5,0));
         EXPECT_FLOAT_EQ(302, (*after0)(5,3,0));
         EXPECT_FLOAT_EQ(310, (*after0)(5,4,0));
         EXPECT_FLOAT_EQ(301, (*after0)(5,6,0));
         EXPECT_FLOAT_EQ(300, (*after0)(5,7,0));
      }
   }
   TEST_F(TestCalibratorKriging, aux) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorKriging cal = CalibratorKriging(mVariable, Options("radius=1 maxElevDiff=1000 efoldDist=1e10 auxVariable=Precip range=0,1 parameters=testing/files/parametersKriging.txt fileType=textSpatial"));

      cal.calibrate(from, parFile);
      FieldPtr after0 = from.getField(mVariable, 0);
      EXPECT_NEAR(307.2,   (*after0)(5,5,0), 1e-3); // Bias should be 4.2, and precip is < 1.
      EXPECT_FLOAT_EQ(310, (*after0)(5,4,0)); // Precip is above 1 so correction should be off
      EXPECT_FLOAT_EQ(304, (*after0)(9,9,0)); // Precip is above 1 so correction should be off
      EXPECT_FLOAT_EQ(304, (*after0)(8,9,0)); // Precip is above 1 so correction should be off
   }
   */
   TEST_F(TestCalibratorKriging, valid) {
      CalibratorKriging(mVariable, Options("radius=100 maxElevDiff=100 efoldDist=2"));
      CalibratorKriging(mVariable, Options("maxElevDiff=100 efoldDist=2"));
      CalibratorKriging(mVariable, Options("radius=100 efoldDist=2"));
      CalibratorKriging(mVariable, Options("radius=100 maxElevDiff=100"));
      CalibratorKriging(mVariable, Options(""));
      CalibratorKriging(mVariable, Options("radius=0"));
      CalibratorKriging(mVariable, Options("maxElevDiff=0"));
      CalibratorKriging(mVariable, Options("efoldDist=0"));
      CalibratorKriging(mVariable, Options("efoldDist=0 operator=add"));
      CalibratorKriging(mVariable, Options("efoldDist=0 operator=multiply"));
   }
   TEST_F(TestCalibratorKriging, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // Negative values
      EXPECT_DEATH(CalibratorKriging(mVariable, Options("radius=-1 maxElevDiff=100 efoldDist=2")), ".*");
      EXPECT_DEATH(CalibratorKriging(mVariable, Options("radius=100 maxElevDiff=-1 efoldDist=2")), ".*");
      EXPECT_DEATH(CalibratorKriging(mVariable, Options("radius=100 maxElevDiff=-100 efoldDist=-1")), ".*");
      EXPECT_DEATH(CalibratorKriging(mVariable, Options("radius=100 maxElevDiff=-100 efoldDist=2")), ".*");

      // Invalid operator
      EXPECT_DEATH(CalibratorKriging(mVariable, Options("radius=100 maxElevDiff=100 efoldDist=2 operator=nonvalidOperator")), ".*");
   }
   TEST_F(TestCalibratorKriging, description) {
      CalibratorKriging::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
