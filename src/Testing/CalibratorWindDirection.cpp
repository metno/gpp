#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Parameters.h"
#include "../Calibrator/WindDirection.h"
#include <gtest/gtest.h>
#include <cmath>

namespace {
   class TestCalibratorWindDirection : public ::testing::Test {
      protected:
         TestCalibratorWindDirection() {
         }
         virtual ~TestCalibratorWindDirection() {
         }
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
         Parameters getParameters(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8, float a9) {
            std::vector<float> parValues(9, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            parValues[2] = a3;
            parValues[3] = a4;
            parValues[4] = a5;
            parValues[5] = a6;
            parValues[6] = a7;
            parValues[7] = a8;
            parValues[8] = a9;
            return Parameters (parValues);
         }
         ParameterFileSimple getParameterFile(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8, float a9) {
            ParameterFileSimple parFile(getParameters(a1, a2, a3, a4, a5, a6, a7, a8, a9));
            return parFile;
         }
   };
   /*
    * TODO: Diagnosing of wind direction is not implemented, therefore this cannot be tested
   TEST_F(TestCalibratorWindDirection, getFactor) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      ParameterFileSimple parFile = getParameterFile(1,0,0,0,0,0,0,0,0);
      CalibratorWindDirection cal(mVariable, Options());
      Parameters par = parFile.getParameters(0);

      EXPECT_FLOAT_EQ(1, cal.getFactor(0, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(45, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(90, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(180, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(270, par));

      parFile = getParameterFile(1,0,1,0,0,0,0,0,0);
      par = parFile.getParameters(0);

      EXPECT_FLOAT_EQ(2, cal.getFactor(0, par));
      EXPECT_FLOAT_EQ(1.7071067, cal.getFactor(45, par));
      EXPECT_FLOAT_EQ(1, cal.getFactor(90, par));
      EXPECT_NEAR(0, cal.getFactor(180, par), 1e-6);
      EXPECT_FLOAT_EQ(1, cal.getFactor(270, par));

      parFile = getParameterFile(1,1,0.5,0,0,0,0,0,0);
      par = parFile.getParameters(0);

      EXPECT_FLOAT_EQ(1.5, cal.getFactor(0, par));
      EXPECT_FLOAT_EQ(2.06066, cal.getFactor(45, par));
      EXPECT_FLOAT_EQ(2, cal.getFactor(90, par));
      EXPECT_FLOAT_EQ(0.5, cal.getFactor(180, par));
      EXPECT_NEAR(0, cal.getFactor(270, par), 1e-6);
   }
   TEST_F(TestCalibratorWindDirection, getFactorInvalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      // Test that getFactor does not give negative values
      ParameterFileSimple parFile = getParameterFile(1,2,3,0,0,0,0,0,0);
      CalibratorWindDirection cal(mVariable, Options());
      Parameters par = parFile.getParameters(0);
      for(int dir = 0; dir <= 360; dir = dir + 15) {
         EXPECT_GE(cal.getFactor(dir, par), 0);
      }

      // Test that getFactor returns invalid factor when parameters are invalid
      parFile = getParameterFile(1,2,3,0,0,0,Util::MV,0,0);
      par = parFile.getParameters(0);
      for(int dir = 0; dir <= 360; dir = dir + 15) {
         EXPECT_FLOAT_EQ(Util::MV, cal.getFactor(dir, par));
      }

      // Invalid values
      // EXPECT_DEATH(cal.getFactor(-1, par), ".*");
   }
   */
   TEST_F(TestCalibratorWindDirection, description) {
      CalibratorWindDirection::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
