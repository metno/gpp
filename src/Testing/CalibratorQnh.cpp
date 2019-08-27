#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Qnh.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorQnh : public ::testing::Test {
      protected:
         TestCalibratorQnh() {
         }
         virtual ~TestCalibratorQnh() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
   };
   TEST_F(TestCalibratorQnh, 10x10) {
      Variable variableQnh = Variable("qnh");
      Variable variableP = Variable("surface_air_pressure");
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorQnh cal = CalibratorQnh(variableQnh, Options("pressureVariable=surface_air_pressure"));
      FieldPtr p = from.getField(variableP, 0);

      for(int i = 0; i < from.getNumTime(); i++)
         from.addField(from.getEmptyField(), variableQnh, i);
      cal.calibrate(from);
      EXPECT_FLOAT_EQ(98334.44, (*p)(5,2,0));
      FieldPtr qnh = from.getField(variableQnh, 0);
      ASSERT_EQ(10, qnh->getNumY());
      ASSERT_EQ(10, qnh->getNumX());
      ASSERT_EQ(1,  qnh->getNumEns());

      // Altitude: 159.6324 Pressure: 98334.44
      EXPECT_FLOAT_EQ(100220.6455, (*qnh)(5,2,0));
   }
   TEST_F(TestCalibratorQnh, calcQnh) {
      EXPECT_FLOAT_EQ(100000, CalibratorQnh::calcQnh(0,100000));
      EXPECT_FLOAT_EQ(0, CalibratorQnh::calcQnh(0,0));

      EXPECT_FLOAT_EQ(100184.6424, CalibratorQnh::calcQnh(100,99000));
      EXPECT_FLOAT_EQ(97826.7259, CalibratorQnh::calcQnh(-100,99000));
      EXPECT_FLOAT_EQ(0, CalibratorQnh::calcQnh(-100,0));
   }
   TEST_F(TestCalibratorQnh, description) {
      CalibratorQnh::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
