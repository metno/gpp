#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Sort.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorSort : public ::testing::Test {
      protected:
         TestCalibratorSort() {
         }
         virtual ~TestCalibratorSort() {
         }
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
         void test(CalibratorSort& cal, File& file, float i0, float i1, float i2, float e0, float e1, float e2) {
            FieldPtr field  = file.getField(mVariable, 0);
            (*field)(0,0,0) = i0;
            (*field)(0,0,1) = i1;
            (*field)(0,0,2) = i2;
            cal.calibrate(file, NULL);
            const FieldPtr after  = file.getField(mVariable, 0);

            EXPECT_FLOAT_EQ(e0, (*after)(0,0,0));
            EXPECT_FLOAT_EQ(e1, (*after)(0,0,1));
            EXPECT_FLOAT_EQ(e2, (*after)(0,0,2));
         }
   };
   TEST_F(TestCalibratorSort, simple) {
      CalibratorSort cal = CalibratorSort(mVariable, Options(""));
      FileFake file(Options("nLat=1 nLon=1 nEns=3 nTime=1"));
      //              Before sort   After sort
      test(cal, file, 3,1,2,        1,2,3);
      test(cal, file, 1,1,2,        1,1,2);
      test(cal, file, 3,1,1,        1,1,3);
      test(cal, file, 3,Util::MV,2, 2,3,Util::MV);
      test(cal, file, 2,Util::MV,2, 2,2,Util::MV);
      test(cal, file, Util::MV,Util::MV,Util::MV, Util::MV,Util::MV,Util::MV);
      test(cal, file, Util::MV,1,Util::MV,   1,Util::MV,Util::MV);
   }
   TEST_F(TestCalibratorSort, description) {
      CalibratorSort::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
