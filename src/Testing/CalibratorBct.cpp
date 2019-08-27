#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Bct.h"
#include <gtest/gtest.h>
#include <boost/math/distributions/students_t.hpp>

namespace {
   class TestCalibratorBct : public ::testing::Test {
      protected:
         TestCalibratorBct() {
         }
         virtual ~TestCalibratorBct() {
         }
         void reset10x10() const {
            Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
         };
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         Variable mVariable = Variable("precipitation_amount");

         void test(CalibratorBct& cal, File& file, ParameterFile& parFile, float i0, float i1, float i2, float e0, float e1, float e2) {
            FieldPtr field  = file.getField(mVariable, 0);
            (*field)(0,0,0) = i0;
            (*field)(0,0,1) = i1;
            (*field)(0,0,2) = i2;
            cal.calibrate(file, &parFile);
            const FieldPtr after  = file.getField(mVariable, 0);

            EXPECT_FLOAT_EQ(e0, (*after)(0,0,0));
            EXPECT_FLOAT_EQ(e1, (*after)(0,0,1));
            EXPECT_FLOAT_EQ(e2, (*after)(0,0,2));
         }
         Parameters getParameters(float a1, float a2, float a3, float a4, float a5, float a6, float a7) {
            std::vector<float> parValues(8, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            parValues[2] = a3;
            parValues[3] = a4;
            parValues[4] = a5;
            parValues[5] = a6;
            parValues[6] = a7;
            return Parameters (parValues);
         }
         ParameterFileSimple getParameterFile(float a1, float a2, float a3, float a4, float a5, float a6, float a7) {
            ParameterFileSimple parFile(getParameters(a1, a2, a3, a4, a5, a6, a7));
            return parFile;
         }
         CalibratorBct getCalibrator(const Options& iOptions=Options("")) {
            return CalibratorBct(mVariable, iOptions);
         }
   };

   TEST_F(TestCalibratorBct, small) {
      /*
      // Set up file
      FileFake file(1, 1, 3, 1);
      // Set up calibrator
      ParameterFileSimple parFile = getParameterFile(-1.1,1.4,0.05,-0.05, 2.03, -0.05, 0.82);
      CalibratorBct cal = getCalibrator(Options(""));

      // High wind case: mean = 5.667 std = 0
      test(cal, file, parFile, 12, 1, 4, 6.6280001, 1.0158194, 3.0753615);
      // Low wind case: mean = 0.3333 std = 0.3333
      test(cal, file, parFile, 0.1, 0, 0.9, 0, 0, 0.59979528);
      */

   }
   // Check that boost can handle large values of the number of degrees of freedom
   // in the t-distribution
   TEST_F(TestCalibratorBct, tdistribution) {
      float par = 10;
      double tau = exp(par);
      boost::math::students_t_distribution<> dist(tau);
      float z = boost::math::quantile(dist, 0.001);
      EXPECT_TRUE(Util::isValid(z));
   }
   TEST_F(TestCalibratorBct, description) {
      CalibratorBct::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
