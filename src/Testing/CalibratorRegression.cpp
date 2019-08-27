#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Regression.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorRegression : public ::testing::Test {
      protected:
         TestCalibratorRegression() {
         }
         virtual ~TestCalibratorRegression() {
         }
         void reset10x10() const {
            Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
         };
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
         std::vector<float> getVector(float iArray[]) {
            return std::vector<float>(iArray, iArray + sizeof(iArray)/sizeof(float));
         }
   };
   // Constant correction to 0.3
   TEST_F(TestCalibratorRegression, 10x10_0order) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/regression0order.txt"));
      CalibratorRegression cal = CalibratorRegression(mVariable, Options());

      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(0.3, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(0.3, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(0.3, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorRegression, 10x10_1order) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/regression1order.txt"));
      CalibratorRegression cal = CalibratorRegression(mVariable, Options());

      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(361.5, (*after)(5,2,0)); // 0.3 + 1.2*301
      EXPECT_FLOAT_EQ(365.1, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(384.3, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorRegression, 10x10_multivariate) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/regression1order.txt"));
      CalibratorRegression cal = CalibratorRegression(mVariable, Options("variables=relative_humidity_2m,air_temperature_2m"));

      FieldPtr before = from.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(301, (*before)(5,2,0));
      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(361.48777506, (*after)(5,2,0)); // 0.3*0.9592502 + 1.2*301
   }
   TEST_F(TestCalibratorRegression, 10x10_2order) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/regression2order.txt"));
      CalibratorRegression cal = CalibratorRegression(mVariable, Options());

      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(-72174.08, (*after)(5,2,0)); // -0.3 + 1.02*301 - 0.8*301^2
      EXPECT_FLOAT_EQ(-73623.02, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(-81593.90, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorRegression, missing_parameters) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileText par(Options("file=testing/files/regressionMissing.txt"));
      CalibratorRegression cal = CalibratorRegression(mVariable, Options());

      cal.calibrate(from, &par);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(Util::MV, (*after)(5,2,0)); // 0.3 + 1.2*301
      EXPECT_FLOAT_EQ(Util::MV, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(Util::MV, (*after)(0,9,0));
   }
   // Incorrect number of data columns
   TEST_F(TestCalibratorRegression, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
	  FileNetcdf from("testing/files/10x10.nc");
      Util::setShowError(false);
      ParameterFileText par(Options("file=testing/files/regressionInvalid1.txt"));
      CalibratorRegression calibrator(mVariable, Options());
      EXPECT_DEATH(calibrator.calibrate(from, &par), ".*");
   }
   // Missing parameter file
   TEST_F(TestCalibratorRegression, invalid2) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
	  FileNetcdf from("testing/files/10x10.nc");
      Util::setShowError(false);
      CalibratorRegression calibrator(mVariable, Options());
      EXPECT_DEATH(calibrator.calibrate(from, NULL), ".*");
   }
   TEST_F(TestCalibratorRegression, training) {
      // R code:
      // data = data.frame(x=c(4,5,9),y=c(3.2,5,14))
      FileNetcdf from("testing/files/10x10.nc");
      std::vector<ObsEns> obsens;
      Ens ens(3,0);
      ens[0] = 4;
      ens[1] = 4;
      ens[2] = 4;
      obsens.push_back(ObsEns(3.2, ens));
      ens[0] = 6;
      ens[1] = 5;
      ens[2] = 4;
      obsens.push_back(ObsEns(5, ens));
      ens[0] = 9;
      ens[1] = 8;
      ens[2] = 10;
      obsens.push_back(ObsEns(14, ens));

      // glm(y ~ x)$coefficients
      CalibratorRegression cal = CalibratorRegression(mVariable, Options("order=1"));
      Parameters par = cal.train(obsens);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(-5.7142825, par[0]);
      EXPECT_FLOAT_EQ(2.185714286, par[1]);

      // glm(y ~ x-1)$coefficients
      cal = CalibratorRegression(mVariable, Options("order=1 intercept=0"));
      par = cal.train(obsens);
      ASSERT_EQ(2, par.size());
      EXPECT_FLOAT_EQ(0, par[0]);
      EXPECT_FLOAT_EQ(1.342623, par[1]);
   }
   TEST_F(TestCalibratorRegression, description) {
      CalibratorRegression::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
