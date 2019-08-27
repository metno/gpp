#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Qq.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorQq : public ::testing::Test {
      protected:
         TestCalibratorQq() {
         }
         virtual ~TestCalibratorQq() {
         }
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
         Parameters getParameters(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8) {
            std::vector<float> parValues(8, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            parValues[2] = a3;
            parValues[3] = a4;
            parValues[4] = a5;
            parValues[5] = a6;
            parValues[6] = a7;
            parValues[7] = a8;
            return Parameters (parValues);
         }
         ParameterFileSimple getParameterFile(float a1, float a2, float a3, float a4, float a5, float a6, float a7, float a8) {
            ParameterFileSimple parFile(getParameters(a1, a2, a3, a4, a5, a6, a7, a8));
            return parFile;
         }
         ObsEns getObsEns(float obs, float ens0, float ens1=Util::MV, float ens2=Util::MV) {
            std::vector<float> ens;
            ens.push_back(ens0);
            if(Util::isValid(ens1))
               ens.push_back(ens1);
            if(Util::isValid(ens2))
               ens.push_back(ens2);
            return ObsEns(obs, ens);
         }
   };
   TEST_F(TestCalibratorQq, 1to1) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,260,300,290,303,300,315);
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1"));

      cal.calibrate(from, &parFile);

      FieldPtr after = from.getField(mVariable, 0);
      FieldPtr after1 = from.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(270, (*after)(5,2,0));      // 301 (row,col)
      EXPECT_FLOAT_EQ(290.8333, (*after)(5,9,0)); // 304
      EXPECT_FLOAT_EQ(305, (*after)(0,9,0));      // 320 // Outside

      EXPECT_FLOAT_EQ(244.1962, (*after1)(5,5,0)); // 284.1962 // Outside
   }
   TEST_F(TestCalibratorQq, meanSlope) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,260,300,290,303,300,315);
      CalibratorQq cal(mVariable ,Options("extrapolation=meanSlope"));

      cal.calibrate(from, &parFile);

      FieldPtr after = from.getField(mVariable, 0);
      FieldPtr after1 = from.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(310, (*after)(0,9,0)); // mean slope = 50/25
      EXPECT_FLOAT_EQ(238.3924, (*after1)(5,5,0)); // 250 - slope*(5.8038)
   }
   TEST_F(TestCalibratorQq, nearestSlope) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,260,300,290,303,300,315);
      CalibratorQq cal(mVariable ,Options("extrapolation=nearestSlope"));

      cal.calibrate(from, &parFile);

      FieldPtr after = from.getField(mVariable, 0);
      FieldPtr after1 = from.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(304.16666667, (*after)(0,9,0)); 
      EXPECT_FLOAT_EQ(244.1962, (*after1)(5,5,0)); // nearest slope = 1
   }
   TEST_F(TestCalibratorQq, zero) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,260,300,290,303,300,315);
      CalibratorQq cal(mVariable ,Options("extrapolation=zero"));

      cal.calibrate(from, &parFile);

      FieldPtr after = from.getField(mVariable, 0);
      FieldPtr after1 = from.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(300, (*after)(0,9,0));
      EXPECT_FLOAT_EQ(250, (*after1)(5,5,0));
   }
   TEST_F(TestCalibratorQq, twoEqual) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,270,301,280,301,320,320);
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1"));

      cal.calibrate(from, &parFile);

      FieldPtr after = from.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(275, (*after)(5,2,0));
   }
   TEST_F(TestCalibratorQq, twoEqualObs) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,280,300,290,302,320,320);
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1"));

      cal.calibrate(from, &parFile);

      FieldPtr after = from.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(285, (*after)(5,2,0)); // raw = 301
   }
   TEST_F(TestCalibratorQq, missingValues) {
      FileNetcdf from("testing/files/10x10.nc");
      ParameterFileSimple parFile = getParameterFile(250,290,260,300,290,303,300,315);
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1"));

      FieldPtr field = from.getField(mVariable, 0);
      (*field)(5,5,0) = Util::MV;
      cal.calibrate(from, &parFile);

      // Can't calibrate
      EXPECT_FLOAT_EQ(Util::MV, (*field)(5,5,0));
      // Unaffected by other missing values
      EXPECT_FLOAT_EQ(270, (*field)(5,2,0));
   }
   TEST_F(TestCalibratorQq, train) {
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1"));
      std::vector<ObsEns> data;
      data.push_back(getObsEns(3, 2));
      data.push_back(getObsEns(4, 3));
      Parameters parameters = cal.train(data);

      ASSERT_EQ(4, parameters.size());
      EXPECT_FLOAT_EQ(3, parameters[0]);
      EXPECT_FLOAT_EQ(2, parameters[1]);
      EXPECT_FLOAT_EQ(4, parameters[2]);
      EXPECT_FLOAT_EQ(3, parameters[3]);
   }
   TEST_F(TestCalibratorQq, trainQuantiles) {
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1 quantiles=0,0.32,0.59,1"));
      std::vector<ObsEns> data;
      data.push_back(getObsEns(4, 2));
      data.push_back(getObsEns(2, 3));
      data.push_back(getObsEns(1, 3));
      data.push_back(getObsEns(5, 6));
      data.push_back(getObsEns(0, 5));
      Parameters parameters = cal.train(data);

      // Sorted obs:  0 1 2 4 5
      // Sorted fcst: 2 3 3 5 6
      ASSERT_EQ(8, parameters.size());
      EXPECT_FLOAT_EQ(0,    parameters[0]); // 0%
      EXPECT_FLOAT_EQ(2,    parameters[1]);
      EXPECT_FLOAT_EQ(1.28, parameters[2]); // 32%
      EXPECT_FLOAT_EQ(3,    parameters[3]);
      EXPECT_FLOAT_EQ(2.72, parameters[4]); // 59%
      EXPECT_FLOAT_EQ(3.72, parameters[5]);
      EXPECT_FLOAT_EQ(5,    parameters[6]); // 100 %
      EXPECT_FLOAT_EQ(6,    parameters[7]);
   }
   TEST_F(TestCalibratorQq, trainQuantilesWithExtra) {
      CalibratorQq cal(mVariable ,Options("extrapolation=1to1 quantiles=0,0.25,1 extraObs=-10,3 extraFcst=-5,1"));
      std::vector<ObsEns> data;
      data.push_back(getObsEns(4, 2));
      data.push_back(getObsEns(2, 3));
      data.push_back(getObsEns(1, 3));
      data.push_back(getObsEns(5, 6));
      data.push_back(getObsEns(0, 5));
      Parameters parameters = cal.train(data);

      // Sorted obs:  0 1 2 4 5
      // Sorted fcst: 2 3 3 5 6
      // With extra obs: -10 0 1 3 5
      // With extra fcst: -5 1 2 3 6
      ASSERT_EQ(10, parameters.size());
      EXPECT_FLOAT_EQ(-10,  parameters[0]);
      EXPECT_FLOAT_EQ(-5,   parameters[1]);
      EXPECT_FLOAT_EQ(0,    parameters[2]);
      EXPECT_FLOAT_EQ(1,    parameters[3]);
      EXPECT_FLOAT_EQ(1,    parameters[4]);
      EXPECT_FLOAT_EQ(2,    parameters[5]);
      EXPECT_FLOAT_EQ(3,    parameters[6]);
      EXPECT_FLOAT_EQ(3,    parameters[7]);
      EXPECT_FLOAT_EQ(5,    parameters[8]);
      EXPECT_FLOAT_EQ(6,    parameters[9]);
   }
   TEST_F(TestCalibratorQq, unevenNumberOfParameters) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      FileNetcdf from("testing/files/10x10.nc");
      std::vector<float> parValues;
      parValues.resize(7,3.2);
      ParameterFileSimple parFile(parValues);

      CalibratorQq calibrator(mVariable, Options("extrapolation=1to1"));
      EXPECT_DEATH(calibrator.calibrate(from, &parFile), ".*");
   }
   TEST_F(TestCalibratorQq, invalidConstruction) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // Invalid extrapolation
      EXPECT_DEATH(CalibratorQq(mVariable ,Options("extrapolation=92090j0f923")), ".*");
   }
   TEST_F(TestCalibratorQq, description) {
      CalibratorQq::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
