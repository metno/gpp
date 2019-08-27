#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Phase.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorPhase : public ::testing::Test {
      protected:
         TestCalibratorPhase() {
         }
         virtual ~TestCalibratorPhase() {
         }
         virtual void SetUp() {
         }
         virtual void TearDown() {
         }
         Parameters getParameters(float a1, float a2) {
            std::vector<float> parValues(8, 0);
            parValues[0] = a1;
            parValues[1] = a2;
            return Parameters (parValues);
         }
         ParameterFileSimple getParameterFile(float a1, float a2) {
            std::vector<float> values;
            values.push_back(a1);
            values.push_back(a2);
            ParameterFileSimple parFile(values);
            return parFile;
         }
         CalibratorPhase getCalibratorWetbulb(Variable iVariable) {
            return CalibratorPhase(iVariable, Options("temperature=air_temperature_2m precipitation=precipitation_amount rh=relative_humidity_2m pressure=surface_air_pressure"));
         }
         CalibratorPhase getCalibrator(Variable iVariable) {
            return CalibratorPhase(iVariable, Options("temperature=air_temperature_2m precipitation=precipitation_amount"));
         }
         void setValues(const File& iFile, float iPrecip, float iTemp) {
            FieldPtr precip   = iFile.getField(Variable("precipitation_amount"), 0);
            FieldPtr temp     = iFile.getField(Variable("air_temperature_2m"), 0);
            (*precip)(0,0,0) = iPrecip;
            (*temp)(0,0,0) = iTemp;
         }
   };

   TEST_F(TestCalibratorPhase, 10x10) {
      FileNetcdf file("testing/files/10x10.nc");
      ParameterFileText parFile(Options("file=testing/files/parametersPhase.txt"));
      Parameters parameters = parFile.getParameters(0);
      ASSERT_EQ(2, parameters.size());
      EXPECT_FLOAT_EQ(273.7, parameters[0]);
      EXPECT_FLOAT_EQ(274.7, parameters[1]);
      CalibratorPhase cal = getCalibratorWetbulb(Variable("phase"));

      cal.calibrate(file, &parFile);
      FieldPtr phase    = file.getField(Variable("phase"), 0);
      FieldPtr precip   = file.getField(Variable("precipitation_amount"), 0);
      FieldPtr temp     = file.getField(Variable("air_temperature_2m"), 0);
      FieldPtr rh       = file.getField(Variable("relative_humidity_2m"), 0);
      FieldPtr pressure = file.getField(Variable("surface_air_pressure"), 0);
      // T      301 K
      // RH     0.95925
      // P      98334 pa
      // Precip 1.38 mm
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseRain, (*phase)(2,5,0));
   }
   TEST_F(TestCalibratorPhase, phases) {
      FileFake file(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      //file.initNewVariable(Variable("phase"));
      //file.initNewVariable(Variable("precipitation_amount"));
      //file.initNewVariable(Variable("air_temperature_2m"));
      //file.initNewVariable(Variable("relative_humidity_2m"));
      //file.initNewVariable(Variable("surface_air_pressure"));
      FieldPtr phase = file.getField(Variable("phase"), 0);
      ParameterFileSimple parFile = getParameterFile(273.7,274.7);
      CalibratorPhase cal = getCalibratorWetbulb(Variable("phase"));
      float wetbulb = CalibratorDiagnoseHumidity::computeWetbulb(270, 101325, 0.95);
      setValues(file, 5, wetbulb);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSnow, (*phase)(0,0,0));

      wetbulb = CalibratorDiagnoseHumidity::computeWetbulb(274, 101325, 0.98);
      setValues(file, 5, wetbulb);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSleet, (*phase)(0,0,0));

      wetbulb = CalibratorDiagnoseHumidity::computeWetbulb(274, 101325, 0.5);
      setValues(file, 0, wetbulb);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseNone, (*phase)(0,0,0));

      wetbulb = CalibratorDiagnoseHumidity::computeWetbulb(275, 101325, 0.5);
      setValues(file, 5, wetbulb);
      cal.calibrate(file, &parFile);
      EXPECT_FLOAT_EQ(CalibratorPhase::PhaseSnow, (*phase)(0,0,0));
   }
   TEST_F(TestCalibratorPhase, description) {
      CalibratorPhase::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
