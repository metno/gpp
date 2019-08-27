#include "../File/Netcdf.h"
#include "../Util.h"
#include "../Calibrator/Calibrator.h"
#include <gtest/gtest.h>

// For each test it is safe to assume that 10x10_copy.nc is identical to 10x10.nc
// After the test is done, it is safe to assume that 10x10_copy.nc is again reverted.
namespace {
   class FileNetcdfTest : public ::testing::Test {
      protected:
         void reset10x10() const {
            Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
         };
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
             reset10x10();
         }
         Variable mVariable;
   };

   TEST_F(FileNetcdfTest, isValid) {
      FileNetcdf file = FileNetcdf("testing/files/validNetcdf1.nc");
   }
   TEST_F(FileNetcdfTest, missingY) {
      FileNetcdf file = FileNetcdf("testing/files/validNetcdf2.nc");
      EXPECT_EQ(1, file.getNumY());
      EXPECT_EQ(10, file.getNumX());
      EXPECT_EQ(10, file.getNumEns());
      EXPECT_EQ(2, file.getNumTime());
   }
   TEST_F(FileNetcdfTest, missingTime) {
      FileNetcdf file = FileNetcdf("testing/files/validNetcdf3.nc");
      EXPECT_EQ(3, file.getNumY());
      EXPECT_EQ(3, file.getNumX());
      EXPECT_EQ(1, file.getNumEns());
      EXPECT_EQ(1, file.getNumTime());
   }
   TEST_F(FileNetcdfTest, missingXandTime) {
      FileNetcdf file = FileNetcdf("testing/files/validNetcdf4.nc");
      EXPECT_EQ(10, file.getNumY());
      EXPECT_EQ(1, file.getNumX());
      EXPECT_EQ(1, file.getNumEns());
      EXPECT_EQ(1, file.getNumTime());
      FieldPtr field = file.getField(Variable("air_temperature_2m"), 0);
      EXPECT_FLOAT_EQ(21, (*field)(0, 0, 0));
      EXPECT_FLOAT_EQ(26, (*field)(5, 0, 0));
   }
   TEST_F(FileNetcdfTest, dimNames) {
      // Check that dim and var options are passed to the parser
      // Altitudes have the x and y dimensions reversed, so test that this works.
      FileNetcdf file = FileNetcdf("testing/files/validNetcdfDimNames.nc", Options("xDim=h2 yDim=h1 timeDim=date ensDim=member latVar=latVar lonVar=lonVar timeVar=date"));
      EXPECT_EQ(3, file.getNumY());
      EXPECT_EQ(2, file.getNumX());
      EXPECT_EQ(2, file.getNumEns());
      EXPECT_EQ(2, file.getNumTime());
      vec2 lats = file.getLats();
      vec2 lons = file.getLons();
      vec2 elevs = file.getElevs();
      for(int i = 0; i < file.getNumY(); i++) {
         for(int j = 0; j < file.getNumX(); j++) {
            EXPECT_FLOAT_EQ(i, lats[i][j]);
            EXPECT_FLOAT_EQ(j, lons[i][j]);
         }
      }
      // Ncview doesn't show the correct elevation field since the field has the y/x
      // dimensions flipped.
      EXPECT_FLOAT_EQ(160, elevs[0][0]);
      EXPECT_FLOAT_EQ(295, elevs[1][0]);
      EXPECT_FLOAT_EQ(11, elevs[2][0]);
      EXPECT_FLOAT_EQ(-13, elevs[0][1]);
      EXPECT_FLOAT_EQ(168, elevs[1][1]);
      EXPECT_FLOAT_EQ(-171, elevs[2][1]);

      FieldPtr field = file.getField(Variable("air_temperature_2m"), 0);
      EXPECT_FLOAT_EQ(1, (*field)(0, 0, 0));
      EXPECT_FLOAT_EQ(27, (*field)(2, 0, 0));

      EXPECT_FLOAT_EQ(28, (*field)(2, 0, 1));
      EXPECT_FLOAT_EQ(32, (*field)(2, 1, 1));
      field = file.getField(Variable("air_temperature_2m"), 1);
      EXPECT_FLOAT_EQ(21, (*field)(1, 0, 0));
      EXPECT_FLOAT_EQ(24, (*field)(1, 1, 0));
      EXPECT_FLOAT_EQ(12, (*field)(0, 1, 1));
      EXPECT_FLOAT_EQ(38, (*field)(2, 1, 1));
   }
   TEST_F(FileNetcdfTest, geopotential) {
      // Test that altitude is computed from surface geopotential
      FileNetcdf file = FileNetcdf("testing/files/validNetcdfGeopotential.nc");
      vec2 elevs = file.getElevs();
      EXPECT_FLOAT_EQ(90/9.81, elevs[0][0]);
      EXPECT_FLOAT_EQ(30/9.81, elevs[1][0]);
      EXPECT_FLOAT_EQ(14/9.81, elevs[2][0]);
      EXPECT_FLOAT_EQ(80/9.81, elevs[0][1]);
      EXPECT_FLOAT_EQ(40/9.81, elevs[1][1]);
      EXPECT_FLOAT_EQ(99/9.81, elevs[2][1]);
   }

   TEST_F(FileNetcdfTest, analysis) {
      // Test that an analysis file is parsed correctly
      FileNetcdf file = FileNetcdf("testing/files/validNetcdfAnalysis.nc");
      EXPECT_EQ(1, file.getNumTime());
      std::vector<double> times = file.getTimes();
      EXPECT_EQ(1, times.size());
      // Use forecast_reference_time
      EXPECT_FLOAT_EQ(1414130400, times[0]);

      FieldPtr field = file.getField(Variable("air_temperature_2m"), 0);
      EXPECT_FLOAT_EQ(300, (*field)(0, 0, 0));
      EXPECT_FLOAT_EQ(303, (*field)(2, 1, 0));
      EXPECT_FLOAT_EQ(307, (*field)(2, 0, 1));
      EXPECT_FLOAT_EQ(Util::MV, (*field)(0, 0, 1));
   }

   TEST_F(FileNetcdfTest, scalarTime) {
      // Test that an analysis file can use a time variable without a dimension
      FileNetcdf file = FileNetcdf("testing/files/validNetcdfAnalysis2.nc");
      EXPECT_EQ(1, file.getNumTime());
      std::vector<double> times = file.getTimes();
      EXPECT_FLOAT_EQ(1414130400, times[0]);
   }

   TEST_F(FileNetcdfTest, overwriteAttribute) {
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history", "test512");
      EXPECT_EQ("test512", file.getGlobalAttribute("history"));
   }
   TEST_F(FileNetcdfTest, addAttribute) {
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history2", "test123");
      EXPECT_EQ("test123", file.getGlobalAttribute("history2"));
   }
   TEST_F(FileNetcdfTest, missingAttribute) {
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      std::string att = file.getGlobalAttribute("qowhoiqfhoiqhdow");
      EXPECT_EQ("", att);
   }
   TEST_F(FileNetcdfTest, appendAttribute) {
      // Check that appending and prepending works
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("history", "empty");
      file.prependGlobalAttribute("history",  "testing");
      file.appendGlobalAttribute("history",  "testing2");
      EXPECT_EQ("testing\nempty\ntesting2", file.getGlobalAttribute("history"));

      // Writing should not thrown an error
      std::vector<Variable> vars;
      file.write(vars);
   }
   TEST_F(FileNetcdfTest, appendAttributeEmpty) {
      // Check that appending and prepending to an empty attribute works
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      file.prependGlobalAttribute("history71623",  "value321");
      file.appendGlobalAttribute("history99311",  "value15");
      EXPECT_EQ("value321", file.getGlobalAttribute("history71623"));
      EXPECT_EQ("value15",  file.getGlobalAttribute("history99311"));
   }
   TEST_F(FileNetcdfTest, setAttribute) {
      // Check that appending and prepending to an empty attribute works
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      file.setGlobalAttribute("att1",     "value93824");
      file.appendGlobalAttribute("att1",  "append");
      file.setGlobalAttribute("att1",     "value321192839819");

      file.setAttribute("air_temperature_2m", "att1", "value71");
      file.setAttribute("air_temperature_2m", "att1", "value72");
      file.setAttribute("air_temperature_2m", "att1", "value73");

      file.setGlobalAttribute("att2",  "value15");
      std::vector<Variable> vars;
      vars.push_back(mVariable);
      file.write(vars);
      FileNetcdf file2 = FileNetcdf("testing/files/10x10_copy.nc");
      EXPECT_EQ("value321192839819", file.getGlobalAttribute("att1"));
      EXPECT_EQ("value15",  file.getGlobalAttribute("att2"));
      EXPECT_EQ("value73",  file.getAttribute("air_temperature_2m", "att1"));
      EXPECT_EQ("",  file.getAttribute("air_temperature_2m", "att2"));
   }
   TEST_F(FileNetcdfTest, inandoutOfDataMode) {
      // Check that we can go in and out of data mode without error
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      // Define attributes
      file.setAttribute("air_temperature_2m", "att1", "value71");
      EXPECT_EQ("value71",  file.getAttribute("air_temperature_2m", "att1"));
      // Read data
      FieldPtr field = file.getField(mVariable, 0);
      std::vector<Variable> vars(1, mVariable);
      // Write data
      file.write(vars);

      // Add more attributes
      file.setAttribute("air_temperature_2m", "att1", "value72");
      EXPECT_EQ("value72",  file.getAttribute("air_temperature_2m", "att1"));
   }
   TEST_F(FileNetcdfTest, setAttributeError) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");

      // Variable do not exist
      EXPECT_DEATH(file.setAttribute("nonvalid_variable", "units", "value93824"), ".*");
      EXPECT_DEATH(file.getAttribute("q", "att1"), ".*");
   }
   TEST_F(FileNetcdfTest, setLongAttribute) {
      // Attempt to create a really long attribute
      {
         FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
         std::stringstream ss;
         for(int i = 0; i < 1e7; i++) {
            ss << "1234567890";
         }
         ss << "1234";
         std::string value = ss.str();
         file.appendGlobalAttribute("history", value);
         std::vector<Variable> vars(1, mVariable);
         file.write(vars);
      }
      // Make sure the attribute hasn't been set to the really long value
      FileNetcdf file = FileNetcdf("testing/files/10x10_copy.nc");
      std::string value = file.getGlobalAttribute("history");
      EXPECT_TRUE(value.size() < 1e8);
   }
   TEST_F(FileNetcdfTest, noTimeDimension) {
      FileNetcdf file("testing/files/validNetcdf3.nc");
      EXPECT_EQ(1, file.getNumTime());
   }
   TEST_F(FileNetcdfTest, createNewVariable) {
      FileNetcdf file("testing/files/10x10_copy.nc");
      std::vector<Variable> vars;
      Variable var("pop");
      vars.push_back(var);
      std::vector<float> pars(8,0);
      Parameters par(pars);
      ParameterFileSimple parFile(par);
      file.initNewVariable(var);
      CalibratorZaga cal(var, Options("popVariable=pop neighbourhoodSize=1 fracThreshold=0.4 popThreshold=0.5 6h=1"));
      cal.calibrate(file, &parFile);
      file.write(vars);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
