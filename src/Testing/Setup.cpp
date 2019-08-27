#include <gtest/gtest.h>
#include "../Setup.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"
typedef Setup MetSetup;

namespace {
   class SetupTest : public ::testing::Test {
      protected:
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };
   TEST_F(SetupTest, test1) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -c zaga -p testing/files/parameters.txt type=text -c accumulate -d smart radius=11"));
      EXPECT_EQ(1, setup.variableConfigurations.size());
      EXPECT_EQ(2, setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ(mVariable, setup.variableConfigurations[0].inputVariable);
   }
   TEST_F(SetupTest, test2) {
      std::vector<std::string> lines;
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -vi air_temperature_2m -v out -vi air_temperature_2m -v out2");
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -vi air_temperature_2m -v out -d nearestNeighbour -vi air_temperature_2m -v out2");
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -vi air_temperature_2m -v out -d nearestNeighbour -vi air_temperature_2m -v out2 -d bilinear");
      for(int i = 0; i < lines.size(); i++) {
         MetSetup setup(Util::split(lines[i]));
         ASSERT_EQ(2, setup.variableConfigurations.size());
         EXPECT_EQ("air_temperature_2m", setup.variableConfigurations[0].inputVariable.name());
         EXPECT_EQ("air_temperature_2m", setup.variableConfigurations[1].inputVariable.name());
         EXPECT_EQ("out", setup.variableConfigurations[0].outputVariable.name());
         EXPECT_EQ("out2", setup.variableConfigurations[1].outputVariable.name());
      }
   }
   TEST_F(SetupTest, test3) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -d nearestNeighbour -vi air_temperature_2m -v out2"));
      ASSERT_EQ(2, setup.variableConfigurations.size());
      EXPECT_EQ("precipitation_amount", setup.variableConfigurations[0].inputVariable.name());
      EXPECT_EQ("air_temperature_2m", setup.variableConfigurations[1].inputVariable.name());
      EXPECT_EQ("precipitation_amount", setup.variableConfigurations[0].outputVariable.name());
      EXPECT_EQ("out2", setup.variableConfigurations[1].outputVariable.name());
   }
   TEST_F(SetupTest, calibratorOptions) {
      // Test that the calibrator picks up the right options
      std::vector<std::string> lines;
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -c neighbourhood radius=2 -p testing/files/parameters.txt type=text opt=1");
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c neighbourhood radius=11 -v air_temperature_2m -c neighbourhood radius=2 -p testing/files/parameters.txt type=text opt=1");
      for(int i = 0; i < lines.size(); i++) {
         // The last calibrator should have a radius of 2
         MetSetup setup(Util::split(lines[i]));
         int last = setup.variableConfigurations.size() - 1;
         EXPECT_EQ(1, setup.variableConfigurations[last].calibrators.size());
         EXPECT_EQ(1, setup.variableConfigurations[last].parameterFileCalibrators.size());
         Options options = setup.variableConfigurations[last].calibrators[0]->getOptions();
         int radius = Util::MV;
         options.getValue("radius", radius);
         EXPECT_EQ(2, radius);

         options = setup.variableConfigurations[last].parameterFileCalibrators[0]->getOptions();
         int opt = Util::MV;
         options.getValue("opt", opt);
         EXPECT_EQ(1, opt);
      }
   }
   TEST_F(SetupTest, calibratorOptionsMultiple) {
      // Check that radis=2 gets assigned to the correct calibrator
      std::vector<std::string> lines;
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -c accumulate -c neighbourhood radius=2 -c deaccumulate");
      for(int i = 0; i < lines.size(); i++) {
         MetSetup setup(Util::split(lines[i]));
         EXPECT_EQ(3, setup.variableConfigurations[0].calibrators.size());
         EXPECT_EQ(3, setup.variableConfigurations[0].parameterFileCalibrators.size());
         for(int c = 0; c < 3; c++) {
            Calibrator* cal = setup.variableConfigurations[0].calibrators[c];
            Options options = cal->getOptions();
            if(cal->name() == "neighbourhood") {
               int radius = Util::MV;
               options.getValue("radius", radius);
               EXPECT_EQ(2, radius);
            }
            else {
               int radius = Util::MV;
               options.getValue("radius", radius);
               EXPECT_EQ(Util::MV, radius);
            }
         }
      }
   }
   TEST_F(SetupTest, calibratorMultipleVariables) {
      // Check that cal options for the first calibrator does not get used for a calibrator
      // for a different variable
      std::vector<std::string> lines;
      lines.push_back("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -c accumulate -v surface_air_pressure -c diagnoseWind x=x_wind_10m y=y_wind_10m compute=speed");
      for(int i = 0; i < lines.size(); i++) {
         MetSetup setup(Util::split(lines[i]));
         EXPECT_EQ(1, setup.variableConfigurations[0].calibrators.size());
         EXPECT_EQ(1, setup.variableConfigurations[1].calibrators.size());
         Calibrator* cal = setup.variableConfigurations[0].calibrators[0];
         EXPECT_EQ(cal->name(), "accumulate");
         cal = setup.variableConfigurations[1].calibrators[0];
         EXPECT_EQ(cal->name(), "diagnoseWind");
      }
   }
   TEST_F(SetupTest, variableOnly) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(mVariable,                setup.variableConfigurations[0].inputVariable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST_F(SetupTest, allowBypass) {
      // Test that undefined variables can be passed down past the downscaler
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v wetbulb -d bypass"));
      ASSERT_EQ(1, setup.variableConfigurations.size());
      EXPECT_EQ(Variable("wetbulb"), setup.variableConfigurations[0].outputVariable.name());
      EXPECT_EQ("bypass", setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST_F(SetupTest, valid) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      EXPECT_EQ(mVariable,  setup.variableConfigurations[0].inputVariable);
      EXPECT_EQ("smart",      setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,            setup.variableConfigurations[0].calibrators.size());
   }
   TEST_F(SetupTest, downscaler_parameters) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -c zaga -p testing/files/parameters.txt type=text"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      EXPECT_EQ(mVariable,  setup.variableConfigurations[0].inputVariable);
      EXPECT_EQ(1,            setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ("zaga",       setup.variableConfigurations[0].calibrators[0]->name());
      EXPECT_EQ("text",       setup.variableConfigurations[0].parameterFileCalibrators[0]->name());
      EXPECT_EQ("testing/files/parameters.txt",       setup.variableConfigurations[0].parameterFileCalibrators[0]->getFilename());
   }
   TEST_F(SetupTest, repeatVariable) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -v air_temperature_2m -d smart -c neighbourhood"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(mVariable,                setup.variableConfigurations[0].inputVariable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST_F(SetupTest, repeatDownscaler) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart -d nearestNeighbour"));
      ASSERT_EQ(1,                  setup.variableConfigurations.size());
      EXPECT_EQ("nearestNeighbour", setup.variableConfigurations[0].downscaler->name());
   }
   TEST_F(SetupTest, complicated) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d nearestNeighbour -d smart -c neighbourhood -c accumulate -c neighbourhood -v precipitation_amount -c zaga -p testing/files/parameters.txt type=text -d gradient"));
      ASSERT_EQ(2,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];
      EXPECT_EQ(mVariable,  varconf.inputVariable);
      EXPECT_EQ("smart",      varconf.downscaler->name());
      ASSERT_EQ(3,            varconf.calibrators.size());
      EXPECT_EQ("neighbourhood",     varconf.calibrators[0]->name());
      EXPECT_EQ("accumulate", varconf.calibrators[1]->name());
      EXPECT_EQ("neighbourhood",     varconf.calibrators[2]->name());

      EXPECT_EQ(Variable("precipitation_amount"), setup.variableConfigurations[1].inputVariable);
      EXPECT_EQ("gradient",   setup.variableConfigurations[1].downscaler->name());
      ASSERT_EQ(1,            setup.variableConfigurations[1].calibrators.size());
      EXPECT_EQ("zaga",       setup.variableConfigurations[1].calibrators[0]->name());
   }
   TEST_F(SetupTest, differentInputOutputVariables) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -vi precipitation_amount -v air_temperature_2m"));
      ASSERT_EQ(1, setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      EXPECT_EQ("precipitation_amount", varconf.inputVariable.name());
      EXPECT_EQ("air_temperature_2m", varconf.outputVariable.name());
   }
   TEST_F(SetupTest, noOutput) {
      // Check that when no output file is specified, the input is used as output
      std::vector<std::string> lines;
      lines.push_back("testing/files/10x10.nc -v air_temperature_2m");
      lines.push_back("testing/files/10x10.nc -v air_temperature_2m -d bilinear");
      for(int i = 0; i < lines.size(); i++) {
         MetSetup setup(Util::split(lines[i]));
         ASSERT_EQ(1, setup.variableConfigurations.size());
         ASSERT_EQ("testing/files/10x10.nc", setup.inputFiles[0]->getFilename());
         ASSERT_EQ("testing/files/10x10.nc", setup.outputFiles[0]->getFilename());
      }
      // Flip the dimensions
      // std::string line("testing/files/10x10.nc xDim=latitude yDim=longitude -v air_temperature_2m -d bilinear");
      std::string line("testing/files/3x3.nc xDim=y yDim=x -v air_temperature_2m -d bilinear");
      MetSetup setup(Util::split(line));
      ASSERT_EQ(1, setup.variableConfigurations.size());
      ASSERT_EQ("testing/files/3x3.nc", setup.inputFiles[0]->getFilename());
      ASSERT_EQ("testing/files/3x3.nc", setup.outputFiles[0]->getFilename());
      vec2 ilats = setup.inputFiles[0]->getLats();
      vec2 olats = setup.outputFiles[0]->getLats();
      ASSERT_EQ(ilats, olats);
   }
   TEST_F(SetupTest, variableOptionsSingle) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m write=0"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      Options vOptions = varconf.outputVariableOptions;
      bool doWrite = true;
      ASSERT_TRUE(vOptions.getValue("write", doWrite));
      ASSERT_FALSE(vOptions.getValue("-d", doWrite));
      EXPECT_EQ(0, doWrite);
   }
   TEST_F(SetupTest, variableOptionsMultiple) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -v surface_air_pressure write=0 -v relative_humidity_2m -v x_wind_10m test=2 -d smart -v y_wind_10m -v precipitation_amount new=2.1 -c neighbourhood"));
      ASSERT_EQ(6, setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      EXPECT_EQ(mVariable, varconf.inputVariable);
      Options vOptions = varconf.outputVariableOptions;
      bool doWrite = true;
      float value = -1;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[1];
      EXPECT_EQ(Variable("surface_air_pressure"), varconf.inputVariable);
      vOptions = varconf.outputVariableOptions;
      ASSERT_TRUE(vOptions.getValue("write", doWrite));
      EXPECT_EQ(0, doWrite);
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[2];
      EXPECT_EQ(Variable("relative_humidity_2m"), varconf.inputVariable);
      vOptions = varconf.outputVariableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[3];
      EXPECT_EQ(Variable("x_wind_10m"), varconf.inputVariable);
      EXPECT_EQ("smart", varconf.downscaler->name());
      vOptions = varconf.outputVariableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      ASSERT_TRUE(vOptions.getValue("test", value));
      EXPECT_FLOAT_EQ(2, value);
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[4];
      EXPECT_EQ(Variable("y_wind_10m"), varconf.inputVariable);
      vOptions = varconf.outputVariableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[5];
      EXPECT_EQ(Variable("precipitation_amount"), varconf.inputVariable);
      vOptions = varconf.outputVariableOptions;
      ASSERT_TRUE(vOptions.getValue("new", value));
      EXPECT_FLOAT_EQ(2.1, value);
      EXPECT_FALSE(vOptions.getValue("-v", value));
      EXPECT_EQ(1, varconf.calibrators.size());
      EXPECT_EQ("neighbourhood", varconf.calibrators[0]->name());
   }
   TEST_F(SetupTest, shouldBeValid) {
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -c neighbourhood -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart -c neighbourhood"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart -c neighbourhood neighbourhood"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart -c neighbourhood neighbourhood -v precipitation_amount -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d nearestNeighbour -d smart -c neighbourhood accumulate neighbourhood -v precipitation_amount -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d nearestNeighbour -v precipitation_amount -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart numSmart=2 -c neighbourhood -v precipitation_amount -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart numSmart=2 -v precipitation_amount -d smart"));
      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -v air_temperature_2m -d smart numSmart=2 -v precipitation_amount -d smart"));

      MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -v precipitation_amount -c zaga -p testing/files/parameters.txt type=text -c zaga -p testing/files/parameters.txt type=text -d nearestNeighbour"));

      // Multiple inputs and outputs
      MetSetup(Util::split("testing/files/10x10.nc,testing/files/10x10_copy.nc testing/files/10x10.nc,testing/files/10x10_copy.nc -v precipitation_amount -c zaga -p testing/files/parameters.txt type=text -c zaga -p testing/files/parameters.txt type=text -d nearestNeighbour"));
   }
   TEST_F(SetupTest, shouldBeInValid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      // No variables
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v -d smart")), ".*");
      // Too many files (should have max two files)
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc testing/files/10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      // Invalid input file
      EXPECT_DEATH(MetSetup(Util::split("testing/files/weoihwoiedoiwe10x10.nc testing/files/10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      // Invalid output file
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/weoihwoiedoiwe10x10.nc -v -d smart testing/files/10x10.nc")), ".*");
      // Nothing after downscaler
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -d")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -d -c neighbourhood")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c neighbourhood -d")), ".*");
      // Nothing after calibrator
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -d nearest -c")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c -d nearest")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c -d nearest")), ".*");

      // Parameters before other schemes
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -p testing/files/parameters.txt type=text -v precipitation_amount -c zaga")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -p testing/files/parameters.txt type=text -c zaga")), ".*");

      // Invalid parameter file
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c zaga -p testing/files/parametersw8e9yhd89hywe89d.txt type=text")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c zaga -p testing/files/parametersw8e9yhd89hywe89d.txt type=text")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c zaga -p testing/files/parametersInvalidTime.txt type=text")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v precipitation_amount -c zaga -p testing/files/parametersw8e9yhd89hywe89d.txt type=netcdf")), ".*");

      // Different number of input and output files
      EXPECT_DEATH(MetSetup(Util::split("\"testing/files/10x10*.nc\" \"testing/files/10x10.nc\" -v precipitation_amount -c zaga -p testing/files/parametersw8e9yhd89hywe89d.txt type=netcdf")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("\"testing/files/10x10.nc\" \"testing/files/10x10*.nc\" -v precipitation_amount -c zaga -p testing/files/parametersw8e9yhd89hywe89d.txt type=netcdf")), ".*");

      // -vi but no -v
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -vi precipitation_amount")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -vi precipitation_amount -d nearestNeighbour")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -vi precipitation_amount -d nearestNeighbour -v out")), ".*");

      // Missing type
      EXPECT_DEATH(MetSetup(Util::split("\"testing/files/10x10*.nc\" \"testing/files/10x10.nc\" -v precipitation_amount -c zaga -p testing/files/parameters.txt")), ".*");

      // Missing dimension
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc xDim=test testing/files/10x10_copy.nc -v precipitation_amount -c zaga -p testing/files/parameters.txt")), ".*");

      // Different input output options on the same file
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc xDim=y testing/files/10x10.nc -v precipitation_amount -c zaga -p testing/files/parameters.txt")), ".*");

      // Variable does not exist in input
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v wetbulb")), ".*");
      EXPECT_DEATH(MetSetup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v wetbulb -d nearestNeighbour")), ".*");

   }
   TEST_F(SetupTest, defaultDownscaler) {
      std::string downscaler = Setup::defaultDownscaler();
      EXPECT_NE("", downscaler);
   }
   TEST_F(SetupTest, variableConfiguration) {
      VariableConfiguration varconf;
   }
   TEST_F(SetupTest, destructor) {
      {
         MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10.nc -v air_temperature_2m -d smart numSmart=2 -v precipitation_amount -d smart"));
      }
   }
   TEST_F(SetupTest, inputoutputOptions) {
      MetSetup setup0(Util::split("testing/files/10x10.nc option1=1 testing/files/10x10_copy.nc option2=2 -v air_temperature_2m write=1 -d smart numSmart=2"));
      int i;
      EXPECT_TRUE(setup0.inputOptions.getValue("option1", i));
      EXPECT_FALSE(setup0.inputOptions.getValue("option2", i));
      EXPECT_FALSE(setup0.inputOptions.getValue("write", i));
      EXPECT_EQ(1, i);
      EXPECT_TRUE(setup0.outputOptions.getValue("option2", i));
      EXPECT_FALSE(setup0.outputOptions.getValue("option1", i));
      EXPECT_FALSE(setup0.outputOptions.getValue("write", i));
      EXPECT_EQ(2, i);
   }
   TEST_F(SetupTest, alias) {
      MetSetup setup(Util::split("testing/files/10x10.nc testing/files/10x10_copy.nc -va tlevel1 name=air_temperature_2m level=1 -v air_temperature_2m -d smart numSmart=2"));
      ASSERT_EQ(1, setup.variableAliases.size());
      ASSERT_TRUE(setup.variableAliases.find("tlevel1") != setup.variableAliases.end());
      Variable var = setup.variableAliases["tlevel1"];
      EXPECT_EQ("air_temperature_2m", var.name());
      EXPECT_EQ(1, var.level());
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
