#include <gtest/gtest.h>
#include "../Setup.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"
typedef Setup MetSetup;

namespace {

   TEST(SetupTest, test1) {
      MetSetup setup(Util::split("input output -v T -c zaga parameters=testing/files/parameters.txt -c accumulate -d smart searchRadius=11"));
      EXPECT_EQ(1,           setup.variableConfigurations.size());
      EXPECT_EQ(2,           setup.variableConfigurations[0].calibrators.size());
      EXPECT_EQ(Variable::T, setup.variableConfigurations[0].variable);
      EXPECT_EQ(11, ((DownscalerSmart*) setup.variableConfigurations[0].downscaler)->getSearchRadius());
   }
   TEST(SetupTest, variableOnly) {
      MetSetup setup(Util::split("input output -v T"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,                setup.variableConfigurations[0].variable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, valid) {
      MetSetup setup(Util::split("input output -v T -d smart"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,  setup.variableConfigurations[0].variable);
      EXPECT_EQ("smart",      setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,            setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, repeatVariable) {
      MetSetup setup(Util::split("input output -v T -v T -d smart -c smooth"));
      ASSERT_EQ(1,                          setup.variableConfigurations.size());
      EXPECT_EQ(Variable::T,                setup.variableConfigurations[0].variable);
      EXPECT_EQ(Setup::defaultDownscaler(), setup.variableConfigurations[0].downscaler->name());
      EXPECT_EQ(0,                          setup.variableConfigurations[0].calibrators.size());
   }
   TEST(SetupTest, repeatDownscaler) {
      MetSetup setup(Util::split("input output -v T -d smart -d nearestNeighbour"));
      ASSERT_EQ(1,                  setup.variableConfigurations.size());
      EXPECT_EQ("nearestNeighbour", setup.variableConfigurations[0].downscaler->name());
   }
   TEST(SetupTest, complicated) {
      MetSetup setup(Util::split("input output -v T -d nearestNeighbour -d smart -c smooth -c accumulate -c smooth -v Precip -c zaga parameters=testing/files/parameters.txt -d gradient"));
      ASSERT_EQ(2,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];
      EXPECT_EQ(Variable::T,  varconf.variable);
      EXPECT_EQ("smart",      varconf.downscaler->name());
      ASSERT_EQ(3,            varconf.calibrators.size());
      EXPECT_EQ("smooth",     varconf.calibrators[0]->name());
      EXPECT_EQ("accumulate", varconf.calibrators[1]->name());
      EXPECT_EQ("smooth",     varconf.calibrators[2]->name());

      EXPECT_EQ(Variable::Precip, setup.variableConfigurations[1].variable);
      EXPECT_EQ("gradient",   setup.variableConfigurations[1].downscaler->name());
      ASSERT_EQ(1,            setup.variableConfigurations[1].calibrators.size());
      EXPECT_EQ("zaga",       setup.variableConfigurations[1].calibrators[0]->name());
   }
   TEST(SetupTest, variableOptionsSingle) {
      MetSetup setup(Util::split("input output -v T write=0"));
      ASSERT_EQ(1,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      Options vOptions = varconf.variableOptions;
      bool doWrite = true;
      ASSERT_TRUE(vOptions.getValue("write", doWrite));
      ASSERT_FALSE(vOptions.getValue("-d", doWrite));
      EXPECT_EQ(0, doWrite);
   }
   TEST(SetupTest, variableOptionsMultiple) {
      MetSetup setup(Util::split("input output -v T -v P write=0 -v RH -v U test=2 -d smart -v V -v Precip new=2.1"));
      ASSERT_EQ(6,            setup.variableConfigurations.size());
      VariableConfiguration varconf = setup.variableConfigurations[0];

      EXPECT_EQ(Variable::T, varconf.variable);
      Options vOptions = varconf.variableOptions;
      bool doWrite = true;
      float value = -1;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[1];
      EXPECT_EQ(Variable::P, varconf.variable);
      vOptions = varconf.variableOptions;
      ASSERT_TRUE(vOptions.getValue("write", doWrite));
      EXPECT_EQ(0, doWrite);
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[2];
      EXPECT_EQ(Variable::RH, varconf.variable);
      vOptions = varconf.variableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[3];
      EXPECT_EQ(Variable::U, varconf.variable);
      EXPECT_EQ("smart", varconf.downscaler->name());
      vOptions = varconf.variableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      ASSERT_TRUE(vOptions.getValue("test", value));
      EXPECT_FLOAT_EQ(2, value);
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[4];
      EXPECT_EQ(Variable::V, varconf.variable);
      vOptions = varconf.variableOptions;
      EXPECT_FALSE(vOptions.getValue("write", doWrite));
      EXPECT_FALSE(vOptions.getValue("-v", value));

      varconf = setup.variableConfigurations[5];
      EXPECT_EQ(Variable::Precip, varconf.variable);
      vOptions = varconf.variableOptions;
      ASSERT_TRUE(vOptions.getValue("new", value));
      EXPECT_FLOAT_EQ(2.1, value);
      EXPECT_FALSE(vOptions.getValue("-v", value));
   }
   TEST(SetupTest, shouldBeValid) {
      MetSetup setup0(Util::split("input output"));
      MetSetup setup1(Util::split("input output -v T -d smart"));
      MetSetup setup2(Util::split("input output -v T -c smooth -d smart"));
      MetSetup setup3(Util::split("input output -v T -d smart -c smooth"));
      MetSetup setup4(Util::split("input output -v T -d smart -c smooth smooth"));
      MetSetup setup5(Util::split("input output -v T -d smart -c smooth smooth -v Precip -d smart"));
      MetSetup setup6(Util::split("input output -v T -d nearestNeighbour -d smart -c smooth accumulate smooth -v Precip -d smart"));
      MetSetup setup7(Util::split("input output -v T -d nearestNeighbour -v Precip -d smart"));
      MetSetup setup8(Util::split("input output -v T -d smart numSmart=2 -c smooth -v Precip -d smart"));
      MetSetup setup9(Util::split("input output -v T -d smart numSmart=2 -v Precip -d smart"));
   }
   TEST(SetupTest, shouldBeInValid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(MetSetup setup2(Util::split("input output -v")), ".*");
      EXPECT_DEATH(MetSetup setup2(Util::split("input output -v -d smart")), ".*");
   }
   TEST(SetupTest, defaultDownscaler) {
      std::string downscaler = Setup::defaultDownscaler();
      EXPECT_NE("", downscaler);
   }
   TEST(SetupTest, variableConfiguration) {
      VariableConfiguration varconf;
   }
   TEST(SetupTest, destructor) {
      {
         MetSetup setup(Util::split("input output -v T -d smart numSmart=2 -v Precip -d smart"));
      }
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
