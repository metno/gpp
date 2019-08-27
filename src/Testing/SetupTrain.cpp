#include <gtest/gtest.h>
#include "../SetupTrain.h"
#include "../Downscaler/Smart.h"
#include "../Calibrator/Calibrator.h"

namespace {

   /*
   TEST(TestSetupTrain, test1) {
      std::vector<std::string> setups;
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -c gaussian -p testing/files/calculatedParameters.txt type=text");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -p testing/files/calculatedParameters.txt type=text -c gaussian");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p testing/files/calculatedParameters.txt type=text -v precipitation_amount -c gaussian");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p testing/files/calculatedParameters.txt type=text -c gaussian -v precipitation_amount");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -p testing/files/calculatedParameters.txt type=text -v precipitation_amount");
      setups.push_back("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -v precipitation_amount -p testing/files/calculatedParameters.txt type=text");
      for(int i = 0; i < setups.size(); i++) {
         SetupTrain setup(Util::split(setups[i]));
         ASSERT_TRUE(setup.observations.size() == 1);
         ASSERT_TRUE(setup.forecasts.size() == 1);
         EXPECT_EQ(setup.observations[0]->getFilename(), "testing/files/trainingDataObs.txt");
         EXPECT_EQ(setup.forecasts[0]->getFilename(), "testing/files/trainingDataFcst.txt");
         EXPECT_EQ("testing/files/calculatedParameters.txt", setup.output->getFilename());
         EXPECT_EQ("text", setup.output->name());
         EXPECT_EQ("gaussian", setup.method->name());
      }
   }
   TEST(TestSetupTrain, shouldBeValid) {
      // The order of -p -c -v should be irrelevant
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -c gaussian -p testing/files/calculatedParameters.txt type=text"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -p testing/files/calculatedParameters.txt type=text -c gaussian"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p testing/files/calculatedParameters.txt type=text -v precipitation_amount -c gaussian"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p testing/files/calculatedParameters.txt type=text -c gaussian -v precipitation_amount"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -p testing/files/calculatedParameters.txt type=text -v precipitation_amount"));
      SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -v precipitation_amount -p testing/files/calculatedParameters.txt type=text"));
   }
   TEST(TestSetupTrain, shouldBeInValid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      // Missing obs and or fcst data
      EXPECT_DEATH(SetupTrain(Util::split("")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("-v precipitation_amount -c zaga -p text file=testing/files/calculatedParameters.txt")), ".*");

      // No variables
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p output.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p output.txt type=text -v")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p output.txt type=text -v -c gaussian")), ".*");
      // Invalid input file
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/weoihwoiedoiwe10x10.nc testing/files/oiwejiojqewqoijoijdwq.nc -p output.txt type=text -v T -c gaussian")), ".*");

      // Nothing after -p
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p")), ".*");


      // Missing -c
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p testing/files/calculatedParameters.txt type=text -v T")), ".*");

      // Missing -p
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c zaga file=testing/files/calculatedParameters.txt -v T")), ".*");

      // Repeat arguments
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -c zaga -p testing/files/calculatedParameters.txt type=text -v T")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -v T -c zaga -p testing/files/calculatedParameters.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -c zaga -p testing/files/calculatedParameters.txt type=text -c gaussian")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -c zaga -c gaussian -p testing/files/calculatedParameters.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -v precipitation_amount -c zaga -p testing/files/calculatedParameters.txt type=text -p testing/files/calculatedParameters.txt type=text")), ".*");
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p testing/files/calculatedParameters.txt type=text -v precipitation_amount -c zaga -p testing/files/calculatedParameters.txt type=text")), ".*");

      // Junk after variable
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -p type=text -v precipitation_amount file=testing/files/calculatedParameters.txt -c zaga")), ".*");

      // Missing type in -o
      EXPECT_DEATH(SetupTrain(Util::split("testing/files/trainingDataObs.txt type=text testing/files/trainingDataFcst.txt type=text -c gaussian -v precipitation_amount -p testing/files/calculatedParameters.txt")), ".*");
   }
*/
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
