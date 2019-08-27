#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileTextTest, singleTime) {
      ParameterFileText file(Options("file=testing/files/parametersSingleTime.txt"));
      ASSERT_EQ(1, file.getTimes().size());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(9, par.size());
      EXPECT_FLOAT_EQ(-1.2021, par[0]);
      EXPECT_FLOAT_EQ(0.0007985, par[8]);

      EXPECT_EQ(file.getParameters(0).getValues(), file.getParameters(10).getValues());
      EXPECT_EQ(file.getParameters(3).getValues(), file.getParameters(25).getValues());
   }

   TEST(ParameterFileTextTest, noTime) {
      ParameterFileText file(Options("file=testing/files/parametersNoTime.txt"));
      ASSERT_EQ(1, file.getTimes().size());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(9, par.size());
      EXPECT_FLOAT_EQ(-1.2021, par[0]);
      EXPECT_FLOAT_EQ(0.0007985, par[8]);

      EXPECT_EQ(file.getParameters(0).getValues(), file.getParameters(10).getValues());
      EXPECT_EQ(file.getParameters(3).getValues(), file.getParameters(25).getValues());
   }

   TEST(ParameterFileTextTest, multipleTime) {
      ParameterFileText file(Options("file=testing/files/parametersMultipleTime.txt"));
      ASSERT_EQ(8, file.getTimes().size());
      Parameters par = file.getParameters(30);
      ASSERT_EQ(8, par.size());
      EXPECT_FLOAT_EQ(0.04198875, par[0]);
      EXPECT_FLOAT_EQ(-0.04039751, par[5]);
   }
   TEST(ParameterFileTextTest, empty) {
      ParameterFileText file(Options("file=testing/files/parametersEmpty.txt"));
      ASSERT_EQ(6, file.getTimes().size());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(0, par.size());
   }
   TEST(ParameterFileTextTest, invalidFiles) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(ParameterFileText(Options("file=testing/files/parametersUnevenRows.txt")), ".*");
      EXPECT_DEATH(ParameterFileText(Options("file=testing/files/parametersInvalidTime.txt")), ".*");
      EXPECT_DEATH(ParameterFileText(Options("file=testing/files/parametersInvalidEntries.txt")), ".*");
      EXPECT_FALSE(ParameterFileText::isValid("testing/files/parametersf98wey8y8y89rwe.txt"));
   }
   TEST(ParameterFileTextTest, invalidTime) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      ParameterFileText file(Options("file=testing/files/parametersMultipleTime.txt"));
      EXPECT_DEATH(file.getParameters(1100), ".*");
      EXPECT_DEATH(file.getParameters(-1), ".*");
   }

   // Spatial
   TEST(ParameterFileTextTest, default) {
      ParameterFileText file(Options("file=testing/files/parametersKriging.txt"));
      std::vector<Location> locations = file.getLocations();
      ASSERT_EQ(4, locations.size());
      std::set<Location> locationsSet(locations.begin(), locations.end());
      std::set<Location> expected;
      expected.insert(Location(9, 9, 800));
      expected.insert(Location(9, 8, 850));
      expected.insert(Location(5, 5, 140));
      expected.insert(Location(0, 0, 150));

      EXPECT_EQ(locationsSet, expected);
   }
   TEST(ParameterFileTextTest, write) {
      ParameterFileText file(Options("file=testing/files/parametersKriging.txt"));
      file.write("testing/files/parametersKriging2.txt");
      ParameterFileText file2(Options("file=testing/files/parametersKriging2.txt"));
      // Should have the same locations
      std::vector<Location> locations1 = file.getLocations();
      std::vector<Location> locations2 = file.getLocations();
      std::sort(locations1.begin(), locations1.end());
      std::sort(locations2.begin(), locations2.end());
      EXPECT_EQ(locations1, locations2);
   }
   TEST(ParameterFileTextTest, missingFirstHour) {
      ParameterFileText file(Options("file=testing/files/parametersKriging.txt"));
      Parameters par = file.getParameters(1, Location(0,0,150));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(4, par[0]);

      // No parameters at time 0, so should find the nearest neighbour
      par = file.getParameters(0, Location(0,0,150));
      ASSERT_EQ(1, par.size());
      EXPECT_FLOAT_EQ(4.2, par[0]);
   }
   TEST(ParameterFileTextTest, description) {
      ParameterFileText::description();
   }
   TEST(ParameterFileTextTest, writeError) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      ParameterFileText file(Options("file=testing/files/parametersKriging.txt"));
      // Shouldn't be able to write to a directory
      EXPECT_DEATH(file.write("testing/files/"), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
