#include "../ParameterFile/ParameterFile.h"
#include "../Util.h"
#include <gtest/gtest.h>
#include <algorithm>

namespace {
   TEST(ParameterFileSimpleTest, valid) {
      std::vector<float> values;
      values.push_back(2.1);
      values.push_back(3.4);
      values.push_back(-2.9);
      ParameterFileSimple file(values);
      ASSERT_EQ(1, file.getTimes().size());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(3, par.size());
      EXPECT_FLOAT_EQ(2.1, par[0]);
      EXPECT_FLOAT_EQ(3.4, par[1]);
      EXPECT_FLOAT_EQ(-2.9, par[2]);
   }
   TEST(ParameterFileSimpleTest, empty) {
      std::vector<float> values;
      ParameterFileSimple file(values);
      ASSERT_EQ(1, file.getTimes().size());
      Parameters par = file.getParameters(0);
      ASSERT_EQ(0, par.size());
      EXPECT_FALSE(file.isLocationDependent());
      std::vector<Location> locations;
      locations = file.getLocations();
      EXPECT_EQ(0, locations.size());
   }

}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
