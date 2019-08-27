#include <gtest/gtest.h>
#include "../Variable.h"
#include "../Util.h"

namespace {

   TEST(VariableTest, test) {
      Variable variable("test", "units");
      EXPECT_EQ("test", variable.name());
      EXPECT_EQ("units", variable.units());

      Variable variable2 = variable;
      EXPECT_EQ("test", variable2.name());
      EXPECT_EQ("units", variable2.units());
      EXPECT_EQ(variable, variable2);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
