#include "../File/Text.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileTextTest : public ::testing::Test {
      protected:
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(FileTextTest, asInput) {
      FileText file("testing/files/validText1.txt", Options());
      FieldPtr field0 = file.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(3.2, (*field0)(0,0,0));
      FieldPtr field1 = file.getField(mVariable, 1);
      EXPECT_FLOAT_EQ(4.1, (*field1)(0,0,0));
   }
   TEST_F(FileTextTest, asEnsemble) {
      FileText file("testing/files/validText2.txt", Options("lat=1 lon=2 elev=3"));
      FieldPtr field0 = file.getField(mVariable, 0);
      ASSERT_EQ(3, field0->getNumEns());
      ASSERT_EQ(2, field0->getNumY());
      ASSERT_EQ(1, field0->getNumX());
      EXPECT_FLOAT_EQ(11, (*field0)(0,0,0));
      EXPECT_FLOAT_EQ(21, (*field0)(0,0,1));
      EXPECT_FLOAT_EQ(-1, (*field0)(0,0,2));
      EXPECT_FLOAT_EQ(3.2, (*field0)(1,0,0));
      EXPECT_FLOAT_EQ(1.5, (*field0)(1,0,1));
      EXPECT_FLOAT_EQ(5.1, (*field0)(1,0,2));
      FieldPtr field1 = file.getField(mVariable, 1);
      ASSERT_EQ(3, field1->getNumEns());
      EXPECT_FLOAT_EQ(Util::MV, (*field1)(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*field1)(0,0,1));
      EXPECT_FLOAT_EQ(Util::MV, (*field1)(0,0,2));
      EXPECT_FLOAT_EQ(4, (*field1)(1,0,0));
      EXPECT_FLOAT_EQ(1, (*field1)(1,0,1));
      EXPECT_FLOAT_EQ(2, (*field1)(1,0,2));
   }
   TEST_F(FileTextTest, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(FileText("testing/files/invalidText1.txt", Options()), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
