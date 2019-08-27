#include "../Field.h"
#include "../Util.h"
#include <gtest/gtest.h>

namespace {
   class FieldTest : public ::testing::Test {
   };

   TEST_F(FieldTest, access) {
      Field field(1, 1, 3, 2);
      ASSERT_EQ(1, field.getNumY());
      ASSERT_EQ(1, field.getNumX());
      ASSERT_EQ(3, field.getNumEns());
      EXPECT_FLOAT_EQ(2, field(0,0,0));
      EXPECT_FLOAT_EQ(2, field(0,0,1));
      EXPECT_FLOAT_EQ(2, field(0,0,2));
   }
   TEST_F(FieldTest, ensembleAccess) {
      float def = 2;
      Field field(3, 2, 3, def);
      field(2,1,0) = 4.1;
      field(2,1,1) = Util::MV;
      field(2,1,2) = 3;
      std::vector<float> vec = field(2,1);
      EXPECT_FLOAT_EQ(4.1, vec[0]);
      EXPECT_FLOAT_EQ(Util::MV, vec[1]);
      EXPECT_FLOAT_EQ(3, vec[2]);
      std::vector<float> vec2 = field(0,1);
      EXPECT_FLOAT_EQ(def, vec2[0]);
      EXPECT_FLOAT_EQ(def, vec2[1]);
      EXPECT_FLOAT_EQ(def, vec2[2]);
   }
   TEST_F(FieldTest, equality) {
      Field field1(3, 2, 3, 3.5);
      Field field2(3, 2, 3, 3.5);
      EXPECT_TRUE(field1  == field2);
      EXPECT_FALSE(field1 != field2);
      field2(0,1,1) = Util::MV;
      EXPECT_FALSE(field1  == field2);
      EXPECT_TRUE(field1 != field2);
      field2(0,1,1) = 3.5;
      EXPECT_TRUE(field1  == field2);
      EXPECT_FALSE(field1 != field2);
      field2(0,1,1) = 2.5;
      EXPECT_FALSE(field1  == field2);
      EXPECT_TRUE(field1 != field2);
   }
   TEST_F(FieldTest, inEquality) {
      Field field1(3, 2, 3, 3.5);
      Field field2(3, 2, 3, 2.2);
      Field field3(2, 1, 3, 4.1);
      EXPECT_FALSE(field1 == field2);
      EXPECT_FALSE(field1 == field3);
      EXPECT_FALSE(field2 == field3);
      EXPECT_TRUE(field1 != field2);
      EXPECT_TRUE(field1 != field3);
      EXPECT_TRUE(field2 != field3);
   }
   TEST_F(FieldTest, assignment) {
      Field field(1, 1, 3, 2);
      field(0,0,0) = 4;
      field(0,0,1) = Util::MV;
      field(0,0,2) = 19.3;
      EXPECT_FLOAT_EQ(4, field(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, field(0,0,1));
      EXPECT_FLOAT_EQ(19.3, field(0,0,2));
      field(0,0,0) = Util::MV;
      field(0,0,1) = 3.1;
      field(0,0,2) = Util::MV;
      EXPECT_FLOAT_EQ(Util::MV, field(0,0,0));
      EXPECT_FLOAT_EQ(3.1, field(0,0,1));
      EXPECT_FLOAT_EQ(Util::MV, field(0,0,2));
   }
   TEST_F(FieldTest, invalidConstruction) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Field(-1, 1, 1, 1), ".*");
      EXPECT_DEATH(Field(1, -1, 1, 1), ".*");
      EXPECT_DEATH(Field(1, 1, -1, 1), ".*");
      EXPECT_DEATH(Field(Util::MV, 1, 1, 1), ".*");
      EXPECT_DEATH(Field(1, Util::MV, 1, 1), ".*");
      EXPECT_DEATH(Field(1, 1, Util::MV, 1), ".*");
      EXPECT_DEATH(Field(-1, Util::MV, 1, 1), ".*");
      EXPECT_DEATH(Field(1, -1, Util::MV, 1), ".*");
      EXPECT_DEATH(Field(Util::MV, 1, -1, 1), ".*");
   }
   TEST_F(FieldTest, invalidAccess) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      Field field(3, 3, 2, 4.0);
      EXPECT_DEATH(field(3,0,0), ".*");
      EXPECT_DEATH(field(0,3,0), ".*");
      EXPECT_DEATH(field(0,0,2), ".*");

      EXPECT_DEATH(field(12,0,0), ".*");
      EXPECT_DEATH(field(0,12,0), ".*");
      EXPECT_DEATH(field(0,0,12), ".*");

      EXPECT_DEATH(field(-1,0,0), ".*");
      EXPECT_DEATH(field(0,-1,0), ".*");
      EXPECT_DEATH(field(0,0,-1), ".*");

      EXPECT_DEATH(field(Util::MV,0,0), ".*");
      EXPECT_DEATH(field(0,Util::MV,0), ".*");
      EXPECT_DEATH(field(0,0,Util::MV), ".*");
   }
   TEST_F(FieldTest, invalidAssignment) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      Field field(3, 3, 2, 4.0);
      EXPECT_DEATH(field(3,0,0) = 1, ".*");
      EXPECT_DEATH(field(0,3,0) = 1, ".*");
      EXPECT_DEATH(field(0,0,2) = 1, ".*");

      EXPECT_DEATH(field(12,0,0) = 1, ".*");
      EXPECT_DEATH(field(0,12,0) = 1, ".*");
      EXPECT_DEATH(field(0,0,12) = 1, ".*");

      EXPECT_DEATH(field(-1,0,0) = 1, ".*");
      EXPECT_DEATH(field(0,-1,0) = 1, ".*");
      EXPECT_DEATH(field(0,0,-1) = 1, ".*");

      EXPECT_DEATH(field(Util::MV,0,0) = 1, ".*");
      EXPECT_DEATH(field(0,Util::MV,0) = 1, ".*");
      EXPECT_DEATH(field(0,0,Util::MV) = 1, ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
