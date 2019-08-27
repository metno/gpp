#include "../Options.h"
#include <gtest/gtest.h>

namespace {
   class OptionsTest : public ::testing::Test {
      protected:
   };

   TEST_F(OptionsTest, accessEmptyOptions) {
      Options options;
      float value = 441;
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_FLOAT_EQ(441, value);

      status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_FLOAT_EQ(441, value);
   }
   TEST_F(OptionsTest, addOptionInt) {
      Options options;
      options.addOption("key",1);
      int value = -1;
      bool status = options.getValue("key", value);
      EXPECT_TRUE(status);
      EXPECT_EQ(1, value);

      value = -1;
      status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_EQ(-1, value);
   }
   TEST_F(OptionsTest, emptyKey) {
      // Empty keys are not allowed
      Options options;
      std::string value = "new";
      bool status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);

      options.addOption("", "test");
      status = options.getValue("", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);
   }
   TEST_F(OptionsTest, emptyValue) {
      // Empty values are not allowed
      Options options;
      options.addOption("test", "");
      std::string value = "new";
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);
   }
   TEST_F(OptionsTest, emptyValueEqualSign) {
      Options options;
      options.addOptions("test=");
      std::string value = "new";
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("new", value);

      int valueInt = 2;
      status = options.getValue("test", valueInt);
      EXPECT_FALSE(status);
      EXPECT_EQ(2, valueInt);
   }
   TEST_F(OptionsTest, addOption) {
      Options options;
      options.addOption("test", 441);

      float value = 2;
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_FLOAT_EQ(441, value);

      value = 2;
      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_FLOAT_EQ(2, value);
   }
   TEST_F(OptionsTest, stringOption) {
      Options options;
      options.addOption("test", "before");

      std::string value = "";
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ("before", value);

      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("before", value);
   }
   TEST_F(OptionsTest, addOptionEqualSign) {
      Options options;
      options.addOptions("test=before");

      std::string value = "";
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ("before", value);

      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("before", value);
   }
   TEST_F(OptionsTest, clear) {
      Options options;
      options.addOptions("test=before");
      options.addOptions("other=some");

      options.clear();
      std::string value = "something";
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("something", value);

      options.addOptions("new=value");
      status = options.getValue("new", value);
      EXPECT_TRUE(status);
      EXPECT_EQ("value", value);

      status = options.getValue("missing", value);
      EXPECT_FALSE(status);
      EXPECT_EQ("value", value);
   }
   TEST_F(OptionsTest, overwrite) {
      Options options;
      options.addOptions("test=before");
      options.addOptions("test=2");

      float value = 3;
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ(2, value);
   }
   TEST_F(OptionsTest, constructor) {
      Options options("test=3 test=2");

      float value = 19;
      bool status = options.getValue("test", value);
      EXPECT_TRUE(status);
      EXPECT_EQ(2, value);
      status = options.getValue("new", value);
      EXPECT_FALSE(status);
      EXPECT_EQ(2, value);
   }
   TEST_F(OptionsTest, constructorEmpty) {
      Options options("");

      float value = 19;
      bool status = options.getValue("test", value);
      EXPECT_FALSE(status);
      EXPECT_EQ(19, value);
   }
   TEST_F(OptionsTest, toString) {
      Options options("test=0   q=4");
      std::string s = options.toString();
      EXPECT_LT(0, s.size());

      Options copyOptions(s);
      int i = -1;
      copyOptions.getValue("test", i);
      EXPECT_EQ(0, i);
      copyOptions.getValue("q", i);
      EXPECT_EQ(4, i);

      options.addOptions("f=3");
      s = options.toString();
      copyOptions = Options(s);
      copyOptions.getValue("f", i);
      EXPECT_EQ(3, i);
   }
   TEST_F(OptionsTest, getValues) {
      Options options("test=1,2,3");
      std::vector<float> values;
      options.getValues("test", values);
      ASSERT_EQ(3, values.size());
      EXPECT_FLOAT_EQ(1, values[0]);
      EXPECT_FLOAT_EQ(2, values[1]);
      EXPECT_FLOAT_EQ(3, values[2]);

      options = Options("test=3,1,2 new=4");
      options.getValues("test", values);
      ASSERT_EQ(3, values.size());
      EXPECT_FLOAT_EQ(3, values[0]);
      EXPECT_FLOAT_EQ(1, values[1]);
      EXPECT_FLOAT_EQ(2, values[2]);
      options.getValues("new", values);
      ASSERT_EQ(1, values.size());
      EXPECT_FLOAT_EQ(4, values[0]);
   }
   TEST_F(OptionsTest, getValuesString) {
      Options options("test=1,2,3");
      std::vector<std::string> values;
      options.getValues("test", values);
      ASSERT_EQ(3, values.size());

      options = Options("test=3,1,2 new=4");
      options.getValues("test", values);
      ASSERT_EQ(3, values.size());
      EXPECT_EQ("3", values[0]);
      EXPECT_EQ("1", values[1]);
      EXPECT_EQ("2", values[2]);
      options.getValues("new", values);
      ASSERT_EQ(1, values.size());
      EXPECT_EQ("4", values[0]);
   }
   TEST_F(OptionsTest, requiredValue) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      Options options = Options("test=3,1,2 new=4");
      std::vector<int> values;
      int value = Util::MV;
      options.getRequiredValues("test", values);
      ASSERT_EQ(3, values.size());
      EXPECT_EQ(3, values[0]);
      EXPECT_EQ(1, values[1]);
      EXPECT_EQ(2, values[2]);
      options.getValue("new", value);
      EXPECT_EQ(4, value);

      EXPECT_DEATH(options.getRequiredValues("missing", values), ".*");
      EXPECT_DEATH(options.getRequiredValue("missing", value), ".*");
   }
   TEST_F(OptionsTest, check) {
      std::vector<int> values;
      int value = Util::MV;
      Options options = Options("test=3,1,2 new=4 att1=1,2,3");
      options.getValues("test", values);
      options.getValue("new", value);
      EXPECT_EQ(4, value);
      EXPECT_FALSE(options.check());
      options.getValue("new", value);
      EXPECT_FALSE(options.check());
      options.getValues("att1", values);
      EXPECT_TRUE(options.check());
   }
   TEST_F(OptionsTest, checkRequired) {
      std::vector<int> values;
      int value = Util::MV;
      Options options = Options("test=3,1,2 new=4 att1=1,2,3");
      options.getRequiredValues("test", values);
      options.getRequiredValue("new", value);
      EXPECT_EQ(4, value);
      EXPECT_FALSE(options.check());
      options.getRequiredValue("new", value);
      EXPECT_FALSE(options.check());
      options.getRequiredValues("att1", values);
      EXPECT_TRUE(options.check());
   }
   TEST_F(OptionsTest, equality) {
      Options options1("test=1 other=2");
      Options options2("other=2 test=1");
      EXPECT_TRUE(options1 == options2);
      EXPECT_TRUE(options2 == options1);
      EXPECT_FALSE(options1 != options2);
      EXPECT_FALSE(options2 != options1);

      Options options3("");
      Options options4("");
      EXPECT_TRUE(options3 == options4);
      EXPECT_TRUE(options4 == options3);
      EXPECT_FALSE(options3 != options4);
      EXPECT_FALSE(options4 != options3);
   }
   TEST_F(OptionsTest, inequality) {
      Options options1("test=1 other=2");
      Options options2("other=2");
      EXPECT_FALSE(options1 == options2);
      EXPECT_FALSE(options2 == options1);
      EXPECT_TRUE(options1 != options2);
      EXPECT_TRUE(options2 != options1);

      Options options3("test=1 other=2");
      Options options4("test=1");
      EXPECT_FALSE(options3 == options4);
      EXPECT_FALSE(options4 == options3);
      EXPECT_TRUE(options3 != options4);
      EXPECT_TRUE(options4 != options3);

      Options options5("test=1 other=2");
      Options options6("test=1 new=2");
      EXPECT_FALSE(options5 == options6);
      EXPECT_FALSE(options6 == options5);
      EXPECT_TRUE(options5 != options6);
      EXPECT_TRUE(options6 != options5);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
