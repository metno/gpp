#include "../File/Arome.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileTest : public ::testing::Test {
   };

   TEST_F(FileTest, 10x10) {
      File* file = File::getScheme("testing/files/10x10.nc");
      EXPECT_EQ("arome", ((FileArome*)file)->name());
   }
   TEST_F(FileTest, 10x10_smart) {
      {
         FileArome from("testing/files/10x10.nc");
         FileArome to("testing/files/10x10_copy.nc");
         EXPECT_TRUE(from.hasVariable(Variable::T));
         DownscalerSmart d(Variable::T);
         std::vector<Variable::Type> variables;
         variables.push_back(Variable::T);
         d.downscale(from, to);

         to.write(variables);
      }
      // Nearest neighbour should give the same values
      FileArome f1("testing/files/10x10.nc");
      FileArome f2("testing/files/10x10_copy.nc");
      FieldPtr p1 = f1.getField(Variable::T, 0);
      FieldPtr p2 = f2.getField(Variable::T, 0);
      EXPECT_NE(*p1, *p2);
   }
   TEST_F(FileTest, hasVariable) {
      FileArome from("testing/files/10x10.nc");
      EXPECT_TRUE(from.hasVariable(Variable::PrecipAcc)); // Derivable
      EXPECT_TRUE(from.hasVariable(Variable::Precip));
      FieldPtr precip = from.getField(Variable::Precip, 0);
      FieldPtr precipAcc = from.getField(Variable::PrecipAcc, 0);
      EXPECT_FLOAT_EQ(0.911191, (*precip)[5][5][0]);
      EXPECT_FLOAT_EQ(0,        (*precipAcc)[5][5][0]);
   }
   TEST_F(FileTest, hasSameDimensions) {
      FileArome f1("testing/files/10x10.nc");
      FileArome f2("testing/files/10x10_copy.nc");
      FileFake f3(3,3,1,1);
      EXPECT_TRUE(f1.hasSameDimensions(f2));
      EXPECT_TRUE(f2.hasSameDimensions(f1));
      EXPECT_FALSE(f1.hasSameDimensions(f3));
      EXPECT_FALSE(f3.hasSameDimensions(f1));
   }
   TEST_F(FileTest, initNewVariable) {
      FileArome f1("testing/files/10x10.nc");
      EXPECT_FALSE(f1.hasVariable(Variable::U));
      f1.initNewVariable(Variable::U);
      EXPECT_TRUE(f1.hasVariable(Variable::U));
      FieldPtr field = f1.getField(Variable::U, 0);
   }
   TEST_F(FileTest, deriveVariables) {
      FileFake file(3, 3, 1, 3);
      ASSERT_TRUE(file.hasVariable(Variable::Precip));
      FieldPtr p0 = file.getField(Variable::Precip, 0);
      FieldPtr p1 = file.getField(Variable::Precip, 1);
      FieldPtr p2 = file.getField(Variable::Precip, 2);
      (*p0)[1][1][0] = 7.4;
      (*p1)[1][1][0] = 3.1;
      (*p2)[1][1][0] = 2.4;

      (*p0)[1][2][0] = 2.9;
      (*p1)[1][2][0] = Util::MV;
      (*p2)[1][2][0] = 6.1;

      (*p0)[0][2][0] = Util::MV;
      (*p1)[0][2][0] = Util::MV;
      (*p2)[0][2][0] = 6.1;

      (*p0)[0][0][0] = Util::MV;
      (*p1)[0][0][0] = 4.6;
      (*p2)[0][0][0] = 6.1;

      FieldPtr acc0 = file.getField(Variable::PrecipAcc, 0);
      FieldPtr acc1 = file.getField(Variable::PrecipAcc, 1);
      FieldPtr acc2 = file.getField(Variable::PrecipAcc, 2);
      EXPECT_FLOAT_EQ(0,        (*acc0)[1][1][0]);
      EXPECT_FLOAT_EQ(3.1,      (*acc1)[1][1][0]);
      EXPECT_FLOAT_EQ(5.5,      (*acc2)[1][1][0]);
      EXPECT_FLOAT_EQ(0,        (*acc0)[1][2][0]);
      EXPECT_FLOAT_EQ(Util::MV, (*acc1)[1][2][0]);
      EXPECT_FLOAT_EQ(Util::MV, (*acc2)[1][2][0]);
      EXPECT_FLOAT_EQ(0,        (*acc0)[0][2][0]);
      EXPECT_FLOAT_EQ(Util::MV, (*acc1)[0][2][0]);
      EXPECT_FLOAT_EQ(Util::MV, (*acc2)[0][2][0]);
      EXPECT_FLOAT_EQ(0,        (*acc0)[0][0][0]);
      EXPECT_FLOAT_EQ(4.6,      (*acc1)[0][0][0]);
      EXPECT_FLOAT_EQ(10.7,     (*acc2)[0][0][0]);
   }
   TEST_F(FileTest, factoryFake) {
      File* f = File::getScheme("fake");
      EXPECT_EQ("fake", f->name());
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
