#include "../File/Netcdf.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileTest : public ::testing::Test {
      protected:
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(FileTest, 10x10) {
      File* file = File::getScheme("testing/files/10x10.nc", Options());
      EXPECT_EQ("netcdf", ((FileNetcdf*)file)->name());
   }
   TEST_F(FileTest, 10x10_smart) {
      {
         FileNetcdf from("testing/files/10x10.nc");
         FileNetcdf to("testing/files/10x10_copy.nc");
         EXPECT_TRUE(from.hasVariable(mVariable));
         DownscalerSmart d(mVariable, mVariable, Options());
         std::vector<Variable> variables;
         variables.push_back(mVariable);
         d.downscale(from, to);

         to.write(variables);
      }
      FileNetcdf f1("testing/files/10x10.nc");
      FileNetcdf f2("testing/files/10x10_copy.nc");
      FieldPtr p1 = f1.getField(mVariable, 0);
      FieldPtr p2 = f2.getField(mVariable, 0);
      EXPECT_NE(*p1, *p2);
   }
   TEST_F(FileTest, hasVariable) {
      FileNetcdf from("testing/files/10x10.nc");
      EXPECT_FALSE(from.hasVariable(Variable("precipitation_amount_acc"))); // Can't derive
      EXPECT_TRUE(from.hasVariable(Variable("precipitation_amount")));
      FieldPtr precip = from.getField(Variable("precipitation_amount"), 0);
      EXPECT_FLOAT_EQ(0.911191, (*precip)(5,5,0));
   }
   TEST_F(FileTest, test) {
      FileNetcdf from("testing/files/10x10.nc");
      EXPECT_FALSE(from.hasVariable(Variable("test")));
      FieldPtr field = from.getField(Variable("test"), 0);
      EXPECT_TRUE(from.hasVariable(Variable("test")));
   }
   TEST_F(FileTest, hasSameDimensions) {
      FileNetcdf f1("testing/files/10x10.nc");
      FileNetcdf f2("testing/files/10x10_copy.nc");
      FileFake f3(Options("nLat=3 nLon=3 nEns=1 nTime=1"));
      EXPECT_TRUE(f1.hasSameDimensions(f2));
      EXPECT_TRUE(f2.hasSameDimensions(f1));
      EXPECT_FALSE(f1.hasSameDimensions(f3));
      EXPECT_FALSE(f3.hasSameDimensions(f1));
   }
   TEST_F(FileTest, initNewVariable) {
      FileNetcdf f1("testing/files/10x10.nc");
      Variable variable("fake");
      EXPECT_FALSE(f1.hasVariable(variable));
      f1.initNewVariable(variable);
      EXPECT_TRUE(f1.hasVariable(variable));
      FieldPtr field = f1.getField(variable, 0);
   }
   TEST_F(FileTest, deriveVariables) {
      FileFake file(Options("nLat=3 nLon=3 nEns=1 nTime=3"));
      Variable variable("precipitation");
      ASSERT_TRUE(file.hasVariable(variable));
      FieldPtr p0 = file.getField(variable, 0);
      FieldPtr p1 = file.getField(variable, 1);
      FieldPtr p2 = file.getField(variable, 2);
      (*p0)(1,1,0) = 7.4;
      (*p1)(1,1,0) = 3.1;
      (*p2)(1,1,0) = 2.4;

      (*p0)(1,2,0) = 2.9;
      (*p1)(1,2,0) = Util::MV;
      (*p2)(1,2,0) = 6.1;

      (*p0)(0,2,0) = Util::MV;
      (*p1)(0,2,0) = Util::MV;
      (*p2)(0,2,0) = 6.1;

      (*p0)(0,0,0) = Util::MV;
      (*p1)(0,0,0) = 4.6;
      (*p2)(0,0,0) = 6.1;

   }
   TEST_F(FileTest, getFieldInvalidTime) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      FileFake f0(Options("nLat=3 nLon=3 nEns=1 nTime=3"));
      EXPECT_DEATH(f0.getField(mVariable, 4), ".*");
      FileNetcdf f1("testing/files/10x10.nc");
      EXPECT_DEATH(f1.getField(mVariable, 100), ".*");
   }
   TEST_F(FileTest, getFieldInvalidTimeAfterValidAccess) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      FileFake f0(Options("nLat=3 nLon=3 nEns=1 nTime=3"));
      f0.getField(mVariable, 0);
      EXPECT_DEATH(f0.getField(mVariable, 4), ".*");
      FileNetcdf f1("testing/files/10x10.nc");
      EXPECT_DEATH(f1.getField(mVariable, 100), ".*");
   }
   TEST_F(FileTest, getFieldInvalidTimePreviouslyRead) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      FileFake f0(Options("nLat=3 nLon=3 nEns=1 nTime=3"));
      f0.getField(mVariable, 1);
      EXPECT_DEATH(f0.getField(mVariable, 4), ".*");
      EXPECT_DEATH(f0.getField(mVariable, 100), ".*");
   }
   TEST_F(FileTest, setgetTimes) {
      FileFake f0(Options("nLat=3 nLon=3 nEns=1 nTime=3"));
      std::vector<double> setTimes(3,0);
      setTimes[0] = 3.123;
      setTimes[1] = 4.624;
      setTimes[2] = 5;
      f0.setTimes(setTimes);

      std::vector<double> getTimes = f0.getTimes();
      EXPECT_EQ(3, getTimes.size());
      EXPECT_DOUBLE_EQ(3.123, getTimes[0]);
      EXPECT_DOUBLE_EQ(4.624, getTimes[1]);
      EXPECT_DOUBLE_EQ(5,     getTimes[2]);
   }
   TEST_F(FileTest, setgetReferenceTime) {
      FileFake f0(Options("nLat=3 nLon=3 nEns=1 nTime=3"));
      f0.setReferenceTime(4.1123);

      double referenceTime = f0.getReferenceTime();
      EXPECT_DOUBLE_EQ(4.1123, referenceTime);
   }
   TEST_F(FileTest, factoryMissing) {
      File* f = File::getScheme("missingfilename", Options());
      EXPECT_EQ(NULL, f);
   }
   TEST_F(FileTest, factoryNorcom) {
      File* f = File::getScheme("file=testing/files/norcom.txt", Options("type=norcomQnh lats=60 lons=9 elevs=100 names=Test numTimes=100 startTime=0 endTime=3"));
      ASSERT_TRUE(f);
      EXPECT_EQ("norcom", f->name());
   }
   /* TODO: Not implemented
   TEST_F(FileTest, deaccumulate) {
      // Create accumulation field
      FileNetcdf from("testing/files/1x1.nc");
      Variable var = Variable("precipitation_amount");

      // Accumulated    0, 3, 4, 4, 5.5,  10, _,12,12,20
      // Deaccumulated  _, 3, 1, 0, 1.5, 4.5, _, _, 0, 8

      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(var, 0))(0,0,0));
      EXPECT_FLOAT_EQ(3, (*from.getField(var, 1))(0,0,0));
      EXPECT_FLOAT_EQ(1, (*from.getField(var, 2))(0,0,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(var, 3))(0,0,0));
      EXPECT_FLOAT_EQ(1.5, (*from.getField(var, 4))(0,0,0));
      EXPECT_FLOAT_EQ(4.5, (*from.getField(var, 5))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(var, 6))(0,0,0));
      EXPECT_FLOAT_EQ(Util::MV, (*from.getField(var, 7))(0,0,0));
      EXPECT_FLOAT_EQ(0, (*from.getField(var, 8))(0,0,0));
      EXPECT_FLOAT_EQ(8, (*from.getField(var, 9))(0,0,0));
   }
   */
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
