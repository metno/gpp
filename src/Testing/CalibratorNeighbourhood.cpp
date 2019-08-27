#include "../File/Fake.h"
#include "../Util.h"
#include "../ParameterFile/ParameterFile.h"
#include "../Calibrator/Neighbourhood.h"
#include <gtest/gtest.h>

namespace {
   class TestCalibratorNeighbourhood : public ::testing::Test {
      protected:
         TestCalibratorNeighbourhood() {
         }
         virtual ~TestCalibratorNeighbourhood() {
         }
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };
   TEST_F(TestCalibratorNeighbourhood, 10x10_double) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 fast=1"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(304.6667, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306.1667, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(303,      (*after)(9,9,0));
      EXPECT_FLOAT_EQ(308.25,   (*after)(0,9,0));
      EXPECT_FLOAT_EQ(302,      (*after)(0,0,0));
      EXPECT_FLOAT_EQ(303,      (*after)(1,0,0));
      EXPECT_FLOAT_EQ(304.6667, (*after)(5,0,0));
      EXPECT_FLOAT_EQ(306.25,   (*after)(9,0,0));
      EXPECT_FLOAT_EQ(305.5,    (*after)(8,0,0));
      EXPECT_FLOAT_EQ(300+61.0/9,    (*after)(8,1,0));

      CalibratorNeighbourhood cal2 = CalibratorNeighbourhood(mVariable ,Options("radius=2 fast=1"));
      cal2.calibrate(from);
      EXPECT_FLOAT_EQ(304.73114, (*after)(5,2,0));
      // Exact equality is difficult to achieve to do accumulated round-off errors
      EXPECT_NEAR(305.355, (*after)(5,9,0), 0.001);
   }
   TEST_F(TestCalibratorNeighbourhood, 10x10_double_slow) {
      // The slow method should give the same results
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 fast=0"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);
      ASSERT_EQ(10, after->getNumY());
      ASSERT_EQ(10, after->getNumX());
      ASSERT_EQ(1,  after->getNumEns());

      EXPECT_FLOAT_EQ(304.6667, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306.1667, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(303,      (*after)(9,9,0));
      EXPECT_FLOAT_EQ(308.25,   (*after)(0,9,0));
      EXPECT_FLOAT_EQ(302,      (*after)(0,0,0));
      EXPECT_FLOAT_EQ(303,      (*after)(1,0,0));
      EXPECT_FLOAT_EQ(304.6667, (*after)(5,0,0));
      EXPECT_FLOAT_EQ(306.25,   (*after)(9,0,0));
      EXPECT_FLOAT_EQ(305.5,    (*after)(8,0,0));
      EXPECT_FLOAT_EQ(300+61.0/9,    (*after)(8,1,0));

      CalibratorNeighbourhood cal2 = CalibratorNeighbourhood(mVariable ,Options("radius=2 fast=0"));
      cal2.calibrate(from);
      EXPECT_FLOAT_EQ(304.73114, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(305.35556, (*after)(5,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, calibrateField) {
      FileNetcdf from("testing/files/10x10.nc");
      std::vector<std::string> stats;
      stats.push_back("mean");
      stats.push_back("median");
      stats.push_back("min");
      stats.push_back("max");
      for(int s = 0; s < stats.size(); s++) {
          std::stringstream ss;
          ss << "radius=1stat= " << stats[s] << " fast=0";
          CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable, Options(ss.str()));
          FieldPtr before = from.getField(mVariable, 0);
          FieldPtr after = from.getEmptyField(0);
          cal.calibrateField(*before, *after);
          cal.calibrateField(*before, *before);
          for(int x = 0; x < before->getNumX(); x++) {
              for(int y = 0; y < before->getNumY(); y++) {
                  EXPECT_EQ((*before)(y, x), (*after)(y, x));
              }
          }
      }
   }
   TEST_F(TestCalibratorNeighbourhood, 10x10_missingValues) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable, Options("radius=1"));
      FieldPtr field = from.getField(mVariable, 0);
      (*field)(4,1,0) = Util::MV;
      (*field)(4,2,0) = Util::MV;
      (*field)(4,3,0) = Util::MV;
      (*field)(5,1,0) = Util::MV;
      (*field)(5,2,0) = Util::MV;
      (*field)(5,3,0) = Util::MV;
      (*field)(6,1,0) = Util::MV;
      (*field)(6,2,0) = Util::MV;
      (*field)(6,3,0) = Util::MV;

      cal.calibrate(from);
      EXPECT_FLOAT_EQ(Util::MV, (*field)(5,2,0));
      EXPECT_FLOAT_EQ(304.6667, (*field)(5,1,0));
   }
   TEST_F(TestCalibratorNeighbourhood, mean) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=mean"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(304.6667, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306.1667, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(308.25,   (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, sum) {
      // The slow method should give the same results
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=sum fast=1"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(304.6667*9, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306.1667*6, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(303*4,      (*after)(9,9,0));
      EXPECT_FLOAT_EQ(308.25*4,   (*after)(0,9,0));
      EXPECT_FLOAT_EQ(302*4,      (*after)(0,0,0));
      EXPECT_FLOAT_EQ(303*6,      (*after)(1,0,0));
      EXPECT_FLOAT_EQ(304.6667*6, (*after)(5,0,0));
      EXPECT_FLOAT_EQ(306.25*4,   (*after)(9,0,0));
      EXPECT_FLOAT_EQ(305.5*6,    (*after)(8,0,0));
      EXPECT_FLOAT_EQ(300*9+61.0, (*after)(8,1,0));
   }
   TEST_F(TestCalibratorNeighbourhood, min) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=min"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(301, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(302, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(302, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, max) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=max"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(310, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(310, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(320, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, std) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=std"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(2.7080128, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(3.2871804, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(7.0133801, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, median) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=median"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(306, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(306, (*after)(5,9,0));
      EXPECT_FLOAT_EQ(305.5, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, quantile) {
      FileNetcdf from("testing/files/10x10.nc");
      CalibratorNeighbourhood cal = CalibratorNeighbourhood(mVariable ,Options("radius=1 stat=quantile quantile=0.9"));

      cal.calibrate(from);
      FieldPtr after = from.getField(mVariable, 0);

      EXPECT_FLOAT_EQ(306.8, (*after)(5,2,0));
      EXPECT_FLOAT_EQ(310,   (*after)(5,9,0));
      EXPECT_FLOAT_EQ(316.1, (*after)(0,9,0));
   }
   TEST_F(TestCalibratorNeighbourhood, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);

      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("stat=invalidStatType")), ".*");

      // Negative radius
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("radius=-1")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("radius=-1 stat=quantile")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("radius=-1 stat=quantile quantile=0.5")), ".*");

      // Missing quantile
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("stat=quantile")), ".*");
      // Invalid quantile
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("stat=quantile quantile=-0.1")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("stat=quantile quantile=1.1")), ".*");
      EXPECT_DEATH(CalibratorNeighbourhood(mVariable, Options("stat=quantile quantile=-999")), ".*");
   }
   TEST_F(TestCalibratorNeighbourhood, description) {
      CalibratorNeighbourhood::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
