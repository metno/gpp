#include "../Util.h"
#include "../File/Netcdf.h"
#include "../Variable.h"
#include "../Calibrator/Calibrator.h"
#include <algorithm>
#include <math.h>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <unistd.h>

namespace {
   class UtilTest : public ::testing::Test {
      protected:
         std::vector<float> getHood(float v1, float v2, float v3) {
            std::vector<float> hood;
            hood.push_back(v1);
            hood.push_back(v2);
            hood.push_back(v3);
            return hood;
         }
         std::vector<float> getHood(float v1, float v2, float v3, float v4) {
            std::vector<float> hood;
            hood.push_back(v1);
            hood.push_back(v2);
            hood.push_back(v3);
            hood.push_back(v4);
            return hood;
         }
         // a00 a01
         // a10 a11
         vec2 get2x2Matrix(float a00, float a01, float a10, float a11) {
            vec2 matrix;
            matrix.resize(2);
            matrix[0].resize(2);
            matrix[1].resize(2);
            matrix[0][0] = a00;
            matrix[1][0] = a10;
            matrix[0][1] = a01;
            matrix[1][1] = a11;
            return matrix;
         }
         // a00 a01 a02
         // a10 a11 a12
         // a20 a21 a22
         vec2 get3x3Matrix(float a00, float a01, float a02, float a10, float a11, float a12, float a20, float a21, float a22) {
            vec2 matrix;
            matrix.resize(3);
            matrix[0].resize(3);
            matrix[1].resize(3);
            matrix[2].resize(3);
            matrix[0][0] = a00;
            matrix[1][0] = a10;
            matrix[2][0] = a20;
            matrix[0][1] = a01;
            matrix[1][1] = a11;
            matrix[2][1] = a21;
            matrix[0][2] = a02;
            matrix[1][2] = a12;
            matrix[2][2] = a22;
            return matrix;
         }
   };

   TEST_F(UtilTest, isValid) {
      EXPECT_FALSE(Util::isValid(Util::MV));
      EXPECT_FALSE(Util::isValid(1.0/0));
      EXPECT_TRUE(Util::isValid(-15));
      EXPECT_TRUE(Util::isValid(0));
      EXPECT_TRUE(Util::isValid(3.14));
   }
   TEST_F(UtilTest, deg2rad) {
      EXPECT_FLOAT_EQ(2*Util::pi, Util::deg2rad(360));
      EXPECT_FLOAT_EQ(-2*Util::pi, Util::deg2rad(-360));
      EXPECT_FLOAT_EQ(0, Util::deg2rad(0));
      EXPECT_FLOAT_EQ(1.0/3*Util::pi, Util::deg2rad(60));
   }
   TEST_F(UtilTest, rad2deg) {
      EXPECT_FLOAT_EQ(360, Util::rad2deg(2*Util::pi));
      EXPECT_FLOAT_EQ(-360, Util::rad2deg(-2*Util::pi));
      EXPECT_FLOAT_EQ(60, Util::rad2deg(1.0/3*Util::pi));
   }
   TEST_F(UtilTest, getDistance) {
      EXPECT_FLOAT_EQ(0, Util::getDistance(60,10,60,10));
      EXPECT_FLOAT_EQ(20037508, Util::getDistance(90,10,-90,10));
      EXPECT_FLOAT_EQ(20037508, Util::getDistance(0,0,0,180));
      EXPECT_FLOAT_EQ(16879114, Util::getDistance(60.5,5.25,-84.75,-101.75));
      EXPECT_FLOAT_EQ(124080.79, Util::getDistance(60,10,61,11));
   }
   TEST_F(UtilTest, getDistance_approx) {
      EXPECT_FLOAT_EQ(0, Util::getDistance(60,10,60,10, true));
      EXPECT_NEAR(20037508, Util::getDistance(90,10,-90,10, true), 100);
      EXPECT_NEAR(20037508, Util::getDistance(0,0,0,180, true), 100);
      // EXPECT_NEAR(16879114, Util::getDistance(60.5,5.25,-84.75,-101.75, true), 100);
      EXPECT_NEAR(124084.21, Util::getDistance(60,10,61,11, true), 100);
   }
   TEST_F(UtilTest, getDistanceInvalid) {
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(Util::MV,5.25,-84.75,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(60.5,Util::MV,-84.75,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(60.5,5.25,Util::MV,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(60.5,5.25,-84.75,Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(Util::MV,5.25,Util::MV,-101.75));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance(Util::MV,Util::MV,Util::MV,Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::getDistance((float) 1/0,5.25,-84.75,-101.75));
   }
   TEST_F(UtilTest, sortFirst) {
      typedef std::pair<int,int> pair;
      std::vector<pair> pairs;
      pairs.push_back(pair(0,2));
      pairs.push_back(pair(1,4));
      pairs.push_back(pair(3,5));
      pairs.push_back(pair(2,3));
      std::sort(pairs.begin(), pairs.end(), Util::sort_pair_first<int,int>());
      EXPECT_EQ(0, pairs[0].first);
      EXPECT_EQ(1, pairs[1].first);
      EXPECT_EQ(2, pairs[2].first);
      EXPECT_EQ(3, pairs[3].first);
      EXPECT_EQ(2, pairs[0].second);
      EXPECT_EQ(4, pairs[1].second);
      EXPECT_EQ(3, pairs[2].second);
      EXPECT_EQ(5, pairs[3].second);
   }
   TEST_F(UtilTest, sortSecond) {
      typedef std::pair<int,int> pair;
      std::vector<pair> pairs;
      pairs.push_back(pair(0,2));
      pairs.push_back(pair(1,4));
      pairs.push_back(pair(3,5));
      pairs.push_back(pair(2,3));
      std::sort(pairs.begin(), pairs.end(), Util::sort_pair_second<int,int>());
      EXPECT_EQ(0, pairs[0].first);
      EXPECT_EQ(2, pairs[1].first);
      EXPECT_EQ(1, pairs[2].first);
      EXPECT_EQ(3, pairs[3].first);
      EXPECT_EQ(2, pairs[0].second);
      EXPECT_EQ(3, pairs[1].second);
      EXPECT_EQ(4, pairs[2].second);
      EXPECT_EQ(5, pairs[3].second);
   }
   TEST_F(UtilTest, calcDate) {
      EXPECT_EQ(20150101, Util::calcDate(20141231, 24));
      EXPECT_EQ(20141231, Util::calcDate(20141231, 0));
      EXPECT_EQ(20141231, Util::calcDate(20141231, 1.5));
      EXPECT_EQ(20150101, Util::calcDate(20141231, 47.9));
      EXPECT_EQ(20141231, Util::calcDate(20150101, -0.1));
      EXPECT_EQ(20141231, Util::calcDate(20150101, -24));
      EXPECT_EQ(20141230, Util::calcDate(20150101, -24.01));
   }
   TEST_F(UtilTest, split) {
      std::string s = "test 1 2 3";
      std::vector<std::string> strings = Util::split(s);
      ASSERT_EQ(4, strings.size());
      EXPECT_EQ("test", strings[0]);
      EXPECT_EQ("1", strings[1]);
      EXPECT_EQ("2", strings[2]);
      EXPECT_EQ("3", strings[3]);
   }
   TEST_F(UtilTest, splitEmpty) {
      std::string s = "";
      std::vector<std::string> strings = Util::split(s);
      EXPECT_EQ(0, strings.size());
   }
   TEST_F(UtilTest, splitEmptyExtraSpace) {
      std::string s = "     ";
      std::vector<std::string> strings = Util::split(s);
      EXPECT_EQ(0, strings.size());
   }
   TEST_F(UtilTest, splitExtraSpace) {
      std::string s = "  test  1 ";
      std::vector<std::string> strings = Util::split(s);
      ASSERT_EQ(2, strings.size());
      EXPECT_EQ("test", strings[0]);
      EXPECT_EQ("1", strings[1]);
   }
   TEST_F(UtilTest, splitDelim) {
      std::string s = "test 1,2, test";
      std::vector<std::string> strings = Util::split(s, ",");
      ASSERT_EQ(3, strings.size());
      EXPECT_EQ("test 1", strings[0]);
      EXPECT_EQ("2", strings[1]);
      EXPECT_EQ(" test", strings[2]);
   }
   TEST_F(UtilTest, splitDelims) {
      std::string s = "test 1,2q3 test,4q4";
      std::vector<std::string> strings = Util::split(s, "q,");
      ASSERT_EQ(5, strings.size());
      EXPECT_EQ("test 1", strings[0]);
      EXPECT_EQ("2", strings[1]);
      EXPECT_EQ("3 test", strings[2]);
      EXPECT_EQ("4", strings[3]);
      EXPECT_EQ("4", strings[4]);
   }
   TEST_F(UtilTest, logit) {
      const float p[] = {0.2,0.9};
      const float exp[] = {-1.386294, 2.197225};
      for(int i = 0; i < sizeof(p)/sizeof(float); i++) {
         float ans = Util::logit(p[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
      EXPECT_FLOAT_EQ(Util::MV, Util::logit((float) 1/0));
      EXPECT_FLOAT_EQ(Util::MV, Util::logit(-1));
      EXPECT_FLOAT_EQ(Util::MV, Util::logit(0));
      EXPECT_FLOAT_EQ(Util::MV, Util::logit(1));
   }
   TEST_F(UtilTest, invlogit) {
      const float x[] = {-2, 3};
      const float exp[] = {0.1192029, 0.9525741};
      for(int i = 0; i < sizeof(x)/sizeof(float); i++) {
         float ans = Util::invLogit(x[i]);
         EXPECT_FLOAT_EQ(exp[i], ans);
      }
      EXPECT_FLOAT_EQ(Util::MV, Util::invLogit((float) 1/0));
   }
   TEST_F(UtilTest, error) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Util::error("test"), ".*");
   }
   TEST_F(UtilTest, warning) {
      Util::setShowWarning(false);
      Util::warning("test");
   }
   TEST_F(UtilTest, status) {
      Util::setShowStatus(false);
      Util::status("test");
      Util::status("test", false);
      Util::status("test", true);
   }
   TEST_F(UtilTest, info) {
      Util::setShowInfo(false);
      Util::info("test");
   }
   TEST_F(UtilTest, getDate) {
      EXPECT_FLOAT_EQ(20150528, Util::getDate(1432816155));
      EXPECT_FLOAT_EQ(20001015, Util::getDate(971623971));
   }
   TEST_F(UtilTest, getTime) {
      EXPECT_FLOAT_EQ(122915, Util::getTime(1432816155));
      EXPECT_FLOAT_EQ(153251, Util::getTime(971623971));
   }
   TEST_F(UtilTest, getUnixTime) {
      EXPECT_FLOAT_EQ(1432816155, Util::getUnixTime(20150528, 122915));
      EXPECT_FLOAT_EQ(971623971,  Util::getUnixTime(20001015, 153251));
      EXPECT_FLOAT_EQ(1432771190, Util::getUnixTime(20150528, -600));
   }
   TEST_F(UtilTest, getCurrentDate) {
      int date = Util::getCurrentDate();
      EXPECT_GT(date, 20000000);
      EXPECT_LT(date, 99999999);
   }
   TEST_F(UtilTest, getCurrentTime) {
      EXPECT_GE(Util::getCurrentTime(), 0);
      EXPECT_LT(Util::getCurrentTime(), 240000);
   }
   TEST_F(UtilTest, getCurrentUnixTime) {
      // A rough check that we are after 2000/01/01 00:00:00
      EXPECT_GT(Util::getCurrentUnixTime(), 946684800);
   }
   TEST_F(UtilTest, getCurrentDateString) {
      std::string date = Util::getCurrentTimeStamp();
      EXPECT_EQ(19, date.size());
   }
   TEST_F(UtilTest, exists) {
      EXPECT_FALSE(Util::exists("testing/files/parametersw89ey8wehf.txt"));
      EXPECT_FALSE(Util::exists("testing/parametersw89ey8wehf.txt"));
      EXPECT_FALSE(Util::exists("testing/q/parametersw89ey8wehf.txt"));
      EXPECT_TRUE(Util::exists("testing/files/10x10.nc"));
   }
   TEST_F(UtilTest, copy) {
      // Check that copying a file works. Do this by modifying a file, and then overwritiing
      // it with known contents.
      bool status = Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
      EXPECT_TRUE(status);
      // Change 10x10_copy.nc
      Variable variable("air_temperature_2m");
      CalibratorNeighbourhood neighbourhood(variable, Options(""));
      FileNetcdf f1("testing/files/10x10.nc");
      FileNetcdf f2("testing/files/10x10_copy.nc");
      FieldPtr field1 = f1.getField(variable, 0);
      FieldPtr field2 = f2.getField(variable, 0);
      EXPECT_EQ((*field2)(1,1,0), (*field1)(1,1,0));
      neighbourhood.calibrate(f2);
      // Values should be different now
      EXPECT_NE((*field2)(1,1,0), (*field1)(1,1,0));
      std::vector<Variable> vars(1,variable);
      f2.write(vars);

      // Use copy. Values should be the same afterwards
      status = Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
      EXPECT_TRUE(status);
      FileNetcdf f3("testing/files/10x10_copy.nc");
      FieldPtr field3 = f3.getField(variable, 0);
      EXPECT_EQ(f1.getNumY(), f2.getNumY());
      EXPECT_EQ(f1.getNumX(), f2.getNumX());
      EXPECT_EQ(f1.getNumEns(), f2.getNumEns());
      EXPECT_EQ((*field3)(1,1,0), (*field1)(1,1,0));
   }
   TEST_F(UtilTest, copyInvalid) {
      bool status = Util::copy("testing/filesDoesNotExist/file1", "testing/filesDoesNotExist/file1");
      EXPECT_FALSE(status);
      status = Util::copy("testing/files/10x10.nc", "testing/filesDoesNotExist/file1");
      EXPECT_FALSE(status);
      status = Util::copy("testing/filesDoesNotExist/file1", "testing/files/10x10_copy.nc");
      EXPECT_FALSE(status);
      status = Util::copy("testing/files/10x10.nc", "testing/files/10x10_copy.nc");
      EXPECT_TRUE(status);
   }
   TEST_F(UtilTest, remove) {
      // Check that removing a file works.
      std::string tempfile = "testing/files/temp123.nc";
      bool status = Util::copy("testing/files/10x10.nc", tempfile);
      status = Util::remove(tempfile);
      EXPECT_TRUE(status);
      EXPECT_FALSE(Util::exists(tempfile));
   }
   TEST_F(UtilTest, removeInvalid) {
      // Check that removing a file works.
      std::string tempfile = "testing/files/temp123.nc";
      bool status = Util::remove(tempfile);
      EXPECT_FALSE(status);
      EXPECT_FALSE(Util::exists(tempfile));
   }
   TEST_F(UtilTest, hasChar) {
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#st", 't'));
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#st", ' '));
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#2t", '2'));
      EXPECT_TRUE(Util::hasChar("te  2384 &$*#2t", '&'));
      EXPECT_FALSE(Util::hasChar("te  2384 &$*#2t", ','));
      EXPECT_FALSE(Util::hasChar("te  2384 &$*#2t", 'q'));

      EXPECT_FALSE(Util::hasChar("", 'q'));
      EXPECT_FALSE(Util::hasChar("", ' '));

      EXPECT_TRUE(Util::hasChar("test\nq2", '2'));
      EXPECT_TRUE(Util::hasChar("test\nq2", 't'));
      EXPECT_TRUE(Util::hasChar("test\nq2", '\n'));
      EXPECT_FALSE(Util::hasChar("test\nq2", ' '));
      EXPECT_FALSE(Util::hasChar("test\nq2", '3'));
   }
   TEST_F(UtilTest, formatDescription) {
      EXPECT_EQ("test      ad qwi qwio\n          wqio dwqion\n          qdwion", Util::formatDescription("test", "ad qwi qwio wqio dwqion qdwion", 10, 22, 0));
      EXPECT_EQ("  test    ad qwi qwio\n          wqio dwqion\n          qdwion", Util::formatDescription("test", "ad qwi qwio wqio dwqion qdwion", 10, 22, 2));
      // Invalid input should not throw an error
      Util::formatDescription("test", "ad qwi qwio wqio dwqion qdwion", 10, 5, 2); // Too narrow message
      Util::formatDescription("test", "ad qwi qwio wqio dwqion qdwion", 10, 11, 2); // Very narrow message
   }
   TEST_F(UtilTest, compute) {
      FileNetcdf from("testing/files/10x10.nc");

      // Odd-sized neighbourhood
      std::vector<float> hood = getHood(3,0,2);
      EXPECT_FLOAT_EQ(2, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(0, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(3, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(1.6666666, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(1.247219, Util::calculateStat(hood, Util::StatTypeStd));

      // Even-sized neighbourhood
      hood = getHood(4,2,0,3);
      EXPECT_FLOAT_EQ(2.5, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(0, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(4, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(2.25, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(1.47902, Util::calculateStat(hood, Util::StatTypeStd));

      // Missing value
      hood = getHood(0,Util::MV,-4,2);
      EXPECT_FLOAT_EQ(0, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(-4, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(2, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(-0.6666666, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(2.494438, Util::calculateStat(hood, Util::StatTypeStd));

      // All same values
      hood = getHood(1,1,1,1);
      EXPECT_FLOAT_EQ(1, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(1, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(1, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(1, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(0, Util::calculateStat(hood, Util::StatTypeStd));

      // Small variance, large mean
      hood = getHood(28123.49,28123.48,28123.49);
      EXPECT_FLOAT_EQ(28123.49, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(28123.48, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(28123.49, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(28123.48666667, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(0.0046035596, Util::calculateStat(hood, Util::StatTypeStd));

      // All missing
      hood = getHood(Util::MV,Util::MV,Util::MV);
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeStd));

      // Empty
      hood = std::vector<float>();
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMedian, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMin, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMax, Util::MV));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeMean));
      EXPECT_FLOAT_EQ(Util::MV, Util::calculateStat(hood, Util::StatTypeStd));
   }
   TEST_F(UtilTest, inverse) {
      // Inverse of identity matrix
      {
         vec2 matrix = get2x2Matrix(1,0,0,1);
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(2, inverse.size());
         ASSERT_EQ(2, inverse[0].size());
         EXPECT_FLOAT_EQ(1, inverse[0][0]);
         EXPECT_FLOAT_EQ(0, inverse[0][1]);
         EXPECT_FLOAT_EQ(0, inverse[1][0]);
         EXPECT_FLOAT_EQ(1, inverse[1][1]);
      }
      // Inverse of a complicated matrix
      {
         vec2 matrix = get2x2Matrix(0.2,3.1,-0.2,4.1);
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(2, inverse.size());
         ASSERT_EQ(2, inverse[0].size());
         EXPECT_FLOAT_EQ(2.84722222, inverse[0][0]);
         EXPECT_FLOAT_EQ(-2.15277778, inverse[0][1]);
         EXPECT_FLOAT_EQ(0.13888889, inverse[1][0]);
         EXPECT_FLOAT_EQ(0.13888889, inverse[1][1]);
      }
      // Inverse of a 1x1 matrix
      {
         vec2 matrix;
         matrix.resize(1);
         matrix[0].resize(1);
         matrix[0][0] = -2.3;
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(1, inverse.size());
         ASSERT_EQ(1, inverse[0].size());
         EXPECT_FLOAT_EQ(-0.4347826, inverse[0][0]);
      }
      // Empty matrix
      {
         vec2 matrix;
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(0, inverse.size());
      }
      // With missing values
      {
         vec2 matrix = get2x2Matrix(0.2,Util::MV,-0.2,4.1);
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(2, inverse.size());
         ASSERT_EQ(2, inverse[0].size());
         EXPECT_FLOAT_EQ(Util::MV, inverse[0][0]);
         EXPECT_FLOAT_EQ(Util::MV, inverse[0][1]);
         EXPECT_FLOAT_EQ(Util::MV, inverse[1][0]);
         EXPECT_FLOAT_EQ(Util::MV, inverse[1][1]);
      }
      // Unstable matrix
      {
         vec2 matrix = get2x2Matrix(0,0,0,0);
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(2, inverse.size());
         ASSERT_EQ(2, inverse[0].size());
         EXPECT_FLOAT_EQ(Util::MV, inverse[0][0]);
         EXPECT_FLOAT_EQ(Util::MV, inverse[0][1]);
         EXPECT_FLOAT_EQ(Util::MV, inverse[1][0]);
         EXPECT_FLOAT_EQ(Util::MV, inverse[1][1]);
      }
      // 3x3
      {
         vec2 matrix = get3x3Matrix(1,0,0.2,0,1,0.2,0.2,0.2,1);
         vec2 inverse = Util::inverse(matrix);
         ASSERT_EQ(3, inverse.size());
         ASSERT_EQ(3, inverse[0].size());
         EXPECT_FLOAT_EQ(1.04347826, inverse[0][0]);
         EXPECT_FLOAT_EQ(0.04347826, inverse[0][1]);
         EXPECT_FLOAT_EQ(-0.2173913, inverse[0][2]);
         EXPECT_FLOAT_EQ(0.04347826, inverse[1][0]);
         EXPECT_FLOAT_EQ(1.04347826, inverse[1][1]);
         EXPECT_FLOAT_EQ(-0.2173913, inverse[1][2]);
         EXPECT_FLOAT_EQ(-0.21739130, inverse[2][0]);
         EXPECT_FLOAT_EQ(-0.21739130, inverse[2][1]);
         EXPECT_FLOAT_EQ(1.0869565, inverse[2][2]);
      }
   }
   TEST_F(UtilTest, gridppVersion) {
      std::string version = Util::gridppVersion();
      EXPECT_NE("", version);
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
