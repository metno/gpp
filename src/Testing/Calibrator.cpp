#include "../Calibrator/Calibrator.h"
#include "../Util.h"
#include "../Options.h"
#include <gtest/gtest.h>
#include <vector>

namespace {
   class TestCalibrator : public ::testing::Test {
      public:
         std::vector<float> vec(const float iVector[], int iSize) {
            std::vector<float> vec(iSize, 0);
            for(int i = 0; i < iSize; i++) {
               vec[i] = iVector[i];
            }
            return vec;
         };
      protected:
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(TestCalibrator, shuffle) {
      const std::vector<float> vec1 = vec((const float[]) {5,1,4,7,6,2,3}, 7);
            std::vector<float> vec2 = vec((const float[]) {32,14,21,0,11,2,5}, 7);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(7, vec1.size());
      ASSERT_EQ(7, vec2.size());
      EXPECT_EQ(14, vec2[0]);
      EXPECT_EQ(0, vec2[1]);
      EXPECT_EQ(11, vec2[2]);
      EXPECT_EQ(32, vec2[3]);
      EXPECT_EQ(21, vec2[4]);
      EXPECT_EQ(2, vec2[5]);
      EXPECT_EQ(5, vec2[6]);
   }
   TEST_F(TestCalibrator, shuffleDifferentSizes) {
      const std::vector<float> vec1 = vec((const float[]) {5,4}, 2);
            std::vector<float> vec2 = vec((const float[]) {1,2,3}, 3);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(2, vec1.size());
      ASSERT_EQ(3, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(3, vec2[2]);
   }
   TEST_F(TestCalibrator, shuffleZeroSize) {
      const std::vector<float> vec1;
            std::vector<float> vec2 = vec((const float[]) {1,2,3}, 3);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(0, vec1.size());
      ASSERT_EQ(3, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(3, vec2[2]);
   }
   TEST_F(TestCalibrator, shuffleZeroSize2) {
      const std::vector<float> vec1;
            std::vector<float> vec2;
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(0, vec1.size());
      ASSERT_EQ(0, vec2.size());
   }
   TEST_F(TestCalibrator, shuffleWithMissing) {
      float ar1[] = {3, Util::MV, 19, 3};
      float ar2[] = {1,2,4,3};
      const std::vector<float> vec1 = vec(ar1, 4);
            std::vector<float> vec2 = vec(ar2, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(4, vec2[2]);
      EXPECT_EQ(3, vec2[3]);
   }
   TEST_F(TestCalibrator, shuffleWithMissing2) {
      float ar1[] = {3, 1, 19,3};
      float ar2[] = {1,2,Util::MV,3};
      const std::vector<float> vec1 = vec(ar1, 4);
            std::vector<float> vec2 = vec(ar2, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      // Should remain unchanged
      EXPECT_EQ(1, vec2[0]);
      EXPECT_EQ(2, vec2[1]);
      EXPECT_EQ(Util::MV, vec2[2]);
      EXPECT_EQ(3, vec2[3]);
   }
   TEST_F(TestCalibrator, shuffleRepeated) {
      const std::vector<float> vec1 = vec((const float[]) {3,1,7,1}, 4);
            std::vector<float> vec2 = vec((const float[]) {1,2,4,3}, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      EXPECT_EQ(3, vec2[0]);
      EXPECT_TRUE(vec2[1]==1 || vec2[1]==2);
      EXPECT_EQ(4, vec2[2]);
      EXPECT_TRUE(vec2[3]==1 || vec2[3]==2);
   }
   TEST_F(TestCalibrator, shuffleRepeated2) {
      const std::vector<float> vec1 = vec((const float[]) {3,2,7,1}, 4);
            std::vector<float> vec2 = vec((const float[]) {1,2,1,3}, 4);
      Calibrator::shuffle(vec1, vec2);
      EXPECT_EQ(4, vec1.size());
      ASSERT_EQ(4, vec2.size());
      EXPECT_EQ(2, vec2[0]);
      EXPECT_TRUE(vec2[1]==1 || vec2[1]==2);
      EXPECT_EQ(3, vec2[2]);
      EXPECT_TRUE(vec2[3]==1 || vec2[3]==2);
   }
   TEST_F(TestCalibrator, factoryZaga) {
      {
         Calibrator* c;
         c = Calibrator::getScheme("zaga", Variable("precipitation_amount"), Options("popThreshold=0.24 fracThreshold=0.34 neighbourhoodSize=24 maxEnsMean=90"));
         EXPECT_TRUE(c);
         EXPECT_EQ("zaga", c->name());
         delete c;
      }
      {
         Calibrator* c;
         c = Calibrator::getScheme("zaga", mVariable, Options("popThreshold=-0.12 outputPop=0 fracThreshold=-0.92 neighbourhoodSize=6 maxEnsMean=40"));
         EXPECT_TRUE(c);
         EXPECT_EQ("zaga", c->name());
         delete c;
      }
   }
   TEST_F(TestCalibrator, factoryNeighbourhood) {
      Calibrator* c;
      c = Calibrator::getScheme("neighbourhood",mVariable, Options("radius=3"));
      EXPECT_TRUE(c);
      EXPECT_EQ("neighbourhood", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryQc) {
      Calibrator* c;
      c = Calibrator::getScheme("qc", mVariable, Options("m3"));
      EXPECT_TRUE(c);
      EXPECT_EQ("qc", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryPhase) {
      Calibrator* c;
      c = Calibrator::getScheme("phase", mVariable, Options("temperature=air_temperature_2m precipitation=precipitation_amount minPrecip=0.771 useWetbulb=0"));
      EXPECT_TRUE(c);
      EXPECT_EQ("phase", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryGaussian) {
      Calibrator* c;
      c = Calibrator::getScheme("gaussian", mVariable, Options(""));
      EXPECT_TRUE(c);
      EXPECT_EQ("gaussian", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryRegression) {
      Calibrator* c;
      c = Calibrator::getScheme("regression", mVariable, Options(""));
      EXPECT_TRUE(c);
      EXPECT_EQ("regression", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryQnh) {
      Calibrator* c;
      c = Calibrator::getScheme("qnh", mVariable, Options("pressureVariable=surface_air_pressure"));
      EXPECT_TRUE(c);
      EXPECT_EQ("qnh", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryWindow) {
      Calibrator* c;
      c = Calibrator::getScheme("window", mVariable, Options("radius=2 stat=quantile quantile=0.5"));
      EXPECT_TRUE(c);
      EXPECT_EQ("window", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryKriging) {
      Calibrator* c;
      c = Calibrator::getScheme("kriging", mVariable, Options("radius=100 maxElevDiff=100 efoldDist=2"));
      EXPECT_TRUE(c);
      EXPECT_EQ("kriging", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryQq) {
      Calibrator* c;
      c = Calibrator::getScheme("qq", mVariable, Options(""));
      EXPECT_TRUE(c);
      EXPECT_EQ("qq", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factorySort) {
      Calibrator* c;
      c = Calibrator::getScheme("sort", mVariable, Options(""));
      EXPECT_TRUE(c);
      EXPECT_EQ("sort", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryBct) {
      Calibrator* c;
      c = Calibrator::getScheme("bct", mVariable, Options(""));
      EXPECT_TRUE(c);
      EXPECT_EQ("bct", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryAltitude) {
      Calibrator* c;
      c = Calibrator::getScheme("altitude", mVariable, Options());
      EXPECT_TRUE(c);
      EXPECT_EQ("altitude", c->name());
      delete c;
   }
   TEST_F(TestCalibrator, factoryValid) {
      Calibrator::getScheme("zaga", mVariable, Options(""));
      Calibrator::getScheme("zaga", mVariable, Options("variable=T"));
      Calibrator::getScheme("neighbourhood", mVariable, Options(""));
   }
   TEST_F(TestCalibrator, descriptions) {
      std::string descriptions = Calibrator::getDescriptions();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
