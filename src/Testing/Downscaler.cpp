#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscaler : public ::testing::Test {
      public:
         vec2 makeVec2(int nLat, int nLon, const std::vector<float>& values) {
            vec2 grid;
            grid.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               grid[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  int index = i*nLon + j;
                  grid[i][j] = values[index];
               }
            }
            return grid;
         };
         void setLatLon(FileFake& iFile, const float iLat[], const float iLon[]) {
            vec2 lat;
            vec2 lon;
            int nLat = iFile.getNumY(); 
            int nLon = iFile.getNumX();
            lat.resize(nLat);
            lon.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               lat[i].resize(nLon);
               lon[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  lat[i][j] = iLat[i];
                  lon[i][j] = iLon[j];
               }
            }
            iFile.setLats(lat);
            iFile.setLons(lon);
         };
      protected:
         virtual void SetUp() {
             mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(TestDownscaler, validDownscalers) {
      Downscaler* d0 = Downscaler::getScheme("nearestNeighbour", mVariable, mVariable, Options());
      Downscaler* d1 = Downscaler::getScheme("smart", mVariable, mVariable, Options("elevRadius=3 numSmart=2 minElevDiff=400"));
      Downscaler* d2 = Downscaler::getScheme("gradient", mVariable, mVariable, Options("elevRadius=5 constantElevGradient=0.04 minElevDiff=213.2"));
      Downscaler* d3 = Downscaler::getScheme("pressure", mVariable, mVariable, Options(""));
      Downscaler* d4 = Downscaler::getScheme("bypass", mVariable, mVariable, Options(""));
      EXPECT_EQ("nearestNeighbour", d0->name());
      EXPECT_EQ("smart", d1->name());
      EXPECT_EQ("gradient", d2->name());
      EXPECT_EQ("pressure", d3->name());
      EXPECT_EQ("bypass", d4->name());
   }
   TEST_F(TestDownscaler, invalidDownscalers) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(Downscaler::getScheme("woehiowciwofew", mVariable, mVariable, Options("searchRadius=5")), ".*");
      EXPECT_DEATH(Downscaler::getScheme("woehiowciwofew", mVariable, mVariable, Options()), ".*");
   }
   TEST_F(TestDownscaler, validNearestNeighbours) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {50,55,60}, (const float[]){0,10});
      setLatLon(to,   (const float[]) {40, 54.99},   (const float[]){-1,9.99});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::clearCache();
      Downscaler::getNearestNeighbourBruteForce(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());

      EXPECT_EQ(0, I[0][0]);
      EXPECT_EQ(0, J[0][0]);
      EXPECT_EQ(0, I[0][1]);
      EXPECT_EQ(1, J[0][1]);
      EXPECT_EQ(1, I[1][0]);
      EXPECT_EQ(0, J[1][0]);
      EXPECT_EQ(1, I[1][1]);
      EXPECT_EQ(1, J[1][1]);
   }
   TEST_F(TestDownscaler, identicalGrid) {
      int nLat = 3;
      int nLon = 4;
      FileFake from(Options("nLat=3 nLon=4 nEns=1 nTime=1"));
      FileFake to(Options("nLat=3 nLon=4 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {50,55,60}, (const float[]){0,3,5,10});
      setLatLon(to,   (const float[]) {50,55,60}, (const float[]){0,3,5,10});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::clearCache();
      Downscaler::getNearestNeighbourBruteForce(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(nLat, I.size());
      ASSERT_EQ(nLon, I[0].size());

      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            EXPECT_EQ(i, I[i][j]);
            EXPECT_EQ(j, J[i][j]);
         }
      }
   }
   TEST_F(TestDownscaler, notIdenticalGridDiffLon) {
      int nLat = 1;
      int nLon = 5;
      FileFake from(Options("nLat=1 nLon=5 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=5 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {1}, (const float[]){0, 1, 2, 3, 4});
      setLatLon(to,   (const float[]) {1}, (const float[]){0, 2, 2.4, 3, 4});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::clearCache();
      Downscaler::getNearestNeighbourBruteForce(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(nLat, I.size());
      ASSERT_EQ(nLon, I[0].size());

      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            EXPECT_EQ(i, I[i][j]);
         }
      }
   }
   TEST_F(TestDownscaler, notIdenticalGridDiffLat) {
      int nLat = 5;
      int nLon = 1;
      FileFake from(Options("nLat=5 nLon=1 nEns=1 nTime=1"));
      FileFake to(Options("nLat=5 nLon=1 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {0, 1, 2, 3, 4}, (const float[]){1});
      setLatLon(to,   (const float[]) {0, 2, 2.4, 3, 4}, (const float[]){1});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::clearCache();
      Downscaler::getNearestNeighbourBruteForce(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(nLat, I.size());
      ASSERT_EQ(nLon, I[0].size());

      for(int i = 0; i < nLat; i++) {
         for(int j = 0; j < nLon; j++) {
            EXPECT_EQ(j, J[i][j]);
         }
      }
   }
   TEST_F(TestDownscaler, unsortedGrid) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {60,50,55}, (const float[]){5,4});
      setLatLon(to,   (const float[]) {56,49},    (const float[]){3,4.6});

      vec2Int I, J, If, Jf;
      Downscaler::getNearestNeighbour(from, to, I, J);
      Downscaler::clearCache();
      Downscaler::getNearestNeighbourBruteForce(from, to, If, Jf);
      EXPECT_EQ(I, If);
      EXPECT_EQ(J, Jf);

      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());

      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(1, J[0][0]);
      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(0, J[0][1]);
      EXPECT_EQ(1, I[1][0]);
      EXPECT_EQ(1, J[1][0]);
      EXPECT_EQ(1, I[1][0]);
      EXPECT_EQ(0, J[1][1]);
   }
   TEST_F(TestDownscaler, cache) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {60,50,55}, (const float[]){5,4});
      setLatLon(to,   (const float[]) {56,49},    (const float[]){3,4.6});

      vec2Int I, J;
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());
      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(1, J[0][0]);

      // Don't change the grid
      Uuid idBefore = from.getUniqueTag();
      setLatLon(from, (const float[]) {60,50,55}, (const float[]){5,4});
      Uuid idAfter = from.getUniqueTag();
      EXPECT_EQ(idBefore, idAfter);
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());
      EXPECT_EQ(2, I[0][0]);
      EXPECT_EQ(1, J[0][0]);

      // Change the grid
      idBefore = from.getUniqueTag();
      setLatLon(from, (const float[]) {60,55,50}, (const float[]){5,4});
      idAfter = from.getUniqueTag();
      EXPECT_NE(idBefore, idAfter);
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(2, I.size());
      ASSERT_EQ(2, I[0].size());
      EXPECT_EQ(1, I[0][0]);
      EXPECT_EQ(1, J[0][0]);
   }
   TEST_F(TestDownscaler, copyConstructor) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to = from;
      vec2Int I, J;
      EXPECT_EQ(from.getUniqueTag(), to.getUniqueTag());
      Downscaler::getNearestNeighbour(from, to, I, J);
      ASSERT_EQ(3, I.size());
      ASSERT_EQ(2, I[0].size());
      for(int i = 0; i < I.size(); i++) {
         for(int j = 0; j < I[0].size(); j++) {
            EXPECT_EQ(i, I[i][j]);
            EXPECT_EQ(j, J[i][j]);
         }
      }
   }
   TEST_F(TestDownscaler, missingLatLon) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      float lat1[] = {50,Util::MV,60};
      float lon1[] = {0,Util::MV};
      float lat2[] = {40, 54};
      float lon2[] = {-1,Util::MV};
      EXPECT_DEATH(setLatLon(from, lat1, lon1), ".*");
      EXPECT_DEATH(setLatLon(to,   lat2, lon2), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
