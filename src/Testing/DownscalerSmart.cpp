#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>
#include <boost/assign/list_of.hpp>

namespace {
   class TestDownscalerSmart : public ::testing::Test {
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
         // iElev must be a flat array with lon varying fastest
         void setLatLonElev(FileFake& iFile, const float iLat[], const float iLon[], const float iElev[]) {
            vec2 lat;
            vec2 lon;
            vec2 elev;
            int nLat = iFile.getNumY(); 
            int nLon = iFile.getNumX();
            lat.resize(nLat);
            lon.resize(nLat);
            elev.resize(nLat);
            for(int i = 0; i < nLat; i++) {
               lat[i].resize(nLon);
               lon[i].resize(nLon);
               elev[i].resize(nLon);
               for(int j = 0; j < nLon; j++) {
                  lat[i][j] = iLat[i];
                  lon[i][j] = iLon[j];
                  elev[i][j] = iElev[j+i*nLon];
               }
            }
            iFile.setLats(lat);
            iFile.setLons(lon);
            iFile.setElevs(elev);
         };
      protected:
         virtual void SetUp() {
             mT = Variable("air_temperature_2m");
             mPrecip = Variable("precipitation_amount");
         }
         virtual void TearDown() {
         }
         Variable mT;
         Variable mPrecip;
   };

   TEST_F(TestDownscalerSmart, isValid) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      setLatLonElev(from, (const float[]) {50,55,60}, (const float[]){0,10}, (const float[]){3, 15, 6, 30, 20, 11});
      setLatLonElev(to,   (const float[]) {54},   (const float[]){9}, (const float[]){10});

      vec3Int I, J;
      DownscalerSmart d(mPrecip, mPrecip, Options("radius=10 num=2"));
      d.getSmartNeighbours(from, to, I, J);

      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(2, I[0][0].size());

      EXPECT_EQ(2, I[0][0][0]);
      EXPECT_EQ(1, J[0][0][0]);
      EXPECT_EQ(1, I[0][0][1]);
      EXPECT_EQ(0, J[0][0][1]);
   }
   TEST_F(TestDownscalerSmart, 10x10) {
      DownscalerSmart d(mT, mT, Options("radius=3 num=2"));
      FileNetcdf from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(mT, 0);
      std::stringstream ss;
      ss << "nLat=1 nLon=5 nEns=1 nTime=" << from.getNumTime();
      FileFake to(Options(ss.str()));
      // Case 1: West boundary outside domain
      // Case 2: Within domain
      // Case 3/4: Nearest neighbour is on the boundary, so only the western half is used
      // Case 5: Elevation is missing, so use the nearest neighbour
      float elev[] = {120, 80, 600, 600, Util::MV};
      setLatLonElev(to, (const float[]) {5}, (const float[]){2,3,12,20,2}, elev);
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(mT, 0);
      ASSERT_EQ(1, toT.getNumY());
      ASSERT_EQ(5, toT.getNumX());
      EXPECT_FLOAT_EQ(303,   toT(0,0,0));
      EXPECT_FLOAT_EQ(304.5, toT(0,1,0));
      EXPECT_FLOAT_EQ(305.5, toT(0,2,0));
      EXPECT_FLOAT_EQ(305.5, toT(0,3,0));
      EXPECT_FLOAT_EQ(301,   toT(0,4,0));
      vec3Int I, J;
      d.getSmartNeighbours(from, to, I, J);
      EXPECT_EQ(4, I[0][2][0]);
      EXPECT_EQ(9, J[0][2][0]);
      EXPECT_EQ(3, I[0][2][1]);
      EXPECT_EQ(8, J[0][2][1]);

      ASSERT_EQ(1, I[0][4].size());
      ASSERT_EQ(1, J[0][4].size());
      EXPECT_EQ(5, I[0][4][0]);
      EXPECT_EQ(2, J[0][4][0]);
   }
   TEST_F(TestDownscalerSmart, 10x10minElevDiff) {
      DownscalerSmart d(mT, mT, Options("radius=3 num=2 minElevDiff=109"));
      FileNetcdf from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(mT, 0);
      std::stringstream ss;
      ss << "nLat=1 nLon=3 nEns=1 nTime=" << from.getNumTime();
      FileFake to(Options(ss.str()));
      float elev[] = {120, 50, Util::MV};
      setLatLonElev(to, (const float[]) {5}, (const float[]){2,2,2}, elev);
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(mT, 0);
      ASSERT_EQ(1, toT.getNumY());
      ASSERT_EQ(3, toT.getNumX());
      EXPECT_FLOAT_EQ(301, toT(0,0,0)); // Nearest neighbour
      EXPECT_FLOAT_EQ(304, toT(0,1,0));
      EXPECT_FLOAT_EQ(301, toT(0,2,0)); // Nearest neighbour
      vec3Int I, J;
      d.getSmartNeighbours(from, to, I, J);
      ASSERT_EQ(2, I[0][1].size());
      ASSERT_EQ(2, J[0][1].size());
      EXPECT_EQ(4, I[0][1][0]);
      EXPECT_EQ(5, J[0][1][0]);
      EXPECT_EQ(2, I[0][1][1]);
      EXPECT_EQ(4, J[0][1][1]);

      // All nearest neighbours
      DownscalerSmart d2(mT, mT, Options("radius=3 num=2 minElevDiff=1200"));
      status = d2.downscale(from, to);
      EXPECT_TRUE(status);
      EXPECT_FLOAT_EQ(301, toT(0,0,0));
      EXPECT_FLOAT_EQ(301, toT(0,1,0));
      EXPECT_FLOAT_EQ(301, toT(0,2,0));
   }
   TEST_F(TestDownscalerSmart, minElevDiff) {
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      setLatLonElev(from, (const float[]) {50,55,60}, (const float[]){0,10}, (const float[]){3, 15, 6, 30, 20, 11});
      setLatLonElev(to,   (const float[]) {55},       (const float[]){10},   (const float[]){22.3});

      DownscalerSmart d(mT, mT, Options("radius=3 num=2 minElevDiff=7.701"));
      vec3Int I, J;

      // Use nearest neighbour (at 30m)
      d.getSmartNeighbours(from, to, I, J);
      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(1, I[0][0].size());
      EXPECT_FLOAT_EQ(1, I[0][0][0]);
      EXPECT_FLOAT_EQ(1, J[0][0][0]);

      // Use best neighbours (at 20m and 15m)
      DownscalerSmart d2(mT, mT, Options("radius=3 num=2 minElevDiff=7.6999"));
      d2.getSmartNeighbours(from, to, I, J);
      ASSERT_EQ(1, I[0].size());
      ASSERT_EQ(2, I[0][0].size());
      EXPECT_FLOAT_EQ(2, I[0][0][0]);
      EXPECT_FLOAT_EQ(0, J[0][0][0]);
      EXPECT_FLOAT_EQ(0, I[0][0][1]);
      EXPECT_FLOAT_EQ(1, J[0][0][1]);

   }
   TEST_F(TestDownscalerSmart, fewerPointsThanSmart) {
      FileFake from(Options("nLat=5 nLon=3 nEns=1 nTime=1"));
      setLatLonElev(from, (const float[]) {4,5,6,7,8}, (const float[]){5,10,15}, (const float[]){70,50,20,80,60,70,50,40,30,20,10,40,50,30,60});
      FileFake to(Options("nLat=1 nLon=4 nEns=1 nTime=1"));
      setLatLonElev(to, (const float[]) {5.5}, (const float[]){2,4, 10,20}, (const float[]){120, 80, 600, 600});
      vec3Int I, J;

      DownscalerSmart d(mPrecip, mPrecip, Options("radius=1 num=20"));
      d.getSmartNeighbours(from, to, I, J);

      ASSERT_EQ(1, I.size());
      ASSERT_EQ(1, J.size());
      ASSERT_EQ(4, I[0].size());
      ASSERT_EQ(4, J[0].size());
      EXPECT_EQ(6, I[0][0].size());
      EXPECT_EQ(6, J[0][0].size());
      EXPECT_EQ(6, I[0][1].size());
      EXPECT_EQ(6, J[0][1].size());
      EXPECT_EQ(9, I[0][2].size());
      EXPECT_EQ(9, J[0][2].size());
      EXPECT_EQ(6, I[0][3].size());
      EXPECT_EQ(6, J[0][3].size());
   }
   TEST_F(TestDownscalerSmart, description) {
      DownscalerSmart::description();
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
