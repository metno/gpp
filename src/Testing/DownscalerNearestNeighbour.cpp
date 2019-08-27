#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class TestDownscalerNearestNeighbour : public ::testing::Test {
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
         TestDownscalerNearestNeighbour() {
         }
         virtual ~TestDownscalerNearestNeighbour() {
         }
         virtual void SetUp() {
            mVariable = Variable("air_temperature_2m");
         }
         virtual void TearDown() {
         }
         Variable mVariable;
   };

   TEST_F(TestDownscalerNearestNeighbour, description) {
      DownscalerNearestNeighbour::description();
   }
   TEST_F(TestDownscalerNearestNeighbour, downscale) {
      DownscalerNearestNeighbour d(mVariable, mVariable, Options());
      FileFake from(Options("nLat=3 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {60,50,55}, (const float[]){5,4});
      setLatLon(to,   (const float[]) {56,49},    (const float[]){3,4.6});
      d.downscale(from, to);
      const Field& fromT = *from.getField(mVariable, 0);
      const Field& toT   = *to.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(fromT(2,1,0), toT(0,0,0));
      EXPECT_FLOAT_EQ(fromT(2,0,0), toT(0,1,0));
      EXPECT_FLOAT_EQ(fromT(1,1,0), toT(1,0,0));
      EXPECT_FLOAT_EQ(fromT(1,0,0), toT(1,1,0));
   }
   TEST_F(TestDownscalerNearestNeighbour, 10x10) {
      DownscalerNearestNeighbour d(mVariable, mVariable, Options());
      FileNetcdf from("testing/files/10x10.nc");
      const Field& fromT  = *from.getField(mVariable, 0);
      std::stringstream ss;
      ss << "nLat=2 nLon=2 nEns=1 nTime=" << from.getNumTime();
      FileFake to(Options(ss.str()));
      setLatLon(to, (const float[]) {5,11}, (const float[]){2,3});
      bool status = d.downscale(from, to);
      EXPECT_TRUE(status);
      const Field& toT   = *to.getField(mVariable, 0);
      ASSERT_EQ(2, toT.getNumY());
      ASSERT_EQ(2, toT.getNumX());
      EXPECT_FLOAT_EQ(301, toT(0,0,0));
      EXPECT_FLOAT_EQ(302, toT(0,1,0));
      EXPECT_FLOAT_EQ(309, toT(1,0,0));
      EXPECT_FLOAT_EQ(301, toT(1,1,0));
      vec2Int I, J;
      d.getNearestNeighbour(from, to, I, J);
      EXPECT_EQ(9, I[1][1]);
      EXPECT_EQ(3, J[1][1]);
      EXPECT_FLOAT_EQ(301, fromT(9,3,0));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
