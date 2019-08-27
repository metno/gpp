#include "../Util.h"
#include "../File/File.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class TestDownscalerBilinear : public ::testing::Test {
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
         vec2 make2x2(int i00, int i01, int i10, int i11) {
            vec2 grid;
            grid.resize(2);
            grid[0].resize(2);
            grid[1].resize(2);
            grid[0][0] = i00;
            grid[0][1] = i01;
            grid[1][0] = i10;
            grid[1][1] = i11;
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

   TEST_F(TestDownscalerBilinear, description) {
      DownscalerBilinear::description();
   }
   TEST_F(TestDownscalerBilinear, rectangular) {
      // Test in python
      // i = scipy.interpolate.interp2d([0,0,2,2], [0,1,0,1],[0,1,2,3]); i(0.123,0.897)
      float x0 = 0;
      float x1 = 0;
      float x2 = 2;
      float x3 = 2;
      float y0 = 0;
      float y1 = 1;
      float y2 = 0;
      float y3 = 1;
      float v0 = 0;
      float v1 = 1;
      float v2 = 2;
      float v3 = 3;
      // 2x1 rectangle
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(1.5, DownscalerBilinear::bilinear(1,0.5, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(1.143, DownscalerBilinear::bilinear(0.246,0.897, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
   }
   TEST_F(TestDownscalerBilinear, parallelogram) {
      float x0 = 0;
      float x1 = 0;
      float x2 = 2;
      float x3 = 2;
      float y0 = 0;
      float y1 = 1;
      float y2 = 0;
      float y3 = 1;
      float v0 = 0;
      float v1 = 1;
      float v2 = 2;
      float v3 = 3;
      x1 = 0.2;
      x3 = 2.2;
      // Parallellogram (half-rectangle)
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(1.5, DownscalerBilinear::bilinear(1.1,0.5, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(0.9636, DownscalerBilinear::bilinear(0.246,0.897, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      y2 = 0.5;
      y3 = 1.5;
      // Parallellogram
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(1.5, DownscalerBilinear::bilinear(1.1,0.75, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(0.94957895, DownscalerBilinear::bilinear(0.246,0.897, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
   }
   TEST_F(TestDownscalerBilinear, diamond) {
      float x0 = -1;
      float x1 = 0;
      float x2 = 0;
      float x3 = 1;
      float y0 = 0;
      float y1 = 1;
      float y2 = -1;
      float y3 = 0;
      float v0 = 0;
      float v1 = 1;
      float v2 = 2;
      float v3 = 3;
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(1.5, DownscalerBilinear::bilinear(0,0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(1.0825, DownscalerBilinear::bilinear(-0.146,0.397, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
   }
   TEST_F(TestDownscalerBilinear, general) {
      float x0 = 0;
      float x1 = 1;
      float x2 = 4;
      float x3 = 5;
      float y0 = 0;
      float y1 = 3;
      float y2 = 0;
      float y3 = 5;
      float v0 = 0;
      float v1 = 1;
      float v2 = 2;
      float v3 = 3;
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      float tol = 0.05;
      // The following do not give identical results, maybe scipy uses a slightly different formulation
      EXPECT_NEAR(1.68666667, DownscalerBilinear::bilinear(3,1.4, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3), tol);
      EXPECT_NEAR(2.29005863, DownscalerBilinear::bilinear(3.87123,2.98321, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3), tol);

      // EXPECT_NEAR(2.74033847, DownscalerBilinear::bilinear(4.87123,2.98321, x2, x3, x0, x1, y2, y3, y0, y1, v2, v3, v0, v1), tol);
   }
   TEST_F(TestDownscalerBilinear, trapezoidVerticalParallel) {
      float x0 = 0;
      float x1 = 0;
      float x2 = 1;
      float x3 = 1;
      float y0 = 0;
      float y1 = 1;
      float y2 = 0.1;
      float y3 = 0.9;
      float v0 = 0;
      float v1 = 1;
      float v2 = 2;
      float v3 = 3;
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      float tol = 0.03;
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_NEAR(1.5, DownscalerBilinear::bilinear(0.5,0.5, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3), tol);
      EXPECT_NEAR(1.95175825, DownscalerBilinear::bilinear(0.903,0.211, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3), tol);
   }
   TEST_F(TestDownscalerBilinear, trapezoidHorizontalParallel) {
      float x0 = 0;
      float x1 = 0.1;
      float x2 = 1;
      float x3 = 0.9;
      float y0 = 0;
      float y1 = 1;
      float y2 = 0;
      float y3 = 1;
      float v0 = 0;
      float v1 = 1;
      float v2 = 2;
      float v3 = 3;
      EXPECT_FLOAT_EQ(v0, DownscalerBilinear::bilinear(x0,y0, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v1, DownscalerBilinear::bilinear(x1,y1, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_FLOAT_EQ(v2, DownscalerBilinear::bilinear(x2,y2, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      float tol = 0.03;
      EXPECT_FLOAT_EQ(v3, DownscalerBilinear::bilinear(x3,y3, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3));
      EXPECT_NEAR(1.5, DownscalerBilinear::bilinear(0.5,0.5, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3), tol);
      EXPECT_NEAR(1.1945165, DownscalerBilinear::bilinear(0.211,0.903, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3), tol);
   }
   // // Now reverse the points, so that v2,v3 are to the left of v0,v1
   // EXPECT_FLOAT_EQ(10.146060606061, DownscalerBilinear::bilinear(1.3,0.6,x1,x1,x0,x0,y0,y1,y0,y1,v2,v3,v0,v1));
   // // Put v1,v3 below v0,v1
   // EXPECT_FLOAT_EQ(10.146060606061, DownscalerBilinear::bilinear(1.3,0.6,x0,x0,x1,x1,y1,y0,y1,y0,v1,v0,v3,v2));
   TEST_F(TestDownscalerBilinear, downscale) {
      // Same as above, but using an input file
      DownscalerBilinear d(mVariable, mVariable, Options());
      FileFake from(Options("nLat=2 nLon=2 nEns=1 nTime=1"));
      FileFake to(Options("nLat=1 nLon=1 nEns=1 nTime=1"));
      setLatLon(from, (const float[]) {0.3, 4.8}, (const float[]){1.1, 2.2});
      setLatLon(to,   (const float[]) {0.6},    (const float[]){1.3});

      Field& fromT = *from.getField(mVariable, 0);
      fromT(0,0,0) = 10;
      fromT(1,0,0) = 10.3;
      fromT(0,1,0) = 10.7;
      fromT(1,1,0) = 10.9;
      d.downscale(from, to);
      const Field& toT   = *to.getField(mVariable, 0);
      EXPECT_FLOAT_EQ(10.14606060606061, toT(0,0,0));
   }
   TEST_F(TestDownscalerBilinear, outside) {
      DownscalerBilinear d(mVariable, mVariable, Options());
      vec2 lats = make2x2(0, 0, 1, 1);
      vec2 lons = make2x2(0, 1, 0, 1);

      int I1, I2, J1, J2;
      int I = 1;
      int J = 1;

      EXPECT_FALSE(d.findCoords(-0.1, -0.1, lats, lons, 0, 0, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords(-0.1 , 0.5, lats, lons, 0, 1, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords(-0.1 , 1.1, lats, lons, 0, 1, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 1.1, -0.1, lats, lons, 1, 0, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 1.1,  0.5, lats, lons, 1, 1, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 1.1,  1.1, lats, lons, 1, 1, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords(-0.1, -0.1, lats, lons, 0, 0, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 0.5, -0.1, lats, lons, 1, 0, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 1.1, -0.1, lats, lons, 1, 0, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords(-0.1,  1.1, lats, lons, 0, 1, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 0.5,  1.1, lats, lons, 1, 1, I1, J1, I2, J2));
      EXPECT_FALSE(d.findCoords( 1.1,  1.1, lats, lons, 1, 1, I1, J1, I2, J2));

      ASSERT_TRUE(d.findCoords(0.9, 0.99, lats, lons, I, J, I1, J1, I2, J2));
      EXPECT_EQ(0, I1);
      EXPECT_EQ(1, I2);
      EXPECT_EQ(0, J1);
      EXPECT_EQ(1, J2);

      // Lats go backwards
      lats = make2x2(1, 1, 0, 0);
      lons = make2x2(0, 1, 0, 1);
      I = 0;
      EXPECT_FALSE(d.findCoords(1.1, 1, lats, lons, I, J, I1, J1, I2, J2));
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
