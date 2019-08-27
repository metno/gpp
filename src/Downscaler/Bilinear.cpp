#include "Bilinear.h"
#include "../File/File.h"
#include "../Util.h"
#include <math.h>

DownscalerBilinear::DownscalerBilinear(const Variable& iInputVariable, const Variable& iOutputVariable, const Options& iOptions) :
      Downscaler(iInputVariable, iOutputVariable, iOptions) {
   iOptions.check();
}

void DownscalerBilinear::downscaleCore(const File& iInput, File& iOutput) const {
   int nTime = iInput.getNumTime();
   vec2 ilats = iInput.getLats();
   vec2 ilons = iInput.getLons();
   vec2 olats = iOutput.getLats();
   vec2 olons = iOutput.getLons();

   // Get nearest neighbour
   vec2Int nearestI, nearestJ;
   Downscaler::getNearestNeighbour(iInput, iOutput, nearestI, nearestJ);

   for(int t = 0; t < nTime; t++) {
      Field& ifield = *iInput.getField(mInputVariable, t);
      Field& ofield = *iOutput.getField(mOutputVariable, t, true);
      downscaleField(ifield, ofield, ilats, ilons, olats, olons, nearestI, nearestJ);
   }
}

float DownscalerBilinear::bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3) {
   float Y1 = y1;
   float Y2 = y3;
   float Y3 = y0;
   float Y4 = y2;
   float X1 = x1;
   float X2 = x3;
   float X3 = x0;
   float X4 = x2;
   float P1 = v1;
   float P2 = v3;
   float P3 = v0;
   float P4 = v2;

   // General method based on: https://stackoverflow.com/questions/23920976/bilinear-interpolation-with-non-aligned-input-points
   // Parallelogram method based on: http://www.ahinson.com/algorithms_general/Sections/InterpolationRegression/InterpolationIrregularBilinear.pdf

   float s = Util::MV, t = Util::MV;
   bool rectangularGrid = (X1 == X3 && X2 == X4 && Y1 == Y2 && Y3 == Y4);
   bool verticalParallel = fabs((X3 - X1)*(Y4 - Y2) - (X4 - X2)*(Y3 - Y1)) <= 1e-6;
   bool horizontalParallel = fabs((X2 - X1)*(Y4 - Y3) - (X4 - X3)*(Y2 - Y1)) <= 1e-6;

   // Compute s and t
   if(verticalParallel && horizontalParallel)
      calcParallelogram(x, y, X1, X2, X3, X4, Y1, Y2, Y3, Y4, t, s);
   else
      calcGeneral(x, y, x0, x1, x2, x3, y0, y1, y2, y3, t, s);

   if(t >= 1 && t <= 1.15)
      t = 1;
   if(t <= 0 && t >= -0.15)
      t = 0;
   if(s >= 1 && s <= 1.15)
      s = 1;
   if(s <= 0 && s >= -0.15)
      s = 0;
   if(!(s >= 0 && s <= 1 && t >= 0 && t <= 1)) {
      std::stringstream ss;
      ss << "Problem with bilinear interpolation. Grid is rotated/distorted in a way that is not supported. s=" << s << " and t=" << t << " are outside [-0.05,1.05].";
      Util::error(ss.str());
   }
   assert(s >= 0 && s <= 1 && t >= 0 && t <= 1);
   float value = P1 * (1 - s) * ( 1 - t) + P2 * s * (1 - t) + P3 * (1 - s) * t + P4 * s * t;

   return value;
}

vec2 DownscalerBilinear::downscaleVec(const vec2& iInput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ) {

   int nLat = iOutputLats.size();
   int nLon = iOutputLats[0].size();

   vec2 output;
   output.resize(nLat);
   for(int i = 0; i < nLat; i++)
      output[i].resize(nLon);

   // Algorithm from here:
   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         // Use the four points surrounding the lookup point. Use the nearest neighbour
         // and then figure out what side of this point the other there points are.
         int I = nearestI[i][j];
         int J = nearestJ[i][j];
         float lat = iOutputLats[i][j];
         float lon = iOutputLons[i][j];
         int I1 = Util::MV;
         int I2 = Util::MV;
         int J1 = Util::MV;
         int J2 = Util::MV;
         bool inside = findCoords(lat, lon, iInputLats, iInputLons, I, J, I1, J1, I2, J2);
         if(inside) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            float v0 = iInput[I1][J1];
            float v1 = iInput[I2][J1];
            float v2 = iInput[I1][J2];
            float v3 = iInput[I2][J2];
            if(Util::isValid(v0) && Util::isValid(v1) &&Util::isValid(v2) &&Util::isValid(v3)) {
               float value = bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
               output[i][j] = value;
            }
            else
               output[i][j] = iInput[I][J];
         }
         else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            output[i][j] = iInput[I][J];
         }
      }
   }
   return output;
}

void DownscalerBilinear::downscaleField(const Field& iInput, Field& iOutput,
            const vec2& iInputLats, const vec2& iInputLons,
            const vec2& iOutputLats, const vec2& iOutputLons,
            const vec2Int& nearestI, const vec2Int& nearestJ) {

   int nLat = iOutput.getNumY();
   int nLon = iOutput.getNumX();
   int nEns = iOutput.getNumEns();

   #pragma omp parallel for
   for(int i = 0; i < nLat; i++) {
      for(int j = 0; j < nLon; j++) {
         int I = nearestI[i][j];
         int J = nearestJ[i][j];
         float lat = iOutputLats[i][j];
         float lon = iOutputLons[i][j];
         // Find out if current lookup point is right/above the nearest neighbour
         int I1 = Util::MV;
         int I2 = Util::MV;
         int J1 = Util::MV;
         int J2 = Util::MV;
         bool inside = findCoords(lat, lon, iInputLats, iInputLons, I, J, I1, J1, I2, J2);
         if(inside) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            for(int ee = 0; ee < nEns; ee++) {
               float v0 = iInput(I1, J1, ee);
               float v1 = iInput(I2, J1, ee);
               float v2 = iInput(I1, J2, ee);
               float v3 = iInput(I2, J2, ee);
               if(Util::isValid(v0) && Util::isValid(v1) &&Util::isValid(v2) &&Util::isValid(v3)) {
                  float value = bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
                  iOutput(i,j,ee) = value;
               }
               else
                  iOutput(i,j,ee) = Util::MV;
            }
         }
         else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            for(int ee = 0; ee < nEns; ee++) {
               iOutput(i,j,ee) = iInput(I, J, ee);
            }
         }
      }
   }
}

bool DownscalerBilinear::findCoords(float iLat, float iLon, const vec2& iLats, const vec2& iLons, int I, int J, int& I1, int& J1, int& I2, int& J2) {
   // Use the four points surrounding the lookup point. Use the nearest neighbour
   // and then figure out what side of this point the other there points are.
   // Find out if current lookup point is right/above the nearest neighbour
   bool isRight = iLon > iLons[I][J];
   bool isAbove = iLat > iLats[I][J];
   int nI = iLats.size();
   int nJ = iLats[0].size();

   // Check if we are outside the grid
   // 0 1 2 3 x
   if(isAbove && I == nI-1 && iLats[I-1][J] < iLats[I][J])
      return false;
   // x 0 1 2 3
   if(!isAbove && I == 0 && iLats[0][J] < iLats[1][J])
      return false;
   // 3 2 1 0 x
   if(isAbove && I == 0 && iLats[0][J] > iLats[1][J])
      return false;
   // x 3 2 1 0
   if(!isAbove && I == nI-1 && iLats[I-1][J] > iLats[I][J])
      return false;

   if(isRight && J == nJ-1 && iLons[I][J-1] < iLons[I][J])
      return false;
   // x 0 1 2 3
   if(!isRight && J == 0 && iLons[I][0] < iLons[I][1])
      return false;
   // 3 2 1 0 x
   if(isRight && J == 0 && iLons[I][0] > iLons[I][1])
      return false;
   // x 3 2 1 0
   if(!isRight && I == nJ-1 && iLons[I][J-1] > iLons[I][J])
      return false;

   // Are the lats/lons increating when the index increases?
   bool Iinc = (I == 0) || (iLats[I][J] > iLats[I-1][J]);
   bool Jinc = (J == 0) || (iLons[I][J] > iLons[I][J-1]);

   // Find out which Is to use
   if(isAbove) {
      if(Iinc) {
         I1 = I;
         I2 = I + 1;
      }
      else {
         I1 = I - 1;
         I2 = I;
      }
   }
   else {
      if(Iinc) {
         I1 = I - 1;
         I2 = I;
      }
      else {
         I1 = I;
         I2 = I + 1;
      }
   }

   // Find out which Js to use
   if(isRight) {
      if(Jinc) {
         J1 = J;
         J2 = J + 1;
      }
      else {
         J1 = J - 1;
         J2 = J;
      }
   }
   else {
      if(Jinc) {
         J1 = J - 1;
         J2 = J;
      }
      else {
         J1 = J;
         J2 = J + 1;
      }
   }
   if((I1 < 0 || I1 >= iLats.size()) ||
      (I2 < 0 || I2 >= iLats.size()) ||
      (J1 < 0 || J1 >= iLats[0].size()) ||
      (J2 < 0 || J2 >= iLats[0].size())) {
      return false;
   }
   return true;
}
bool DownscalerBilinear::calcParallelogram(float x, float y, float X1, float X2, float X3, float X4, float Y1, float Y2, float Y3, float Y4, float &t, float &s) {
   // std::cout << "Method 3: Parallelogram" << std::endl;
   float X31 = X3 - X1;
   float X21 = X2 - X1;
   float Y42 = Y4 - Y2;
   float Y21 = Y2 - Y1;
   float Y31 = Y3 - Y1;
   float Y43 = Y4 - Y3;
   float X42 = X4 - X2;
   float X43 = X4 - X3;

   float A = X21;
   float B = X31;
   float C = Y21;
   float D = Y31;
   assert(A * D != B * C);
   float det = 1 / (A * D - B * C);
   s = det * ((x - X1) * (D) + (y - Y1) * (-B));
   t = det * ((x - X1) * (-C) + (y - Y1) * (A));
   return true;
}

bool DownscalerBilinear::calcGeneral(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float &t, float &s) {
   // std::cout << "Method general" << std::endl;
   float a = -x0 + x2;
   float b = -x0 + x1;
   float c = x0 - x1 - x2 + x3;
   float d = x - x0;
   float e = -y0 + y2;
   float f = -y0 + y1;
   float g = y0 - y1 - y2 + y3;
   float h = y - y0;
   float alpha=Util::MV, beta=Util::MV;
   float Y1 = y1;
   float Y2 = y3;
   float Y3 = y0;
   float Y4 = y2;
   float X1 = x1;
   float X2 = x3;
   float X3 = x0;
   float X4 = x2;
   float X31 = X3 - X1;
   float X21 = X2 - X1;
   float Y42 = Y4 - Y2;
   float Y21 = Y2 - Y1;
   float Y31 = Y3 - Y1;
   float Y43 = Y4 - Y3;
   float X42 = X4 - X2;
   float X43 = X4 - X3;
   if(2 * c * e - 2 * a * g != 0 && 2 * c * f - 2 * b * g != 0) {
      alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
      beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
   }
   else if(2 * c * f - 2 * b * g == 0) {
      // std::cout << "Need to diagnose t" << std::endl;
      alpha = -(b * e - a * f + d * g - c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h, 2)))/(2 * c * e - 2 * a * g);
      float s = alpha;
      float t;
      if(Y3 + Y43 * s - Y1 - Y21 * s == 0) {
         assert(X3 + X43 * s - X1 - X21 * s != 0);
         t = (x - X1 - X21 * s) / (X3 + X43 * s - X1 - X21 * s);
      }
      else {
         t = (y - Y1 - Y21 * s) / (Y3 + Y43 * s - Y1 - Y21 * s);
      }
      beta = 1 - t;
   }
   else if(2 * c * e - 2 * a * g == 0) {
      // std::cout << "Need to diagnose s" << std::endl;
      beta  = (b * e - a * f - d * g + c * h + sqrt(-4 * (c * e - a * g) * (d * f - b * h) + pow(b * e - a * f + d * g - c * h,2)))/(2 * c * f - 2 * b * g);
      float t =  1 - beta;
      float s;
      if(Y2 + Y42 * t - Y1 - Y31 * t == 0) {
         assert(X2 + X42 * t - X1 - X31 * t != 0);
         s = (x - X1 - X31 * t) / (X2 + X42 * t - X1 - X31 * t);
      }
      else {
         s = (y - Y1 - Y31 * t) / (Y2 + Y42 * t - Y1 - Y31 * t);
      }
      alpha = s;
   }
   s = alpha;
   t = 1 - beta;
   return true;
}

std::string DownscalerBilinear::description(bool full) {
   std::stringstream ss;
   if(full)
      ss << Util::formatDescription("-d bilinear", "Bilinear interpolation using the 4 points surrounding the lookup point. If the lookup point is outside the input domain, then the nearest neighbour is used.") << std::endl;
   else
      ss << Util::formatDescription("-d bilinear", "Bilinear interpolation") << std::endl;
   return ss.str();
}
