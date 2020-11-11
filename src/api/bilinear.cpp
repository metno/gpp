#include "gridpp.h"
#include <math.h>
#include <iostream>

using namespace gridpp;

namespace {
    // Bilinear interpolation for a given point
    float calc(const Grid& grid, const vec2& ivalues, float lat, float lon);
    // Bilinear interpolation based on 4 surrounding points with coordinates (x0,y0), (x1,y1), etc
    // and values v0, v1, etc
    float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3);
    // Compute s,t when the points form a parallelogram
    bool calcParallelogram(float x, float y, float X1, float X2, float X3, float X4, float Y1, float Y2, float Y3, float Y4, float &t, float &s);
    // Compute s,t when the points do not form a parallelogram
    bool calcGeneral(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float &t, float &s);
}
vec2 gridpp::bilinear(const Grid& igrid, const Grid& ogrid, vec2 ivalues) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();
    vec2 iInputLats = igrid.get_lats();
    vec2 iInputLons = igrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output;
    output.resize(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    // Algorithm from here:
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            // Use the four points surrounding the lookup point. Use the nearest neighbour
            // and then figure out what side of this point the other there points are.
            float lat = iOutputLats[i][j];
            float lon = iOutputLons[i][j];
            output[i][j] = ::calc(igrid, ivalues, lat, lon);
        }
    }
    return output;
}

vec gridpp::bilinear(const Grid& igrid, const Points& opoints, vec2 ivalues) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();
    vec2 iInputLats = igrid.get_lats();
    vec2 iInputLons = igrid.get_lons();

    int nPoints = iOutputLats.size();

    vec output(nPoints, gridpp::MV);

    // Algorithm from here:
    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        // Use the four points surrounding the lookup point. Use the nearest neighbour
        // and then figure out what side of this point the other there points are.
        float lat = iOutputLats[i];
        float lon = iOutputLons[i];
        output[i] = ::calc(igrid, ivalues, lat, lon);
    }
    return output;
}

namespace {
    bool calcParallelogram(float x, float y, float X1, float X2, float X3, float X4, float Y1, float Y2, float Y3, float Y4, float &t, float &s) {
        // std::cout << "Method 3: Parallelogram" << std::endl;
        float X31 = X3 - X1;
        float X21 = X2 - X1;
        float Y21 = Y2 - Y1;
        float Y31 = Y3 - Y1;

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

    bool calcGeneral(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float &t, float &s) {
        // std::cout << "Method general" << std::endl;
        float a = -x0 + x2;
        float b = -x0 + x1;
        float c = x0 - x1 - x2 + x3;
        float d = x - x0;
        float e = -y0 + y2;
        float f = -y0 + y1;
        float g = y0 - y1 - y2 + y3;
        float h = y - y0;
        float alpha=gridpp::MV, beta=gridpp::MV;
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
    float bilinear(float x, float y, float x0, float x1, float x2, float x3, float y0, float y1, float y2, float y3, float v0, float v1, float v2, float v3) {
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

       float s = gridpp::MV, t = gridpp::MV;
       bool rectangularGrid = (X1 == X3 && X2 == X4 && Y1 == Y2 && Y3 == Y4);
       // TODO: Why are the tolerances so high?
       bool verticalParallel = fabs((X3 - X1)*(Y4 - Y2) - (X4 - X2)*(Y3 - Y1)) <= 1e-4;
       bool horizontalParallel = fabs((X2 - X1)*(Y4 - Y3) - (X4 - X3)*(Y2 - Y1)) <= 1e-4;

       // std::cout << rectangularGrid << " " << fabs((X3 - X1)*(Y4 - Y2) - (X4 - X2)*(Y3 - Y1)) << " " << fabs((X2 - X1)*(Y4 - Y3) - (X4 - X3)*(Y2 - Y1)) << std::endl;

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
          throw std::runtime_error(ss.str());
       }
       assert(s >= 0 && s <= 1 && t >= 0 && t <= 1);
       float value = P1 * (1 - s) * ( 1 - t) + P2 * s * (1 - t) + P3 * (1 - s) * t + P4 * s * t;

       return value;
    }

    float calc(const Grid& grid, const vec2& ivalues, float lat, float lon) {
        int I1 = gridpp::MV;
        int I2 = gridpp::MV;
        int J1 = gridpp::MV;
        int J2 = gridpp::MV;
        const vec2 iInputLats = grid.get_lats();
        const vec2 iInputLons = grid.get_lons();
        float output = gridpp::MV;
        bool inside = grid.get_box(lat, lon, I1, J1, I2, J2);
        // std::cout << "Coords: " << I1 << " " << J1 << " " << I2 << " " << J2 << " " << inside << std::endl;
        if(inside) {
            float x0 = iInputLons[I1][J1];
            float x1 = iInputLons[I2][J1];
            float x2 = iInputLons[I1][J2];
            float x3 = iInputLons[I2][J2];
            float y0 = iInputLats[I1][J1];
            float y1 = iInputLats[I2][J1];
            float y2 = iInputLats[I1][J2];
            float y3 = iInputLats[I2][J2];

            float v0 = ivalues[I1][J1];
            float v1 = ivalues[I2][J1];
            float v2 = ivalues[I1][J2];
            float v3 = ivalues[I2][J2];
            if(gridpp::is_valid(v0) && gridpp::is_valid(v1) && gridpp::is_valid(v2) && gridpp::is_valid(v3)) {
                float value = ::bilinear(lon, lat, x0, x1, x2, x3, y0, y1, y2, y3, v0, v1, v2, v3);
                output = value;
            }
            else {
                ivec nn = grid.get_nearest_neighbour(lat, lon);
                output = ivalues[nn[0]][nn[1]];
            }
        }
        else {
            // The point is outside the input domain. Revert to nearest neighbour
            // Util::warning("Point is outside domain, cannot bilinearly interpolate");
            ivec nn = grid.get_nearest_neighbour(lat, lon);
            output = ivalues[nn[0]][nn[1]];
        }
        return output;
    }
}
