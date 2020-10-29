#include "gridpp.h"
#include <cmath>

using namespace gridpp;

float gridpp::wind_speed(float xwind, float ywind) {
    return sqrt(xwind * xwind + ywind * ywind);
}
vec gridpp::wind_speed(const vec& xwind, const vec& ywind) {
    if(xwind.size() != ywind.size())
        throw std::invalid_argument("xwind and ywind must be of the same size");
    int N = xwind.size();
    vec values(N, gridpp::MV);
    for(int n = 0; n < N; n++) {
        values[n] = wind_speed(xwind[n], ywind[n]);
    }
    return values;
}
float gridpp::wind_direction(float xwind, float ywind) {
    float dir = std::atan2(-xwind, -ywind) * 180 / gridpp::pi;
    if(dir < 0)
        dir += 360;

    return dir;
}
vec gridpp::wind_direction(const vec& xwind, const vec& ywind) {
    if(xwind.size() != ywind.size())
        throw std::invalid_argument("xwind and ywind must be of the same size");
    int N = xwind.size();
    vec values(N, gridpp::MV);
    for(int n = 0; n < N; n++) {
        values[n] = wind_direction(xwind[n], ywind[n]);
    }
    return values;
}
