#include "gridpp.h"

using namespace gridpp;

float gridpp::pressure(float ielev, float oelev, float ipressure, float itemperature) {
    float g0 = 9.80665;
    float M = 0.0289644;
    float R = 8.3144598;
    float value = gridpp::MV;
    if(gridpp::is_valid(ielev) && gridpp::is_valid(oelev) && gridpp::is_valid(ipressure) && gridpp::is_valid(itemperature))
        value = ipressure * exp(-g0 * M * (oelev-ielev) / (R * itemperature));
    return value;
}

vec gridpp::pressure(const vec& ielev, const vec& oelev, const vec& ipressure, const vec& itemperature) {
    int N = ielev.size();
    if(oelev.size() != N || ipressure.size() != N || itemperature.size() != N)
        throw std::invalid_argument("pressure: Input arguments must be of the same size");

    vec values(N);
    #pragma omp parallel for
    for(int i = 0; i < N; i++) {
        values[i] = gridpp::pressure(ielev[i], oelev[i], ipressure[i], itemperature[i]);
    }
    return values;
}
