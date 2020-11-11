#include "gridpp.h"
#include <iostream>
#include <math.h>

using namespace gridpp;

vec2 gridpp::nearest(const Grid& igrid, const Grid& ogrid, vec2 ivalues) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            int I = indices[0];
            int J = indices[1];
            output[i][j] = ivalues[I][J];
        }
    }
    return output;
}
vec3 gridpp::nearest(const Grid& igrid, const Grid& ogrid, vec3 ivalues) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nTime = ivalues.size();
    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec3 output(nTime);
    for(int t = 0; t < nTime; t++) {
        output[t].resize(nLat);
        for(int i = 0; i < nLat; i++)
            output[t][i].resize(nLon);
    }

    #pragma omp parallel for collapse(2)
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(iOutputLats[i][j], iOutputLons[i][j]);
            int I = indices[0];
            int J = indices[1];
            for(int t = 0; t < nTime; t++) {
                output[t][i][j] = ivalues[t][I][J];
            }
        }
    }
    return output;
}
vec gridpp::nearest(const Grid& igrid, const Points& opoints, vec2 ivalues) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();

    vec output(nPoints);

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        ivec indices = igrid.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        int I = indices[0];
        int J = indices[1];
        output[i] = ivalues[I][J];
    }
    return output;
}
vec2 gridpp::nearest(const Grid& igrid, const Points& opoints, vec3 ivalues) {
    if(gridpp::compatible_size(igrid, ivalues))
       throw std::invalid_argument("Grid size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();
    int nTime = ivalues.size();

    vec2 output(nTime);
    for(int t = 0; t < nTime; t++) {
        output[t].resize(nPoints, gridpp::MV);
    }

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        ivec indices = igrid.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        int I = indices[0];
        int J = indices[1];
        for(int t = 0; t < nTime; t++) {
            output[t][i] = ivalues[t][I][J];
        }
    }
    return output;
}
vec gridpp::nearest(const Points& ipoints, const Points& opoints, vec ivalues) {
    if(ipoints.size() != ivalues.size())
        throw std::invalid_argument("Points size is not the same as values");
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();

    vec output(nPoints);

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        int index = ipoints.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        output[i] = ivalues[index];
    }
    return output;
}
vec2 gridpp::nearest(const Points& ipoints, const Points& opoints, vec2 ivalues) {
    if(gridpp::compatible_size(ipoints, ivalues)) {
        std::stringstream ss;
        ss << "Input points (" << ipoints.size() << ") and values (" << ivalues.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    vec iOutputLats = opoints.get_lats();
    vec iOutputLons = opoints.get_lons();

    int nPoints = iOutputLats.size();
    int nTime = ivalues.size();

    vec2 output(nTime);
    for(int t = 0; t < nTime; t++) {
        output[t].resize(nPoints, gridpp::MV);
    }

    #pragma omp parallel for
    for(int i = 0; i < nPoints; i++) {
        int index = ipoints.get_nearest_neighbour(iOutputLats[i], iOutputLons[i]);
        for(int t = 0; t < nTime; t++) {
            output[t][i] = ivalues[t][index];
        }
    }
    return output;
}
