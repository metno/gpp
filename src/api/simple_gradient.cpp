#include "gridpp.h"

using namespace gridpp;

vec gridpp::simple_gradient(const Grid& igrid, const Points& opoints, vec2 ivalues, float elev_gradient) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 ielevs = igrid.get_elevs();
    vec oelevs = opoints.get_elevs();
    vec olats = opoints.get_lats();
    vec olons = opoints.get_lons();

    vec output(opoints.size(), gridpp::MV);

    // #pragma omp parallel for
    for(int i = 0; i < opoints.size(); i++) {
        ivec indices = igrid.get_nearest_neighbour(olats[i], olons[i]);
        float oelev = oelevs[i];
        float ielev = ielevs[indices[0]][indices[1]];
        float elev_diff = oelev - ielev;
        float elev_corr = elev_diff * elev_gradient;
        output[i] = ivalues[indices[0]][indices[1]] + elev_corr;
    }
    return output;
}

vec2 gridpp::simple_gradient(const Grid& igrid, const Grid& ogrid, vec2 ivalues, float elev_gradient) {
    if(gridpp::compatible_size(igrid, ivalues))
        throw std::invalid_argument("Grid size is not the same as values");
    vec2 ielevs = igrid.get_elevs();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();

    int nLat = ogrid.size()[0];
    int nLon = ogrid.size()[1];

    vec2 output;
    output.resize(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    // #pragma omp parallel for
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec indices = igrid.get_nearest_neighbour(olats[i][j], olons[i][j]);
            float oelev = oelevs[i][j];
            float ielev = ielevs[indices[0]][indices[1]];
            float elev_diff = oelev - ielev;
            float elev_corr = elev_diff * elev_gradient;
            output[i][j] = ivalues[indices[0]][indices[1]] + elev_corr;
        }
    }
    return output;
}
