#include "gridpp.h"
#include <math.h>
#include <algorithm>
#include <armadillo>
#include <assert.h>
#include <exception>

using namespace gridpp;

namespace {
    typedef arma::mat mattype;
    typedef arma::vec vectype;
    typedef arma::cx_mat cxtype;

    void check_vec(vec2 input, int Y, int X);
    void check_vec(vec input, int S);

    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
}

vec2 gridpp::optimal_interpolation(const gridpp::Grid& bgrid,
        const vec2& background,
        const gridpp::Points& points,
        const vec& pobs,
        const vec& pratios,
        const vec& pbackground,
        const gridpp::StructureFunction& structure,
        int max_points) {
    double s_time = gridpp::clock();

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");

    if(bgrid.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background grid and observations points must be of same coordinate type (lat/lon or x/y)");
    }
    if(background.size() != bgrid.size()[0] || background[0].size() != bgrid.size()[1]) {
        std::stringstream ss;
        ss << "input field (" << bgrid.size()[0] << "," << bgrid.size()[1] << ") is not the same size as the grid (" << background.size() << "," << background[0].size() << ")";
        throw std::invalid_argument(ss.str());
    }
    if(pobs.size() != points.size()) {
        std::stringstream ss;
        ss << "Observations (" << pobs.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(pratios.size() != points.size()) {
        std::stringstream ss;
        ss << "Ratios (" << pratios.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(pbackground.size() != points.size()) {
        std::stringstream ss;
        ss << "Background (" << pbackground.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }

    int nY = background.size();
    int nX = background[0].size();

    gridpp::Points bpoints = bgrid.to_points();
    vec background1(nY * nX);
    int count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            background1[count] = background[y][x];
            count++;
        }
    }
    vec output1 = optimal_interpolation(bpoints, background1, points, pobs, pratios, pbackground, structure, max_points);
    vec2 output = gridpp::init_vec2(nY, nX);
    count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            output[y][x] = output1[count];
            count++;
        }
    }
    return output;
}

vec gridpp::optimal_interpolation(const gridpp::Points& bpoints,
        const vec& background,
        const gridpp::Points& points,
        const vec& pobs,  // gObs
        const vec& pratios,   // gCi
        const vec& pbackground,
        const gridpp::StructureFunction& structure,
        int max_points) {
    double s_time = gridpp::clock();

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");

    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background points and observations points must be of same coordinate type (lat/lon or x/y)");
    }
    if(background.size() != bpoints.size()) {
        std::stringstream ss;
        ss << "Input field (" << bpoints.size() << ") is not the same size as the grid (" << background.size() << ")";
        throw std::invalid_argument(ss.str());
    }
    if(pobs.size() != points.size()) {
        std::stringstream ss;
        ss << "Observations (" << pobs.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(pratios.size() != points.size()) {
        std::stringstream ss;
        ss << "Ratios (" << pratios.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }

    int nY = background.size();

    // Prepare output matrix
    vec output(nY);

    vec blats = bpoints.get_lats();
    vec blons = bpoints.get_lons();
    vec belevs = bpoints.get_elevs();
    vec blafs = bpoints.get_lafs();

    vec plats = points.get_lats();
    vec plons = points.get_lons();
    vec pelevs = points.get_elevs();
    vec plafs = points.get_lafs();
    int nS = plats.size();

    float localizationRadius = structure.localization_distance();

    // Compute the background value at observation points (Y)
    vec gY = pbackground;

    double curr_time = gridpp::clock();
    CoordinateType coordinate_type = points.get_coordinate_type();
    #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];

        // Find observations within localization radius
        // TODO: Check that the chosen ones have elevation
        ivec lLocIndices0 = points.get_neighbours(lat, lon, localizationRadius);
        if(lLocIndices0.size() == 0) {
            // If we have too few observations though, then use the background
            output[y] = background[y];
            continue;
        }
        std::vector<int> lLocIndices;
        lLocIndices.reserve(lLocIndices0.size());
        std::vector<std::pair<float,int> > lRhos0;

        // Calculate gridpoint to observation rhos
        lRhos0.reserve(lLocIndices0.size());
        for(int i = 0; i < lLocIndices0.size(); i++) {
            int index = lLocIndices0[i];
            if(gridpp::is_valid(pobs[index])) {
                Point p1(plats[index], plons[index], pelevs[index], plafs[index], coordinate_type);
                Point p2(lat, lon, elev, laf, coordinate_type);
                float rho = structure.corr_background(p1, p2);
                if(rho > 0) {
                    lRhos0.push_back(std::pair<float,int>(rho, i));
                }
            }
        }

        // Make sure we don't use too many observations
        arma::vec lRhos;
        if(max_points > 0 && lRhos0.size() > max_points) {
            // If we have too many locations, then only keep the best ones based on rho.
            // Otherwise, just use the last locations added
            lRhos = arma::vec(max_points);
            std::sort(lRhos0.begin(), lRhos0.end(), ::sort_pair_first<float,int>());
            for(int i = 0; i < max_points; i++) {
                // The best values start at the end of the array
                int index = lRhos0[lRhos0.size() - 1 - i].second;
                lLocIndices.push_back(lLocIndices0[index]);
                lRhos(i) = lRhos0[lRhos0.size() - 1 - i].first;
            }
        }
        else {
            lRhos = arma::vec(lRhos0.size());
            for(int i = 0; i < lRhos0.size(); i++) {
                int index = lRhos0[i].second;
                lLocIndices.push_back(lLocIndices0[index]);
                lRhos(i) = lRhos0[i].first;
            }
        }

        int lS = lLocIndices.size();
        if(lS == 0) {
            // If we have too few observations though, then use the background
            output[y] = background[y];
            continue;
        }

        vectype lObs(lS);
        // Compute Y (model at obs-locations)
        vectype lY(lS);
        // Current grid-point to station error covariance matrix
        mattype lG(1, lS, arma::fill::zeros);
        // Station to station error covariance matrix
        mattype lP(lS, lS, arma::fill::zeros);
        // Station variance
        mattype lR(lS, lS, arma::fill::zeros);
        for(int i = 0; i < lS; i++) {
            int index = lLocIndices[i];
            lObs(i) = pobs[index];
            lY(i) = gY[index];
            lR(i, i) = pratios[index];
            lG(0, i) = lRhos(i);
            Point p1(plats[index], plons[index], pelevs[index], plafs[index]);
            for(int j = 0; j < lS; j++) {
                int index_j = lLocIndices[j];
                Point p2(plats[index_j], plons[index_j], pelevs[index_j], plafs[index_j]);
                lP(i, j) = structure.corr(p1, p2);
            }
        }
        mattype lGSR = lG * arma::inv(lP + lR);
        vectype dx = lGSR * (lObs - lY);
        output[y] = background[y] + dx[0];
    }

    return output;
}

vec gridpp::optimal_interpolation_transform(const gridpp::Points& bpoints,
        const vec& background,
        float bsigma,
        const gridpp::Points& points,
        const vec& pobs,
        const vec& psigmas,
        const vec& pbackground,
        const gridpp::StructureFunction& structure,
        int max_points,
        const gridpp::Transform& transform) {

    // Check input data
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");

    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background points and observations points must be of same coordinate type (lat/lon or x/y)");
    }
    if(background.size() != bpoints.size()) {
        std::stringstream ss;
        ss << "Input field (" << bpoints.size() << ") is not the same size as the grid (" << background.size() << ")";
        throw std::invalid_argument(ss.str());
    }
    if(pobs.size() != points.size()) {
        std::stringstream ss;
        ss << "Observations (" << pobs.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }
    if(psigmas.size() != points.size()) {
        std::stringstream ss;
        ss << "Sigmas (" << psigmas.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }

    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coordinate type (lat/lon or x/y)");
    }
    int nY = background.size();
    int nS = points.size();

    // Transform the background
    vec background_transformed = background;
    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float value = background_transformed[y];
        if(gridpp::is_valid(value))
            background_transformed[y] = transform.forward(value);
    }

    vec pbackground_transformed = pbackground;
    vec pobs_transformed = pobs;
    vec pratios(nS);
    for(int s = 0; s < nS; s++) {
        float value = pbackground_transformed[s];
        if(gridpp::is_valid(value))
            pbackground_transformed[s] = transform.forward(value);

        value = pobs_transformed[s];
        if(gridpp::is_valid(value))
            pobs_transformed[s] = transform.forward(value);

        pratios[s] = psigmas[s] * psigmas[s] / bsigma / bsigma;
    }

    vec analysis = optimal_interpolation(bpoints, background_transformed, points, pobs_transformed, pratios, pbackground_transformed, structure, max_points);

    // Transform the background
    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float value = analysis[y];
        if(gridpp::is_valid(value))
            analysis[y] = transform.backward(value);
    }

    return analysis;
}
vec2 gridpp::optimal_interpolation_transform(const gridpp::Grid& bgrid,
        const vec2& background,
        float bsigma,
        const gridpp::Points& points,
        const vec& pobs,
        const vec& psigmas,
        const vec& pbackground,
        const gridpp::StructureFunction& structure,
        int max_points,
        const gridpp::Transform& transform) {

    if(bgrid.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background grid and observations points must be of same coordinate type (lat/lon or x/y)");
    }
    int nY = background.size();
    int nX = background[0].size();
    int nS = points.size();

    // Transform the background
    vec2 background_transformed = background;
    // #pragma omp parallel for
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            float value = background_transformed[y][x];
            if(gridpp::is_valid(value))
                background_transformed[y][x] = transform.forward(value);
        }
    }

    vec pbackground_transformed = pbackground;
    vec pobs_transformed = pobs;
    vec pratios(nS);
    for(int s = 0; s < nS; s++) {
        float value = pbackground_transformed[s];
        if(gridpp::is_valid(value))
            pbackground_transformed[s] = transform.forward(value);

        value = pobs_transformed[s];
        if(gridpp::is_valid(value))
            pobs_transformed[s] = transform.forward(value);

        pratios[s] = psigmas[s] * psigmas[s] / bsigma / bsigma;
    }

    vec2 analysis = optimal_interpolation(bgrid, background_transformed, points, pobs_transformed, pratios, pbackground_transformed, structure, max_points);

    // Transform the background
    // #pragma omp parallel for
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            float value = analysis[y][x];
            if(gridpp::is_valid(value))
                analysis[y][x] = transform.backward(value);
        }
    }

    return analysis;
}
