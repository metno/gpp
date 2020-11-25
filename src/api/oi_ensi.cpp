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

template<class Matrix>
void print_matrix(Matrix matrix) {
       matrix.print(std::cout);
}

template void print_matrix< ::mattype>(::mattype matrix);
template void print_matrix< ::cxtype>(::cxtype matrix);

vec3 gridpp::optimal_interpolation_ensi(const gridpp::Grid& bgrid,
        const vec3& background,
        const gridpp::Points& points,
        const vec& pobs,
        const vec& psigmas,
        const vec2& pbackground,
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
    if(psigmas.size() != points.size()) {
        std::stringstream ss;
        ss << "Sigmas (" << psigmas.size() << ") and points (" << points.size() << ") size mismatch";
        throw std::invalid_argument(ss.str());
    }

    int nY = background.size();
    int nX = background[0].size();
    int nE = background[0][0].size();

    gridpp::Points bpoints = bgrid.to_points();
    vec2 background1 = gridpp::init_vec2(nY * nX, nE);
    int count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nE; e++) {
                background1[count][e] = background[y][x][e];
            }
            count++;
        }
    }
    vec2 output1 = optimal_interpolation_ensi(bpoints, background1, points, pobs, psigmas, pbackground, structure, max_points);
    vec3 output = gridpp::init_vec3(nY, nX, nE);
    count = 0;
    for(int y = 0; y < nY; y++) {
        for(int x = 0; x < nX; x++) {
            for(int e = 0; e < nE; e++) {
                output[y][x][e] = output1[count][e];
            }
            count++;
        }
    }
    return output;
}
vec2 gridpp::optimal_interpolation_ensi(const gridpp::Points& bpoints,
        const vec2& background,
        const gridpp::Points& points,
        const vec& pobs,  // gObs
        const vec& psigmas,   // pci
        const vec2& pbackground,
        const gridpp::StructureFunction& structure,
        int max_points) {
    if(max_points < 0)
        throw std::invalid_argument("max_points must be >= 0");
    if(bpoints.get_coordinate_type() != points.get_coordinate_type()) {
        throw std::invalid_argument("Both background and observations points must be of same coorindate type (lat/lon or x/y)");
    }
    if(background.size() != bpoints.size())
        throw std::invalid_argument("Input field is not the same size as the grid");
    if(pobs.size() != points.size())
        throw std::invalid_argument("Observations and points exception mismatch");
    if(psigmas.size() != points.size())
        throw std::invalid_argument("Sigmas and points size mismatch");
    if(pbackground.size() != points.size())
        throw std::invalid_argument("Background and points size mismatch");

    double s_time = gridpp::clock();
    int mX = -1;
    int mY = -1;
    int mMinValidEns = 5;
    int numParameters = 2;
    float sigmac = 0.5;
    float delta = 1;
    bool mExtrapolate = false;
    bool mDiagnose = false;

    int nY = background.size();
    int nEns = background[0].size();

    // Prepare output matrix
    vec2 output;
    output.resize(nY);
    for(int y = 0; y < nY; y++) {
        output[y].resize(nEns);
    }
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

    // Compute Y
    vec2 gY = pbackground;
    vec gYhat(nS);
    for(int i = 0; i < nS; i++) {
        float mean = gridpp::calc_statistic(gY[i], gridpp::Mean);
        for(int e = 0; e < nEns; e++) {
            float value = gY[i][e];
            if(gridpp::is_valid(value) && gridpp::is_valid(mean)) {
                gY[i][e] -= mean;
            }
        }
        gYhat[i] = mean;
    }

    // Calculate number of valid members
    int nValidEns = 0;
    std::vector<int> validEns;
    for(int e = 0; e < nEns; e++) {
        int numInvalid = 0;
        for(int y = 0; y < nY; y++) {
            float value = background[y][e];
            if(!gridpp::is_valid(value))
                numInvalid++;
        }
        if(numInvalid == 0) {
            validEns.push_back(e);
            nValidEns++;
        }
    }
    // TODO: Deal with single member mode
    bool singleMemberMode = nValidEns < mMinValidEns;

    // #pragma omp parallel for
    for(int y = 0; y < nY; y++) {
        float lat = blats[y];
        float lon = blons[y];
        float elev = belevs[y];
        float laf = blafs[y];
        Point p1 = bpoints.get_point(y);

        // Create list of locations for this gridpoint
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
            Point p2 = points.get_point(index);
            float rho = structure.corr_background(p1, p2);
            if(rho > 0) {
                lRhos0.push_back(std::pair<float,int>(rho, i));
            }
        }

        // Make sure we don't use too many observations
        arma::vec lRhos;
        if(max_points > 0 && lRhos0.size() > max_points) {
            // If sorting is enabled and we have too many locations, then only keep the best ones based on rho.
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
        vectype lElevs(lS);
        vectype lLafs(lS);
        for(int i = 0; i < lLocIndices.size(); i++) {
            int index = lLocIndices[i];
            lObs[i] = pobs[index];
            lElevs[i] = pelevs[index];
            lLafs[i] = plafs[index];
        }

        // Compute Y (model at obs-locations)
        mattype lY(lS, nValidEns);
        vectype lYhat(lS);

        for(int i = 0; i < lS; i++) {
            // Use the nearest neighbour for this location
            int index = lLocIndices[i];
            for(int e = 0; e < nValidEns; e++) {
                int ei = validEns[e];
                lY(i, e) = gY[index][ei];
            }
            lYhat[i] = gYhat[index];
        }

        // Compute Rinv
        mattype Rinv(lS, lS, arma::fill::zeros);
        if(numParameters == 2) {
            for(int i = 0; i < lS; i++) {
                int index = lLocIndices[i];
                Rinv(i, i) = lRhos[i] / (psigmas[index] * psigmas[index]);
            }
        }
        else if(numParameters == 3) {
            /*
            // Inverting the matrix is more complicated, since the radar observations
            // have covariances. Therefore invert the covariance matrix for the radar part and
            // insert the values into the big inverse matrix.
            // std::cout << "Computing R matrix" << std::endl;
            // R = get_precipitation_r(gRadarL, pci, lLocIndices, lRhos);
            // Compute little R
            std::vector<int> gRadarIndices;
            gRadarIndices.reserve(lS);
            std::vector<int> lRadarIndices;
            lRadarIndices.reserve(lS);
            for(int i = 0; i < lS; i++) {
                int index = lLocIndices[i];
                if(gRadarL[index] > 0) {
                    gRadarIndices.push_back(index);
                    lRadarIndices.push_back(i);
                }
            }
            int lNumRadar = gRadarIndices.size();

            // Compute R tilde r
            mattype radarR(lNumRadar, lNumRadar, arma::fill::zeros);
            for(int i = 0; i < lNumRadar; i++) {
                for(int j = 0; j < lNumRadar; j++) {
                    int gIndex_i = gRadarIndices[i];
                    int gIndex_j = gRadarIndices[j];
                    int lIndex_i = lRadarIndices[i];
                    int lIndex_j = lRadarIndices[j];
                    if(i == j) {
                        radarR(i, i) = 1;
                    }
                    else {
                        // Equation 5
                        float dist = Util::getDistance(gLocations[gIndex_i].lat(), gLocations[gIndex_i].lon(), gLocations[gIndex_j].lat(), gLocations[gIndex_j].lon(), true);
                        float h = dist / mHLengthC;
                        float rho = (1 + h) * exp(-h);
                        radarR(i, j) = rho;
                    }
                }
            }

            float cond = arma::rcond(radarR);
            if(cond <= 0) {
                std::stringstream ss;
                ss << "Condition number of " << cond << " for radar values. Using raw values";
                gridpp::warning(ss.str());
                for(int e = 0; e < nEns; e++) {
                    (*output)(y, x, e) = (*field)(y, x, e); // Util::MV;
                }
                continue;
            }

            mattype radarRinv(lNumRadar, lNumRadar, arma::fill::zeros);
            radarRinv = arma::inv(radarR);

            for(int i = 0; i < lS; i++) {
                int index = lLocIndices[i];
                Rinv(i, i) = lRhos[i] / (sigma * sigma * pci[index]);
            }
            // Overwrite where we have radar pixels
            for(int i = 0; i < lNumRadar; i++) {
                int ii = lRadarIndices[i];
                for(int j = 0; j < lNumRadar; j++) {
                    int jj = lRadarIndices[j];
                    Rinv(ii, jj) = sqrt(lRhos[ii] * lRhos[jj]) / (sigmac * sigmac) * radarRinv(i, j);
                }
            }
            */
        }
        else {
            abort();
        }

        // Compute C matrix
        // k x nS * nS x nS
        mattype C(nValidEns, lS);
        C = lY.t() * Rinv;

        mattype Pinv(nValidEns, nValidEns);
        float diag = 1 / delta * (nValidEns - 1);

        Pinv = C * lY + diag * arma::eye<mattype>(nValidEns, nValidEns);
        float cond = arma::rcond(Pinv);
        if(cond <= 0) {
            std::stringstream ss;
            ss << "Condition number of " << cond << ". Using raw values";
            gridpp::warning(ss.str());
            for(int e = 0; e < nEns; e++) {
                output[y][e] = background[y][e]; // Util::MV;
            }
            continue;
        }

        // Compute sqrt of matrix. Armadillo 6.6 has this function, but on many systems this
        // is not available. Therefore, compute sqrt using the method found in 6.6
        // cxtype Wcx(nValidEns, nValidEns);
        // status = arma::sqrtmat(Wcx, (nValidEns - 1) * P);
        // mattype W = arma::real(Wcx);

        mattype P = arma::inv(Pinv);
        vectype eigval;
        mattype eigvec;
        bool status = arma::eig_sym(eigval, eigvec, (nValidEns - 1) * P);
        if(!status) {
            std::cout << "Cannot find eigenvector:" << std::endl;
            std::cout << "Lat: " << lat << std::endl;
            std::cout << "Lon: " << lon << std::endl;
            std::cout << "Elev: " << elev << std::endl;
            std::cout << "Laf: " << laf << std::endl;
            std::cout << "Pinv" << std::endl;
            print_matrix<mattype>(Pinv);
            std::cout << "P" << std::endl;
            print_matrix<mattype>(P);
            std::cout << "Y:" << std::endl;
            print_matrix<mattype>(lY);
            std::cout << "lObs:" << std::endl;
            print_matrix<mattype>(lObs);
            std::cout << "Yhat" << std::endl;
            print_matrix<mattype>(lYhat);
        }
        eigval = sqrt(eigval);
        mattype Wcx = eigvec * arma::diagmat(eigval) * eigvec.t();
        mattype W = arma::real(Wcx);

        if(W.n_rows == 0) {
            std::stringstream ss;
            ss << "Could not find the real part of W. Using raw values.";
            gridpp::warning(ss.str());
            for(int e = 0; e < nEns; e++) {
                output[y][e] = background[y][e];
            }
            continue;
        }

        // Compute PC
        mattype PC(nValidEns, lS);
        PC = P * C;

        // Compute w
        vectype w(nValidEns);
        if(mDiagnose)
            w = PC * (arma::ones<vectype>(lS));
        else
            w = PC * (lObs - lYhat);

        // Add w to W
        for(int e = 0; e < nValidEns; e++) {
            for(int e2 = 0; e2 < nValidEns; e2 ++) {
                W(e, e2) = W(e, e2) + w(e) ;
            }
        }

        // Compute X (perturbations about model mean)
        vectype X(nValidEns);
        float total = 0;
        int count = 0;
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            float value = background[y][ei];
            if(gridpp::is_valid(value)) {
                X(e) = value;
                total += value;
                count++;
            }
        }
        float ensMean = total / count;
        for(int e = 0; e < nValidEns; e++) {
            X(e) -= ensMean;
        }

        // Write debugging information
        if(y == mY) {
            std::cout << "Lat: " << lat << std::endl;
            std::cout << "Lon: " << lon << " " << lat << " " << std::endl;
            std::cout << "Elev: " << elev << std::endl;
            std::cout << "Laf: " << laf << std::endl;
            std::cout << "Num obs: " << lS << std::endl;
            std::cout << "Num ens: " << nValidEns << std::endl;
            std::cout << "rhos" << std::endl;
            print_matrix<mattype>(lRhos);
            std::cout << "P" << std::endl;
            print_matrix<mattype>(P);
            std::cout << "C" << std::endl;
            print_matrix<mattype>(C);
            std::cout << "C * lY" << std::endl;
            print_matrix<mattype>(C * lY);
            std::cout << "PC" << std::endl;
            print_matrix<mattype>(PC);
            std::cout << "W" << std::endl;
            print_matrix<mattype>(W);
            std::cout << "w" << std::endl;
            print_matrix<mattype>(w);
            std::cout << "Y:" << std::endl;
            print_matrix<mattype>(lY);
            std::cout << "Yhat" << std::endl;
            print_matrix<mattype>(lYhat);
            std::cout << "lObs" << std::endl;
            print_matrix<mattype>(lObs);
            std::cout << "lObs - Yhat" << std::endl;
            print_matrix<mattype>(lObs - lYhat);
            std::cout << "X" << std::endl;
            print_matrix<mattype>(X);
            std::cout << "elevs" << std::endl;
            print_matrix<mattype>(lElevs);
            std::cout << "lafs" << std::endl;
            print_matrix<mattype>(lLafs);
            std::cout << "Analysis increment:" << std::endl;
            print_matrix<mattype>(X.t() * W);
            std::cout << "My: " << arma::mean(arma::dot(lObs - lYhat, lRhos) / lS) << std::endl;
        }

        // Compute analysis
        for(int e = 0; e < nValidEns; e++) {
            int ei = validEns[e];
            float total = 0;
            for(int k = 0; k < nValidEns; k++) {
                total += X(k) * W(k, e);
            }

            float currIncrement = total;

            float raw = ensMean;

            ///////////////////////////////
            // Anti-extrapolation filter //
            ///////////////////////////////
            if(!mExtrapolate) {
                // Don't allow a final increment that is larger than any increment
                // at station points
                float maxInc = arma::max(lObs - (lY[e] + lYhat));
                float minInc = arma::min(lObs - (lY[e] + lYhat));
                if(y == mY) {
                    std::cout << "Increments: " << maxInc << " " << minInc << " " << currIncrement << std::endl;
                }

                // The increment for this member. currIncrement is the increment relative to
                // ensemble mean
                float memberIncrement = currIncrement - X(e);
                // Adjust increment if it gives a member increment that is outside the range
                // of the observation increments
                if(y == mY) {
                    std::cout << "Analysis increment: " << memberIncrement << " " << ensMean << " " << currIncrement << " " << X(e) << std::endl;
                }
                if(maxInc > 0 && memberIncrement > maxInc) {
                    currIncrement = maxInc + X(e);
                }
                else if(maxInc < 0 && memberIncrement > 0) {
                    currIncrement = 0 + X(e);
                }
                else if(minInc < 0 && memberIncrement < minInc) {
                    currIncrement = minInc + X(e);
                }
                else if(minInc > 0 && memberIncrement < 0) {
                    currIncrement = 0 + X(e);
                }
                if(y == mY) {
                    std::cout << "Final increment: " << currIncrement << " " << currIncrement - X(e) << std::endl;
                }
            }
            output[y][ei] = ensMean + currIncrement;
        }
    }

    // Back-transform
    /*
    #pragma omp parallel for
    for(int x = 0; x < nX; x++) {
        for(int y = 0; y < nY; y++) {
            float value = (*output)(y, x, e);
            if(gridpp::is_valid(value))
                (*output)(y, x, e) = invTransform(value);
        }
    }
    */
    return output;
}
