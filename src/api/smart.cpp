#include "gridpp.h"

using namespace gridpp;

namespace {
    template<class T1, class T2> struct sort_pair_first {
        bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) {
            return left.first < right.first;
        };
    };
}
vec2 gridpp::smart(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, int num, const gridpp::StructureFunction& structure) {
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    vec2 ielevs = igrid.get_elevs();
    vec2 ilafs = igrid.get_lafs();

    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    vec2 oelevs = ogrid.get_elevs();
    vec2 olafs = ogrid.get_lafs();

    int nLat = ogrid.size()[0];
    int nLon = ogrid.size()[1];

    vec2 output;
    output.resize(nLat);
    for(int i = 0; i < nLat; i++)
        output[i].resize(nLon);

    float dist = structure.localization_distance();

    // #pragma omp parallel for
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec2 indices = igrid.get_neighbours(olats[i][j], olons[i][j], dist);
            std::vector<std::pair<float,int> > rhos;
            rhos.reserve(indices.size());
            Point p1 = ogrid.get_point(i, j);
            for(int s = 0; s < indices.size(); s++) {
                int ii = indices[s][0];
                int jj = indices[s][1];
                Point p2 = igrid.get_point(ii, jj);
                float rho = structure.corr(p1, p2);
                rhos.push_back(std::pair<float,int>(rho, s));
            }
            std::sort(rhos.begin(), rhos.end(), ::sort_pair_first<float,int>());
            assert(rhos[0].first <= rhos[rhos.size()-1].first);
            float sum = 0;
            int count = 0;
            for(int s = 0; s < std::min(num, int(indices.size())); s++) {
                int index = rhos[rhos.size() - 1 - s].second;
                int ii = indices[index][0];
                int jj = indices[index][1];
                float value = ivalues[ii][jj];
                sum += value;
                count++;
            }
            float value = gridpp::MV;
            if(count > 0)
                value = sum / count;
            output[i][j] = value;
        }
    }
    return output;
}
#if 0
vec2 gridpp::smart(const Grid& igrid, const Grid& ogrid, const vec2 ivalues, int radius, float num, float min_elev_diff) {
    vec2 iOutputLats = ogrid.get_lats();
    vec2 iOutputLons = ogrid.get_lons();

    int nLat = iOutputLats.size();
    int nLon = iOutputLats[0].size();

    vec2 output(nLat);
    for(int i = 0; i < nLat; i++) {
        output[i].resize(nLon, 0);
    }

    // #pragma omp parallel for
    for(int i = 0; i < nLat; i++) {
        for(int j = 0; j < nLon; j++) {
            ivec2 indices = ::get_smart_neighbour(igrid, iOutputLats[i][j], iOutputLons[i][j]);
            int N = indices.size();
            for(int n = 0; n < N; n++) {
                int I = indices[n][0];
                int J = indices[n][1];
                output[i][j] += ivalues[I][J];
            }
            output[i][j] / N;
        }
    }
    return output;
namespace {
def get_smart_neighbour(const Grid& grid, float lat, float lon, float altitude, int radius, float num, float min_elev_diff) {
    util::not_implemented_error();
    /*
    ivec2 indices = grid.get_neighbours(lat, lon, radius);
    int N = indices.size();
    vec
    for(int n = 0; n < N; n++) {
        int I = indices[n][0];
        int J = indices[n][1];
        float curr_altitude = grid.get_elevs()[I][J];
        bool isWithinMinElev = util::is_valid(altitude) && util::is_valid(curr_altitude) && util::is_valid(mmin_elev_diff) && fabs(altitude - curr_altitude) <= mmin_elev_diff;
        // Compute elevation differences on the stencil surrounding the current point
        std::vector<std::pair<int, float> > elevDiff; // elevation difference of points in stencil
        // but placed in 1D array, to simplify sort
        elevDiff.reserve(numSearch);
        // Keep track of which I/J corresponds to indices in elevDiff
        std::vector<int> Ilookup;
        std::vector<int> Jlookup;
        Ilookup.reserve(numSearch);
        Jlookup.reserve(numSearch);

        int index = 0;
        for(int ii = std::max(0, Ic-mRadius); ii <= std::min(iFrom.getNumY()-1, Ic+mRadius); ii++) {
            for(int jj = std::max(0, Jc-mRadius); jj <= std::min(iFrom.getNumX()-1, Jc+mRadius); jj++) {
                float ielev = ielevs[ii][jj];
                float diff = 1e10;
                if(Util::is_valid(ielev) && Util::is_valid(oelev))
                    diff = abs(ielev - oelev);
                elevDiff.push_back(std::pair<int, float>(index, diff));
                Ilookup.push_back(ii);
                Jlookup.push_back(jj);
                index++;
            }
        }

        std::sort(elevDiff.begin(), elevDiff.end(), Util::sort_pair_second<int,float>());

        std::stringstream ss;
        ss << "Smart neighbours for " << i << " " << j << " " << oelev;
        Util::info(ss.str());

        // Use nearest neighbour if all fails
        if(elevDiff.size() == 0) {
            iI[i][j].push_back(Ic);
            iJ[i][j].push_back(Jc);
        }
        else {
            int N = std::min((int) elevDiff.size(), mNum);
            iI[i][j].resize(N, Util::MV);
            iJ[i][j].resize(N, Util::MV);

            for(int n = 0; n < N; n++) {
                int index = elevDiff[n].first;
                iI[i][j][n] = Ilookup[index];
                iJ[i][j][n] = Jlookup[index];
                std::stringstream ss;
                ss << "   " << iI[i][j][n] << " " << iJ[i][j][n] << " " << elevDiff[n].second;
                Util::info(ss.str());
            }
        }
    }
    */
}
}
#endif
