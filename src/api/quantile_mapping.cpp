#include "gridpp.h"

using namespace gridpp;

vec2 gridpp::quantile_mapping_curve(const vec& ref, const vec& fcst, vec quantiles) {
    if(ref.size() != fcst.size())
        throw std::invalid_argument("ref and fcst must be of the same size");

    if(quantiles.size() > 0) {
        for(int i = 0; i < quantiles.size(); i++) {
            float curr = quantiles[i];
            if(!gridpp::is_valid(curr) || curr > 1 || curr < 0)
                throw std::invalid_argument("Quantiles must be >= 0 and <= 1");
        }
    }
    vec2 curve(2);
    if(ref.size() == 0)
        return curve;
    else if(ref.size() == 1) {
        curve[0] = fcst;
        curve[1] = ref;
        return curve;
    }

    vec ref_sort = ref;
    std::sort(ref_sort.begin(), ref_sort.end());
    vec fcst_sort = fcst;
    std::sort(fcst_sort.begin(), fcst_sort.end());
    int N = quantiles.size();
    int S = fcst_sort.size();
    if(N == 0) {
        curve[0] = fcst_sort;
        curve[1] = ref_sort;
    }
    else {
        curve[0].resize(N);
        curve[1].resize(N);
        for(int i = 0; i < N; i++) {
            int index = quantiles[i] * (S - 1);
            curve[0][i] = fcst[index];
            curve[1][i] = ref[index];
        }
    }
    return curve;
}
