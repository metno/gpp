#include "gridpp.h"
#include <iostream>
#include <algorithm>
#include <boost/math/tools/minima.hpp>
using boost::math::tools::brent_find_minima;

using namespace gridpp;

namespace {
    float partial(float ref_diff, float fcst_diff, float sigma) {
        float total = 0;
        float ref_total = 0;
        float fcst_total = 1;
        if(ref_diff > 0)
            ref_total = 1;
        else
            ref_total = 1 + (ref_diff) / sigma;
        if(fcst_diff > 0)
            fcst_total = 1;
        else
            fcst_total = 1 + (fcst_diff) / sigma;

        return ref_total * fcst_total;
    }
    class Data {
        public:
            Data(const vec& ref, const vec& fcst, float threshold, gridpp::Metric metric) {
                this->ref = ref;
                this->fcst = fcst;
                this->metric = metric;
                this->threshold = threshold;
                int N = ref.size();
                min = gridpp::MV;
                max = gridpp::MV;
                for(int i = 0; i < N; i++) {
                    if(!gridpp::is_valid(min) || fcst[i] < min) {
                        min = fcst[i];
                    }
                    if(!gridpp::is_valid(max) || fcst[i] > max) {
                        max = fcst[i];
                    }
                }
            };
            vec ref;
            vec fcst;
            float threshold;
            gridpp::Metric metric;
            float min;
            float max;
    };
    class Score {
        public:
            Score(float a, float b, float c, float d, gridpp::Metric metric) {
                this->a = a;
                this->b = b;
                this->c = c;
                this->d = d;
                this->metric = metric;
            };
            float a;
            float b;
            float c;
            float d;
            gridpp::Metric metric;
    };
    class Score_func {
        public:
            Score_func(const Data& idata) : data(idata) {
            };
            float operator()(double threshold) {
                float s = - gridpp::calc_score(data.ref, data.fcst, data.threshold, threshold, data.metric);
                // if(s > -0.0001)
                //     s = 0;
                // std::cout << threshold << " " << s << std::endl;
                return s;
            }
            float calc_uncertainty(double threshold) {
                float a = 0, b = 0, c = 0, d = 0;
                int N = data.fcst.size();
                for(int t = 0; t < N; t++) {
                    if((data.ref[t] > data.threshold) && (data.fcst[t] > threshold))
                        a++;
                    else if((data.ref[t] <= data.threshold) && (data.fcst[t] > threshold))
                        b++;
                    else if((data.ref[t] > data.threshold) && (data.fcst[t] <= threshold))
                        c++;
                    else if((data.ref[t] <= data.threshold) && (data.fcst[t] <= threshold))
                        d++;
                }
                float s = - gridpp::calc_score(a, b, c, d, data.metric);

                float zero = 0;
                float s0 = - gridpp::calc_score(std::max(zero, a-1), std::max(zero, b-1), std::max(zero, c-1), std::max(zero, d-1), data.metric);
                float diff = fabs(fabs(s - s0) / std::min(s, s0));
                if(diff == -gridpp::MV)
                    diff = gridpp::MV;
                // if(a + b < 50)
                //     diff = 1;
                return diff;
            }
            const Data& data;
    };
}

vec2 gridpp::metric_optimizer_curve(const vec& ref, const vec& fcst, const vec& thresholds, gridpp::Metric metric) {
    if(ref.size() != fcst.size())
        throw std::invalid_argument("ref and fcst not the same size");

    int N = thresholds.size();
    vec x(N);
    vec y(N);
    float min = 0;
    float max = 10;
    vec2 results(2);
    results[0].reserve(N);
    results[1].reserve(N);

    float last_diff = 0;
    for(int i = 0; i < N; i++) {
        float value = gridpp::get_optimal_threshold(ref, fcst, thresholds[i], metric);
        if(gridpp::is_valid(value)) {
            results[0].push_back(value);
            results[1].push_back(thresholds[i]);
        }
    }
    return results;
}

float gridpp::get_optimal_threshold(const vec& ref, const vec& fcst, float threshold, gridpp::Metric metric) {
    bool remove_near_zero = true;
    bool remove_uncertain = false;
    bool remove_at_boundary = true;
    // int bits = std::numeric_limits<double>::digits;
    int bits = 32;
    Data data(ref, fcst, threshold, metric);
    Score_func func(data);

    /* Compute values at even intervals, to find where to start optimization search */
    int B = 10;
    vec bins(B);
    vec values(B);
    int min_index = 0;
    float min_value = func(data.min);
    bins[0] = data.min;
    for(int b = 1; b < B; b++) {
        bins[b] = data.min + (data.max - data.min) / (B - 1) * b;
        float curr_x = data.min + (data.max - data.min) / (B - 1) * b;
        float curr_value = func(curr_x);
        // std::cout << curr_x << " " << curr_value << " " << min_index << " " << min_value << std::endl;
        if(curr_value < min_value) {
            min_index = b;
            min_value = curr_value;
        }
    }
    int left_index = min_index > 0 ? min_index - 1 : 0;
    int right_index = min_index < B - 1 ? min_index + 1 : B - 1;
    // std::cout << "MIN: " << data.min << " " << "MAX: " << data.max << std::endl;
    // std::cout << "MIN: " << bins[left_index] << " " << "MAX: " << bins[right_index] << std::endl;
    std::pair<double, double> r = brent_find_minima<Score_func, double>(func, bins[left_index], bins[right_index], bits);
    float x = r.first;
    float score = - r.second;
    if(!gridpp::is_valid(score))
        x = gridpp::MV;
    if(remove_near_zero && score <= 0.0001)
        x = gridpp::MV;
    if(remove_uncertain) {
        float diff = func.calc_uncertainty(r.first);
        if(!gridpp::is_valid(diff))
            x = gridpp::MV;
        else if(diff > 1e-3)
            x = gridpp::MV;
        // std::cout << threshold << " -> " << r.first << " " << score << " " << diff << " " << x << std::endl;
    }
    if(remove_at_boundary) {
        float s0 = - func(data.min);
        float s1 = - func(data.max);
        if (fabs(r.second - s0) < 0.001 || fabs(r.second - s1) < 0.001)
            x = gridpp::MV;
    }
    return x;
}
float gridpp::calc_score(const vec& ref, const vec& fcst, float threshold, gridpp::Metric metric) {
    return calc_score(ref, fcst, threshold, threshold, metric);
}
float gridpp::calc_score(const vec& ref, const vec& fcst, float threshold, float fthreshold, gridpp::Metric metric) {
    float a = 0, b = 0, c = 0, d = 0;
    int N = fcst.size();
    for(int t = 0; t < N; t++) {
        if((ref[t] > threshold) && (fcst[t] > fthreshold))
            a++;
        else if((ref[t] <= threshold) && (fcst[t] > fthreshold))
            b++;
        else if((ref[t] > threshold) && (fcst[t] <= fthreshold))
            c++;
        else if((ref[t] <= threshold) && (fcst[t] <= fthreshold))
            d++;
    }
    return calc_score(a, b, c, d, metric);
}
float gridpp::calc_score(float a, float b, float c, float d, gridpp::Metric metric) {
    if(metric == gridpp::Ets) {
        float N = a + b + c + d;
        float ar = (a + b) / 1.0 / N * (a + c);
        if(a + b + c - ar == 0)
            return gridpp::MV;
        return (a - ar) / 1.0 / (a + b + c - ar);
    }
    else if(metric == gridpp::Ts) {
        float N = a + b + c + d;
        return a / 1.0 / (a + b + c);
    }
    else if(metric == gridpp::Pc) {
        float N = a + b + c + d;
        return (a + d) / N;
    }
    else if(metric == gridpp::Kss) {
        if((a + c) * (b + d) == 0)
            return gridpp::MV;
        return (a * d - b * c) * 1.0 / ((a + c) * (b + d));
    }
    else if(metric == gridpp::Bias) {
        if(b + c == 0)
            return gridpp::MV;
        // std::cout << a << " " << b << " " << c << " " << d << std::endl;
        float value = 1-float(fabs(b - c))/(b + c);
        return value;
    }
    else if(metric == gridpp::Hss) {
        float denom = ((a + c) * (c + d) + (a + b) * (b + d));
        if(denom == 0)
            return gridpp::MV;
        return 2.0 * (a * d - b * c) / denom;
    }
}
