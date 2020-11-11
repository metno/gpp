#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec gridpp::apply_curve(const vec& fcst, const vec2& curve, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    if(curve.size() != 2)
        throw std::invalid_argument("Curve must have a first dimension size of 2");
    if(curve[0].size() == 0 || curve[1].size() == 0)
        throw std::invalid_argument("x and y vectors in curve cannot have size 0");
    if(curve[0].size() != curve[1].size())
        throw std::invalid_argument("x and y vectors in curve not the same size");

    int N = fcst.size();
    int C = curve[0].size();
    vec output(N, gridpp::MV);
    const vec& curve_fcst = curve[0];
    const vec& curve_ref = curve[1];
    float smallestObs  = curve_ref[0];
    float smallestFcst = curve_fcst[0];
    float largestObs   = curve_ref[C-1];
    float largestFcst  = curve_fcst[C-1];
    for(int i = 0; i < N; i++) {
        float input = fcst[i];
        // Linear interpolation within curve
        if(input > smallestFcst && input < largestFcst) {
            output[i] = gridpp::interpolate(input, curve_fcst, curve_ref);
        }
        // Extrapolate outside curve
        else {
            float nearestObs;
            float nearestFcst;
            gridpp::Extrapolation policy;
            if(input <= smallestFcst) {
                nearestObs  = smallestObs;
                nearestFcst = smallestFcst;
                policy = policy_below;
            }
            else {
                nearestObs  = largestObs;
                nearestFcst = largestFcst;
                policy = policy_above;
            }
            float slope = 1;
            if(policy == gridpp::Zero) {
                slope = 0;
            }
            if(policy == gridpp::OneToOne || C <= 1) {
                slope = 1;
            }
            else if(policy == gridpp::MeanSlope) {
                float dObs  = largestObs - smallestObs;
                float dFcst = largestFcst - smallestFcst;
                slope = dObs / dFcst;
            }
            else if(policy == gridpp::NearestSlope) {
                float dObs;
                float dFcst;
                if(input <= smallestFcst) {
                    dObs  = curve_ref[1] - curve_ref[0];
                    dFcst = curve_fcst[1] - curve_fcst[0];
                }
                else {
                    dObs  = curve_ref[C-1] - curve_ref[C-2];
                    dFcst = curve_fcst[C-1] - curve_fcst[C-2];
                }
                slope = dObs / dFcst;
            }
            output[i] = nearestObs + slope * (input - nearestFcst);
        }
    }
    return output;
}
vec2 gridpp::apply_curve(const vec2& fcst, const vec2& curve, gridpp::Extrapolation policy_below, gridpp::Extrapolation policy_above) {
    int nY = fcst.size();
    int nX = fcst[0].size();
    vec2 output(nY);
    for(int y = 0; y < nY; y++) {
        output[y] = apply_curve(fcst[y], curve, policy_below, policy_above);
    }
    return output;
}
vec2 gridpp::monotonize_curve(vec2 curve) {
    if(curve.size() != 2)
        throw std::invalid_argument("Curve must have a first dimension size of 2");
    if(curve[0].size() == 0 || curve[1].size() == 0)
        throw std::invalid_argument("x and y vectors in curve cannot have size 0");
    if(curve[0].size() != curve[1].size())
        throw std::invalid_argument("x and y vectors in curve not the same size");

    bool debug = false;

    // Remove missing
    int N = curve[0].size();
    ivec keep;
    keep.reserve(N); // Array of indices to keep (both x and y are valid)
    for(int i = 0; i < N; i++) {
        if(gridpp::is_valid(curve[0][i]) && gridpp::is_valid(curve[1][i]))
            keep.push_back(i);
    }
    if(keep.size() != N) {
        vec2 curve_copy = curve;
        curve[0].clear();
        curve[1].clear();
        curve[0].resize(keep.size(), 0);
        curve[1].resize(keep.size(), 0);
        for(int i = 0; i < keep.size(); i++) {
            curve[0][i] = curve_copy[0][keep[i]];
            curve[1][i] = curve_copy[1][keep[i]];
        }
    }
    N = curve[0].size();
    ivec new_indices;
    new_indices.reserve(N);

    int start = 1;
    float prev = curve[0][0];

    /* A deviation is a part of the curve is not monotonic. This occurs when a set of x values are
       not in increasing order. Remove all points that are inside the non-increasing setion.
    */
    bool deviation = false;
    int deviation_index = 0;
    float x_min = curve[0][0]; // The lowest x value in the deviation section
    float x_max = curve[0][0]; // The highest x value in the deivation section

    new_indices.push_back(0);
    float tol = 0.1;
    for(int i = start; i < N; i++) {
        float x = curve[0][i];
        float y = curve[1][i];
        assert(gridpp::is_valid(x));
        assert(gridpp::is_valid(y));
        if(deviation) {
            // Currently inside a deviation section
            if(x < x_min) {
                x_min = x;
            }
            if(x > x_max + tol) {
                // We are past the deviation section
                x_max = x;
                if(debug)
                    std::cout << "End of loop at " << i << " x=" << x << "(" << x_min << "," << x_max << ")" << std::endl;
                // Remove points inside the deviation section, that is all points above x_min
                for(int j = new_indices.size() - 1; j >= 0; j--) {
                    int index = new_indices[j];
                    if(curve[0][index] < x_min - tol)
                        break;
                    else {
                        if(debug)
                            std::cout << "Removing index=" << index << " x=" << curve[0][index] << std::endl;
                        new_indices.pop_back();
                    }
                }
                new_indices.push_back(i);
                deviation = false;
                prev = x;
                x_max = x;
            }
        }
        else {
            if(x <= prev + tol) {
                // We have a new deviation
                deviation = true;
                deviation_index = i;
                x_min = x;
                if(debug)
                    std::cout << "Detecting loop from " << i << " x=" << x << " " << x_max << std::endl;
            }
            else {
                new_indices.push_back(i);
                prev = x;
                x_max = x;
            }
        }
    }
    if(deviation) {
        if(debug)
            std::cout << "Finished in a deviation " << deviation_index << std::endl;
        // Find latest x point
        for(int j = new_indices.size() - 1; j >= 0; j--) {
            if(debug)
                std::cout << new_indices[j] << " " << deviation_index << std::endl;
            int index = new_indices[j];
            if(curve[0][index] >= x_min)
                new_indices.pop_back();
        }
    }

    vec2 new_curve(2);
    new_curve[0].resize(new_indices.size());
    new_curve[1].resize(new_indices.size());
    for(int i = 0; i < new_indices.size(); i++) {
        assert(curve[0].size() > new_indices[i]);
        new_curve[0][i] = curve[0][new_indices[i]];
        new_curve[1][i] = curve[1][new_indices[i]];
    }
    return new_curve;
}
