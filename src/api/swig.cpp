#include "gridpp.h"
#include <iostream>

using namespace gridpp;

float* gridpp::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
float gridpp::test_vec_input(const vec& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++)
        total += input[i];
    return total;
}
int gridpp::test_ivec_input(const ivec& input) {
    int total = 0;
    for(int i = 0; i < input.size(); i++) {
        total += int(input[i]);
    }
    return total;
}
float gridpp::test_vec2_input(const vec2& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++) {
        for(int j = 0; j < input[i].size(); j++) {
            total += input[i][j];
        }
    }
    return total;
}
float gridpp::test_vec3_input(const vec3& input) {
    float total = 0;
    for(int i = 0; i < input.size(); i++) {
        for(int j = 0; j < input[i].size(); j++) {
            for(int k = 0; k < input[i][j].size(); k++) {
                total += input[i][j][k];
            }
        }
    }
    return total;
}
vec gridpp::test_vec_output() {
    vec output(3, swig_default_value);
    return output;
}
vec2 gridpp::test_vec2_output() {
    vec2 output(3);
    for(int i = 0; i < 3; i++)
        output[i].resize(3, swig_default_value);
    return output;
}
vec3 gridpp::test_vec3_output() {
    vec3 output(3);
    for(int i = 0; i < 3; i++) {
        output[i].resize(3);
        for(int j = 0; j < 3; j++)
            output[i][j].resize(3, swig_default_value);
    }
    return output;
}

ivec gridpp::test_ivec_output() {
    ivec output(3, swig_default_value);
    return output;
}
ivec2 gridpp::test_ivec2_output() {
    ivec2 output(3);
    for(int i = 0; i < 3; i++)
        output[i].resize(3, swig_default_value);
    return output;
}

ivec3 gridpp::test_ivec3_output() {
    ivec3 output(3);
    for(int i = 0; i < 3; i++) {
        output[i].resize(3);
        for(int j = 0; j < 3; j++)
            output[i][j].resize(3, swig_default_value);
    }
    return output;
}

float gridpp::test_vec_argout(vec& distances) {
    distances.clear();
    distances.resize(10, swig_default_value);
    return 0;
}
float gridpp::test_vec2_argout(vec2& distances) {
    distances.clear();
    distances.resize(10);
    for(int i = 0; i < 10; i++)
        distances[i].resize(10, swig_default_value);
    return 0;
}
void gridpp::test_not_implemented_exception() {
    throw gridpp::not_implemented_exception();
}
