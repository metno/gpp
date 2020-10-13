#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec2 gridpp::fill(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool outside) {
    if(gridpp::compatible_size(igrid, input))
        throw std::invalid_argument("Grid size is not the same as values");
    if(points.size() != radii.size())
        throw std::invalid_argument("Points size is not the same as radii size");
    for(int i = 0; i < radii.size(); i++) {
        if(radii[i] < 0)
            throw std::invalid_argument("All radius sizes must be 0 or greater");
    }

    double s_time = gridpp::clock();
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    vec2 output;
    if(outside) {
        output.resize(input.size());
        for(int i = 0; i < input.size(); i++) {
            output[i].resize(input[i].size(), value);
        }
    }
    else {
        output = input;
    }
    for(int i = 0; i < points.size(); i++) {
        ivec2 I = igrid.get_neighbours(lats[i], lons[i], radii[i]);
        // std::cout << " " << i << " " << I.size() << " " << lats[i] << " " << lons[i] << " " << radii[i] << std::endl;
        for(int j = 0; j < I.size(); j++) {
            if(outside)
                output[I[j][0]][I[j][1]] = input[I[j][0]][I[j][1]];
            else
                output[I[j][0]][I[j][1]] = value;
        }
    }
    // std::cout << "Fill time: " << gridpp::clock() - s_time << std::endl;
    return output;
}

vec2 gridpp::fill_missing(const vec2& values) {
    int Y = values.size();
    int X = values[0].size();
    vec2 results_y(Y);
    vec2 results_x(Y);
    int count_y = 0;
    for(int y = 0; y < Y; y++) {
        results_y[y].resize(X, gridpp::MV);
        int last = 0;
        int next = -1;
        for(int x = 0; x < X; x++) {
            float curr = values[y][x];
            if(!gridpp::is_valid(curr)) {
                // Find next
                if(next < x) {
                    for(next = x; next < X; next++) {
                        if(gridpp::is_valid(values[y][next])) {
                            break;
                        }
                    }
                }
                if(next >= X) {
                    continue;
                }
                // Interpolate
                float value_last = values[y][last];
                float value_next = values[y][next];
                results_y[y][x] = (value_last) + (value_next - value_last) * (x - last) / (next - last);
                count_y++;
            }
            else {
                last = x;
                results_y[y][x] = values[y][x];
            }
        }
    }

    int count_x = 0;
    for(int y = 0; y < Y; y++) {
        results_x[y].resize(X, gridpp::MV);
    }

    for(int x = 0; x < X; x++) {
        int last = 0;
        int next = -1;
        for(int y = 0; y < Y; y++) {
            float curr = values[y][x];
            if(!gridpp::is_valid(curr)) {
                // Find next
                if(next < y) {
                    for(next = y; next < Y; next++) {
                        if(gridpp::is_valid(values[next][x])) {
                            break;
                        }
                    }
                }
                // Interpolate
                if(next >= X) {
                    continue;
                }
                float value_last = values[last][x];
                float value_next = values[next][x];
                results_x[y][x] = (value_last) + (value_next - value_last) * (y - last) / (next - last);
                count_x++;
            }
            else {
                last = y;
                results_x[y][x] = values[y][x];
            }
        }
    }
    vec2 results(Y);
    for(int y = 0; y < Y; y++) {
        results[y].resize(X, gridpp::MV);
        for(int x = 0; x < X; x++) {
            int count = 0;
            float total = 0;
            if(gridpp::is_valid(results_y[y][x])) {
                total += results_y[y][x];
                count++;
            }
            if(gridpp::is_valid(results_x[y][x])) {
                total += results_x[y][x];
                count++;
            }
            if(count > 0)
                results[y][x] = total / count;
        }
    }
    // std::cout << "Number of infilled values: " << count_y << " " << count_x << std::endl;
    return results;
}
