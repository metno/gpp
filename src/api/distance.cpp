#include "gridpp.h"
#include <iostream>

using namespace gridpp;

vec gridpp::distance(const Grid& grid, const Points& points, int num) {
    if(grid.get_coordinate_type() != points.get_coordinate_type())
        throw std::invalid_argument("Incompatible coordinate types");

    int size = points.size();
    vec output(size);
    vec2 ilats = grid.get_lats();
    vec2 ilons = grid.get_lons();
    CoordinateType coordinate_type = grid.get_coordinate_type();
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    for(int i = 0; i < size; i++) {
        ivec2 indices = grid.get_closest_neighbours(lats[i], lons[i], num);
        float max_dist = 0;
        for(int k = 0; k < indices.size(); k++) {
            int y_index = indices[k][0];
            int x_index = indices[k][1];
            float dist = gridpp::KDTree::calc_distance(lats[i], lons[i], ilats[y_index][x_index], ilons[y_index][x_index], coordinate_type);
            if(dist > max_dist)
                max_dist = dist;
        }
        output[i] = max_dist;
    }
    return output;
}

vec2 gridpp::distance(const Grid& igrid, const Grid& ogrid, int num) {
    if(igrid.get_coordinate_type() != ogrid.get_coordinate_type())
        throw std::invalid_argument("Incompatible coordinate types");

    ivec size = ogrid.size();
    vec2 output(size[0]);
    vec2 ilats = igrid.get_lats();
    vec2 ilons = igrid.get_lons();
    CoordinateType coordinate_type = igrid.get_coordinate_type();
    vec2 olats = ogrid.get_lats();
    vec2 olons = ogrid.get_lons();
    for(int i = 0; i < size[0]; i++) {
        output[i].resize(size[1], 0);
        for(int j = 0; j < size[1]; j++) {
            ivec2 indices = igrid.get_closest_neighbours(olats[i][j], olons[i][j], num);
            float max_dist = 0;
            for(int k = 0; k < indices.size(); k++) {
                ivec index = indices[k];
                // TODO: Deal with coordinate type
                float dist = gridpp::KDTree::calc_distance(ilats[index[0]][index[1]], ilons[index[0]][index[1]], olats[i][j], olons[i][j], coordinate_type);
                if(dist > max_dist)
                    max_dist = dist;
            }
            output[i][j] = max_dist;
        }
    }
    return output;
}

vec2 gridpp::distance(const Points& points, const Grid& grid, int num) {
    if(points.get_coordinate_type() != grid.get_coordinate_type())
        throw std::invalid_argument("Incompatible coordinate types");

    ivec size = grid.size();
    vec2 output(size[0]);
    vec lats = points.get_lats();
    vec lons = points.get_lons();
    CoordinateType coordinate_type = points.get_coordinate_type();
    vec2 ilats = grid.get_lats();
    vec2 ilons = grid.get_lons();
    for(int i = 0; i < size[0]; i++) {
        output[i].resize(size[1], 0);
        for(int j = 0; j < size[1]; j++) {
            ivec indices = points.get_closest_neighbours(ilats[i][j], ilons[i][j], num);
            float max_dist = 0;
            for(int k = 0; k < indices.size(); k++) {
                int index = indices[k];
                float dist = gridpp::KDTree::calc_distance(lats[index], lons[index], ilats[i][j], ilons[i][j], coordinate_type);
                if(dist > max_dist)
                    max_dist = dist;
            }
            output[i][j] = max_dist;
        }
    }
    return output;
}
