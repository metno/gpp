#include "gridpp.h"
#include <iostream>

using namespace gridpp;

gridpp::Grid::Grid() {
    vec lats;
    vec lons;
    mTree = KDTree(lats, lons);

}
gridpp::Grid::Grid(vec2 lats, vec2 lons, vec2 elevs, vec2 lafs, CoordinateType type) {
    mLats = lats;
    mLons = lons;
    mElevs = elevs;
    mLafs = lafs;
    mX = lats[0].size();
    int N = lats.size() * lats[0].size();
    vec lats0(N);
    vec lons0(N);
    int count = 0;
    for(int i = 0; i < lats.size(); i++) {
        for(int j = 0; j < lats[0].size(); j++) {
            lats0[count] = lats[i][j];
            lons0[count] = lons[i][j];
            count++;
        }
    }
    KDTree test = KDTree(lats0, lons0, type);
    mTree = test;

    if(mElevs.size() != lats.size() or mElevs[0].size() != lats[0].size()) {
        mElevs.clear();
        mElevs.resize(lats.size());
        for(int i = 0; i < lats.size(); i++) {
            mElevs[i].resize(lats[0].size(), gridpp::MV);
        }
    }
    if(mLafs.size() != lats.size() or mLafs[0].size() != lats[0].size()) {
        mLafs.clear();
        mLafs.resize(lats.size());
        for(int i = 0; i < lats.size(); i++) {
            mLafs[i].resize(lats[0].size(), gridpp::MV);
        }
    }
}

int gridpp::Grid::get_num_neighbours(float lat, float lon, float radius) const {
    ivec indices = mTree.get_neighbours(lat, lon, radius);
    return indices.size();
}

ivec2 gridpp::Grid::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const {
    ivec indices = mTree.get_neighbours_with_distance(lat, lon, radius, distances);
    return get_indices(indices);
}

ivec2 gridpp::Grid::get_neighbours(float lat, float lon, float radius) const {
    ivec indices = mTree.get_neighbours(lat, lon, radius);
    return get_indices(indices);
}

ivec2 gridpp::Grid::get_closest_neighbours(float lat, float lon, int num) const {
    ivec indices = mTree.get_closest_neighbours(lat, lon, num);
    return get_indices(indices);
}
ivec gridpp::Grid::get_nearest_neighbour(float lat, float lon) const {
    ivec2 I = get_closest_neighbours(lat, lon, 1);
    if(I.size() > 0)
        return I[0];
    else
        return ivec();
}
vec2 gridpp::Grid::get_lats() const {
    return mLats;
}
vec2 gridpp::Grid::get_lons() const {
    return mLons;
}
vec2 gridpp::Grid::get_elevs() const {
    return mElevs;
}
vec2 gridpp::Grid::get_lafs() const {
    return mLafs;
}
vec2 gridpp::Grid::get_2d(vec input) const {
    int Y = input.size() / mX;
    vec2 output(Y);
    int count = 0;
    for(int i = 0; i < Y; i++) {
        output[i].resize(mX, 0);
        for(int j = 0; j < mX; j++) {
            output[i][j] = input[count];
            count++;
        }
    }
    return output;
}
ivec gridpp::Grid::get_indices(int index) const {
    ivec results(2, 0);
    assert(index < mTree.size());
    results[0] = index / mX;
    results[1] = index % mX;
    return results;
}
ivec2 gridpp::Grid::get_indices(ivec indices) const {
    ivec2 results(indices.size());
    for(int i = 0; i < indices.size(); i++) {
        results[i] = get_indices(indices[i]);
    }
    return results;
}
ivec gridpp::Grid::size() const {
    ivec indices(2, 0);
    if(mTree.size() > 0) {
        indices = get_indices(mTree.size() - 1);
        indices[0]++;
        indices[1]++;
    }
    return indices;
}
Points gridpp::Grid::to_points() const {
    int Y = size()[0];
    int X = size()[1];
    vec elevs(Y * X);
    vec lafs(Y * X);
    int count = 0;
    for(int y = 0; y < Y; y++) {
        for(int x = 0; x < X; x++) {
            elevs[count] = mElevs[y][x];
            lafs[count] = mLafs[y][x];
            count++;
        }
    }
    return gridpp::Points(mTree, elevs, lafs);
}
CoordinateType gridpp::Grid::get_coordinate_type() const {
    return mTree.get_coordinate_type();
}
bool gridpp::Grid::get_box(float lat, float lon, int& Y1, int& X1, int& Y2, int& X2) const {
    int xdir = 1;
    int ydir = 1;
    ivec nn = get_nearest_neighbour(lat, lon);
    if(nn.size() != 2)
        return false;
    int Y = nn[0];
    int X = nn[1];
    int nY = size()[0];
    int nX = size()[1];
    Y1 = -1;
    Y2 = -1;
    X1 = -1;
    X2 = -1;
    if(nX <= 1 || nY <= 1)
        return false;
    if(Y == 0) {
        Y1 = 0;
        Y2 = 1;
    }
    else if(Y == nY - 1) {
        Y1 = nY - 2;
        Y2 = nY - 1;
    }
    if(X == 0) {
        X1 = 0;
        X2 = 1;
    }
    else if(X == nX - 1) {
        X1 = nX - 2;
        X2 = nX - 1;
    }

    int it = 0;
    bool found = false;
    while(it < 4) {
        xdir = -1 + 2 * (it %2);
        ydir = -1 + 2 * (it < 2);
        it++;
        if(Y == 0 && ydir == -1)
            continue;
        else if(Y == nY - 1 && ydir == 1)
            continue;
        else if(X == 0 && xdir == -1)
            continue;
        else if(X == nX - 1 && xdir == 1)
            continue;
        // std::cout << "Iteration " << it << " " << ydir << " " << xdir << std::endl;
        Point A(mLats[Y][X], mLons[Y][X]);
        Point B(mLats[Y + ydir][X], mLons[Y + ydir][X]);
        Point C(mLats[Y + ydir][X + xdir], mLons[Y + ydir][X + xdir]);
        Point D(mLats[Y][X + xdir], mLons[Y][X + xdir]);
        Point curr(lat, lon);
        // std::cout << A.lat << "," << B.lat << "," << C.lat << "," << D.lat << "," << curr.lat << std::endl;
        // std::cout << A.lon << "," << B.lon << "," << C.lon << "," << D.lon << "," << curr.lon << std::endl;
        bool inrectangle = point_in_rectangle(A, B, C, D, curr);
        // std::cout << " inrectangle: " << inrectangle << std::endl;
        if(inrectangle) {
            found = true;
            break;
        }
    }
    if(!found)
        return false;

    if(xdir == 1) {
        X1 = X;
        X2 = X + 1;
    }
    else if(xdir == -1) {
        X1 = X - 1;
        X2 = X;
    }
    if(ydir == 1) {
        Y1 = Y;
        Y2 = Y + 1;
    }
    else if(ydir == -1) {
        Y1 = Y - 1;
        Y2 = Y;
    }
    return true;
}
