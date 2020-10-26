#include "gridpp.h"

using namespace gridpp;

gridpp::Points::Points() {
    // TODO: Deal with the empty case. Can't really find nearest neighbours then

}
gridpp::Points::Points(vec lats, vec lons, vec elevs, vec lafs, CoordinateType type) {
    int N = lats.size();
    if(lons.size() != N)
        throw std::invalid_argument("Cannot create points with unequal lat and lon sizes");
    if(elevs.size() != 0 && elevs.size() != N)
        throw std::invalid_argument("'elevs' must either be size 0 or the same size at lats/lons");
    if(lafs.size() != 0 && lafs.size() != N)
        throw std::invalid_argument("'lafs' must either be size 0 or the same size at lats/lons");
    mLats = lats;
    mLons = lons;
    mElevs = elevs;
    mLafs = lafs;
    KDTree tree = KDTree(lats, lons, type);
    mTree = tree;
    if(mElevs.size() != N) {
        mElevs.clear();
        mElevs.resize(N, gridpp::MV);
    }
    if(mLafs.size() != N) {
        mLafs.clear();
        mLafs.resize(N, gridpp::MV);
    }
}
gridpp::Points::Points(KDTree tree, vec elevs, vec lafs) {
    mElevs = elevs;
    mLafs = lafs;
    mTree = tree;
    mLats = tree.get_lats();
    mLons = tree.get_lons();
}

int gridpp::Points::get_num_neighbours(float lat, float lon, float radius) const {
    return mTree.get_num_neighbours(lat, lon, radius);
}

ivec gridpp::Points::get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const {
    return mTree.get_neighbours_with_distance(lat, lon, radius, distances);
}

ivec gridpp::Points::get_neighbours(float lat, float lon, float radius) const {
    return mTree.get_neighbours(lat, lon, radius);
}

ivec gridpp::Points::get_closest_neighbours(float lat, float lon, int num) const {
    return mTree.get_closest_neighbours(lat, lon, num);
}
int gridpp::Points::get_nearest_neighbour(float lat, float lon) const {
    ivec I = get_closest_neighbours(lat, lon, 1);
    if(I.size() > 0)
        return I[0];
    else
        return -1;
}
vec gridpp::Points::get_lats() const {
    return mLats;
}
vec gridpp::Points::get_lons() const {
    return mLons;
}
vec gridpp::Points::get_elevs() const {
    return mElevs;
}
vec gridpp::Points::get_lafs() const {
    return mLafs;
}
int gridpp::Points::size() const {
    return mLats.size();
}
ivec gridpp::Points::get_in_domain_indices(const gridpp::Grid& grid) const {
    ivec indices;
    indices.reserve(size());
    ivec gsize = grid.size();
    int Y = gsize[0];
    int X = gsize[1];

    for(int s = 0; s < size(); s++) {
        int Y1, X1, Y2, X2;
        bool inside  = grid.get_box(mLats[s], mLons[s], Y1, X1, Y2, X2);
        if(inside) {
            indices.push_back(s);
        }
    }
    return indices;
}
gridpp::Points gridpp::Points::get_in_domain(const gridpp::Grid& grid) const {
    vec lats, lons, elevs, lafs;
    lats.reserve(size());
    lons.reserve(size());
    elevs.reserve(size());
    lafs.reserve(size());
    ivec indices = get_in_domain_indices(grid);
    int S = indices.size();
    for(int s = 0; s < S; s++) {
        lats.push_back(mLats[indices[s]]);
        lons.push_back(mLons[indices[s]]);
        elevs.push_back(mElevs[indices[s]]);
        lafs.push_back(mLafs[indices[s]]);
    }
    const gridpp::Points new_points(lats, lons, elevs, lafs);
    return new_points;
}
gridpp::Points& gridpp::Points::operator=(gridpp::Points other) {
    std::swap(mLats, other.mLats);
    std::swap(mLons, other.mLons);
    std::swap(mElevs, other.mElevs);
    std::swap(mLafs, other.mLafs);
    std::swap(mTree, other.mTree);
    return *this;
}
gridpp::Points::Points(const gridpp::Points& other) {
    mLats = other.mLats;
    mLons = other.mLons;
    mElevs = other.mElevs;
    mLafs = other.mLafs;
    mTree = KDTree(mLats, mLons, mTree.get_coordinate_type());
}
CoordinateType gridpp::Points::get_coordinate_type() const {
    return mTree.get_coordinate_type();
}
