#include "gridpp.h"

using namespace gridpp;

gridpp::StructureFunction::StructureFunction(float localization_distance) {
    if(!gridpp::is_valid(localization_distance) || localization_distance < 0)
        throw std::invalid_argument("Invalid 'h' in structure");

    mLocalizationDistance = localization_distance;
}
float gridpp::StructureFunction::corr_background(const Point& p1, const Point& p2) const {
    return corr(p1, p2);
}
float gridpp::StructureFunction::barnes_rho(float dist, float length) const {
    if(!gridpp::is_valid(length) || length == 0)
        // Disabled
        return 1;
    if(!gridpp::is_valid(dist))
        return 0;
    float v = dist / length;
    return exp(-0.5 * v * v);
}
float gridpp::StructureFunction::cressman_rho(float dist, float length) const {
    if(!gridpp::is_valid(length) || length == 0)
        // Disabled
        return 1;
    if(!gridpp::is_valid(dist))
        return 0;
    if(dist >= length)
        return 0;
    return (length * length - dist * dist) / (length * length + dist * dist);
}
float gridpp::StructureFunction::localization_distance() const {
    return mLocalizationDistance;
}
gridpp::MultipleStructure::MultipleStructure(const StructureFunction& structure_h, const StructureFunction& structure_v, const StructureFunction& structure_w):
    gridpp::StructureFunction(structure_h.localization_distance()) {
    m_structure_h = structure_h.clone();
    m_structure_v = structure_v.clone();
    m_structure_w = structure_w.clone();
}
float gridpp::MultipleStructure::corr(const Point& p1, const Point& p2) const {
    Point p1_h(p1.lat, p1.lon, p1.elev, p1.laf, p1.type);
    Point p2_h(p2.lat, p2.lon, p1.elev, p1.laf, p1.type);
    Point p1_v(p1.lat, p1.lon, p1.elev, p1.laf, p1.type);
    Point p2_v(p1.lat, p1.lon, p2.elev, p1.laf, p1.type);
    Point p1_w(p1.lat, p1.lon, p1.elev, p1.laf, p1.type);
    Point p2_w(p1.lat, p1.lon, p1.elev, p2.laf, p1.type);
    float corr_h = m_structure_h->corr(p1_h, p2_h);
    float corr_v = m_structure_v->corr(p1_v, p2_v);
    float corr_w = m_structure_w->corr(p1_w, p2_w);
    return corr_h * corr_v * corr_w;
}
gridpp::StructureFunction* gridpp::MultipleStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::MultipleStructure(*m_structure_h, *m_structure_v, *m_structure_w);
    return val;
}

/** Barnes */
gridpp::BarnesStructure::BarnesStructure(float h, float v, float w, float hmax) :
    gridpp::StructureFunction(h) {
    if(gridpp::is_valid(hmax) && hmax < 0)
        throw std::invalid_argument("hmax must be >= 0");
    if(!gridpp::is_valid(v) || v < 0)
        throw std::invalid_argument("v must be >= 0");
    if(!gridpp::is_valid(w) || w < 0)
        throw std::invalid_argument("w must be >= 0");
    if(gridpp::is_valid(hmax))
        mLocalizationDistance = hmax;
    else {
        // Calculate the horizontal localization radius. For min_rho=0.0013, the factor is 3.64
        float default_min_rho = 0.0013;
        mLocalizationDistance = sqrt(-2*log(default_min_rho)) * h;
    }
    mH = h;
    mV = v;
    mW = w;
}
float gridpp::BarnesStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    if(hdist > localization_distance())
        return 0;
    float rho = gridpp::StructureFunction::barnes_rho(hdist, mH);
    if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
        float vdist = p1.elev - p2.elev;
        rho *= gridpp::StructureFunction::barnes_rho(vdist, mV);
    }
    if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
        float lafdist = p1.laf - p2.laf;
        rho *= gridpp::StructureFunction::barnes_rho(lafdist, mW);
    }
    return rho;
}
gridpp::StructureFunction* gridpp::BarnesStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::BarnesStructure(mH, mV, mW, mLocalizationDistance);
    return val;
}

/** Cressman */
gridpp::CressmanStructure::CressmanStructure(float h, float v, float w) :
    gridpp::StructureFunction(h) {
    if(!gridpp::is_valid(v) || v < 0)
        throw std::invalid_argument("v must be >= 0");
    if(!gridpp::is_valid(w) || w < 0)
        throw std::invalid_argument("w must be >= 0");
    mH = h;
    mV = v;
    mW = w;
}
float gridpp::CressmanStructure::corr(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    float rho = gridpp::StructureFunction::cressman_rho(hdist, mH);
    if(gridpp::is_valid(p1.elev) && gridpp::is_valid(p2.elev)) {
        float vdist = p1.elev - p2.elev;
        rho *= gridpp::StructureFunction::cressman_rho(vdist, mV);
    }
    if(gridpp::is_valid(p1.laf) && gridpp::is_valid(p2.laf)) {
        float lafdist = p1.laf - p2.laf;
        rho *= gridpp::StructureFunction::cressman_rho(lafdist, mW);
    }
    return rho;
}
gridpp::StructureFunction* gridpp::CressmanStructure::clone() const {
    gridpp::StructureFunction* val = new gridpp::CressmanStructure(mH, mV, mW);
    return val;
}

/** CrossValidation */
gridpp::CrossValidation::CrossValidation(StructureFunction& structure, float dist) :
        StructureFunction(0){
    if(!gridpp::is_valid(dist) || dist < 0)
        throw std::invalid_argument("Invalid 'dist' in CrossValidation structure");
    mLocalizationDistance = structure.localization_distance();
    m_structure = structure.clone();
    m_dist = dist;
}
float gridpp::CrossValidation::corr(const Point& p1, const Point& p2) const {
    return m_structure->corr(p1, p2);
}
float gridpp::CrossValidation::corr_background(const Point& p1, const Point& p2) const {
    float hdist = gridpp::KDTree::calc_distance_fast(p1, p2);
    if(gridpp::is_valid(m_dist)) {
        if(hdist <= m_dist)
            return 0;
    }
    return m_structure->corr_background(p1, p2);
}
gridpp::StructureFunction* gridpp::CrossValidation::clone() const {
    gridpp::StructureFunction* val = new gridpp::CrossValidation(*m_structure);
    return val;
}

