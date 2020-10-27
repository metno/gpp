#ifndef GRIDPP_API_H
#define GRIDPP_API_H
#include <vector>
#include <string>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#ifdef _OPENMP
    #include <omp.h>
#endif

#define GRIDPP_VERSION "0.5.0a2"
#define __version__ GRIDPP_VERSION

namespace gridpp {
    /** **************************************
     * @name Short-hand notation for vectors of different dimensions sizes
     * ***************************************/ /**@{*/
    // Preferred vector types
    typedef std::vector<float> vec;
    typedef std::vector<vec> vec2;
    typedef std::vector<vec2> vec3;
    typedef std::vector<int> ivec;
    typedef std::vector<ivec> ivec2;
    typedef std::vector<ivec2> ivec3;

    // Only use when double is required
    typedef std::vector<double> dvec;
    typedef std::vector<dvec> dvec2;
    /**@}*/

    /** **************************************
     * @name Constants
     * Functions that assimilate observations onto a gridded background
     * ***************************************/ /**@{*/
    /** Missing value indicator */
    static const float MV = NAN;
    /** Missing value indicator in gridpp command-line tool */
    static const float MV_CML = -999;
    /** Mathematical constant pi */
    static const float pi = 3.14159265;
    /** Radius of the earth [m] */
    static const double radius_earth = 6.378137e6;
    /**@}*/

    class KDTree;
    class Points;
    class Grid;
    class Point;
    class Nearest;
    class StructureFunction;
    class Transform;

    /** Methods for extrapolating outside a curve */
    enum Extrapolation {
            OneToOne = 0,      /**< Continue past the end-points using a slope of 1 */
            MeanSlope = 10,    /**< Continue past the end-points using the mean slope of the curve*/
            NearestSlope = 20, /**< Continue past the end-points using the slope of the two lowermost or uppermost points in the curve */
            Zero = 30,         /**< Continue past the end-points using a slope of 0 */
        };

    /** Statistical operations to reduce a vector to a scalar */
    enum Statistic {
        Mean      = 0,  /**< Mean of values */
        Min       = 10, /**< Minimum of values */
        Median    = 20, /**< Mean of values */
        Max       = 30, /**< Maximum of values */
        Quantile  = 40, /**< A quantile from values */
        Std       = 50, /**< Standard deviation of values */
        Variance  = 60, /**< Population variance of values */
        Sum       = 70, /**< Sum of values */
        Unknown   = -1  /**< Unknown statistic */
    };

    /** Binary verification metrics */
    enum Metric {
        Ets       = 0,  /**< Equitable threat score */
        Ts        = 1,  /**< Threat score */
        Kss       = 20, /**< Hannsen-Kuiper skill score */
        Pc        = 30, /**< Proportion correct */
        Bias      = 40, /**< Bias */
        Hss       = 50, /**< Heidke skill score */
    };

    /** Method for statistical correction */
    enum CorrectionType {
        Qq        = 0,        /**< Quantile mapping */
        Multiplicative = 10,  /**< Multiplicative */
        Additive  = 20,       /**< Additive */
    };

    /** Types of coordinates for position of points */
    enum CoordinateType {
        Geodetic = 0,      /**< Latitude and longitude */
        Cartesian = 1,     /**< X and Y */
    };

    /** **************************************
     * @name Data assimilation methods
     * Functions that merge observations with a background field
     * ***************************************/ /**@{*/

    /** Optimal interpolation for a deterministic gridded field
      * @param bgrid Grid of background field
      * @param background 2D field of background values
      * @param points Points of observations
      * @param pobs Vector of observations
      * @param pratios Vector of ratio of observation error variance to background variance
      * @param pbackground Background with observation operator
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param cross_validation_distance Exclude observations within this distance of the gridpoint [m]
    */
    vec2 optimal_interpolation(const Grid& bgrid,
            const vec2& background,
            const Points& points,
            const vec& pobs,
            const vec& pratios,
            const vec& pbackground,
            const StructureFunction& structure,
            int max_points,
            float cross_validation_distance=MV);

    /** Optimal interpolation for a deterministic vector of points
      * @param bpoints Points of background field
      * @param background 1D field of background values
      * @param points Points of observations
      * @param pobs Vector of observations
      * @param pratios Vector of ratio of observation error variance to background variance
      * @param pbackground Background with observation operator
      * @param structure Structure function
      * @param max_points Maximum number of observations to use inside localization zone; Use 0 to disable
      * @param cross_validation_distance Exclude observations within this distance of the gridpoint [m]
    */
    vec optimal_interpolation(const Points& bpoints,
            const vec& background,
            const Points& points,
            const vec& pobs,
            const vec& pratios,
            const vec& pbackground,
            const StructureFunction& structure,
            int max_points,
            float cross_validation_distance=MV);

    /** This is its own function because the variance parameterization is different than in the
      * non-transformed case */
    vec2 optimal_interpolation_transform(const Grid& bgrid,
            const vec2& background,
            float bsigma,
            const Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec& pbackground,
            const StructureFunction& structure,
            int max_points,
            const Transform& transform);

    vec optimal_interpolation_transform(const Points& bpoints,
            const vec& background,
            float bsigma,
            const Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec& pbackground,
            const StructureFunction& structure,
            int max_points,
            const Transform& transform,
            float cross_validation_distance=MV);

    /** Optimal interpolation using a structure function based on an ensemble 
      * See Lussana et al 2019 (DOI: 10.1002/qj.3646)
      * @param input 3D field of background values (Y, X, E)
      * @param bgrid grid corresponding to input
      * @param pobs vector of observations
      * @param pci vector of ci values
      * @param points observation points
    */
    vec3 optimal_interpolation_ensi(const Grid& bgrid,
            const vec3& background,
            const Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec2& pbackground,
            const StructureFunction& structure,
            int max_points,
            float cross_validation_distance=MV);

    vec2 optimal_interpolation_ensi(const Points& bpoints,
            const vec2& background,
            const Points& points,
            const vec& pobs,
            const vec& psigmas,
            const vec2& pbackground,
            const StructureFunction& structure,
            int max_points,
            float cross_validation_distance=MV);

    /** Fill in values inside or outside a set of circles
      * @param input Deterministic values with dimensions Y, X
      * @param radii Circle radii for each point
      * @param value Fill in this value
      * @param outside if True, fill outside circles, if False, fill inside circles
    */
    vec2 fill(const Grid& igrid, const vec2& input, const Points& points, const vec& radii, float value, bool outside);
    /**@}*/

    /** **************************************
     * @name Spatial neighbourhood filters
     * Functions that apply neighbourhood filters on a gridded field
     * ***************************************/ /**@{*/

    /** Spatial neighbourhood filter, computing a statistic for a sliding square window
      * @param input 2D grid of values
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
    */
    vec2 neighbourhood(const vec2& input, int halfwidth, Statistic statistic);

    /** Spatial neighbourhood filter for an ensemble of fields
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param statistic Statistic to compute
    */
    vec2 neighbourhood(const vec3& input, int halfwidth, Statistic statistic);

    /** Computes a quantile in a sliding square neighbourhood
      * @param input 2D grid of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
    */
    vec2 neighbourhood_quantile(const vec2& input, float quantile, int halfwidth);

    /** Computes a quantile in a sliding square neighbourhood for an ensemble of fields
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
    */
    vec2 neighbourhood_quantile(const vec3& input, float quantile, int halfwidth);

    /** Fast and approximate neighbourhood quantile
      * @param input 2D grid of values
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, float quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile for ensemble of fields
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param quantile Quantile to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec3& input, float quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile, with spatially varying quantile
      * @param input 2D grid of values
      * @param quantile 2D grid quantiles to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec2& input, const vec2& quantile, int halfwidth, const vec& thresholds);

    /** Fast and approximate neighbourhood quantile for ensemble of fields, with spatially varying quantile
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param quantile 2D grid quantiles to compute (between 0 and 1)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param thresholds Vector of thresholds to use to approximate value
    */
    vec2 neighbourhood_quantile_fast(const vec3& input, const vec2& quantile, int halfwidth, const vec& thresholds);

    /** Spatial neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
      * @param input 2D grid of values
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param operation one of min, mean, median, max
    */
    vec2 neighbourhood_brute_force(const vec2& input, int halfwidth, Statistic statistic);

    /** Spatial neighbourhood filter without any shortcuts. This is quite slow and is only useful for testing.
      * @param input 3D grid of values with dimensions (Y, X, E)
      * @param halfwidth Filter halfwidth in number of gridpoints
      * @param operation one of min, mean, median, max
    */
    vec2 neighbourhood_brute_force(const vec3& input, int halfwidth, Statistic statistic);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile
      * @param input 2D (Y, X) array of values
      * @param num_thresholds Number of thresholds
    */
    vec get_neighbourhood_thresholds(const vec2& input, int num_thresholds);

    /** Calculate appropriate approximation thresholds for neighbourhood quantile based on an * ensemble
      * @param input 3D (Y, X, T) array of values
      * @param num_thresholds Number of thresholds
    */
    vec get_neighbourhood_thresholds(const vec3& input, int num_thresholds);

    /** Deprecated: Compute neighbourhood statistic on ensemble field
      * @deprecated Use neighbourhood() function */
    vec2 neighbourhood_ens(const vec3& input, int halfwidth, Statistic statistic);
    /** Deprecated: Compute neighbourhood quantiles on ensemble field
      * @deprecated Use neighbourhood_quantile() function */
    vec2 neighbourhood_quantile_ens(const vec3& input, float quantile, int halfwidth);
    /** Deprecated: Compute neighbourhood quantiles fast on ensemble field
      * @deprecated Use neighbourhood_quantile_fast() function */
    vec2 neighbourhood_quantile_ens_fast(const vec3& input, float quantile, int radius, const vec& thresholds);
    /**@}*/

    /** **************************************
     * @name Calibration methods              
     * Functions that apply statistical methods to data
     * ***************************************/ /**@{*/

    /** Create quantile mapping calibration curve
     *  @param ref Reference values (observations)
     *  @param fcst Forecast values
     *  @param quantiles Vector of quantiles to extract. If empty, use all values.
     *  @return Calibration curve
    */
    vec2 quantile_mapping_curve(const vec& ref, const vec& fcst, vec quantiles=vec());

    /** Create calibration curve that optimizes a metric
     *  @param ref Reference values (observations)
     *  @param fcst Forecast values
     *  @param thresholds Thresholds for computing optimal values for
     *  @param metric A Metric to optimize for
     *  @return Calibration curve
    */

    vec2 metric_optimizer_curve(const vec& ref, const vec& fcst, const vec& thresholds, Metric metric);
    /** Apply arbitrary calibration curve to 1D forecasts
     *  @param fcst 1D vector of forecast values
     *  @param curve Calibration curve
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @return Calibrated forecasts
    */

    vec apply_curve(const vec& fcst, const vec2& curve, Extrapolation policy_below, Extrapolation policy_above);

    /** Apply arbitrary calibration curve to 2D forecasts
     *  @param fcst 2D grid of forecast values
     *  @param curve Calibration curve
     *  @param policy_below Extrapolation policy below curve
     *  @param policy_above Extrapolation policy above curve
     *  @return Calibrated forecasts
    */
    vec2 apply_curve(const vec2& fcst, const vec2& curve, Extrapolation policy_below, Extrapolation policy_above);

    /** Ensure calibration curve is monotonic, by removing points
     *  @param curve Calibration curve
     *  @returns Monotonic calibration curve
    */
    vec2 monotonize_curve(vec2 curve);

    float get_optimal_threshold(const vec& ref, const vec& fcst, float threshold, Metric metric);

    float calc_score(float a, float b, float c, float d, Metric metric);
    float calc_score(const vec& ref, const vec& fcst, float threshold, Metric metric);
    float calc_score(const vec& ref, const vec& fcst, float threshold, float fthreshold, Metric metric);

    void advection_implicit(const vec2& y_dist, const vec2& x_dist, float dt, ivec2& y_coord, ivec2& x_coord);

    vec2 correction(const Grid& rgrid, const vec2& rvalues, const Points& npoints, const vec& nvalues, float mean_radius, float outer_radius, float inner_radius, int min_num, int max_num, CorrectionType type, ivec2& count);
    // Apply correction based on multiple timesteps
    vec2 correction(const Grid& rgrid, const vec3& rvalues, const Points& npoints, const vec2& nvalues, const vec2& apply_values, float mean_radius, float outer_radius, float inner_radius, int min_num, int max_num, CorrectionType type, ivec2& count);

    /**@}*/

    /** **************************************
     * @name Downscaling methods
     * Functions that interpolate data from one grid to another
     * ***************************************/ /**@{*/

    /** Nearest neighbour dowscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values on the output grid
    */
    vec2 nearest(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);
    vec3 nearest(const Grid& igrid, const Grid& ogrid, const vec3 ivalues);

    /** Nearest neighbour dowscaling grid to point
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec nearest(const Grid& igrid, const Points& opoints, const vec2 ivalues);
    vec2 nearest(const Grid& igrid, const Points& opoints, const vec3 ivalues);

    /** Bilinear downscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output grid to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values on the output grid
    */
    vec2 bilinear(const Grid& igrid, const Grid& ogrid, const vec2 ivalues);

    /** Bilinear downscaling grid to points
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @return Values for the output points
    */
    vec bilinear(const Grid& igrid, const Points& opoints, const vec2 ivalues);

    vec2 simple_gradient(const Grid& igrid, const Grid& ogrid, const vec2 ivalues, float elev_gradient);
    vec simple_gradient(const Grid& igrid, const Points& opoints, const vec2 ivalues, float elev_gradient);

    /** Smart neighbour downscaling grid to grid
      * @param igrid Input grid
      * @param ogrid Output points to downscale to
      * @param ivalues 2D vector of values on the input grid
      * @param num Number of neighbours to average
      * @param structure Structure function for determining similarity
      * @return Values for the output points
    */
    vec2 smart(const Grid& igrid, const Grid& ogrid, const vec2& ivalues, int num, const StructureFunction& structure);
    /**@}*/

    /** **************************************
     * @name Grid calculations
     * Functions that calculate statistics on a grid
     * ***************************************/ /**@{*/

    /** For each point, counts the number of gridpoints within the radius
     *  @param grid Grid
     *  @param points Points
     *  @param radius Radius [m]
     *  @return Number of gridpoints
    */
    vec count(const Grid& grid, const Points& points, float radius);

    /** For each gridpoint, counts the number of points within the radius
     *  @param grid Grid
     *  @param points Points
     *  @param radius Radius [m]
     *  @return Number of points
    */
    vec2 count(const Points& points, const Grid& grid, float radius);

    /** For each point, calculates the distance to nearest gridpoint
     *  @param grid Grid
     *  @param points Points
     *  @param num Number of points
     *  @return Distance [m] to nearest gridpoint for each point
    */
    vec distance(const Grid& grid, const Points& points, int num=1);

    /** For each output gridpoint, calculate the distance to nearest input gridpoint
     *  @param grid Grid
     *  @param ogrid Output grid
     *  @param num Number of points
     *  @return Distance [m] to nearest gridpoint for each gridpoint
    */
    vec2 distance(const Grid& igrid, const Grid& ogrid, int num=1);

    /** For each gridpoint, calculates the distance to nearest point
     *  @param points Points
     *  @param grid Grid
     *  @param num Number of points
     *  @return Distance [m] to nearest gridpoint for each point
    */
    vec2 distance(const Points& points, const Grid& grid, int num=1);

    /** Fill in missing values based on nearby values
      * @param values 2D array of values
      * @return 2D array of values without any missing values
    */
    vec2 fill_missing(const vec2& values);

    /** Aggregate points onto a grid. Writes MV where there are not enough observations
      * @param grid Grid to aggregate to
      * @param points Points with values
      * @param values Values at points
      * @param radius Circle radius for aggregate points [m]
      * @param min_num Minimum number of points in radius to create a value
      * @param statistic Statistic to compute on points within radius
    */
    vec2 gridding(const Grid& grid, const Points& points, const vec& values, float radius, int min_num, Statistic statistic);

    /**@}*/

    /** ****************************************
     * @name Diagnosing meteorological variables
     * Functions that diagnose a meteorological variable based on other variables
     * *****************************************/ /**@{*/

    /** Calculate dewpoint temperature from temperature and relative humidity
     *  @param temperature Temperature [K]
     *  @param relative_humidity Relative humidity [1]
     *  @returns Dewpoint temperature [K]
    */
    float dewpoint(float temperature, float relative_humidity);

    /** Vector version of dewpoint calculation
     *  @param temperature Temperatures [K]
     *  @param relative_humidity Relative humidities [1]
     *  @returns Dewpoint temperatures [K]
    */
    vec dewpoint(const vec& temperature, const vec& relative_humidity);

    /** Calculate pressure at a new elevation
     *  @param ielev Elevation at start point
     *  @param oelev Elevation at new point
     *  @param ipressure Pressure at start point
     *  @param itemperature Temperature at start point
     *  @return Pressure at new point
     */
    float pressure(float ielev, float oelev, float ipressure, float itemperature=288.15);

    /** Calculate Vector version of pressure calculation
     *  @param ielev Elevations at start point
     *  @param oelev Elevations at new point
     *  @param ipressure Pressures at start point
     *  @param itemperature Temperatures at start point
     *  @return Pressures at new points
     */
    vec pressure(const vec& ielev, const vec& oelev, const vec& ipressure, const vec& itemperature);

    /** Diagnose QNH from pressure and altitude
     *  @param pressure Pressure at point [pa]
     *  @param altitude Altitude of point [m]
     *  @returns QNH [pa]
    */
    float qnh(float pressure, float altitude);

    /** Vector version of QNH calculation
     *  @param pressure Pressures at points [pa]
     *  @param altitude Altitudes of points [m]
     *  @returns QNH [pa]
    */
    vec qnh(const vec& pressure, const vec& altitude);

    /** Calculate relative humidity from temperature and dewpoint temperature
     *  @param temperature Temperature [K]
     *  @param dewpoint Dewpoint temperature [K]
     *  @returns Relative humidity [1]
    */
    float relative_humidity(float temperature, float dewpoint);

    /** Vector version of relative humidity calculation
     *  @param temperature Temperatures [K]
     *  @param dewpoint Dewpoint temperatures [K]
     *  @returns Relative humidities [1]
    */
    vec relative_humidity(const vec& temperature, const vec& dewpoint);

    /** Calculate wetbulb temperature from temperature, pressure, and relative humidity
     *  @param temperature Temperature [K]
     *  @param pressure Air pressure [pa]
     *  @param Relative humidity [1]
     *  @returns Wetbulb temperature [K]
    */
    float wetbulb(float temperature, float pressure, float relative_humidity);

    /** Vector version of wetbulb calculation
     *  @param temperature Temperatures [K]
     *  @param pressure Air pressures [pa]
     *  @param Relative humidities [1]
     *  @returns Wetbulb temperatures [K]
    */
    vec wetbulb(const vec& temperature, const vec& pressure, const vec& relative_humidity);

    /** Diagnose wind speed from its components
     *  @param xwind X-component of wind [any unit]
     *  @param ywind Y-component of wind [any unit]
     *  @return Wind speed [any unit]
     * */
    float wind_speed(float xwind, float ywind);

    /** Vector version of wind speed calculation
     *  @param xwind X-components of wind [any unit]
     *  @param ywind Y-components of wind [any unit]
     *  @return Wind speeds [any unit]
     * */
    vec wind_speed(const vec& xwind, const vec& ywind);

    /**@}*/

    /** ****************************************
     * @name OpenMP settings
     * Functions that configure OpenMP
     * *****************************************/ /**@{*/
    /** Set the number of OpenMP threads to use. Overwrides OMP_NUM_THREAD env variable. */
    void set_omp_threads(int num);

    /** Sets the number of OpenMP threads to 1 if OMP_NUM_THREADS undefined */
    void initialize_omp();

    /** ****************************************
     * @name Utilities
     * Helper functions
     * *****************************************/ /**@{*/
    // vec2 calc_gradient(const vec2& values, const vec2& aux, int radius);
    // ivec regression(const vec& x, const vec& y);

    /** Convert name of a statistic enum */
    Statistic get_statistic(std::string name);

    /** The gridpp version
     * @return The gridpp version
    */
    std::string version();

    double clock();
    void debug(std::string string);
    void warning(std::string string);
    void error(std::string string);
    void future_deprecation_warning(std::string function, std::string other="");
    bool is_valid(float value);
    float calc_statistic(const vec& array, Statistic statistic);
    float calc_quantile(const vec& array, float quantile);
    vec calc_statistic(const vec2& array, Statistic statistic);
    vec calc_quantile(const vec2& array, float quantile=MV);
    int num_missing_values(const vec2& iArray);

    /** Find the index in a vector that is equal or just below a value
     *  @param iX Lookup value
     *  @param iValues Lookup vector. Must be sorted.
     *  @return The index into iValues that is equal or just below iX
    */
    int get_lower_index(float iX, const std::vector<float>& iValues);

    /** Find the index in a vector that is equal or just above a value
     *  @param iX Lookup value
     *  @param iValues Lookup vector. Must be sorted.
     *  @return The index into iValues that is equal or just above iX
    */
    int get_upper_index(float iX, const std::vector<float>& iValues);

    /** Piecewise linear interpolation
     *  If x is outside the range of iX, then the min/max value of iY is used
     *  @param x Interpolation to this point
     *  @param iX Vector of x-axis values. Vector must be sorted.
     *  @param iY Vector of y-axis values corresponding to iX.
     *  @return Y value corresponding to x
    */
    float interpolate(float x, const std::vector<float>& iX, const std::vector<float>& iY);
    void not_implemented_error();

    /** Initialize a vector of size Y, X, with a given value */
    ivec2 init_ivec2(int Y, int X, int value);
    vec2 init_vec2(int Y, int X, float value=MV);

    /** Initialize a vector of size Y, X, E, with a given value */
    ivec3 init_ivec3(int Y, int X, int E, int value);
    vec3 init_vec3(int Y, int X, int E, float value=MV);

    /** Get reasonably spaced quantiles from a vector of values, ignoring duplicate values
      *  but including the first number after duplicated values. Include the lowest and highest
      *  values.
      *  @param values vector of values (unsorted, and no invalid values)
      *  @param num number of thresholds to get
    */
    vec calc_even_quantiles(const vec& values, int num);

    /** Computes gradients based on values in neighbourhood
     *  @param grid Grid
     *  @param base Dependent variable. Missing values are not used.
     *  @param values Independent variable. Missing values are not used.
     *  @param radius Neighbourhood radius in number of gridpoints
     *  @param min_nim Minimum number of points required to compute gradient
     *  @param min_range Minimum range of base to compute gradient
     *  @param default_gradient Use this gradient if minimum number is not met
    */
    vec2 calc_gradient(const Grid& grid, const vec2& base, const vec2& values, int radius, int min_num=2, float min_range=0, float default_gradient=0);

    /** Check if the grid is the same size as the 2D vector */
    bool compatible_size(const Grid& grid, const vec2& v);
    bool compatible_size(const Grid& grid, const vec3& v);

    /** Checks if a point is located inside a rectangle formed by 4 points. The 4 points must be
      * provided in an order that draws out a rectangle (either clockwise or counter-clockwise)
      * @param A: A point in the rectangle
      * @param B: A point in the rectangle
      * @param C: A point in the rectangle
      * @param D: A point in the rectangle
      * @param m: The point to test if it is inside
      * @return True if the point is inside, False otherwise
    */
    bool point_in_rectangle(const Point& A, const Point& B, const Point& C, const Point& D, const Point& m );

    /** ****************************************
     * @name SWIG testing functions
     * Functions for testing the SWIG interface. Not useful for any other purpose.
     * *****************************************/ /**@{*/
    /** Special function whose presense is needed for SWIG */
    float* test_array(float* v, int n);
    /** Testing function for 1D input vector */
    float test_vec_input(const vec& input);
    /** Testing function for 1D input vector */
    int test_ivec_input(const ivec& input);
    /** Testing function for 2D input vector */
    float test_vec2_input(const vec2& input);
    /** Testing function for 3D input vector */
    float test_vec3_input(const vec3& input);
    /** Testing function for 1D output vector */
    vec test_vec_output();
    ivec test_ivec_output();
    /** Testing function for 2D output vector */
    vec2 test_vec2_output();
    ivec2 test_ivec2_output();
    /** Testing function for 3D output vector */
    vec3 test_vec3_output();
    ivec3 test_ivec3_output();

    /** Testing function for 1D vector treated as output */
    float test_vec_argout(vec& distances);
    /** Testing function for 2D vector treated as output */
    float test_vec2_argout(vec2& distances);

    /** Default value used to fill array in SWIG testing functions. Not useful for any other purpose. */
    static const float swig_default_value = -1;

    /**@}*/

    /** Represents a single point in some coordinate system */
    class Point {
        public:
            /** Constructor
              * @param lat: Latitude coordinate
              * @param lon: Longitude coordinate
              * @param elev: Elevation
              * @param laf: Land area fraction (between 0 and 1)
              * @param type: Coordinate type for lat and lon
            */
            Point(float lat, float lon, float elev=MV, float laf=MV, CoordinateType type=Geodetic);
            float lat;
            float lon;
            float elev;
            float laf;
            CoordinateType type;
    };
    /** Covariance structure function */
    class StructureFunction {
        public:
            StructureFunction(float localization_distance);
            /** Correlation between two points */
            virtual float corr(const Point& p1, const Point& p2) const = 0;
            virtual float corr_background(const Point& p1, const Point& p2) const;
            /** Maximum distance for which an observation can have an impact (localization)
              * @return Distance [m]
            */
            float localization_distance() const;
            virtual StructureFunction* clone() const = 0;
        protected:
            float barnes_rho(float dist, float length) const;
            float cressman_rho(float dist, float length) const;
            float mLocalizationDistance;
    };
    /** Simple structure function based on distance, elevation, and land area fraction */
    class BarnesStructure: public StructureFunction {
        public:
            /** Exponential structure function
              * @param h: Horizontal decorrelation length [m]
              * @param v: Vertical decorrelation length [m]. If 0, disable decorrelation.
              * @param w: Land/sea decorrelation length [1]. If 0, disable decorrelation.
              * @param hmax: Truncate horizontal correlation at this length [m]. If undefined, 3.64 * h.
            */
            BarnesStructure(float h, float v=0, float w=0, float hmax=MV);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
        private:
            float mH;
            float mV;
            float mW;
    };

    /** Simple structure function based on distance, elevation, and land area fraction */
    class CressmanStructure: public StructureFunction {
        public:
            CressmanStructure(float h, float v=0, float w=0);
            float corr(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
        private:
            float mH;
            float mV;
            float mW;
    };
    class CrossValidation: public StructureFunction {
        public:
            CrossValidation(StructureFunction& structure, float dist=0);
            float corr(const Point& p1, const Point& p2) const;
            float corr_background(const Point& p1, const Point& p2) const;
            StructureFunction* clone() const;
        private:
            StructureFunction* m_structure;
            float m_dist;
    };

    class Transform {
        public:
            virtual float forward(float value) const = 0;
            virtual float backward(float value) const = 0;
            virtual vec2 forward(const vec2& input) const;
            virtual vec2 backward(const vec2& input) const;
    };
    class Identity : public Transform {
        public:
            float forward(float value) const;
            float backward(float value) const;
    };
    class Log : public Transform {
        public:
            float forward(float value) const;
            float backward(float value) const;
    };
    class BoxCox : public Transform {
        public:
            BoxCox(float threshold);
            float forward(float value) const;
            float backward(float value) const;
        private:
            float mThreshold;
    };

    /** Helper class for Grid and Points */
    class KDTree {
        public:
            KDTree(vec lats, vec lons, CoordinateType type=Geodetic);
            KDTree& operator=(KDTree other);
            KDTree(const KDTree& other);
            KDTree() {};

            /** Find single nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             * */
            int get_nearest_neighbour(float lat, float lon) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            ivec get_neighbours(float lat, float lon, float radius) const;

            /** Find all points with a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             *  @param distances Vector to store separation distances [m]
             * */
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const;

            /** Find the number of points within a radius
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param radius Lookup radius [m]
             * */
            int get_num_neighbours(float lat, float lon, float radius) const;

            /** Find a set of nearest points
             *  @param lat Latitude of lookup-point
             *  @param lon Longitude of lookup-point
             *  @param num Number of points to find
             * */
            ivec get_closest_neighbours(float lat, float lon, int num) const;


            /** Convert lat/lons to 3D cartesian coordinates with the centre of the earth as the origin
             *  @param lats vector of latitudes [deg]
             *  @param lons vector of longitudes [deg]
             *  @param x_coords vector of x-coordinates [m]
             *  @param y_coords vector of y-coordinates [m]
             *  @param z_coords vector of z-coordinates [m]
             * */
            bool convert_coordinates(const vec& lats, const vec& lons, vec& x_coords, vec& y_coords, vec& z_coords) const;

            /** Same as above, but convert a single lat/lon to 3D cartesian coordinates
             *  @param lat latitude [deg]
             *  @param lon longitude [deg]
             *  @param x_coord x-coordinate [m]
             *  @param y_coord y-coordinate [m]
             *  @param z_coord z-coordinate [m]
             * */
            bool convert_coordinates(float lat, float lon, float& x_coord, float& y_coord, float& z_coord) const;
            static float deg2rad(float deg);
            static float rad2deg(float deg);
            static float calc_distance(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);
            static float calc_distance(float x0, float y0, float z0, float x1, float y1, float z1);
            static float calc_distance_fast(float lat1, float lon1, float lat2, float lon2, CoordinateType type=Geodetic);
            static float calc_distance_fast(const Point& p1, const Point& p2);
            vec get_lats() const;
            vec get_lons() const;
            int size() const;
            CoordinateType get_coordinate_type() const;
        protected:
            typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
            typedef std::pair<point, unsigned> value;
            typedef boost::geometry::model::box<point> box;
            boost::geometry::index::rtree< value, boost::geometry::index::quadratic<16> > mTree;
            vec mLats;
            vec mLons;
            CoordinateType mType;
    };

    /** Represents a vector of locations and their metadata */
    class Points  {
        public:
            Points();
            /** Initialize a new grid
             *  @param lats: vector of latitudes [degrees]
             *  @param lons: vector of longitudes [degrees]
             *  @param elevs: vector of elevations [m]
             *  @param lafs: vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Points(vec lats, vec lons, vec elevs=vec(), vec lafs=vec(), CoordinateType type=Geodetic);
            Points(KDTree tree, vec elevs=vec(), vec lafs=vec());
            Points& operator=(Points other);
            Points(const Points& other);
            // Returns -1 if there are no neighbours
            int get_nearest_neighbour(float lat, float lon) const;
            ivec get_neighbours(float lat, float lon, float radius) const;
            ivec get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const;
            int get_num_neighbours(float lat, float lon, float radius) const;
            ivec get_closest_neighbours(float lat, float lon, int num) const;

            vec get_lats() const;
            vec get_lons() const;
            vec get_elevs() const;
            vec get_lafs() const;
            int size() const;
            ivec get_in_domain_indices(const Grid& grid) const;
            Points get_in_domain(const Grid& grid) const;
            CoordinateType get_coordinate_type() const;
        private:
            KDTree mTree;
            vec mLats;
            vec mLons;
            vec mElevs;
            vec mLafs;
    };

    /** Represents a 2D grid of locations and their metadata */
    class Grid {
        public:
            Grid();

            /** Initialize a new grid
             *  @param lats: 2D vector of latitudes [degrees]
             *  @param lons: 2D vector of longitudes [degrees]
             *  @param elevs: 2D vector of elevations [m]
             *  @param lafs: 2D vector of land area fractions [1]
             *  @param type: Coordinate type
            */
            Grid(vec2 lats, vec2 lons, vec2 elevs=vec2(), vec2 lafs=vec2(), CoordinateType type=Geodetic);
            ivec get_nearest_neighbour(float lat, float lon) const;
            ivec2 get_neighbours(float lat, float lon, float radius) const;
            ivec2 get_neighbours_with_distance(float lat, float lon, float radius, vec& distances) const;
            int get_num_neighbours(float lat, float lon, float radius) const;
            ivec2 get_closest_neighbours(float lat, float lon, int num) const;

            bool get_box(float lat, float lon, int& Y1_out, int& X1_out, int& Y2_out, int& X2_out) const;

            /** Convert grid to a vector of points */
            Points to_points() const;

            vec2 get_lats() const;
            vec2 get_lons() const;
            vec2 get_elevs() const;
            vec2 get_lafs() const;
            ivec size() const;
            CoordinateType get_coordinate_type() const;
        private:
            KDTree mTree;
            int mX;
            vec2 get_2d(vec input) const;
            ivec get_indices(int index) const;
            ivec2 get_indices(ivec indices) const;
            vec2 mLats;
            vec2 mLons;
            vec2 mElevs;
            vec2 mLafs;
    };
};
#endif
