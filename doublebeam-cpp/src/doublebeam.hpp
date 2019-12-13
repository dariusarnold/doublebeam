#ifndef DOUBLEBEAM_CPP_DOUBLEBEAM_H
#define DOUBLEBEAM_CPP_DOUBLEBEAM_H

#include <utility>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "model.hpp"
#include "raytracing.hpp"
#include "seismodata.hpp"
#include "twopoint.hpp"
#include "utils.hpp"


class FractureParameters {
public:
    FractureParameters(double depth, double phi_hat_x, double phi_hat_y,
                       int num_fracture_orientations, double spacing_min, double spacing_max,
                       int num_fracture_spacings);

    /**
     * Top depth of fracture plane.
     */
    double depth;
    /**
     * X component of central fracture orientation (unit vector perpendicular to fracture plane).
     */
    double phi_hat_x;
    /**
     * Y component of central fracture orientation (unit vector perpendicular to fracture plane).
     */
    double phi_hat_y;
    /**
     * Array of possible fracture orientations. Shape N, 2 where N is the number of fracture
     * orientations. Second axis contains x, y components of fracture orientations.
     * DB algorithm scans all possible fracture orientations and calculates the scatteres slowness.
     */
    std::vector<math::Vector2> orientations;
    /**
     * Array of possible fracture spacings.
     * DB algorithm scans over this array and calculates scattered slowness for every entry.
     */
    std::vector<double> spacings;
};

struct DoubleBeamResult {
    /**
     * Construct stacking amplitude sigma as a matrix with all results for one fixed fracture
     * spacing in one row and all results for one fixed fracture orientation in one column.
     * @param num_of_fracture_spacings Used as number of rows in the matrix.
     * @param num_of_fracture_orientations Used as number of columns in the matrix.
     */
    DoubleBeamResult(size_t num_of_fracture_spacings, size_t num_of_fracture_orientations);
    /**
     * Result of double beam algorithm: Matrix with stacking amplitude sigma.
     * TODO explain how fracture spacing/orientation is represented in this, eg. spacing values are
     *  column wise and orientations are row wise.
     */
    Eigen::ArrayXXcd data;
};

class DoubleBeam {
public:
    DoubleBeam(const VelocityModel& model);

    DoubleBeamResult algorithm(std::vector<position_t> source_geometry,
                               std::vector<position_t> target_geometry, SeismoData data,
                               FractureParameters fracture_info, Meter beam_width,
                               AngularFrequency beam_frequency, double window_length,
                               double max_stacking_distance);

private:
    const VelocityModel& model;
    TwoPointRayTracing twopoint;
    RayTracer tracer;
};

#endif // DOUBLEBEAM_CPP_DOUBLEBEAM_H
