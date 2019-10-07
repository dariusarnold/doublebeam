#ifndef DOUBLEBEAM_CPP_DOUBLEBEAM_H
#define DOUBLEBEAM_CPP_DOUBLEBEAM_H

#include <utility>

#include <xtensor/xtensor.hpp>

#include "model.hpp"
#include "raytracing.hpp"
#include "twopoint.hpp"
#include "utils.hpp"


/**
 * Generate coordinates of a evenly spaced grid.
 * @param x0 X start coordinate.
 * @param x1 X end coordinate.
 * @param y0 Y start coordinate.
 * @param y1 Y end coordinate.
 * @param depth Depth of all grid points.
 * @param num_x Number of points along the x axis.
 * @param num_y Number of points along the y axis.
 * @return Container of all grid locations with shape (num_x, num_y, 3) where the first axis
 * is the x axis, the second axis is the y axis and the third axis contains the z values.
 */
std::vector<position_t> grid_coordinates(double x0, double x1, double y0, double y1, double depth,
                                       int num_x, int num_y);


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

struct DoubleBeamParameters {};

class DoubleBeam {
public:
    DoubleBeam(const VelocityModel& model);

    void algorithm(std::vector<position_t> source_geometry, std::vector<position_t> target_geometry,
                   FractureParameters fracture_info, double beam_width, double beam_frequency,
                   double window_length);

private:
    const VelocityModel& model;
    TwoPointRayTracing twopoint;
    RayTracer tracer;
};

#endif // DOUBLEBEAM_CPP_DOUBLEBEAM_H