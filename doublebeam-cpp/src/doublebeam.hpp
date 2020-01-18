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
    FractureParameters(math::Vector2 phi_hat, int num_fracture_orientations, Meter spacing_min,
                       Meter spacing_max, int num_fracture_spacings);
    /**
     * Central fracture orientation vector (unit vector perpendicular to fracture plane).
     * This vector is horizontal (has only x, y component) since only vertical fractures are
     * allowed.
     */
    math::Vector2 phi_hat;

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
    std::vector<Meter> spacings;
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
     * The matrix has num_of_fracture_orientations columns and num_of_fracture_spacings rows.
     * This means the results are stored as different orientations in columns and different spacings
     * in columns. The 0,0 index of the matrix corresponds to the first fracture orientation and the
     * first fracture spacing. Downwards the spacing index increases and to the right the
     * orientations index increases.
     */
    Eigen::ArrayXXcd data;
};

class DoubleBeam {
public:
    explicit DoubleBeam(const VelocityModel& model);

    /**
     * Performs doublebeam algorithm and calculates scattering coefficient sigma
     * @param source_geometry Vector of source positions
     * @param target Target position in the subsurface.
     * @param data Contains seismograms.
     * @param fracture_info Fracture information
     * @param beam_width Width of Gauss beam
     * @param beam_frequency Frequency for Gauss beam
     * @param window_length Window length for cutting seismic data. The real part of the added
     * complex traveltime from source and receiver beam gives the center point where the seismograms
     * are cut. The window will be centered around that value.
     * @param max_stacking_distance If a source or receiver beam is farther than this distance from
     * a source/receiver, it will not be evaluated in the stacking process.
     * @return
     */
    DoubleBeamResult algorithm(const std::vector<Position>& source_geometry, Position target,
                               const SeismoData& data, const FractureParameters& fracture_info,
                               Meter beam_width, AngularFrequency beam_frequency,
                               Second window_length, Meter max_stacking_distance);

private:
    Eigen::ArrayXXcd calc_sigma_for_sbc(const Position& source_beam_center, const Position& target,
                                        const FractureParameters& fracture_info,
                                        const SeismoData& data, Meter beam_width,
                                        AngularFrequency beam_frequency,
                                        const std::vector<WaveType>& ray_code, Second window_length,
                                        Meter max_stacking_distance);

    const VelocityModel& model;
    TwoPointRayTracing twopoint;
    RayTracer tracer;
};

#endif // DOUBLEBEAM_CPP_DOUBLEBEAM_H
