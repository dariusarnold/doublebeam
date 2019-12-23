#ifndef DOUBLEBEAM_CPP_MODEL_HPP
#define DOUBLEBEAM_CPP_MODEL_HPP

#include <filesystem>
#include <optional>
#include <tuple>
#include <vector>

#include "raytracing_types.hpp"
#include "units.hpp"


struct Layer {
    Meter top_depth;
    Meter bot_depth;
    Velocity velocity;
};

Meter layer_height(const Layer& layer);

bool operator==(const Layer& l1, const Layer& l2);


class VelocityModel {
    /**
     * Depth of all interfaces, including the top and the bottom one
     * where the layer ends.
     */
    std::vector<Meter> m_interface_depths;

    /**
     * Array of layers.
     */
    std::vector<Layer> layers_m;

    /**
     * Horizontal dimension of velocity model.
     */
    Meter x0_, x1_;
    Meter y0_, y1_;

public:
    /**
     * Create velocity model.
     * Model has a valid dimension from x0 to x1 along x axis and y0 to y1 along y axis.
     * @param layers Sequence of layers. First layer in sequence should be the top layer, with
     * below layers following.
     * @param x1 Second x coordinate of model width.
     * @param y1 Second y coordinate of model width.
     * @param x0 First x coordinate of model width.
     * @param y0 First y coordinate of model width.
     */
    explicit VelocityModel(const std::vector<Layer>& layers, Meter x1, Meter y1, Meter x0 = 0_meter,
                           Meter y0 = 0_meter);

    /**
     * Get layer by index.
     * @param index Index of layer.
     * @return Layer at index.
     */
    const Layer& operator[](size_t index) const;

    /**
     * Get index of layer at position x, y, z.
     * Top of layer is inclusive, bottom is exclusive. The bottom most
     * layer is a special case, since there the bottom of the layer also
     * belongs to the layer.
     * @param x X coordinate in m.
     * @param y Y coordinate in m.
     * @param z Depth in m.
     * @return Index of layer in velocity model, 0 indexed. Empty when point (x, y, z) is  outside
     * of velocity model.
     */
    [[nodiscard]] std::optional<size_t> layer_index(Meter x, Meter y, Meter z) const;

    /**
     * Get index of layer at depth z.
     * Top of layer is inclusive, bottom is exclusive. The bottom most
     * layer is a special case, since there the bottom of the layer also
     * belongs to the layer.
     * @param z Depth in m.
     * @return Index of layer in velocity model, 0 indexed. Empty when depth z is outside velocity
     * model.
     */
    [[nodiscard]] std::optional<std::size_t> layer_index(Meter z) const;

    /**
     * Evaluate model at a certain position.
     * @param x X coordinate.
     * @param y Y coordinate.
     * @param z Depth in meter.
     * @return Velocity in m/s.
     */
    std::optional<Velocity> eval_at(Meter x, Meter y, Meter z) const;

    /**
     * Evaluate model at position.
     * @param position
     * @return
     */
    std::optional<Velocity> eval_at(const Position& position) const;

    /**
     * Return true if point given by depth, z coordinates is within the velocity model.
     * @param z Depth of point in m.
     * @return
     */
    bool in_model(Meter x, Meter y, Meter z) const;

    struct InterfaceVelocities {
        Velocity above;
        Velocity below;
    };

    /**
     * Return velocities above and below the closest interface.
     * 0 is returned for the velocity outside of the model, e.g. when the
     * interface between the top layer and the one above is requested.
     * @param z Depth in meter.
     * @return velocity at closest interface to depth z in m/s.
     */
    [[nodiscard]] InterfaceVelocities interface_velocities(Meter z) const;

    /**
     * Return top and bottom depth of model in m.
     * @return Pair of (top, bottom).
     */
    std::pair<Meter, Meter> get_top_bottom() const;

    /**
     * Get x0, x1 pair.
     */
    std::pair<Meter, Meter> get_x_extent() const;

    /**
     * Get y0, y1 pair.
     */
    std::pair<Meter, Meter> get_y_extent() const;

    Meter get_x0() const;

    Meter get_x1() const;

    Meter get_y0() const;

    Meter get_y1() const;

    /**
     * Get width of model along x direction
     * @return
     */
    Meter x_width() const;

    /**
     * Get width of model along y direction
     * @return
     */
    Meter y_width() const;

    /**
     * Compare two velocity models. Return true if they are the same.
     * @param other Other VelocityModel to compare this instance to.
     */
    bool operator==(const VelocityModel& other) const;

    /**
     * Put model information (dimensions) into stream.
     * @param os
     * @param model
     * @return
     */
    friend std::ostream& operator<<(std::ostream& os, const VelocityModel& model);

    /**
     * Define begin and end to allow range-for iteration over layers.
     * @return
     */
    std::vector<Layer>::const_iterator begin() const;
    std::vector<Layer>::const_iterator end() const;

    /**
     * Return number of layers in model.
     * @return
     */
    std::size_t num_layers() const;

    const std::vector<Meter>& interface_depths() const;

    /**
     * Get Layer at depth z.
     */
    Layer get_layer(Meter x, Meter y, Meter z) const;

    /**
     * Get number of interfaces between two points.
     * Depths do not have to have any order.
     * @param z1
     * @param z2
     * @return
     */
    size_t number_of_interfaces_between(Meter z1, Meter z2) const;

    bool in_horizontal_extent(Meter x, Meter y) const;
};


/**
 * Create a velocity model from a file. For file format, see README.
 */
VelocityModel read_velocity_file(const std::filesystem::path& filepath);

/**
 * Find highest velocity between two points in the model.
 * Depths do not have to be sorted.
 * @param depth1 z coordinate of point 1.
 * @param depth2 z coordinate of point 2.
 * @param model Velocity model.
 * @return Highest velocity (m/s) between depth 1 and depth 2 for a direct ray.
 */
Velocity highest_velocity_between(Meter depth1, Meter depth2, const VelocityModel& model);

#endif // DOUBLEBEAM_CPP_MODEL_HPP
