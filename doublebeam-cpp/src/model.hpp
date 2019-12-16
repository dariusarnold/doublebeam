#ifndef DOUBLEBEAM_CPP_MODEL_HPP
#define DOUBLEBEAM_CPP_MODEL_HPP

#include <filesystem>
#include <optional>
#include <tuple>
#include <vector>


struct Layer {
    double top_depth;
    double bot_depth;
    double velocity;
};

double layer_height(const Layer& layer);

bool operator==(const Layer& l1, const Layer& l2);


class VelocityModel {
    /**
     * Depth of all interfaces, including the top and the bottom one
     * where the layer ends.
     */
    std::vector<double> m_interface_depths;

    /**
     * Array of layers.
     */
    std::vector<Layer> layers;

    /**
     * Horizontal dimension of velocity model.
     */
    double x0_, x1_;
    double y0_, y1_;

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
    explicit VelocityModel(const std::vector<Layer>& layers, double x1, double y1, double x0 = 0,
                           double y0 = 0);

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
    [[nodiscard]] std::optional<size_t> layer_index(double x, double y, double z) const;

    /**
     * Get index of layer at depth z.
     * Top of layer is inclusive, bottom is exclusive. The bottom most
     * layer is a special case, since there the bottom of the layer also
     * belongs to the layer.
     * @param z Depth in m.
     * @return Index of layer in velocity model, 0 indexed. Empty when depth z is outside velocity
     * model.
     */
    [[nodiscard]] std::optional<std::size_t> layer_index(double z) const;

    /**
     * Evaluate model at a certain position.
     * @param x X coordinate.
     * @param y Y coordinate.
     * @param z Depth in meter.
     * @return Velocity in m/s.
     */
    std::optional<double> eval_at(double x, double y, double z) const;

    /**
     * Return true if point given by depth, z coordinates is within the velocity model.
     * @param z Depth of point in m.
     * @return
     */
    bool in_model(double x, double y, double z) const;

    struct InterfaceVelocities {
        double above;
        double below;
    };

    /**
     * Return velocities above and below the closest interface.
     * 0 is returned for the velocity outside of the model, e.g. when the
     * interface between the top layer and the one above is requested.
     * @param z Depth in meter.
     * @return velocity at closest interface to depth z in m/s.
     */
    InterfaceVelocities interface_velocities(double z) const;

    using iterator = std::vector<double>::const_iterator;
    /**
     * Return interface velocities between depths z1 and z2.
     * Interfaces have two velocities associated with them, one above and one below the interface.
     * The depths do not have to be sorted.
     * If a depth lies on an interface, that interface is not included in the returned range.
     * @param z1 First depth, has to be inside model.
     * @param z2 Second depth, has to be inside model.
     * @return Pair of iterators pointing to the first interface velocity and one past the last
     * interface velocity.
     */
    std::pair<iterator, iterator> interface_velocities(double z1, double z2) const;

    /**
     * Return top and bottom depth of model in m.
     * @return Pair of (top, bottom).
     */
    std::pair<double, double> get_top_bottom() const;

    /**
     * Get x0, x1 pair.
     */
    std::pair<double, double> get_x_extent() const;

    /**
     * Get y0, y1 pair.
     */
    std::pair<double, double> get_y_extent() const;

    double x0() const;

    double x1() const;

    double y0() const;

    double y1() const;

    /**
     * Get width of model along x direction
     * @return
     */
    double x_width() const;

    /**
     * Get width of model along y direction
     * @return
     */
    double y_width() const;

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

    const std::vector<double>& interface_depths() const;

    /**
     * Get Layer at depth z.
     */
    Layer get_layer(double x, double y, double z) const;

    /**
     * Get number of interfaces between two points.
     * Depths do not have to have any order.
     * @param z1
     * @param z2
     * @return
     */
    size_t number_of_interfaces_between(double z1, double z2) const;

    bool in_horizontal_extent(double x, double y) const;
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
double highest_velocity_between(double depth1, double depth2, const VelocityModel& model);

#endif // DOUBLEBEAM_CPP_MODEL_HPP
