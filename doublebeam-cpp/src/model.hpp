#ifndef DOUBLEBEAM_CPP_MODEL_HPP
#define DOUBLEBEAM_CPP_MODEL_HPP

#include <tuple>
#include <vector>

struct Layer {
    double top_depth;
    double bot_depth;
    double intercept;
    double gradient;
};

class VelocityModel {
    /**
     * Depth of all interfaces, including the top and the bottom one
     * where the layer ends.
     */
    std::vector<double> interface_depths;
    /**
     * Pair of above, below velocity of every interface (in m/s).
     * Includes interfaces at top and bottom where the model ends.
     */
    std::vector<std::pair<double, double>> _interface_velocities;
    /**
     * Velocity intercepts of all layers.
     */
    std::vector<double> intercepts;
    /**
     * Velocity gradient of all layers.
     */
    std::vector<double> gradients;
    /**
     * Array of layers.
     */
    std::vector<Layer> layers;

  public:
    explicit VelocityModel(const std::vector<Layer>& layers);

    /**
     * Get layer by index.
     * @param index Index of layer.
     * @return Layer at index.
     */
    Layer operator[](size_t index) const;

    /**
     * Get index of layer at depth z.
     * Top of layer is inclusive, bottom is exclusive. The bottom most
     * layer is a special case, since there the bottom of the layer also
     * belongs to the layer.
     * @param z Depth in m.
     * @return Index of layer in velocity model, 0 indexed.
     * @throw std::domain_error is thrown when depth is out of model range
     */
    size_t layer_index(double z) const;

    /**
     * Evaluate model at a certain depth.
     * @param z Depth in meter.
     * @return Velocity in m/s.
     */
    double eval_at(double z) const;

    /**
     * Return velocities above and below the closest interface.
     * 0 is returned for the velocity outside of the model, e.g. when the
     * interface between the top layer and the one above is requested.
     * @param z Depth in meter.
     * @return velocity at closest interface to depth z in m/s.
     */
    std::pair<double, double> interface_velocities(double z);

    /**
     * Return top and bottom depth of model in m.
     * @return Pair of (top, bottom).
     */
    std::pair<double, double> get_top_bottom();
};

#endif // DOUBLEBEAM_CPP_MODEL_HPP
