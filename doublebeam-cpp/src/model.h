#ifndef DOUBLEBEAM_CPP_MODEL_H
#define DOUBLEBEAM_CPP_MODEL_H

#include <vector>
#include <tuple>


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
     */
    size_t layer_index(double z) const;

    /**
     * Evaluate model at a certain depth.
     * @param z Depth in meter.
     * @return Velocity in m/s.
     */
    double eval_at(double z) const;
};


#endif //DOUBLEBEAM_CPP_MODEL_H
