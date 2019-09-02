#include <algorithm>

#include "model.h"


/**
 * Helper function to calculate velocity in layer.
 * @return Velocity (m/s) in layer at depth. No bounds checking is done.
 */
double layer_velocity(Layer layer, double depth) {
    return layer.intercept + layer.gradient * depth;
}


VelocityModel::VelocityModel(const std::vector<Layer>& layers) : layers(layers) {
    if (layers.empty()) {
        throw std::invalid_argument("Velocity model has to contain layers");
    }
    // First get interface depths
    auto top_layer = layers.front();
    interface_depths.push_back(top_layer.top_depth);
    std::for_each(layers.begin(), layers.end(), [&](const auto& layer){interface_depths.push_back(layer.bot_depth);});
    // Second get interface velocities with two special cases for first and last layer
    double vel_above = 0., vel_below;
    for (auto& l: layers) {
        vel_below = layer_velocity(l, l.top_depth);
        _interface_velocities.emplace_back(vel_above, vel_below);
        vel_above = layer_velocity(l, l.bot_depth);
    }
    _interface_velocities.emplace_back(vel_above, 0);
}

Layer VelocityModel::operator[](size_t index) const {
    return layers[index];
}

size_t VelocityModel::layer_index(double z) const {
    if (z < interface_depths.front() or z > interface_depths.back()){
        throw std::domain_error("Evaluating model outside of its depth range: " + std::to_string(z));
    }
    auto greater = std::upper_bound(interface_depths.begin(), interface_depths.end(), z);
    return std::min(std::distance(interface_depths.begin(), greater) - 1, static_cast<long int>(layers.size()) - 1);
}

double VelocityModel::eval_at(double z) const {
    Layer layer;
    try {
        layer = layers[layer_index(z)];
    } catch (const std::domain_error&) {
        return -1;
    }

    return layer.gradient * z + layer.intercept;
}

std::pair<double, double> VelocityModel::interface_velocities(double z) {
    auto index = layer_index(z);
    auto half_depth = layers[index].top_depth + 0.5 * (layers[index].bot_depth - layers[index].top_depth);
    if (z < half_depth) {
        return _interface_velocities[index];
    }
    return _interface_velocities[index + 1];
}

std::pair<double, double> VelocityModel::get_top_bottom() {
    return {interface_depths.front(), interface_depths.back()};
}
