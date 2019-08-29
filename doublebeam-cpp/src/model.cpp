#include <algorithm>

#include "model.h"


VelocityModel::VelocityModel(const std::vector<Layer>& layers): layers(layers) {
    interface_depths.push_back(layers[0].top_depth);
    for (auto& l: layers){
        interface_depths.push_back(l.bot_depth);
    }
}

Layer VelocityModel::operator[](size_t index) const{
    return layers[index];
}

size_t VelocityModel::layer_index(double z) const{
    auto greater = std::upper_bound(interface_depths.begin(), interface_depths.end(), z);
    return std::min(std::distance(interface_depths.begin(), greater) - 1, static_cast<long int>(layers.size()) - 1);
}

double VelocityModel::eval_at(double z) const{
    auto layer = layers[layer_index(z)];
    return layer.gradient * z + layer.intercept;
}