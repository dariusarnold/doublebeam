struct Layer{
    double top_depth;
    double bot_depth;
    double intercept;
    double gradient;
};


class VelocityModel{
    std::vector<double> interface_depths;
    std::vector<double> intercepts;
    std::vector<double> gradients;
    std::vector<Layer> layers;
public:
    explicit VelocityModel(const std::vector<Layer>& layers): layers(layers) {
        interface_depths.push_back(layers[0].top_depth);
        for (auto& l: layers){
            interface_depths.push_back(l.bot_depth);
        }
    }

    Layer operator[](size_t index) const{
        return layers[index];
    }

    int64_t layer_index(double z) const{
        auto greater = std::upper_bound(interface_depths.begin(), interface_depths.end(), z);
        return std::min(std::distance(interface_depths.begin(), greater) - 1, static_cast<long int>(layers.size()) - 1);
    }

    double eval_at(double z) const{
        auto layer = layers[layer_index(z)];
        return layer.gradient * z + layer.intercept;
    }
};