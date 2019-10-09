#include <algorithm>
#include <fstream>

#include "model.hpp"
#include "utils.hpp"


namespace fs = std::filesystem;


/**
 * Helper function to calculate velocity in layer.
 * @return Velocity (m/s) in layer at depth. No bounds checking is done.
 */
double layer_velocity(Layer layer, double depth) {
    return layer.intercept + layer.gradient * depth;
}


VelocityModel::VelocityModel(const std::vector<Layer>& layers, double x1, double y1, double x0,
                             double y0) :
        layers(layers),
        x0_(x0),
        x1_(x1),
        y0_(y0),
        y1_(y1) {
    if (layers.empty()) {
        throw std::invalid_argument("Velocity model has to contain layers.");
    }
    if (x_width() != y_width()) {
        throw std::invalid_argument(
            impl::Formatter()
            << "Different model widths along x and y axis not allowed. Widths are: x " << x_width()
            << " y " << y_width());
    }
    // First get interface depths
    auto top_layer = layers.front();
    m_interface_depths.push_back(top_layer.top_depth);
    std::for_each(layers.begin(), layers.end(),
                  [&](const auto& layer) { m_interface_depths.push_back(layer.bot_depth); });
    // Second get interface velocities with two special cases for first and last layer
    double velocity = 0.;
    _interface_velocities.push_back(velocity);
    for (auto& l : layers) {
        velocity = layer_velocity(l, l.top_depth);
        _interface_velocities.push_back(velocity);
        velocity = layer_velocity(l, l.bot_depth);
        _interface_velocities.push_back(velocity);
    }
    _interface_velocities.emplace_back(0);
}

VelocityModel read_velocity_file(const fs::path& filepath) {
    if (filepath.empty()) {
        throw std::invalid_argument("Can't read velocity model from empty path " +
                                    filepath.string());
    }
    if (not fs::exists(filepath)) {
        throw std::invalid_argument("Can't find file " + fs::absolute(filepath).string());
    }
    std::ifstream file(filepath);
    std::string line;
    std::istringstream iss;
    std::vector<double> widths;
    double depth_top, depth_bot, velocity_top, velocity_bottom, intercept, gradient;
    // extract widths
    int number_of_widths = 4;
    while (number_of_widths > 0) {
        std::getline(file, line);
        if (line.empty() or line.find('#') != std::string::npos) {
            continue;
        }
        iss = std::istringstream(line);
        iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        double x;
        iss >> x;
        widths.push_back(x);
        --number_of_widths;
    }
    // extract layers
    std::vector<Layer> ls;
    while (std::getline(file, line)) {
        // skip empty line and comments
        if (line.empty() or line.find('#') != std::string::npos) {
            continue;
        }
        iss = std::istringstream(line);
        iss >> depth_top;
        iss.ignore(1, ',');
        iss >> depth_bot;
        iss.ignore(1, ',');
        iss >> velocity_top;
        iss.ignore(1, ',');
        iss >> velocity_bottom;
        gradient = (velocity_bottom - velocity_top) / (depth_bot - depth_top);
        intercept = velocity_top - gradient * depth_top;
        ls.push_back({depth_top, depth_bot, intercept, gradient});
    }
    // widths from file are given as x0 x1 y0 y1, but constructor takes them as x1 y1 x0 y0
    return VelocityModel{ls, widths[1], widths[3], widths[0], widths[2]};
}

bool operator==(const Layer& l1, const Layer& l2) {
    return l1.gradient == l2.gradient and l1.top_depth == l2.top_depth and
           l1.bot_depth == l2.bot_depth and l1.intercept == l2.intercept;
}


const Layer& VelocityModel::operator[](size_t index) const {
    return layers[index];
}

std::optional<size_t> VelocityModel::layer_index(double x, double y, double z) const {
    if (not in_model(x, y, z)) {
        return std::nullopt;
    }
    auto greater = std::upper_bound(m_interface_depths.begin(), m_interface_depths.end(), z);
    return std::min(std::distance(m_interface_depths.begin(), greater) - 1,
                    static_cast<long int>(layers.size()) - 1);
}


std::optional<std::ptrdiff_t> VelocityModel::layer_index(double z) const {
    return layer_index(x0_, y0_, z);
}


std::optional<double> VelocityModel::eval_at(double x, double y, double z) const {
    auto index = layer_index(x, y, z);
    if (not index) {
        return {};
    }
    return layer_velocity(layers[index.value()], z);
}

std::pair<double, double> VelocityModel::interface_velocities(double z) const {
    auto index = layer_index(z).value();
    auto half_depth =
        layers[index].top_depth + 0.5 * (layers[index].bot_depth - layers[index].top_depth);
    // index is layer index, which goes from 0 to n-1. Indices for interface velocities go from
    // 0 to 2*(n+1)+1 since the first layer has 2 interface velocity pairs and every layer after
    // that adds one interface velocity pair.
    if (z < half_depth) {
        return {_interface_velocities[2 * index], _interface_velocities[2 * index + 1]};
    }
    return {_interface_velocities[2 * (index + 1)], _interface_velocities[2 * (index + 1) + 1]};
}

std::pair<double, double> VelocityModel::get_top_bottom() const {
    return {m_interface_depths.front(), m_interface_depths.back()};
}

bool VelocityModel::operator==(const VelocityModel& other) const {
    return layers == other.layers and x0_ == other.x0_ and x1_ == other.x1_ and y0_ == other.y0_ and
           y1_ == other.y1_;
}

bool VelocityModel::in_model(double x, double y, double z) const {
    bool in_depth = z >= m_interface_depths.front() and z <= m_interface_depths.back();
    bool in_widths = x0_ <= x and x <= x1_ and y0_ <= y and y <= y1_;
    return in_depth and in_widths;
}

std::pair<VelocityModel::iterator, VelocityModel::iterator>
VelocityModel::interface_velocities(double z1, double z2) const {
    auto [z_low, z_high] = std::minmax(z1, z2);
    auto index_low = layer_index(z_low).value();
    auto index_high = layer_index(z_high).value();
    return {_interface_velocities.begin() + 2 * (index_low + 1),
            _interface_velocities.begin() + 2 * (index_high + 1)};
}

std::vector<Layer>::const_iterator VelocityModel::begin() const {
    return layers.begin();
}

std::vector<Layer>::const_iterator VelocityModel::end() const {
    return layers.end();
}

std::ptrdiff_t VelocityModel::size() const {
    return layers.size();
}

const std::vector<double>& VelocityModel::interface_depths() const {
    return m_interface_depths;
}

Layer VelocityModel::get_layer(double x, double y, double z) const {
    auto index = layer_index(x, y, z);
    return layers[index.value()];
}

std::pair<double, double> VelocityModel::get_x_extent() const {
    return {x0_, x1_};
}

std::pair<double, double> VelocityModel::get_y_extent() const {
    return {y0_, y1_};
}

double VelocityModel::x0() const {
    return x0_;
}

double VelocityModel::x1() const {
    return x1_;
}

double VelocityModel::y0() const {
    return y0_;
}

double VelocityModel::y1() const {
    return y1_;
}

double VelocityModel::x_width() const {
    return x1_ - x0_;
}

double VelocityModel::y_width() const {
    return y1_ - y0_;
}

std::ostream& operator<<(std::ostream& os, const VelocityModel& model) {
    os << "VelocityModel(x0 = " << model.x0_ << ", x1 = " << model.x1_ << ", y0 = " << model.y0_
       << ", y1 = " << model.y1_ << ", z0 = " << model.interface_depths().front()
       << ", z1 = " << model.interface_depths().back() << ")";
    return os;
}

size_t VelocityModel::number_of_interfaces_between(double z1, double z2) const {
    // sort by depth since index is unsigned so we can only subtract the smaller one from the larger
    // one otherwise
    auto [min, max] = std::minmax(z1, z2);
    return layer_index(max).value() - layer_index(min).value();
}

bool VelocityModel::in_horizontal_extent(double x, double y) const {
    return x0_ <= x and x <= x1_ and y0_ <= y and y <= y1_;
}

double highest_velocity_between(double source_depth, double receiver_depth,
                                const VelocityModel& model) {
    auto receiver_index = model.layer_index(receiver_depth).value();
    auto source_index = model.layer_index(source_depth).value();
    auto receiver_layer = model[receiver_index];
    auto source_layer = model[source_index];
    // if points in same layer, maximum has to be at either point for linear velocity gradient
    if (source_index == receiver_index) {
        return std::max(layer_velocity(receiver_layer, receiver_depth),
                        layer_velocity(source_layer, source_depth));
    }
    auto [iterator_begin, iterator_end] = model.interface_velocities(receiver_depth, source_depth);
    auto max_vel = std::max_element(iterator_begin, iterator_end);
    return std::max({*max_vel, layer_velocity(source_layer, source_depth),
                     layer_velocity(receiver_layer, receiver_depth)});
}
