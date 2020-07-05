/*
 * Copyright (C) 2019-2020  Darius Arnold
 *
 * This file is part of doublebeam.
 *
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include <algorithm>
#include <fstream>

#include "model.hpp"
#include "utils.hpp"


namespace fs = std::filesystem;


std::vector<Meter> get_interface_depths(const std::vector<Layer>& layers) {
    if (layers.empty()) {
        throw std::invalid_argument("Velocity model has to contain layers.");
    }
    std::vector<Meter> interface_depths;
    interface_depths.reserve(layers.size() + 1);
    auto top_layer = layers.front();
    interface_depths.push_back(top_layer.top_depth);
    std::for_each(layers.begin(), layers.end(),
                  [&](const auto& layer) { interface_depths.push_back(layer.bot_depth); });
    return interface_depths;
}


VelocityModel::VelocityModel(const std::vector<Layer>& layers, Meter x1, Meter y1, Meter x0,
                             Meter y0) :
        m_interface_depths(get_interface_depths(layers)),
        layers_m(layers),
        x0_(x0),
        x1_(x1),
        y0_(y0),
        y1_(y1) {}

VelocityModel read_velocity_file(const fs::path& filepath) {
    if (filepath.empty()) {
        throw std::invalid_argument("Can't read velocity model from empty path " +
                                    filepath.string());
    }
    if (not fs::exists(filepath)) {
        throw std::invalid_argument("Can't find file " + filepath.string());
    }
    std::ifstream file(filepath);
    std::string line;
    std::istringstream iss;
    std::vector<Meter> widths;
    Meter depth_top, depth_bot;
    Velocity velocity;
    // extract widths
    int number_of_widths = 4;
    while (number_of_widths > 0) {
        std::getline(file, line);
        if (line.empty() or line.find('#') != std::string::npos) {
            continue;
        }
        iss = std::istringstream(line);
        iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        Meter x;
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
        iss >> velocity;
        ls.push_back({depth_top, depth_bot, velocity});
    }
    // widths from file are given as x0 x1 y0 y1, but constructor takes them as x1 y1 x0 y0
    return VelocityModel{ls, widths[1], widths[3], widths[0], widths[2]};
}

bool operator==(const Layer& l1, const Layer& l2) {
    return l1.top_depth == l2.top_depth and l1.bot_depth == l2.bot_depth and
           l1.velocity == l2.velocity;
}


const Layer& VelocityModel::operator[](size_t index) const {
    return layers_m[index];
}

std::optional<size_t> VelocityModel::layer_index(Meter x, Meter y, Meter z) const {
    if (not in_model(x, y, z)) {
        return std::nullopt;
    }
    auto greater = std::upper_bound(m_interface_depths.begin(), m_interface_depths.end(), z);
    return std::min(std::distance(m_interface_depths.begin(), greater) - 1,
                    static_cast<long int>(layers_m.size()) - 1);
}


std::optional<std::size_t> VelocityModel::layer_index(Meter z) const {
    return layer_index(x0_, y0_, z);
}


std::optional<Velocity> VelocityModel::eval_at(Meter x, Meter y, Meter z) const {
    auto index = layer_index(x, y, z);
    if (not index) {
        return {};
    }
    return layers_m[index.value()].velocity;
}

std::optional<Velocity> VelocityModel::eval_at(const Position& position) const {
    return eval_at(position.x, position.y, position.z);
}

Meter layer_height(const Layer& layer) {
    return layer.bot_depth - layer.top_depth;
}

VelocityModel::InterfaceVelocities VelocityModel::interface_velocities(Meter z) const {
    size_t index = layer_index(z).value();
    const Velocity outside_velocity(0);
    auto half_depth = layers_m[index].top_depth + 0.5 * layer_height(layers_m[index]);
    if (index == 0 and z < half_depth) {
        return {outside_velocity, layers_m[index].velocity};
    }
    if (index == num_layers() - 1 and z > half_depth) {
        return {layers_m[index].velocity, outside_velocity};
    }
    if (z < half_depth) {
        // interface above current layer
        return {layers_m[index - 1].velocity, layers_m[index].velocity};
    }
    // interface below current layer
    return {layers_m[index].velocity, layers_m[index + 1].velocity};
}

std::pair<Meter, Meter> VelocityModel::get_top_bottom() const {
    return {m_interface_depths.front(), m_interface_depths.back()};
}

bool VelocityModel::operator==(const VelocityModel& other) const {
    return layers_m == other.layers_m and x0_ == other.x0_ and x1_ == other.x1_ and y0_ == other.y0_ and
           y1_ == other.y1_;
}

bool VelocityModel::in_model(Meter x, Meter y, Meter z) const {
    bool in_depth = z >= m_interface_depths.front() and z <= m_interface_depths.back();
    bool in_widths = x0_ <= x and x <= x1_ and y0_ <= y and y <= y1_;
    return in_depth and in_widths;
}

std::vector<Layer>::const_iterator VelocityModel::begin() const {
    return layers_m.begin();
}

std::vector<Layer>::const_iterator VelocityModel::end() const {
    return layers_m.end();
}

std::size_t VelocityModel::num_layers() const {
    return layers_m.size();
}

const std::vector<Meter>& VelocityModel::interface_depths() const {
    return m_interface_depths;
}

Layer VelocityModel::get_layer(Meter x, Meter y, Meter z) const {
    auto index = layer_index(x, y, z);
    return layers_m[index.value()];
}

std::pair<Meter, Meter> VelocityModel::get_x_extent() const {
    return {x0_, x1_};
}

std::pair<Meter, Meter> VelocityModel::get_y_extent() const {
    return {y0_, y1_};
}

Meter VelocityModel::get_x0() const {
    return x0_;
}

Meter VelocityModel::get_x1() const {
    return x1_;
}

Meter VelocityModel::get_y0() const {
    return y0_;
}

Meter VelocityModel::get_y1() const {
    return y1_;
}

Meter VelocityModel::x_width() const {
    return x1_ - x0_;
}

Meter VelocityModel::y_width() const {
    return y1_ - y0_;
}

std::ostream& operator<<(std::ostream& os, const VelocityModel& model) {
    os << "VelocityModel(x0 = " << model.x0_ << ", x1 = " << model.x1_ << ", y0 = " << model.y0_
       << ", y1 = " << model.y1_ << ", z0 = " << model.interface_depths().front()
       << ", z1 = " << model.interface_depths().back() << ")";
    return os;
}

size_t VelocityModel::number_of_interfaces_between(Meter z1, Meter z2) const {
    // sort by depth since index is unsigned so we can only subtract the smaller one from the larger
    // one otherwise
    auto [min, max] = std::minmax(z1, z2);
    return layer_index(max).value() - layer_index(min).value();
}

bool VelocityModel::in_horizontal_extent(Meter x, Meter y) const {
    return x0_ <= x and x <= x1_ and y0_ <= y and y <= y1_;
}

Velocity highest_velocity_between(Meter depth1, Meter depth2, const VelocityModel& model) {
    auto receiver_index = model.layer_index(depth2).value();
    auto source_index = model.layer_index(depth1).value();
    // if points in same layer, maximum has to be at either point for linear velocity gradient
    Velocity max_velocity(0);
    for (auto i = std::min(receiver_index, source_index);
         i <= std::max(receiver_index, source_index); ++i) {
        if (model[i].velocity > max_velocity) {
            max_velocity = model[i].velocity;
        }
    }
    return max_velocity;
}
