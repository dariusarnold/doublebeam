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
#include <cmath>

#include "model.hpp"
#include "raytracing_helpers.hpp"
#include "utils.hpp"


Slowness snells_law(const Slowness& slowness, Velocity v_above, Velocity v_below,
                    WaveType wave_type) {
    if (wave_type == WaveType::Reflected) {
        return {slowness.px, slowness.py, -slowness.pz};
    }
    double minus_plus = -1.;
    if (not seismo::ray_direction_down(slowness)) {
        // for transmitted upgoing ray
        std::swap(v_above, v_below);
    }
    // handle only special case of horizontal interface where normal is vertical.
    // n should be oriented to the side the transmitted wave propagates for the
    // minus_plus relation to work,
    double nz = std::copysign(1., slowness.pz.get());
    // dot product can be simplified since bx and ny are 0, nz is -1
    double p_dot_n = nz * slowness.pz.get();
    double eps = std::copysign(1., p_dot_n);
    // since in the original formula everything subtracted from p is multiplied by n
    // only the pz component changes for horizontal interfaces.
    InverseVelocity pz(
        slowness.pz.get() -
        (p_dot_n + minus_plus * eps *
                       std::sqrt(1. / (v_below.get() * v_below.get()) -
                                 1. / (v_above.get() * v_above.get()) + p_dot_n * p_dot_n)) *
            nz);
    return {slowness.px, slowness.py, pz};
}

Slowness snells_law(const RayState& old_state, const VelocityModel& model, WaveType wave_type) {
    auto [position, slowness, travel_time, arclength] = old_state;
    auto [v_above, v_below] = model.interface_velocities(position.z);
    return snells_law(slowness, v_above, v_below, wave_type);
}