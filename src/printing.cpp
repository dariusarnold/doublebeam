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
#include "printing.hpp"


std::ostream& operator<<(std::ostream& os, const WaveType& wave_type) {
    return os << to_char(wave_type);
}

std::ostream& operator<<(std::ostream& os, std::vector<WaveType> v) {
    for (auto wave_type : v) {
        os << wave_type;
    }
    return os;
}
