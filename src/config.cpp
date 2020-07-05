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
#include "config.hpp"

namespace config {

    std::string_view get_source_filename() {
        return "sources.txt";
    }

    std::string_view get_receiver_filename() {
        return "receivers.txt";
    }

    std::string_view shotdata_foldername() {
        return "shotdata";
    }

    std::string_view get_binary_seismogram_filename() {
        return "data.bin";
    }

    std::string_view get_binary_seismogram_extension() {
        auto index_of_dot = get_binary_seismogram_filename().rfind('.');
        return get_binary_seismogram_filename().substr(index_of_dot);
    }

    std::regex get_seismogram_file_regex() {
        return std::regex{get_seismogram_file_regex_str()};
    }

    std::string get_seismogram_file_regex_str() {
        return R"(receiver_[0-9]*.txt)";
    }

    UnitVector2 get_phi_hat() {
        return UnitVector2(1, 0);
    }
} // namespace config