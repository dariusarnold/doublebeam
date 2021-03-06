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
#ifndef DOUBLEBEAM_CPP_CONFIG_HPP
#define DOUBLEBEAM_CPP_CONFIG_HPP

#include <cmath>
#include <string_view>
#include <regex>

/**
 * This file contains constants, filenames, paths which are used by the main
 * program. Changing the values, recompiling the library and relinking with
 * the main program to change these settings.
 */


namespace config {

    /**
     * Name of source file, example: sources.txt
     */
    std::string_view get_source_filename();

    /**
     * Name of receiver file, example: receivers.txt
     */
    std::string_view get_receiver_filename();

    /**
     * Filename for binary seismograms.
     * WITH file ending, eg. data.bin
     */
    std::string_view get_binary_seismogram_filename();

    /**
     * Extension for binary seismograms.
     * This will be automatically determined by get_binary_seismogram_filename.
     */
    std::string_view get_binary_seismogram_extension();

    /**
     * Name of shotdata folder, which contains source subfolders.
     */
    std::string_view shotdata_foldername();


    struct UnitVector2 {
        UnitVector2(double xx, double yy) :
                x(xx / std::sqrt(xx * xx + yy * yy)), y(yy / std::sqrt(xx * xx + yy * yy)) {}

        double x, y;
    };

    /**
     * Set central direction of fracture normal (phi_hat) in scattering equation.
     * The algorithm will scan
     * @return
     */
    UnitVector2 get_phi_hat();

    /**
     * Regular expression used to find/select seismogram files.
     * All files in a source folder matching this expression will be loaded. This expression
     * therefore should only match files of one component (x, y, z), since only one component
     * can be used by the doublebeam algorithm.
     * The function returning std::regex is autocreated by the function returning std::string,
     * so only the string function has to be changed.
     */
    std::regex get_seismogram_file_regex();
    std::string get_seismogram_file_regex_str();

} // namespace config

#endif // DOUBLEBEAM_CPP_CONFIG_HPP
