#ifndef DOUBLEBEAM_CPP_CONFIG_HPP
#define DOUBLEBEAM_CPP_CONFIG_HPP

#include <cmath>
#include <regex>
#include <string_view>

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
     * Filename for binary seismograms WITH file ending, eg. data.bin.
     * Binary seismograms are loaded before text files, since they are a lot
     * faster. This means if the source folder contains a file with this name,
     * seismograms will be loaded from it instead of the text files.
     * These three versions are for the three components (x, y, z).
     */
    std::string_view get_binary_seismogram_filename_x();
    std::string_view get_binary_seismogram_filename_y();
    std::string_view get_binary_seismogram_filename_z();

    /**
     * Regular expression used to find/select seismogram files.
     * All files in a source folder matching this expression will be loaded. This expression
     * therefore should only match files of one component (x, y, z), since only one component
     * can be used by the doublebeam algorithm.
     * @return
     */
    std::regex get_seismogram_file_regex_x();
    std::regex get_seismogram_file_regex_y();
    std::regex get_seismogram_file_regex_z();

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
    std::regex get_seismogram_file_regex_x();
    std::regex get_seismogram_file_regex_y();
    std::regex get_seismogram_file_regex_z();
    std::string get_seismogram_file_regex_str_x();
    std::string get_seismogram_file_regex_str_y();
    std::string get_seismogram_file_regex_str_z();

} // namespace config

#endif // DOUBLEBEAM_CPP_CONFIG_HPP
