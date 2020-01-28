#ifndef DOUBLEBEAM_CPP_CONFIG_HPP
#define DOUBLEBEAM_CPP_CONFIG_HPP

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
     * Filename for binary seismograms.
     */
    std::string_view get_binary_seismogram_filename();

    /**
     * Name of shotdata folder, which contains source subfolders.
     */
    std::string_view shotdata_foldername();

} // namespace config

#endif // DOUBLEBEAM_CPP_CONFIG_HPP
