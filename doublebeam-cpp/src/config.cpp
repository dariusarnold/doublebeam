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

    std::string_view get_binary_seismogram_filename_x() {
        return "datax.bin";
    }

    std::string_view get_binary_seismogram_filename_y() {
        return "datay.bin";
    }

    std::string_view get_binary_seismogram_filename_z() {
        return "dataz.bin";
    }

    std::string_view get_binary_seismogram_extension() {
        auto index_of_dot = get_binary_seismogram_filename_x().rfind('.');
        return get_binary_seismogram_filename_x().substr(index_of_dot);
    }

    std::regex get_seismogram_file_regex_x() {
        return std::regex{get_seismogram_file_regex_x()};
    }

     std::regex get_seismogram_file_regex_y() {
        return std::regex{get_seismogram_file_regex_y()};
    }

     std::regex get_seismogram_file_regex_z() {
        return std::regex{get_seismogram_file_regex_z()};
    }

    std::string get_seismogram_file_regex_str_x() {
        return R"(seismo.x.*.sdu)";
    }

    std::string get_seismogram_file_regex_str_y() {
        return R"(seismo.y.*.sdu)";
    }

    std::string get_seismogram_file_regex_str_z() {
        return R"(seismo.z.*.sdu)";
    }

    UnitVector2 get_phi_hat() {
        return UnitVector2(1, 0);
    }
} // namespace config