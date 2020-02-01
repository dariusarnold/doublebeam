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
        return "dataz.bin";
    }

    std::string_view get_binary_seismogram_extension() {
        auto index_of_dot = get_binary_seismogram_filename().rfind('.');
        return get_binary_seismogram_filename().substr(index_of_dot);
    }

    std::regex get_seismogram_file_regex() {
        return std::regex{get_seismogram_file_regex_str()};
    }

    std::string get_seismogram_file_regex_str() {
        return R"(seismo.z.*.sdu)";
    }

    UnitVector2 get_phi_hat() {
        return UnitVector2(1, 0);
    }
} // namespace config