#include "config.hpp"

namespace config {

    std::string_view get_source_filename() {
        return "sources.txt";
    }

    std::string_view get_receiver_filename() {
        return "receivers.txt";
    }

    std::string_view get_binary_filename() {
        return "data.bin";
    }

    std::string_view shotdata_foldername() {
        return "shotdata";
    }

    std::string_view get_binary_seismogram_filename() {
        return "data.bin";
    }
} // namespace config