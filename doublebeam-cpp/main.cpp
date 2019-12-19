#include <chrono>
#include <complex>

#include "doublebeam.hpp"
#include "eigen_helpers.hpp"
#include "io.hpp"
#include "printing.hpp"
#include "raytracing.hpp"
#include "twopoint.hpp"
#include "units.hpp"
#include "utils.hpp"


int main(int argc, char* argv[]) {
    if (argc != 2) {
        throw std::runtime_error("Specify config file path as positional argument");
    }
    auto options = read_config_file(std::filesystem::path(argv[1]));
    std::ios_base::sync_with_stdio(false);
    auto vm = read_velocity_file(options.model_params.path);
    auto db = DoubleBeam(vm);
    auto source_beam_centres = seismo::grid_coordinates(
        options.sbc_params.x0, options.sbc_params.x1, options.sbc_params.y0, options.sbc_params.y1,
        options.sbc_params.z, options.sbc_params.num_x, options.sbc_params.num_y);
    // TODO fracture depth is unused
    FractureParameters fractures(
        options.fracture_params.phi_hat, options.fracture_params.num_orientations,
        options.fracture_params.spacings_min, options.fracture_params.spacings_max,
        options.fracture_params.num_spacings);
    auto data = SeismoData(options.seismo_data_params.path);
    auto a = std::chrono::high_resolution_clock::now();
    auto result =
        db.algorithm(source_beam_centres, options.target.position, data, fractures,
                     options.seismo_data_params.source_frequency, Meter(options.beam_params.width),
                     hertz_to_angular(Frequency(options.beam_params.frequency)),
                     options.beam_params.window_length, options.beam_params.max_stacking_distance);
    auto b = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime db : " << std::chrono::duration_cast<std::chrono::seconds>(b - a).count()
              << " s" << std::endl;
    std::string result_filename = "result";
    std::filesystem::path result_path(result_filename + "0.txt");
    int index = 0;
    while (std::filesystem::exists(result_path)) {
        ++index;
        result_path = std::filesystem::path(result_filename + (std::to_string(index)) + ".txt");
    }
    std::ofstream file{result_path};
    if (file.is_open()) {
        file << result.data;
    }
    return 0;
}