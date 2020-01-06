#include <chrono>

#include "doublebeam.hpp"
#include "io.hpp"


int main(int argc, char* argv[]) {
    if (argc != 2) {
        throw std::runtime_error("Specify config file path as positional argument");
    }
    auto options = read_config_file(std::filesystem::path(argv[1]));
    std::ios_base::sync_with_stdio(false);
    auto vm = read_velocity_file(options.seismo_data_params.velocity_model_path);
    auto db = DoubleBeam(vm);
    auto source_beam_centres = seismo::grid_coordinates(
        options.sbc_params.x0, options.sbc_params.x1, options.sbc_params.y0, options.sbc_params.y1,
        options.sbc_params.z, options.sbc_params.num_x, options.sbc_params.num_y);
    FractureParameters fractures(
        options.fracture_params.phi_hat, options.fracture_params.num_orientations,
        options.fracture_params.spacings_min, options.fracture_params.spacings_max,
        options.fracture_params.num_spacings);
    auto data = SeismoData(options.seismo_data_params.data_path);
    auto a = std::chrono::high_resolution_clock::now();
    auto result =
        db.algorithm(source_beam_centres, options.target.position, data, fractures,
                     options.beam_params.width, options.beam_params.frequency,
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
        file << options;
        file << "\n[result]\n";
        file << result.data;
    }
    return 0;
}