#include <chrono>

#include <boost/program_options.hpp>
#include <fmt/format.h>

#include "doublebeam.hpp"
#include "io.hpp"


/**
 * Find first unused filename with the scheme "resultX.txt" where X is an increasing integer
 * starting at 0.
 */
std::filesystem::path find_unused_result_filename() {
    std::string result_filename = "result";
    std::filesystem::path result_path(result_filename + "0.txt");
    int index = 0;
    while (std::filesystem::exists(result_path)) {
        ++index;
        result_path = std::filesystem::path(result_filename + (std::to_string(index)) + ".txt");
    }
    return result_path;
}


int main(int argc, char* argv[]) {
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    try {
        DoubleBeamOptions opts;
        // positional options for config file
        po::positional_options_description positional;
        positional.add("config_file", 1);

        // store config file on own description so the keyword argument is not shown in usage
        po::options_description hidden("Hidden options");
        // clang-format off
        hidden.add_options()
            ("config_file,c", po::value<fs::path>()->value_name("PATH"), "Path to config file.");
        // clang-format on

        // options only allowed on command line
        po::options_description cmd("Commandline options");
        // clang-format off
        cmd.add_options()
            ("help,h", "Produce help message.");
        // clang-format on

        // add general options, can be specified on command line and/or in config file
        po::options_description general("Config/keyword arguments");
        // clang-format off
        general.add_options()
            ("data.model,m", po::value(&opts.seismo_data_params.velocity_model_path)->value_name("PATH")->required(),
                "Path to velocity model file")
            ("data.path,d", po::value(&opts.seismo_data_params.data_path)->value_name("PATH")->required(),
                "Path to folder containing seismic data.")
            ("target.x", po::value(&opts.target.position.x)->required(),
                "x (meter) coordinate of target position.")
            ("target.y", po::value(&opts.target.position.y)->required(),
                "y (meter) coordinate of target position.")
            ("target.z", po::value(&opts.target.position.z)->required(),
                "z (meter) coordinate of target position.")
            ("source_beam_centers.x0", po::value(&opts.sbc_params.x0)->required(),
                "Min x coordinate of source beam center rectangle.")
            ("source_beam_centers.x1", po::value(&opts.sbc_params.x1)->required(),
                "Max x coordinate of source beam center rectangle.")
            ("source_beam_centers.y0", po::value(&opts.sbc_params.y0)->required(),
                "Min y coordinate of source beam center rectangle.")
            ("source_beam_centers.y1", po::value(&opts.sbc_params.y1)->required(),
                "Max y coordinate of source beam center rectangle.")
            ("source_beam_centers.z", po::value(&opts.sbc_params.z)->default_value(0_meter),
                "Depth in meter of source beam centers. Usually at surface.")
            ("source_beam_centers.num_x", po::value(&opts.sbc_params.num_x)->value_name("NX")->required(),
                "Number of source beam centers along x axis.")
            ("source_beam_centers.num_y", po::value(&opts.sbc_params.num_y)->value_name("NY")->required(),
                "Number of source beam centers along y axis.")
            ("fractures.spacing_min", po::value(&opts.fracture_params.spacings_min)->required(),
                "Min fracture spacing in meter.")
            ("fractures.spacing_max", po::value(&opts.fracture_params.spacings_max)->required(),
                "Max fracture spacing in meter.")
            ("fractures.num_spacings", po::value(&opts.fracture_params.num_spacings)->required(),
                "Number of fracture spacings to scan between spacing_min and spacing_max..")
            ("fractures.num_orientations", po::value(&opts.fracture_params.num_orientations)->required(),
                "Number of fracture orientations to scan.")
            ("beam.width,w", po::value(&opts.beam_params.width)->required(),
                "Initial width of Gauss beam in meter.")
            ("beam.frequency,f", po::value(&opts.beam_params.frequency)->required(),
                "Initial frequency of Gauss beam in hertz.")
            ("beam.window_length,l", po::value(&opts.beam_params.window_length)->required(),
                "Window length (seconds) of cut window around traveltime of beam.")
            ("beam.max_stacking_distance", po::value(&opts.beam_params.max_stacking_distance)->required(),
                "Maximum stacking distance at beam end point at the surface. Only sources/receivers"
                " within this distance will be used for stacking.")
            ;
        // clang-format on

        po::options_description all_options, options_shown_in_usage;
        all_options.add(cmd).add(general).add(hidden);
        options_shown_in_usage.add(cmd).add(general);
        po::variables_map options;
        po::store(
            po::command_line_parser(argc, argv).options(all_options).positional(positional).run(),
            options);
        // if help is given, ignore other options, print help and exit
        if (options.count("help")) {
            fmt::print("Gaussian double beam fracture parameter determination.\n\n");
            fmt::print("Usage:\n{} [config_file_path] [keyword_arguments]\n\n",
                       fs::path(argv[0]).filename().string());
            // TODO this should have the same width as the options description (80 columns)
            fmt::print("config_file_path can give a text file in the ini-like Boost.program_options"
                       " format. In this file options for the program can be specified. Command "
                       "line options override the values in the config file.\n\n");
            fmt::print("{}\n", options_shown_in_usage);
            exit(0);
        }
        // if path to config file is given, read options from it. Since options given first are
        // preferred, this means cmd options override the config file parameters.
        if (options.count("config_file")) {
            auto config_path = options["config_file"].as<fs::path>();
            std::ifstream config_file{config_path};
            if (config_file) {
                po::store(po::parse_config_file(config_file, all_options), options);
            } else {
                throw po::error(fmt::format("Failed to open config file {}", config_path));
            }
        }
        po::notify(options);

        auto vm = read_velocity_file(opts.seismo_data_params.velocity_model_path);
        auto db = DoubleBeam(vm);
        auto source_beam_centres = seismo::grid_coordinates(
            opts.sbc_params.x0, opts.sbc_params.x1, opts.sbc_params.y0, opts.sbc_params.y1,
            opts.sbc_params.z, opts.sbc_params.num_x, opts.sbc_params.num_y);
        FractureParameters fractures(math::Vector2{1, 0}, opts.fracture_params.num_orientations,
                                     opts.fracture_params.spacings_min,
                                     opts.fracture_params.spacings_max,
                                     opts.fracture_params.num_spacings);
        auto data = SeismoData(opts.seismo_data_params.data_path);
        auto a = std::chrono::high_resolution_clock::now();
        auto result =
            db.algorithm(source_beam_centres, opts.target.position, data, fractures,
                         opts.beam_params.width, hertz_to_angular(opts.beam_params.frequency),
                         opts.beam_params.window_length, opts.beam_params.max_stacking_distance);
        auto b = std::chrono::high_resolution_clock::now();
        std::cout << "Runtime db : "
                  << std::chrono::duration_cast<std::chrono::seconds>(b - a).count() << " s"
                  << std::endl;
        auto result_path = find_unused_result_filename();
        std::ofstream file{result_path};
        if (file.is_open()) {
            file << opts;
            file << "\n[result]\n";
            file << result.data;
        }
        fmt::print("Saved result in {}\n", result_path.string());
        return 0;
    } catch (const po::unknown_option& er) {
        fmt::print(std::cerr, "{}\nMaybe a typo?", er.what());
    } catch (const po::error& er) {
        fmt::print(std::cerr, "{}\n", er.what());
        exit(-1);
    }
}