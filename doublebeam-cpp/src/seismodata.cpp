#include <algorithm>
#include <iterator>
#include <regex>

#include <boost/range/adaptor/indexed.hpp>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include "config.hpp"
#include "io.hpp"
#include "seismodata.hpp"
#include "utils.hpp"


SeismoData::SeismoData(const std::filesystem::path& project_folder,
                       std::string_view source_file_name, std::string_view receiver_file_name) :
        seismograms(project_folder, source_file_name, receiver_file_name) {}

bool PositionWithIndex::operator==(const PositionWithIndex& other) const {
    return x == other.x and y == other.y and z == other.z and index == other.index;
}

std::ostream& operator<<(std::ostream& os, const Source& s) {
    os << "Source_" << s.index << "(" << s.x << ", " << s.y << ", " << s.z << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Receiver& r) {
    os << "Receiver_" << r.index << "(" << r.x << ", " << r.y << ", " << r.z << ")";
    return os;
}

const std::vector<Source>& SeismoData::sources() const {
    return seismograms.sources;
}

const std::vector<Receiver>& SeismoData::receivers() const {
    return seismograms.receivers;
}

Seismograms::Seismograms(const std::filesystem::path& project_folder,
                         std::string_view source_file_name, std::string_view receiver_file_name) :
        sources(read_sourcefile(project_folder / source_file_name)),
        receivers(read_receiverfile(project_folder / receiver_file_name)),
        source_kd_tree(sources),
        receiver_kd_tree(receivers),
        datax(),
        datay(),
        dataz(),
        timesteps(),
        common_timestep(-1) {
    read_all_seismograms(project_folder);
    common_timestep = Second(timesteps[1] - timesteps[0]);
    std::for_each(
        sources.begin(), sources.end(),
        [start_index = sources[0].index - 1](Source& source) { source.index -= start_index; });
    std::for_each(receivers.begin(), receivers.end(),
                  [start_index = receivers[0].index - 1](Receiver& receiver) {
                      receiver.index -= start_index;
                  });
}


/**
 * Read timesteps from the first normal file with ending .txt in source_path.
 * @param source_path
 * @return
 */
std::vector<double> read_timesteps_from_some_seismogram(std::filesystem::path& source_path) {
    for (const auto& p : std::filesystem::directory_iterator(source_path)) {
        if (std::filesystem::is_regular_file(p) and
            p.path().extension() != config::get_binary_seismogram_extension()) {
            return read_timesteps(p);
        }
    }
    throw std::runtime_error("Could not find seismogram in folder " + source_path.string());
}


/**
 * Get sorted vector of seismograms in a folder.
 */
std::vector<std::filesystem::path>
get_sorted_seismogram_files(const std::filesystem::path& sourcepath,
                            const std::regex& seismo_file_regex) {
    namespace fs = std::filesystem;
    std::vector<fs::path> seismo_files;
    std::copy_if(fs::directory_iterator(sourcepath), fs::directory_iterator(),
                 std::back_inserter(seismo_files), [&](const auto& dir_entry) {
                     return fs::is_regular_file(dir_entry) and
                            std::regex_search(dir_entry.path().filename().string(),
                                              seismo_file_regex);
                 });
    std::sort(seismo_files.begin(), seismo_files.end());
    return seismo_files;
}


/**
 * Check if binary data is available for source folder and load that first,
 * else load data from text seismograms.
 * @param sourcepath Path to source folder.
 * @param num_receivers Number of receivers per source.
 * @param subarray Array large enough to hold num_receivers * timesteps values.
 * Should be part of the larger seismogram array, where data for all sources is stored.
 */
void load_seismogram_data(const std::filesystem::path& sourcepath, size_t num_receivers,
                          const gsl::span<double> subarray, const std::regex& seismo_file_regex,
                          std::string_view binary_seismo_filename) {
    namespace fs = std::filesystem;
    if (auto binary_file = sourcepath / binary_seismo_filename;
        // read binary data if it exists,
        fs::exists(binary_file)) {
        fmt::print("Loading {} (binary)\n", sourcepath.string());
        load_binary_seismograms(binary_file, num_receivers, subarray);
    } else {
        // fall back to text data
        auto seismo_files = get_sorted_seismogram_files(sourcepath, seismo_file_regex);
        if (seismo_files.empty()) {
            throw std::invalid_argument(fmt::format("Found no seismograms in {}\n"
//                                                    "Was searching for {}.\n",
//                                                    sourcepath,
//                                                    config::get_seismogram_file_regex_str()
                                                    ));
        }
        fmt::print("Loading {} (text)\n", sourcepath.string());
        for (const auto& seismo_file : seismo_files | boost::adaptors::indexed()) {
            auto ampl = read_amplitude(seismo_file.value());
            std::copy(ampl.begin(), ampl.end(),
                      subarray.data() + seismo_file.index() * ampl.size());
        }
    }
}


/**
 * Finds shotdata folder and builds a sorted vector of all source folders.
 */
std::vector<std::filesystem::path>
get_sorted_source_folders(const std::filesystem::path& project_folder) {
    namespace fs = std::filesystem;
    auto p = project_folder / config::shotdata_foldername();
    // error if shotdata folder missing
    if (not fs::is_directory(p) or not fs::exists(p)) {
        throw std::runtime_error(impl::Formatter()
                                 << "No directory shotdata in " << project_folder);
    }
    // build sorted list of source directories
    // sort because the folder names are assumed to be related to the source index, eg source #1
    // will have folder name "source_1" or similar.
    std::vector<fs::path> source_paths;
    std::copy_if(fs::directory_iterator(p), fs::directory_iterator(),
                 std::back_inserter(source_paths),
                 [](const auto& dir_entry) { return dir_entry.is_directory(); });
    std::sort(source_paths.begin(), source_paths.end());
    return source_paths;
}


void Seismograms::read_all_seismograms(const std::filesystem::path& project_folder) {
    auto source_paths = get_sorted_source_folders(project_folder);
    timesteps = read_timesteps_from_some_seismogram(source_paths[0]);
    datax.resize(sources.size() * receivers.size() * timesteps.size());
    datay.resize(sources.size() * receivers.size() * timesteps.size());
    dataz.resize(sources.size() * receivers.size() * timesteps.size());
    const size_t number_of_datapoints_per_source = receivers.size() * timesteps.size();
    // iterate over source directories and read seismograms
    for (const auto& sourcepath : source_paths | boost::adaptors::indexed()) {
        load_seismogram_data(
            sourcepath.value(), receivers.size(),
            gsl::span<double>(datax.data() + sourcepath.index() * number_of_datapoints_per_source,
                              number_of_datapoints_per_source),
            config::get_seismogram_file_regex_x(), config::get_binary_seismogram_filename_x());
        load_seismogram_data(
            sourcepath.value(), receivers.size(),
            gsl::span<double>(datay.data() + sourcepath.index() * number_of_datapoints_per_source,
                              number_of_datapoints_per_source),
            config::get_seismogram_file_regex_y(), config::get_binary_seismogram_filename_y());
        load_seismogram_data(
            sourcepath.value(), receivers.size(),
            gsl::span<double>(dataz.data() + sourcepath.index() * number_of_datapoints_per_source,
                              number_of_datapoints_per_source),
            config::get_seismogram_file_regex_y(), config::get_binary_seismogram_filename_z());
    }
}

KDTreeSearchResults<Source> SeismoData::get_sources(const Position& position, Meter radius) const {
    return seismograms.source_kd_tree.get_positions(position, radius);
}

KDTreeSearchResults<Receiver> SeismoData::get_receivers(const Position& position,
                                                        Meter radius) const {
    return seismograms.receiver_kd_tree.get_positions(position, radius);
}

size_t SeismoData::num_receivers() const {
    return seismograms.receivers.size();
}

size_t SeismoData::num_sources() const {
    return seismograms.sources.size();
}

AngularFrequency SeismoData::sampling_frequency() const {
    return AngularFrequency(2 * M_PI / timestep().get());
}

size_t SeismoData::num_samples() const {
    return seismograms.timesteps.size();
}

Second SeismoData::timestep() const {
    return seismograms.common_timestep;
}

Second SeismoData::time_length() const {
    return Second{seismograms.timesteps.back() - seismograms.timesteps.front()};
}
