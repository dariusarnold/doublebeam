#include <algorithm>
#include <iterator>

#include "io.hpp"
#include "seismodata.hpp"
#include "utils.hpp"

SeismoData::SeismoData(const std::filesystem::path& project_folder,
                       const std::string& source_file_name, const std::string& receiver_file_name) :
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

Seismogram<double> SeismoData::get_seismogram(const Source& s, const Receiver& r) {
    // subtract 1 since files use 1 based indexing while vector uses zero based indexing
    auto seismogram_index = (s.index - 1) * num_receivers() + (r.index - 1);
    return Seismogram(seismograms.data.data() + (seismogram_index * num_samples()), num_samples(),
                      seismograms.timesteps.data(), num_samples());
}

Seismogram<const double> SeismoData::get_seismogram(const Source& s, const Receiver& r) const {
    // subtract 1 since files use 1 based indexing while vector uses zero based indexing
    auto seismogram_index = (s.index - 1) * num_receivers() + (r.index - 1);
    return Seismogram(seismograms.data.data() + seismogram_index * num_samples(), num_samples(),
                      seismograms.timesteps.data(), num_samples());
}

Seismogram<double> SeismoData::get_seismogram(const Source& s, const Receiver& r, double t0,
                                              double t1) {
    auto seismo = get_seismogram(s, r);
    ptrdiff_t begin_offset = std::ceil(t0 / timestep());
    // +1 because end should point to one past the end.
    ptrdiff_t end_offset = std::floor(t1 / timestep()) + 1;
    begin_offset = std::clamp(begin_offset, 0L, static_cast<ptrdiff_t>(seismo.size()) - 1);
    end_offset = std::clamp(end_offset, 0L, static_cast<ptrdiff_t>(seismo.size()));
    return Seismogram(seismo.data.data() + begin_offset, seismo.data.data() + end_offset,
                      seismo.timesteps.data() + begin_offset, seismo.timesteps.data() + end_offset);
}

Seismogram<const double> SeismoData::get_seismogram(const Source& s, const Receiver& r, double t0,
                                                    double t1) const {
    auto seismo = get_seismogram(s, r);
    ptrdiff_t begin_offset = std::ceil(t0 / timestep());
    // +1 because end should point to one past the end.
    ptrdiff_t end_offset = std::floor(t1 / timestep()) + 1;
    begin_offset = std::clamp(begin_offset, 0L, static_cast<ptrdiff_t>(seismo.size()) - 1);
    end_offset = std::clamp(end_offset, 0L, static_cast<ptrdiff_t>(seismo.size()));
    return Seismogram(seismo.data.data() + begin_offset, seismo.data.data() + end_offset,
                      seismo.timesteps.data() + begin_offset, seismo.timesteps.data() + end_offset);
}


const std::vector<Source>& SeismoData::sources() const {
    return seismograms.sources;
}

const std::vector<Receiver>& SeismoData::receivers() const {
    return seismograms.receivers;
}

Seismograms::Seismograms(const std::filesystem::path& project_folder,
                         const std::string& source_file_name,
                         const std::string& receiver_file_name) :
        sources(read_sourcefile(project_folder / source_file_name)),
        receivers(read_receiverfile(project_folder / receiver_file_name)),
        source_kd_tree(sources),
        receiver_kd_tree(receivers),
        data(),
        timesteps(),
        common_timestep(-1) {
    read_all_seismograms(project_folder);
    common_timestep = timesteps[1] - timesteps[0];
}


/**
 * Read timesteps from the first normal file with ending .txt in source_path.
 * @param source_path
 * @return
 */
std::vector<double> read_timesteps_from_some_seismogram(std::filesystem::path& source_path) {
    for (const auto& p : std::filesystem::directory_iterator(source_path)) {
        if (std::filesystem::is_regular_file(p) and p.path().extension() != ".bin") {
            return read_timesteps(p);
        }
    }
    throw std::runtime_error("Could not find seismogram in folder " + source_path.string());
}


void Seismograms::read_all_seismograms(const std::filesystem::path& project_folder) {
    namespace fs = std::filesystem;
    auto p = project_folder / "shotdata";
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
    timesteps = read_timesteps_from_some_seismogram(source_paths[0]);
    data.resize(sources.size() * receivers.size() * timesteps.size());
    // iterate over source directories and read seismograms
    auto source_index = 0;
    for (const auto& sourcepath : source_paths) {
        // read binary data if it exists, else fall back to text data
        if (auto binary_file = sourcepath / "data.bin"; fs::exists(binary_file)) {
            load_binary_seismograms2(
                binary_file, receivers.size(),
                gsl::span<double>(data.data() + source_index * timesteps.size() * receivers.size(),
                                  data.data() +
                                      (source_index + 1) * timesteps.size() * receivers.size()));
        } else {
            std::vector<fs::path> seismo_files;
            std::copy_if(fs::directory_iterator(sourcepath), fs::directory_iterator(),
                         std::back_inserter(seismo_files),
                         [](const auto& dir_entry) { return fs::is_regular_file(dir_entry); });
            std::sort(seismo_files.begin(), seismo_files.end());
            auto receiver_index = 0;
            for (const auto& seismo_file : seismo_files) {
                auto ampl = read_amplitude(seismo_file);
                std::copy(ampl.begin(), ampl.end(),
                          data.data() + receiver_index * timesteps.size() +
                              source_index * receivers.size() * timesteps.size());
                ++receiver_index;
            }
        }
        ++source_index;
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
    return AngularFrequency(2 * M_PI / timestep());
}

size_t SeismoData::num_samples() const {
    return seismograms.timesteps.size();
}

double SeismoData::timestep() const {
    return seismograms.common_timestep;
}