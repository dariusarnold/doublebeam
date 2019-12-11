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

Seismogram cut(const Seismogram& seismogram, double t0, double t1) {
    auto timestep = seismogram.timesteps[1] - seismogram.timesteps[0];
    size_t start_index = std::ceil(t0 / timestep);
    // shift by half a timestep to make t1 inclusive, meaning if t1 is equal to a time value of the
    // seismogram, this value will be included in the range.
    size_t end_index = std::ceil((t1 + timestep * 0.5) / timestep);
    end_index = std::min(seismogram.timesteps.size(), end_index);
    if (start_index > seismogram.timesteps.size()) {
        return {};
    }
    return {seismogram.data.begin() + start_index, seismogram.data.begin() + end_index};
}

Seismogram& SeismoData::get_seismogram(const Source& s, const Receiver& r) {
    // subtract 1 since files use 1 based indexing while vector uses zero based indexing
    auto seismogram_index = (s.index - 1) * seismograms.receivers.size() + (r.index - 1);
    return Seismogram(seismograms.seismograms.begin() + );
}

const Seismogram& SeismoData::get_seismogram(const Source& s, const Receiver& r) const {
    return seismograms.seismograms[(s.index - 1) * num_receivers() + (r.index - 1)];
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
        receivers(read_receiverfile(project_folder / receiver_file_name)) {
    read_all_seismograms(project_folder);
}


/**
 * Read timesteps from the first normal file with ending .txt in source_path.
 * @param source_path
 * @return
 */
std::vector<double> read_timesteps_from_some_seismogram(std::filesystem::path& source_path) {
    for (const auto& p : source_path) {
        if (std::filesystem::is_regular_file(p) and p.extension() == ".txt") {
            return read_timesteps(p);
        }
    }
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
    seismograms.reserve(sources.size() * receivers.size() * timesteps.size());
    // iterate over source directories and read seismograms
    auto i = 0;
    for (const auto& sourcepath : source_paths) {
        // read binary data if it exists, else fall back to text data
        if (auto binary_file = sourcepath / "data.bin"; fs::exists(binary_file)) {
            load_binary_seismograms2(
                binary_file, receivers.size(),
                gsl::span<double>(seismograms.data() + i * timesteps.size(),
                                  seismograms.data() + (i + 1) * timesteps.size()));
        } else {
            std::vector<fs::path> seismo_files;
            std::copy_if(fs::directory_iterator(sourcepath), fs::directory_iterator(),
                         std::back_inserter(seismo_files),
                         [](const auto& dir_entry) { return fs::is_regular_file(dir_entry); });
            std::sort(seismo_files.begin(), seismo_files.end());
            for (const auto& seismo_file : seismo_files) {
                auto ampl = read_amplitude(seismo_file);
                std::copy(ampl.begin(), ampl.end(), seismograms.data() + i * timesteps.size());
            }
        }
        ++i;
    }
}

size_t SeismoData::num_receivers() const {
    return seismograms.receivers.size();
}

size_t SeismoData::num_sources() const {
    return seismograms.sources.size();
}

AngularFrequency SeismoData::sampling_frequency() const {
    auto sample_timestep =
        seismograms.seismograms.front().timesteps[1] - seismograms.seismograms.front().timesteps[0];
    return AngularFrequency(2 * M_PI / sample_timestep);
}

size_t SeismogramPart::size() const {
    return std::distance(begin, end);
}
