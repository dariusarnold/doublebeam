#include <algorithm>
#include <iterator>

#include "io.hpp"
#include "seismodata.hpp"
#include "utils.hpp"

SeismoData::SeismoData(const std::filesystem::path& project_folder,
                       const std::string& source_file_name, const std::string& receiver_file_name) :
        seismograms(project_folder, source_file_name, receiver_file_name) {}

bool Source::operator==(const Source& other) const {
    return x == other.x and y == other.y and z == other.z and index == other.index;
}

std::ostream& operator<<(std::ostream& os, const Source& s) {
    os << "Source_" << s.index << "(" << s.x << ", " << s.y << ", " << s.z << ")";
    return os;
}

bool Receiver::operator==(const Receiver& other) const {
    return x == other.x and y == other.y and z == other.z and index == other.index;
}

std::ostream& operator<<(std::ostream& os, const Receiver& r) {
    os << "Receiver_" << r.index << "(" << r.x << ", " << r.y << ", " << r.z << ")";
    return os;
}

Seismogram cut(const Seismogram& seismogram, double t0, double t1) {
    // TODO test/specify what happens for t0/t1 outside of timesteps range
    if (seismogram.data.size() != seismogram.timesteps.size()) {
        throw std::invalid_argument(
            impl::Formatter() << "Got seismogram of size " << seismogram.data.size() << " but only "
                              << seismogram.timesteps.size() << " time steps.");
    }
    // If I could assume that all time samples are evenly spaced, it would be easy to calculate the
    // index of the start/end cut. Since I can't make this assumption for all input, a binary search
    // is applied to find the indices of t0 and t1.
    auto i1 = std::lower_bound(seismogram.timesteps.begin(), seismogram.timesteps.end(), t0);
    // t1 can't be before t0, so decrease search range
    auto i2 = std::upper_bound(i1, seismogram.timesteps.end(), t1);
    auto start_index = std::distance(seismogram.timesteps.begin(), i1);
    auto end_index = std::distance(seismogram.timesteps.begin(), i2);
    Seismogram out{std::vector<double>{seismogram.timesteps.begin() + start_index,
                                       seismogram.timesteps.begin() + end_index},
                   std::vector<double>{seismogram.data.begin() + start_index,
                                       seismogram.data.begin() + end_index}};
    return out;
}

Seismogram& SeismoData::operator()(const Source& s, const Receiver& r) {
    // subtract 1 since files use 1 based indexing while vector uses zero based indexing
    return seismograms.seismograms[(s.index - 1) * seismograms.receivers.size() + (r.index - 1)];
}

const Seismogram& SeismoData::operator()(const Source& s, const Receiver& r) const {
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
    // iterate over source directories and read seismograms
    for (const auto& sourcepath : source_paths) {
        // read binary data if it exists, else fall back to text data
        if (auto binary_file = sourcepath / "data.bin"; fs::exists(binary_file)) {
            auto data = load_binary_seismograms(receivers.size(), binary_file);
            for (auto seismogram : data) {
                seismograms.emplace_back(std::move(seismogram.first), std::move(seismogram.second));
            }
        } else {
            std::vector<fs::path> seismo_files;
            std::copy_if(fs::directory_iterator(sourcepath), fs::directory_iterator(),
                         std::back_inserter(seismo_files),
                         [](const auto& dir_entry) { return fs::is_regular_file(dir_entry); });
            std::sort(seismo_files.begin(), seismo_files.end());
            for (const auto& seismo_file : seismo_files) {
                auto seismogram = read_seismogram(seismo_file);
                seismograms.emplace_back(std::move(seismogram.first), (seismogram.second));
            }
        }
    }
}

size_t SeismoData::num_receivers() const {
    return seismograms.receivers.size();
}

size_t SeismoData::num_sources() const {
    return seismograms.sources.size();
}