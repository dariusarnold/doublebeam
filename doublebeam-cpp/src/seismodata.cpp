#include "seismodata.hpp"
#include "io.hpp"
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

Seismogram& SeismoData::operator()(const Source& s, const Receiver& r) {
    // subtract 1 since files use 1 based indexing while vector uses zero based indexing
    return seismograms.data[(s.index - 1) * seismograms.receivers.size() + (r.index - 1)];
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
    std::vector<fs::path> source_paths;
    // build sorted list of source directories
    for (auto source_dir : fs::directory_iterator(p)) {
        if (not fs::is_directory(source_dir)) {
            // skip normal files without error
            continue;
        }
        source_paths.push_back(source_dir);
    }
    // sort because the folder names are assumed to be related to the source index, eg source #1
    // will have folder name "source_1" or similar.
    std::sort(source_paths.begin(), source_paths.end());
    // iterate over source directories and read seismograms
    for (const auto& sourcepath : source_paths) {
        std::vector<fs::path> seismo_files;
        for (auto seismo_file : fs::directory_iterator(sourcepath)) {
            if (fs::is_regular_file(seismo_file)) {
                seismo_files.push_back(seismo_file);
            }
        }
        std::sort(seismo_files.begin(), seismo_files.end());
        for (const auto& seismo_file : seismo_files) {
            data.emplace_back(read_amplitude(seismo_file));
        }
    }
}