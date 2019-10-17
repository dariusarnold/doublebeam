#ifndef DOUBLEBEAM_CPP_SEISMODATA_HPP
#define DOUBLEBEAM_CPP_SEISMODATA_HPP

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>

struct Source {
    double x, y, z;
    // index of source in source file and project directory.
    size_t index;
    bool operator==(const Source& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Source& s);
};

struct Receiver {
    double x, y, z;
    size_t index;
    // index of receiver in receiver file and project directory
    bool operator==(const Receiver& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Receiver& r);
};

struct Seismograms {
    /**
     * Read all seismograms, sources and receivers from folder.
     * @param project_folder
     */
    Seismograms(const std::filesystem::path& project_folder,
                const std::string& source_file_name = "sources.txt",
                const std::string& receiver_file_name = "receivers.txt");

    std::vector<Source> sources;
    std::vector<Receiver> receivers;
    // common time steps of all seismograms
    std::vector<double> times;
    std::vector<std::vector<double>> data;

private:
    void read_all_seismograms(const std::filesystem::path& project_folder);
};

class SeismoData {
public:
    SeismoData(const std::filesystem::path& project_folder,
               const std::string& source_file_name = "sources.txt",
               const std::string& receiver_file_name = "receivers.txt");
    /**
     * Return specific seismogram given from shot at source recorded at receiver.
     * Won't check if source or receiver are at the correct positions. Use only sources and
     * receivers as read from source/receiver file that belong to the same project.
     * @param s
     * @param r
     * @return
     */
    std::vector<double>& operator()(const Source& s, const Receiver& r);
    /**
     * Get access to list of sources.
     */
    const std::vector<Source>& sources() const;
    /**
     * Get access to list of receivers.
     * @return
     */
    const std::vector<Receiver>& receivers() const;

private:
    Seismograms seismograms;
};


#endif // DOUBLEBEAM_CPP_SEISMODATA_HPP
