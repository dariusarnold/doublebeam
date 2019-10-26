#ifndef DOUBLEBEAM_CPP_SEISMODATA_HPP
#define DOUBLEBEAM_CPP_SEISMODATA_HPP

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>

struct Source {
    // coordinates of source
    double x, y, z;
    // index of source in source file and project directory.
    size_t index;
    bool operator==(const Source& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Source& s);
};

struct Receiver {
    // coordinates of receiver
    double x, y, z;
    // index of receiver in directory structure
    size_t index;
    // index of receiver in receiver file and project directory
    bool operator==(const Receiver& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Receiver& r);
};

struct Seismogram {
    Seismogram() = default;
    Seismogram(std::vector<double>&& d) : data(d) {}

    // amplitude data
    std::vector<double> data;
};


/**
 * Cut seismogram to time samples between time t0 and time t1.
 * Both t0 and t1 are inclusive.
 * Seismogram and t have to have the same size.
 * @param seismogram
 * @param t
 * @param t0
 * @param t1
 * @return New seismogram containing only the amplitude samples between t0 and t1.
 */
Seismogram cut(const Seismogram& seismogram, const std::vector<double>& t, double t0, double t1);


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
    std::vector<Seismogram> data;

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
    Seismogram& operator()(const Source& s, const Receiver& r);
    /**
     * Get access to list of sources.
     */
    const std::vector<Source>& sources() const;
    /**
     * Get access to list of receivers.
     * @return
     */
    const std::vector<Receiver>& receivers() const;

    /**
     * Get number of receivers.
     * @return
     */
    size_t num_receivers() const;

    /**
     * Get number of sources.
     * @return
     */
    size_t num_sources() const;

    /**
     * Get access to common timesteps of all seismograms
     * @return
     */
    const std::vector<double>& timesteps() const;

private:
    Seismograms seismograms;
};


#endif // DOUBLEBEAM_CPP_SEISMODATA_HPP
