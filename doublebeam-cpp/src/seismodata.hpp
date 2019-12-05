#ifndef DOUBLEBEAM_CPP_SEISMODATA_HPP
#define DOUBLEBEAM_CPP_SEISMODATA_HPP

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>

#include "units.hpp"


struct PositionWithIndex {
    double x, y, z;
    // index of source in source file or of receiver in receiver file.
    size_t index;
    /**
     * Compare position and index for equality.
     * @param other
     * @return
     */
    bool operator==(const PositionWithIndex& other) const;
};

struct Source : public PositionWithIndex{
    friend std::ostream& operator<<(std::ostream& os, const Source& s);
};

struct Receiver : public PositionWithIndex{
    friend std::ostream& operator<<(std::ostream& os, const Receiver& r);
};

struct Seismogram {
    Seismogram() {}

    Seismogram(std::vector<double> t, std::vector<double> d) :
            timesteps(std::move(t)), data(std::move(d)) {}

    // time data
    std::vector<double> timesteps{};
    // amplitude data
    std::vector<double> data{};
};


/**
 * Cut seismogram to time samples between time t0 and time t1.
 * Both t0 and t1 are inclusive.
 * Seismogram and t have to have the same size.
 * @param seismogram
 * @param t0
 * @param t1
 * @return New seismogram containing only the amplitude samples between t0 and t1.
 */
Seismogram cut(const Seismogram& seismogram, double t0, double t1);


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
    std::vector<Seismogram> seismograms{};

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

    const Seismogram& operator()(const Source& s, const Receiver& r) const;

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
     * Return sampling frequency in rad/s.
     * This assumes all seismograms are sampled with the same timestep.
     */
     AngularFrequency sampling_frequency() const;

private:
    Seismograms seismograms;
};


#endif // DOUBLEBEAM_CPP_SEISMODATA_HPP
