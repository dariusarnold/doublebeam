#ifndef DOUBLEBEAM_CPP_SEISMODATA_HPP
#define DOUBLEBEAM_CPP_SEISMODATA_HPP

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>

#include <gsl/span>

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

struct Source : public PositionWithIndex {
    friend std::ostream& operator<<(std::ostream& os, const Source& s);
};

struct Receiver : public PositionWithIndex {
    friend std::ostream& operator<<(std::ostream& os, const Receiver& r);
};

/**
 * Seismogram class can represent full seismogram or partial (cut out).
 */
struct Seismogram {
public:
    Seismogram(double* data_begin, size_t data_size, double* timesteps_begin,
               size_t timesteps_size) :
            data(data_begin, data_size), timesteps(timesteps_begin, timesteps_size) {};

    gsl::span<double> data;
    gsl::span<double> timesteps;
};


/**
 * Cut seismogram to time samples between time t0 and time t1.
 * Both t0 and t1 are inclusive.
 * @param seismogram
 * @param t0 start time
 * @param t1 end time
 * @return Part of seismogram containing only the amplitude samples between t0 and t1.
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
    // All amplitude data in one vector, ordered first by source, then by receiver.
    // For S sources and R receivers, will contain S*R seismograms.
    // First R seismograms will belong to source 1, the following R seismograms to source 2 and so
    // on.
    std::vector<double> seismograms;
    // Common timesteps of all seismograms.
    std::vector<double> timesteps;

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
    Seismogram& get_seismogram(const Source& s, const Receiver& r);

    const Seismogram& get_seismogram(const Source& s, const Receiver& r) const;

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
