#ifndef DOUBLEBEAM_CPP_SEISMODATA_HPP
#define DOUBLEBEAM_CPP_SEISMODATA_HPP

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <vector>

#include <gsl/span>
#include <boost/container_hash/hash.hpp>

#include "kdtree.hpp"
#include "units.hpp"


struct PositionWithIndex : Position {
    PositionWithIndex(Meter xx, Meter yy, Meter zz, size_t ind) :
            Position(xx, yy, zz), index(ind) {}

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
    Source(Meter xx, Meter yy, Meter zz, size_t ind) : PositionWithIndex(xx, yy, zz, ind) {}

    friend std::ostream& operator<<(std::ostream& os, const Source& s);
};

struct Receiver : public PositionWithIndex {
    Receiver(Meter xx, Meter yy, Meter zz, size_t ind) : PositionWithIndex(xx, yy, zz, ind) {}

    friend std::ostream& operator<<(std::ostream& os, const Receiver& r);
};

/**
 * Seismogram class can represent full seismogram or partial (cut out).
 */
template <typename T>
struct Seismogram {
public:
    template <typename Container>
    Seismogram(Container amplitudes, Container times) :
            Seismogram(amplitudes.data(), amplitudes.size(), times.data(), times.size()) {}

    Seismogram(T* data_begin, size_t data_size, T* timesteps_begin, size_t timesteps_size) :
            data(data_begin, data_size), timesteps(timesteps_begin, timesteps_size) {}

    Seismogram(T* data_begin, T* data_end, T* timesteps_begin, T* timesteps_end) :
            data(data_begin, data_end), timesteps(timesteps_begin, timesteps_end) {}

    [[nodiscard]] size_t size() const {
        return data.size();
    }

    const gsl::span<T> data;
    const gsl::span<T> timesteps;

    bool operator==(const Seismogram& other) const {
        return data.data() == other.data.data() and timesteps.data() == other.timesteps.data() and
               size() == other.size();
    }

    friend std::size_t hash_value(const Seismogram& seismogram) {
        std::size_t seed = 0;
        boost::hash_combine(seed, seismogram.data.data());
        boost::hash_combine(seed, seismogram.timesteps.data());
        boost::hash_combine(seed, seismogram.size());
        return seed;
    }
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
    KDTree<Source> source_kd_tree;
    KDTree<Receiver> receiver_kd_tree;

    // All amplitude data in one vector, ordered first by source, then by receiver.
    // For S sources and R receivers, will contain S*R seismograms.
    // First R seismograms will belong to source 1, the following R seismograms to source 2 and so
    // on.
    std::vector<double> data;
    // Common timesteps of all seismograms.
    std::vector<double> timesteps;

    // Distance between two time samples.
    double common_timestep;

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
    Seismogram<double> get_seismogram(const Source& s, const Receiver& r);

    [[nodiscard]] Seismogram<const double> get_seismogram(const Source& s, const Receiver& r) const;

    /**
     * Cut seismogram to time samples between time t0 and time t1.
     * Both t0 and t1 are inclusive.
     * @param s Source position
     * @param r Receiver position
     * @param t0 start time
     * @param t1 end time
     * @return Part of seismogram containing only the amplitude samples between t0 and t1.
     */
    Seismogram<double> get_seismogram(const Source& s, const Receiver& r, double t0, double t1);

    [[nodiscard]] Seismogram<const double> get_seismogram(const Source& s, const Receiver& r,
                                                          double t0, double t1) const;

    /**
     * Get access to list of sources.
     */
    [[nodiscard]] const std::vector<Source>& sources() const;
    /**
     * Get access to list of receivers.
     * @return
     */
    [[nodiscard]] const std::vector<Receiver>& receivers() const;

    /**
     * Get number of receivers.
     * @return
     */
    [[nodiscard]] size_t num_receivers() const;

    /**
     * Get number of sources.
     * @return
     */
    [[nodiscard]] size_t num_sources() const;

    /**
     * Get number of samples in a seismogram.
     * @return
     */
    [[nodiscard]] size_t num_samples() const;

    /**
     * Get time difference between two samples, ie. the sampling period.
     * @return
     */
    [[nodiscard]] double timestep() const;

    /**
     * Return sampling frequency in rad/s.
     * This assumes all seismograms are sampled with the same timestep.
     */
    [[nodiscard]] AngularFrequency sampling_frequency() const;

    /**
     * Get Sources within radius around position.
     */
    KDTreeSearchResults<Source> get_sources(const Position& position, Meter radius) const;
    /**
     * Get Receivers within radius around position.
     */
    KDTreeSearchResults<Receiver> get_receivers(const Position& position, Meter radius) const;


private:
    Seismograms seismograms;
};


#endif // DOUBLEBEAM_CPP_SEISMODATA_HPP
