#ifndef DOUBLEBEAM_CPP_SEISMODATA_HPP
#define DOUBLEBEAM_CPP_SEISMODATA_HPP

struct Source {
    double x, y, z;
    size_t index;
    bool operator==(const Source& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Source& s);
};

struct Receiver {
    double x, y, z;
    size_t index;
    bool operator==(const Receiver& other) const;
    friend std::ostream& operator<<(std::ostream& os, const Receiver& r);
};
#endif // DOUBLEBEAM_CPP_SEISMODATA_HPP
