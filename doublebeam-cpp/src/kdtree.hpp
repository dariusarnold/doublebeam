#ifndef DOUBLEBEAM_CPP_KDTREE_HPP
#define DOUBLEBEAM_CPP_KDTREE_HPP

#include <array>
#include <vector>

#include <fmt/format.h>
#include <nanoflann.hpp>

#include "units.hpp"
#include "raytracing_types.hpp"


/**
 * Search result class that improves access to the results of a radius search in the KDTree.
 * Instead of just the indices
 * Allows range based for iteration over the results of a radiusSearch.
 * Allows [] access to results.
 * @tparam PosType
 */
template <typename PosType>
class KDTreeSearchResults {
    using search_result = std::pair<size_t, double>;

    /**
     * Iterator for access to positions from the indices returned from a search.
     */
    class SearchResultsIterator {
    public:
        SearchResultsIterator(std::vector<search_result>::const_iterator it,
                              const std::vector<PosType>& pos) :
                matches_iterator(it), positions(pos) {}

        SearchResultsIterator& operator++() {
            ++matches_iterator;
            return *this;
        }

        bool operator!=(const SearchResultsIterator& other) {
            return matches_iterator != other.matches_iterator;
        }

        /**
         * Dereference iterator to position.
         * @return Current position in matches.
         */
        const PosType& operator*() {
            return positions[matches_iterator->first];
        }

    private:
        std::vector<search_result>::const_iterator matches_iterator;
        const std::vector<PosType>& positions;
    };

public:
    explicit KDTreeSearchResults(const std::vector<PosType>& pos) : positions(pos) {}

    // implicitly convert to reference so this class can be passed to radiusSearch instead of
    // the member.
    operator std::vector<search_result>&() {
        return matches;
    }

    const PosType& operator[](size_t index) {
        return positions[matches[index].first];
    }

    SearchResultsIterator begin() {
        return SearchResultsIterator(matches.begin(), positions);
    }

    SearchResultsIterator end() {
        return SearchResultsIterator(matches.end(), positions);
    }

    [[nodiscard]] size_t size() const {
        return matches.size();
    }

private:
    std::vector<search_result> matches{};
    const std::vector<PosType>& positions;
};

/**
 * KDTree class that wraps the adaptor and provides a search method to iterate over results
 * instead of returning indices.
 * @tparam PosType Type of Position object, can be a derived class.
 */
template <typename PosType>
class KDTree {
    /**
     * Adaptor to treat a std::vector<Position> (or a derived class, then given as template
     * parameter) as data for a kd tree. Functions prefixed with kdtree_ are there to implement the
     * DatasetAdaptor interface.
     * @tparam PositionType Position struct with x, y, z members.
     */
    class KDTreeAdaptor {
    public:
        explicit KDTreeAdaptor(const std::vector<PosType>& positions) : positions_m(positions) {}

        [[nodiscard]] size_t kdtree_get_point_count() const {
            return positions_m.size();
        }

        [[nodiscard]] double kdtree_get_pt(const size_t index, int dimension) const {
            switch (dimension) {
            case 0:
                return positions_m[index].x.get();
            case 1:
                return positions_m[index].y.get();
            case 2:
                return positions_m[index].z.get();
            default:
                throw std::invalid_argument(
                    fmt::format("Wrong dimension {} passed, maximum 2.", dimension));
            }
        }

        template <typename BBox>
        bool kdtree_get_bbox(BBox&) const {
            // return false to default to a standard bbox computation loop
            return false;
        }

        /**
         * Return reference to underlying position data set.
         */
        [[nodiscard]] const std::vector<PosType>& positions() const {
            return positions_m;
        }

    private:
        const std::vector<PosType>& positions_m;
    };


public:
    explicit KDTree(const std::vector<PosType>& positions_) :
            positions(KDTreeAdaptor(positions_)), tree(3, KDTreeAdaptor(positions_)) {
        tree.buildIndex();
    }

    /**
     * Find positions around center within radius.
     * @param center Point to find nearest positions around.
     * @param radius Positions with radius less than this are returned.
     */
    KDTreeSearchResults<PosType> get_positions(const Position& center, Meter radius) const {
        // Transfer to array since nanoflann expects [] indexable type
        std::array<double, 3> c{center.x.get(), center.y.get(), center.z.get()};
        KDTreeSearchResults results(tree.dataset.positions());
        tree.radiusSearch(c.data(), radius.get(), results, search_unsorted);
        return results;
    }

private:
    // set sorted to false since I dont care about sorting results by distance.
    nanoflann::SearchParams search_unsorted{32, 0, false};
    using kd_tree =
        nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDTreeAdaptor>,
                                            KDTreeAdaptor, 3>;
    KDTreeAdaptor positions;
    kd_tree tree;
};


#endif // DOUBLEBEAM_CPP_KDTREE_HPP
