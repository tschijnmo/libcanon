/** Canonicalization of edge-labelled directed acyclic graphs.
 *
 * Edge-labelled directly acyclic graphs (ELDAGs) are generalization of the
 * normal DAGs with a label for each edge.  And the labels for the edges can
 * possibly be permuted.  The algorithm here is based on the spirit of the
 * algorithm by B McKay for normal graphs.
 */

#ifndef LIBCANON_ELDAG_H
#define LIBCANON_ELDAG_H

#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>

namespace libcanon {

/** Data type for an Eldag.
 *
 * This is basically a CSR format to stored the children of each node.  This
 * data structure should have very high cache friendliness.
 *
 * Note that variant information like the colour and symmetry at each node is
 * not stored here.
 *
 * It is implemented here as a struct since the structure should be transparent
 * to developers.
 */

struct Eldag {

    /** The edges from each node to other nodes.
     *
     * Stored as list of indices for the index of the nodes that they are
     * connected to.
     */

    std::vector<size_t> edges;

    /** The start index for the connections of each node.
     *
     * This is basically the same as the IA array for CSR.
     */

    std::vector<size_t> ia;

    /** Gets the number of nodes in the Eldag. */

    size_t size() const { return ia.size() - 1; }
};

/** A partition of points.
 *
 * Here each point is identified by a natural number.  Here the partition is
 * stored is a slightly redundant way.  But it supports fast query of both a
 * permutation and the colour of any point.
 */

class Partition {
public:
    /** Array always means a vector of the primary integer type in this class.
     */
    using Array = std::vector<size_t>;

    /** Initializes a trivial partition. */

    Partition(size_t size)
        : perm(size)
        , begins(size)
        , ends(size)
    {
        std::iota(perm.begin(), perm.end(), 0);
        std::fill(begins.begin(), begins.end(), 0);
        std::fill(ends.begin(), ends.end(), size);
    }

    /** Copy constructs a partition from a given one. */

    Partition(const Partition& base) = default;

    /** Constructs a new partition by individualizing a given point.
     *
     * The given point should be in a non-singleton cell.  And it will be
     * placed before other points in the same cell.
     */

    Partition(const Partition& base, Point point)
        : Partition(base)
    {
        auto begin_idx = begins[point];
        auto end_idx = ends[point];
        if (end_idx - begin_idx < 2)
            return; // Just to copy for singleton cells.

        for (size_t i = begin_idx; i < end_idx; ++i) {
            auto& curr_point = perm[i];
            if (curr_point == point) {
                std::swap(curr_point, perm[begin_idx]);
                ends[point] = begin_idx + 1;
            } else {
                begins[curr_point] = begin_idx + 1;
            }
        }
    }

    /** Splits a given partition by the given keys.
     *
     * The keys should be given by an random access container holding the keys
     * in the same order of the points in the cell as given by iteration by
     * `begin_cell` and `end_cell`.
     *
     * Whether splitting actually occurred will be returned as a boolean.
     */

    template <typename T> bool split_by_key(Point point, const T& keys)
    {
        auto cell_size = get_cell_size(point);
        if (cell_size == 1)
            return false;

        // The original cell, it should be in correspondence with the keys.
        std::vector<Point> orig_cell(cell_begin(point), cell_end(point));

        // Sorted indices within the cell.
        std::vector<size_t> sorted_idxes(cell_size);
        std::iota(sorted_idxes.begin(), sorted_idxes.end(), 0);
        std::sort(sorted_idxes.begin(), sorted_idxes.end(),
            [&](auto x, auto y) { return keys[x] < keys[y]; });

        size_t base_idx = begins[point];
        size_t curr_begin = base_idx; // Begin index of the current new cell.
        size_t n_groups = 0; // Number of groups found.

        for (size_t i = 0; i < cell_size; ++i) {

            auto dest_idx = base_idx + i;
            perm[dest_idx] = orig_cell[sorted_idxes[i]];

            if (i == 0 || keys[sorted_idxes[i]] != keys[sorted_idxes[i - 1]]) {
                // For a new group.

                ++n_groups;
                for (size_t i = curr_begin; i < dest_idx; ++i) {
                    ends[perm[i]] = dest_idx;
                }
                curr_begin = dest_idx;
            }

            begins[perm[i]] = curr_begin;
        }

        return n_groups > 1;
    }

    /** Tests if a partition is a leaf partition with only singletons. */

    bool is_leaf() const
    {
        for (size_t i = 0; i < size(); ++i) {
            if (get_cell_size(i) != 1)
                return false;
        }

        return true;
    }

    /** Gets the size of the cell of a point. */

    size_t get_cell_size(Point point) const
    {
        return ends[point] - begins[point];
    }

    /** Gets an iterator for the points of a cell. */

    Array::const_iterator cell_begin(Point point) const
    {
        auto it = perm.cbegin();
        std::advance(it, begins[point]);
        return it;
    }

    /** Gets the sentinel for the iterator over points in a cell. */

    Array::const_iterator cell_end(Point point) const
    {
        auto it = perm.cbegin();
        std::advance(it, ends[point]);
        return it;
    }

    /** Gets a point in the first cell of the partition. */

    Point get_first() const { return perm[0]; }

    /** Gets a point in the next cell of the given cell.
     *
     * The size of the partition will be returned if the given point is in the
     * last cell.
     */

    Point next_cell(Point point) const
    {
        auto end_idx = ends[point];
        if (end_idx < size()) {
            return perm[end_idx];
        } else {
            return size();
        }
    }

    /** Gets the size of the entire permutation domain. */

    size_t size() const { return perm.size(); }

private:
    /** Permutation of the given points */
    Array perm;

    /** The begin and end index of the cell of all points. */
    Array begins;
    Array ends;
};

}

#endif // LIBCANON_ELDAG_H
