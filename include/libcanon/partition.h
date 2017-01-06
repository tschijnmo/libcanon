/** Partitions of point.
 *
 * Partitions are special types of cosets of setwise stabilizer subgroups.  It
 * is mostly useful for the canonicalization of graphs, but the code here is
 * fully general and can be used on different problems.
 *
 */

#ifndef LIBCANON_PARTITION_H
#define LIBCANON_PARTITION_H

#include <algorithm>
#include <vector>

#include <libcanon/perm.h>

namespace libcanon {

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

    /** Initializes an empty partition.
     */
    Partition()
        : perm()
        , begins()
        , ends()
    {
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
     * The keys should be given by a functor taking a point to return an
     * totally ordered value or reference.
     *
     * Whether splitting actually occurred will be returned as a boolean.
     */

    template <typename T> bool split_by_key(Point point, T get_key)
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
            [&](auto x, auto y) { return get_key(x) < get_key(y); });

        size_t base_idx = begins[point];
        size_t curr_begin = base_idx; // Begin index of the current new cell.
        size_t n_groups = 0; // Number of groups found.

        for (size_t i = 0; i < cell_size; ++i) {

            auto dest_idx = base_idx + i;
            perm[dest_idx] = orig_cell[sorted_idxes[i]];

            if (i == 0
                || get_key(sorted_idxes[i]) != get_key(sorted_idxes[i - 1])) {
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

    /** Tests if a partition is a discrete partition with only singletons. */

    bool is_discrete() const
    {
        return std::all_of(begin(), end(),
            [&](auto point) { return get_cell_size(point) == 1; });
    }

    /** Gets the size of the cell of a point. */

    size_t get_cell_size(Point point) const
    {
        return ends[point] - begins[point];
    }

    /** Gets the colour of a point.
     *
     * The assigned colour will just be an integer that is equal for all points
     * in the same cell and is ordered the same.
     */

    size_t get_colour(Point point) const { return ends[point]; }

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

    /** Gets the pre-image array of a permutation for the partition. */

    const Array& get_pre_imgs() const { return perm; }

    /** Gets the pre-image of a given point.
     */

    Point operator>>(Point point) const { return perm[point]; }

    /** Makes a simple permutation object for the partition.
     */

    Simple_perm make_perm() const { return { perm.begin(), perm.end() }; }

    /** Iterator type for cells in the partition.
     *
     * Dereferencing the iterator will give a point in the cell.
     */

    class Cell_it {
    public:
        /** Constructs a cell iterator.
         */

        Cell_it(const Partition& partition, Point curr)
            : partition{ partition }
            , curr{ curr }
        {
        }

        /** Increments a cell iterator.
         */

        Cell_it& operator++()
        {
            curr = partition.next_cell(curr);

            return *this;
        }

        /** Dereferences a cell iterator.
         */

        Point operator*() const { return curr; }

        /** Compares two iterators for equality.
         */

        bool operator==(const Cell_it& other)
        {
            return this->curr == other.curr;
        }

        /** Compares two iterators for inequality.
         */

        bool operator!=(const Cell_it& other) { return !(*this == other); }

    private:
        Point curr;
        const Partition& partition;
    };

    Cell_it begin() const { return { *this, get_first() }; }

    Cell_it end() const { return { *this, size() }; }

private:
    /** Permutation of the given points */
    Array perm;

    /** The begin and end index of the cell of all points. */
    Array begins;
    Array ends;
};

} // End namespace libcanon

#endif
