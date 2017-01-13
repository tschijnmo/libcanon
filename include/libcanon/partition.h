/** Partitions of points.
 *
 * Partitions are special types of cosets of setwise stabilizer subgroups.  It
 * is mostly useful for the canonicalization of graphs, but the code here is
 * fully general and can be used on different problems.
 */

#ifndef LIBCANON_PARTITION_H
#define LIBCANON_PARTITION_H

#include <algorithm>
#include <cassert>
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
    /** Initializes a trivial partition.
     */

    explicit Partition(size_t size)
        : perm_(size)
        , begins_(size)
        , ends_(size)
    {
        std::iota(perm_.begin(), perm_.end(), 0);
        std::fill(begins_.begin(), begins_.end(), 0);
        std::fill(ends_.begin(), ends_.end(), size);
    }

    /** Initializes an empty partition.
     *
     * This is helpful for holding the memory temporarily.
     */

    Partition()
        : perm_()
        , begins_()
        , ends_()
    {
    }

    /** Copy constructs a partition from a given one.
     */

    Partition(const Partition& base) = default;

    /** Move constructs a partition.
     */

    Partition(Partition&& base) = default;

    /** Copy assigns a partition from another one.
     */

    Partition& operator=(const Partition& op) = default;

    /** Move assigns a partition.
     */

    Partition& operator=(Partition&& op) = default;

    /** Constructs a new partition by individualizing a given point.
     *
     * The given point should be in a non-singleton cell.  And it will be
     * placed before other points in the same cell.
     *
     * Note that a new partition is created, without the base partition mutated
     * in any way.
     */

    Partition(const Partition& base, Point point)
        : Partition(base)
    {
        auto begin_idx = begins_[point];
        auto end_idx = ends_[point];
        if (end_idx - begin_idx < 2)
            return; // Just to copy for singleton cells.

        for (size_t i = begin_idx; i < end_idx; ++i) {
            auto& curr_point = perm_[i];
            if (curr_point == point) {
                std::swap(curr_point, perm_[begin_idx]);
                ends_[point] = begin_idx + 1;
            } else {
                begins_[curr_point] = begin_idx + 1;
            }
        }
    }

    /** Splits a given cell by the given keys.
     *
     * The keys should be given by a functor taking a point to return an
     * totally ordered value or reference.
     *
     * Whether splitting actually occurred will be returned as a boolean.
     *
     * The cell can be given by any point in the cell.
     */

    template <typename T> bool split_by_key(Point point, T get_key)
    {
        auto cell_size = get_cell_size(point);
        if (cell_size == 1)
            return false;

        // Sorted the point within the cell according to the given key
        // function.
        std::sort(cell_begin(point), cell_end(point),
            [&](auto x, auto y) { return get_key(x) < get_key(y); });

        size_t base_idx = begins_[point];
        size_t curr_begin = base_idx; // Begin index of the current new cell.
        size_t n_groups = 0; // Number of groups found.

        for (size_t i = 0; i < cell_size; ++i) {

            auto idx = base_idx + i;
            auto curr_point = perm_[idx];

            if (i == 0 || get_key(curr_point) != get_key(perm_[idx - 1])) {
                // For a new group.

                ++n_groups;
                for (size_t j = curr_begin; j < idx; ++j) {
                    ends_[perm_[j]] = idx; // Set ends_
                }
                curr_begin = dest_idx;
            }

            begins_[curr_point] = curr_begin; // Set begins_
        }

        return n_groups > 1;
    }

    /** Tests if a partition is a discrete partition with only singletons.
     */

    bool is_discrete() const
    {
        return std::all_of(begin(), end(),
            [&](auto point) { return this->get_cell_size(point) == 1; });
    }

    /** Gets the size of the cell of a point. */

    size_t get_cell_size(Point point) const
    {
        return ends_[point] - begins_[point];
    }

    /** Gets the colour of a point.
     *
     * The assigned colour will just be an integer that is equal for all points
     * in the same cell and is ordered the same.
     */

    size_t get_colour(Point point) const { return ends_[point]; }

    /** Gets an iterator for the points of a cell.
     */

    Point_vec::const_iterator cell_begin(Point point) const
    {
        auto it = perm_.cbegin();
        std::advance(it, begins_[point]);
        return it;
    }

    /** Gets the sentinel for the iterator over points in a cell.
     */

    Point_vec::const_iterator cell_end(Point point) const
    {
        auto it = perm_.cbegin();
        std::advance(it, ends_[point]);
        return it;
    }

    /** Gets a point in the first cell of the partition.
     *
     * Note that this method cannot be called on an empty partition.
     */

    Point get_first() const
    {
        assert(size() > 0);
        return perm_[0];
    }

    /** Gets a point in the next cell of the given cell.
     *
     * The size of the partition will be returned if the given point is in the
     * last cell.
     */

    Point next_cell(Point point) const
    {
        auto end_idx = ends_[point];
        if (end_idx < size()) {
            return perm_[end_idx];
        } else {
            return size();
        }
    }

    /** Gets the size of the entire permutation domain.
     */

    size_t size() const { return perm_.size(); }

    /** Gets the pre-image array of a permutation for the partition.
     */

    const Point_vec& get_pre_imgs() const { return perm_; }

    /** Gets the pre-image of a given point.
     */

    Point operator>>(Point point) const { return perm_[point]; }

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
            : partition_(&partition)
            , curr_(curr)
        {
        }

        /** Increments a cell iterator.
         */

        Cell_it& operator++()
        {
            curr_ = partition_->next_cell(curr_);

            return *this;
        }

        /** Dereferences a cell iterator.
         */

        Point operator*() const { return curr_; }

        /** Compares two iterators for equality.
         */

        bool operator==(const Cell_it& other)
        {
            assert(this.partition_ == other.partition_);
            return this->curr_ == other.curr_;
        }

        /** Compares two iterators for inequality.
         */

        bool operator!=(const Cell_it& other) { return !(*this == other); }

    private:
        /** The current point.
         */

        Point curr_;

        /** Pointer to the partition to be looped over.
         */

        const Partition* partition_;
    };

    Cell_it begin() const { return { *this, get_first() }; }

    Cell_it end() const { return { *this, size() }; }

private:
    /** Permutation of the given points
     */

    Point_vec perm_;

    /** The begin index of the cell of all points.
     */

    std::vector<size_t> begins_;

    /** The end index of the cell of all points.
     */

    std::vector<size_t> ends_;
};

} // End namespace libcanon

#endif
