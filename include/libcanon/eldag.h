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
#include <memory>
#include <numeric>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <libcanon/perm.h>
#include <libcanon/string.h>

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

    /** Creates an Eldag of a given size.
     *
     * Note that here the vectors are only reserved space.  The actual sizes
     * are still empty.
     */

    Eldag(size_t n_nodes, size_t n_edges)
        : edges{}
        , ia{}
    {
        ia.reserve(n_nodes + 1);
        edges.reserve(n_edges);
    }
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

/** The data type to given symmetries for the nodes in an Eldag.
 *
 * Here the container is fixed to vector for the ease of coding.  Null values
 * are considered to be the absence of any symmetries.
 *
 * This type is going to be used for the public driver functions for Eldag
 * canonicalization.
 */

template <typename A> using Symms = std::vector<const Sims_transv<Perm<A>>*>;

/** The data type for the permutations to be applied to nodes in an Eldag.
 *
 * This is going to be in the return value of the public driver function for
 * Eldag canonicalization.  Note that here the value owns the reference to the
 * permutations.
 */

template <typename A> using Perms = std::vector<std::unique_ptr<Perm<A>>>;

/** Acts permutations on an Eldag.
 *
 * Here the global permutation of the nodes should be given in a
 * permutation-like object supporting operators `<<` and `>>`, the local
 * permutation acting on each node should be given as an indexable container
 * containing pointer-like objects pointing to the actual permutations on each
 * node.  Null values indicates the absence of any permutation.
 */

template <typename G, typename L>
Eldag act(const Eldag& eldag, const G& gl_perm, const L& perms)
{
    size_t size = eldag.size();
    Eldag res{ size, eldag.edges.size() };

    // The index that the next connection is going to be written to.
    size_t curr_idx = 0;

    // Constructs the nodes in the result one-by-one.
    for (size_t curr = 0; curr < size; ++curr) {

        Point pre_img = gl_perm >> curr;
        size_t prev_base = eldag.ia[pre_img];
        size_t valence = eldag.ia[pre_img + 1] - prev_base;

        for (size_t i = 0; i < valence; ++i) {
            size_t offset = i;

            if (perms[curr]) {
                offset = *perms[curr] >> offset;
            }

            Point conn = gl_perm << eldag.edges[prev_base + offset];
            res.edges[curr_idx] = conn;
            ++curr_idx;
        }

        res.ia.push_back(curr_idx);
    }

    res.ia.push_back(curr_idx);

    return res;
}

namespace internal {

    /** Data type for owned references to symmetries. */

    template <typename A>
    using Owned_symms = std::vector<std::unique_ptr<Sims_transv<Perm<A>>>>;

    /** Data type for borrowed references to permutations. */

    template <typename A> using Borrowed_perms = std::vector<const Perm<A>*>;

    /**
     * Data type for a coset in the canonicalization of an Eldag.
     */

    template <typename A> class Eldag_coset {
    public:
        /** Constructs a coset based on initial partition and symmetries.
         *
         * This is useful for the creation of root coset.
         */

        Eldag_coset(const Partition& init_part, Symms<A> init_symms)
            : partition{ init_part }
            , symms(init_symms)
            , perms(init_part.size(), nullptr)
            , refined_symms()
            , refined_perms()
        {
        }

        /** Constructs a coset by acting a permutation.
         *
         * This is useful for the pruning of the refinement tree.  Note that
         * the created coset is not actually fully valid one.  Just the
         * individualized point is correctly set for the sake of removing the
         * actual coset.
         */

        Eldag_coset(const Eldag_coset& base, const Simple_perm& perm)
            : individualized{ perm >> base.get_individualized() }
            , partition()
            , perms()
            , symms()
            , refined_perms()
            , refined_symms()
        {
        }

        /** Constructs a coset by individualizing a point.
         *
         * This is useful for the individualization step.  Note that this
         * construction is better used on cosets after refinement.
         */

        Eldag_coset(const Eldag_coset& base, Point point)
            : individualized(point)
            , partition(base.partition, point)
            , symms(base.symms)
            , perms(base.perms)
            , refined_perms()
            , refined_symms()
        {
        }

        /** Individualizes the refined partition.
         *
         * This method can be called after the testing of leaf, where the
         * partition is going to be refined.
         */

        std::vector<Eldag_coset> individualize_first_nonsingleton()
        {
            auto cell = std::find_if(partition.begin(), partition.end(),
                [&](auto cell) { return partition.get_cell_size(cell) > 1; });

            // When used correctly, we should find a non-singleton cell here.

            std::vector<Eldag_coset> res;
            std::transform(partition.cell_begin(*cell),
                partition.cell_end(*cell), std::back_inserter(res),
                [&](auto point) { return Eldag_coset(*this, point); });

            return res;
        }

        /** Tests if the result is a leaf.
         *
         * In addition to testing, most importantly, the Weisfeiler-Lehman
         * refinement will be performed here.
         */

        bool is_leaf()
        {
            refine();
            return partition.is_discrete();
        }

        /** Gets the point individualized during its construction.
         */

        Point get_individualized() const { return individualized; }

        /** Gets the size of the graph.
         */

        size_t size() const { return partition.size(); }

        /** Compares two cosets for equality.
         *
         * Here only the point that is individualized is compared.  So it is
         * actually only applicable to peers.
         */

        bool operator==(const Eldag_coset& other) const
        {
            return get_individualized() == other.get_individualized();
        }

        /** Evaluates the hash of a coset.
         *
         * Note that only the individualized point is used.  So it is only
         * useful for retrieval amongst it peers.
         */

        size_t hash() const
        {
            return get_individualized();
        }

    private:
        /** Refines the currently holding partition and symmetry.
         *
         * Most of the actual work of the canonicalization of Eldag actually
         * comes here.
         */

        void refine() 
        {
            // TODO: Add the actual refinement here.
        }

        //
        // Data fields
        //

        /** The partition of graph nodes.
         *
         * It can be coarse initially, and then refined.
         */

        Partition partition;

        /** The current permutations applied to the nodes.
         */

        Borrowed_perms<A> perms;

        /** The current symmetries for each of the nodes.
         */

        Symms<A> symms;

        /** The permutation on each node after refinement.
         */

        Perms<A> refined_perms;

        /** Refined symmetries for each node.
         */

        Owned_symms<A> refined_symms;

        /**
         * The point that is individualized in the construction of this coset.
         */

        Point individualized;
    };

    /** Data type for a permutation on an Eldag.
     *
     * Inside it, we have not only the global permutation of the nodes, but the
     * permutation on valences of each nodes as well.
     *
     * Note that the permutation stored in the graph automorphism are not of
     * this type, but rather simple permutations of the nodes.  This is
     * achieved by delegate to proxy types for the expression `g | ~h`.
     */

    template <typename A> struct Eldag_perm {

        /** The permutations on each of the nodes. */

        Perms<A> perms;

        /** The global partition of the nodes.
         *
         * It should be a singleton partition already.
         */

        Partition partition;

        /** The class for the inversion of an Eldag permutation. */

        struct Inv_eldag_perm {
            const Eldag_perm& operand;
        };

        /** Inverses a permutation.
         *
         * No actual inversion is done.  Just a proxy is returned.
         */

        Inv_eldag_perm operator~() const { return { *this }; }

        /** Forms an automorphism.
         *
         * This function will be called by the generic canonicalization
         * function by the evaluation of the expression `p | ~q`.
         */

        Simple_perm operator|(const Inv_eldag_perm& inv)
        {
            size_t size = partition.size();
            std::vector<Point> pre_imgs(size);

            const auto& p = partition.get_pre_imgs();
            const auto& q = inv.operand.get_pre_imgs();

            for (size_t i = 0; i < size; ++i) {
                // By ~q, q[i] <- i;
                pre_imgs[q[i]] = p[i];
            }

            return { pre_imgs.begin(), pre_imgs.end() };
        }

        /** Gets the pre-image of a given point.
         */

        Point operator>>(Point point) const { return partition >> point; }
    };

    /**
     * The actual refiner for Eldag canonicalization.
     */

    template <typename A> class Eldag_refiner {
    public:
        //
        // Types required by the refiner protocol.
        //

        using Coset = Eldag_coset<A>;
        using Perm = Eldag_perm<A>;
        using Transv = Sims_transv<Simple_perm>;
        using Structure = Eldag;
        using Container = std::unordered_map<Eldag, Perm>;

        /** Refines the given coset.
         */

        std::vector<Coset> refine(const Coset& coset)
        {
            return coset.individualize_first_nonsingleton();
        }

        /** Decides if a coset is a leaf.
         */

        bool is_leaf(Coset& coset) { return coset.is_leaf(); }

        /** Gets a permutation from a coset.
         */

        Perm get_a_perm(const Coset& coset) { return coset.get_a_perm(); }

        /** Acts a permutation on an Eldag.
         */

        Eldag act(const Eldag& eldag, Perm perm)
        {
            return act(eldag, perm.partition.make_perm(), perm.perms);
        }

        /** Left multiplies an automorphism onto a coset.
         */

        Coset left_mult(const Simple_perm& perm, const Coset& coset)
        {
            return { coset, perm };
        }

        /** Creates a transversal system for automorphisms.
         */

        auto create_transv(const Coset& upper, const Coset& lower)
        {
            return std::move(std::make_unique<Transv>(
                upper.size(), lower.get_individualized()));
        }
    };
} // End namespace libcanon::internal.

template <typename A, typename F>
std::tuple<Simple_perm, Perms<A>, Sims_transv<Perm<A>>> canon_eldag(
    const Eldag& eldag, const Symms<A>& symms, F init_colour)
{
    Partition init_part{ eldag.size() };
    init_part.split_by_key(0, init_colour);

    using Refiner = internal::Eldag_refiner<A>;
    Refiner refiner{};
    internal::Eldag_coset<A> root_coset{ init_part, symms };
    typename Refiner::Container container{};

    auto aut = add_all_candidates(refiner, eldag, root_coset, container);

    const auto& canon_form
        = std::min_element(container.begin(), container.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

    return { canon_form.second.make_perm(), std::move(canon_form.second.perms),
        std::move(aut) };
}

} // End namespace libcanon.

#endif // LIBCANON_ELDAG_H
