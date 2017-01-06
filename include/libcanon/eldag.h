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
#include <libcanon/sims.h>
#include <libcanon/partition.h>
#include <libcanon/utils.h>

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

    /** Evaluates the hash of an Eldag.
     *
     * This is a very simple-minded hash function, just the values in the two
     * arrays are combined.
     */

    size_t hash() const
    {
        size_t seed = 0;
        for (size_t i : edges) {
            combine_hash(seed, i);
        }
        for (size_t i : ia) {
            combine_hash(seed, i);
        }
        return seed;
    }
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

        size_t hash() const { return get_individualized(); }

    private:
        /** The orbit label of valences of a node.
         */

        using Orbit = std::vector<Point>;

        /** An edge in detail.
         *
         * It is a sorted list of all the orbit labels of all the edges from
         * the parent node to the children.
         */

        using Detailed_edge = std::vector<size_t>;

        /** Detailed description of the in/out edges between a node and a cell.
         *
         * The result should be sorted for set semantics.  It is implemented as
         * a struct here to have better control over the comparison order among
         * them.
         */

        class Detailed_edges : public std::vector<Detailed_edge> {

            /** Compares two for order.
             *
             * Here larger number of connections is considered smaller.  In
             * this way, when we sweep the cells from beginning to end, nodes
             * with more connection with nodes of earlier colour will be put
             * earlier.
             */

            bool operator<(const Detailed_edges& other) const
            {
                return this->size() < other.size() || *this < other;
            }
        };

        /** Detailed description of the connection between a node and a cell.
         *
         * Here we have information about both the in and the out edges.
         * Default lexicographical order will be used.
         *
         * Normally for Eldag, either the in edges or the out edges are empty.
         */

        using Conns = std::pair<Detailed_edges, Detailed_edges>;

        /** Refines the currently holding partition and symmetry.
         *
         * Most of the actual work of the canonicalization of Eldag actually
         * comes here.
         */

        void refine(const Eldag& eldag)
        {
            // Refine until fixed point.

            while (true) {

                // First refine the automorphism at each node.
                refine_nodes(eldag);

                //
                // Unary split.
                //

                auto orbits = form_orbits();
                bool split = false;
                // Back up the current partition for unary split.
                std::vector<size_t> curr_partition(
                    partition.begin(), partition.end());

                for (auto i : curr_partition) {
                    split |= partition.split_by_key(
                        i, [&](auto point) -> const Orbit& {
                            return orbits[point];
                        });
                }
                if (split) {
                    // Binary split is expansive.  Refine as much as possible
                    // before carrying it out.
                    continue;
                }

                //
                // Binary split.
                //

                split = false;
                std::vector<Conns> conns(partition.size());

                // Here for each loop, we always take advantage of the latest
                // refine due to the special semantics of looping over cells in
                // a partition.

                for (auto splittee : partition) {
                    for (auto splitter : partition) {

                        update_conns4cell(
                            conns, eldag, orbits, splittee, splitter);

                        split |= partition.split_by_key(
                            splittee, [&](auto point) -> const Conns& {
                                return conns[point];
                            });
                    }
                }

                if (split) {
                    continue;
                } else {
                    // Now we reached a fixed point.
                    break;
                }
            };
        }

        /** The class representing the valences of a node.
         *
         * This class is designed to be directly interoperable with the string
         * canonicalization facilities.
         */

        class Valence {
            Valence(const Eldag& eldag, const Partition& partition, Point node,
                const Perm<A>& perm)
                : eldag{ eldag }
                , partition{ partition }
                , node{ node }
                , perm{ perm }
            {
            }

            /** Gets the colour of the node connected to a given slot.
             */

            Point operator[](size_t slot) const
            {
                size_t offset = perm >> slot;
                return partition.get_colour(
                    eldag.edges[eldag.ia[node + offset]]);
            }

        private:
            const Eldag& eldag;
            const Partition& partition;
            Point node;
            const Perm<A>& perm;
        };

        /** Refines the automorphism of all nodes.
         *
         * In this function, the permutations and symmetries for each node will
         * be refined as much as possible based on the current partition of the
         * nodes.
         */

        void refine_nodes(const Eldag& eldag)
        {
            // Refinement for nodes are independent.
            for (size_t i = 0; i < partition.size(); ++i) {

                auto curr_symm = symms[i];
                if (!curr_symm) {
                    // Nodes without symmetry.
                    continue;
                }

                Valence valence{ eldag, i, perms[i] };

                auto canon_res = canon_string(valence, *curr_symm);

                refined_perms[i]
                    = std::make_unique<Perm<A>>(std::move(canon_res.first));
                refined_symms[i] = std::move(canon_res.second);
                perms[i] = refined_perms[i].get();
                symms[i] = refined_symms[i].get();
            }
        }

        /** Forms the orbit label array for all nodes.
         */

        std::vector<Orbit> form_orbits() const
        {
            std::vector<Orbit> res{};

            std::transform(symms.begin(), symms.end(), std::back_inserter(res),
                [&](const auto aut) { return aut.get_orbits(); });
            return res;
        }

        /** Updates the connection information.
         */

        void update_conns4cell(std::vector<Conns>& conns, const Eldag& eldag,
            const std::vector<Orbit>& orbits, Point splittee, Point splitter)
        {
            // Short-cut singleton partitions.  They will be automatically
            // short-cut in split by key function any way.
            if (partition.get_cell_size(splittee) < 2)
                return;

            std::for_each(partition.cell_begin(splittee),
                partition.cell_end(splittee), [&](Point base) {
                    Detailed_edges& ins = conns[base].first;
                    Detailed_edges& outs = conns[base].second;

                    ins.clear();
                    outs.clear();

                    std::for_each(partition.cell_begin(splitter),
                        partition.cell_end(splitter), [&](Point curr) {
                            // Add connection with the current point.
                            add_detailed_edges(ins, eldag, orbits, curr, base);
                            add_detailed_edges(outs, eldag, orbits, base, curr);
                        });

                    std::sort(ins.begin(), ins.end());
                    std::sort(outs.begin(), outs.end());
                });
        }

        /** Adds detailed edge from a point to another if there is any.
         *
         * Empty connection will not be added.
         */

        void add_detailed_edges(Detailed_edges& dest, const Eldag& eldag,
            const std::vector<Orbit>& orbits, Point from, Point to)
        {
            Detailed_edge edge{};
            for (size_t i = eldag.ia[from]; i < eldag.ia[from + 1]; ++i) {
                if (eldag.edges[i] == to) {
                    Point index = perms[from] ? *perms[from] << i : i;
                    edge.push_back(orbits[from][index]);
                }
            }
            if (edge.size() == 0)
                return;
            std::sort(edge.begin(), edge.end());
            dest.push_back(std::move(edge));

            return;
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

//
// Std namespace injection for default hashing.
//

namespace std {

/** Hasher for Eldag.
 */

template <> struct hash<libcanon::Eldag> {
    size_t operator()(const libcanon::Eldag& eldag) const
    {
        return eldag.hash();
    }
};
}

#endif // LIBCANON_ELDAG_H
