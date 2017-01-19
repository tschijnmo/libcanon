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
#include <cassert>
#include <iterator>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include <libcanon/partition.h>
#include <libcanon/perm.h>
#include <libcanon/sims.h>
#include <libcanon/string.h>
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

    /** Gets the number of nodes in the Eldag.
     */

    size_t size() const { return ia.size() - 1; }

    /** Gets the number of edges in the Eldag.
     */

    size_t n_edges() const { return edges.size(); }

    /** Gets the number of valences of a node.
     */

    size_t n_valences(Point node) const
    {
        assert(node >= 0 && node < size());
        return ia[node + 1] - ia[node];
    }

    /** Creates an Eldag of a given size.
     *
     * Note that here the vectors are only set to the correct size without any
     * valid content.
     */

    Eldag(size_t n_nodes, size_t n_edges)
        : edges(n_edges)
        , ia(n_nodes + 1)
    {
    }

    /** Creates an empty Eldag.
     *
     * This empty Eldag is not completely empty, a zero is already added to the
     * ia array.
     */

    Eldag()
        : edges()
        , ia{ 0 }
    {
    }

    /** Copy constructs an eldag.
     */

    Eldag(const Eldag& eldag) = default;

    /** Move constructs an eldag.
     */

    Eldag(Eldag&& eldag) = default;

    /** Pushes the current number of edges into ia.
     *
     * This method can be helpful after finish adding edges from a node.
     */

    void update_ia() { ia.push_back(edges.size()); }

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
 * This type is going to be mainly used as input for the symmetries of nodes in
 * the public driver functions for Eldag canonicalization.
 */

template <typename P> using Node_symms = std::vector<const Sims_transv<P>*>;

/** The data type for the permutations to be applied to nodes in an Eldag.
 *
 * A vector of unique pointer to the given permutation type.  This is going to
 * be mainly in the return value of the public driver function for Eldag
 * canonicalization.  Note that here the value owns the reference to the
 * permutations.
 */

template <typename P> using Node_perms = std::vector<std::unique_ptr<P>>;

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
    Eldag res(size, eldag.n_edges());
    res.ia.clear();

    // The index that the next connection is going to be written to.
    size_t curr_idx = 0;
    res.ia.push_back(curr_idx);

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

            assert(offset >= 0 && offset < valence);

            Point conn = gl_perm << eldag.edges[prev_base + offset];
            res.edges[curr_idx] = conn;
            ++curr_idx;
        }

        res.ia.push_back(curr_idx);
    }

    return res;
}

namespace internal {

    /** Data type for owned references to symmetries.
     */

    template <typename P>
    using Owned_symms = std::vector<std::unique_ptr<Sims_transv<P>>>;

    /** Data type for borrowed references to permutations.
     */

    template <typename P> using Borrowed_perms = std::vector<const P*>;

} // End namespace internal

/** Data type for a coset in the canonicalization of an Eldag.
 *
 * Basically, it has a global partition of the nodes, the symmetries and
 * applied permutations for each of the nodes, which can be from the parent
 * coset or from the refined permutation and symmetries owned by itself.
 */

template <typename P> class Eldag_coset {
public:
    /** Constructs a coset based on initial partition and symmetries.
     *
     * This is useful for the creation of root coset.
     *
     * The templating here is just for perfect forwarding.  The initial
     * partition is assumed to be a Partition object.
     */

    template <typename T>
    Eldag_coset(T&& init_part, const Node_symms<P>& init_symms)
        : partition_(std::forward<T>(init_part))
        , symms_(init_symms)
        , perms_(init_part.size(), nullptr)
        , refined_symms_()
        , refined_perms_()
        , individualized_(init_part.size())
    {
        assert(init_part.size() == init_symms.size());
    }

    /** Constructs a coset by individualizing a point.
     *
     * This is useful for the individualization step.  Note that this
     * construction is better used on cosets after refinement.
     */

    Eldag_coset(const Eldag_coset& base, Point point)
        : individualized_(point)
        , partition_(base.partition_, point)
        , symms_(base.symms_)
        , perms_(base.perms_)
        , refined_perms_()
        , refined_symms_()
    {
    }

    /** Constructs a coset by acting a permutation.
     *
     * This is useful for the pruning of the refinement tree.  Note that the
     * created coset is not actually fully valid one.  Just the individualized
     * point is correctly set for the sake of removing the actual coset in a
     * hash set/map.
     */

    Eldag_coset(const Eldag_coset& base, const Simple_perm& perm)
        : individualized_(perm >> base.individualized())
        , partition_()
        , perms_()
        , symms_()
        , refined_perms_()
        , refined_symms_()
    {
    }

    /** Individualizes the refined partition.
     *
     * This method can be called after the testing of leaf, where the partition
     * is going to be refined.
     *
     * Currently we eagerly put the result in a vector, it should better be
     * done lazily.
     */

    std::vector<Eldag_coset> individualize_first_nonsingleton()
    {
        auto cell = std::find_if(partition_.begin(), partition_.end(),
            [&](auto cell) { return partition_.get_cell_size(cell) > 1; });

        // When used correctly, we should find a non-singleton cell here.
        assert(cell != partition_.end());

        std::vector<Eldag_coset> res;
        std::transform(partition.cell_begin(*cell), partition.cell_end(*cell),
            std::back_inserter(res),
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
        refine(); // Core work.
        return partition_.is_discrete();
    }

    /** Gets the point individualized during its construction.
     */

    Point individualized() const { return individualized_; }

    /** Gets the size of the graph.
     */

    size_t size() const { return partition_.size(); }

    /** Compares two cosets for equality.
     *
     * Here only the point that is individualized is compared.  So it is
     * actually only applicable to peers.
     */

    bool operator==(const Eldag_coset& other) const
    {
        return individualized() == other.individualized();
    }

    /** Evaluates the hash of a coset.
     *
     * Note that only the individualized point is used.  So it is only
     * useful for retrieval amongst it peers.
     */

    size_t hash() const { return individualized(); }

private:
    //
    // Core refinement facilities
    //

    // Data types for node distinguishing.

    /** The orbit label of valences of a node.
     *
     * Values of this type can be used to map a valence of a node to its orbit
     * label.
     */

    using Orbit = std::vector<Point>;

    /** The orbits of nodes in an Eldag.
     *
     * The orbit vector for any node should be able to be queried by indexing
     * it.
     */

    using Orbits = std::vector<Orbit>;

    /** An edge in detail.
     *
     * It is a sorted list of all the orbit labels of all the edges from
     * the parent node to the children.
     */

    using Detailed_edge = std::vector<size_t>;

    /** Detailed description of the in/out edges between a node and a cell.
     *
     * The result should be sorted for set semantics.  It is implemented as a
     * struct here to have better control over the comparison order among them.
     * Since the values are mostly only used for lexicographical comparison,
     * without need to be queried for a specific edge, here it is based on
     * vector rather than multi-set for better cache friendliness with the same
     * complexity.
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
            return this->size() > other.size() || *this < other;
        }
    };

    /** Detailed description of the connection between a node and a cell.
     *
     * Here we have information about both the in and the out edges.  Default
     * lexicographical order will be used.
     *
     * Normally for Eldag, either the in edges or the out edges are just empty.
     */

    using Conn = std::pair<Detailed_edges, Detailed_edges>;

    /** The connection information of all nodes.
     *
     * This vector is designed to be directly indexed by node labels for the
     * connection description.  Note that normally only the entries for nodes
     * in the same cell are comparable.
     */

    using Conns = std::vector<Conn>;

    // Driver function.

    /** Refines the global partition and local symmetries.
     *
     * Most of the actual work of the canonicalization of Eldag actually
     * comes here.
     */

    void refine(const Eldag& eldag)
    {
        // Refine until fixed point.

        while (true) {

            // First refine the automorphism at each node and form the orbits.
            refine_nodes(eldag);
            Orbits orbits = form_orbits(eldag);

            //
            // Unary split based on orbits.
            //

            bool split = false;
            // Back up the current partition for unary split.
            std::vector<size_t> curr_partition(
                partition_.begin(), partition_.end());

            for (auto i : curr_partition) {
                split |= partition_.split_by_key(i,
                    [&](auto point) -> const Orbit& { return orbits[point]; });
            }
            if (split) {
                // Binary split is expensive.  Refine as much as possible
                // before carrying it out.
                continue;
            }

            // Now all the nodes in the same cell must have identical orbit
            // labelling.  So that their comparison by connection labels start
            // to make sense.

            //
            // Binary split.
            //

            split = false;
            Conns conns(size());

            // Here for each loop, we always take advantage of the latest
            // refine due to the special semantics of looping over cells in a
            // partition.

            std::for_each(
                partition_.rbegin(), partition_.rend(), [&](auto splittee) {
                    for (auto splitter : partition) {

                        update_conns4cell(
                            conns, eldag, orbits, splittee, splitter);

                        split |= partition.split_by_key(
                            splittee, [&](auto point) -> const Conn& {
                                return conns[point];
                            });
                    }
                });

            if (split) {
                continue;
            } else {
                // Now we reached a fixed point.
                break;
            }
        }; // End main loop.
    }

    // Local refinement of the symmetries.

    /** The class representing the valences of a node.
     *
     * This class is designed to be directly interoperable with the string
     * canonicalization facilities.
     */

    class Valence {
        Valence(const Eldag& eldag, const Partition& partition, Point node,
            const P* perm)
            : eldag_(&eldag)
            , partition_(&partition)
            , node_(node)
            , perm_(perm)
        {
            assert(eldag.size() == partition.size());
            assert(0 <= node && node < partition.size());
        }

        /** Gets the colour of the node connected to a given slot.
         */

        Point operator[](size_t slot) const
        {
            size_t offset = slot;
            if (perm_) {
                offset = *perm_ >> slot;
            }

            auto base_idx = eldag_.ia[node_];
            assert(offset >= 0 && offset < eldag_.n_valences(node_));

            return partition_.get_colour(eldag_->edges[base_idx + offset]);
        }

    private:
        /** The Eldag that the valence is for.
         */

        const Eldag* eldag_;

        /** The current partition of the nodes in the Eldag.
         */

        const Partition* partition_;

        /** The node label for this.
         */

        Point node_;

        /** The current permutation of the node.
         */

        const P* perm_;
    };

    /** Refines the automorphism of all nodes.
     *
     * In this function, the permutations and symmetries for each node will be
     * refined as much as possible based on the current partition of the nodes.
     */

    void refine_nodes(const Eldag& eldag)
    {
        // Refinement for nodes are independent.
        for (size_t i = 0; i < size(); ++i) {

            auto curr_symm = symms[i];
            if (!curr_symm) {
                // Nodes without symmetry, without need for *further*
                // refinement.
                continue;
            }

            Valence valence(eldag, partition_, i, perms_[i]);

            auto canon_res = canon_string(valence, *curr_symm);

            if (!perms_[i]) {
                // When nothing has ever been applied to the node.

                refined_perms_[i]
                    = std::make_unique<P>(std::move(canon_res.first));
            } else {
                // Just to be safe about the order of evaluation.
                auto new_perm
                    = std::make_unique<P>(*perms_[i] | canon_res.first);
                refined_perms_[i] = std::move(new_perm);
            }
            refined_symms[i] = std::move(canon_res.second);
            perms[i] = refined_perms[i].get();
            symms[i] = refined_symms[i].get();
        }
    }

    /** Forms the orbit label array for all nodes.
     */

    Orbits form_orbits(const Eldag& eldag) const
    {
        Orbits res{};

        for (size_t node = 0; node < eldag.size(); ++node) {
            size_t n_valences = eldag.n_valences(node);

            Orbit orbit(n_valences);
            std::iota(orbit.begin(), orbit.end(), 0);

            for (const Sims_transv<P>* transv = symms[node]; transv != nullptr;
                 transv = transv->get_next()) {
                Point target = transv->target();
                for (const auto& perm : *transv) {
                    orbit[perm >> target] = orbit[target];
                }
            }

            res.push_back(std::move(orbit));
        }

        return res;
    }

    // Binary refinements.

    /** Updates the connection information.
     */

    void update_conns4cell(Conns& conns, const Eldag& eldag,
        const Orbits& orbits, Point splittee, Point splitter)
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
        const Orbits& orbits, Point from, Point to)
    {
        Detailed_edge edge{};
        size_t base_idx = eldag.ia[from];
        size_t n_valences = eldag.n_valences(from);

        for (size_t i = 0; i < n_valences; ++i) {
            Point conn_node = eldag.edges[base_idx + i];

            if (conn_node == to) {
                size_t index = perms_[from] ? *perms_[from] << i : i;
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

    Partition partition_;

    /** The current permutations applied to the nodes.
     */

    internal::Borrowed_perms<P> perms_;

    /** The current symmetries for each of the nodes.
     */

    Node_symms<P> symms_;

    /** The permutation on each node after refinement.
     */

    Node_Perms<P> refined_perms_;

    /** Refined symmetries for each node.
     */

    internal::Owned_symms<P> refined_symms_;

    /**
     * The point that is individualized in the construction of this coset.
     */

    Point individualized_;
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

template <typename P> struct Eldag_perm {

    /** The permutations on each of the nodes. */

    Perms<P> perms;

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

        return std::move(pre_imgs);
    }

    /** Gets the pre-image of a given point.
     */

    Point operator>>(Point point) const { return partition >> point; }

    /** Gets the permutations for each node.
     *
     * The returned vector will have borrowed reference for the permutations.
     */

    internal::Borrowed_perms<P> get_perms() const
    {
        internal::Borrowed_perms res{};

        for (const auto& i : perms) {
            res.push_back(i.get());
        }

        return res;
    }
};

/** Acts an eldag with a given permutation.
 *
 * This is a convenience wrapper function for acting the permutations
 * encapsulated in an Eldag_perm object onto an Eldag.
 */

template <typename P> Eldag act(const Eldag& eldag, const Eldag_perm<P>& perm)
{
    return act(eldag, perm.partition.make_perm(), perm.get_perms());
}

/** The actual refiner for Eldag canonicalization.
 */

template <typename P> class Eldag_refiner {
public:
    //
    // Types required by the refiner protocol.
    //

    using Coset = Eldag_coset<P>;
    using Structure = Eldag;

    //
    // Types not required by the refiner protocol.
    //
    // They are included here just for convenience.

    using Perm = Eldag_perm<P>;

    /** Refines the given coset.
     */

    std::vector<Coset> refine(const Eldag& eldag, const Coset& coset)
    {
        return coset.individualize_first_nonsingleton();
    }

    /** Decides if a coset is a leaf.
     */

    bool is_leaf(const Eldag& eldag, Coset& coset)
    {
        return coset.is_leaf(eldag);
    }

    /** Gets a permutation from a coset.
     */

    Perm get_a_perm(const Coset& coset) { return coset.get_a_perm(); }

    /** Acts a permutation on an Eldag.
     */

    Eldag act(const Perm& perm, const Eldag& eldag) { return act(eldag, perm); }

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
        assert(upper.size() == lower.size());

        return std::make_unique<Sims_transv<P>>(
            upper.size(), lower.individualized());
    }
};

/** Canonicalizes the given Eldag.
 *
 * Similar to the case of string canonicalization, here the permutation bring
 * the Eldag into canonical form is returned.  But the automorphism group
 * returned is with respect to the original graph rather than the canonical
 * form.
 */

template <typename P, typename F>
std::pair<Eldag_perm<P>, Sims_transv<P>> canon_eldag(
    const Eldag& eldag, const Symms<P>& symms, F init_colour)
{
    Partition init_part(eldag.size());
    init_part.split_by_key(0, init_colour);

    Eldag_refiner<P> refiner{};
    Eldag_coset<P> root_coset(init_part, symms);

    using Container = std::unorderd_map<Eldag, Eldag_perm<P>>;
    Container container{};

    auto aut = add_all_candidates(refiner, eldag, root_coset, container);

    const auto& canon_form
        = std::min_element(container.begin(), container.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

    return { std::move(canon_form.second), std::move(aut) };
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

/** Hash for Eldag coset.
 */

template <typename P> struct hash<libcanon::Eldag_coset<P>> {
    size_t operator()(const libcanon::Eldag_coset<P>& coset) const
    {
        return coset.hash();
    }
};

} // End namespace std

#endif // LIBCANON_ELDAG_H
