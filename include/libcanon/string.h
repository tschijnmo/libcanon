/** Canonicalization of strings.
 *
 * In this module, the generic canonicalization algorithm is applied to string
 * canonicalization problem.  Data types required for the generic
 * canonicalization algorithm to work are defined.  And the hashing of some of
 * them are injected into the `std` namespace for ease of use.  All these
 * culminates in the \ref canon_string function, which canonicalizes any given
 * string combinatorial object.
 *
 */

#ifndef LIBCANON_STRING_H
#define LIBCANON_STRING_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <libcanon/canon.h>
#include <libcanon/perm.h>
#include <libcanon/sims.h>

namespace libcanon {

/** Cosets for Sims transversal system.
 *
 * A left coset for one level of  Sims transversal system is basically a
 * pointer to a transversal system with a pointer to the selected point to be
 * moved to the target.  Here we also have a pointer to the previous level of
 * coset that the current coset is in.  The product of all these permutations
 * gives a permutation in the coset as in the whole group.  For graceful
 * treatment of the whole group as a coset, a pointer to the next level of
 * transversal is also included.
 */

template <typename P> class Sims_coset {
public:
    /** Constructs a Sims coset by expanding from a previous one.
     *
     * \param prev The previous level of Sims coset.
     * \param selected The point selected to be moved into the target.  Note
     * that it needs to be a point in the orbit.
     */

    Sims_coset(const Sims_coset& prev, Point selected)
        : prev_(&prev)
        , curr_(prev.next())
        , selected_(selected)
        , perm_(curr_->get_repr(selected))
        , next_(curr_->next())
    {
        // We either have a coset representative or we have the implicit
        // identity.

        assert(perm_ || selected_ == curr_->target());
    }

    /** Constructs a root coset for the group given as a Sims transversal.
     */

    Sims_coset(const Sims_transv<P>& group)
        : prev_(nullptr)
        , curr_(nullptr)
        , selected_(group.size()) // Not going to be used.
        , perm_(nullptr)
        , next_(&group)
    {
    }

    /** Gets the pre-image of a point by the permutation labelling this coset.
     *
     * This function is more suitable for the query of just a few points.  For
     * full construction of a permutation, use `get_a_perm`.
     */

    friend Point operator>>(const Sims_coset& coset, Point point)
    {
        for (const Sims_coset* i = &coset; !i->is_root(); i = i->prev_) {
            if (i->perm_) {
                point = *i->perm_ >> point;
            } else {
                // Else we assume we chose identity.
                assert(i->selected() == i->curr_->target());
            }
        }

        return point;
    }

    /** Gets the image of a point by a permutation in the coset.
     */

    friend Point operator<<(const Sims_coset& coset, Point point)
    {
        auto perms = coset.gather_perms();
        for (auto i = perms.crbegin(); i != perms.crend(); ++i) {
            point = **i << point;
        }

        return point;
    }

    /** Gets a permutation in the coset.
     */

    P get_a_perm() const
    {
        auto perms = gather_perms();

        return chain<P>(size(), perms.crbegin(), perms.crend());
    }

    /** Tests if two cosets are equal.
     *
     * For performance reason, here we just compare the selected point and
     * assume that all the rest matches.
     */

    bool operator==(const Sims_coset& other) const
    {
        assert(curr_ == other.curr_);
        return this->selected_ == other.selected_;
    }

    /** Computes the hash of the coset.
     *
     *  In the same vein as the equality comparison, here we just use the point
     *  moved into the target as the hash.
     */

    size_t hash() const { return selected_; }

    /** Gets the selected point to be put into the target.
     *
     * This method should *not* be called on root coset.
     */

    Point selected() const
    {
        assert(!is_root());
        return selected_;
    }

    /** Gets the size of the permutation domain.
     */

    size_t size() const
    {
        // Root are guaranteed to have next, others are guaranteed to have
        // current.
        return is_root() ? next_->size() : curr_->size();
    }

    /** Gets the previous coset.
     *
     * It can only be called on non-root coset.
     */

    const Sims_coset* prev() const
    {
        assert(!is_root());
        return prev_;
    }

    /** Gets pointer to the next level of transversal system.
     */

    const Sims_transv<P>* next() const { return next_; }

    /** Decides if a coset is the full group.
     */

    bool is_root() const { return !curr_; }

private:
    //
    // Utility methods
    //

    /** Gather all the permutation steps.
     *
     * The result will be a vector of points, with the last entry being the
     * first permutation.  Note that null values indicating implicit identity
     * will not be added.
     */

    std::vector<const P*> gather_perms() const
    {
        std::vector<const P*> perms{};

        for (const Sims_coset* i = this; !i->is_root(); i = i->prev_) {
            if (i->perm_)
                perms.push_back(i->perm_);
        }

        return perms;
    }

    //
    // Data fields
    //

    /** Pointer to the previous level of coset, nullptr for first level.
     */

    const Sims_coset* prev_;

    /** The current level of transversal.
     *
     * It will be set to null for the root coset.
     */

    const Sims_transv<P>* curr_;

    /** The label for the coset.
     */

    Point selected_;

    /** Reference to the permutation chosen by this level of coset.
     */

    const P* perm_;

    /** Pointer to the next level of subgroup.
     */

    const Sims_transv<P>* next_;
};

/** The type for the alphabet of a string combinatorial structure.
 *
 * The type is set to the type obtained by indexing the structure by integral
 * values.
 *
 * Note that this may lead the problems for cases like a boolean vector.
 */

template <typename S>
using Alphabet_of = std::decay_t<decltype(std::declval<S>()[(size_t)0])>;

/** The result of acting a permutation on a string combinatorial object.
 *
 * For convenience of implementation, here the result is set to a vector formed
 * from the alphabet of the corresponding combinatorial structure.  Actual
 * users can actually be isolated from this.
 */

template <typename S>
class Sims_act_res : public std::vector<Alphabet_of<S>>,
                     private std::hash<Alphabet_of<S>> {
public:
    /** Constructs an empty result.
     */

    Sims_act_res(size_t size)
        : std::vector<Alphabet_of<S>>(size)
        , std::hash<Alphabet_of<S>>()
    {
    }

    /** Computes the hash of the result.
     */

    size_t hash() const
    {
        size_t hash = 0;
        for (const auto& i : *this) {
            combine_hash(
                hash, static_cast<const std::hash<Alphabet_of<S>>&>(*this)(i));
        }

        return hash;
    }
};

/** Acts a permutation on a string combinatorial object.
 *
 * Here the result is not required to be the same as the given object.  Here
 * the result type is required to be constructible from a given size, and its
 * indexed result is assignable.
 */

template <typename R, typename P, typename S>
R act_string(const P& perm, const S& orig)
{
    size_t size = perm.size();
    R result(size);

    for (size_t i = 0; i < size; ++i) {
        result[i] = orig[perm >> i];
    }

    return result;
}

/** Container for the candidates
 *
 * Here we use a simple unordered map from the action result to permutation.
 * The hashing for the action result are injected into the std namespace.
 */

template <typename S, typename P>
using Sims_candidates = std::unordered_map<Sims_act_res<S>, P>;

/** Refiner for string canonicalization problem based on Sims transversal.
 *
 * For a problem, in addition to the type of permutation `P`, we also need the
 * type of the combinatorial structure `S`, which needs to be indexable to give
 * comparable and hashable alphabet.
 */

template <typename S, typename P> struct Sims_refiner {
public:
    //
    // Types required by the refiner protocol.
    //

    using Coset = Sims_coset<P>;
    using Structure = S;

    //
    // Methods required by the refiner protocol
    //

    /** Refines a non-leaf coset.
     *
     * Currently the refined cosets are all put into a vector eagerly.  It
     * might be better done lazily.
     */

    std::vector<Coset> refine(const S& obj, const Coset& coset) const
    {

        // Only valid for non-leaf cosets.
        assert(coset.next());
        const Sims_transv<P>& transv(*coset.next());
        Point target = transv.target();

        // Find all points of minimum colour in the orbit.
        std::vector<Point> points_w_min{};

        // First we take the target itself as the min colour.
        auto min = obj[coset >> target];
        points_w_min.push_back(target);

        // Then we loop over other points in the orbit and update.
        for (const auto& perm : transv) {
            Point src = perm >> target;
            Point orig = coset >> src;
            auto colour = obj[orig];

            if (colour == min) {
                points_w_min.push_back(src);
            } else if (colour < min) {
                min = colour;
                points_w_min.clear();
                points_w_min.push_back(src);
            }
            // Do nothing when a greater point is found.
        }

        std::vector<Coset> children{}; // Named return value.
        std::transform(points_w_min.cbegin(), points_w_min.cend(),
            std::back_inserter(children),
            [&](auto src) { return Sims_coset<P>(coset, src); });

        return children;
    }

    /** Tests if a given coset is a leaf coset.
     */

    bool is_leaf(const S& obj, const Coset& coset) const
    {
        return !coset.next();
    }

    /** Gets a permutation in a leaf coset.
     *
     * This should only be called on leaf cosets.
     */

    P get_a_perm(const Coset& coset) const
    {
        // It is hard to assert leaf here.
        return coset.get_a_perm();
    }

    /** Acts a given permutation to a combinatorial object.
     *
     * The result is an instance of the action result, not of the same type as
     * the given structure.
     */

    Sims_act_res<S> act(const P& perm, const S& obj) const
    {
        return act_string<Sims_act_res<S>>(perm, obj);
    }

    /** Left multiplies a coset by a permutation.
     */

    Coset left_mult(const P& perm, const Coset& coset) const
    {
        const auto& prev = *coset.prev();
        return { prev, prev << (perm >> (prev >> coset.selected())) };
    }

    /** Creates a transversal system.
     */

    auto create_transv(const Coset& upper, const Coset& lower) const
    {
        size_t size = upper.size();
        assert(lower.size() == size);

        return std::make_unique<Sims_transv<P>>(
            upper >> lower.selected(), size);
    }
};

/** Canonicalize the given string.
 *
 * This is the main driver function for string canonicalization problem.  Note
 * that the canonical form itself is not returned.  But rather, we return a
 * permutation bringing the combinatorial object to the canonical form and the
 * automorphism of the *canonicalized* object.  This can be understood as the
 * left coset of permutation canonicalizing the given object.
 *
 * The transversal system for the automorphism group is minimized.  Null
 * pointers indicates the absence of any symmetry.
 */

template <typename S, typename P>
std::pair<P, std::unique_ptr<Sims_transv<P>>> canon_string(
    const S& input, const Sims_transv<P>& group)
{
    using Refiner = Sims_refiner<S, P>;
    Refiner refiner{};
    Sims_coset<P> whole_group(group);
    Sims_candidates<S, P> candidates{};

    auto aut = add_all_candidates(refiner, input, whole_group, candidates);

    auto canon_form = std::min_element(candidates.begin(), candidates.end(),
        [](const auto& a, const auto& b) { return a.first < b.first; });

    P& canon_perm = canon_form->second;

    // Conjugate and minimize the automorphism group.
    aut->conj(canon_perm);
    auto min_aut = min_transv(std::move(aut));

    return { std::move(canon_perm), std::move(min_aut) };
}

} // End namespace libcanon

//
// Injection into the standard namespace for hashing.
//

namespace std {

/** Hasher for string combinatorial structures */

template <typename S> struct hash<libcanon::Sims_act_res<S>> {
    size_t operator()(const libcanon::Sims_act_res<S>& obj) const
    {
        return obj.hash();
    }
};

/** Hasher for Sims transversal cosets */

template <typename P> struct hash<libcanon::Sims_coset<P>> {
    size_t operator()(const libcanon::Sims_coset<P>& coset) const
    {
        return coset.hash();
    }
};
}

#endif // LIBCANON_STRING_H
