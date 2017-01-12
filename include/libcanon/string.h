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
        , curr_(prev->next())
        , selected_(selected)
        , perm_(curr->get_repr(selected))
        , next_(curr->next())
    {
        // We either have a coset representative or we have the implicit
        // identity.

        assert(perm_ || selected_ == curr->target);
    }

    /** Constructs a root coset for the group given as a Sims transversal.
     */

    Sims_coset(const Sims_transv<P>& group)
        : prev_(nullptr)
        , curr_(nullptr)
        , selected(0) // Not going to be used.
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
        for (const Sims_coset* i = &coset; !i->is_root(); i = curr->prev_) {
            if (i->perm_) {
                point = *i->perm >> point;
            } else {
                // Else we assume we chose identity.
                assert(i->selected() == i->curr_->target());
            }
        }

        return point;
    }

    /** Gets a permutation in the coset.
     */

    P get_a_perm() const
    {
        std::vector<const P*> perms{};

        for (const Sims_coset* i = this; !i->is_root(); i = coset->prev_) {
            perms.push_back(i->perm_); // The permutation can be null.
        }

        return chain(size(), perms.crbegin(), perms.crend());
    }

    /** Tests if two cosets are equal.
     *
     * For performance reason, here we just compare the selected point and
     * assume that all the rest matches.
     */

    bool operator==(const Sims_coset& other)
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

// clang-format off

/** The type for the alphabet of a string combinatorial structure.
 *
 * The type is set to the type obtained by indexing the structure by integral
 * values.
 *
 * Note that this may lead the problems for cases like a boolean vector.
 */

template <typename S>
using Alphabet_of<S> = std::remove_reference_t<std::result_of<
    decltype(&S::operator[])(size_t)
>>;

// clang-format on

/** Refiner for string canonicalization problem based on Sims transversal.
 *
 * For a problem, in addition to the type of permutation `P`, we also need the
 * type of the combinatorial structure `S` and the comparator for the alphabet
 * in the string `C`
 *
 * `S` needs to be indexable and able to be constructed from a size.  The
 * comparator needs to implement `equals` and `is_less` methods for comparing
 * the contents in the string.
 */

template <typename S, typename P> class Sims_refiner {
public:
    //
    // Types required by the refiner protocol.
    //

    using Coset = Sims_coset<P>;
    using Perm = P;
    using Transv = Sims_transv<P>;
    using Structure = S;
    using Container = std::unordered_map<S, P>;

    template <typename T>
    Sims_refiner(size_t size)
        : size{ size }
    {
    }

    //
    // Methods required by the refiner protocol
    //

    /** Tests if a given coset is a leaf coset. */
    bool is_leaf(const Coset& coset) const { return !coset.get_next(); }

    /** Refines a non-leaf coset */
    std::vector<Coset> refine(const Structure& obj, const Coset& coset) const
    {
        std::vector<Coset> children{};

        Point target = coset.get_next()->get_target();
        // Only valid for non-leaf cosets.

        // Find all points of minimum colour in the orbit.
        auto min = obj(coset >> target);
        children.emplace_back(coset, target);

        for (const auto& perm : coset.get_curr()) {
            Point src = perm >> target;
            Point orig = coset >> src;
            auto colour = obj(orig);
            if (colour == min) {
                children.emplace_back(coset, src);
            } else if (colour < min) {
                min = colour;
                children.clear();
                children.emplace_back(coset, src);
            }
            // Do nothing when a greater point is found.
        }

        return children;
    }

    /** Gets a permutation in a leaf coset */
    Perm get_a_perm(const Coset& coset) const
    {
        std::unique_ptr<P> perm = coset.get_a_perm();
        if (perm) {
            return std::move(*perm);
        } else {
            return Perm{ size }; // Create the identity permutation.
        }
    }

    /** Acts a given permutation to a combinatorial object. */
    Structure act(const Perm& perm, const Structure& obj) const
    {
        Structure result(size);

        for (size_t i = 0; i < size; ++i) {
            result[perm << i] = obj[i];
        }

        return result;
    }

    /** Left multiplies a coset by a permutation */
    Coset left_mult(const Perm& perm, const Coset& coset) const
    {
        return { coset.get_prev(), perm >> coset.get_selected() };
    }

    auto get_transv(const Coset& upper, const Coset& lower) const
    {
        return std::move(
            std::make_unique<Transv>(size, lower >> lower.get_target()));
    }

private:
    size_t size;
};

/** Canonicalize the given string.
 *
 * This is the main driver function for string canonicalization problem.  Note
 * that in addition to the other constrains, the combinatorial object needs to
 * satisfy the concept required in the template \ref String_structure.
 */

template <typename S, typename P>
std::pair<S, std::unique_ptr<Sims_transv<P>>> canon_string(
    const S& input, const Sims_transv<P>& group)
{
    using Refiner = Sims_refiner<S, P>;
    Refiner refiner{ input.size() };
    Sims_coset<P> whole_group(group);
    typename Refiner::Container container{};

    auto aut = add_all_candidates(refiner, input, whole_group, container);

    const auto& canon_form
        = std::min_element(container.begin(), container.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

    aut->conj(canon_form.second);
    return { std::move(canon_form.first), std::move(aut) };
}

} // End namespace libcanon

//
// Injection into the standard namespace for hashing.
//

namespace std {

/** Hasher for string combinatorial structures */

template <typename S> struct hash<libcanon::String_structure<S>> {
    size_t operator()(const libcanon::String_structure<S>& obj) const
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
