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

/** Base type for string combinatorial structures.
 *
 * A string combinatorial structure need to,
 *
 * 1. Have method `size` to get the size of the structure.
 *
 * 2. Able to be constructed from a size.
 *
 * 3. Indexable by `[]`, to get the actual alphabet at the given index, also
 * indexable by `()`, to get the colour of the alphabet at the given index.
 * Note that the colour will be used for equality comparison by `==` and order
 * comparison by `<`, and hashing by std::hash.
 *
 * Then by sub classing the specialization of this template with the actual
 * type will make the equality and lexicographical comparison based on the
 * actual colour available.  Also available are the iterators for colours and
 * hashing.
 */

template <typename S> class String_structure {
public:
    /** The actual data type for string structures. */

    using Structure = S;
    using Colour = std::remove_reference_t<std::result_of_t<S(Point)>>;
    using Alphabet = std::remove_reference_t<decltype(std::declval<S>()[0])>;

    /** Gets the size of the string */

    size_t size() const { return static_cast<const S&>(*this).size(); }

    /** Gets the alphabet at the given index */

    auto& operator[](Point idx) const
    {
        return static_cast<const S&>(*this)[idx];
    }

    /** Gets the colour of the alphabet at the given index */

    auto operator()(Point idx) const
    {
        return static_cast<const S&>(*this)(idx);
    }

    /** Lexicographical quality comparison based on colour. */

    bool operator==(const String_structure& other) const
    {
        return std::equal(this->begin_colour(), this->end_colour(),
            other.begin_colour(), other.end_colour());
    }

    /** Lexicographical ordering based on colour */

    bool operator<(const String_structure& other) const
    {
        return std::lexicographical_compare(this->begin_colour(),
            this->end_colour(), other.begin_colour(), other.end_colour());
    }

    size_t hash() const
    {
        size_t seed = 0;
        std::hash<Colour> hasher{};
        std::for_each(this->begin_colour(), this->end_colour(),
            [&](auto colour) { combine_hash(seed, hasher(colour)); });
        return seed;
    }

    //
    // Easy iteration over the colours at each index.
    //

    /** The iterator class for colour iteration */

    class Colour_it {

    public:
        /** Creates a colour iterator. */

        Colour_it(const S& obj, Point idx)
            : obj{ obj }
            , idx{ idx }
        {
        }

        /** Increments the iterator. */

        Colour_it& operator++()
        {
            ++idx;
            return *this;
        }

        /** Equality comparison.
         *
         * Here for performance reasons, we just compare the index.
         */

        bool operator==(const Colour_it& other) { return idx == other.idx; }

        /** Inequality comparison. */

        bool operator!=(const Colour_it& other)
        {
            return !(this->operator==(other));
        }

        /** Dereferencing. */

        auto operator*() const { return obj(idx); }

    private:
        const S& obj;
        Point idx;
    };

    Colour_it begin_colour() const
    {
        return { static_cast<const S&>(*this), 0 };
    }

    Colour_it end_colour() const
    {
        return { static_cast<const S&>(*this), this->size() };
    }
};

/** Cosets for Sims transversal system.
 *
 * A left coset for one level of  Sims transversal system is basically a
 * pointer to a permutation with a reference to the next level of Sims
 * transversal.  Here we also have a pointer to the previous level of coset
 * that the current coset is in.  The product of all these permutations gives a
 * permutation in the actual coset.
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
        : prev{ *prev }
        , curr{ prev->next }
        , selected{ selected }
        , perm{ curr->get_repr(selected) }
        , next{ curr->next->get() }
    {
    }

    /**
     * Constructs a root coset for the group given as a Sims transversal.
     */

    Sims_coset(const Sims_transv<P>& group)
        : curr{ nullptr }
        , next{ *group }
    {
    }

    /** Gets the pre-image of a point by the permutation labelling this coset.
     *
     * This function is more suitable for the query of just a few points.  For
     * full construction of a permutation, use `get_a_perm`.
     */

    friend Point operator>>(const Sims_coset* coset, Point point)
    {
        for (; coset->curr != nullptr; coset = coset->prev) {
            if (coset->perm) {
                point = *coset->perm >> point;
            } // Else we assume we chose identity.
        }
        return point;
    }

    /**
     * Gets a permutation in the coset.
     *
     * Null pointer for identity permutation.
     */

    std::unique_ptr<P> get_a_perm() const
    {
        std::unique_ptr<P> res;
        std::unique_ptr<P> prod;

        for (const Sims_coset* coset = this; coset->curr != nullptr;
             coset = coset->prev) {
            if (coset->perm) {
                if (res) {
                    prod = std::make_unique<P>(*coset->perm | *res);
                    res.swap(prod);
                } else {
                    res = std::make_unique<P>(*coset->perm); // Copy.
                }
            }
        }
        return std::move(res);
    }

    /** Tests if two cosets are equal.
     *
     * For performance reason, here we just compare the selected point and
     * assume that all the rest matches.
     */

    bool operator==(const Sims_coset& other)
    {
        return this->selected == other.selected;
    }

    /** Computes the hash of the coset.
     *
     *  In the same vein as the equality comparison, here we just use the point
     *  moved into the target as the hash.
     */

    size_t hash() const { return selected; }

    /** Gets the selected point to be put into the target.
     *
     * Result for root coset is undefined.
     */

    Point get_selected() const { return selected; }

    /** Gets the previous coset.
     *
     * It can only be called on non-root coset.
     */

    const Sims_coset& get_prev() const { return *prev; }

    /** Gets the current transversal system.
     *
     * It can only be called on non-root coset.
     */

    const Sims_transv<P>& get_curr() const { return *curr; }

    /** Gets the next level of transversal system */

    const Sims_transv<P>* get_next() const { return next; }

private:
    /** Pointer to the previous level of coset, nullptr for first level. */
    const Sims_coset* prev;

    /** The current level of transversal. */
    const Sims_transv<P>* curr;

    /** The label for the coset. */
    Point selected;

    /** Reference to the permutation chosen by this level of coset. */
    const P* perm;

    /** Pointer to the next level of subgroup. */
    const Sims_transv<P>* next;
};

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
