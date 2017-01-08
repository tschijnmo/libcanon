/** \file
 *
 * General utilities for permutation groups.
 *
 * This module contains some general data structure and manipulations for
 * permutation groups.  Most importantly, \ref Perm is a fairly generic data
 * type for permutations.  And \ref adapt_transv is a generic function to adapt
 * transversal systems.
 *
 */

#ifndef LIBCANON_PERM_H
#define LIBCANON_PERM_H

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

#include <libcanon/utils.h>

/** Libcanon public namespace.
 */

namespace libcanon {

//
// Permutation data type
// ---------------------
//

/** Type for a point to be permuted by a permutation.
 *
 * Here points are labelled simply by the first few natural numbers.  The
 * built-in `size_t` is used for easy indexing of the array for pre-image and
 * image arrays.
 */

using Point = size_t;

/** Base type for permutation expressions.
 *
 * This base type is for all kinds of expressions formed from permutations.
 * The type parameter should be set to the actual type.
 *
 * For permutations $g$, $h$ and a point $\alpha$, operation supported are
 *
 * - `g << alpha` for the image of $\alpha$ under this permutation.
 *
 * - `g >> alpha` for the pre-image of $\alpha$.
 *
 * These two notations are motivated by putting the resulted point to the left
 * of the expression, which is common in assignment expressions in imperative
 * programming languages.
 *
 * Also supported are `~` operator for inversion and `|` operator for
 * multiplication.  For the multiplication, it needs to be noted that
 * permutations always acts from left mathematically.  So we actually have
 *
 * ```
 * (p | q) << a == q << (p << a)
 * ```
 *
 *
 */

template <typename T> class Perm_expr {

public:
    /** Gets the pre-image of a point.
     */

    friend Point operator>>(const Perm_expr& perm, Point point)
    {
        return static_cast<const T&>(perm) >> point;
    }

    /** Gets the image of a point.
     */

    friend Point operator<<(const Perm_expr& perm, Point point)
    {
        return static_cast<const T&>(perm) << point;
    }

    /** Gets the accompanied action of a permutation.
     */

    auto acc() const { return static_cast<const T&>(*this).acc(); }

    /** Gets the size of the permutation domain.
     */

    size_t size() const { return static_cast<const T&>(*this).size(); }
};

/** Atomic permutation type.
 *
 * Here a permutation is stored redundantly by two arrays for both the preimage
 * of all points and the post image of all points.  Also a permutation can be
 * accompanied by an action, which is going to be composed by `^` for
 * multiplication and `|` with itself for inversion.  In this way, when the
 * accompanied action is a product of the $\mathbb{Z}_2$ group encoded as
 * integral types, nothing needs to be added.  Also the type for accompanied
 * action needs to be able to be constructed from integral zero for the
 * identity action.
 */

template <typename A> class Perm : public Perm_expr<Perm<A>> {
public:
    /** Constructs a permutation from a preimage array. */
    template <typename Input_it>
    Perm(Input_it begin, Input_it end, A acc_input = 0)
        : image{}
        , pre_image{}
        , acc{ acc_input }
    {
        // To avoid the ambiguity in vector constructor when the iterator type
        // is integral.
        std::copy(begin, end, std::back_inserter(pre_image));

        size = pre_image.size();
        set_image();
    }

    /** Creates an identity permutation */
    Perm(size_t size)
        : image(size)
        , pre_image(size)
        , acc{ 0 }
    {
        std::iota(image.begin(), image.end(), 0);
        std::iota(pre_image.begin(), pre_image.end(), 0);
    }

    //  Here we need to put some default explicitly so that they will not be
    //  override by the generic constructor from expressions.

    /** Copy-constructs a permutation. */
    Perm(const Perm& perm) = default;

    /** Move-constructs a permutation */
    Perm(Perm&& perm) = default;

    /** Constructs a permutation from an expression. */
    template <typename T>
    Perm(const Perm_expr<T>& expr)
        : size{ expr.size() }
        , image{ size }
        , pre_image{ size }
        , acc{ expr.get_acc() }
    {
        // Duplicate loops for cache friendliness.
        size_t i;

        for (i = 0; i < size; ++i) {
            image[i] = expr << i;
        }

        for (i = 0; i < size; ++i) {
            pre_image[i] = expr >> i;
        }
    }

    /** Constructs a permutation from a chain of permutations.
     *
     * Here the result will be the multiplication of all permutation in the
     * chain, which is given as an iterator over *pointers* to permutations.
     *
     * Compared with implementation based on expression templates, this
     * specialized function is a lot more cache-friendly.
     */

    template <typename Input_it>
    Perm(size_t size, Input_it begin, Input_it end)
        : size{ size }
        , image{}
        , pre_image{}
    {
        Point_vec pts1(size);
        Point_vec pts2(size);
        std::iota(pts1.begin(), pts1.end(), 0);

        Point_vec* src = &pts1;
        Point_vec* dest = &pts2;

        bool got_acc = false;

        std::for_each(begin, end, [&](auto element) {
            const Perm& perm{ *element };
            // Invalid iterator given if error happens here.

            if (got_acc) {
                acc = acc ^ perm.get_acc();
            } else {
                acc = perm.get_acc();
            }

            for (size_t i = 0; i < size; ++i) {
                (*dest)[perm << i] = (*src)[i];
            }

            std::swap(src, dest);
        });

        pre_image = std::move(*src);
        set_image();
    }

    //
    // Permutation expression operations.
    //

    /** Gets the pre-image of a point. */
    friend Point operator>>(const Perm& perm, Point point)
    {
        return perm.pre_image[point];
    }

    /** Gets the image of a point. */
    friend Point operator<<(const Perm& perm, Point point)
    {
        return perm.image[point];
    }

    /** Gets the accompanied action of a permutation. */
    auto get_acc() const { return acc; }

    /** Gets the size of the permutation domain */
    size_t get_size() const { return image.size(); }

    // Permutations are immutable.
    Perm& operator=(const Perm& perm) = delete;
    Perm& operator=(Perm&& perm) = delete;

    /** Gets the earliest point moved by the permutation.
     *
     * A point `size` will be returned if the permutation is identity.
     */

    Point get_earliest_moved() const
    {
        for (Point i = 0; i < size; ++i) {
            if (i != *this << i)
                return i;
        }
        return size;
    }

private:
    /** Sets the image array from the pre-image array */
    void set_image()
    {
        image.resize(size);
        for (size_t i = 0; i < size; i++) {
            image[pre_image[i]] = i;
        }
    }

    size_t size;
    using Point_vec = std::vector<Point>;
    Point_vec image;
    Point_vec pre_image;
    A acc;
};

/** Expression for inverted permutations.
 *
 * This is a lazy expression for inverted permutations.  A constant reference
 * is stored for its operand.
 */

template <typename T> class Inv_perm : public Perm_expr<Inv_perm<T>> {
public:
    Inv_perm(const T& operand)
        : operand{ operand }
    {
    }

    //
    // Permutation expression operations.
    //

    /** Gets the pre-image of a point. */
    friend Point operator>>(const Inv_perm& perm, Point point)
    {
        return perm.operand << point;
    }

    /** Gets the image of a point. */
    friend Point operator<<(const Inv_perm& perm, Point point)
    {
        return perm.operand >> point;
    }

    /** Gets the accompanied action of a permutation. */
    auto get_acc() const { return operand.get_acc() | operand.get_acc(); }

    /** Gets the size of the permutation domain */
    size_t get_size() const { return operand.size(); }

private:
    const T& operand;
};

/** Inverts a permutation expression.
 *
 * The familiar notation of `~` is used.  The result is a lazily evaluated
 * inversion expression.
 */

template <typename T> auto operator~(const Perm_expr<T>& expr)
{
    return Inv_perm<T>(expr);
}

/** Expression for the product of two permutation expressions.
 *
 * The two expressions does not have to be atomic.  Any permutation expression
 * works.
 */

template <typename T1, typename T2>
class Perm_prod : public Perm_expr<Perm_prod<T1, T2>> {
public:
    Perm_prod(const T1& left, const T2& right)
        : left{ left }
        , right{ right }
    {
    }

    //
    // Permutation expression operations.
    //

    /** Gets the pre-image of a point. */
    friend Point operator>>(const Perm_prod& perm, Point point)
    {
        return perm.left >> (perm.right >> point);
    }

    /** Gets the image of a point. */
    friend Point operator<<(const Perm_prod& perm, Point point)
    {
        return perm.left << (perm.right << point);
    }

    /** Gets the accompanied action of a permutation. */
    auto get_acc() const { return left.get_acc() ^ right.get_acc(); }

    /** Gets the size of the permutation domain
     *
     * No checking is done here.  Checking are performed during (before) the
     * construction.
     */

    size_t get_size() const { return left.size(); }
private:
    const T1& left;
    const T2& right;
};

/** Multiplies to permutation expression.
 *
 * The result is a lazily-evaluated permutation expression.  Note that the two
 * operands need to have the same size and accompanied action type.
 */

template <typename T1, typename T2>
auto operator|(const Perm_expr<T1>& left, const Perm_expr<T2>& right)
{
    static_assert(std::is_same<decltype(left.get_acc()),
                      decltype(right.get_acc())>::value,
        "Two operands need to have the same accompanied action type.");
    assert(left.get_size() == right.get_size()
        && "Two operands need act on the same domain");

    return Perm_prod<T1, T2>(left, right);
}

/** Type for simple permutations.
 *
 * This should be sufficient for most of the purposes, when the permutation are
 * either accompanied by no other action, and when the action is simply
 * isomorphic to the tensor product of a few $\mathbb{Z}_2$ groups.
 */

using Simple_perm = Perm<char>;

//
// Transversal adaptation
// ----------------------
//

#ifdef __cpp_concepts
// clang-format off

template <typename T>
concept bool Transv = requires () {
    typename T::Perm;
} && Simple_iterable<Transv, typename T::Perm> && requires (
    T transv, typename T::Perm perm) {
    { transv.next } -> std::unique_ptr<T>;
    { transv.has(perm) } -> bool;
    { transv.get_repr(perm) } -> const T::Perm*;
    { transv.insert(perm) };
};

// clang-format on
#endif

//
// Utilities for the transversal adaptation.
//

namespace internal {

    /**
     * Tests if a transversal is a leaf.
     *
     * A leaf transversal is one without a next transversal.
     */

    template <typename T> bool is_leaf_transv(const T& transv)
    {
        return transv.next == nullptr;
    }

    /** Process a permutation for a transversal.
     *
     * This is the core function for transversal system adaptation.  If the
     * given permutation is in the subgroup next level, it will be added to the
     * container for permutations to pass.  Or it will be added to the given
     * transversal container if no representative is already present for its
     * coset.
     */

    template <typename P, typename T, typename V>
    void proc_perm_for_transv(P&& perm, T& transv, V& perms_to_pass)
    {
        if (transv.next->has(perm)) {
            perms_to_pass.push_back(std::forward<P>(perm));
        } else {
            auto repr = transv.get_repr(perm);
            if (repr == nullptr) {
                transv.insert(std::forward<P>(perm));
            }
        }
    }

} // End namespace internal.

/** Adapts a transversal system into another.
 *
 * The input and output transversal systems are assumed to be for the same
 * group up to the same leaf.  Note that this function takes the ownership of
 * the input transversal and deallocates it after the adaptation.
 *
 * The transversals can be any data type satisfying the concept of a
 * transversal system `Transv`.  Most importantly,
 *
 * - A type `Perm` needs to be defined for the permutation type it contains.
 *
 * - It needs to be an iterable of permutations in the transversal system.
 *   Note that identity should be skipped in the iteration, since it is always
 *   present.
 *
 * - It needs to have an attribute named `next` for a unique pointer to the
 *   next level of subgroup.
 *
 * - It needs to support methods `has`, `get_repr`, and `insert`, to test if a
 *   permutation is in the subgroup of the transversal system, to get a
 *   constant pointer to a coset representation (NULL if not present), and to
 *   insert a permutation into the transversal.
 *
 */

template <typename T> void adapt_trasv(std::unique_ptr<T> input, T& output)
{
    // Gather all transversals to loop over them in reverse order.
    std::vector<T*> inputs{};
    T* curr = input.get();
    while (!internal::is_leaf_transv(*curr)) {
        inputs.push_back(curr);
        curr = curr->next.get();
    }

    for (auto curr_input = inputs.rbegin(); curr_input != input.rend();
         ++curr_input) {

        using Perm_vector = std::vector<typename T::Perm>;
        Perm_vector passed_perms{};
        for (auto& i : *curr_input) {
            passed_perms.push_back(std::move(i));
        }
        Perm_vector perms_to_pass{};

        T* curr_output = &output;
        while (!internal::is_leaf_transv(*curr_output)) {

            std::vector<const typename T::Perm*> existing{};
            std::transform(begin(*curr_output), end(*curr_output),
                std::back_inserter(existing), [](auto& perm) { return &perm; });

            for (const auto& passed_perm : passed_perms) {
                for (auto existing_perm : existing) {
                    internal::proc_perm_for_transv(passed_perm | *existing_perm,
                        *curr_output, perms_to_pass);
                }
                internal::proc_perm_for_transv(
                    std::move(passed_perm), *curr_output, perms_to_pass);
            }

            curr_output = curr_output->next.get();
            passed_perms.clear();
            passed_perms.swap(perms_to_pass);

        } // End loop output level.
    } // End loop input level.

    return;
}

} // End namespace libcanon.

#endif // LIBCANON_PERM_H
