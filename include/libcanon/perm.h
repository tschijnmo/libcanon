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

/** The type for vectors of points.
 *
 * Since the points are identified by integral values, this vector can be used
 * as a mapping from points to points, which is useful for things like image
 * arrays and pre-image arrays.
 */

using Point_vec = std::vector<Point>;

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
 * of all points and the image of all points.  Also a permutation can be
 * accompanied by an action, which is going to be composed by `^` for
 * multiplication and `|` with itself for inversion.  In this way, when the
 * accompanied action is a product of the $\mathbb{Z}_2$ group encoded as
 * integral types, nothing needs to be added.  Also the type for accompanied
 * action needs to be able to be constructed from integral zero for the
 * identity action.
 */

template <typename A> class Perm : public Perm_expr<Perm<A>> {
public:
    /** Constructs a permutation from a preimage array.
     */

    template <typename Input_it>
    Perm(Input_it begin, Input_it end, A acc = 0)
        : images_()
        , pre_images_()
        , acc_(acc)
    {
        // To avoid the ambiguity in vector constructor when the iterator type
        // is integral.
        std::copy(begin, end, std::back_inserter(pre_images_));
        set_images();
    }

    /** Creates an identity permutation.
     */

    Perm(size_t size)
        : images_(size)
        , pre_images_(size)
        , acc_(0)
    {
        std::iota(images_.begin(), images_.end(), 0);
        std::iota(pre_images_.begin(), pre_images_.end(), 0);
    }

    /** Creates an empty permutation
     *
     * This constructor might be useful for cases where a placeholder is
     * needed.  Its values does not actually constitute a valid permutation of
     * any finite domain and their usage in any way might lead to undefined
     * behaviour.
     */

    Perm()
        : images_()
        , pre_images_()
        , acc_(0)
    {
    }

    //  Here we need to put some default explicitly so that they will not be
    //  override by the generic constructor from expressions.

    /** Copy-constructs a permutation.
     */

    Perm(const Perm& perm) = default;

    /** Move-constructs a permutation.
     */

    Perm(Perm&& perm) = default;

    /** Constructs a permutation from an expression.
     */

    template <typename T>
    Perm(const Perm_expr<T>& expr)
        : images_(expr.size())
        , pre_images_(expr.size())
        , acc_(expr.acc())
    {
        // Duplicate loops for cache friendliness.
        size_t i;
        size_t size = expr.size();

        for (i = 0; i < size; ++i) {
            images_[i] = expr << i;
        }

        for (i = 0; i < size; ++i) {
            pre_images_[i] = expr >> i;
        }
    }

    /** Constructs a permutation from a chain of permutations.
     *
     * Here the result will be the multiplication of all permutation in the
     * chain, which is given as an iterator over *pointers* to permutations.
     *
     * Compared with implementation based on expression templates, this
     * specialized function is a lot more cache-friendly.
     *
     * When an empty iterator is given, the result will be the identity
     * permutation.
     */

    template <typename Input_it>
    Perm(size_t size, Input_it begin, Input_it end)
        : images_()
        , pre_images_()
        , acc_(0)
    {
        Point_vec pts1(size);
        Point_vec pts2(size);
        std::iota(pts1.begin(), pts1.end(), 0);

        // Src always points to the pre-image array.
        Point_vec* src = &pts1;
        Point_vec* dest = &pts2;

        std::for_each(begin, end, [&](auto element) {
            const Perm& perm{ *element };
            // Invalid iterator given if error happens here.

            acc_ = acc_ ^ perm.acc();

            for (size_t i = 0; i < size; ++i) {
                (*dest)[i] = (*src)[perm >> i];
            }

            std::swap(src, dest);
        });

        pre_images_ = std::move(*src);
        set_images();
    }

    //
    // Permutation expression operations.
    //

    /** Gets the pre-image of a point.
     */

    friend Point operator>>(const Perm& perm, Point point)
    {
        return perm.pre_images_[point];
    }

    /** Gets the image of a point.
     */

    friend Point operator<<(const Perm& perm, Point point)
    {
        return perm.images_[point];
    }

    /** Gets the accompanied action of a permutation.
     */

    A acc() const { return acc_; }

    /** Gets the size of the permutation domain.
     */

    size_t size() const { return images_.size(); }

    // Default assignment behaviour.

    Perm& operator=(const Perm& perm) = default;
    Perm& operator=(Perm&& perm) = default;

    /** Gets the earliest point moved by the permutation.
     *
     * A point `size` will be returned if the permutation is identity.
     */

    Point get_earliest_moved() const
    {
        Point i;
        for (i = 0; i < size(); ++i) {
            if (i != *this >> i)
                break
        }
        return i;
    }

private:
    /** Sets the image array from the pre-image array.
     */

    void set_images()
    {
        images_.resize(size());
        for (size_t i = 0; i < size(); ++i) {
            images_[pre_images_[i]] = i;
        }
    }

    Point_vec images_;
    Point_vec pre_images_;
    A acc_;
};

/** Expression for inverted permutations.
 *
 * This is a lazy expression for inverted permutations.  A constant reference
 * is stored for its operand.  So copy, move, and assignment are all
 * disallowed.
 */

template <typename T> class Inv_perm : public Perm_expr<Inv_perm<T>> {
public:
    Inv_perm(const T& operand)
        : operand_(operand)
    {
    }

    //
    // Permutation expression operations.
    //

    /** Gets the pre-image of a point.
     */

    friend Point operator>>(const Inv_perm& perm, Point point)
    {
        return perm.operand_ << point;
    }

    /** Gets the image of a point.
     */

    friend Point operator<<(const Inv_perm& perm, Point point)
    {
        return perm.operand_ >> point;
    }

    /** Gets the accompanied action of a permutation.
     */

    auto acc() const { return operand_.acc() | operand_.acc(); }

    /** Gets the size of the permutation domain
     */

    size_t size() const { return operand_.size(); }

private:
    const T& operand_;
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
    /** Constructs the expression from the operands.
     */

    Perm_prod(const T1& left, const T2& right)
        : left_(left)
        , right_(right)
    {
    }

    //
    // Permutation expression operations.
    //

    /** Gets the pre-image of a point.
     */

    friend Point operator>>(const Perm_prod& perm, Point point)
    {
        return perm.left_ >> (perm.right_ >> point);
    }

    /** Gets the image of a point.
     *
     * Despite our notation, the permutations still action from right, as in
     * the convention of computational group theory.
     */

    friend Point operator<<(const Perm_prod& perm, Point point)
    {
        return perm.right_ << (perm.left_ << point);
    }

    /** Gets the accompanied action of a permutation.
     */

    auto acc() const { return left_.acc() ^ right_.acc(); }

    /** Gets the size of the permutation domain
     *
     * No checking is done here.  Checking are performed during (before) the
     * construction.
     */

    size_t size() const { return left_.size(); }

private:
    const T1& left_;
    const T2& right_;
};

/** Multiplies two permutation expressions.
 *
 * The result is a lazily-evaluated permutation expression.  Note that the two
 * operands need to have the same size and accompanied action type.
 */

template <typename T1, typename T2>
auto operator|(const Perm_expr<T1>& left, const Perm_expr<T2>& right)
{
    static_assert(
        std::is_same<decltype(left.acc()), decltype(right.acc())>::value,
        "Two operands need to have the same accompanied action type.");
    assert(left.size() == right.size()
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
    { transv.next() } -> T*;
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

    /** Process a permutation for a transversal.
     *
     * This is the core function for transversal system adaptation.  If the
     * given permutation is in the subgroup of the transversal, it will be
     * added to the container for permutations to pass.  Or it will be added to
     * the given transversal container if no representative is already present
     * for its coset.
     */

    template <typename P, typename T, typename V>
    void proc_perm_for_transv(P&& perm, T& transv, V& perms_to_pass)
    {
        if (transv.has(perm)) {
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
 * group up to the same leaf.  Note that this function moves all permutations
 * out of the input transversal system.  So the input transversal will become
 * an empty one after the adaptation.
 *
 * The transversals can be any data type satisfying the concept of a
 * transversal system `Transv`.  Most importantly,
 *
 * - A type `Perm` needs to be defined for the permutation type it contains.
 *
 * - It needs to be an iterable of permutations in the transversal system.
 *   Note that identity *should* be skipped in the iteration, since it is
 *   always present.
 *
 * - It needs to have a method named `next` for a pointer to the next level of
 *   subgroup.  End of the chain should be given by a NULL value.
 *
 * - It needs to support methods `has`, `get_repr`, and `insert`, to test if a
 *   permutation is in the subgroup of the transversal system, to get a
 *   constant pointer to a coset representative (NULL if not present), and to
 *   insert a permutation into the transversal.
 */

template <typename T> void adapt_trasv(T& input, T& output)
{
    // Gather all transversals to loop over them in reverse order.
    std::vector<T*> inputs{};
    for (T* curr = &input; curr; curr = curr->next()) {
        inputs.push_back(curr);
    }

    std::for_each(inputs.rbegin(), input.rend(), [&](T* curr_input) {

        using Perm_vector = std::vector<typename T::Perm>;
        Perm_vector passed_perms{};
        std::move(begin(*curr_input), end(*curr_input),
            std::back_inserter(passed_perms));

        for (T* curr_output = &output; curr_output;
             curr_output = curr_output->next()) {

            std::vector<const typename T::Perm*> existing{};
            std::transform(begin(*curr_output), end(*curr_output),
                std::back_inserter(existing),
                [](const auto& perm) { return &perm; });

            // Passed permutation times identity. It is treated before the
            // actual products so that we can use move whenever it is
            // possible.
            //
            // TODO: fix the moving semantics here.

            for (const auto& passed_perm : passed_perms) {
                internal::proc_perm_for_transv(
                    std::move(passed_perm), *curr_output, perms_to_pass);
            }

            // Other products, passed times existing.

            for (const auto& passed_perm : passed_perms) {
                for (auto existing_perm : existing) {
                    internal::proc_perm_for_transv(passed_perm | *existing_perm,
                        *curr_output, perms_to_pass);
                }
            }

            passed_perms.clear();
            passed_perms.swap(perms_to_pass);

        } // End loop output level.

    }); // End loop input level.

    // If the two transversal are indeed for the same core kernel, we should
    // have nothing left now.
    assert(passed_perms.empty());
    return;
}

} // End namespace libcanon.

#endif // LIBCANON_PERM_H
