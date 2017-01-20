/**
 * \file
 *
 * Small utilities to be reused.
 */

#ifndef LIBCANON_UTILS_H
#define LIBCANON_UTILS_H

namespace libcanon {

#ifdef __cpp_concepts
// clang-format off

//
// Simple iterable concept.
//
// This is certainly going to be removed when the range standard library goes
// forward.
//

template<typename I, typename T>
concept bool Simple_iterable = requires (I iterable) {
    { begin(i) };
    { end(i) };
    { ++begin(i) };
    { begin(i) != end(i) } -> bool;
    { *begin(i) } -> T&&;
}

// clang-format on
#endif

/** Combines the given hash values.
 *
 * The algorithm is adapted from the boost hash library.
 *
 * Normally the data type for the hash can just be `size_t`.
 */

template <typename T> void combine_hash(T& seed, T value)
{
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/** Ensures the type to be a unique pointer.
 *
 * Similar to the standard `remove_reference` meta function, the type
 * inside the unique pointer will be defined on the `type` attribute.
 */

template <typename T> struct Ensure_unique_ptr {
};

template <typename T> struct Ensure_unique_ptr<std::unique_ptr<T>> {
    using type = T;
};

template <typename T>
using Ensure_unique_ptr_t = typename Ensure_unique_ptr<T>::type;

/** Compares by degrev ordering wrapping over the original ordering.
 *
 * The two objects needs to be from classes with `size` method.  The one with
 * greater size is considered less, as in the degrevlex ordering commonly used
 * in Groebner basis theory.  When the sizes are equal, the two quantities will
 * be cast to the given base type and compared.
 *
 * Note that the given base type need to be given as the base type itself
 * without ref or cv quantification.
 */

template <typename Base, typename T>
bool is_degrev_less(const T& left, const T& right)
{
    size_t size_l = left.size();
    size_t size_r = right.size();

    using Base_const_ref = const Base&;

    return size_l > size_r || (size_l == size_r
                                  && static_cast<Base_const_ref>(left)
                                      < static_cast<Base_const_ref>(right));
}

} // End namespace libcanon.

#endif // LIBCANON_UTILS_H
