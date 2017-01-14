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

} // End namespace libcanon.

#endif // LIBCANON_UTILS_H
