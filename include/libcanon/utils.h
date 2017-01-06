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

} // End namespace libcanon.

#endif // LIBCANON_UTILS_H
