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
}
#endif // LIBCANON_UTILS_H
