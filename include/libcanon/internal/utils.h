/**
 * Small utilities.
 */

#ifndef LIBCANON_UTILS_H
#define LIBCANON_UTILS_H

#include <type_traits>
#include <utility>

namespace libcanon::internal {

//
// Concepts related to iterator and range.
//
// They are sure going to be removed when the C++ standard moves forward.

template<typename I, typename E, typename V>
concept bool Simple_iterator = requires (I iter, E end) {
    { ++iter } -> &I;
    { *iter } -> V;
    { iter != end } -> bool;
};

template<typename R, typename V>
concept bool Range_of = requires (R range) {
    { begin(range) };
    { end(range) };
    requires Simple_iterator<
        std::remove_reference_t<decltype(begin(range))>,
        std::remove_reference_t<decltype(end(range))>, V
    >;
};

}
#endif //LIBCANON_UTILS_H
