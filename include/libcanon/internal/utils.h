/**
 * Small utilities.
 */

#ifndef LIBCANON_UTILS_H
#define LIBCANON_UTILS_H

#include <type_traits>
#include <tuple>
#include <utility>

#include<boost/iterator/transform_iterator.hpp>

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


//
// Utility functions.
//

/**
 * Get the iterator for keys from the given container.
 *
 * @param container A map-style container.
 * @param method The container method to get the original value iterator.
 * @return An iterator looping over the key only.
 */

template<typename T, typename R>
auto get_key_iter(const T& container, R (T::*method)())
{
    auto iter = (container.*method)();
    return boost::make_transform_iterator(iter, [](auto&& i) {
        return std::get<0>(std::forward<decltype(i)>(i));
    });
}


/**
 * Get the constant begin iterator for the keys.
 */

template <typename T>
auto get_key_cbegin(const T& container)
{
    return get_key_iter(container, T::cbegin);
}


/**
 * Get the constant end iterator for the keys.
 */
template <typename T>
auto get_key_cend(const T& container)
{
    return get_key_iter(container, T::cend);
}

}
#endif //LIBCANON_UTILS_H
