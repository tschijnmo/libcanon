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
// Concepts related to iterator and iterable.
//
// They are mainly mimicked and simplified version of ReadableIterator and Range
// tailored to the special needs of libcanon.  Sure they are going to be removed
// when the C++ standard moves forward.

template<typename I, typename S, typename V>
concept bool Simple_iterator = requires (I iter, S sentinel) {
    { ++iter } -> I&;
    { *iter } -> const V&;
    { iter != sentinel } -> bool;
};

template<typename T>
using iterator_t = decltype(begin(std::declval<T&>()));
template<typename T>
using sentinel_t = decltype(end(std::declval<T&>()));

template<typename I, typename V>
concept bool Simple_iterable = requires {
    typename iterator_t<I>;
    typename sentinel_t<I>;
    requires Simple_iterator<iterator_t<I>, sentinel_t<I>, V>;
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
auto get_key_iter(const T& container, R (T::*method)() const)
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
    return get_key_iter(container, &T::cbegin);
}


/**
 * Get the constant end iterator for the keys.
 */
template <typename T>
auto get_key_cend(const T& container)
{
    return get_key_iter(container, &T::cend);
}

}
#endif //LIBCANON_UTILS_H
