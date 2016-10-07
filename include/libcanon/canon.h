/**
 * General canonicalization of combinatorial objects.
 */

#ifndef LIBCANON_CANON_H
#define LIBCANON_CANON_H

#include<unordered_map>

#include<libcanon/internal/utils.h>

namespace libcanon {

//
// Concepts of the generic types
// -----------------------------
//

/**
 * The concept of a type of object being permutable by some other object.
 */
template<typename O, typename A>
concept bool Permutable = requires (O obj, A action) {
    {obj ^ action} -> O;
};


namespace internal {
template<typename O, typename C>
using Candidates_container = std::unordered_map<O, const C&>;
}


/**
 * The concept of a type of object being canonicalizable by a canonicalizer.
 *
 * The canonicalizer should hold all the information that is needed to
 * canonicalize the objects of the given type.
 */

template<typename O, typename C>
concept bool Canonicalizable = requires {
    typename C::Group_type;
    typename C::Coset_type;
    typename C::Perm_type;

    requires Permutable<O, C::Perm_type>;
} && requires (
        O obj, C caner, typename C::Coset_type coset,
        internal::Candidates_container<O, C::Coset_type> cont
) {
    { caner.refine(obj, coset) } -> internal::Range_of<C::Coset_type>;

    { caner.choose(
        internal::get_key_cbegin(cont), internal::get_key_cend(cont)
    ) } -> O;
};


/**
 * The concept that a symmetry specification is tractable by a canonicalizer.
 */
template<typename G, typename C>
concept bool Tractable = requires (G group, C caner) {
    { caner.get_full_coset(group) } -> typename C::Coset_type;
};

}  // End namespace libcanon.

#endif //LIBCANON_CANON_H
