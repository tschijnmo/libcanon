/**
 * General canonicalization of combinatorial objects.
 */

#ifndef LIBCANON_CANON_H
#define LIBCANON_CANON_H

#include<tuple>
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


/**
 * Canonicalize the combinatorial object.
 *
 * The given combinatorial object will be canonicalized by successive refinement
 * of its isomorphism group.
 *
 */

template<typename O, typename G, typename C>
auto canonicalize(const O& obj, const G& isom_grp, C& caner)
requires Canonicalizable<O, C> && Tractable<G, C>
{
    // Initialize the container for the automorphism group.
    auto full_group = caner.get_full_coset(isom_grp);

    // Create the lazy generator of the candidates.
    internal::Candidates<Obj, Refiner> candidates{obj, full_group, caner};

    // Choose the canonical isomorph and get the automorphism group.
    return std::make_tuple(
        caner.choose(candidates.cbegin(), candidates.cend()),
        std::move(candidates.get_aut())
    );
};

}  // End namespace libcanon.

#endif //LIBCANON_CANON_H
