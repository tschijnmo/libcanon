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
#include <iterator>
#include <memory>
#include <vector>

#include <libcanon/utils.h>

namespace libcanon {

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
    { transv.next } -> std::unique_ptr<T>;
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

    /**
     * Tests if a transversal is a leaf.
     *
     * A leaf transversal is one without a next transversal.
     */

    template <typename T> bool is_leaf_transv(const T& transv)
    {
        return transv.next == nullptr;
    }

    /** Process a permutation for a transversal.
     *
     * This is the core function for transversal system adaptation.  If the
     * given permutation is in the subgroup next level, it will be added to the
     * container for permutations to pass.  Or it will be added to the given
     * transversal container if no representative is already present for its
     * coset.
     */

    template <typename P, typename T, typename V>
    void proc_perm_for_transv(P&& perm, T& transv, V& perms_to_pass)
    {
        if (transv.next->has(perm)) {
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
 * group up to the same leaf.  Note that this function takes the ownership of
 * the input transversal and deallocates it after the adaptation.
 *
 * The transversals can be any data type satisfying the concept of a
 * transversal system `Transv`.  Most importantly,
 *
 * - A type `Perm` needs to be defined for the permutation type it contains.
 *
 * - It needs to be an iterable of permutations in the transversal system.
 *   Note that identity should be skipped in the iteration, since it is always
 *   resent.
 *
 * - It needs to have an attribute named `next` for a unique pointer to the
 *   next level of subgroup.
 *
 * - It needs to support methods `has`, `get_repr`, and `insert`, to test if a
 *   permutation is in the subgroup of the transversal system, to get a
 *   constant pointer to a coset representation (NULL if not present), and to
 *   insert a permutation into the transversal.
 *
 */

template <typename T> void adapt_trasv(std::unique_ptr<T> input, T& output)
{
    // Gather all transversals to loop over them in reverse order.
    std::vector<T*> inputs{};
    T* curr = input.get();
    while (!internal::is_leaf_transv(*curr)) {
        inputs.push_back(curr);
        curr = curr->next.get();
    }

    for (auto curr_input = inputs.rbegin(); curr_input != input.rend();
         ++curr_input) {

        using Perm_vector = std::vector<typename T::Perm>;
        Perm_vector passed_perms{};
        for (auto& i : *curr_input) {
            passed_perms.push_back(std::move(i));
        }
        Perm_vector perms_to_pass{};

        T* curr_output = &output;
        while (!internal::is_leaf_transv(*curr_output)) {

            std::vector<const typename T::Perm*> existing{};
            std::transform(begin(*curr_output), end(*curr_output),
                std::back_inserter(existing), [](auto& perm) { return &perm; });

            for (const auto& passed_perm : passed_perms) {
                for (auto existing_perm : existing) {
                    internal::proc_perm_for_transv(passed_perm | *existing_perm,
                        *curr_output, perms_to_pass);
                }
                internal::proc_perm_for_transv(
                    std::move(passed_perm), *curr_output, perms_to_pass);
            }

            curr_output = curr_output->next.get();
            passed_perms.clear();
            passed_perms.swap(perms_to_pass);

        } // End loop output level.
    } // End loop input level.

    return;
}

} // End namespace libcanon.

#endif // LIBCANON_PERM_H
