/** \file
 *
 * General canonicalization of combinatorial objects.
 *
 * Synopsis
 * --------
 *
 * This file contains generic driver functions for the canonicalization of
 * combinatorial objects.  The problem specific information should be mostly
 * put into the refiner class, whose type must implement the refiner protocol,
 * which is defined as follows.
 *
 *
 * Refiner protocol
 * ----------------
 *
 * - The class for refiners needs to define the following types,
 *
 *   + `Coset` for the type of a coset in the isomorphism group.
 *
 *   + `Structure` for the type of the combinatorial structure.
 *
 * - It has a method named `refine` capable of returning an iterable of
 *   `Coset`s when it is given one combinatorial object and a coset.
 *
 * - It has a method named `is_leaf`, to be called with an combinatorial object
 *   and a coset, to decide if a given coset is a leaf.
 *
 * - It has a method named `get_a_perm` to get a permutation from a coset.
 *
 * - The permutations need to be able to act on objects of the combinatorial
 *   structure and cosets by methods `act(perm, obj)` or `left_mult(perm,
 *   coset)` of the refiner.
 *
 * - The refiner need to have a method named `create_transv` capable of taking
 *   a coset $aH$ and another coset $bK$ to return a unique pointer to a
 *   transversal container for $b K b^{-1}$ in $a H a^{-1}$.  The transversal
 *   must satisfy the concept for a transversal container defined in `Perm.h`.
 *   And it needs to accept permutations given by  `perm1 | ~perm2` by method
 *   `insert`.
 *
 * All of these requirements are put in the `Refiner` concept, which is
 * automatically used on compilers supporting the C++ Concept-Lite technical
 * report.
 *
 */

#ifndef LIBCANON_CANON_H
#define LIBCANON_CANON_H

#include <memory>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <libcanon/utils.h>

namespace libcanon {

//
// Meta functions for types in a Refiner
//
// All these meta functions are named with ending `_of` to deduce the
// corresponding type of a given refiner.  Failure of reduction indicates
// dissatisfaction of the refiner concept.
//

/** The coset type for a refiner.
 *
 * This meta function ensures that a type `Coset` is defined inside it.
 */

template <typename R> using Coset_of<R> = typename R::Coset;

/** The structure type for a refiner.
 *
 * The type is obtained by just read from the `Structure` attribute from the
 * refiner class.
 */

template <typename R> using Structure_of<R> = typename R::Structure;

/** Permutation type for a refiner.
 *
 * This meta function also requires the presence of a `get_a_perm` method to
 * return a permutation of a leaf coset when it is called with the coset.
 */

template <typename R>
using Perm_of<R> = std::result_of_t<decltype (&R::get_a_perm)(Coset_of<R>)>;

/** The result of acting a permutation on a structure.
 *
 * Sometimes the result of acting a permutation on a structure is an object of
 * the same type by itself.  However, there are cases where this flexibility is
 * very useful.  This meta function also ensures the presence of the `act`
 * method.
 */

template <typename R>
using Act_res_of
    = std::result_of_t<decltype (&R::act)(Perm_of<R>, Structure_of<R>)>;

/** The transversal type used by a refiner.
 *
 * This meta function also ensures that the `create_transv` function is present
 * and returns a unique pointer to something.
 */

// clang-format off

template <typename R>
using Transv_of = internal::Ensure_unique_ptr_t<std::result_of_t<
    decltype(&R::create_transv)(Coset_of<R>, Coset_of<R>)
>>;

// clang-format on

//
// The refiner protocol.
//

#ifdef __cpp_concepts
// clang-format off

//
// The refiner concept.
// 

template<typename R>
concept bool Refiner = requires () {
    typename R::Coset;
    typename R::Perm;
    typename R::Transv;
    typename R::Structure;
    typename R::Container;
} && requires (R refiner, typename R::Coset coset, typename R::Structure obj, 
        std::unique_ptr<typename R::Transv> transv, 
        typename R::Transv* target_transv, typename R::Container container, 
        typename R::Perm perm, typename R::Perm perm2) {
    // For refiner itself.
    { refiner.refine(coset) } -> Simple_iterable<typename R::Coset>;
    { refiner.is_leaf(coset) } -> bool;
    { refiner.get_a_perm(coset) } -> typename R::Perm;

    // For the action.
    { obj >> perm } -> R::Structure;
    { coset >> perm } -> R::Coset;

    // For the container.
    { container.emplace(obj, perm) };
    { container.emplace(obj, perm | ~perm2) };

    // For the transversal container.
    {
        refiner.create_transv(coset, transv) 
    } -> std::unique_ptr<typename R::Transv>;
    { adapt_transv(transv, target_transv) }
};

// clang-format on
#endif

//
// Some internal data structure and algorithms.
//

namespace internal {

    /** Applies given action for each children of the given coset.
     *
     * Basically the same semantics as the `for_each` algorithm in the standard
     * library, just here the callback is called with an extra boolean argument
     * giving if the children is the first one.
     *
     * Note that the given coset is assumed to be non-leaf.  And the result is
     * assumed to be non-empty.
     */

    template <typename R, typename H>
    void for_each_children(
        R& refiner, const typename R::Coset& coset, H handler)
    {
        auto children = refiner.refine(coset);
        auto child_iter = begin(children);
        auto child_sentinel = end(children);

        handler(*child_iter, true);
        for (++child_iter; child_iter != child_sentinel; ++child_iter) {
            handler(*child_iter, false);
        }
    }

    /** Experimental path from a given coset.
     *
     * A candidate and its permutation from the successive refinement of the
     * given coset can be queried.  And all candidates can be added if it is
     * needed, given an automorphism group.
     */

    template <typename R> class Exp_path {
    public:
        /**
         * Initialize experimental path object.
         */

        Exp_path(R& refiner, const typename R::Structure& obj,
            const typename R::Coset& coset)
            : refiner{ refiner }
            , obj{ obj }
            , base{ &coset }
        {
            if (refiner.is_leaf(coset)) {
                perm = std::make_unique<R::Perm>(refiner.get_a_perm(coset));
                form = std::make_unique<R::Structure>(obj >> *perm);
            } else {
                internal::for_each_children(
                    refiner, coset, [this](auto&& children, bool if_first) {
                        this->children.insert(
                            std::forward<decltype(children)>(children));
                    });
                curr = std::make_unique<Exp_path>(
                    refiner, obj, *(children.begin()));
            }
        }

        /**
         * Sets the base point of the experimental path node.
         */

        void set_base(const typename R::Coset& coset) { base = &coset; }

        /**
         * Gets the form in a leaf node.
         */
        const auto& get_form() { return *form; }

        /**
         * Gets the permutation in a leaf node.
         */

        const auto& get_perm() { return *perm; }

        /**
         * Gets the leaf experimental path node.
         */

        const Exp_path& get_a_leaf() const
        {
            return curr ? curr->get_a_leaf() : *this;
        }

        std::unique_ptr<typename R::Transv> prepare_transv() const
        {
            return refiner.create_transv(
                *base, curr ? curr->prepare_transv() : nullptr);
        }

        /** Adds all candidates from the given coset.
         *
         * All candidates from the coset `base` will be added to the given
         * container.
         *
         * \param aut A pointer to a transversal system giving the automorphism
         * group of the object in the conjugated subgroup of the base.  It is
         * assumed that the transversal system is already adapted to the
         * current subgroup chain.
         *
         * \return A pointer to a new transversal system for the same subgroup
         * of the automorphism group.
         */

        std::unique_ptr<typename R::Transv> add_all_candidates(
            std::unique_ptr<typename R::Transv> aut,
            typename R::Container& container)
        {
            if (form) {
                container.emplace(std::move(*form), std::move(*perm));
            }

            while (curr) {
                aut->next
                    = curr->add_all_candidates(std::move(aut.next), container);

                const auto& anchor = *(curr->base);
                for (const auto& i : aut) {
                    children.erase(anchor >> i);
                }
                children.erase(anchor);

                if (!children.empty()) {
                    curr = std::make_unique<Exp_path>(
                        refiner, obj, children.begin());
                    auto new_aut = curr.prepare_transv();
                    adapt_transv(std::move(aut), new_aut.get());
                    aut = std::move(new_aut);
                } else {
                    curr = nullptr;
                }
            }

            return std::move(aut);
        }

    private:
        R& refiner;
        const typename R::Structure& obj;
        const typename R::Coset* base;

        std::unordered_set<typename R::Coset> children;
        std::unique_ptr<Exp_path> curr;

        std::unique_ptr<typename R::Perm> perm;
        std::unique_ptr<typename R::Structure> form;
    };

} // End namespace internal

/**
 * Adds all candidates from the successive refinement of the given coset.
 *
 * This is the core function of the canonicalization module.  All candidates
 * from the successive refinement of the given coset is going to be added to the
 * given container.
 */

template <typename R>
auto add_all_candidates(R& refiner, const typename R::Structure& obj,
    const typename R::Coset& coset, typename R::Container& container)
{
    std::unique_ptr<typename R::Transv> aut{}; // Named return value.

    if (refiner.is_leaf(coset)) {
        auto perm = refiner.get_a_perm(coset);
        container.emplace(obj >> perm, std::move(perm));
        aut = refiner.create_transv(coset, nullptr);
    }

    using Exp_path_type = internal::Exp_path<R>;
    using Exp_path_container
        = std::unordered_map<typename R::Coset, std::unique_ptr<Exp_path_type>>;
    Exp_path_container pending{};

    internal::for_each_children(refiner, coset, [&](auto child, bool if_first) {
        if (if_first) {
            std::unique_ptr<typename R::Transv> inner_aut;
            inner_aut = add_all_candidates(refiner, obj, child, container);
            aut = refiner.create_transv(coset, std::move(inner_aut));
        } else {
            auto exp_path
                = std::make_unique<Exp_path_type>(refiner, obj, child);
            const auto& leaf = exp_path.get_a_leaf();

            auto prev_perm = container.find(leaf.get_form());
            if (prev_perm == container.end()) {
                auto res
                    = pending.emplace(std::move(child), std::move(exp_path));
                // Reset the pointer to base.
                auto& kv_pair = *res.first;
                kv_pair.second->set_base(kv_pair.first);
            } else {
                aut.emplace(prev_perm->second | ~leaf.ger_perm());
            }
        }
    });

    while (!pending.empty()) {
        auto& trial = *pending.begin();
        auto& child = trial.first;
        auto& exp_path = *trial.second; 
        auto new_aut = exp_path.prepare_transv();
        adapt_transv(std::move(aut), new_aut.get());
        aut = std::move(new_aut);
        aut.next = exp_path.add_all_candidates(std::move(aut.next), container);
        for (const auto& i : aut) {
            pending.erase(child >> i);
        }
        pending.erase(child);
    }

    return aut;
}

} // End namespace libcanon.

#endif // LIBCANON_CANON_H
