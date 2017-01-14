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

#include <algorithm>
#include <cassert>
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

template <typename R> using Coset_of = typename R::Coset;

/** The structure type for a refiner.
 *
 * The type is obtained by just read from the `Structure` attribute from the
 * refiner class.
 */

template <typename R> using Structure_of = typename R::Structure;

/** Permutation type for a refiner.
 *
 * This meta function also requires the presence of a `get_a_perm` method to
 * return a permutation of a leaf coset when it is called with the coset.
 */

template <typename R>
using Perm_of = std::result_of_t<decltype (&R::get_a_perm)(Coset_of<R>)>;

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

/** The concept for a refiner.
 *
 * Here basically the kinds of operations needed for a refiner are all listed,
 * mostly for documentation purpose.  So all expressions are added even if they
 * are already checked in the type meta functions for refiners.
 */

template<typename R>
concept bool Refiner = requires () {
    typename Coset_of<R>;
    typename Structure_of<R>;
    typename Perm_of<R>;
    typename Action_res_of<R>;
    typename Transv_of<R>;
} && requires (R refiner,
        Coset_of<R> coset, Structure_of<R> obj,
        Perm_of_R<R> perm, Perm_of_R<R> perm2,
        Transv_of_R<R> transv, Transv_of_R<R> target_transv) {
    // For refiner itself.
    { refiner.refine(obj, coset) } -> Simple_iterable<Coset_of<R>>;
    { refiner.is_leaf(obj, coset) } -> bool;

    // For the action.
    { refiner.act(perm, obj) } -> Act_res_of<R>;
    { refiner.left_mult(perm, coset) } -> Coset_of<R>;

    // For the transversal container.
    { transv.insert(perm | ~perm2) };
    { adapt_transv(transv, target_transv) }
};

/** Concept for a container able to work with a refiner.
 *
 * Here we require that the container is able to hold pairs of an action result
 * and a permutation.  Also the container can be queried for the presence of
 * any form, to yield a pair composed from an action result and a permutation.
 */

template<typename R, typename C>
concept bool Refiner_container = requires (
        Container container, Act_res_of<R> res, Perm_of<R> perm) {
    { container.emplace(res, perm) };
    { *container.find(res) } -> std::pair<const Act_res_of<R>, Perm_of<R>>;
    typename Container::reference;
};

// clang-format on
#endif

//
// Some internal data structure and algorithms.
//

/** Experimental path from a given coset.
 *
 * The central class for the generic canonicalization algorithm.  All core work
 * are performed here.
 *
 * When given leaf cosets, here a permutation and the corresponding form is
 * stored.  For non-leaf cosets, here we store all its children, one of which
 * is guaranteed to have a path extended from it.
 *
 * A candidate and its permutation from the successive refinement of the given
 * coset can be queried.  And all candidates descending from it can be added if
 * it is needed, with or without a given automorphism group.
 *
 * Note that this class is neither copyable or assignable.
 */

template <typename R> class Exp_path {
public:
    //
    // Data types for the refiner.
    //

    using Coset = Coset_of<R>;
    using Structure = Structure_of<R>;
    using Transv = Transv_of<R>;
    using Act_res = Act_res_of<R>;
    using Perm = Perm_of<R>;

    /**
     * Initialize experimental path object.
     */

    Exp_path(R& refiner, const Structure& obj, const Coset& coset)
        : refiner_(refiner)
        , obj_(obj)
        , base_(&coset)
    {
        if (refiner.is_leaf(obj, coset)) {
            // For leaf states, we store the permutation and the action result.

            perm_ = std::make_unique<Perm>(refiner_.get_a_perm(coset));
            form_ = std::make_unique<Act_res>(refiner_.act(*perm, obj));
        } else {
            // For non-leaf states, we store all the refinements and extend a
            // path from one of the children.

            auto children = refiner.refine(coset);
            auto child_iter = begin(children);
            auto child_sentinel = end(children);

            // Non-leaf cosets are guaranteed to have non-empty refinement.
            assert(child_iter != child_sentinel);

            std::for_each(child_iter, child_sentinel, [&](auto&& child) {
                children_.emplace(
                    std::forward<decltype(child)>(child), nullptr);
            });

            set_curr();
        }
    }

    //
    // Disabling of some unnecessary operations.
    //

    Exp_path(const Exp_path& exp_path) = delete;
    Exp_path(Exp_path&& exp_path) = delete;
    Exp_path& operator=(const Exp_path& exp_path) = delete;
    Exp_path& operator=(Exp_path&& exp_path) = delete;

    /** Gets the form in a leaf node.
     */

    Act_res& form() const
    {
        assert(is_leaf());
        return *form_;
    }

    /** Gets the permutation in a leaf node.
     */

    Perm& perm() const
    {
        assert(is_leaf());
        return *perm_;
    }

    /** Gets the leaf experimental path node.
     */

    const Exp_path& get_a_leaf() const
    {
        return is_leaf() ? *this : curr_exp_path_->get_a_leaf();
    }

    /** Prepares a transversal system best suited for the current path.
     *
     * This transversal system has the current nodes at each level as the
     * anchor points.
     */

    std::unique_ptr<Transv> prepare_transv() const
    {
        if (is_leaf()) {
            return nullptr;
        } else {
            // For non-leaf nodes.
            auto transv = refiner.create_transv(base_, *curr_coset_);
            transv.set_next(curr_exp_path_->prepare_transv());
            return transv;
        }
    }

    /** Adds all candidates from the given coset with automorphism group.
     *
     * All candidates from the base coset will be added to the given container.
     *
     * \param aut A pointer to a transversal system giving the automorphism
     * group of the object in the conjugated subgroup of the base.  It is
     * assumed that the transversal system is already adapted to the current
     * subgroup chain.
     *
     * \return A pointer to a new transversal system for the same subgroup of
     * the automorphism group.
     */

    template <typename C>
    std::unique_ptr<Transv> add_all_candidates(
        std::unique_ptr<Transv> aut, C& container)
    {
        if (is_leaf()) {
            // Both the form and the permutation can be safely moved since
            // they are not going to be used after this method is called.
            container.emplace(std::move(*form_), std::move(*perm_));
            return nullptr;
        }

        while (curr_coset_) {
            // Add all candidates from the current children.
            //
            // Here care is taken for the undefined evaluation order.
            auto new_inner = curr_exp_path_->add_all_candidates(
                aut.release_next(), container);
            aut.set_next(std::move(new_inner));

            // Remove all identical siblings of the current children.
            const auto& anchor = *curr_coset_;
            for (const auto& i : aut) {
                auto n_erased = children_.erase(refiner_.left_mult(i, anchor));
                assert(n_erased == 1);
            }
            auto n_erased = children_.erase(anchor);
            assert(n_erased == 1);

            if (!children.empty()) {
                set_curr();
                auto new_aut = prepare_transv();
                adapt_transv(*aut, *new_aut);
                aut = std::move(new_aut);
            } else {
                curr_coset_ = nullptr;
                curr_exp_path_ = nullptr;
            }
        }

        return aut;
    }

    /** Adds all candidates from without known automorphism.
     *
     * All candidates from the coset `base` will be added to the given
     * container.
     *
     * A pointer to a new transversal system from the conjugate group for the
     * base coset is returned.
     */

    template <typename C>
    std::unique_ptr<Transv> add_all_candidates(C& container)
    {
        std::unique_ptr<Transv> aut{}; // Named return value.

        if (is_leaf()) {
            container.emplace(std::move(*form_), std::move(*perm_));
            return aut;
        }

        aut = create_transv(base_, *curr_coset_);
        aut->set_next(curr_exp_path_->add_all_candidates(container));

        // Create experimental paths for all refinement to find all the
        // identical siblings.
        set_paths();

        // Loop over all the experimental path by iterator for unordered_map
        // for faster removal without hashing and equality comparison involved.

        auto child_it = children_.begin();
        while (child_it != children.end()) {
            Exp_path* curr = child_it->second.get();

            const auto& leaf = curr->get_a_leaf();

            auto existing = container.find(leaf.form());
            if (existing == container.end()) {
                // Non identical siblings will be treated later.
                ++i;
                continue;
            } else {
                // When we find an identical sibling.
                if (curr != curr_exp_path_) {
                    // Skip identity permutation.
                    aut->insert(existing->second | ~leaf.perm());
                }
                auto to_remove = i; // Make a copy of the iterator.
                ++i;
                children_.erase(to_remove);
                // Note that the pointers curr_exp_path_ and curr_coset_
                // will be invalidated after this loop.
            }
        }

        // Now we have treated an identical class of siblings and got a
        // complete description of the automorphism group.

        if (set_curr()) {
            auto new_aut = create_transv(base_, *curr_coset_);
            adapt_transv(*aut, *new_aut);
            auto final_aut
                = this->add_all_candidates(std::move(new_aut), container);
            return final_aut;
        } else {
            return aut;
        }
    }

private:
    //
    // Internal utility methods.
    //

    /** Sets the current coset and experimental path to begin of children.
     *
     * False will be returned if the children container is empty.
     */

    bool set_curr()
    {
        auto begin = children_.begin();
        if (begin == children_.end())
            return false;

        curr_coset_ = &begin->first;
        set_path(*begin);
        curr_exp_path_ = begin->second.get();

        return true;
    }

    /** Sets the paths for all children.
     */

    void set_paths()
    {
        for (auto i& : children_) {
            set_path(i);
        }
    }

    /** Sets the path of a given child.
     *
     * The path is created only when there is no one already.
     */

    void set_path(Children::reference child)
    {
        if (!child.second) {
            child.second
                = std::make_unique<Exp_path>(refiner_, obj_, child.first);
        }
    }

    /** Tests if an experimental path node is a leaf state.
     */

    bool is_leaf() const { return perm_; }

    //
    // Data fields.
    //

    // References to basic information for convenience.
    R& refiner_;
    const Structure& obj_;
    const Coset& base_;

    // For non-leaf nodes
    using Children = std::unorderd_map<Coset, std::unique_ptr<Exp_path>>;
    Children children_;
    Coset* curr_coset_;
    Exp_path* curr_exp_path_;

    // For leaf nodes.
    std::unique_ptr<Perm> perm_;
    std::unique_ptr<Structure> form_;
};

/** Adds all candidates from the successive refinement of the given coset.
 *
 * This is a shallow wrapper over the actual methods in the Exp_path class.  In
 * this way, a pure functional interface for the facility is available without
 * having to touch the data structure.
 */

template <typename R, typename C>
auto add_all_candidates(R& refiner, const Structure_of<R>& obj,
    const Coset_of<R>& coset, C& container)
{
    Exp_path<R> exp_path(refiner, obj, coset);
    return exp_path.add_all_candidates(container);
}

} // End namespace libcanon.

#endif // LIBCANON_CANON_H
