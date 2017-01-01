/** Canonicalization of strings.
 *
 * In this module, the generic canonicalization algorithm is applied to string
 * canonicalization problem.  It is primarily comprised of two parts, the
 * compilation of any list of generators for a permutation group to a
 * transversal system for pointwise stabilizers, and a refiner based on this
 * transversal system.
 */

#ifndef LIBCANON_STRING_CANON_H
#define LIBCANON_STRING_CANON_H

#include <memory>
#include <vector>

#include <libcanon/perm.h>

namespace libcanon {

/** Transversal systems for pointwise stabilisers.
 *
 * This class is designed to hold transversal systems for subgroup chains of
 * pointwise stabilizers.  Here each level of subgroup will have one point
 * stabilized, called the target, and all the cosets are labelled by the point
 * that is moved into the target.
 */

template <typename P> class SimsTransv {
public:
    //
    // Required by transversal interface
    //

    /** The type of permutations contained. */
    using Perm = P;

    /** Unique pointer to the next level of transversal */
    std::unique_ptr<SimsTransv> next;

    /** Tests if the permutation could be in this subgroup. */
    bool has(const P& perm) { return target == (perm << target); }

    /** Gets the coset representative of a permutation. */
    const P* get_repr(const P& perm)
    {
        const auto& found = transv[perm >> target];

        // Slight redundancy for clarity.
        if (found) {
            return found.get();
        } else {
            return nullptr;
        }
    }

    /** Inserts a permutation into the transversal system
     *
     * No checking of previous coset representative is done.  The given
     * permutation will just overwrite the previous one.  Permutations inside
     * the subgroup are ignored.
     */

    template <typename T> void insert(T&& perm)
    {
        Point label = perm >> target;
        if (label == target)
            return;
        transv[label] = std::make_unique<P>(std::forward<T>(perm));
    }

    //
    // Problem specific
    //

    /**
     * Creates a transversal system.
     */

    SimsTransv(size_t target, size_t size)
        : target{ target }
        , transv{ size }
    {
    }

private:
    const Point target;
    std::vector<std::unique_ptr<P>> transv;
};
}

#endif // LIBCANON_STRING_CANON_H
