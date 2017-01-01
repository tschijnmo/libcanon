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
     * The given permutation is inserted only when there is no representative
     * for that coset yet.
     */

    template <typename T> void insert(T&& perm)
    {
        Point label = perm >> target;
        if (label == target || transv[label])
            return;
        transv[label] = std::make_unique<P>(std::forward<T>(perm));
    }

    //
    // Problem specific
    //

    /** Inserts a permutation given by a unique pointer.
     *
     * Different from general insertion of a permutation, here no checking of
     * previous coset representative is done.  The given permutation will just
     * overwrite the previous one.  Permutations inside the subgroup are
     * ignored.
     */

    void insert(std::unique_ptr<P> perm)
    {
        Point label = *perm >> target;
        if (label == target)
            return;
        transv[label] = std::move(perm);
    }

    /**
     * Gets a representative for the coset of moving the given point to target.
     */

    const P* get_repr(Point point)
    {
        return transv[point].get();
    }

    /**
     * Creates a transversal system.
     */

    SimsTransv(size_t target, size_t size)
        : target{ target }
        , transv{ size }
    {
    }

    Point get_target() const { return target; }

private:
    const Point target;
    std::vector<std::unique_ptr<P>> transv;
};

//
// Schreier-Sims algorithm for forming a transversal system
//

namespace internal {

    /** Builds one level of Sims transversal system.
     *
     * The earliest non-stabilized point will be found and its orbit is used to
     * build a transversal system of the pointwise stabilizer subgroup of that
     * point.  Null will be returned when all points up to the given size are
     * stabilized.
     */

    template <typename P>
    std::unique_ptr<SimsTransv<P>> build_sims_transv(
        Point begin, size_t size, std::vector<P>& gens)
    {

        // First we find the orbit and build a (right) transversal system by
        // the standard algorithm.

        std::vector<Point> orbit{};
        orbit.reserve(size);
        std::unique_ptr<SimsTransv<P>> tentative;

        for (size_t target = begin; target < size; ++target) {
            orbit.clear();
            orbit.push_back(target);

            tentative = std::make_unique<SimsTransv<P>>(target, size);

            size_t curr_idx = 0;
            while (curr_idx < orbit.size()) {
                auto curr_point = orbit[curr_idx];
                auto curr_perm = tentative.get_repr(curr_point);

                for (const auto& i : gens) {
                    auto loc = i << curr_point;
                    auto repr = tentative.get_repr(loc);
                    if (loc == target || repr)
                        continue;
                    // Now we have a new point in orbit.
                    orbit.push_back(loc);

                    if (curr_point == target) {
                        tentative.insert(~i);
                    } else {
                        tentative.insert(~i | *repr);
                    }
                }
            } // End while looping over orbit.

            if (orbit.size() == 1) {
                tentative.release();
            } else {
                // We found a point not stabilized.
                break;
            }
        } // End loop over point to find a non-stabilized.

        if (!tentative)
            return nullptr;

        std::vector<P> schreier_gens
            = form_schreier_gens(*tentative.get(), gens);
        gens.swap(schreier_gens);

        return std::move(tentative);
    } // End function build_sims_transv
}

/** Builds a Sims transversal system.
 *
 * For performance reason, here we takes a vector of generators by value.  The
 * Sims structured is always built on the earliest point not stabilized by the
 * group.
 */

template <typename P>
std::unique_ptr<SimsTransv<P>> build_sims_sys(size_t size, std::vector<P> gens)
{
    using Ptr2transv = std::unique_ptr<SimsTransv<P>>;

    Ptr2transv head;
    Ptr2transv* dest = &head;
    Ptr2transv curr;
    Point begin = 0;

    while ((curr = internal::build_sims_transv(begin, size, gens))) {
        begin = curr->get_target();
        *dest = std::move(curr);
        dest = &dest->next;
    }

    return std::move(head);
}

} // End namespace libcanon
#endif // LIBCANON_STRING_CANON_H
