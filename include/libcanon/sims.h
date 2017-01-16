/** \file
 *
 * Generic code for Sims transversal systems.
 *
 * The code in this section can be used not only for string canonicalization
 * problem, it can be used for other problem as well, since it is fully
 * general.
 */

#ifndef LIBCANON_SIMS_H
#define LIBCANON_SIMS_H

#include <cassert>
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

template <typename P> class Sims_transv {
public:
    //
    // Required by transversal interface
    //

    /** The type of permutations contained.
     */

    using Perm = P;

    /** Tests if the permutation could be in this subgroup.
     */

    bool has(const P& perm) const { return target_ == (perm << target_); }

    /** Gets the coset representative of a permutation.
     */

    const P* get_repr(const P& perm) const { return get_repr(perm >> target_); }

    /** Inserts a permutation into the transversal system.
     *
     * The given permutation is inserted only when there is no representative
     * for that coset yet.  The pointer to the inserted permutation will be
     * returned if the permutation is really inserted.
     */

    template <typename T> const P* insert(T&& perm)
    {
        Point label = perm >> target_;

        if (label == target_) {
            // Subset is always already represented by the implicit identity.
            return nullptr;
        }

        auto& slot = transv_[label];

        if (!slot) {
            slot = std::make_unique<P>(std::forward<T>(perm));
            return slot.get();
        } else {
            return nullptr;
        }
    }

    /** Gets the pointer to the next level of transversal system.
     */

    Sims_transv* next() const { return next_.get(); }

    /** Gets a new next level of the transversal system.
     */

    void set_next(std::unique_ptr<Sims_transv> new_next)
    {
        next_.swap(new_next);
        return;
    }

    /** Releases pointer to the next level of the transversal system.
     */

    std::unique_ptr<Sims_transv> release_next() { return std::move(next_); }

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
        Point label = *perm >> target_;
        if (label == target_)
            return;
        transv_[label] = std::move(perm);
    }

    /** Gets a representative for the coset moving the given point to target.
     */

    const P* get_repr(Point point) const
    {
        const auto& found = transv_[point];

        // Slight redundancy for clarity.
        if (found) {
            return found.get();
        } else {
            return nullptr;
        }
    }

    /** Creates a transversal system.
     */

    Sims_transv(size_t target, size_t size)
        : target_(target)
        , transv_(size)
    {
    }

    /** Gets the target of the Sims transversal.
     */

    Point target() const { return target_; }

    /** Gets the size of the permutation domain for the transversal.
     */

    size_t size() const { return transv_.size(); }

    /** Conjugates the group formed by the transversal system.
     *
     * After calling this method, the group formed by the transversal system
     * will be transformed from the original $G$ to $p^- G p$, where $p$ is the
     * given permutation.  Note that the placement of inversion is slightly
     * non-conventional.  This is for the convenient of application in
     * canonicalization problem.
     */

    template <typename E> void conj(const E& perm)
    {
        Point new_target = perm << target_;
        size_t size = transv_.size();
        Transv_container new_transv(size);

        for (size_t i = 0; i < size; ++i) {
            auto& repr = transv_[i];
            if (repr) {
                auto new_perm = std::make_unique<P>(~perm | *repr | perm);
                Point new_label = *new_perm >> target_;

                // new_perm is guaranteed to be non-identity.
                new_transv[new_label] = std::move(new_perm);
            }
        }

        transv_.swap(new_transv);
        if (next_)
            next_->conj(perm);
    }

private:
    /** The target point for the transversal system.
     *
     * The permutations are organized by the point that is moved into the
     * target.
     */

    const Point target_;

    using Transv_container = std::vector<std::unique_ptr<P>>;

    /** The actual container for transversals.
     */

    Transv_container transv_;

    /** Unique pointer to the next level of transversal.
     */

    std::unique_ptr<Sims_transv> next_;

    //
    // Iteration over a transversal.
    //

    /** Iterator type for the sims transversal system.
     *
     * This is a simple input iterator.
     */

    class Sims_transv_it {
    public:
        /** Constructs an iterator.
         */

        Sims_transv_it(size_t curr, const Sims_transv& transv)
            : curr_(curr)
            , transv_(&transv)
        {
            // Forward to a valid index unless we are building an end.
            if (curr_ < transv.size())
                find_next_present();
        }

        /** Increments the iterator.
         */

        Sims_transv_it& operator++()
        {
            ++curr_;
            find_next_present();
            return *this;
        }

        /** Deference the iterator.
         *
         * A constant reference to a permutation will be returned.
         */

        const P& operator*() { return *transv_->transv_[curr_]; }

        /** Compares iterators for equality.
         */

        bool operator==(const Sims_transv_it& sentinel)
        {
            return (curr_ == sentinel.curr_);
        }

        /** Compares for inequality.
         */

        bool operator!=(const Sims_transv_it& sentinel)
        {
            return (curr_ != sentinel.curr_);
        }

    private:
        /** The current index in the transversal vector.
         */

        size_t curr_;

        /** Pointer to the transversal system.
         */

        const Sims_transv* transv_;

        /** Increments the current index to one with a permutation.
         *
         * If the current index has a permutation, nothing will be done.  If no
         * more permutation is available, the index will be set to the end.
         */

        void find_next_present()
        {
            while (!transv_->get_repr(curr_) && curr_ < transv_->size()) {
                ++curr_;
            }
        }
    };

public:
    friend auto begin(const Sims_transv& container)
    {
        return Sims_transv_it(0, container);
    }

    friend auto end(const Sims_transv& container)
    {
        return Sims_transv_it(container.size(), container);
    }
};

//
// Schreier-Sims algorithm for forming a transversal system
//

namespace internal {

    /** Special container for generating set of permutations.
     *
     * In this container, the Jerrum filter is applied to all given
     * permutations so that the number of permutations is guaranteed to be at
     * most $n - 1$, where $n$ denotes the size of the permutation domain.
     *
     * This class subclasses normal vector over permutations, so that all the
     * iteration facility is immediately available.
     */

    template <typename P> class Jerrum_container : public std::vector<P> {
    public:
        /** Constructs an empty Jerrum container.
         */

        Jerrum_container(size_t size)
            : std::vector<P>()
            , graph_(size)
        {
            this->reserve(size - 1);
        }

        /**
         * Inserts a new generator.
         *
         * The new generator will be processed by Jerrum filter.
         */

        void insert_gen(P gen)
        {
            size_t size = graph_.size(); // Size of the permutation domain.
            assert(gen.size() == size);

            while (true) {
                Point earliest = gen.get_earliest_moved();
                if (earliest == size)
                    return; // Identity.
                Point dest = gen << earliest;

                std::vector<const P*> stack{};
                bool found_path = find_path(earliest, dest, stack);

                if (!found_path) {
                    // When no loop is found.
                    size_t idx = this->size(); // Number of perms we have.
                    this->push_back(std::move(gen));
                    graph_[earliest].push_back({ dest, idx });
                    break;
                } else {
                    // When a loop if formed by the new perm.
                    //
                    // The product of all permutations along the path.
                    auto prod = chain<P>(size, stack.begin(), stack.end());
                    // Here construction from an iterator of pointers are
                    // required.

                    P transfed(prod | ~gen);
                    gen = std::move(transfed);
                    continue;
                }
            } // End main loop.
        }

    private:
        /** Finds a path from the given source to the given destination.
         *
         * Note that this depth-first algorithm assumes that the graph cannot
         * be cyclic, which is always the case here.
         */

        bool find_path(
            Point src, Point dest, std::vector<const P*>& stack) const
        {
            if (src == dest)
                return true;

            for (const auto& i : graph_[src]) {
                stack.push_back(&(*this)[i.second]);
                if (find_path(i.first, dest, stack)) {
                    return true;
                } else {
                    stack.pop_back();
                }
            }

            return false;
        }

        //
        // Data fields
        //

        using Row = std::vector<std::pair<Point, size_t>>;
        using Graph = std::vector<Row>;

        /** The Jerrum connection graph.
         *
         * For each point, we store its destinations and index in the base
         * vector for the corresponding permutation.
         */

        Graph graph_;
    };

    /** Forms the Schreier generators.
     *
     * The given transversal is assumed to hold representative for each left
     * coset of the group generated by `gens`.
     */

    template <typename P>
    std::vector<P> form_schreier_gens(size_t size,
        const Sims_transv<P>& tentative, const std::vector<P>& gens)
    {
        Jerrum_container<P> result(size);

        // Add perm \bar{perm}^- to the result.
        //
        // TODO: Make sure the correctness of the theorem.
        auto add2result = [&](const auto& perm) {
            if (tentative.has(perm)) {
                // Special treatment since identity is never explicitly stored.
                result.insert_gen(perm);
            } else {
                result.insert_gen(perm | ~(*tentative.get_repr(perm)));
            }
        };

        for (const auto& i : gens) {
            add2result(i); // Identity is not in tentative.
            for (const auto& j : tentative) {
                add2result(j | i);
            }
        }

        return static_cast<std::vector<P>&&>(result);
    }

    /** Builds one level of Sims transversal system.
     *
     * The earliest non-stabilized point will be found and its orbit is used to
     * build a transversal system of the pointwise stabilizer subgroup of that
     * point.  Null will be returned when all points up to the given size are
     * stabilized.
     *
     * The transversal will be returned, and the generators will be modified to
     * the Schreier generators generating the stabilizer subgroup of the
     * selected point.
     */

    template <typename P>
    std::unique_ptr<Sims_transv<P>> build_sims_transv(
        Point begin, size_t size, std::vector<P>& gens)
    {

        // First we find the orbit and build a (right) transversal system by
        // the standard algorithm.

        std::vector<Point> orbit{};
        orbit.reserve(size);
        std::unique_ptr<Sims_transv<P>> tentative;

        for (size_t target = begin; target < size; ++target) {
            orbit.clear();
            orbit.push_back(target);

            tentative = std::make_unique<Sims_transv<P>>(target, size);

            //
            // Standard algorithm for orbits and coset representatives.
            //
            // Just here we need to be careful that we need the permutations to
            // permute the points in the orbit into the target.
            //

            for (size_t curr_idx = 0; curr_idx < orbit.size(); ++curr_idx) {
                auto curr_point = orbit[curr_idx];
                auto curr_perm = tentative->get_repr(curr_point);

                for (const auto& i : gens) {
                    auto new_point = i << curr_point;
                    auto repr = tentative->get_repr(new_point);
                    if (new_point == target || repr)
                        continue;

                    // Now we have a new point in orbit.
                    orbit.push_back(new_point);

                    if (curr_point == target) {
                        tentative->insert(~i);
                    } else {
                        tentative->insert(~i | *curr_perm);
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
            = form_schreier_gens(size, *tentative, gens);
        gens.swap(schreier_gens);

        return tentative;
    } // End function build_sims_transv

} // End namespace internal

/** Builds a Sims transversal system.
 *
 * For performance reason, here we takes a vector of generators by value.  The
 * Sims structured is always built on the earliest point not stabilized by the
 * group.
 *
 * A null pointer will be returned if the generators only generate the trivial
 * group.
 */

template <typename P>
std::unique_ptr<Sims_transv<P>> build_sims_sys(size_t size, std::vector<P> gens)
{
    using Transv = Sims_transv<P>;

    // The dummy head for the transversal system.
    Transv head(0, size);

    Transv* curr = &head;
    std::unique_ptr<Transv> new_transv;

    Point begin = 0;

    while ((new_transv = internal::build_sims_transv(begin, size, gens))) {
        begin = new_transv->target();
        curr->set_next(std::move(new_transv));
        curr = curr->next();
    }

    return head.release_next();
}

} // End namespace libcanon

#endif
