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

#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include <libcanon/canon.h>
#include <libcanon/perm.h>

namespace libcanon {

//
// Generic code for Sims transversal systems
// -----------------------------------------
//
// The code in this section can be used not only for string canonicalization
// problem, it can be used for other problem as well, since it is fully
// general.
//

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

    /** The type of permutations contained. */
    using Perm = P;

    /** Unique pointer to the next level of transversal */
    std::unique_ptr<Sims_transv> next;

    /** Tests if the permutation could be in this subgroup. */
    bool has(const P& perm) { return target == (perm << target); }

    /** Gets the coset representative of a permutation. */
    const P* get_repr(const P& perm) { return get_repr(perm >> target); }

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
        const auto& found = transv[point];

        // Slight redundancy for clarity.
        if (found) {
            return found.get();
        } else {
            return nullptr;
        }
    }

    /**
     * Creates a transversal system.
     */

    Sims_transv(size_t target, size_t size)
        : target{ target }
        , transv{ size }
    {
    }

    Point get_target() const { return target; }

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
        Point new_target = perm << target;
        size_t size = transv.size();
        Transv_container new_transv(size);

        for (size_t i = 0; i < size; ++i) {
            auto& repr = transv[i];
            if (repr) {
                auto new_perm = std::make_unique<P>(~perm | *repr | perm);
                Point new_label = *new_perm >> target;

                // new_perm is guaranteed to be non-identity.
                new_transv[new_label] = std::move(new_perm);
            }
        }

        transv.swap(new_transv);
        if (next)
            next->conj(perm);
    }

private:
    const Point target;
    using Transv_container = std::vector<std::unique_ptr<P>>;
    Transv_container transv;

    //
    // Iteration over a transversal.
    //

    /**
     * Iterator type for the sims transversal system.
     *
     * This is a simple forward iterator.
     */

    class Sims_transv_it {

        /** Constructs an iterator. */
        Sims_transv_it(size_t curr, const Sims_transv& transv)
            : curr{ curr }
            , transv{ transv }
        {
            if (curr == 0)
                find_next_present();
        }

        Sims_transv_it& operator++()
        {
            ++curr;
            find_next_present();
            return *this;
        }

        P& operator*() { return *(transv.transv[curr]); }

        bool operator==(const Sims_transv_it& sentinel)
        {
            return (curr == sentinel.curr);
        }

        bool operator!=(const Sims_transv_it& sentinel)
        {
            return (curr != sentinel.curr);
        }

    private:
        size_t curr;
        const Sims_transv& transv;

        /** Increments the current index to a one with permutation.
         *
         * If the current index has a permutation, nothing will be done.  If no
         * more permutation is available, the index will be set to the end.
         */

        void find_next_present()
        {
            while (!transv.get_repr(curr) && curr < transv.get_size()) {
                ++curr;
            }
        }
    };

public:
    friend Sims_transv_it begin(const Sims_transv& container)
    {
        return { 0, container };
    }

    friend Sims_transv_it end(const Sims_transv& container)
    {
        return { container.transv.size(), container };
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
     */

    template <typename P> class Jerrum_container : public std::vector<P> {
    public:
        /** Constructs the jerrum container */
        Jerrum_container(size_t size)
            : std::vector<P>{}
            , size{ size }
            , graph(size) // For C++11 compatibility.
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
            while (true) {
                Point earliest = gen.get_earliest_moved();
                if (earliest == size)
                    return; // Identity.
                Point dest = gen << earliest;

                std::vector<const P*> stack{};
                bool found_path = find_path(earliest, dest, stack);
                if (!found_path) {
                    // When no loop is found.
                    size_t idx = this->size();
                    this->push_back(gen);
                    graph[earliest].push_back({ dest, idx });
                    break;
                } else {
                    P prod{ size, stack.begin(), stack.end() };
                    P transfed{ prod | ~gen };
                    gen = std::move(transfed);
                    continue;
                }
            } // End main loop.
        }

    private:
        /**
         * Finds a path from the given source to the given destination.
         */

        bool find_path(
            Point src, Point dest, std::vector<const P*>& stack) const
        {
            if (src == dest)
                return true;

            for (const auto& i : graph[src]) {
                stack.push_back(&((*this)[i.second]));
                if (find_path(i.first, dest, stack)) {
                    return true;
                } else {
                    stack.pop_back();
                }
            }

            return false;
        }

        // For each point, we store its destinations and index in the base
        // vector for the corresponding permutation.
        using Row = std::vector<std::pair<Point, size_t>>;
        using Graph = std::vector<Row>;

        size_t size;
        Graph graph;
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
        Jerrum_container<P> result{ size };

        // Add perm \bar{perm}^- to the result.
        auto add2result = [&](const auto& perm) {
            if (tentative.has(perm)) {
                // Special treatment since identity is never explicitly stored.
                result.insert_gen(perm);
            } else {
                result.insert_gen(perm | ~tentative.get_repr(perm));
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
            = form_schreier_gens(size, *tentative, gens);
        gens.swap(schreier_gens);

        return std::move(tentative);
    } // End function build_sims_transv

} // End namespace internal

/** Builds a Sims transversal system.
 *
 * For performance reason, here we takes a vector of generators by value.  The
 * Sims structured is always built on the earliest point not stabilized by the
 * group.
 */

template <typename P>
std::unique_ptr<Sims_transv<P>> build_sims_sys(size_t size, std::vector<P> gens)
{
    using Ptr2transv = std::unique_ptr<Sims_transv<P>>;

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

//
// String canonicalization problem
// -------------------------------
//
// By using the code for Sims transversal system, the code in this section can
// be used for the canonicalization of strings under a given permutation group.
//

/** Cosets for Sims transversal system.
 *
 * A left coset for one level of  Sims transversal system is basically a
 * pointer to a permutation with a reference to the next level of Sims
 * transversal.  Here we also have a pointer to the previous level of coset
 * that the current coset is in.  The product of all these permutations gives a
 * permutation in the actual coset.
 */

template <typename P> class Sims_coset {
public:
    /** Constructs a Sims coset by expanding from a previous one.
     *
     * \param prev The previous level of Sims coset.
     * \param selected The point selected to be moved into the target.  Note
     * that it needs to be a point in the orbit.
     */

    Sims_coset(const Sims_coset& prev, Point selected)
        : prev{ *prev }
        , curr{ prev->next }
        , selected{ selected }
        , perm{ curr->get_repr(selected) }
        , next{ curr->next->get() }
    {
    }

    /**
     * Constructs a root coset for the group given as a Sims transversal.
     */

    Sims_coset(const Sims_transv<P>& group)
        : curr{ nullptr }
        , next{ *group }
    {
    }

    /** Gets the pre-image of a point by the permutation labelling this coset.
     *
     * This function is more suitable for the query of just a few points.  For
     * full construction of a permutation, use `get_a_perm`.
     */

    friend Point operator>>(const Sims_coset* coset, Point point)
    {
        for (; coset->curr != nullptr; coset = coset->prev) {
            if (coset->perm) {
                point = *coset->perm >> point;
            } // Else we assume we chose identity.
        }
        return point;
    }

    /**
     * Gets a permutation in the coset.
     *
     * Null pointer for identity permutation.
     */

    std::unique_ptr<P> get_a_perm() const
    {
        std::unique_ptr<P> res;
        std::unique_ptr<P> prod;

        for (const Sims_coset* coset = this; coset->curr != nullptr;
             coset = coset->prev) {
            if (coset->perm) {
                if (res) {
                    prod = std::make_unique<P>(*coset->perm | *res);
                    res.swap(prod);
                } else {
                    res = std::make_unique<P>(*coset->perm); // Copy.
                }
            }
        }
        return std::move(res);
    }

    /** Tests if two cosets are equal.
     *
     * For performance reason, here we just compare the selected point and
     * assume that all the rest matches.
     */

    bool operator==(const Sims_coset& other)
    {
        return this->selected == other.selected;
    }

    /** Computes the hash of the coset.
     *
     *  In the same vein as the equality comparison, here we just use the point
     *  moved into the target as the hash.
     */

    size_t hash() const { return selected; }

    /** Gets the selected point to be put into the target.
     *
     * Result for root coset is undefined.
     */

    Point get_selected() const { return selected; }

    /** Gets the previous coset.
     *
     * It can only be called on non-root coset.
     */

    const Sims_coset& get_prev() const { return *prev; }

    /** Gets the current transversal system.
     *
     * It can only be called on non-root coset.
     */

    const Sims_transv<P>& get_curr() const { return *curr; }

    /** Gets the next level of transversal system */

    const Sims_transv<P>* get_next() const { return next; }

private:
    /** Pointer to the previous level of coset, nullptr for first level. */
    const Sims_coset* prev;

    /** The current level of transversal. */
    const Sims_transv<P>* curr;

    /** The label for the coset. */
    Point selected;

    /** Reference to the permutation chosen by this level of coset. */
    const P* perm;

    /** Pointer to the next level of subgroup. */
    const Sims_transv<P>* next;
};

/** Refiner for string canonicalization problem based on Sims transversal.
 *
 * For a problem, in addition to the type of permutation `P`, we also need the
 * type of the combinatorial structure `S` and the comparator for the alphabet
 * in the string `C`
 *
 * `S` needs to be indexable and able to be constructed from a size.  The
 * comparator needs to implement `equals` and `is_less` methods for comparing
 * the contents in the string.
 */

template <typename S, typename P, typename C> class Sims_refiner {
public:
    //
    // Types required by the refiner protocol.
    //

    using Coset = Sims_coset<P>;
    using Perm = P;
    using Transv = Sims_transv<P>;
    using Structure = S;
    using Container = std::unordered_map<S, P>;

    /** Comparator for the alphabet.  */

    using Alpha_comp = C;

    template <typename T>
    Sims_refiner(size_t size, T&& comp)
        : size{ size }
        , comp{ std::forward<T>(comp) }
    {
    }

    //
    // Methods required by the refiner protocol
    //

    /** Tests if a given coset is a leaf coset. */
    bool is_leaf(const Coset& coset) const { return !coset.get_next(); }

    /** Refines a non-leaf coset */
    std::vector<Coset> refine(const Structure& obj, const Coset& coset) const
    {
        std::vector<Coset> children{};

        Point target = coset.get_next()->get_target();
        // Only valid for non-leaf cosets.

        // Find the points of minimum alphabet in the orbit.
        Point min = coset >> target;
        children.emplace_back(coset, target);

        for (const auto& perm : coset.get_curr()) {
            Point src = perm >> target;
            Point orig = coset >> src;
            if (comp.equals(orig, min)) {
                children.emplace_back(coset, src);
            } else if (comp.less(orig, min)) {
                min = orig;
                children.clear();
                children.emplace_back(coset, src);
            }
            // Do nothing when a greater point is found.
        }

        return children;
    }

    /** Gets a permutation in a leaf coset */
    Perm get_a_perm(const Coset& coset) const
    {
        std::unique_ptr<P> perm = coset.get_a_perm();
        if (perm) {
            return std::move(*perm);
        } else {
            return Perm{ size }; // Create the identity permutation.
        }
    }

    /** Acts a given permutation to a combinatorial object. */
    Structure act(const Perm& perm, const Structure& obj) const
    {
        Structure result(size);

        for (size_t i = 0; i < size; ++i) {
            result[perm << i] = obj[i];
        }

        return result;
    }

    /** Left multiplies a coset by a permutation */
    Coset left_mult(const Perm& perm, const Coset& coset) const
    {
        return { coset.get_prev(), perm >> coset.get_selected() };
    }

    auto get_transv(const Coset& upper, const Coset& lower) const
    {
        return std::move(
            std::make_unique<Transv>(size, lower >> lower.get_target()));
    }

private:
    size_t size;
    C comp;
};

/** Canonicalize the given string.
 *
 * This is the main driver function for string canonicalization problem.  Note
 * that in addition to the other constrains, the combinatorial object needs to
 * be hashable, equality comparable, and has order `<` defined, which is better
 * compatible with the given comparator for the alphabet is used.
 */

template <typename S, typename P, typename C>
std::pair<S, std::unique_ptr<Sims_transv<P>>> canon_string(
    const S& input, const Sims_transv<P>& group, C&& comp)
{
    using Refiner = Sims_refiner<S, P, C>;
    Refiner refiner{ input.size(), std::forward<C>(comp) };
    Sims_coset<P> whole_group(group);
    typename Refiner::Container container{};

    auto aut = add_all_candidates(refiner, input, whole_group, container);

    const auto& canon_form
        = std::min_element(container.begin(), container.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });

    aut->conj(canon_form.second);
    return { std::move(canon_form.first), std::move(aut) };
}

} // End namespace libcanon
#endif // LIBCANON_STRING_CANON_H
