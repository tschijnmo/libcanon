/**
 * General utilities for permutation groups.
 */


#ifndef LIBCANON_PERM_H
#define LIBCANON_PERM_H

#include<vector>
#include<unordered_map>


//
// The container of DAG for storing permutation groups
// ===================================================
//


template<typename G, typename CP>
class Perm_decomp;  // Forward declaration.


/**
 * Permutation groups given as directed cyclic graph.
 *
 * This DAG is designed to be the container to hold permutation groups given by
 * multiple towers of subgroups and coset transversals.  Essentially, the nodes
 * in this DAG are the subgroups of the given permutation group.  And each node
 * could contain possibly multiple subgroups of it along with a transversal over
 * all of its cosets.  Here a pair of a subgroup and the corresponding
 * transversal is called a decomposition.  This is a generalization of the tower
 * of subgroup data structure commonly seen in computational group theory.
 *
 * For generality, this container is fully generic.  The subgroup can be
 * specified in any suitable way, and the transversals can contain either
 * permutations or some description of the cosets.
 *
 */

template<typename G, typename T>
struct Perm_DAG {
    using Group_type = G;
    using Transv_type = T;
    using Decomp_type = Perm_decomp<G, T>;
    using Decomps_type = std::vector<Decomp_type>;
    using Node_type = decltype(nodes)::value_type;
    using Nodes_type = std::unordered_map<G, Decomps_type>;

    Nodes_type nodes;

    // Copying would invalidate all the references, hence disabled.
    Perm_DAG(Perm_DAG&) = delete;
    Perm_DAG& operator=(Perm_DAG&) = delete;
};


template<typename G, typename T>
struct Perm_decomp {
    const Perm_DAG<G, T>::Node_type& child;
    std::vector<T> transv;
};

#endif //LIBCANON_PERM_H
