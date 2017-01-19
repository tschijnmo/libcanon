/** Tests of the Eldag canonicalization facility.
 *
 * In this module, we have test cases for the Eldag canonicalization facility
 * for different  kinds of eldags, like those with or without symmetry on its
 * nodes, and highly symmetrical ones and highly unsymmetrical ones.
 */

#include <vector>
#include <memory>

#include <gtest/gtest.h>

#include <libcanon/eldag.h>

using namespace libcanon;

/** Tests the canonicalization facility on a simple triangle.
 *
 * Here all the nodes are unary, and the automorphism group should be the
 * cyclic C3 group.
 */

TEST(triangle_graph, can_be_canonicalized)
{
    // Build the graph.
    const size_t n_nodes = 3;
    Eldag triangle{};
    for (size_t i = 0; i < n_nodes; ++i) {
        triangle.edges.push_back((i + 1) % n_nodes);
        triangle.update_ia();
    }

    // Canonicalize the graph.
    Node_symms<Simple_perm> symms(n_nodes, nullptr);
    auto res = canon_eldag(triangle, symms, [](auto point) { return 0; });

    // TODO: Add checking.
}
