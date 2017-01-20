/** Tests of the Eldag canonicalization facility.
 *
 * In this module, we have test cases for the Eldag canonicalization facility
 * for different  kinds of eldags, like those with or without symmetry on its
 * nodes, and highly symmetrical ones and highly unsymmetrical ones.
 */

#include <algorithm>
#include <memory>
#include <vector>

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

    // Here we just test the canonicalizing perm can be acted on the triangle.
    Eldag canon_form = act_eldag(res.first, triangle);

    // The automorphism group is the core of the symmetry testing.

    ASSERT_TRUE(res.second);
    const auto& transv = *res.second;
    EXPECT_FALSE(transv.next());
    Point target = transv.target();

    std::vector<Point> points{ 0, 1, 2 };
    points.erase(std::find(points.begin(), points.end(), target));
    ASSERT_EQ(points.size(), 2);

    size_t n_aut = 0;
    for (const auto& i : transv) {
        ++n_aut;

        // Make sure that they are all cyclic permutation.
        Point sel_point = i >> target;
        Point other_point = points[0] == sel_point ? points[1] : points[0];
        EXPECT_EQ(other_point, i >> sel_point);
        EXPECT_EQ(target, i >> other_point);
    }
    EXPECT_EQ(n_aut, 2);
}
