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

TEST(Test_triangle_graph, can_be_canonicalized)
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

/** Forms a star graph.
 *
 * A star graph has a root node that is connected to other nodes in the graph.
 */

template <typename It> Eldag make_star_graph(It begin, It end)
{
    Eldag star{};
    star.edges.insert(star.edges.end(), begin, end);
    for (size_t i = 0; i < star.edges.size() + 1; ++i)
        star.update_ia();
    return star;
}

/** Tests the canonicalization of a non-symmetric star graph.
 *
 * Here in the star graph, we have one root node and three child nodes.  There
 * is no symmetry in the valences of the root node.
 */

TEST(Test_non_symm_star_graph, can_be_canonicalized)
{

    std::vector<Point> children{ 1, 2, 3 };

    // The expected form of the canonical form.
    //
    // We strive to make the canonical form look like this since it is the most
    // natural form.

    Eldag expected_canon = make_star_graph(children.begin(), children.end());

    // Now we loop over some other possible forms.
    do {

        auto star = make_star_graph(children.begin(), children.end());

        Node_symms<Simple_perm> symms(star.size(), nullptr);
        auto res = canon_eldag(star, symms, [](auto i) { return 0; });

        // This graph should have no symmetry.
        EXPECT_FALSE(res.second);

        auto canon_form = act_eldag(res.first, star);
        EXPECT_EQ(canon_form, expected_canon);

    } while (std::next_permutation(children.begin(), children.end()));
}

/** Tests the canonicalization of a symmetric star graph.
 *
 * This test is very similar to the test of non-symmetric star graph, just here
 * we have symmetry in the valences of the root node.
 *
 * Significant code duplication can be seen.  However, due to the different
 * code segments that is covered by the two tests, it is put here as a separate
 * test.
 */

TEST(Test_symm_star_graph, can_be_canonicalized)
{

    std::vector<Point> children{ 1, 2, 3 };

    Eldag expected_canon = make_star_graph(children.begin(), children.end());

    // The symmetries.
    Node_symms<Simple_perm> symms(children.size() + 1, nullptr);
    auto node_symm = build_sims_sys<Simple_perm>(
        children.size(), { { 2, 0, 1 }, { 1, 0, 2 } });
    symms[0] = node_symm.get();

    do {

        auto star = make_star_graph(children.begin(), children.end());

        auto res = canon_eldag(star, symms, [](auto i) { return 0; });

        auto canon_form = act_eldag(res.first, star);
        EXPECT_EQ(canon_form, expected_canon);

        // We have automorphism in this graph now, it should be S3 among the
        // three child nodes.

        ASSERT_TRUE(res.second);
        const auto& transv1 = *res.second;
        ASSERT_TRUE(transv1.next());
        const auto& transv2 = *transv1.next();
        EXPECT_FALSE(transv2.next());

        auto target1 = transv1.target();
        Point_vec orbit1{ target1 };
        for (const auto& i : transv1) {
            orbit1.push_back(i >> target1);
        }
        std::sort(orbit1.begin(), orbit1.end());
        Point_vec expected_orbit1{ 1, 2, 3 };
        EXPECT_EQ(orbit1, expected_orbit1);

        auto target2 = transv2.target();
        EXPECT_NE(target1, target2);
        Point_vec orbit2{ target2 };
        for (const auto& i : transv2) {
            orbit2.push_back(i >> target2);
        }
        EXPECT_EQ(orbit2.size(), 2);

        // The second orbit should be the same as the first when the target of
        // the first orbit is added.
        orbit2.push_back(target1);
        std::sort(orbit2.begin(), orbit2.end());
        EXPECT_EQ(orbit2, expected_orbit1);

    } while (std::next_permutation(children.begin(), children.end()));
}
