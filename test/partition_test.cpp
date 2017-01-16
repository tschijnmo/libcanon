/** Tests the facilities for partition.
 *
 * Here the basic facilities for partition, like refinement and
 * individualization, are tested.  This is a preparation for the actual tests
 * for eldag canonicalization.
 */

#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include <libcanon/partition.h>

using namespace libcanon;

/** Test fixture for the partition tests.
 *
 * Here we simply start with a trivial partition.
 */

class Partition_test : public ::testing::Test {
public:
    /** Sets up the test fixture.
     */

    Partition_test()
        : size(4)
        , trivial(size)
    {
    }

    /** The size of the trivial partition.
     */

    size_t size;

    /** The trivial partition.
     */

    Partition trivial;
};

/** Tests the trivial partition.
 *
 * Here we just have some have simple tests.  Most of the facilities are going
 * to be tested on non-trivial partitions.
 */

TEST_F(Partition_test, test_trivial)
{
    Point_vec trivial_vec(size);
    std::iota(trivial_vec.begin(), trivial_vec.end(), 0);
    Partition::Normal_form normal{ trivial_vec };

    EXPECT_EQ(trivial.get_normal_form(), normal);
    EXPECT_EQ(trivial.make_perm(), Simple_perm(trivial_vec));
    EXPECT_EQ(trivial, Partition(normal));
}

/** Tests refinement facility.
 *
 * In this test case, the refinement facility and other things are tested.
 */

TEST_F(Partition_test, test_refinement)
{

    // Swap the first two and the second two points.

    Partition refined(trivial);
    refined.split_by_key(trivial.get_first(), [](auto point) {
        if (point == 0 || point == 1) {
            return 1;
        } else {
            return 0;
        }
    });

    // First check the correctness of the partition.

    Partition::Normal_form expected_nf{ { 2, 3 }, { 0, 1 } };
    Partition expected_part(expected_nf);
    EXPECT_EQ(refined.get_normal_form(), expected_nf);
    EXPECT_EQ(refined, Partition(expected_nf));

    // Test some simple queries of the partition as a whole.

    EXPECT_FALSE(refined.is_discrete());
    EXPECT_EQ(refined.size(), size);
    const auto& pre_imgs = refined.get_pre_imgs();
    auto perm = refined.make_perm();
    for (Point i = 0; i < size; ++i) {
        for (Point j : { pre_imgs[i], perm >> i, refined >> i }) {
            if (i == 0 || i == 1) {
                EXPECT_TRUE(j == 2 || j == 3);
            } else {
                EXPECT_TRUE(j == 0 || j == 1);
            }
        }
    }

    // Test the looping and point properties.

    // Vector of cells.
    Point_vec forward_cells(refined.begin(), refined.end());
    Point_vec backward_cells(refined.rbegin(), refined.rend());

    // Expand the cells into all points.
    std::vector<std::vector<Point_vec>> cells{};

    for (const auto& i : { forward_cells, backward_cells }) {
        std::vector<Point_vec> curr{};
        for (auto j : i) {
            curr.emplace_back(refined.cell_begin(j), refined.cell_end(j));
        }
        ASSERT_EQ(curr.size(), 2);

        for (const auto& j : curr) {
            ASSERT_EQ(j.size(), 2);
            // All points in the same cell has the same colour.
            EXPECT_EQ(refined.get_colour(j[0]), refined.get_colour(j[1]));
            EXPECT_NE(j[0], j[1]);
        }

        cells.push_back(std::move(curr));
    }

    // The first cell should have lower colour for forward iteration.

    for (auto i : cells[0][0]) {
        for (auto j : cells[0][1]) {
            EXPECT_NE(i, j);
            EXPECT_LT(refined.get_colour(i), refined.get_colour(j));
        }
    }

    // The opposite for reverse iteration.

    for (auto i : cells[1][0]) {
        for (auto j : cells[1][1]) {
            EXPECT_NE(i, j);
            EXPECT_GT(refined.get_colour(i), refined.get_colour(j));
        }
    }
}
