/** Tests the facilities for partition.
 *
 * Here the basic facilities for partition, like refinement and
 * individualization, are tested.  This is a preparation for the actual tests
 * for eldag canonicalization.
 */

#include <numeric>

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
}
