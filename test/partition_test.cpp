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

