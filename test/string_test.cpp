/** Tests for string canonicalization.
 *
 * Here we are not only testing the string canonicalization problem, the code
 * for basic permutation facilities in perm.h, as well as the code for Sims
 * transversal systems in sims.h, are both tested here.
 */

#include <gtest/gtest.h>

#include <libcanon/perm.h>
#include <libcanon/sims.h>
#include <libcanon/string.h>

using namespace libcanon;

//
// Test fixture
// ------------
//

/** Test fixture for S3 group.
 *
 * This is a slight complication of the S3 group for better test coverage of
 * different facilities in the code.
 *
 * 1. The permutation domain has six points, the last half of the points must
 * follow the first half, as in 3-body interaction tensor in many-body physics.
 *
 * 2. The permutation parity of the first three points are included as the
 * accompanied action.
 *
 */

class S3_test : public ::testing::Test {
public:
    /** Sets up the test fixture.
     */

    S3_test()
        : cyclic({ 1, 2, 0, 4, 5, 3 })
        , transpose(std::vector<size_t>({ 1, 0, 2, 4, 3, 5 }),
              1) // Test another constructor.
        , gens{ cyclic, transpose }
        , corresp{ 3, 4, 5 }
    {
    }

    /** The cyclic permutation of the three points.
     */

    Simple_perm cyclic;

    /** Transposition of the first two points.
     */

    Simple_perm transpose;

    /** Generators of the full group.
     */

    std::vector<Simple_perm> gens;

    /** Mapping from the handle points to the slave points.
     */

    std::vector<size_t> corresp;
};

//
// Tests of basic permutation facility
// -----------------------------------
//

/** Tests of the basic methods of the perm class.
 */

TEST_F(S3_test, perm_methods)
{

    // Modular arithmetic of the point labels.
    //
    // Here we hard code the result so that we can avoid the problems with
    // unsigned integer arithmetic.

    auto plus1 = [](Point point) -> Point {
        if (point == 2) {
            return 0;
        } else if (point == 5) {
            return 3;
        } else {
            return point + 1;
        }
    };

    auto minus1 = [](Point point) -> Point {
        if (point == 0) {
            return 2;
        } else if (point == 3) {
            return 5;
        } else {
            return point - 1;
        }
    };

    // Test the pre-image and image operator by using the cyclic permutation.
    for (size_t i = 0; i < 6; ++i) {
        EXPECT_EQ(cyclic >> i, plus1(i));
        EXPECT_EQ(cyclic << i, minus1(i));
    }

    // Test the accompanied action.
    EXPECT_EQ(cyclic.acc(), 0);
    EXPECT_EQ(transpose.acc(), 1);

    // Test the size query.
    for (const auto& i : { cyclic, transpose }) {
        EXPECT_EQ(i.size(), 6);
    }

    // Test the detection of the first moved point.
    for (const auto& i : { cyclic, transpose }) {
        EXPECT_EQ(i.get_earliest_moved(), 0);
    }

    // Another test for identity permutations.
    size_t test_size = 6;
    Simple_perm identity(test_size);
    EXPECT_EQ(identity.size(), test_size);
    EXPECT_EQ(identity.get_earliest_moved(), test_size);
}
