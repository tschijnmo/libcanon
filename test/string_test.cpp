/** Tests for string canonicalization.
 *
 * Here we are not only testing the string canonicalization problem, the code
 * for basic permutation facilities in perm.h, as well as the code for Sims
 * transversal systems in sims.h, are both tested here.
 */

#include <vector>

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
        : size(6)
        , cyclic({ 1, 2, 0, 4, 5, 3 })
        , transpose(std::vector<size_t>({ 1, 0, 2, 4, 3, 5 }),
              1) // Test another constructor.
        , gens{ cyclic, transpose }
        , corresp{ 3, 4, 5 }
    {
    }

    /** Size of the permutation domain.
     */

    size_t size;

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
    Simple_perm identity(size);
    EXPECT_EQ(identity.size(), size);
    EXPECT_EQ(identity.get_earliest_moved(), size);
}

/** Tests of the inversion of permutations.
 *
 * Here we use the involutory transposition for the testing.
 */

TEST_F(S3_test, inv_perm)
{

    // The unevaluated inversion.
    auto unevaled = ~transpose;
    EXPECT_EQ(unevaled.size(), size);
    EXPECT_EQ(unevaled, transpose);

    // Evaluations of the evaluated inversion.
    Simple_perm evaled(unevaled);
    EXPECT_EQ(evaled.size(), size);
    EXPECT_EQ(evaled, transpose);
    EXPECT_EQ(evaled, unevaled);

    // Since the equality operator actually mostly tests the pre-image array,
    // here we do additional tests for the image.

    auto test_image = [&](const auto& expr) {
        for (size_t i = 0; i < size; ++i) {
            EXPECT_EQ(expr << i, transpose << i);
        }
    };

    test_image(unevaled);
    test_image(evaled);
}

/** Tests of the product of permutations.
 *
 * Here we not only test simple products, the composition with inversion is
 * also thoroughly tested.
 */

TEST_F(S3_test, mult_perm)
{
    // Identity, for convenience.
    Simple_perm identity(size);

    // Identity test facility.
    auto expect_id = [&](const auto& expr) {
        Simple_perm evaled(expr);
        EXPECT_EQ(expr, identity);
        EXPECT_EQ(evaled, identity);
        EXPECT_EQ(expr, evaled);
        EXPECT_EQ(evaled.get_earliest_moved(), size);
    };

    // t is for the transposition, here we test t t = 1.
    expect_id(transpose | transpose);

    // c is for the cyclic generator, here we test c c c = 1.
    expect_id(cyclic | cyclic | cyclic);

    // Another test for t t t = t.
    auto ttt_unevaled = transpose | transpose | transpose;
    Simple_perm ttt_evaled(ttt_unevaled);
    EXPECT_EQ(ttt_unevaled, transpose);
    EXPECT_EQ(ttt_evaled, transpose);

    // Tests for the composition of multiplication and inverse.

    auto test_prod_inv = [&](const auto& op) {
        expect_id(op | ~op);
        expect_id(~op | op);
    };

    // Atomic ones.
    test_prod_inv(transpose);
    test_prod_inv(cyclic);

    // Non-atomic ones.
    test_prod_inv(cyclic | cyclic);
    test_prod_inv(~cyclic);
}

/** Tests chaining of permutations.
 *
 * Here we test the chaining of an indefinite number of permutations.
 */

TEST_F(S3_test, chain_perm)
{
    Simple_perm identity(size);

    // The cyclic generator, order 3.
    std::vector<const Simple_perm*> cyclics(4, &cyclic);
    auto c3_res = chain<Simple_perm>(size, cyclics.begin(), cyclics.end() - 1);
    EXPECT_EQ(c3_res, identity);
    EXPECT_EQ(c3_res.get_earliest_moved(), size);

    auto c4_res = chain<Simple_perm>(size, cyclics.begin(), cyclics.end());
    EXPECT_EQ(c4_res, cyclic);

    // The transpose generator, order 2.
    std::vector<const Simple_perm*> transposes(6, &transpose);
    auto t5_res
        = chain<Simple_perm>(size, transposes.begin(), transposes.end() - 1);
    EXPECT_EQ(t5_res, transpose);

    auto t6_res
        = chain<Simple_perm>(size, transposes.begin(), transposes.end());
    EXPECT_EQ(t6_res, identity);
    EXPECT_EQ(t6_res.get_earliest_moved(), size);
}
