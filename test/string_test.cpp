/** Tests for string canonicalization.
 *
 * Here we are not only testing the string canonicalization problem, the code
 * for basic permutation facilities in perm.h, as well as the code for Sims
 * transversal systems in sims.h, are both tested here.
 */

#include <algorithm>
#include <memory>
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

class Test_S3 : public ::testing::Test {
public:
    //
    // Data and setting up
    //

    /** Sets up the test fixture.
     */

    Test_S3()
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

    //
    // Utility methods
    //

    /** Acts a permutation on the first half on the whole string.
     *
     * The second half is assumed to follow the first half.
     */

    template <typename T>
    T act_by_win(const std::vector<size_t>& perm_win, const T& orig)
    {
        size_t half_size = size / 2;

        T result(size);
        for (size_t i = 0; i < half_size; ++i) {
            result[i] = orig[perm_win[i]];
            result[i + half_size] = orig[perm_win[i] + half_size];
        }

        return result;
    }
};

//
// Tests of basic permutation facility
// -----------------------------------
//

/** Tests of the basic methods of the perm class.
 */

TEST_F(Test_S3, has_basic_perm_methods)
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

TEST_F(Test_S3, can_be_inverted)
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

TEST_F(Test_S3, can_be_multiplied)
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

TEST_F(Test_S3, can_be_chained)
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

//
// Tests of the Schreier-Sims algorithm
// ------------------------------------
//

/** Tests the Schreier-Sims algorithm.
 *
 * Here we test the building of a Sims transversal system by using the
 * Schreier-Sims algorithm.
 */

TEST_F(Test_S3, can_be_formed_into_sims_system)
{
    auto sys = build_sims_sys(size, gens);

    // Flatten the transversals out for easy checking.
    using Transv = Sims_transv<Simple_perm>;
    std::vector<const Transv*> transvs{};
    for (Transv* i = sys.get(); i; i = i->next()) {
        transvs.push_back(i);
    }
    ASSERT_EQ(transvs.size(), 2);

    // Examine the first level of transversal system.
    const Transv& first = *transvs.front();
    Point target = 0;
    EXPECT_EQ(first.target(), target);
    Point expect = 1;
    for (const auto& i : first) {
        EXPECT_EQ(i >> target, expect);
        ++expect;
    }
    EXPECT_EQ(expect, 3);

    // Examine the second level.
    const Transv& second = *transvs.back();
    target = 1;
    EXPECT_EQ(second.target(), target);
    expect = 2;
    for (const auto& i : second) {
        EXPECT_EQ(i >> target, expect);
        ++expect;
    }
    EXPECT_EQ(expect, 3);
}

/** Tests the transversal adaptation algorithm.
 *
 * Here we test the transversal adaptation algorithm by adapting the Sims
 * transversal system from the Schreier-Sims algorithm into another subgroup
 * chain.
 */

TEST_F(Test_S3, can_be_adapted_to_another_subgroup_chain)
{
    using Transv = Sims_transv<Simple_perm>;

    auto orig = build_sims_sys(size, gens);

    // Here we simply use the subgroup chain of stabilizer of 2 > the
    // stabilizer of 1 and 2.
    Transv new_first(2, size);
    new_first.set_next(std::make_unique<Transv>(1, size));
    const Transv& new_second = *new_first.next();

    adapt_transv(*orig, new_first);

    // Examine the first level.
    Point target = 2;
    Point expect = 0;
    EXPECT_EQ(new_first.target(), target);
    for (const auto& i : new_first) {
        EXPECT_EQ(i >> target, expect);
        ++expect;
    }
    EXPECT_EQ(expect, target);

    // Examine the second level.
    target = 1;
    expect = 0;
    EXPECT_EQ(new_second.target(), target);
    for (const auto& i : new_second) {
        EXPECT_EQ(i >> target, expect);
        ++expect;
    }
    EXPECT_EQ(expect, target);
}

//
// Test of a canonicalization of strings
// -------------------------------------
//

/** Tests the canonicalization of a non-symmetrical string.
 *
 * Here the string to be canonicalized does not really contain any non-trivial
 * automorphism.  This is more of a test for the refinement function for
 * strings.
 */

TEST_F(Test_S3, canonicalizes_non_symmetric_string)
{
    using Structure = std::vector<char>;

    // The original form of the string.  The string should always be permuted
    // back to this form.
    Structure orig{ 'a', 'b', 'c', 'x', 'y', 'z' };

    // The current permutation of the first three point.
    std::vector<size_t> curr_perm{ 0, 1, 2 };

    // The isomorphism group.
    auto iso = build_sims_sys(size, gens);

    // Loop over all forms of the structure.
    do {

        // Assemble the input structure.
        auto input_form = act_by_win(curr_perm, orig);

        auto res = canon_string(input_form, *iso);

        const auto& canon_perm = res.first;
        const auto& aut = res.second;

        auto canon_form = act_string<Structure>(canon_perm, input_form);
        EXPECT_EQ(canon_form, orig);

        // No automorphism for this string.
        EXPECT_FALSE(aut);

    } while (std::next_permutation(curr_perm.begin(), curr_perm.end()));
}

/** Tests the canonicalization of a string with automorphism.
 *
 * Here the string to be canonicalized contains some automorphism.   This thus
 * checks the discovery of automorphisms.
 */

TEST_F(Test_S3, canonicalizes_symmetric_string)
{
    using Structure = std::vector<char>;

    Structure orig{ 'a', 'a', 'a', 'x', 'y', 'y' };

    std::vector<size_t> curr_perm{ 0, 1, 2 };

    // The isomorphism group.
    auto iso = build_sims_sys(size, gens);

    // Loop over all forms of the structure.
    do {

        // Assemble the input structure.
        auto input_form = act_by_win(curr_perm, orig);

        auto res = canon_string(input_form, *iso);

        const auto& canon_perm = res.first;
        const auto& aut = res.second;

        auto canon_form = act_string<Structure>(canon_perm, input_form);
        EXPECT_EQ(canon_form, orig);

        EXPECT_TRUE(aut);
        EXPECT_FALSE(aut->next());

        Point target = aut->target();
        auto aut_it = begin(*aut);
        auto aut_sentin = end(*aut);

        // Note that we have the automorphism of the canonical rather than the
        // original form for string canonicalization.

        if (target == 1) {
            EXPECT_EQ(*aut_it >> target, 2);
        } else if (target == 2) {
            EXPECT_EQ(*aut_it >> target, 1);
        } else {
            EXPECT_TRUE(false);
        }

        ++aut_it;
        EXPECT_EQ(aut_it, aut_sentin);

    } while (std::next_permutation(curr_perm.begin(), curr_perm.end()));
}
