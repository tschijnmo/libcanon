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

class C3_test : public ::testing::Test {
public:
    /** Sets up the test fixture.
     */

    C3_test()
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
