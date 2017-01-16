/** Tests the facilities for partition.
 *
 * Here the basic facilities for partition, like refinement and
 * individualization, are tested.  This is a preparation for the actual tests
 * for eldag canonicalization.
 */

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

