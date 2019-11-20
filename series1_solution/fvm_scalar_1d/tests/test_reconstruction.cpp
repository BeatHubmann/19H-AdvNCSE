#include <gtest/gtest.h>

#include <ancse/reconstruction.hpp>

//// ANCSE_CUT_START_TEMPLATE
TEST(TestMinMod, Example) {
    ASSERT_DOUBLE_EQ(minmod(-1.0, -2.0), -1.0);
    ASSERT_DOUBLE_EQ(minmod(2.0, 1.0), 1.0);

    ASSERT_DOUBLE_EQ(minmod(-1.0, 3.0), 0.0);

    ASSERT_DOUBLE_EQ(minmod(0.0, -1.0), 0.0);
    ASSERT_DOUBLE_EQ(minmod(1.0, 0.0), 0.0);
}
//// ANCSE_END_TEMPLATE

// Test by checking an easy example.
TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    double ua = 1.0, ub = 2.0;
    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}

//// ANCSE_CUT_START_TEMPLATE
// Testing some simple examples is surprisingly effective. It helps to choose
// the example carefully such that the situation is 'generic'.
TEST(TestPWLinear, Example) {
    auto rc = PWLinearReconstruction{MinMod{}};

    double ua = 1.5, ub = 2.0, uc = 3.0, ud = 3.5;
    auto [uL, uR] = rc(ua, ub, uc, ud);

    ASSERT_DOUBLE_EQ(uL, 2.25);
    ASSERT_DOUBLE_EQ(uR, 2.75);
}
//// ANCSE_END_TEMPLATE
