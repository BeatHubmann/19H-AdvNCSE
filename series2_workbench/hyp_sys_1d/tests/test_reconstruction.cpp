#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <ancse/reconstruction.hpp>


TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    Eigen::VectorXd ua(3), ub(3);
    ua << 1.0, 0, 1.50;
    ub << 0.1, 0, 0.15;


    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}

TEST(TestPWLinear, Example) {
    auto rc = PWLinearReconstruction{MinMod{}};

    Eigen::VectorXd ua(3), ub(3), uc(3), ud(3), uL_true(3), uR_true(3);
    
    ua << 1.5,
          0.0,
          0.0;

    ub << 2.0,
          0.0,
          0.0;

    uc << 3.0,
          0.0,
          0.0;
    
    ud << 3.5,
          0.0,
          0.0;

    auto [uL, uR] = rc(ua, ub, uc, ud);

    uL_true << 2.25,
               0.0,
               0.0;

    uR_true << 2.75,
               0.0,
               0.0;

    ASSERT_EQ(uL, uL_true);
    ASSERT_EQ(uR, uR_true);
}