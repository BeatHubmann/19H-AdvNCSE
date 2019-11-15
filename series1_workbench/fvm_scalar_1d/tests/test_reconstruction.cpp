#include <gtest/gtest.h>

#include <ancse/reconstruction.hpp>


// Test by checking an easy example.
TEST(TestPWConstant, Example) {
    auto rc = PWConstantReconstruction{};

    double ua = 1.0, ub = 2.0;
    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub);
}

TEST(TestMinmod, Example)
{
    auto rc= minmodReconstruction{};

    Eigen::VectorXd a;
    double desired, result;
    int minPos;

    a= Eigen::VectorXd::Ones(10);
    desired= 1.0;
    result= rc.minmod(a);
    ASSERT_EQ(desired, result);
    
    a[0]= -1;
    desired= 0.0;
    result= rc.minmod(a);
    ASSERT_EQ(desired, result); 

    a= Eigen::VectorXd::Zero(10);
    desired= 0.0;
    result= rc.minmod(a);
    ASSERT_EQ(desired, result);

    a= Eigen::VectorXd::Random(10).cwiseAbs();
    a.cwiseAbs().minCoeff(&minPos);
    desired= a[minPos];
    result= rc.minmod(a);
    ASSERT_EQ(desired, result);

    double ua = 1.0, ub = 2.0;
    auto [uL, uR] = rc(ua, ub);

    ASSERT_EQ(uL, ua);
    ASSERT_EQ(uR, ub); 
}