#include <gtest/gtest.h>

#include <ancse/boundary_condition.hpp>

TEST(TestBoundaryCondition, Periodic) {
    int n_cells = 10;
    int n_ghost = 2;

    Eigen::VectorXd u(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        u[i] = i;
    }

    //// ANCSE_CUT_START_TEMPLATE
    //// ANCSE_COMMENT Test you boundary condition on a small example.
    //// ANCSE_COMMENT Off-by-one errors are very common here.

    auto bc = PeriodicBC(n_ghost);
    bc(u);

    ASSERT_DOUBLE_EQ(u[0], 6.0);
    ASSERT_DOUBLE_EQ(u[1], 7.0);

    for (int i = n_ghost; i < n_cells - n_ghost; ++i) {
        ASSERT_DOUBLE_EQ(u[i], i) << "Failed on i = " << i;
    }

    ASSERT_DOUBLE_EQ(u[8], 2.0);
    ASSERT_DOUBLE_EQ(u[9], 3.0);
    //// ANCSE_END_TEMPLATE
}

TEST(TestBoundaryCondition, Outflow) {
    int n_cells = 10;
    int n_ghost = 2;

    Eigen::VectorXd u(n_cells);
    for (int i = 0; i < n_cells; ++i) {
        u[i] = i;
    }

    //// ANCSE_CUT_START_TEMPLATE
    //// ANCSE_COMMENT Test you boundary condition on a small example.
    //// ANCSE_COMMENT Off-by-one errors are very common here.
    auto bc = OutflowBC(n_ghost);
    bc(u);

    ASSERT_DOUBLE_EQ(u[0], 2.0);
    ASSERT_DOUBLE_EQ(u[1], 2.0);

    for (int i = n_ghost; i < n_cells - n_ghost; ++i) {
        ASSERT_DOUBLE_EQ(u[i], i) << "Failed on i = " << i;
    }

    ASSERT_DOUBLE_EQ(u[8], 7.0);
    ASSERT_DOUBLE_EQ(u[9], 7.0);
    //// ANCSE_END_TEMPLATE
}
