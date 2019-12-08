#include <gtest/gtest.h>

#include <ancse/grid.hpp>

TEST(Grid, cell_center)
{
    double a = 2.0, b = 4.0;
    double dx = 2.0 / 10.0;

    int n_cells = 14;
    int n_ghost = 2;

    auto grid = Grid({a, b}, n_cells, n_ghost);

    // The cell-center of the first interior cell is
    //    a + 0.5 *dx
    ASSERT_DOUBLE_EQ(cell_center(grid, n_ghost), a + 0.5*dx);
}

TEST(Grid, ReferenceConversions)
{
    double a= -1.0, b= 1.0;

    int n_cells= 14;
    int n_ghost= 2;

    double dx= (b - a) / (n_cells - 2 * n_ghost);

    auto grid= Grid({a, b}, n_cells, n_ghost);

    double x1= -0.25;
    double x2=  2.0 / 3.0;
    double x3= 0.999999999;
    double x4= 0.0;

    // x1 is 75% into respective cell
    ASSERT_NEAR(reference_point(grid, x1), 0.75, 1e-10);

    // x2 is 33% into respective cell 
    ASSERT_NEAR(reference_point(grid, x2), 1.0/3.0, 1e-10);

    // x3 is at rightmost border of respective cell:
    ASSERT_NEAR(reference_point(grid, x3), 0.999999, 1e-6);

    // x4 is at leftmost border of respective cell:
    ASSERT_NEAR(reference_point(grid, x4), 0.0, 1e-10);

    // compound test with cell_center
     ASSERT_NEAR(reference_point(grid, cell_center(grid, n_ghost)), 0.5, 1e-10);

    // compound test with cell_center, cell_idx
     ASSERT_NEAR(reference_point(grid, cell_center(grid, cell_idx(grid, a + 0.5 * dx))), 0.5, 1e-10);
}

TEST(Grid, cell_idx)
{
    double a = 2.0, b = 4.0;
    double dx = 2.0 / 10.0;

    int n_cells = 14;
    int n_ghost = 2;

    auto grid = Grid({a, b}, n_cells, n_ghost);

    // x= 2.3 is in the 4th overall (including ghosts) grid cell with index 3
    double x= 2.3;
    ASSERT_EQ(cell_idx(grid, x), 3);

    // compound test with cell_center
    ASSERT_DOUBLE_EQ(cell_center(grid, cell_idx(grid, a + 0.5 * dx)), a + 0.5 * dx);
}