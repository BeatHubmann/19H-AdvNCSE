#include <gtest/gtest.h>

#include <ancse/numerical_flux.hpp>

template <class NumericalFlux>
void check_consistency_euler(const NumericalFlux &nf) {

    auto model = Euler();
    int n_vars = model.get_nvars();

    Eigen::MatrixXd u(n_vars, 3);
    u.col(0) << 1, 0, 1.5;
    u.col(1) << 0.125, 0.25, 0.40;
    u.col(2) << 0.11432, -0.11432, 0.26432;

    double TOL = 1E-10;
    for (int k=0; k<u.cols(); k++) {
        ASSERT_LE( (model.flux(u.col(k)) - nf(u.col(k), u.col(k))).norm(), TOL);
    }
}

TEST(TestCentralFlux, consistency) {
    std::shared_ptr<Model> model = std::make_shared<Euler>();
    auto central_flux = CentralFlux(model);

    check_consistency_euler(central_flux);
}

TEST(TestRusanov, consistency)
{
    std::shared_ptr<Model> model= std::make_shared<Euler>();
    auto rusanov_flux= Rusanov(model);

    check_consistency_euler(rusanov_flux);
}

TEST(TestLaxFriedrichs, consistency)
{
    std::shared_ptr<Model> model= std::make_shared<Euler>();
    auto grid = Grid({0.0, 1.0}, 14, 2);
    auto sim_time = std::make_shared<SimulationTime>(0.0, 0.1, 0.2, 0);
    auto lxf_flux= LaxFriedrichs(grid, model, sim_time);

    check_consistency_euler(lxf_flux);
}

TEST(TestRoe, consistency)
{
    std::shared_ptr<Model> model= std::make_shared<Euler>();
    auto roe_flux= Roe(model);

    check_consistency_euler(roe_flux);
}

TEST(TestHLL, consistency)
{
    std::shared_ptr<Model> model= std::make_shared<Euler>();
    auto HLL_flux= HLL(model);

    check_consistency_euler(HLL_flux);
}

TEST(TestHLLc, consistency)
{
    std::shared_ptr<Model> model= std::make_shared<Euler>();
    auto HLLc_flux= HLLc(model);

    check_consistency_euler(HLLc_flux);
}