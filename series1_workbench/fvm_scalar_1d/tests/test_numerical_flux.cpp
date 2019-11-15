#include <gtest/gtest.h>

#include <ancse/numerical_flux.hpp>

template <class NumericalFlux>
void check_consistency(const NumericalFlux &nf) {
    auto model = make_dummy_model();

    auto to_check = std::vector<double>{1.0, 2.0, 3.3251, -1.332};
    for (auto u : to_check) {
        ASSERT_DOUBLE_EQ(model.flux(u), nf(u, u));
    }
}

TEST(TestCentralFlux, consistency) {
    auto model = make_dummy_model();
    auto central_flux = CentralFlux(model);

    check_consistency(central_flux);
}

TEST(TestRusanovFlux, consistency)
{
    auto model= make_dummy_model();
    auto rusanov_flux= RusanovFlux(model);

    check_consistency(rusanov_flux);
}

TEST(TestLxFFlux, consistency)
{
    auto model= make_dummy_model();
    const double dx{0.1};
    const double dt{0.1};
    auto lxf_flux= LxFFlux(model, dx, dt);
    check_consistency(lxf_flux);
}
