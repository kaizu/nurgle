#pragma once

#include <cmath>
#include <ostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/odeint.hpp>

#include <nurgle/csv.hpp>
#include <nurgle/utility.hpp>


namespace nurgle
{

namespace ode
{

namespace odeint = boost::numeric::odeint;

struct michaelis_menten
{
    template <typename Treaction_>
    double operator()(
        std::vector<double> const& reactants,
        std::vector<double> const& products,
        std::vector<double> const& enzymes,
        double const volume, double const t,
        Treaction_ const& reaction)
    {
        double const forward = utils::product(reactants.begin(), reactants.end(), reaction.forward);
        double const reverse = utils::product(products.begin(), products.end(), reaction.reverse);
    
        double const Km = 1.0;
        double denom1 = 1.0;
        for (size_t i = 0; i < reactants.size(); ++i)
        {
            assert(reaction.reactant_coefficients[i] > 0);
            denom1 *= std::pow(1 + reactants[i] / Km, static_cast<int>(reaction.reactant_coefficients[i]));
        }
        double denom2 = 1.0;
        for (size_t i = 0; i < products.size(); ++i)
        {
            assert(reaction.product_coefficients[i] > 0);
            denom2 *= std::pow(1 + products[i] / Km, static_cast<int>(reaction.product_coefficients[i]));
        }
        assert(denom1 + denom2 - 1.0 != 0.0);
        return (forward - reverse) / (denom1 + denom2 - 1);
    }
};

struct mass_action
{
    template <typename Treaction_>
    double operator()(
        std::vector<double> const& reactants,
        std::vector<double> const& products,
        std::vector<double> const& enzymes,
        double const volume, double const t,
        Treaction_ const& reaction)
    {
        double const forward = utils::product(reactants.begin(), reactants.end(), reaction.forward);
        double const reverse = utils::product(products.begin(), products.end(), reaction.reverse);
        return forward - reverse;
    }
};

// template <typename Treaction_>
// double evaluate_reaction(
//     std::vector<double> const& reactants,
//     std::vector<double> const& products,
//     std::vector<double> const& enzymes,
//     double const volume, double const t,
//     Treaction_ const& reaction)
// {
//     // if (reaction.reverse == 0.0)
//     // {
//     //     if (reactants.size() != 0)
//     //         return utils::product(reactants.begin(), reactants.end(), reaction.forward, [](double const x){ return x / (x + 0.25); });
//     //     else
//     //         return utils::product(products.begin(), products.end(), reaction.forward, [](double const x){ return 1.0 / (x + 0.25); });
//     // }
//     // else if (reaction.forward == 0.0)
//     // {
//     //     if (products.size() != 0)
//     //        return utils::product(products.begin(), products.end(), -reaction.reverse, [](double const x){ return x / (x + 0.25); });
//     //     else
//     //        return utils::product(reactants.begin(), reactants.end(), -reaction.reverse, [](double const x){ return 1.0 / (x + 0.25); });
//     // }
//     // return 0.0;
// 
//     // double const forward = utils::product(reactants.begin(), reactants.end(), reaction.forward, [](double const x){ return x / (x + 0.25); });
//     // double const reverse = utils::product(products.begin(), products.end(), reaction.reverse, [](double const x){ return x / (x + 0.25); });
//     double const forward = utils::product(reactants.begin(), reactants.end(), reaction.forward);
//     double const reverse = utils::product(products.begin(), products.end(), reaction.reverse);
// 
//     double const Km = 1.0;
//     double denom1 = 1.0;
//     for (size_t i = 0; i < reactants.size(); ++i)
//     {
//         assert(reaction.reactant_coefficients[i] > 0);
//         // if (reactants[i] < 0.0) std::cout << "[t=" << t << "] reactants[i] == " << reactants[i] << std::endl;
//         // assert(1 + reactants[i] / Km >= 1.0);
//         // assert(std::pow(1 + reactants[i] / Km, static_cast<int>(reaction.reactant_coefficients[i])) >= 1.0);
//         denom1 *= std::pow(1 + reactants[i] / Km, static_cast<int>(reaction.reactant_coefficients[i]));
//     }
//     double denom2 = 1.0;
//     for (size_t i = 0; i < products.size(); ++i)
//     {
//         assert(reaction.product_coefficients[i] > 0);
//         // if (products[i] < 0.0) std::cout << "[t=" << t << "] products[i] == " << products[i] << std::endl;
//         denom2 *= std::pow(1 + products[i] / Km, static_cast<int>(reaction.product_coefficients[i]));
//     }
//     // assert(denom1 >= 1.0);
//     // assert(denom2 >= 1.0);
//     // assert(denom1 + denom2 - 1.0 > 0.0);
//     assert(denom1 + denom2 - 1.0 != 0.0);
//     return (forward - reverse) / (denom1 + denom2 - 1);
// 
//     // double forward = reaction.forward;
//     // if (forward > 0.0)
//     // {
//     //     for (size_t i = 0; i < reactants.size(); ++i)
//     //     {
//     //         assert(reaction.reactant_coefficients[i] > 0);
//     //         forward *= std::pow(reactants[i], reaction.reactant_coefficients[i]);
//     //     }
//     // }
// 
//     // double reverse = reaction.reverse;
//     // if (reverse > 0.0)
//     // {
//     //     for (size_t i = 0; i < products.size(); ++i)
//     //     {
//     //         assert(reaction.product_coefficients[i] > 0);
//     //         reverse *= std::pow(products[i], reaction.product_coefficients[i]);
//     //     }
//     // }
// 
//     // return forward - reverse;
// }

template <typename Tratelaw>
struct ODESystem
{
    typedef Tratelaw ratelaw_type;

    typedef enum
        {
            RUNGE_KUTTA_CASH_KARP54 = 1,
            ROSENBROCK4_CONTROLLER = 2,
            EULER = 3,
        } solver_type;

    typedef boost::numeric::ublas::vector<double> state_type;
    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef std::vector<state_type::size_type> index_container_type;
    typedef std::vector<double> coefficient_container_type;

    struct reaction_type
    {
        std::string name;
        index_container_type reactants;
        coefficient_container_type reactant_coefficients;
        index_container_type products;
        coefficient_container_type product_coefficients;
        index_container_type enzymes;
        double forward;
        double reverse;
    };

    typedef std::vector<reaction_type> reaction_container_type;
    typedef std::vector<double> state_container_type;

    typedef Pool pool_type;
    typedef pool_type::id_type id_type;
    typedef pool_type::value_type value_type;
    typedef pool_type::bool_container_type bool_container_type;

    struct deriv_func
    {
        reaction_container_type const& reactions;
        double const volume;
        bool_container_type const& is_constant;

        deriv_func(reaction_container_type const& reactions, double const& volume,
                   bool_container_type const& is_constant)
            : reactions(reactions), volume(volume), is_constant(is_constant)
        {
            ;
        }

        void operator()(state_type const& x, state_type& dxdt, double const& t) const
        {
            // for (size_t i = 0; i < x.size(); ++i)
            // {
            //     auto const val = x[i];
            //     if (isnan(val))
            //     {
            //         std::cout << i << " " << val << std::endl;
            //     }
            //     assert(isnan(val) == false);
            // }

            std::fill(dxdt.begin(), dxdt.end(), 0.0);

            for (auto const& reaction : reactions)
            {
                state_container_type reactant_state(reaction.reactants.size());
                state_container_type product_state(reaction.products.size());
                state_container_type enzyme_state(reaction.enzymes.size());

                std::size_t cnt;

                cnt = 0;
                for (auto const& idx : reaction.reactants)
                {
                    reactant_state[cnt] = x[idx];
                    cnt++;
                }

                cnt = 0;
                for (auto const& idx : reaction.products)
                {
                    product_state[cnt] = x[idx];
                    cnt++;
                }

                cnt = 0;
                for (auto const& idx : reaction.enzymes)
                {
                    enzyme_state[cnt] = x[idx];
                    cnt++;
                }

                double flux = ratelaw_type()(reactant_state, product_state, enzyme_state, volume, t, reaction);

                cnt = 0;
                for (auto const& idx : reaction.reactants)
                {
                    if (!is_constant[idx])
                    {
                        dxdt[idx] -= (flux * reaction.reactant_coefficients[cnt]);
                    }
                    cnt++;
                }

                cnt = 0;
                for (auto const& idx : reaction.products)
                {
                    if (!is_constant[idx])
                    {
                        dxdt[idx] += (flux * reaction.product_coefficients[cnt]);
                    }
                    cnt++;
                }
            }
        }
    };

    struct jacobi_func
    {
        reaction_container_type const& reactions;
        double const volume;
        bool_container_type const& is_constant;
        double const abs_tol, rel_tol;

        jacobi_func(reaction_container_type const& reactions, double const& volume,
                    bool_container_type const& is_constant,
                    double const abs_tol, double const rel_tol)
            : reactions(reactions), volume(volume), is_constant(is_constant),
            abs_tol(abs_tol), rel_tol(rel_tol)
        {
            ;
        }

        void operator()(
            state_type const& x, matrix_type& jacobi, double const& t, state_type& dfdt) const
        {
            std::fill(dfdt.begin(), dfdt.end(), 0.0);
            std::fill(jacobi.data().begin(), jacobi.data().end(), 0.0);

            double const SQRTETA(1.4901161193847656e-08);
            double const r0(1.0);
            // double fac(0.0);
            // for (std::size_t k(0); k < dfdt.size(); ++k)
            // {
            //     double const ewtk(atol + rtol * x[k]);
            //     fac = std::max(fac, dfdt[k] * ewtk);
            // }
            // double const r0(1000.0 * h * ETA * dfdt.size() * fac);  //XXX: h means the step interval
            // {
            //     double const ewtj(atol + rtol * x[j]);
            //     double const dyj(std::max(SQRTETA * abs(x[j]), r0 * ewtj));
            // }

            double const ht(1.0e-10);

            for (auto const& reaction : reactions)
            {
                index_container_type::size_type reactants_size(reaction.reactants.size());
                index_container_type::size_type products_size(reaction.products.size());
                state_container_type reactant_state(reactants_size);
                state_container_type product_state(products_size);
                state_container_type enzyme_state(reaction.enzymes.size());

                std::size_t cnt;

                cnt = 0;
                for (auto const& idx : reaction.reactants)
                {
                    reactant_state[cnt] = x[idx];
                    cnt++;
                }

                cnt = 0;
                for (auto const& idx : reaction.products)
                {
                    product_state[cnt] = x[idx];
                    cnt++;
                }

                cnt = 0;
                for (auto const& idx : reaction.enzymes)
                {
                    enzyme_state[cnt] = x[idx];
                    cnt++;
                }

                double const flux0 = ratelaw_type()(reactant_state, product_state, enzyme_state, volume, t, reaction);

                {
                    double const flux = ratelaw_type()(reactant_state, product_state, enzyme_state, volume, t + ht, reaction);
                    double const flux_deriv = (flux - flux0) / ht;

                    if (flux_deriv != 0.0)
                    {
                        for(std::size_t k(0); k < reaction.reactants.size(); k++)
                        {
                            matrix_type::size_type const row = reaction.reactants[k];
                            if (!is_constant[row])
                            {
                                double const coeff = reaction.reactant_coefficients[k];
                                dfdt[row] -= coeff * flux_deriv;
                            }
                        }

                        for(std::size_t k(0); k < reaction.products.size(); k++)
                        {
                            matrix_type::size_type const row = reaction.products[k];
                            if (!is_constant[row])
                            {
                                double const coeff = reaction.product_coefficients[k];
                                dfdt[row] += coeff * flux_deriv;
                            }
                        }
                    }
                }

                for (std::size_t j(0); j < reactant_state.size(); j++)
                {
                    double const ewt = abs_tol + rel_tol * std::abs(reactant_state[j]);
                    double const h = std::max(SQRTETA * std::abs(reactant_state[j]), r0 * ewt);
                    state_container_type h_shift(reactant_state);
                    h_shift[j] += h;
                    double const flux = ratelaw_type()(h_shift, product_state, enzyme_state, volume, t, reaction);
                    double const flux_deriv = (flux - flux0) / h;
                    matrix_type::size_type const col = reaction.reactants[j];

                    for(std::size_t k(0); k < reaction.reactants.size(); k++)
                    {
                        matrix_type::size_type const row = reaction.reactants[k];
                        if (!is_constant[row])
                        {
                            double const coeff = reaction.reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                    }

                    for(std::size_t k(0); k < reaction.products.size(); k++)
                    {
                        matrix_type::size_type const row = reaction.products[k];
                        if (!is_constant[row])
                        {
                            double const coeff = reaction.product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                }

                for (std::size_t j(0); j < product_state.size(); j++)
                {
                    double const ewt = abs_tol + rel_tol * std::abs(product_state[j]);
                    double const h = std::max(SQRTETA * std::abs(product_state[j]), r0 * ewt);
                    state_container_type h_shift(product_state);
                    h_shift[j] += h;
                    double const flux = ratelaw_type()(reactant_state, h_shift, enzyme_state, volume, t, reaction);
                    double const flux_deriv = (flux - flux0) / h;
                    matrix_type::size_type const col = reaction.products[j];

                    for(std::size_t k(0); k < reaction.reactants.size(); k++)
                    {
                        matrix_type::size_type const row = reaction.reactants[k];
                        if (!is_constant[row])
                        {
                            double const coeff = reaction.reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                    }

                    for(std::size_t k(0); k < reaction.products.size(); k++)
                    {
                        matrix_type::size_type const row = reaction.products[k];
                        if (!is_constant[row])
                        {
                            double const coeff = reaction.product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                }

                for (std::size_t j(0); j < enzyme_state.size(); j++)
                {
                    double const ewt = abs_tol + rel_tol * std::abs(enzyme_state[j]);
                    double const h = std::max(SQRTETA * std::abs(enzyme_state[j]), r0 * ewt);
                    state_container_type h_shift(enzyme_state);
                    h_shift[j] += h;
                    double const flux = ratelaw_type()(reactant_state, product_state, h_shift, volume, t, reaction);
                    double const flux_deriv = (flux - flux0) / h;
                    matrix_type::size_type const col = reaction.enzymes[j];

                    for(std::size_t k(0); k < reaction.reactants.size(); k++)
                    {
                        matrix_type::size_type const row = reaction.reactants[k];
                        if (!is_constant[row])
                        {
                            double const coeff = reaction.reactant_coefficients[k];
                            jacobi(row, col) -= coeff * flux_deriv;
                        }
                    }

                    for(std::size_t k(0); k < reaction.products.size(); k++)
                    {
                        matrix_type::size_type const row = reaction.products[k];
                        if (!is_constant[row])
                        {
                            double const coeff = reaction.product_coefficients[k];
                            jacobi(row, col) += coeff * flux_deriv;
                        }
                    }
                }
            }
        }
    };

    state_type state_init;
    reaction_container_type reactions;
    double volume;
    double abs_tol, rel_tol;

    ODESystem(double const abs_tol = 1e-12, double const rel_tol = 1e-8)
        : reactions(), volume(1.0), abs_tol(abs_tol), rel_tol(rel_tol)
    {
        ;
    }

    virtual ~ODESystem()
    {
        ;
    }

    template <typename Treaction_>
    void generate_reactions(std::vector<Treaction_> const& inputs, pool_type& pool)
    {
        std::unordered_map<id_type, state_type::size_type> index_map;
        state_type::size_type idx = 0;
        for (auto const& name : pool.variables)
        {
            index_map[name] = idx;
            idx++;
        }

        reactions.clear();
        reactions.reserve(inputs.size());

        for (Treaction_ const& r0 : inputs)
        {
            reaction_type r;
            r.name = r0.name;
            r.forward = r0.forward;
            r.reverse = r0.reverse;

            r.reactants.reserve(r0.left.size());
            r.products.reserve(r0.right.size());
            r.enzymes.reserve(r0.enzymes.size());

            for (auto const& cmp : r0.left)
            {
                auto const& it = index_map.find(std::get<0>(cmp));
                if (it != index_map.end())
                {
                    r.reactants.push_back(index_map[std::get<0>(cmp)]);
                }
                else
                {
                    size_t i = pool.update(std::get<0>(cmp));
                    r.reactants.push_back(i);
                    index_map[std::get<0>(cmp)] = i;
                }
                r.reactant_coefficients.push_back(std::get<1>(cmp));
            }

            for (auto const& cmp : r0.right)
            {
                auto const& it = index_map.find(std::get<0>(cmp));
                if (it != index_map.end())
                {
                    r.products.push_back(index_map[std::get<0>(cmp)]);
                }
                else
                {
                    size_t i = pool.update(std::get<0>(cmp));
                    r.products.push_back(i);
                    index_map[std::get<0>(cmp)] = i;
                }
                r.product_coefficients.push_back(std::get<1>(cmp));
            }

            for (auto const& name : r0.enzymes)
            {
                auto const& it = index_map.find(name);
                if (it != index_map.end())
                {
                    r.enzymes.push_back(index_map[name]);
                }
                else
                {
                    size_t i = pool.update(name);
                    r.enzymes.push_back(i);
                    index_map[name] = i;
                }
            }

            reactions.push_back(r);
        }
    }

    void synchronize(pool_type const& pool)
    {
        state_init.resize(pool.size());

        for (size_t idx = 0; idx < pool.size(); idx++)
        {
            state_init[idx] = static_cast<double>(pool.values[idx]);
        }
    }

    void integrate(pool_type& pool, double const t, double const dt)
    {
        if (pool.size() != state_init.size())
        {
            if (pool.size() > state_init.size())
            {
                state_init.resize(pool.size(), 0.0);
            }
            else
            {
                throw std::runtime_error("The size of the pool given is smaller than expected.");
            }
        }

        std::vector<double> timelog;
        std::vector<state_type> statelog;

        // synchronize(pool);

        size_t steps;

        double const subdt = dt * 0.1;

        // steps = odeint::integrate(
        //     deriv_func(reactions, volume, pool.is_constant), state_init, t, t + dt, subdt,
        //     [&](const state_type& state, const double t) {
        //         timelog.push_back(t);
        //         statelog.push_back(state);
        //         }
        //     );

        solver_type solver = RUNGE_KUTTA_CASH_KARP54;
        double const abs_tol(1e-6), rel_tol(1e-6);

        switch (solver)
        {
        case RUNGE_KUTTA_CASH_KARP54:
            {
                /* This solver doesn't need the jacobian */
                typedef odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
                steps = odeint::integrate_adaptive(
                    odeint::make_controlled<error_stepper_type>(abs_tol, rel_tol),
                    deriv_func(reactions, volume, pool.is_constant), state_init, t, t + dt, subdt,
                    [&](const state_type& state, const double t) {
                        timelog.push_back(t);
                        statelog.push_back(state);
                        }
                    );
            }
            break;
        case ROSENBROCK4_CONTROLLER:
            {
                typedef odeint::rosenbrock4<state_type::value_type> error_stepper_type;
                steps = odeint::integrate_adaptive(
                    odeint::make_controlled<error_stepper_type>(abs_tol, rel_tol),
                    std::make_pair(deriv_func(reactions, volume, pool.is_constant), jacobi_func(reactions, volume, pool.is_constant, abs_tol, rel_tol)),
                    state_init, t, t + dt, subdt,
                    [&](const state_type& state, const double t) {
                        timelog.push_back(t);
                        statelog.push_back(state);
                        }
                    );
            }
            break;
        case EULER:
            {
                typedef odeint::euler<state_type> stepper_type;
                steps = odeint::integrate_const(
                    stepper_type(), deriv_func(reactions, volume, pool.is_constant), state_init, t, t + dt, subdt,
                    [&](const state_type& state, const double t) {
                        timelog.push_back(t);
                        statelog.push_back(state);
                        }
                    );
            }
            break;
        default:
            throw std::runtime_error("The given solver is not supported\n");
        };

        assert(steps > 0);
        for (size_t idx = 0; idx < pool.size(); idx++)
        {
            // pool.values[idx] = statelog[steps][idx];
            pool.values[idx] += statelog[steps][idx] - statelog[0][idx]; //XXX: This may lead negative values.
            // if (pool.values[idx] < 0)
            // {
            //     std::cerr << idx << ": " << pool.variables[idx] << ": " << pool.values[idx] << " += " << statelog[steps][idx] << " - " << statelog[0][idx] << std::endl;
            // }
            // assert(pool.values[idx] >= 0);
        }
    }

    double evaluate(reaction_type const& reaction, double const t = 0.0) const
    {
        ODESystem::state_container_type reactant_state(reaction.reactants.size());
        ODESystem::state_container_type product_state(reaction.products.size());
        ODESystem::state_container_type enzyme_state(reaction.enzymes.size());

        std::size_t cnt;

        cnt = 0;
        for (auto const& idx : reaction.reactants)
        {
            reactant_state[cnt] = state_init[idx];
            cnt++;
        }

        cnt = 0;
        for (auto const& idx : reaction.products)
        {
            product_state[cnt] = state_init[idx];
            cnt++;
        }

        cnt = 0;
        for (auto const& idx : reaction.enzymes)
        {
            enzyme_state[cnt] = state_init[idx];
            cnt++;
        }

        return ratelaw_type()(reactant_state, product_state, enzyme_state, volume, t, reaction);
    }

    void dump_fluxes(std::ostream& out, double const t = 0.0) const
    {
        for (auto const& reaction : reactions)
        {
            out << reaction.name << "," << evaluate(reaction, t) << std::endl;
        }
    }
};

}; // ode

}; // nurgle
