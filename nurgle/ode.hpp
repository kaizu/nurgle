#pragma once

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

double evaluate_reaction(
    std::vector<double> const& left,
    std::vector<double> const& right,
    std::vector<double> const& enzymes,
    double const volume, double const t, bool const reversible)
{
    // // if (enzymes.size() == 0)
    // // {
    // //     return 0.0;
    // // }
    // if (utils::product(enzymes.begin(), enzymes.end(), 1.0) == 0.0)
    // {
    //     return 0.0;
    // }

    // double const k = pow(0.1, left.size());
    // double const flux = utils::product(left.begin(), left.end(), k);
    // if (!reversible)
    // {
    //     return flux;
    // }
    // return flux - utils::product(right.begin(), right.end(), k);

    double const k = 1.0;
    double const flux = utils::product(
        left.begin(), left.end(), k, [](double const& conc) {
            double const Km = 1.0;
            return conc / (Km + conc);
            });
    if (!reversible)
    {
        return flux;
    }
    return flux - utils::product(
        right.begin(), right.end(), k, [](double const& conc) {
            double const Km = 1.0;
            return conc / (Km + conc);
            });
}

struct ODESystem
{
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
        index_container_type reactants;
        coefficient_container_type reactant_coefficients;
        index_container_type products;
        coefficient_container_type product_coefficients;
        index_container_type enzymes;
        bool reversible;
        std::string name;
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

                double flux = evaluate_reaction(reactant_state, product_state, enzyme_state, volume, t, reaction.reversible);

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

        jacobi_func(reaction_container_type const& reactions, double const& volume,
                    bool_container_type const& is_constant)
            : reactions(reactions), volume(volume), is_constant(is_constant)
        {
            ;
        }

        void operator()(
            state_type const& x, matrix_type& jacobi, double const& t, state_type& dfdt) const
        {
            std::fill(dfdt.begin(), dfdt.end(), 0.0);
            std::fill(jacobi.data().begin(), jacobi.data().end(), 0.0);

            const double h(1.0e-8);
            const double ht(1.0e-10);

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

                double const flux0 = evaluate_reaction(reactant_state, product_state, enzyme_state, volume, t, reaction.reversible);

                {
                    double const flux = evaluate_reaction(reactant_state, product_state, enzyme_state, volume, t + ht, reaction.reversible);
                    double const flux_deriv = (flux - flux0) / h;

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
                    state_container_type h_shift(reactant_state);
                    h_shift[j] += h;
                    double const flux = evaluate_reaction(h_shift, product_state, enzyme_state, volume, t, reaction.reversible);
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
                    state_container_type h_shift(product_state);
                    h_shift[j] += h;
                    double const flux = evaluate_reaction(reactant_state, h_shift, enzyme_state, volume, t, reaction.reversible);
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
                    state_container_type h_shift(enzyme_state);
                    h_shift[j] += h;
                    double const flux = evaluate_reaction(reactant_state, product_state, h_shift, volume, t, reaction.reversible);
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

    ODESystem()
        : reactions(), volume(1.0)
    {
        ;
    }

    virtual ~ODESystem()
    {
        ;
    }

    template <typename Treaction_>
    void generate_reactions(pool_type& pool, std::vector<Treaction_> const& org)
    {
        std::unordered_map<id_type, state_type::size_type> index_map;
        state_type::size_type idx = 0;
        for (auto const& name : pool.variables)
        {
            index_map[name] = idx;
            idx++;
        }

        reactions.clear();
        reactions.reserve(org.size());

        for (Treaction_ const& r0 : org)
        {
            reaction_type r;
            r.reversible = r0.reversible;
            r.name = r0.name;

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

    void set_state(pool_type const& pool)
    {
        state_init.resize(pool.size());

        for (size_t idx = 0; idx < pool.size(); idx++)
        {
            state_init[idx] = static_cast<double>(pool.values[idx]);
        }
    }

    void integrate(pool_type& pool, double const t, double const dt)
    {
        std::vector<double> timelog;
        std::vector<state_type> statelog;

        // set_state(pool);

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
                    std::make_pair(deriv_func(reactions, volume, pool.is_constant), jacobi_func(reactions, volume, pool.is_constant)),
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

    void dump_variables(std::string const& filename, pool_type const& pool, double const t) const
    {
        std::ofstream ofs(filename, std::ios::out);
        assert(ofs.is_open());

        ofs << "#t=" << t << std::endl;

        for (size_t idx = 0; idx < pool.size(); idx++)
        {
            ofs << pool.variables[idx] << "," << pool.values[idx] << "," << (pool.is_constant[idx] ? 1 : 0) << std::endl;
        }
    }

    void dump_fluxes(std::string const& filename, pool_type const& pool, double const t) const
    {
        ODESystem::state_type x(pool.size());
        for (size_t idx = 0; idx < pool.size(); idx++)
        {
            x[idx] = static_cast<double>(pool.values[idx]);
        }

        std::ofstream ofs(filename, std::ios::out);
        assert(ofs.is_open());

        ofs << "#t=" << t << std::endl;

        for (auto const& reaction : reactions)
        {
            ODESystem::state_container_type reactant_state(reaction.reactants.size());
            ODESystem::state_container_type product_state(reaction.products.size());
            ODESystem::state_container_type enzyme_state(reaction.enzymes.size());

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

            double flux = evaluate_reaction(reactant_state, product_state, enzyme_state, volume, t, reaction.reversible);
            // if (flux != 0.0)
            {
                ofs << reaction.name << "," << flux << std::endl;
            }
        }
    }
};

}; // ode

}; // nurgle
