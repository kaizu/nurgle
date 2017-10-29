#pragma once

#include <nurgle/defs.hpp>
#include <nurgle/logger.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/ode.hpp>

// #define NDEBUG

namespace nurgle
{

struct ChemicalReaction
{
    std::string name;
    std::vector<std::tuple<std::string, double>> left;
    std::vector<std::tuple<std::string, double>> right;
    double forward;
    double reverse;
};

std::vector<ChemicalReaction>
    read_chemical_reactions(std::string const filename)
{
    typedef utils::csv<std::string, std::string, std::string, std::string, std::string> csv_type;
    typedef ChemicalReaction ret_type;
    return csv_type::read<ret_type>(
        filename,
        [](csv_type::row_type&& x) {
            ret_type reaction;
            reaction.name = std::get<0>(x);
            reaction.left = utils::_csv<std::tuple<std::string, double>, ':'>::read(std::istringstream(std::get<1>(x)), ';');
            reaction.right = utils::_csv<std::tuple<std::string, double>, ':'>::read(std::istringstream(std::get<2>(x)), ';');
            reaction.forward = std::stod(std::get<3>(x));
            reaction.reverse = std::stod(std::get<4>(x));
            return reaction;
            });
}

struct EnzymaticChemicalReactionEvent: public Event<World>
{
    typedef Event<World> base_type;
    typedef typename base_type::world_type world_type;
    typedef ChemicalReaction reaction_type;
    typedef ode::ODESystem system_type;
    typedef system_type::pool_type pool_type;

    double dt;
    double t;
    std::vector<reaction_type> reactions;
    system_type system;

    EnzymaticChemicalReactionEvent(
        double const dt,
        std::vector<reaction_type> const& reactions)
        : dt(dt), t(0.0), reactions(reactions), system()
    {
        ;
    }

    EnzymaticChemicalReactionEvent(
        double const dt,
        std::string const& filename)
        : dt(dt), t(0.0), reactions(read_chemical_reactions(filename)), system()
    {
        ;
    }

    virtual ~EnzymaticChemicalReactionEvent()
    {
        ;
    }

    void synchronize(world_type& w)
    {
        // for (auto const& reaction : reactions)
        // {
        //     for (auto const& enzyme : reaction.enzymes)
        //     {
        //         w.pool.update(enzyme, w.entities.count(enzyme));
        //     }
        // }
        system.synchronize(w.pool);
    }

    double draw_next_time(world_type& w) override
    {
        if (system.reactions.size() == 0)
        {
            system.generate_reactions(reactions, w.pool);
            assert(system.reactions.size() > 0);
            synchronize(w);
        }
        return t + dt;
    }

    void interrupt(world_type& w) override
    {
        if (w.t <= t)
        {
            return;
        }

        _fire(w, w.t - t);
        t = w.t;
    }

    std::vector<std::string> fire(world_type& w) override
    {
        LOG_DEBUG("EnzymaticChemicalReactionEvent => %1%", w.t);
        return _fire(w, dt);
    }

    std::vector<std::string> _fire(world_type& w, double const dt_)
    {
        //XXX: workaround
        if (dt_ < 1e-100)
        {
            t += dt_;
            return {};
        }

        // dump_populations(w.pool);
        // system.dump_fluxes("fluxes.csv", w.pool, t);
        // system.dump_variables("compounds.csv", w.pool, t);

        system.integrate(w.pool, t, dt_);
        t += dt_;
        synchronize(w);

        return {};
    }

    // static inline void dump_populations(pool_type const& pool, std::vector<Pool::id_type> const& names)
    // {
    //     for (auto const& name : names)
    //     {
    //         LOG_ERROR("%1% => %2% %3%", name, pool.get(name), (pool.check_constant(name) ? "(const)" : ""));
    //     }
    // }

    // static inline void dump_populations(pool_type const& pool)
    // {
    //     std::vector<Pool::id_type>::const_iterator it1(pool.variables.begin());
    //     std::vector<Pool::value_type>::const_iterator it2(pool.values.begin());
    //     Pool::bool_container_type::const_iterator it3(pool.is_constant.begin());

    //     for (; it1 != pool.variables.end(); ++it1, ++it2, ++it3)
    //     {
    //         Pool::value_type const& value = (*it2);
    //         if (value != 0.0 && (*it1).find("@") != std::string::npos)
    //         {
    //             LOG_ERROR("%1% => %2% %3%", (*it1), value, (*it3 ? "(const)" : ""));
    //         }
    //     }
    // }
};

template <typename... Trest_>
std::unique_ptr<Event<World>> generate_enzymatic_chemical_reaction_event(Trest_ const& ... rest)
{
    return std::unique_ptr<Event<World>>(new EnzymaticChemicalReactionEvent(rest...));
}

}; // nurgle