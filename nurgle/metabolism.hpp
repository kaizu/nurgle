#pragma once

#include <nurgle/defs.hpp>
#include <nurgle/logger.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/ode.hpp>
#include <nurgle/utility.hpp>

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
    std::vector<std::string> enzymes;
};

std::vector<ChemicalReaction>
    read_chemical_reactions(std::string const filename)
{
    typedef utils::csv<std::string, std::string, std::string, std::string, std::string, std::string> csv_type;
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
            reaction.enzymes = utils::split(std::get<5>(x), ';');
            return reaction;
            });
}

template <typename Tratelaw>
struct EnzymaticChemicalReactionEvent: public Event<World>
{
    typedef Tratelaw ratelaw_type;
    typedef Event<World> base_type;
    typedef typename base_type::world_type world_type;
    typedef ChemicalReaction reaction_type;
    typedef ode::ODESystem<ratelaw_type> system_type;
    typedef typename system_type::pool_type pool_type;

    double dt;
    double t;
    std::vector<reaction_type> reactions;
    system_type system;
    bool dirty;

    EnzymaticChemicalReactionEvent(
        double const dt,
        std::vector<reaction_type> const& reactions)
        : dt(dt), t(0.0), reactions(reactions), system(), dirty(true)
    {
        ;
    }

    EnzymaticChemicalReactionEvent(
        double const dt,
        std::string const& filename)
        : dt(dt), t(0.0), reactions(read_chemical_reactions(filename)), system(), dirty(true)
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
        if (dirty)
        {
            system.generate_reactions(reactions, w.pool);
            assert(system.reactions.size() > 0);
            // synchronize(w);
            dirty = false;
        }
        assert(w.t == t);
        synchronize(w);
        return t + dt;
    }

    void interrupt(world_type& w) override
    {
        if (w.t <= t)
        {
            return;
        }

        _fire(w, w.t);
        t = w.t;
    }

    std::vector<std::string> fire(world_type& w) override
    {
        LOG_DEBUG("EnzymaticChemicalReactionEvent => %1%", w.t);
        return _fire(w, t + dt);
    }

    std::vector<std::string> _fire(world_type& w, double const next_time_)
    {
        assert(next_time_ >= t);

        //XXX: workaround
        if (next_time_ - t < 1e-100)
        {
            t = next_time_;
            return {};
        }

        // dump_pool(std::cout, w.pool);
        // system.dump_fluxes("fluxes.csv", w.pool, t);
        // system.dump_variables("compounds.csv", w.pool, t);

        system.integrate(w.pool, t, next_time_ - t);
        t = next_time_;
        synchronize(w);

        return {};
    }
};

template <typename Tratelaw, typename... Trest_>
std::unique_ptr<Event<World>> generate_enzymatic_chemical_reaction_event(Trest_ const& ... rest)
{
    return std::unique_ptr<Event<World>>(new EnzymaticChemicalReactionEvent<Tratelaw>(rest...));
}

}; // nurgle
