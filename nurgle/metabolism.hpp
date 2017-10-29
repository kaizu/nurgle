#pragma once

#include <nurgle/defs.hpp>
#include <nurgle/logger.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/ode.hpp>

// #define NDEBUG

namespace nurgle
{

struct EnzymaticChemicalReactionEvent: public Event<World>
{
    typedef Event<World> base_type;
    typedef typename base_type::world_type world_type;
    typedef ode::ChemicalReaction reaction_type;
    typedef ode::ODESystem system_type;
    typedef system_type::pool_type pool_type;

    static inline void dump_populations(pool_type const& pool, std::vector<Pool::id_type> const& names)
    {
        for (auto const& name : names)
        {
            LOG_ERROR("%1% => %2% %3%", name, pool.get(name), (pool.check_constant(name) ? "(const)" : ""));
        }
    }

    static inline void dump_populations(pool_type const& pool)
    {
        std::vector<Pool::id_type>::const_iterator it1(pool.variables.begin());
        std::vector<Pool::value_type>::const_iterator it2(pool.values.begin());
        Pool::bool_container_type::const_iterator it3(pool.is_constant.begin());

        for (; it1 != pool.variables.end(); ++it1, ++it2, ++it3)
        {
            Pool::value_type const& value = (*it2);
            if (value != 0.0 && (*it1).find("@") != std::string::npos)
            {
                LOG_ERROR("%1% => %2% %3%", (*it1), value, (*it3 ? "(const)" : ""));
            }
        }
    }

    double dt;
    double t;
    std::vector<reaction_type> reactions;
    system_type system;

    EnzymaticChemicalReactionEvent(
        double const dt,
        std::vector<reaction_type> const& reactions)
        : dt(dt), t(0.0), reactions(reactions), system()
    {
        // system.generate_reactions(reactions);
    }

    virtual ~EnzymaticChemicalReactionEvent()
    {
        ;
    }

    void set_state(world_type& w)
    {
        // for (auto const& reaction : reactions)
        // {
        //     for (auto const& enzyme : reaction.enzymes)
        //     {
        //         w.pool.update(enzyme, w.entities.count(enzyme));
        //     }
        // }
        system.set_state(w.pool);
    }

    double draw_next_time(world_type& w) override
    {
        if (system.reactions.size() == 0)
        {
            system.generate_reactions(w.pool, reactions);
            assert(system.reactions.size() > 0);
            set_state(w);
        }
        return t + dt;
    }

    // bool is_mutated(std::string const type) const override
    // {
    //     return false;
    // }

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

        // std::string const filename = "metabolites.csv";
        // std::ofstream ofs(filename, std::ios::app);
        // assert(ofs.is_open());
        // ofs << t
        //     << "," << system.pool.get("GLC-6-P@CCO-PERI-BAC")
        //     << "," << system.pool.get("GLC-6-P@CCO-CYTOSOL")
        //     << "," << system.pool.get("FRUCTOSE-6P@CCO-CYTOSOL")
        //     << "," << system.pool.get("FRUCTOSE-16-DIPHOSPHATE@CCO-CYTOSOL")
        //     << "," << system.pool.get("DIHYDROXY-ACETONE-PHOSPHATE@CCO-CYTOSOL")
        //     << "," << system.pool.get("GAP@CCO-CYTOSOL")
        //     << "," << system.pool.get("DPG@CCO-CYTOSOL")
        //     << "," << system.pool.get("G3P@CCO-CYTOSOL")
        //     << "," << system.pool.get("2-PG@CCO-CYTOSOL")
        //     << "," << system.pool.get("PHOSPHO-ENOL-PYRUVATE@CCO-CYTOSOL")
        //     << "," << system.pool.get("PYRUVATE@CCO-CYTOSOL")
        //     << "," << system.pool.get("PYRUVATE@CCO-PERI-BAC")
        //     << std::endl;

        system.dump_fluxes("fluxes.csv", w.pool, t);
        system.dump_variables("compounds.csv", w.pool, t);

        system.integrate(w.pool, t, dt_);
        t += dt_;
        set_state(w);

        // {
        //     std::ofstream ofs("timecourse.csv", std::ios::out | std::ios::app);
        //     assert(ofs.is_open());
        //     ofs << t << ",";
        //     ofs << w.pool.get("APS@CCO-CYTOSOL") << ",";
        //     ofs << w.pool.get("PAPS@CCO-CYTOSOL") << ",";
        //     ofs << w.pool.get("SO3@CCO-CYTOSOL") << ",";
        //     ofs << w.entities.count("PAPSSULFOTRANS-CPLX") << ",";
        //     ofs << w.entities.count("RED-THIOREDOXIN-MONOMER") << ",";
        //     ofs << w.entities.count("OX-THIOREDOXIN-MONOMER") << ",";
        //     ofs << w.entities.count("ADENYLYLSULFKIN-CPLX") << std::endl;
        // }

        return {};
    }
};

template <typename... Trest_>
std::unique_ptr<Event<World>> generate_enzymatic_chemical_reaction_event(Trest_ const& ... rest)
{
    return std::unique_ptr<Event<World>>(new EnzymaticChemicalReactionEvent(rest...));
}

}; // nurgle
