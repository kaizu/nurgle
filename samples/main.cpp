#include <iostream>

#include <nurgle/csv.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/logger.hpp>

#include <nurgle/metabolism.hpp>


void timecourse(
    std::string const& filename, nurgle::World const& w, nurgle::ode::ODESystem const& system,
    std::ios_base::openmode mode = std::ios_base::app)
{
    std::ofstream ofs(filename, mode);
    assert(ofs.is_open());

    if (mode == std::ios_base::out)
    {
        ofs << "Time";
        for (auto const& name : w.pool.variables)
        {
            ofs << "," << name;
        }
        for (auto const& reaction : system.reactions)
        {
            ofs << "," << reaction.name;
        }
        ofs << std::endl;
    }

    ofs << w.t;
    for (auto const& val : w.pool.values)
    {
        ofs << "," << val;
    }
    for (auto const& reaction : system.reactions)
    {
        ofs << "," << system.evaluate(reaction, w.t);
    }
    ofs << std::endl;
}

int main(int argc, char* argv[])
{
    using namespace nurgle;
    LOG_LEVEL(Logger::L_DEBUG);

    std::string const sep("/");
    std::string const pathto(argc > 1 ? argv[1] : ".");

    World w;
    unsigned int const seed = 0;
    w.rng.seed(seed);

    read_pool(pathto + sep + "compounds.csv", w.pool);

    double const dt = 3.0;
    double const duration = (argc > 2 ? std::atof(argv[2]) : dt * 100);
    std::string target("Acetoacetyl_ACPs_c");

    EventScheduler<World> scheduler;

    auto const event_id1 = scheduler.insert(
        generate_enzymatic_chemical_reaction_event(dt, pathto + sep + "metabolism.csv"));
    auto const& system = scheduler.as<EnzymaticChemicalReactionEvent>(event_id1).system;

    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", dt, [&](World& w) -> std::vector<EventScheduler<World>::token_type> {
            timecourse("timecourse.csv", w, system);
            std::cout << "The current time is " << w.t << "." << std::endl;
            std::cout << target << " = " << w.pool.get(target) << std::endl;
            return {};
            }), -1);

    scheduler.update(w);
    timecourse("timecourse.csv", w, system, std::ios_base::out);

    w.pool.update(target, 1.1);

    scheduler.run(w, duration);

    // dump_pool(std::cout, w.pool);

    return 0;
}
