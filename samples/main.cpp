#include <iostream>

#include <nurgle/csv.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/logger.hpp>

#include <nurgle/metabolism.hpp>


void timecourse(std::string const& filename, nurgle::World const& w)
{
    std::ofstream ofs(filename, (w.t == 0 ? std::ios::out : std::ios::app));
    assert(ofs.is_open());

    if (w.t == 0)
    {
        for (auto const& name : w.pool.variables)
        {
            ofs << "," << name;
        }
        ofs << std::endl;
    }

    ofs << w.t;
    for (auto const& val : w.pool.values)
    {
        ofs << "," << val;
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
    w.pool.update("Acetoacetyl_ACPs_c", 1.1);

    EventScheduler<World> scheduler;

    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", 1.0, [&](World& w) -> std::vector<EventScheduler<World>::token_type> {
            timecourse("timecourse.csv", w);
            std::cout << "The current time is " << w.t << "." << std::endl;
            std::cout << "Acetoacetyl_ACPs_c" << " = " << w.pool.get("Acetoacetyl_ACPs_c") << std::endl;
            return {};
            }), -1);

    scheduler.insert(
        generate_enzymatic_chemical_reaction_event(1.0, pathto + sep + "metabolism.csv"));

    double const duration = (argc > 2 ? std::atof(argv[2]) : 300);
    scheduler.run(w, duration);

    // dump_pool(std::cout, w.pool);

    return 0;
}
