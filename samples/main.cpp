#include <iostream>

#include <nurgle/csv.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/logger.hpp>

#include <nurgle/metabolism.hpp>


int main(int argc, char* argv[])
{
    using namespace nurgle;
    LOG_LEVEL(Logger::L_DEBUG);

    std::string const sep("/");
    std::string const pathto(argc > 0 ? argv[1] : ".");

    World w;
    unsigned int const seed = 0;
    w.rng.seed(seed);

    read_pool(pathto + sep + "compounds.csv", w.pool);
    w.pool.update("Acetoacetyl_ACPs_c", 0.5);
    std::cout << "Acetoacetyl_ACPs_c" << " = " << w.pool.get("Acetoacetyl_ACPs_c") << std::endl;

    EventScheduler<World> scheduler;

    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", 10.0, [&](World& w) -> std::vector<EventScheduler<World>::token_type> {
            std::cout << "The current time is " << w.t << "." << std::endl;
            std::cout << "Acetoacetyl_ACPs_c" << " = " << w.pool.get("Acetoacetyl_ACPs_c") << std::endl;
            return {};
            }));

    scheduler.insert(
        generate_enzymatic_chemical_reaction_event(10.0, pathto + sep + "metabolism.csv"));

    scheduler.run(w, 100.0);

    // dump_pool(std::cout, w.pool);

    return 0;
}
