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

    EventScheduler<World> scheduler;

    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", 1.0, [&](World& w) -> std::vector<EventScheduler<World>::token_type> {
            std::cout << "The current time is " << w.t << "." << std::endl;
            return {};
            }));

    scheduler.insert(
        generate_enzymatic_chemical_reaction_event(1.0, pathto + sep + "metabolism.csv"));

    scheduler.run(w, 10.0);

    dump_pool(std::cout, w.pool);

    return 0;
}
