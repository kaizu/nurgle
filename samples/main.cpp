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

    for (auto const& elem : utils::csv<std::string, double, bool>::read(pathto + sep + "compounds.csv"))
    {
        LOG_DEBUG("%1% : %2% : %3%", std::get<0>(elem), std::get<1>(elem), std::get<2>(elem));
        w.pool.update(std::get<0>(elem), std::get<1>(elem), std::get<2>(elem));
    }

    {
        // LOG_INFO("%1% reactions were read from [%2%].", metabolism.size(), pathto + sep + "metabolism.csv");
        auto event1 = generate_enzymatic_chemical_reaction_event(1.0, read_chemical_reactions(pathto + sep + "metabolism.csv"));
        double const tnext = event1->draw_next_time(w);
        LOG_DEBUG("tnext = %1%", tnext);
        dynamic_cast<EnzymaticChemicalReactionEvent&>(*event1).system.dump_fluxes(std::cout);
    }

    EventScheduler<World> scheduler;
    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", 10.0, [&](World& w) -> std::vector<EventScheduler<World>::token_type> {
            std::cout << "The current time is " << w.t << "." << std::endl;
            return {};
            }));

    scheduler.run(w, 100.0);

    return 0;
}
