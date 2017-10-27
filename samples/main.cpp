#include <iostream>

#include <nurgle/csv.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/logger.hpp>


int main(int argc, char* argv[])
{
    using namespace nurgle;
    LOG_LEVEL(Logger::L_INFO);

    std::string const sep("/");
    std::string const pathto(argc > 0 ? argv[1] : ".");

    World w;

    for (auto elem : utils::csv<std::string, double>::read(pathto + sep + "compounds.csv"))
    {
        LOG_DEBUG("%s : %e", std::get<0>(elem), std::get<1>(elem));
        w.pool.update(std::get<0>(elem), std::get<1>(elem));
    }

    std::vector<std::tuple<std::string, std::string, double>> const metabolism
        = utils::csv<std::string, std::string, double>::read(pathto + sep + "metabolism.csv");
    for (auto elem : metabolism)
    {
        LOG_DEBUG("%s : %s : %e", std::get<0>(elem), std::get<1>(elem), std::get<2>(elem));
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
