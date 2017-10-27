#include <iostream>

#include <nurgle/csv.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>

int main(int argc, char* argv[])
{
    using namespace nurgle;

    std::string const sep("/");
    std::string const pathto(argc > 0 ? argv[1] : ".");

    std::string const filename1("compounds.csv");
    std::vector<std::tuple<std::string, double>> const compounds
        = utils::csv<std::string, double>::read(pathto + sep + filename1);
    for (auto elem : compounds)
    {
        std::cout << std::get<0>(elem) << " : " << std::get<1>(elem) << std::endl;
    }

    std::string const filename2("metabolism.csv");
    std::vector<std::tuple<std::string, std::string, double>> const metabolism
        = utils::csv<std::string, std::string, double>::read(pathto + sep + filename2);
    for (auto elem : metabolism)
    {
        std::cout << std::get<0>(elem) << " : " << std::get<1>(elem) << " : " << std::get<2>(elem) << std::endl;
    }

    EventScheduler<World> scheduler;
    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", 10.0, [&](World& w) -> std::vector<std::string> {
            std::cout << "The current time is " << w.t << "." << std::endl;
            return {};
            }));

    World w;
    scheduler.run(w, 100.0);

    return 0;
}
