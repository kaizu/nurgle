#include <iostream>

#include <nurgle/csv.hpp>
#include <nurgle/event_scheduler.hpp>
#include <nurgle/world.hpp>
#include <nurgle/logger.hpp>

#include <nurgle/metabolism.hpp>


void timecourse(
    std::string const& filename,
    nurgle::World const& w,
    nurgle::ode::ODESystem<nurgle::ode::michaelis_menten> const& system,
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

std::unique_ptr<nurgle::Event<nurgle::World>> generate_protein_event(double const dt, std::string const& pathto, std::string const& sep)
{
    using namespace nurgle;
    typedef EnzymaticChemicalReactionEvent<ode::mass_action> event_type;

    auto event = std::make_unique<event_type>(dt);

    {
        typedef utils::csv<std::string, std::string> csv_type;
        typedef ChemicalReaction ret_type;
        csv_type::for_each(pathto + sep + "gene_product_map.csv", [&event](csv_type::row_type&& x) {
            ret_type reaction;
            reaction.name = std::get<0>(x);
            reaction.left = decltype(reaction.left)(0);
            reaction.right = decltype(reaction.right)(1, std::make_tuple(std::get<1>(x), 1.0));
            reaction.forward = 1.0;  // synthesis
            reaction.reverse = 1.0;  // degradation
            reaction.enzymes = decltype(reaction.enzymes)(0);
            event->add_reaction(reaction);
            });
    }
    {
        typedef utils::csv<std::string, std::string> csv_type;
        typedef ChemicalReaction ret_type;
        csv_type::for_each(pathto + sep + "plexes.csv", [&event](csv_type::row_type&& x) {
                {
                    ret_type reaction;
                    reaction.name = std::get<0>(x);
                    reaction.left = utils::_csv<std::tuple<std::string, double>, ':'>::read(std::istringstream(std::get<1>(x)), ';');
                    reaction.right = decltype(reaction.right)(1, std::make_tuple(std::get<0>(x), 1.0));
                    reaction.forward = 1.0;  // binding
                    reaction.reverse = 0.0;  // unbinding
                    reaction.enzymes = decltype(reaction.enzymes)(0);
                    event->add_reaction(reaction);
                }
                {
                    ret_type reaction;
                    reaction.name = std::get<0>(x);
                    reaction.left = decltype(reaction.left)(1, std::make_tuple(std::get<0>(x), 1.0));
                    reaction.right = decltype(reaction.right)(0);
                    reaction.forward = 1.0;  // degradation
                    reaction.reverse = 0.0;
                    reaction.enzymes = decltype(reaction.enzymes)(0);
                    event->add_reaction(reaction);
                }
            });
    }

    return std::unique_ptr<Event<World>>(event.release());
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

    double const dt = 100.0;
    double const duration = (argc > 2 ? std::atof(argv[2]) : dt * 100);
    // std::string target("GLC_e");
    std::string target("LEU_c");
    // std::string target("Acetoacetyl_ACPs_c");

    EventScheduler<World> scheduler;

    auto const event_id1 = scheduler.insert(
        generate_enzymatic_chemical_reaction_event<ode::michaelis_menten>(dt, pathto + sep + "metabolism.csv"));
    auto const& system = scheduler.as<EnzymaticChemicalReactionEvent<ode::michaelis_menten>>(event_id1).system;

    auto const event_id2 = scheduler.insert(
        generate_protein_event(dt, pathto, sep));
    // auto const event_id2 = scheduler.insert(
    //     generate_enzymatic_chemical_reaction_event<ode::mass_action>(
    //         dt, read_gene_expressions(pathto + sep + "gene_product_map.csv")));

    scheduler.insert(generate_fixed_interval_callback_event<World>(
        "TimerEvent", dt, [&](World& w) -> std::vector<EventScheduler<World>::token_type> {
            timecourse("timecourse.csv", w, system);
            std::cout << "The current time is " << w.t << "." << std::endl;
            std::cout << target << " = " << w.pool.get(target) - 1.0 << std::endl;
            std::cout << "EG10048-MONOMER = " << w.pool.get("EG10048-MONOMER") << std::endl;
            std::cout << "ABC-56-CPLX = " << w.pool.get("ABC-56-CPLX") << std::endl;
            return {};
            }), -1);

    scheduler.update(w);
    timecourse("timecourse.csv", w, system, std::ios_base::out);

    // w.pool.update(target, 1.1);
    w.pool.update(target, 1.0 + 1e-1);

    scheduler.run(w, duration);

    // dump_pool(std::cout, w.pool);

    return 0;
}
