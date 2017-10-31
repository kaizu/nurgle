#pragma once

#include <utility>
#include <unordered_map>
#include <boost/operators.hpp>

#include <nurgle/defs.hpp>

namespace nurgle
{

// struct World;
// rng_type& getrng(World&);
// double gett(World const&);
// void sett(World&, double);

template <typename _Tworld>
struct Event
{
    typedef _Tworld world_type;
    typedef std::string token_type;

    unsigned long long num_steps;

    Event()
        : num_steps(0)
    {
        ;
    }

    virtual ~Event()
    {
        ; // do nothing
    }

    virtual std::vector<token_type> fire(world_type& w) = 0;

    virtual std::vector<token_type> const accessors() const
    {
        return {};
    }

    virtual double draw_next_time(world_type& w)
    {
        return Inf;
    }

    virtual void interrupt(world_type& w)
    {
        ; // do nothing
    }
};

template <typename _Tworld>
struct StochasticEvent: public Event<_Tworld>
{
    typedef Event<_Tworld> base_type;
    typedef typename base_type::world_type world_type;

    virtual ~StochasticEvent()
    {
        ; // do nothing
    }

    virtual double propensity(world_type& w) = 0;

    double draw_next_time(world_type& w) final override
    {
        double const a = propensity(w);
        if (a <= 0.0)
        {
            return Inf;
        }
        double const rnd = std::uniform_real_distribution<double>(0.0, 1.0)(getrng(w));
        double const dt = log(1.0 / rnd) / a;
        return gett(w) + dt;
    }
};

template <typename _Tworld>
struct FixedIntervalCallbackEvent: public Event<_Tworld>
{
    typedef Event<_Tworld> base_type;
    typedef typename base_type::world_type world_type;
    typedef typename base_type::token_type token_type;

    typedef std::function<std::vector<token_type>(world_type&)> operator_type;
    std::string name;
    double next_time, dt;
    operator_type op;

    FixedIntervalCallbackEvent(std::string const name, double const dt, operator_type op)
        : name(name), next_time(0.0), dt(dt), op(op)
    {
        ;
    }

    ~FixedIntervalCallbackEvent()
    {
        ;
    }

    std::vector<token_type> fire(world_type& w) override
    {
        auto retval = op(w);
        next_time += dt;
        return retval;
    }

    double draw_next_time(world_type& w) override
    {
        return next_time;
    }
};

template <typename _Tworld>
std::unique_ptr<Event<_Tworld>> generate_fixed_interval_callback_event(
    std::string const name, double const dt,
    typename FixedIntervalCallbackEvent<_Tworld>::operator_type op)
{
    return std::unique_ptr<Event<_Tworld>>(
        new FixedIntervalCallbackEvent<_Tworld>(name, dt, op));
}

template <class ForwardIterator>
ForwardIterator min_element(
    ForwardIterator const& first, ForwardIterator const& last, rng_type& rng)
{
    if (first == last)
    {
        return first;
    }

    ForwardIterator it = first;
    ForwardIterator result = it++;
    std::vector<typename std::iterator_traits<ForwardIterator>::difference_type>
        indices(1, 0);
    for (; it != last; ++it)
    {
        if (*it < *result)
        {
            result = it;
            indices.clear();
            indices.push_back(std::distance(first, it));
        }
        else if (*it == *result)
        {
            indices.push_back(std::distance(first, it));
        }
    }

    if (indices.size() == 1)
    {
        return result;
    }

    size_t const idx = std::uniform_int_distribution<size_t>(0, indices.size() - 1)(rng);
    return std::next(first, indices[idx]);
}

template <typename _Tworld>
struct EventScheduler
{
    typedef Event<_Tworld> event_type;
    typedef typename event_type::world_type world_type;
    typedef typename event_type::token_type token_type;

    typedef unsigned int event_id_type;

    typedef struct EventQueue
        : private boost::less_than_comparable<
            boost::equality_comparable<EventQueue>>
    {
        event_id_type id;
        std::unique_ptr<event_type> event;
        int priority;
        double next_time;

        EventQueue(event_id_type id, std::unique_ptr<event_type>&& event, int priority, double next_time)
            : id(id), event(std::move(event)), priority(priority), next_time(next_time)
        {}

        bool operator<(EventQueue const& x) const
        {
            return (next_time == x.next_time ? priority > x.priority : next_time < x.next_time);
        }

        bool operator==(EventQueue const& x) const
        {
            return priority == x.priority && next_time == x.next_time;
        }
    } event_queue_type;

    std::vector<event_queue_type> events;

    std::unordered_map<token_type, std::vector<size_t>> dependencies_;

    size_t num_events() const
    {
        return events.size();
    }

    double next_time() const
    {
        if (events.size() == 0)
        {
            return Inf;
        }
        return (*min_element(events.begin(), events.end())).next_time;
    }

    void insert(std::unique_ptr<event_type> event, int priority = 0)
    {
        size_t const idx = events.size();

        for (auto const& acc : event->accessors())
        {
            auto it = dependencies_.find(acc);
            if (it == dependencies_.end())
            {
                dependencies_.insert(std::make_pair(acc, std::vector<size_t>{idx}));
            }
            else
            {
                (*it).second.push_back(idx);
            }
        }

        events.push_back(event_queue_type(0, std::move(event), priority, 0.0));
    }

    void run(world_type& w, double const duration)
    {
        double const upto = gett(w) + duration;
        update(w);
        while (step(w, upto))
        {
            ;
        }
    }

    bool step(world_type& w, double const upto)
    {
        if (gett(w) >= upto)
        {
            return false;
        }
        else if (events.size() == 0)
        {
            sett(w, upto);
            return false;
        }

        auto it = min_element(events.begin(), events.end(), getrng(w));
        if ((*it).next_time <= upto)
        {
            sett(w, (*it).next_time);
            std::vector<token_type> const mutators = (*it).event->fire(w);
            (*it).event->num_steps++;
            update(w, (*it), mutators);
            return true;
        }
        else
        {
            sett(w, upto);
            return false;
        }
    }

    void update(world_type& w)
    {
        for (size_t i = 0; i < events.size(); ++i)
        {
            events[i].event->interrupt(w);
            events[i].next_time = events[i].event->draw_next_time(w);
        }
    }

    void update(world_type& w, event_queue_type& queue, std::vector<token_type> const& mutators)
    {
        queue.event->interrupt(w);
        queue.next_time = queue.event->draw_next_time(w);

        for (token_type const& type : mutators)
        {
            // std::cout << "mutataor = '" << type << "'" << std::endl;
            auto it = dependencies_.find(type);
            if (it == dependencies_.end())
            {
                continue;
            }

            for (size_t const idx : (*it).second)
            {
                events[idx].event->interrupt(w);
                events[idx].next_time = events[idx].event->draw_next_time(w);
                // std::cout << "[" << idx << "] was updated. => " << events[idx].next_time << std::endl;
            }
        }

        // for (size_t i = 0; i < events.size(); ++i)
        // {
        //     auto pred = [&](std::string const type) -> bool {
        //         auto const accessors = events[i].event->accessors();
        //         return std::find(accessors.begin(), accessors.end(), type) != accessors.end(); };

        //     if (i == pos || std::any_of(mutators.begin(), mutators.end(), pred))
        //     {
        //         events[i].event->interrupt(w);
        //         events[i].next_time = events[i].event->draw_next_time(w);
        //     }
        // }
    }
};

} // nurgle
