#pragma once

#include <nurgle/defs.hpp>
#include <nurgle/pool.hpp>

namespace nurgle
{

struct World
{
    double t = 0.0;
    Pool pool;
    rng_type rng;

    bool paranoiac = false;
};

rng_type& getrng(World& w)
{
    return w.rng;
}

double gett(World const& w)
{
    return w.t;
}

void sett(World& w, double t)
{
    w.t = t;
}

void diagnosis(World const& w)
{
    assert(gett(w) >= 0.0);
}

}; // nurgle
