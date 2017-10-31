#pragma once

#include <ostream>

#include <nurgle/defs.hpp>

namespace nurgle
{

struct Pool
{
    typedef std::string id_type;
    typedef double value_type;
    typedef std::vector<bool> bool_container_type;  //XXX: be careful

    std::vector<id_type> variables;
    std::vector<value_type> values;
    bool_container_type is_constant;

    Pool()
        : variables(), values(), is_constant()
    {
        ;
    }

    virtual ~Pool()
    {
        ;
    }

    size_t size() const
    {
        return values.size();
    }

    size_t update(id_type const& name, value_type const& value = 0, bool const constant = false)
    {
        auto const& it = std::find(variables.begin(), variables.end(), name);
        if (it != variables.end())
        {
            size_t const idx = std::distance(variables.begin(), it);
            assert(idx >= 0 && idx < values.size());
            values[idx] = value;
            is_constant[idx] = constant;
            return idx;
        }
        else
        {
            size_t const idx = variables.size();
            assert(values.size() == idx);
            assert(is_constant.size() == idx);
            variables.push_back(name);
            values.push_back(value);
            is_constant.push_back(constant);
            return idx;
        }
        throw std::runtime_error("never get here.");
    }

    value_type get(id_type const& name) const
    {
        auto const& it = std::find(variables.begin(), variables.end(), name);
        if (it == variables.end())
        {
            return 0.0;
        }
        return values[std::distance(variables.begin(), it)];
    }

    bool check_constant(id_type const& name) const
    {
        auto const& it = std::find(variables.begin(), variables.end(), name);
        if (it == variables.end())
        {
            return false;
        }
        return is_constant[std::distance(variables.begin(), it)];
    }
};

void read_pool(std::string const& filename, Pool& pool)
{
    for (auto const& elem : utils::csv<std::string, double, bool>::read(filename))
    {
        // LOG_DEBUG("%1% : %2% : %3%", std::get<0>(elem), std::get<1>(elem), std::get<2>(elem));
        pool.update(std::get<0>(elem), std::get<1>(elem), std::get<2>(elem));
    }
}

void dump_pool(std::ostream& out, Pool const& pool)
{
    for (size_t i = 0; i < pool.size(); ++i)
    {
        unsigned int c = (pool.is_constant[i] ? 1 : 0);
        out << pool.variables[i] << "," << pool.values[i] << "," << c << std::endl;
    }
}

template <typename Tcont_>
void dump_pool(std::ostream& out, Pool const& pool, Tcont_ const& names)
{
    for (Pool::id_type const& name : names)
    {
        unsigned int c = (pool.check_constant(name) ? 1 : 0);
        out << name << "," << pool.get(name) << "," << c << std::endl;
    }
}

} // nurgle
