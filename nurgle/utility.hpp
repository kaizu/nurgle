#pragma once

#include <sstream>

namespace nurgle
{

namespace utils
{

std::vector<std::string> split(std::string const& s, char delim, bool ignore = false)
{
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim))
    {
        if (not (ignore || item.empty()))
        {
            elems.push_back(std::move(item));
        }
    }
    return std::move(elems);
}

template<class InputIt, class T, class UnaryOperation>
inline T product(InputIt first, InputIt last, T init, UnaryOperation op)
{
    return std::accumulate(first, last, init, [&op](T const& a, T const& b) -> T {
        return std::multiplies<T>()(a, op(b));
        });
}

template<class InputIt, class T>
inline T product(InputIt first, InputIt last, T init)
{
    return std::accumulate(first, last, init, std::multiplies<T>());
}

} // utils

} // nurgle
