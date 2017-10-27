#pragma once

#include <assert.h>
#include <vector>
#include <tuple>
#include <regex>
#include <string>
#include <fstream>
#include <sstream>

namespace nurgle
{

namespace utils
{

template <typename T_, std::size_t I = 0>
struct match_pattern {};

template <std::size_t I>
struct match_pattern<std::string, I>
{
    static std::string str(char delim)
    {
        return R"aaa(((?:"[^"]*")|(?:[^)aaa" + std::string(1, delim) + R"aaa(]*)))aaa";
    }

    static std::string get(std::match_results<std::string::const_iterator> const& sm)
    {
        return sm.str(I + 1);
    }
};

template <std::size_t I>
struct match_pattern<long, I>
{
    static std::string str(char delim)
    {
        return std::string(R"(([+-]?\d+))");
    }

    static long get(std::match_results<std::string::const_iterator> const& sm)
    {
        return std::stol(sm.str(I + 1));
    }
};

template <std::size_t I>
struct match_pattern<double, I>
{
    static std::string str(char delim)
    {
        return std::string(R"(([+-]?\d+[.]?\d*(?:[eE][+-])?\d*))");
    }

    static double get(std::match_results<std::string::const_iterator> const& sm)
    {
        return std::stod(sm.str(I + 1));
    }
};

template <std::size_t I>
struct match_pattern<bool, I>
{
    static std::string str(char delim)
    {
        return std::string(R"(([01]))");
    }

    static double get(std::match_results<std::string::const_iterator> const& sm)
    {
        return sm.str(I + 1) != "0";
    }
};

template <std::size_t I, typename T1_, typename T2_>
struct match_pattern<std::pair<T1_, T2_>, I>
{
    static std::string str(char delim)
    {
        return match_pattern<T1_, I>::str(delim) + delim + match_pattern<T2_, I + 1>::str(delim);
    }

    static std::pair<T1_, T2_> get(std::match_results<std::string::const_iterator> const& sm)
    {
        return std::make_pair(
            match_pattern<T1_, I>::get(sm),
            match_pattern<T2_, I + 1>::get(sm));
    }
};

template <std::size_t I>
struct match_pattern<std::tuple<>, I>
{
    static std::string str(char delim)
    {
        return "";
    }

    static std::tuple<> get(std::match_results<std::string::const_iterator> const& sm)
    {
        return std::make_tuple();
    }
};

template <std::size_t I, typename _Head, typename... _Tail>
struct match_pattern<std::tuple<_Head, _Tail...>, I>
{
    static std::string str(char delim)
    {
        std::string const s =
            match_pattern<_Head, I>::str(delim)
            + match_pattern<std::tuple<_Tail...>, I + 1>::str(delim);
        if (I == 0)
        {
            return s;
        }
        return delim + s;
    }

    static std::tuple<_Head, _Tail...> get(std::match_results<std::string::const_iterator> const& sm)
    {
        return std::tuple_cat(
            std::make_tuple(match_pattern<_Head, I>::get(sm)),
            match_pattern<std::tuple<_Tail...>, I + 1>::get(sm));
    }
};

template <typename T_, char delimiter = ','>
struct _csv
{
    typedef T_ row_type;

    static std::regex pattern;

    static std::vector<row_type> read(std::string const& filename, char lineterminator = '\n')
    {
        std::ifstream is(filename);
        assert(!is.fail());
        return read(std::move(is), lineterminator);
    }

    template <typename Tret_>
    static std::vector<Tret_> read(std::string const& filename, std::function<Tret_ (row_type&&)> const& op, char lineterminator = '\n')
    {
        std::ifstream is(filename);
        assert(!is.fail());
        return read(std::move(is), op, lineterminator);
    }

    static void for_each(std::string const& filename, std::function<void (row_type&&)> const& op, char lineterminator = '\n')
    {
        std::ifstream is(filename);
        assert(!is.fail());
        for_each(std::move(is), op, lineterminator);
    }

    static std::vector<row_type> read(std::istream&& is, char lineterminator = '\n')
    {
        return read<row_type>(std::move(is), [](row_type&& x) { return x; }, lineterminator);
    }

    template <typename Tret_>
    static std::vector<Tret_> read(std::istream&& is, std::function<Tret_ (row_type&&)> const& op, char lineterminator = '\n')
    {
        typedef Tret_ ret_type;
        std::vector<ret_type> retval;
        for_each(std::move(is), [&retval, &op](row_type&& x) { retval.push_back(op(std::move(x))); }, lineterminator);
        return std::move(retval);
    }

    static void for_each(std::istream&& is, std::function<void (row_type&&)> const& op, char lineterminator = '\n')
    {
        std::smatch sm;
        std::string line;

        std::size_t i = 0;
        while (getline(is, line, lineterminator))
        {
            ++i;

            if (!line.empty() && line[0] == '#')
            {
                continue;
            }
            else if (std::regex_match(line, sm, pattern))
            {
                op(match_pattern<row_type>::get(sm));
            }
            else
            {
                std::stringstream ss;
                ss << "An ill-formed line was given [Line " << i << "]";
                throw std::runtime_error(ss.str());
            }
        }
    }
};

template<typename T_, char delimiter>
std::regex _csv<T_, delimiter>::pattern = std::regex("^" + match_pattern<row_type>::str(delimiter) + R"([\n\r]*$)");

template <typename... Args_>
struct csv: public _csv<std::tuple<Args_...>>
{
    typedef _csv<std::tuple<Args_...>> base_type;
    typedef typename base_type::row_type row_type;
};

} // utils

} // nurgle
