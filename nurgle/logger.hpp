#pragma once

#include <ostream>
#include <string>
#include <unordered_map>
#include <iostream>

#include <boost/format.hpp>

namespace nurgle
{

namespace extras
{

inline void applyvars(boost::format& fmt)
{
    ;  // do nothing
}

template<typename Tfirst_, typename... Trest_>
void applyvars(boost::format& fmt, Tfirst_ const& first, Trest_ const&... rest)
{
    fmt = fmt % first;
    applyvars(fmt, rest...);
}

template<typename Tarray_, std::size_t SIZE>
constexpr std::size_t lengthof(Tarray_ const (&array)[SIZE])
{
    return SIZE;
}

} // extras

class Logger
{
public:

    enum level
    {
        L_OFF = 0,
        L_DEBUG = 1,
        L_INFO = 2,
        L_WARNING = 3,
        L_ERROR = 4,
        L_FATAL = 5
    };

public:

    Logger(std::string const& name);
    Logger(Logger const&) = delete;  //XXX: noncopyable
    ~Logger();

    void level(enum level lv);
    enum level level() const;
    std::string stringize_error_level(enum level lv);

    std::string const name() const
    {
        return name_;
    }

    template <typename... Trest_>
    void log(enum level lv, char const* format, Trest_ const&... rest)
    {
        if (lv < level_)
        {
            return;
        }

        std::cerr << boost::format("%s: %-8s ") % name_ % stringize_error_level(lv);

        {
            boost::format fmt(format);
            extras::applyvars(fmt, rest...);
            std::cerr << fmt;
        }

        std::cerr << std::endl;
    }

    void flush()
    {
        std::cerr << std::flush;
    }

    template <typename... Trest_>
    void debug(char const* format, Trest_ const&... rest)
    {
        log(L_DEBUG, format, rest...);
    }

    template <typename... Trest_>
    void info(char const* format, Trest_ const&... rest)
    {
        log(L_INFO, format, rest...);
    }

    template <typename... Trest_>
    void warn(char const* format, Trest_ const&... rest)
    {
        log(L_WARNING, format, rest...);
    }

    template <typename... Trest_>
    void error(char const* format, Trest_ const&... rest)
    {
        log(L_ERROR, format, rest...);
    }

    template <typename... Trest_>
    void fatal(char const* format, Trest_ const&... rest)
    {
        log(L_FATAL, format, rest...);
    }

    static Logger& get_logger(std::string const& name);
    // static Logger& get_logger(std::string const& name)
    // {
    //     typedef std::unordered_map<std::string, std::unique_ptr<Logger>> loggers_type;
    //     static loggers_type loggers;
    //     std::pair<loggers_type::iterator, bool> i(
    //             loggers.insert(loggers_type::value_type(name, nullptr)));
    //     if (i.second)
    //     {
    //         std::unique_ptr<Logger> tmp(new Logger(name));
    //         (*i.first).second.swap(tmp);
    //     }
    //     return *(*i.first).second;
    // }

protected:

    std::string const name_;
    enum level level_;
};

inline void LOG_LEVEL(enum Logger::level lv)
{
    Logger::get_logger("nurgle").level(lv);
}

template <typename... Trest_>
inline void LOG_DEBUG(char const* format, Trest_ const&... rest)
{
    Logger::get_logger("nurgle").debug(format, rest...);
}

template <typename... Trest_>
inline void LOG_INFO(char const* format, Trest_ const&... rest)
{
    Logger::get_logger("nurgle").info(format, rest...);
}

template <typename... Trest_>
inline void LOG_WARN(char const* format, Trest_ const&... rest)
{
    Logger::get_logger("nurgle").warn(format, rest...);
}

template <typename... Trest_>
inline void LOG_ERROR(char const* format, Trest_ const&... rest)
{
    Logger::get_logger("nurgle").error(format, rest...);
}

template <typename... Trest_>
inline void LOG_FATAL(char const* format, Trest_ const&... rest)
{
    Logger::get_logger("nurgle").fatal(format, rest...);
}

} // nurgle
