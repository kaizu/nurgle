#include "logger.hpp"

namespace nurgle
{

Logger::Logger(std::string const& name)
    : name_(name)
{
    // std::cout << "Logger named [" << name_ << "]" << " was created." << std::endl;
}

Logger::~Logger()
{
    // std::cout << "Logger named [" << name_ << "]" << " was destroyed." << std::endl;
}

void Logger::level(enum Logger::level lv)
{
    level_ = lv;
}

enum Logger::level Logger::level() const
{
    return level_;
}

std::string Logger::stringize_error_level(enum level lv)
{
    static std::string names[] = {
        "OFF",
        "DEBUG",
        "INFO",
        "WARN",
        "ERROR",
        "FATAL"
    };
    return static_cast<std::size_t>(lv) >= extras::lengthof(names) ? "???": names[lv];
}

Logger& Logger::get_logger(std::string const& name)
{
    typedef std::unordered_map<std::string, std::unique_ptr<Logger>> loggers_type;
    static loggers_type loggers;
    std::pair<loggers_type::iterator, bool> i(
            loggers.insert(loggers_type::value_type(name, nullptr)));
    if (i.second)
    {
        std::unique_ptr<Logger> tmp(new Logger(name));
        (*i.first).second.swap(tmp);
    }
    return *(*i.first).second;
}

} // nurgle
