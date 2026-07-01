// logger.ixx
module;
#include <unistd.h>   // for isatty
#include <cstdio>     // for fileno (C++ alternative)
#include <print>
#include <string>
#include <string_view>
#include <fstream>
#include <mutex>
#include <chrono>
#include <source_location>
#include <format>
export module logger;

// ANSI color codes
export{

namespace Color {
    constexpr std::string_view RESET   = "\033[0m";
    constexpr std::string_view RED     = "\033[31m";
    constexpr std::string_view GREEN   = "\033[32m";
    constexpr std::string_view YELLOW  = "\033[33m";
    constexpr std::string_view BLUE    = "\033[34m";
    constexpr std::string_view MAGENTA = "\033[35m";
    constexpr std::string_view CYAN    = "\033[36m";

    constexpr std::string_view WHITE   = "\033[37m";
    constexpr std::string_view BOLD    = "\033[1m";
    constexpr std::string_view DIM     = "\033[2m";
}

enum class LogLevel {
    TRACE, DEBUG, INFO, WARNING, ERROR, FATAL
};

class Logger {
private:
    std::ofstream file_;
    LogLevel min_level_ = LogLevel::DEBUG;
    bool console_output_ = true;
    bool file_output_ = false;
    std::mutex mutex_;

    // Determine if output is a terminal (supports colors)
    static bool is_terminal() {
        #ifdef _WIN32
            return false; // Windows console needs special handling
        #else
            return ::isatty(fileno(stdout));
        #endif
    }

    std::string_view level_to_string(LogLevel level) const {
        switch(level) {
            case LogLevel::TRACE:   return "-trace-";
            case LogLevel::DEBUG:   return "-debug-";
            case LogLevel::INFO:    return "--info-";
            case LogLevel::WARNING: return "--warn-";
            case LogLevel::ERROR:   return "-error-";
            case LogLevel::FATAL:   return "-fatal-";
        }
        return "UNKNOWN";
    }

    std::string_view level_to_color(LogLevel level) const {
        switch(level) {
            case LogLevel::TRACE:   return Color::BLUE;
            case LogLevel::DEBUG:   return Color::CYAN;
            case LogLevel::INFO:    return Color::GREEN;
            case LogLevel::WARNING: return Color::YELLOW;
            case LogLevel::ERROR:   return Color::RED;
            case LogLevel::FATAL:   return Color::MAGENTA;
        }
        return Color::RESET;
    }

    std::string format_timestamp() const {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;

        std::tm tm;
        #ifdef _WIN32
            localtime_s(&tm, &time);
        #else
            localtime_r(&time, &tm);
        #endif

        return std::format("{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
            tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
            tm.tm_hour, tm.tm_min, tm.tm_sec);
    }

    void write_to_console(LogLevel level, std::string_view message, bool color) const {
        if (!console_output_) return;

        if (color && is_terminal()) {
            std::println("{}[{:7s}]{} {}",
                level_to_color(level),
                level_to_string(level),
                Color::RESET,
                message);
        } else {
            std::println("[{:7s}] {}", level_to_string(level), message);
        }
    }

    void write_to_file(LogLevel level, std::string_view message) {
        if (!file_output_ || !file_.is_open()) return;

        file_ << std::format("[{}] [{:7s}] {}\n",
            format_timestamp(),
            level_to_string(level),
            message);
        file_.flush();
    }

public:
    Logger() = default;

    ~Logger() {
        if (file_.is_open()) {
            file_.close();
        }
    }

    void set_min_level(LogLevel level) { min_level_ = level; }
    void enable_console(bool enable) { console_output_ = enable; }

    void enable_file(const std::string& filename, bool truncate = true) {
        std::lock_guard lock(mutex_);
        if (file_.is_open()) file_.close();

        auto mode = truncate ? std::ios::out : std::ios::app;
        file_.open(filename, mode);
        file_output_ = file_.is_open();
    }

    void log(LogLevel level, std::string_view message,
             const std::source_location& loc = std::source_location::current()) {
        if (level < min_level_) return;

        std::lock_guard lock(mutex_);

        // Format with location info for DEBUG and above
        std::string formatted_message;
        formatted_message = std::string(message);

        write_to_console(level, formatted_message, true);
        write_to_file(level, formatted_message);

        if (level == LogLevel::FATAL) {
            std::terminate();
        }
    }

    template<typename... Args>
    void trace(std::format_string<Args...> fmt, Args&&... args) {
        log(LogLevel::TRACE, std::format(fmt, std::forward<Args>(args)...));
    }

    template<typename... Args>
    void debug(std::format_string<Args...> fmt, Args&&... args) {
        log(LogLevel::DEBUG, std::format(fmt, std::forward<Args>(args)...));
    }

    template<typename... Args>
    void info(std::format_string<Args...> fmt, Args&&... args) {
        log(LogLevel::INFO, std::format(fmt, std::forward<Args>(args)...));
    }

    template<typename... Args>
    void warn(std::format_string<Args...> fmt, Args&&... args) {
        log(LogLevel::WARNING, std::format(fmt, std::forward<Args>(args)...));
    }

    template<typename... Args>
    void error(std::format_string<Args...> fmt, Args&&... args) {
        log(LogLevel::ERROR, std::format(fmt, std::forward<Args>(args)...));
    }

    template<typename... Args>
    void fatal(std::format_string<Args...> fmt, Args&&... args) {
        log(LogLevel::FATAL, std::format(fmt, std::forward<Args>(args)...));
    }
};
}
// Global logger instance
export inline Logger pklog;

