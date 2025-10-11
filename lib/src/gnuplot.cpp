#include "gnuplot.hpp"

#include <cstdio>
#include <stdexcept>
#include <vector>

GnuplotCommand::GnuplotCommand(const std::string& gnuplotPath)
    : gnuplotPath_(gnuplotPath), fp_(nullptr, pclose) {
    openGnuplotProcess();
}

GnuplotCommand::~GnuplotCommand() { closeGnuplotProcess(); }

GnuplotCommand::GnuplotCommand(GnuplotCommand&& other) noexcept
    : gnuplotPath_(std::move(other.gnuplotPath_)),
      fp_(std::move(other.fp_)),
      mutex_(),
      commandQueue_(std::move(other.commandQueue_)) {}

GnuplotCommand& GnuplotCommand::operator=(GnuplotCommand&& other) noexcept {
    if (this != &other) {
        gnuplotPath_ = std::move(other.gnuplotPath_);
        fp_ = std::move(other.fp_);
        commandQueue_ = std::move(other.commandQueue_);
    }
    return *this;
}

template <typename... Args>
void GnuplotCommand::command(const char* format, Args&&... args) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!fp_) {
        throw std::runtime_error("Gnuplot process is not initialized.");
    }

    std::string cmd = formatString(format, std::forward<Args>(args)...);
    commandQueue_.push(std::move(cmd));
}

void GnuplotCommand::flush() {
    std::lock_guard<std::mutex> lock(mutex_);
    while (!commandQueue_.empty()) {
        const auto& cmd = commandQueue_.front();
        if (fprintf(fp_.get(), "%s\n", cmd.c_str()) < 0) {
            throw std::runtime_error("Failed to write to gnuplot.");
        }
        commandQueue_.pop();
    }

    if (fflush(fp_.get()) != 0) {
        throw std::runtime_error("Failed to flush the file pointer.");
    }
}

void GnuplotCommand::openGnuplotProcess() {
    fp_.reset(popen(gnuplotPath_.c_str(), "w"));
    if (!fp_) {
        throw std::runtime_error("Failed to start gnuplot process.");
    }
}

void GnuplotCommand::closeGnuplotProcess() {
    if (fp_) {
        fp_.reset();
    }
}

// // Explicit instantiation of the template function
// template void GnuplotCommand::command<>(const char* format);
std::string GnuplotCommand::formatString(const char* format, Args&&... args) {
    int size = snprintf(nullptr, 0, format, std::forward<Args>(args)...);
    if (size < 0) {
        throw std::runtime_error("Error during formatting the string.");
    }
    std::vector<char> buf(size + 1);
    snprintf(buf.data(), buf.size(), format, std::forward<Args>(args)...);
    return std::string(buf.data(), buf.data() + size);
}
