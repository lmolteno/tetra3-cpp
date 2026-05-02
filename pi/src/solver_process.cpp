#include "solver_process.h"

#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

double json_float(const std::string &s, const char *key, double def = 0) {
    std::string needle = std::string("\"") + key + "\"";
    auto pos = s.find(needle);
    if (pos == std::string::npos) return def;
    pos = s.find(':', pos);
    if (pos == std::string::npos) return def;
    return std::strtod(s.c_str() + pos + 1, nullptr);
}

bool json_bool(const std::string &s, const char *key, bool def = false) {
    std::string needle = std::string("\"") + key + "\"";
    auto pos = s.find(needle);
    if (pos == std::string::npos) return def;
    pos = s.find(':', pos);
    if (pos == std::string::npos) return def;
    pos++;
    while (pos < s.size() && s[pos] == ' ') pos++;
    return s.compare(pos, 4, "true") == 0;
}

}

SolverProcess::SolverProcess(Config cfg) : cfg_(std::move(cfg)) {}

SolverProcess::~SolverProcess() {
    shutdown();
}

bool SolverProcess::start() {
    int in_pipe[2];   // parent writes [1], child reads [0]
    int out_pipe[2];  // child writes [1], parent reads [0]
    if (pipe(in_pipe) < 0 || pipe(out_pipe) < 0) {
        std::cerr << "SolverProcess: pipe failed: " << std::strerror(errno) << "\n";
        return false;
    }

    pid_t pid = fork();
    if (pid < 0) {
        std::cerr << "SolverProcess: fork failed: " << std::strerror(errno) << "\n";
        return false;
    }

    if (pid == 0) {
        // Child: wire pipes onto stdin/stdout, exec solve_cli.
        dup2(in_pipe[0], STDIN_FILENO);
        dup2(out_pipe[1], STDOUT_FILENO);
        close(in_pipe[0]);
        close(in_pipe[1]);
        close(out_pipe[0]);
        close(out_pipe[1]);

        std::vector<std::string> args = {cfg_.binary_path, "--batch"};
        if (!cfg_.db_stars.empty())    { args.push_back("--db-stars");    args.push_back(cfg_.db_stars); }
        if (!cfg_.db_patterns.empty()) { args.push_back("--db-patterns"); args.push_back(cfg_.db_patterns); }
        char fov_buf[32], sigma_buf[32], crop_buf[32], bg_buf[32];
        std::snprintf(fov_buf,   sizeof(fov_buf),   "%g", cfg_.fov_deg);
        std::snprintf(sigma_buf, sizeof(sigma_buf), "%g", cfg_.detection_sigma);
        std::snprintf(crop_buf,  sizeof(crop_buf),  "%d", cfg_.crop_size);
        std::snprintf(bg_buf,    sizeof(bg_buf),    "%d", cfg_.bg_tile);
        args.push_back("--fov");             args.push_back(fov_buf);
        args.push_back("--detection-sigma"); args.push_back(sigma_buf);
        args.push_back("--crop");            args.push_back(crop_buf);
        args.push_back("--bg-tile");         args.push_back(bg_buf);

        std::vector<char *> argv;
        argv.reserve(args.size() + 1);
        for (auto &a : args) argv.push_back(a.data());
        argv.push_back(nullptr);

        execvp(cfg_.binary_path.c_str(), argv.data());
        std::fprintf(stderr, "SolverProcess: execvp(%s) failed: %s\n",
                     cfg_.binary_path.c_str(), std::strerror(errno));
        std::_Exit(127);
    }

    // Parent
    close(in_pipe[0]);
    close(out_pipe[1]);
    pid_ = pid;
    stdin_fd_ = in_pipe[1];
    stdout_fd_ = out_pipe[0];
    running_ = true;
    reader_ = std::thread(&SolverProcess::reader_loop, this);
    return true;
}

bool SolverProcess::is_busy() const {
    std::lock_guard<std::mutex> lock(mu_);
    return in_flight_;
}

bool SolverProcess::submit(const std::string &path) {
    {
        std::lock_guard<std::mutex> lock(mu_);
        if (in_flight_) return false;
        in_flight_ = true;
        pending_path_ = path;
    }

    std::string line = path + "\n";
    const char *p = line.data();
    std::size_t left = line.size();
    while (left > 0) {
        ssize_t n = write(stdin_fd_, p, left);
        if (n < 0) {
            if (errno == EINTR) continue;
            std::cerr << "SolverProcess: write failed: " << std::strerror(errno) << "\n";
            std::lock_guard<std::mutex> lock(mu_);
            in_flight_ = false;
            return false;
        }
        p += n;
        left -= static_cast<std::size_t>(n);
    }
    return true;
}

std::optional<SolverProcess::ResultStamp> SolverProcess::poll_fresh() {
    std::lock_guard<std::mutex> lock(mu_);
    if (!latest_ || latest_->seq <= last_polled_seq_) return std::nullopt;
    last_polled_seq_ = latest_->seq;
    return latest_;
}

void SolverProcess::reader_loop() {
    std::string buf;
    char chunk[4096];
    while (running_) {
        ssize_t n = read(stdout_fd_, chunk, sizeof(chunk));
        if (n < 0) {
            if (errno == EINTR) continue;
            break;
        }
        if (n == 0) break;  // EOF
        buf.append(chunk, static_cast<std::size_t>(n));

        // Process complete lines
        std::size_t start = 0;
        for (std::size_t i = 0; i < buf.size(); i++) {
            if (buf[i] != '\n') continue;
            std::string line = buf.substr(start, i - start);
            start = i + 1;
            if (line.empty()) continue;

            SolveResult r = parse_result(line);
            std::lock_guard<std::mutex> lock(mu_);
            latest_ = ResultStamp{r, pending_path_, next_seq_++};
            in_flight_ = false;
            pending_path_.clear();
        }
        if (start > 0) buf.erase(0, start);
    }
    running_ = false;
}

SolveResult SolverProcess::parse_result(const std::string &json) {
    SolveResult r{};
    r.solved = json_bool(json, "solved");
    if (r.solved) {
        r.ra            = static_cast<float>(json_float(json, "ra_rad"));
        r.dec           = static_cast<float>(json_float(json, "dec_rad"));
        r.roll          = static_cast<float>(json_float(json, "roll_rad"));
        r.fov           = static_cast<float>(json_float(json, "fov_rad"));
        r.rmse          = static_cast<float>(json_float(json, "rmse"));
        r.num_matches   = static_cast<int>(json_float(json, "num_matches"));
        r.solve_time_ms = static_cast<float>(json_float(json, "solve_time_ms"));
    } else {
        r.solve_time_ms = static_cast<float>(json_float(json, "solve_time_ms"));
    }
    return r;
}

void SolverProcess::shutdown() {
    if (!running_ && pid_ < 0) return;
    running_ = false;
    if (stdin_fd_ >= 0) {
        close(stdin_fd_);
        stdin_fd_ = -1;
    }
    if (reader_.joinable()) reader_.join();
    if (stdout_fd_ >= 0) {
        close(stdout_fd_);
        stdout_fd_ = -1;
    }
    if (pid_ > 0) {
        // Give the child a moment, then SIGTERM if still alive.
        for (int i = 0; i < 10; i++) {
            int status;
            pid_t r = waitpid(pid_, &status, WNOHANG);
            if (r == pid_) { pid_ = -1; return; }
            if (r < 0) break;
            usleep(20000);
        }
        kill(pid_, SIGTERM);
        int status;
        waitpid(pid_, &status, 0);
        pid_ = -1;
    }
}
