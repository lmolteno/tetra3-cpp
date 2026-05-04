#pragma once

#include <tetra3/types.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <sys/types.h>

// Owns a `solve_cli --batch` subprocess. The pi_tracker daemon writes frame
// paths to it and consumes JSON results on a reader thread. Only one solve
// runs at a time — submit() returns false (drops the frame) while the solver
// is still working on the previous frame.
class SolverProcess {
public:
    struct Config {
        std::string binary_path = "solve_cli";
        std::string db_stars;          // empty -> let solve_cli use its default
        std::string db_patterns;
        float fov_deg = 11.0f;
        float detection_sigma = 5.0f;
        int crop_size = 720;
        int bg_tile = 128;
        // Single-coefficient lens distortion correction. Calibrated on the
        // 260504_pi_pictures dataset (Camera Module 3 NoIR Standard, 4.74mm
        // lens). Set to NaN to disable. Default cuts Dec scatter ~3x and
        // RMSE ~5x on full-frame solves; small effect on the crop pass.
        float distortion_k = 0.024f;
    };

    struct ResultStamp {
        SolveResult result{};
        std::string frame_path;                     // which frame produced it
        std::uint64_t seq = 0;
        // Stages reported by solve_cli (-1 = absent in this result).
        float load_ms = -1.0f;
        float detect_ms = -1.0f;
        float bg_ms = -1.0f;
        float mask_ms = -1.0f;
        float cluster_ms = -1.0f;
        // Wall-clock interval between submit and result-arrived in pi_tracker.
        // Includes IPC + everything solve_cli did.
        float wall_ms = 0.0f;
        // When the result was received (steady_clock). Lets the caller compute
        // how stale this result is.
        std::chrono::steady_clock::time_point received_at{};
    };

    explicit SolverProcess(Config cfg);
    ~SolverProcess();

    SolverProcess(const SolverProcess &) = delete;
    SolverProcess &operator=(const SolverProcess &) = delete;

    // Spawn the child. Returns false if exec fails.
    bool start();

    // Push a frame path to be solved. Returns false if a previous solve is
    // still in flight (the frame is dropped).
    bool submit(const std::string &path);

    bool is_busy() const;

    // Returns the latest result if it hasn't been polled before, otherwise
    // std::nullopt. The same result is never returned twice.
    std::optional<ResultStamp> poll_fresh();

private:
    Config cfg_;
    pid_t pid_ = -1;
    int stdin_fd_ = -1;
    int stdout_fd_ = -1;

    std::thread reader_;
    std::atomic<bool> running_{false};

    mutable std::mutex mu_;
    bool in_flight_ = false;
    std::string pending_path_;
    std::chrono::steady_clock::time_point pending_submitted_at_{};
    std::optional<ResultStamp> latest_;
    std::uint64_t last_polled_seq_ = 0;
    std::uint64_t next_seq_ = 1;

    void reader_loop();
    void shutdown();

    struct ParsedResult {
        SolveResult result{};
        float load_ms = -1.0f;
        float detect_ms = -1.0f;
        float bg_ms = -1.0f;
        float mask_ms = -1.0f;
        float cluster_ms = -1.0f;
    };
    static ParsedResult parse_result(const std::string &json);
};
