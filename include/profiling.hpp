#pragma once

/* Throwaway memory profiler for the construction-profiling branch.
 *
 * Samples process RSS / allocator-live-heap and every registered named gauge
 * every `interval` ms, appending long-format rows to <dir>/samples.tsv:
 *
 *   t_ms  phase  rss  hwm  heap_live  heap_peak  gauge  bytes
 *
 * A named gauge is just an atomic byte counter that a structure updates when
 * it grows/shrinks (cheap, race-free).  The sampler also always writes one
 * `__rss__` row per tick so there is a timeline even with no gauges.
 *
 * On stop() it writes <dir>/summary.txt with the max of every gauge and when
 * it occurred, plus overall RSS peaks and the -g budget verdict.
 *
 * Linux is the supported platform (VmRSS + mallinfo2).  macOS is implemented
 * too (task_info + getrusage) so the output is not silent when developing
 * locally; other platforms return 0 for the live values.
 *
 * Controls (env):
 *   FULGOR_PROF=0        disable completely
 *   FULGOR_PROF_MS=500   sampling interval (min 50)
 *   FULGOR_PROF_DIR=...  output directory (default: <out dir>/prof/<stem>-<epoch_ms>)
 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <shared_mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <sys/resource.h>

#if defined(__GLIBC__)
#include <malloc.h>
#endif
#if defined(__linux__)
#include <fcntl.h>
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <mach/mach.h>
#endif

namespace fulgor::prof {

inline uint64_t peak_rss_bytes() {
    struct rusage ru;
    if (::getrusage(RUSAGE_SELF, &ru) != 0) return 0;
#if defined(__APPLE__)
    return static_cast<uint64_t>(ru.ru_maxrss);
#else
    return static_cast<uint64_t>(ru.ru_maxrss) * 1024ULL;
#endif
}

inline uint64_t rss_bytes() {
#if defined(__linux__)
    int fd = ::open("/proc/self/status", O_RDONLY | O_CLOEXEC);
    if (fd < 0) return 0;
    char buf[4096];
    ssize_t total = 0;
    for (;;) {
        ssize_t n = ::read(fd, buf + total, sizeof(buf) - 1 - static_cast<size_t>(total));
        if (n <= 0) break;
        total += n;
        if (static_cast<size_t>(total) >= sizeof(buf) - 1) break;
    }
    ::close(fd);
    if (total <= 0) return 0;
    buf[total] = 0;
    char const* p = std::strstr(buf, "VmRSS:");
    if (!p) return 0;
    p += 6;
    while (*p == ' ' or *p == '\t') ++p;
    char* end = nullptr;
    unsigned long kb = std::strtoul(p, &end, 10);
    if (end == p) return 0;
    return static_cast<uint64_t>(kb) * 1024ULL;
#elif defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info),
                  &count) != KERN_SUCCESS) {
        return 0;
    }
    return static_cast<uint64_t>(info.resident_size);
#else
    return 0;
#endif
}

inline uint64_t heap_live_bytes() {
#if defined(__GLIBC__) && defined(__GLIBC_PREREQ)
#if __GLIBC_PREREQ(2, 33)
    struct mallinfo2 mi = ::mallinfo2();
    return static_cast<uint64_t>(mi.uordblks);
#else
    struct mallinfo mi = ::mallinfo();
    return static_cast<uint64_t>(mi.uordblks);
#endif
#elif defined(__APPLE__)
    // malloc_zone_statistics would need a specific zone; not worth it for a
    // branch whose runs happen on Linux.
    return 0;
#else
    return 0;
#endif
}

struct gauge;

inline void register_gauge(gauge*);
inline void unregister_gauge(gauge*);
inline bool enabled();
inline int64_t now_ms();
/* A named atomic byte counter.  Structures call set()/add(); the sampler only
   reads it, so gauge updates never block.  Each gauge also keeps its own peak
   (and the time it occurred) so even a run shorter than the sampling interval
   still shows up in the summary when the gauge is destroyed. */
struct gauge {
    gauge() = default;
    explicit gauge(std::string name) : m_name(std::move(name)) { register_gauge(this); }
    ~gauge() { unregister_gauge(this); }

    gauge(gauge const&) = delete;
    gauge& operator=(gauge const&) = delete;
    gauge(gauge&&) = delete;
    gauge& operator=(gauge&&) = delete;

    void set(uint64_t v) {
        m_value.store(v, std::memory_order_relaxed);
        note_peak(v);
    }
    void add(int64_t d) {
        const uint64_t v =
            m_value.fetch_add(static_cast<uint64_t>(d), std::memory_order_relaxed) +
            static_cast<uint64_t>(d);
        if (d > 0) note_peak(v);
    }
    uint64_t get() const { return m_value.load(std::memory_order_relaxed); }
    uint64_t peak() const { return m_peak.load(std::memory_order_relaxed); }
    int64_t peak_at_ms() const { return m_peak_ms.load(std::memory_order_relaxed); }
    std::string const& name() const { return m_name; }

private:
    void note_peak(uint64_t v) {
        uint64_t prev = m_peak.load(std::memory_order_relaxed);
        if (v <= prev) return;
        const int64_t t = enabled() ? now_ms() : 0;
        while (v > prev && !m_peak.compare_exchange_weak(prev, v, std::memory_order_relaxed)) {
        }
        if (v == m_peak.load(std::memory_order_relaxed)) {
            m_peak_ms.store(t, std::memory_order_relaxed);
        }
    }

    friend void unregister_gauge(gauge*);

    std::string m_name;
    std::atomic<uint64_t> m_value{0};
    std::atomic<uint64_t> m_peak{0};
    std::atomic<int64_t> m_peak_ms{0};
};

struct gauge_record {
    std::string name;
    uint64_t max_bytes = 0;
    int64_t max_at_ms = 0;
};

struct sampler_state {
    std::shared_mutex reg_mu;              // guards `gauges`
    std::vector<gauge*> gauges;            // guarded by reg_mu
    std::mutex mu;                         // guards phase + records + streams
    std::string phase{"(none)"};
    std::vector<gauge_record> records;
    std::ofstream samples;
    std::ofstream phases;
    std::filesystem::path dir;
    std::atomic<bool> running{false};
    std::thread thread;
    std::chrono::steady_clock::time_point t0;
    uint32_t interval_ms = 500;
    uint64_t budget_bytes = 0;
    uint64_t rss_peak = 0;
    int64_t rss_peak_ms = 0;
    uint64_t heap_peak = 0;
    int64_t heap_peak_ms = 0;
    int64_t tick = 0;
};

inline sampler_state& state() {
    // Leaked on purpose: the sampler may outlive any static destructor order.
    static sampler_state* s = new sampler_state();
    return *s;
}

inline bool enabled() { return state().running.load(std::memory_order_relaxed); }

inline int64_t now_ms() {
    auto& s = state();
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() -
                                                                 s.t0)
        .count();
}

inline void update_record(std::string const& name, uint64_t bytes, int64_t t_ms) {
    auto& s = state();
    for (auto& r : s.records) {
        if (r.name == name) {
            if (bytes > r.max_bytes) {
                r.max_bytes = bytes;
                r.max_at_ms = t_ms;
            }
            return;
        }
    }
    s.records.push_back({name, bytes, t_ms});
}

inline void write_phase_row(char const* what, std::string const& name) {
    auto& s = state();
    s.phases << now_ms() << '\t' << what << '\t' << name << '\t' << rss_bytes() << '\t'
             << peak_rss_bytes() << '\t' << heap_live_bytes() << '\n';
    s.phases.flush();
}

inline void sampler_loop() {
    auto& s = state();
    while (s.running.load(std::memory_order_relaxed)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(s.interval_ms));
        if (!s.running.load(std::memory_order_relaxed)) break;

        const int64_t t_ms = now_ms();
        const uint64_t rss = rss_bytes();
        const uint64_t hwm = peak_rss_bytes();
        const uint64_t heap = heap_live_bytes();

        std::string phase;
        {
            std::lock_guard<std::mutex> lk(s.mu);
            phase = s.phase;
            if (rss > s.rss_peak) {
                s.rss_peak = rss;
                s.rss_peak_ms = t_ms;
            }
            if (heap > s.heap_peak) {
                s.heap_peak = heap;
                s.heap_peak_ms = t_ms;
            }
            ++s.tick;
        }

        {
            std::lock_guard<std::mutex> lk(s.mu);
            s.samples << t_ms << '\t' << phase << '\t' << rss << '\t' << hwm << '\t' << heap
                      << '\t' << s.heap_peak << '\t' << "__rss__" << '\t' << 0 << '\n';
        }

        std::shared_lock<std::shared_mutex> reg(s.reg_mu);
        for (gauge* g : s.gauges) {
            const uint64_t v = g->get();
            {
                std::lock_guard<std::mutex> lk(s.mu);
                update_record(g->name(), v, t_ms);
                s.samples << t_ms << '\t' << phase << '\t' << rss << '\t' << hwm << '\t' << heap
                          << '\t' << s.heap_peak << '\t' << g->name() << '\t' << v << '\n';
            }
        }
        {
            std::lock_guard<std::mutex> lk(s.mu);
            s.samples.flush();
        }
    }
}

inline void start(std::filesystem::path dir, uint32_t interval_ms = 500,
                  uint64_t budget_bytes = 0) {
    auto& s = state();
    if (s.running.exchange(true)) return;  // already running
    if (interval_ms < 50) interval_ms = 50;
    s.interval_ms = interval_ms;
    s.budget_bytes = budget_bytes;
    s.dir = std::move(dir);
    s.t0 = std::chrono::steady_clock::now();

    std::error_code ec;
    std::filesystem::create_directories(s.dir, ec);

    s.samples.open(s.dir / "samples.tsv", std::ios::out | std::ios::trunc);
    s.phases.open(s.dir / "phases.tsv", std::ios::out | std::ios::trunc);
    s.samples << "t_ms\tphase\trss\thwm\theap_live\theap_peak\tgauge\tbytes\n";
    s.phases << "t_ms\tevent\tphase\trss\thwm\theap_live\n";
    s.samples.flush();
    s.phases.flush();

    if (!s.samples.is_open()) {
        s.running.store(false, std::memory_order_relaxed);
        return;
    }

    std::cout << "[prof] sampling every " << s.interval_ms << " ms -> "
              << std::filesystem::absolute(s.dir).string() << std::endl;
    s.thread = std::thread(sampler_loop);
}

inline void stop() {
    auto& s = state();
    if (!s.running.exchange(false)) return;
    if (s.thread.joinable()) s.thread.join();

    {
        std::lock_guard<std::mutex> lk(s.mu);
        const uint64_t rss = rss_bytes();
        const uint64_t heap = heap_live_bytes();
        if (rss > s.rss_peak) {
            s.rss_peak = rss;
            s.rss_peak_ms = now_ms();
        }
        if (heap > s.heap_peak) {
            s.heap_peak = heap;
            s.heap_peak_ms = now_ms();
        }
    }

    std::ofstream out(s.dir / "summary.txt", std::ios::out | std::ios::trunc);
    if (!out.is_open()) return;

    out << "=== construction profiling summary ===\n";
    out << "samples: " << s.tick << "\n";
    out << "peak rss: " << s.rss_peak << " bytes at " << s.rss_peak_ms << " ms\n";
    out << "peak hwm: " << peak_rss_bytes() << " bytes\n";
    out << "peak heap_live: " << s.heap_peak << " bytes at " << s.heap_peak_ms << " ms\n";
    if (s.budget_bytes > 0) {
        out << "budget (-g): " << s.budget_bytes << " bytes\n";
        if (s.rss_peak > s.budget_bytes) {
            out << "VERDICT: RSS OVER budget by " << (s.rss_peak - s.budget_bytes) << " bytes\n";
        } else {
            out << "VERDICT: within budget\n";
        }
    }
    out << "\ngauge maxima:\n";
    out << "bytes\tat_ms\tname\n";
    std::sort(s.records.begin(), s.records.end(),
              [](gauge_record const& a, gauge_record const& b) { return a.max_bytes > b.max_bytes; });
    for (auto const& r : s.records) {
        out << r.max_bytes << '\t' << r.max_at_ms << '\t' << r.name << '\n';
    }
    if (s.budget_bytes > 0) {
        for (auto const& r : s.records) {
            if (r.max_bytes > s.budget_bytes) {
                out << "OVER -g: " << r.name << " reached " << r.max_bytes << " bytes (+"
                    << (r.max_bytes - s.budget_bytes) << ")\n";
            }
        }
    }
    out.flush();

    std::cout << "[prof] summary: " << std::filesystem::absolute(s.dir / "summary.txt").string()
              << std::endl;
}

inline void register_gauge(gauge* g) {
    if (g == nullptr) return;
    auto& s = state();
    std::unique_lock<std::shared_mutex> lk(s.reg_mu);
    s.gauges.push_back(g);
}

inline void unregister_gauge(gauge* g) {
    if (g == nullptr) return;
    auto& s = state();
    /* Same lock order as the sampler (reg -> mu) so there is no deadlock. */
    std::unique_lock<std::shared_mutex> lk(s.reg_mu);
    {
        std::lock_guard<std::mutex> mlk(s.mu);
        update_record(g->name(), g->peak(), g->peak_at_ms());
    }
    auto it = std::find(s.gauges.begin(), s.gauges.end(), g);
    if (it != s.gauges.end()) s.gauges.erase(it);
}

/* RAII phase label: nested guards restore the enclosing phase. */
struct phase_guard {
    explicit phase_guard(std::string name) {
        if (!enabled()) return;
        auto& s = state();
        std::lock_guard<std::mutex> lk(s.mu);
        m_prev = s.phase;
        s.phase = std::move(name);
        m_name = s.phase;
        m_active = true;
        write_phase_row("enter", m_name);
    }
    ~phase_guard() {
        if (!m_active) return;
        auto& s = state();
        std::lock_guard<std::mutex> lk(s.mu);
        write_phase_row("exit", m_name);
        s.phase = m_prev;
    }
    phase_guard(phase_guard const&) = delete;
    phase_guard& operator=(phase_guard const&) = delete;

private:
    bool m_active = false;
    std::string m_prev;
    std::string m_name;
};

/* Starts the sampler for the duration of a build.  Env-controlled. */
struct session {
    session(std::filesystem::path output_filename, uint64_t budget_bytes = 0) {
        if (char const* e = std::getenv("FULGOR_PROF"); e != nullptr and std::strcmp(e, "0") == 0) {
            return;
        }
        uint32_t ms = 500;
        if (char const* e = std::getenv("FULGOR_PROF_MS")) {
            long v = std::atol(e);
            if (v > 0) ms = static_cast<uint32_t>(v);
        }
        std::filesystem::path dir;
        if (char const* e = std::getenv("FULGOR_PROF_DIR")) {
            dir = e;
        } else {
            auto parent = output_filename.parent_path();
            if (parent.empty()) parent = ".";
            auto stem = output_filename.stem().string();
            if (stem.empty()) stem = "run";
            auto ts = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::system_clock::now().time_since_epoch())
                          .count();
            dir = parent / "prof" / (stem + "-" + std::to_string(ts));
        }
        start(std::move(dir), ms, budget_bytes);
        m_started = enabled();
    }
    ~session() {
        if (m_started) stop();
    }
    session(session const&) = delete;
    session& operator=(session const&) = delete;

private:
    bool m_started = false;
};

}  // namespace fulgor::prof
