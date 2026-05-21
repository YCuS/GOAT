#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct Pwm {
    std::string id;
    std::vector<std::array<double, 4>> columns; // A, C, G, T
};

struct Threshold {
    std::string motif_id;
    double full_score = 0.0;
    double core_score = 0.0;
    int core_start = 0; // 0-based
    int motif_length = 0;
};

static std::string trim(const std::string &text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

static int base_index(char base) {
    switch (base) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't':
    case 'U': case 'u': return 3;
    default: return -1;
    }
}

static std::array<double, 4> normalized_row(double a, double c, double g, double t) {
    std::array<double, 4> row{a, c, g, t};
    double sum = std::accumulate(row.begin(), row.end(), 0.0);
    if (sum <= 0.0) return {0.25, 0.25, 0.25, 0.25};
    for (double &value : row) value /= sum;
    return row;
}

static std::vector<Pwm> load_pwms(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open PWM file: " + path);

    std::vector<Pwm> motifs;
    Pwm current;
    std::string line;

    auto flush_current = [&]() {
        if (!current.id.empty() && !current.columns.empty()) {
            motifs.push_back(std::move(current));
        }
        current = Pwm{};
    };

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        if (line.rfind("Motif:", 0) == 0) {
            flush_current();
            current.id = trim(line.substr(6));
            continue;
        }

        std::istringstream row_stream(line);
        double a = 0.0, c = 0.0, g = 0.0, t = 0.0;
        if (row_stream >> a >> c >> g >> t) {
            current.columns.push_back(normalized_row(a, c, g, t));
        }
    }
    flush_current();
    return motifs;
}

static double column_entropy(const std::array<double, 4> &column) {
    double entropy = 0.0;
    for (double p : column) {
        if (p > 0.0) entropy -= p * std::log(p);
    }
    return entropy;
}

static int choose_core_start(const Pwm &pwm, int core_len) {
    if (core_len <= 0 || core_len > static_cast<int>(pwm.columns.size())) return 0;

    double best_entropy = std::numeric_limits<double>::infinity();
    int best_start = 0;
    for (int start = 0; start + core_len <= static_cast<int>(pwm.columns.size()); ++start) {
        double total_entropy = 0.0;
        for (int i = 0; i < core_len; ++i) {
            total_entropy += column_entropy(pwm.columns[start + i]);
        }
        if (total_entropy < best_entropy) {
            best_entropy = total_entropy;
            best_start = start;
        }
    }
    return best_start;
}

static double score_region(const std::string &sequence,
                           int sequence_start,
                           const Pwm &pwm,
                           int pwm_start,
                           int length) {
    static constexpr double eps = 1e-12;
    double score = 0.0;

    for (int i = 0; i < length; ++i) {
        const int base = base_index(sequence[sequence_start + i]);
        const double probability = (base >= 0) ? pwm.columns[pwm_start + i][base] : 0.25;
        score += -std::log(std::max(probability, eps));
    }
    return score;
}

static uint32_t stable_seed(const std::string &text) {
    uint32_t hash = 2166136261u;
    for (unsigned char ch : text) {
        hash ^= ch;
        hash *= 16777619u;
    }
    return hash;
}

static std::string sample_sequence_from_pwm(const Pwm &pwm, std::mt19937 &rng) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string sequence;
    sequence.reserve(pwm.columns.size());

    for (const auto &column : pwm.columns) {
        std::discrete_distribution<int> distribution(column.begin(), column.end());
        sequence.push_back(bases[distribution(rng)]);
    }
    return sequence;
}

static double percentile(std::vector<double> values, double pct) {
    if (values.empty()) return 0.0;
    pct = std::clamp(pct, 0.0, 100.0);
    std::sort(values.begin(), values.end());
    const double rank = (pct / 100.0) * static_cast<double>(values.size());
    size_t index = static_cast<size_t>(std::ceil(rank));
    if (index == 0) index = 1;
    index = std::min(index - 1, values.size() - 1);
    return values[index];
}

static Threshold calculate_threshold(const Pwm &pwm,
                                     double full_pct,
                                     double core_pct,
                                     int iterations,
                                     int core_len) {
    const int motif_len = static_cast<int>(pwm.columns.size());
    if (core_len <= 0 || core_len > motif_len) {
        throw std::runtime_error("Invalid core length for motif: " + pwm.id);
    }

    const int core_start = choose_core_start(pwm, core_len);
    std::vector<double> full_scores;
    std::vector<double> core_scores;
    full_scores.reserve(iterations);
    core_scores.reserve(iterations);

    std::mt19937 rng(stable_seed(pwm.id));
    for (int i = 0; i < iterations; ++i) {
        const std::string sequence = sample_sequence_from_pwm(pwm, rng);
        full_scores.push_back(score_region(sequence, 0, pwm, 0, motif_len));
        core_scores.push_back(score_region(sequence, core_start, pwm, core_start, core_len));
    }

    return {pwm.id,
            percentile(std::move(full_scores), full_pct),
            percentile(std::move(core_scores), core_pct),
            core_start,
            motif_len};
}

int main(int argc, char *argv[]) {
    if (argc != 8) {
        std::cerr
            << "Usage: " << argv[0]
            << " <pwm.txt> <thresholds.tsv> <full_percentile> <core_percentile>"
            << " <iterations> <core_length> <precision>\n";
        return 1;
    }

    try {
        const std::string pwm_path = argv[1];
        const std::string out_path = argv[2];
        const double full_pct = std::stod(argv[3]);
        const double core_pct = std::stod(argv[4]);
        const int iterations = std::stoi(argv[5]);
        const int core_len = std::stoi(argv[6]);
        const int precision = std::stoi(argv[7]);

        if (iterations <= 0) throw std::runtime_error("iterations must be > 0");

        const auto started = std::chrono::steady_clock::now();
        const auto motifs = load_pwms(pwm_path);
        if (motifs.empty()) throw std::runtime_error("No motifs found in PWM file");

        std::ofstream out(out_path);
        if (!out) throw std::runtime_error("Cannot create threshold file: " + out_path);

        out << "# motif_id\tfull_threshold\tcore_threshold\tcore_start\tmotif_length\n";
        int written = 0;
        int skipped = 0;

        for (const auto &pwm : motifs) {
            try {
                const Threshold threshold =
                    calculate_threshold(pwm, full_pct, core_pct, iterations, core_len);
                out << threshold.motif_id << '\t'
                    << std::fixed << std::setprecision(precision)
                    << threshold.full_score << '\t'
                    << threshold.core_score << '\t'
                    << threshold.core_start << '\t'
                    << threshold.motif_length << '\n';
                ++written;
            } catch (const std::exception &error) {
                ++skipped;
                std::cerr << "Skipping " << pwm.id << ": " << error.what() << '\n';
            }
        }

        const auto finished = std::chrono::steady_clock::now();
        const auto ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(finished - started).count();

        std::cerr << "Thresholds written: " << written << " motifs";
        if (skipped > 0) std::cerr << " (" << skipped << " skipped)";
        std::cerr << "\nOutput: " << out_path << "\nElapsed: " << ms << " ms\n";
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
