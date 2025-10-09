#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <thread>
#include <future>
#include <mutex>
#include <regex>
#include <chrono>
#include <limits> // 新增：choose_core_start_by_min_sum_entropy 依赖

// Position Weight Matrix structure for threshold calculation
struct PositionWeightMatrix
{
    std::map<char, std::vector<double>> nucleotide_probabilities;
    int get_length() const { return nucleotide_probabilities.at('A').size(); }
};

// 新增：计算motif的GC含量（按列归一化后取 (G+C)/sum，再对列求平均）
static double compute_motif_gc(const PositionWeightMatrix &pwm)
{
    const int L = pwm.get_length();
    if (L <= 0)
        return 0.0;
    double gc_sum = 0.0;
    for (int i = 0; i < L; ++i)
    {
        double a = pwm.nucleotide_probabilities.at('A')[i];
        double c = pwm.nucleotide_probabilities.at('C')[i];
        double g = pwm.nucleotide_probabilities.at('G')[i];
        double t = pwm.nucleotide_probabilities.at('T')[i];
        double s = a + c + g + t;
        if (s <= 0.0)
        {
            gc_sum += 0.5; // 退化列，按0.5处理以避免偏置
        }
        else
        {
            gc_sum += (g + c) / s;
        }
    }
    return gc_sum / static_cast<double>(L);
}

// Results container for calculated thresholds
struct MotifThresholds
{
    std::string motif_identifier;
    double full_motif_threshold;
    double core_region_threshold;
    int core_start_position;
    bool calculation_successful = false;
    // 新增字段
    double gc_content = 0.0;
    double full_motif_threshold_unadjusted = 0.0; // unadjusted full threshold (for reference)
    // strict thresholds (used when consecutive hits are detected in search)
    double strict_full_threshold = 0.0;
    double strict_core_threshold = 0.0;
    // relaxed thresholds (explicit values for reference/debugging)
    double relaxed_full_threshold = 0.0;
    double relaxed_core_threshold = 0.0;
    // cluster index (1-based); 0 means not clustered
    int cluster_index = 0;
};

// Parse PWM file and extract probability matrices with memory optimization
std::map<std::string, std::vector<std::vector<double>>>
load_probability_matrices(const std::string &pwm_file)
{
    std::ifstream file(pwm_file);
    if (!file)
        throw std::runtime_error("Failed to open PWM file: " + pwm_file);

    std::map<std::string, std::vector<std::vector<double>>> motif_matrices;
    std::string line, current_motif_id;
    std::vector<std::vector<double>> current_matrix;
    current_matrix.reserve(30); // Reserve for typical motif length

    // Pattern to match probability matrix lines (4 decimal numbers)
    std::regex probability_line_pattern(R"(^\s*[\d.eE+-]+(\s+[\d.eE+-]+){3}\s*$)");

    while (std::getline(file, line))
    {
        // Remove carriage returns
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        // Detect motif header
        if (line.compare(0, 6, "Motif:") == 0)
        {
            // Store previous motif if complete
            if (!current_motif_id.empty() && !current_matrix.empty())
            {
                motif_matrices[current_motif_id] = std::move(current_matrix);
            }

            // Initialize new motif
            current_motif_id = line.substr(6);
            current_matrix.clear();
            current_matrix.reserve(30);
        }
        // Parse probability values
        else if (std::regex_match(line, probability_line_pattern))
        {
            std::istringstream stream(line);
            std::vector<double> position_probabilities;
            position_probabilities.reserve(4);

            double probability;
            while (stream >> probability)
            {
                position_probabilities.push_back(probability);
            }

            if (position_probabilities.size() == 4)
            {
                current_matrix.push_back(std::move(position_probabilities));
            }
        }
    }

    // Store final motif
    if (!current_motif_id.empty() && !current_matrix.empty())
    {
        motif_matrices[current_motif_id] = std::move(current_matrix);
    }

    return motif_matrices;
}

// Convert matrix format to PWM structure for calculations
PositionWeightMatrix create_pwm_from_matrix(const std::vector<std::vector<double>> &probability_matrix)
{
    PositionWeightMatrix pwm;

    // Pre-allocate vectors for nucleotides
    for (char nucleotide : {'A', 'C', 'G', 'T'})
    {
        pwm.nucleotide_probabilities[nucleotide].reserve(probability_matrix.size());
    }

    // Fill PWM structure: rows become columns
    for (const auto &position_probs : probability_matrix)
    {
        pwm.nucleotide_probabilities['A'].push_back(position_probs[0]);
        pwm.nucleotide_probabilities['C'].push_back(position_probs[1]);
        pwm.nucleotide_probabilities['G'].push_back(position_probs[2]);
        pwm.nucleotide_probabilities['T'].push_back(position_probs[3]);
    }
    return pwm;
}

// Generate random sequence according to PWM probabilities
std::string generate_sequence_from_pwm(const PositionWeightMatrix &pwm,
                                       int sequence_length,
                                       std::mt19937 &random_generator)
{
    const std::string nucleotides = "ACGT";
    std::string generated_sequence;
    generated_sequence.reserve(sequence_length);

    for (int position = 0; position < sequence_length; ++position)
    {
        std::vector<double> position_probabilities;
        position_probabilities.reserve(4);
        double probability_sum = 0;

        // Collect probabilities for this position
        for (char nucleotide : nucleotides)
        {
            position_probabilities.push_back(pwm.nucleotide_probabilities.at(nucleotide)[position]);
            probability_sum += position_probabilities.back();
        }

        // Generate nucleotide based on probabilities
        if (probability_sum > 0)
        {
            // Normalize probabilities
            for (auto &prob : position_probabilities)
                prob /= probability_sum;

            std::discrete_distribution<> nucleotide_distribution(
                position_probabilities.begin(), position_probabilities.end());
            generated_sequence += nucleotides[nucleotide_distribution(random_generator)];
        }
        else
        {
            // Fallback for zero probabilities
            generated_sequence += nucleotides[random_generator() % 4];
        }
    }
    return generated_sequence;
}

// Calculate log-likelihood score for sequence against PWM
double calculate_log_likelihood_score(const std::string &sequence, const PositionWeightMatrix &pwm)
{
    double total_score = 0;
    constexpr double zero_probability_penalty = 10.0; // Penalty for zero probabilities

    for (size_t position = 0; position < sequence.size(); ++position)
    {
        double nucleotide_probability = pwm.nucleotide_probabilities.at(sequence[position])[position];
        total_score += (nucleotide_probability > 0) ? -std::log(nucleotide_probability) : zero_probability_penalty;
    }
    return total_score;
}

// Calculate per-column entropy H = -sum_b p_b log p_b (nats) with small epsilon for stability
static inline double column_entropy_at(const PositionWeightMatrix &pwm, int pos)
{
    static constexpr double eps = 1e-12;
    double H = 0.0;
    for (char nuc : {'A', 'C', 'G', 'T'})
    {
        double p = pwm.nucleotide_probabilities.at(nuc)[pos];
        if (p <= 0.0)
            p = eps;
        H += -p * std::log(p);
    }
    return H;
}

// Choose core start index by minimizing the sum of entropies over a window of length core_len
static int choose_core_start_by_min_sum_entropy(const PositionWeightMatrix &pwm, int core_len)
{
    const int L = pwm.get_length();
    if (core_len <= 0 || core_len > L)
        return 0;
    double best_sum = std::numeric_limits<double>::infinity();
    int best_start = 0;
    for (int s = 0; s + core_len <= L; ++s)
    {
        double sumH = 0.0;
        for (int j = s; j < s + core_len; ++j)
        {
            sumH += column_entropy_at(pwm, j);
        }
        if (sumH < best_sum)
        {
            best_sum = sumH;
            best_start = s;
        }
    }
    return best_start;
}

// Calculate thresholds for a single motif using parallel Monte Carlo simulation
MotifThresholds calculate_motif_thresholds_parallel(const std::string &motif_id,
                                                    const PositionWeightMatrix &full_pwm,
                                                    double full_threshold_percentile,
                                                    double core_threshold_percentile,
                                                    int simulation_iterations,
                                                    int core_region_length,
                                                    // 保留 relaxed/strict 百分位，但取消 GC 插值逻辑
                                                    double /*gc_threshold_unused*/,
                                                    double relaxed_full_pct,
                                                    double relaxed_core_pct,
                                                    double strict_full_pct,
                                                    double strict_core_pct)
{
    MotifThresholds result{motif_id, 0, 0, 0, false};
    int motif_length = full_pwm.get_length();
    if (motif_length < core_region_length)
        return result;

    // 新增：计算GC含量
    double gc = compute_motif_gc(full_pwm);
    result.gc_content = gc;

    // Use thread-local random generator for thread safety
    thread_local std::mt19937 random_generator(std::random_device{}());

    std::vector<double> full_motif_scores, core_region_scores;
    full_motif_scores.reserve(simulation_iterations);
    core_region_scores.reserve(simulation_iterations);

    // Select core as the contiguous window with minimal sum of (-log p), i.e., minimal entropy sum
    int core_start_position = choose_core_start_by_min_sum_entropy(full_pwm, core_region_length);

    // Build core PWM from the selected window
    PositionWeightMatrix core_pwm;
    for (char nucleotide : {'A', 'C', 'G', 'T'})
    {
        core_pwm.nucleotide_probabilities[nucleotide].reserve(core_region_length);
        for (int pos = core_start_position;
             pos < core_start_position + core_region_length && pos < motif_length; ++pos)
        {
            core_pwm.nucleotide_probabilities[nucleotide].push_back(
                full_pwm.nucleotide_probabilities.at(nucleotide)[pos]);
        }
    }

    // Monte Carlo simulation
    for (int iteration = 0; iteration < simulation_iterations; ++iteration)
    {
        // Generate sequence from full PWM
        auto simulated_sequence = generate_sequence_from_pwm(full_pwm, motif_length, random_generator);
        full_motif_scores.push_back(calculate_log_likelihood_score(simulated_sequence, full_pwm));

        // Score the selected core region (fixed start determined by PWM)
        if (simulated_sequence.size() >= core_region_length)
        {
            auto core_sequence = simulated_sequence.substr(core_start_position, core_region_length);
            core_region_scores.push_back(calculate_log_likelihood_score(core_sequence, core_pwm));
        }
    }

    // Calculate threshold values from percentiles
    std::sort(full_motif_scores.begin(), full_motif_scores.end());
    std::sort(core_region_scores.begin(), core_region_scores.end());

    // 取消 GC 插值：直接使用原始百分位作为常规阈值（cluster != 3 使用）
    auto pct_full_norm = full_threshold_percentile;
    auto pct_core_norm = core_threshold_percentile;

    auto idx_from_percent = [](size_t n, double pct)
    {
        if (n == 0)
            return 0;
        int idx = static_cast<int>(std::floor(n * pct / 100.0)) - 1;
        if (idx < 0)
            idx = 0;
        if (idx >= static_cast<int>(n))
            idx = static_cast<int>(n) - 1;
        return idx;
    };

    // indices for normal thresholds（不再做 GC 调整）
    int full_threshold_index_eff = idx_from_percent(full_motif_scores.size(), pct_full_norm);
    int full_threshold_index_orig = idx_from_percent(full_motif_scores.size(), full_threshold_percentile);
    int core_threshold_index = idx_from_percent(core_region_scores.size(), pct_core_norm);

    // indices for strict thresholds (always computed regardless of GC)
    int strict_full_index = idx_from_percent(full_motif_scores.size(), strict_full_pct);
    int strict_core_index = idx_from_percent(core_region_scores.size(), strict_core_pct);

    // indices for explicit relaxed thresholds (independent of GC condition)
    int relaxed_full_index = idx_from_percent(full_motif_scores.size(), relaxed_full_pct);
    int relaxed_core_index = idx_from_percent(core_region_scores.size(), relaxed_core_pct);

    // full_motif_threshold 与未调整阈值保持一致（原始百分位）
    double thr_orig = full_motif_scores[full_threshold_index_orig];
    result.full_motif_threshold_unadjusted = thr_orig;
    result.full_motif_threshold = thr_orig;
    result.core_region_threshold = core_region_scores[core_threshold_index];
    // strict thresholds to be used in search for consecutive hits
    result.strict_full_threshold = full_motif_scores[strict_full_index];
    result.strict_core_threshold = core_region_scores[strict_core_index];
    // relaxed thresholds (explicit)
    result.relaxed_full_threshold = full_motif_scores[relaxed_full_index];
    result.relaxed_core_threshold = core_region_scores[relaxed_core_index];
    result.core_start_position = core_start_position;
    result.calculation_successful = true;
    return result;
}

int main(int argc, char *argv[])
{
    // keep backward-compat positional args, plus optional flags for GC/relaxed/strict
    if (argc < 8)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <input_pwm_file> <output_threshold_file> <full_threshold_percentile> "
                  << "<core_threshold_percentile> <simulation_iterations> <core_length> <output_precision>\n";
        std::cerr << "Optional flags:\n"
                  << "  --gc-thresh=X           (default: 0.30)\n"
                  << "  --relaxed-full=PCT      (default: equal to full pct)\n"
                  << "  --relaxed-core=PCT      (default: equal to core pct)\n"
                  << "  --strict-full=PCT       (default: equal to full pct)\n"
                  << "  --strict-core=PCT       (default: equal to core pct)\n";
        return 1;
    }

    try
    {
        std::string input_pwm_file = argv[1];
        std::string output_threshold_file = argv[2];
        double full_threshold_percentile = std::stod(argv[3]);
        double core_threshold_percentile = std::stod(argv[4]);
        int simulation_iterations = std::stoi(argv[5]);
        int core_region_length = std::stoi(argv[6]);
        int output_precision = std::stoi(argv[7]);

        // defaults for optional flags
        double gc_threshold = 0.30;
        double relaxed_full_pct = full_threshold_percentile;
        double relaxed_core_pct = core_threshold_percentile;
        double strict_full_pct = full_threshold_percentile;
        double strict_core_pct = core_threshold_percentile;
        // clustering related flags
        bool clustering_enabled = false;
        int cluster_k = 4;
        int cluster_relaxed_index = 3; // which cluster uses relaxed thresholds (1-based)
        std::string cluster_plot_path; // optional SVG output

        // parse optional flags
        for (int i = 8; i < argc; ++i)
        {
            std::string a = argv[i];
            auto kv = [&](const std::string &key, double &dst)
            {
                if (a.rfind(key + "=", 0) == 0) { dst = std::stod(a.substr(key.size() + 1)); return true; }
                return false; };
            if (kv("--gc-thresh", gc_threshold))
            {
            }
            else if (kv("--relaxed-full", relaxed_full_pct))
            {
            }
            else if (kv("--relaxed-core", relaxed_core_pct))
            {
            }
            else if (kv("--strict-full", strict_full_pct))
            {
            }
            else if (kv("--strict-core", strict_core_pct))
            {
            }
            else if (a.rfind("--cluster-enable=", 0) == 0)
            {
                std::string v = a.substr(17);
                clustering_enabled = (v == "1" || v == "true" || v == "TRUE");
            }
            else if (a.rfind("--cluster-k=", 0) == 0)
            {
                cluster_k = std::max(1, std::stoi(a.substr(12)));
            }
            else if (a.rfind("--cluster-relaxed=", 0) == 0)
            {
                cluster_relaxed_index = std::max(1, std::stoi(a.substr(18)));
            }
            else if (a.rfind("--cluster-plot=", 0) == 0)
            {
                cluster_plot_path = a.substr(15);
            }
            else
            {
                std::cerr << "Warning: unknown flag ignored: " << a << "\n";
            }
        }

        // std::cerr << "Loading PWM probability matrices..." << std::endl;
        auto motif_matrices = load_probability_matrices(input_pwm_file);
        std::cerr << "Found " << motif_matrices.size() << " motifs for processing" << std::endl;

        // std::cerr << "Starting parallel threshold calculations..." << std::endl;
        auto calculation_start_time = std::chrono::high_resolution_clock::now();

        // Pre-compute GC mean/variance features for clustering
        std::vector<std::string> motif_ids;
        motif_ids.reserve(motif_matrices.size());
        std::vector<double> gc_means, gc_vars;
        gc_means.reserve(motif_matrices.size());
        gc_vars.reserve(motif_matrices.size());
        for (const auto &kv : motif_matrices)
        {
            PositionWeightMatrix pwm = create_pwm_from_matrix(kv.second);
            int L = pwm.get_length();
            if (L <= 0)
                continue;
            double sum_gc = 0.0;
            double sum_var_gc = 0.0; // 按列 p_gc * (1 - p_gc) 求和 (不再用与均值差平方的统计方差)
            int valid_cols = 0;
            for (int i = 0; i < L; ++i)
            {
                double a = pwm.nucleotide_probabilities['A'][i];
                double c = pwm.nucleotide_probabilities['C'][i];
                double g = pwm.nucleotide_probabilities['G'][i];
                double t = pwm.nucleotide_probabilities['T'][i];
                double s = a + c + g + t;
                if (s <= 0)
                    continue;
                double pgc = (c + g) / s;
                sum_gc += pgc;
                sum_var_gc += pgc * (1.0 - pgc); // Bernoulli 方差 p(1-p) 按位累加
                valid_cols++;
            }
            if (valid_cols == 0)
                continue;
            double mean_gc = sum_gc / valid_cols; // motif GC 均值
            // sum_var_gc 已是所需“方差”度量（按列求和），直接用于后续 min-max 归一化
            motif_ids.push_back(kv.first);
            gc_means.push_back(mean_gc);
            gc_vars.push_back(sum_var_gc);
        }
        std::map<std::string, int> cluster_assign; // motif_id -> cluster
        // norm_feats: (raw_gc_mean, normalized_variance) 用于聚类（gc 不归一化，与 notebook 一致）
        std::vector<std::pair<double, double>> norm_feats;
        // 记录绘图需要的范围（gc 原始范围 + 方差归一化后 0-1）
        double min_m = 0.0, max_m = 1.0; // gc mean 原始范围用于 SVG 映射
        if (clustering_enabled && !motif_ids.empty())
        {
            min_m = *std::min_element(gc_means.begin(), gc_means.end());
            max_m = *std::max_element(gc_means.begin(), gc_means.end());
            double min_v = *std::min_element(gc_vars.begin(), gc_vars.end());
            double max_v = *std::max_element(gc_vars.begin(), gc_vars.end());
            double range_v = (max_v > min_v) ? (max_v - min_v) : 1.0;
            norm_feats.reserve(motif_ids.size());
            for (size_t i = 0; i < motif_ids.size(); ++i)
            {
                double raw_gc = gc_means[i];
                double nv = (gc_vars[i] - min_v) / range_v; // 仅方差做 min-max 归一化
                norm_feats.emplace_back(raw_gc, nv);
            }
            // Simple KMeans on (raw_gc, nv)
            int k = std::min(cluster_k, (int)norm_feats.size());
            std::vector<std::pair<double, double>> centers(norm_feats.begin(), norm_feats.begin() + k);
            std::vector<int> labels(norm_feats.size(), 0);
            for (int iter = 0; iter < 60; ++iter)
            {
                bool changed = false;
                for (size_t i = 0; i < norm_feats.size(); ++i)
                {
                    double best = 1e300;
                    int bc = 0;
                    for (int c = 0; c < k; ++c)
                    {
                        double dx = norm_feats[i].first - centers[c].first;
                        double dy = norm_feats[i].second - centers[c].second;
                        double d = dx * dx + dy * dy;
                        if (d < best)
                        {
                            best = d;
                            bc = c;
                        }
                    }
                    if (labels[i] != bc + 1)
                    {
                        labels[i] = bc + 1;
                        changed = true;
                    }
                }
                std::vector<double> sx(k, 0.0), sy(k, 0.0);
                std::vector<int> cnt(k, 0);
                for (size_t i = 0; i < norm_feats.size(); ++i)
                {
                    int c = labels[i] - 1;
                    sx[c] += norm_feats[i].first;
                    sy[c] += norm_feats[i].second;
                    cnt[c]++;
                }
                for (int c = 0; c < k; ++c)
                {
                    if (cnt[c] > 0)
                    {
                        centers[c].first = sx[c] / cnt[c];
                        centers[c].second = sy[c] / cnt[c];
                    }
                }
                if (!changed)
                    break;
            }
            for (size_t i = 0; i < motif_ids.size(); ++i)
                cluster_assign[motif_ids[i]] = labels[i];
            std::cerr << "Clustering (raw GC mean + normalized variance): motifs=" << motif_ids.size() << ", k=" << cluster_k << ", relaxed_cluster=" << cluster_relaxed_index << "\n";
        }

        // Launch parallel calculations for all motifs (after determining cluster assignments)
        std::vector<std::future<MotifThresholds>> calculation_futures;
        calculation_futures.reserve(motif_matrices.size());
        for (const auto &kv : motif_matrices)
        {
            calculation_futures.push_back(std::async(std::launch::async, [&, id = kv.first, mat = kv.second]()
                                                     {
                auto pwm = create_pwm_from_matrix(mat);
                auto thr = calculate_motif_thresholds_parallel(id, pwm, full_threshold_percentile, core_threshold_percentile,
                                                               simulation_iterations, core_region_length, gc_threshold,
                                                               relaxed_full_pct, relaxed_core_pct, strict_full_pct, strict_core_pct);
                // apply cluster relaxation
                if (thr.calculation_successful && clustering_enabled) {
                    auto itc = cluster_assign.find(id);
                    if (itc != cluster_assign.end()) {
                        thr.cluster_index = itc->second;
                        if (thr.cluster_index == cluster_relaxed_index) {
                            if (thr.relaxed_full_threshold > 0) thr.full_motif_threshold = thr.relaxed_full_threshold;
                            if (thr.relaxed_core_threshold > 0) thr.core_region_threshold = thr.relaxed_core_threshold;
                        }
                    }
                }
                return thr; }));
        }

        // Collect results and write to output file
        std::ofstream output_file(output_threshold_file);
        if (!output_file)
        {
            throw std::runtime_error("Failed to create output file: " + output_threshold_file);
        }

        int completed_calculations = 0;
        for (auto &future : calculation_futures)
        {
            auto tr = future.get();
            if (!tr.calculation_successful)
                continue;
            output_file << tr.motif_identifier << '\t'
                        << std::fixed << std::setprecision(output_precision)
                        << tr.full_motif_threshold << '\t'
                        << tr.core_region_threshold << '\t'
                        << tr.core_start_position << '\t'
                        << tr.gc_content << '\t'
                        << tr.relaxed_full_threshold << '\t'
                        << tr.relaxed_core_threshold << '\t'
                        << tr.strict_full_threshold << '\t'
                        << tr.strict_core_threshold << '\t'
                        << tr.cluster_index << '\n';
            completed_calculations++;
        }

        // Optional cluster plot SVG
        if (clustering_enabled && !cluster_plot_path.empty() && !norm_feats.empty())
        {
            std::ofstream svg(cluster_plot_path);
            if (svg)
            {
                int W = 800, H = 600;
                double pad = 50;
                double plotW = W - 2 * pad, plotH = H - 2 * pad;
                svg << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W << "' height='" << H << "'>\n";
                svg << "<rect x='0' y='0' width='" << W << "' height='" << H << "' fill='white' stroke='none'/>\n";
                // axes
                svg << "<line x1='" << pad << "' y1='" << H - pad << "' x2='" << W - pad << "' y2='" << H - pad << "' stroke='black' stroke-width='1'/>\n";
                svg << "<line x1='" << pad << "' y1='" << H - pad << "' x2='" << pad << "' y2='" << pad << "' stroke='black' stroke-width='1'/>\n";
                svg << "<text x='" << (W / 2) << "' y='" << (H - 10) << "' font-size='14' text-anchor='middle'>GC Mean</text>\n";
                svg << "<text transform='translate(15," << (H / 2) << ") rotate(-90)' font-size='14' text-anchor='middle'>Normalized GC Variance</text>\n";
                static const char *colors[] = {"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"};
                // draw points
                double range_m_plot = (max_m > min_m) ? (max_m - min_m) : 1.0;
                for (size_t i = 0; i < norm_feats.size(); ++i)
                {
                    auto itc = cluster_assign.find(motif_ids[i]);
                    int cid = (itc != cluster_assign.end()) ? itc->second : 1;
                    cid = (cid - 1) % 10;
                    // x 使用原始 GC 均值按范围线性映射，y 使用已归一化方差
                    double x = pad + ((norm_feats[i].first - min_m) / range_m_plot) * plotW;
                    double y = H - pad - norm_feats[i].second * plotH;
                    bool relaxedC = (cluster_assign[motif_ids[i]] == cluster_relaxed_index);
                    svg << "<circle cx='" << x << "' cy='" << y << "' r='4' fill='" << colors[cid] << "' stroke='" << (relaxedC ? "black" : "none") << "' stroke-width='" << (relaxedC ? 1.5 : 0) << "' opacity='0.75'/>\n";
                }
                // legend
                double lx = W - pad - 140;
                double ly = pad + 10;
                double ldy = 18;
                int shown = std::min(cluster_k, 10);
                for (int c = 1; c <= shown; ++c)
                {
                    svg << "<rect x='" << lx << "' y='" << (ly + (c - 1) * ldy - 10) << "' width='12' height='12' fill='" << colors[(c - 1) % 10] << "' stroke='" << (c == cluster_relaxed_index ? "black" : "none") << "'/>";
                    svg << "<text x='" << (lx + 18) << "' y='" << (ly + (c - 1) * ldy) << "' font-size='12'>cluster " << c << (c == cluster_relaxed_index ? " (relaxed)" : "") << "</text>\n";
                }
                svg << "</svg>\n";
                std::cerr << "Cluster plot saved: " << cluster_plot_path << "\n";
            }
            else
            {
                std::cerr << "Warning: cannot write cluster plot to " << cluster_plot_path << "\n";
            }
        }

        auto calculation_end_time = std::chrono::high_resolution_clock::now();
        auto total_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            calculation_end_time - calculation_start_time);

        std::cerr << "Threshold calculation completed successfully!" << std::endl;
        std::cerr << "Total execution time: " << total_execution_time.count() << " ms" << std::endl;
    }
    catch (const std::exception &error)
    {
        std::cerr << "ERROR: " << error.what() << std::endl;
        return 1;
    }

    return 0;
}