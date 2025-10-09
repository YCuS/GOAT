#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <thread>
#include <future>
#include <atomic>
#include <mutex>
#include <chrono>
#include <filesystem>
#include <array>

// 全局开关：是否允许gap（对齐算法）
static bool g_allow_gaps = false;

// 背景分布与插入/缺失罚分（可调）
static constexpr std::array<double, 4> BG_UNIFORM = {0.25, 0.25, 0.25, 0.25};
static constexpr double LAMBDA_INS = -1;       // 序列侧插入每个位点罚分（nats）
static constexpr double TAU_GAP_MOT = -1;      // motif 侧缺失每列的基线罚分 ln(4)
static constexpr double MIN_MATCH_FRAC = 0.80; // 最少匹配列占比（含 core）
static constexpr int ALIGN_BAND = 128;         // band width for DP (|i-j| <= ALIGN_BAND)
static constexpr int MAX_EXT_WINDOW = 512;     // max bases to extend on each side from core

// 新增：放在任何用到 base_idx 的函数（如 logp_col）之前
static inline int base_idx(char c)
{
    switch (c)
    {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
    case 'U':
    case 'u':
        return 3;
    default:
        return -1; // 未知碱基
    }
}

// Position Weight Matrix with cache-optimized storage
struct PositionWeightMatrix
{
    std::string motif_identifier;
    std::vector<std::array<double, 4>> position_probabilities; // P
    std::vector<std::array<double, 4>> position_log_probs;     // ln(P), precomputed

    void reserve_positions(size_t position_count)
    {
        position_probabilities.reserve(position_count);
        position_log_probs.reserve(position_count);
    }
    size_t get_length() const { return position_probabilities.size(); }
};

// Threshold parameters for motif detection
struct MotifDetectionThresholds
{
    std::string motif_identifier;
    double full_motif_threshold;
    double core_region_threshold;
    int core_start_position;
    // Extended fields (optional in threshold file)
    double gc_content = 0.0;
    double full_motif_threshold_unadjusted = 0.0;
    double strict_full_threshold = 0.0;
    double strict_core_threshold = 0.0;
    double relaxed_full_threshold = 0.0;
    double relaxed_core_threshold = 0.0;
    int cluster_index = 0; // 新增：解析含 cluster_index 的扩展列
};

// DNA sequence with precomputed reverse complement
struct DNASequence
{
    std::string sequence_header;
    std::string forward_sequence;
    std::string reverse_complement;

    DNASequence(std::string header, std::string sequence)
        : sequence_header(std::move(header)), forward_sequence(std::move(sequence))
    {
        reverse_complement = compute_reverse_complement(forward_sequence);
    }

private:
    // Compute reverse complement using lookup table for performance
    static std::string compute_reverse_complement(const std::string &sequence)
    {
        std::string complement;
        complement.reserve(sequence.size());

        // Initialize complement lookup table once
        static char complement_lookup[256];
        static bool lookup_initialized = false;

        if (!lookup_initialized)
        {
            std::fill(complement_lookup, complement_lookup + 256, 0);
            complement_lookup['A'] = 'T';
            complement_lookup['T'] = 'A';
            complement_lookup['C'] = 'G';
            complement_lookup['G'] = 'C';
            complement_lookup['a'] = 't';
            complement_lookup['t'] = 'a';
            complement_lookup['c'] = 'g';
            complement_lookup['g'] = 'c';
            lookup_initialized = true;
        }

        // Build reverse complement
        for (auto nucleotide_iter = sequence.rbegin(); nucleotide_iter != sequence.rend(); ++nucleotide_iter)
        {
            char complement_base = complement_lookup[static_cast<unsigned char>(*nucleotide_iter)];
            complement += (complement_base != 0) ? complement_base : *nucleotide_iter;
        }
        return complement;
    }
};

// Motif hit result: sequence_name, motif_id, position, strand, score, matched_sequence
using MotifSearchResult = std::tuple<std::string, std::string, int, std::string, double, std::string>;

// Thread-safe collector for search results
class SearchResultCollector
{
private:
    std::vector<MotifSearchResult> collected_results;
    std::mutex collection_mutex;

public:
    void add_search_results(std::vector<MotifSearchResult> &&new_results)
    {
        std::lock_guard<std::mutex> lock(collection_mutex);
        collected_results.insert(collected_results.end(),
                                 std::make_move_iterator(new_results.begin()),
                                 std::make_move_iterator(new_results.end()));
    }

    void write_results_to_file(const std::string &output_filename)
    {
        std::ofstream output_file(output_filename);
        if (!output_file)
        {
            throw std::runtime_error("Failed to create output file: " + output_filename);
        }

        for (const auto &result : collected_results)
        {
            output_file << std::get<0>(result) << '\t'  // sequence_name
                        << std::get<1>(result) << '\t'  // motif_id
                        << std::get<2>(result) << '\t'  // position
                        << std::get<3>(result) << '\t'  // strand
                        << std::get<4>(result) << '\t'  // score
                        << std::get<5>(result) << '\n'; // matched_sequence
        }
    }

    size_t get_result_count() const { return collected_results.size(); }
};

// Load PWM definitions from file with memory preallocation
std::vector<PositionWeightMatrix> load_position_weight_matrices(const std::string &pwm_file)
{
    std::ifstream file(pwm_file);
    if (!file)
        throw std::runtime_error("Failed to open PWM file: " + pwm_file);

    std::vector<PositionWeightMatrix> pwm_collection;
    pwm_collection.reserve(100); // Reserve space for typical motif count

    PositionWeightMatrix current_pwm;
    std::string line;

    while (std::getline(file, line))
    {
        // Clean input line
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty())
            continue;

        // Remove leading/trailing whitespace
        auto content_start = line.find_first_not_of(" \t");
        auto content_end = line.find_last_not_of(" \t");
        if (content_start == std::string::npos)
            continue;
        line = line.substr(content_start, content_end - content_start + 1);

        // Process motif header
        if (line.compare(0, 6, "Motif:") == 0)
        {
            // Save previous motif if complete
            if (!current_pwm.position_probabilities.empty())
            {
                pwm_collection.push_back(std::move(current_pwm));
            }

            // Initialize new motif
            current_pwm = PositionWeightMatrix();
            current_pwm.motif_identifier = line.substr(6);
            current_pwm.reserve_positions(30); // Reserve for typical motif length
        }
        // Process probability data
        else
        {
            std::istringstream stream(line);
            double prob_A, prob_C, prob_G, prob_T;
            if (stream >> prob_A >> prob_C >> prob_G >> prob_T)
            {
                current_pwm.position_probabilities.push_back({prob_A, prob_C, prob_G, prob_T});
                const double eps = 1e-12;
                current_pwm.position_log_probs.push_back({std::log(prob_A > 0.0 ? prob_A : eps),
                                                          std::log(prob_C > 0.0 ? prob_C : eps),
                                                          std::log(prob_G > 0.0 ? prob_G : eps),
                                                          std::log(prob_T > 0.0 ? prob_T : eps)});
            }
        }
    }

    // Add final motif
    if (!current_pwm.position_probabilities.empty())
    {
        pwm_collection.push_back(std::move(current_pwm));
    }

    std::cerr << "Loaded " << pwm_collection.size() << " PWMs" << std::endl;
    return pwm_collection;
}

// Load detection thresholds from file (supports extended columns)
std::map<std::string, MotifDetectionThresholds> load_detection_thresholds(const std::string &threshold_file)
{
    std::ifstream file(threshold_file);
    if (!file)
        throw std::runtime_error("Failed to open threshold file: " + threshold_file);

    std::map<std::string, MotifDetectionThresholds> threshold_map;
    std::string line;

    while (std::getline(file, line))
    {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty())
            continue;

        std::istringstream stream(line);
        MotifDetectionThresholds thresholds;
        // 支持多种列格式：
        // A) 9列 (旧): motif_id full core start gc relaxed_full relaxed_core strict_full strict_core
        // B) 10列 (当前): 在 A 基础上追加 cluster_index
        // C) 其它历史扩展: 可能包含 full_unadjusted 等
        if (stream >> thresholds.motif_identifier >> thresholds.full_motif_threshold >>
            thresholds.core_region_threshold >> thresholds.core_start_position)
        {
            // Read the remainder of the line as tokens to support multiple formats
            std::vector<double> tokens;
            {
                double x;
                while (stream >> x)
                    tokens.push_back(x);
            }
            // Defaults
            thresholds.gc_content = 0.0;
            thresholds.relaxed_full_threshold = thresholds.full_motif_threshold;
            thresholds.relaxed_core_threshold = thresholds.core_region_threshold;
            thresholds.strict_full_threshold = thresholds.full_motif_threshold;
            thresholds.strict_core_threshold = thresholds.core_region_threshold;
            // Parse based on token count
            if (!tokens.empty())
            {
                thresholds.gc_content = tokens[0];
                // tokens size mapping:
                // size==5: gc + relaxed_full + relaxed_core + strict_full + strict_core
                // size==6: gc + relaxed_full + relaxed_core + strict_full + strict_core + cluster_index
                // size==6 (old variant with full_unadjusted) ambiguity: detect by ordering heuristics
                // size==6 legacy (gc, full_unadj, strict_full, strict_core, relaxed_full, relaxed_core) -> pattern where tokens[1] > tokens[4] typically (full_unadj > relaxed_full). We'll attempt heuristic.
                if (tokens.size() == 5)
                {
                    thresholds.relaxed_full_threshold = tokens[1];
                    thresholds.relaxed_core_threshold = tokens[2];
                    thresholds.strict_full_threshold = tokens[3];
                    thresholds.strict_core_threshold = tokens[4];
                }
                else if (tokens.size() == 6)
                {
                    // Heuristic to distinguish legacy vs new with cluster_index
                    bool looks_legacy = (tokens[1] > tokens[4]) && (tokens[2] > tokens[5]);
                    if (looks_legacy)
                    {
                        // legacy ordering: gc, full_unadj, strict_full, strict_core, relaxed_full, relaxed_core
                        thresholds.full_motif_threshold_unadjusted = tokens[1];
                        thresholds.strict_full_threshold = tokens[2];
                        thresholds.strict_core_threshold = tokens[3];
                        thresholds.relaxed_full_threshold = tokens[4];
                        thresholds.relaxed_core_threshold = tokens[5];
                    }
                    else
                    {
                        // new ordering with cluster index: gc, relaxed_full, relaxed_core, strict_full, strict_core, cluster_index
                        thresholds.relaxed_full_threshold = tokens[1];
                        thresholds.relaxed_core_threshold = tokens[2];
                        thresholds.strict_full_threshold = tokens[3];
                        thresholds.strict_core_threshold = tokens[4];
                        thresholds.cluster_index = static_cast<int>(std::round(tokens[5]));
                    }
                }
                else if (tokens.size() == 4)
                {
                    // gc + relaxed_full + relaxed_core + strict_full (missing strict_core)
                    thresholds.relaxed_full_threshold = tokens[1];
                    thresholds.relaxed_core_threshold = tokens[2];
                    thresholds.strict_full_threshold = tokens[3];
                }
                else if (tokens.size() == 3)
                {
                    // gc + relaxed_full + relaxed_core
                    thresholds.relaxed_full_threshold = tokens[1];
                    thresholds.relaxed_core_threshold = tokens[2];
                }
                else if (tokens.size() >= 7)
                {
                    // If future extended format appears, take first 5 meaningful after gc per new spec
                    thresholds.relaxed_full_threshold = tokens[1];
                    thresholds.relaxed_core_threshold = tokens[2];
                    thresholds.strict_full_threshold = tokens[3];
                    thresholds.strict_core_threshold = tokens[4];
                    thresholds.cluster_index = static_cast<int>(std::round(tokens[5]));
                }
            }
            threshold_map[thresholds.motif_identifier] = thresholds;
        }
    }

    std::cerr << "Loaded " << threshold_map.size() << " threshold sets" << std::endl;
    return threshold_map;
}

// Load DNA sequences from FASTA file
std::vector<DNASequence> load_dna_sequences_from_fasta(const std::string &fasta_file)
{
    std::ifstream file(fasta_file);
    if (!file)
        throw std::runtime_error("Failed to open sequence file: " + fasta_file);

    std::vector<DNASequence> sequence_collection;
    std::string line, current_header, current_sequence;

    while (std::getline(file, line))
    {
        if (line.empty())
            continue;

        if (line[0] == '>')
        {
            // Save previous sequence if complete
            if (!current_header.empty() && !current_sequence.empty())
            {
                sequence_collection.emplace_back(std::move(current_header), std::move(current_sequence));
            }

            // Start new sequence
            current_header = line.substr(1);
            current_sequence.clear();
        }
        else
        {
            current_sequence += line;
        }
    }

    // Add final sequence
    if (!current_header.empty() && !current_sequence.empty())
    {
        sequence_collection.emplace_back(std::move(current_header), std::move(current_sequence));
    }

    std::cerr << "Loaded " << sequence_collection.size() << " DNA sequences" << std::endl;
    return sequence_collection;
}

// Calculate log-likelihood score using optimized lookup table
double calculate_motif_score_optimized(const std::string &sequence,
                                       int start_position,
                                       const PositionWeightMatrix &pwm)
{
    double total_score = 0.0;
    // Nucleotide to array index mapping: A=0, C=1, G=2, T=3
    static int nucleotide_to_index[256];
    static bool lookup_initialized = false;

    if (!lookup_initialized)
    {
        std::fill(nucleotide_to_index, nucleotide_to_index + 256, -1);
        nucleotide_to_index['A'] = 0;
        nucleotide_to_index['a'] = 0;
        nucleotide_to_index['C'] = 1;
        nucleotide_to_index['c'] = 1;
        nucleotide_to_index['G'] = 2;
        nucleotide_to_index['g'] = 2;
        nucleotide_to_index['T'] = 3;
        nucleotide_to_index['t'] = 3;
        nucleotide_to_index['U'] = 3;
        nucleotide_to_index['u'] = 3;
        lookup_initialized = true;
    }

    // Calculate log-likelihood score for each position
    for (size_t position = 0; position < pwm.position_probabilities.size(); ++position)
    {
        int nucleotide_index = nucleotide_to_index[static_cast<unsigned char>(sequence[start_position + position])];
        if (nucleotide_index >= 0)
        {
            double probability = pwm.position_probabilities[position][nucleotide_index];
            if (probability > 0.0)
            {
                total_score += -std::log(probability);
            }
            else
            {
                std::cerr << "Warning: PWM for motif '"
                          << pwm.motif_identifier
                          << "' contains a zero probability at position "
                          << position
                          << ". Consider using Bayesian mode for more robust scoring."
                          << std::endl;
                total_score += -std::log(0.000001); // Use a small fallback value
            }
        }
    }

    return total_score;
}

// 取某列 PWM 对某碱基的 log 概率（已保证无 0，直接取 log）
static inline double logp_col(const std::array<double, 4> &row, char base)
{
    int k = base_idx(base);
    if (k < 0)
        return std::log(0.25);
    return std::log(row[k]);
}

// 使用已为该列预先取对数后的行（row_ln）版本，更快
static inline double logp_col_ln(const std::array<double, 4> &row_ln, char base)
{
    int k = base_idx(base);
    if (k < 0)
        return std::log(0.25);
    return row_ln[k];
}

// Return PWM column probability for a base. Unknown base -> 0.25 fallback.
static inline double prob_col(const std::array<double, 4> &row, char base)
{
    int k = base_idx(base);
    if (k < 0)
        return 0.25; // unknown base fallback
    double p = row[k];
    if (p < 0.0)
        p = 0.0; // guard (shouldn't happen)
    return p;
}

// 对齐阶段：最大化得分（匹配加 log p，gap 不加分不扣分）；终止取 dp[i][M] 最大
struct AlignResult
{
    double cost;               // 此处语义变为“得分”（sum of log p），保留字段名不改
    std::string aligned_motif; // 长度 == motif_seg.size()，motif gap 用 '-'
    int motif_gaps;
    int seq_gaps;
};

static AlignResult
align_free_gaps(const std::vector<std::array<double, 4>> &motif_seg,
                const std::string &seq_seg)
{
    const int M = static_cast<int>(motif_seg.size());
    const int N = static_cast<int>(seq_seg.size());
    if (M == 0)
        return {0.0, std::string(), 0, 0};

    const size_t W = static_cast<size_t>(M) + 1;
    auto ID = [W](int i, int j) -> size_t
    { return static_cast<size_t>(i) * W + static_cast<size_t>(j); };

    const double NEG_INF = -1e300;
    std::vector<double> dp((static_cast<size_t>(N) + 1) * W, NEG_INF);
    std::vector<char> op((static_cast<size_t>(N) + 1) * W, 0);
    dp[ID(0, 0)] = 0.0;

    for (int j = 1; j <= M; ++j)
    {
        dp[ID(0, j)] = dp[ID(0, j - 1)] + TAU_GAP_MOT;
        op[ID(0, j)] = 'G';
    }
    for (int i = 1; i <= N; ++i)
    {
        dp[ID(i, 0)] = dp[ID(i - 1, 0)] + LAMBDA_INS;
        op[ID(i, 0)] = 'S';
    }

    for (int i = 1; i <= N; ++i)
    {
        // band: only compute j in [i-ALIGN_BAND, i+ALIGN_BAND]
        const int jmin = std::max(1, i - ALIGN_BAND);
        const int jmax = std::min(M, i + ALIGN_BAND);
        const char base = seq_seg[i - 1];

        for (int j = jmin; j <= jmax; ++j)
        {
            double s_match = dp[ID(i - 1, j - 1)] + prob_col(motif_seg[j - 1], base);
            double s_mgap = dp[ID(i, j - 1)] + TAU_GAP_MOT;
            double s_sgap = dp[ID(i - 1, j)] + LAMBDA_INS;

            double best = s_match;
            char bop = 'M';
            if (s_mgap >= best)
            {
                best = s_mgap;
                bop = 'G';
            }
            if (s_sgap > best)
            {
                best = s_sgap;
                bop = 'S';
            }

            dp[ID(i, j)] = best;
            op[ID(i, j)] = bop;
        }
    }

    double best = NEG_INF;
    int i_best = 0;
    for (int i = 0; i <= N; ++i)
    {
        double v = dp[ID(i, M)];
        if (v > best)
        {
            best = v;
            i_best = i;
        }
    }
    if (best <= NEG_INF / 2)
        return {best, std::string(), 0, 0};

    std::string aligned(M, '-');
    int motif_gaps = 0, seq_gaps = 0;
    int i = i_best, j = M;
    while (j > 0 || i > 0)
    {
        char bop = (i >= 0 && j >= 0) ? op[ID(i, j)] : 0;
        if (j > 0 && i > 0 && bop == 'M')
        {
            aligned[j - 1] = seq_seg[i - 1];
            --i;
            --j;
        }
        else if (j > 0 && bop == 'G')
        {
            ++motif_gaps;
            --j;
        }
        else if (i > 0 && bop == 'S')
        {
            ++seq_gaps;
            --i;
        }
        else
        {
            if (j > 0)
            {
                ++motif_gaps;
                --j;
            }
            else if (i > 0)
            {
                ++seq_gaps;
                --i;
            }
            else
                break;
        }
    }
    return {best, aligned, motif_gaps, seq_gaps};
}

// 基于 core 锚定的双侧对齐；回溯后再计算最终分数：
// final = Σ(-log p_j(b_j)) + 2 * (motif_gaps + seq_gaps)
static std::pair<double, std::string>
score_with_core_anchored_alignment(const std::string &seq,
                                   int p_core,
                                   const PositionWeightMatrix &pwm,
                                   int core_start, int core_len)
{
    const int L = (int)pwm.get_length();
    const int n = (int)seq.size();
    const int core_end = core_start + core_len;
    if (p_core < 0 || p_core + core_len > n)
        return {INFINITY, std::string()};

    // Left side window (limit extension)
    std::string left_aln;
    int gaps_left = 0;
    if (core_start > 0)
    {
        int left_len = std::min(p_core, std::min(core_start + 2 * core_len, MAX_EXT_WINDOW));
        int left_start = p_core - left_len;
        std::string left_seq = seq.substr(left_start, left_len);
        std::string left_seq_rev(left_seq.rbegin(), left_seq.rend());

        std::vector<std::array<double, 4>> motif_left_rev;
        motif_left_rev.reserve(core_start);
        for (int j = core_start - 1; j >= 0; --j)
            motif_left_rev.push_back(pwm.position_probabilities[j]);

        auto lr = align_free_gaps(motif_left_rev, left_seq_rev);
        if (!lr.aligned_motif.empty())
        {
            left_aln.assign(lr.aligned_motif.rbegin(), lr.aligned_motif.rend());
            gaps_left = lr.motif_gaps + lr.seq_gaps;
        }
        else
            return {INFINITY, std::string()};
    }

    // Core (contiguous)
    std::string core_aln;
    core_aln.reserve(core_len);
    for (int t = 0; t < core_len; ++t)
    {
        core_aln.push_back(seq[p_core + t]);
    }

    // Right side window (limit extension)
    std::string right_aln;
    int gaps_right = 0;
    if (core_end < L)
    {
        int right_avail = n - (p_core + core_len);
        int right_len = std::min(right_avail, std::min(L - core_end + 2 * core_len, MAX_EXT_WINDOW));
        std::string right_seq = seq.substr(p_core + core_len, right_len);

        std::vector<std::array<double, 4>> motif_right;
        motif_right.reserve(L - core_end);
        for (int j = core_end; j < L; ++j)
            motif_right.push_back(pwm.position_probabilities[j]);

        auto rr = align_free_gaps(motif_right, right_seq);
        if (!rr.aligned_motif.empty())
        {
            right_aln = std::move(rr.aligned_motif);
            gaps_right = rr.motif_gaps + rr.seq_gaps;
        }
        else
            return {INFINITY, std::string()};
    }

    // Join and minimum matched columns
    std::string aligned = left_aln + core_aln + right_aln;
    int matched_cols = 0;
    for (char c : aligned)
        if (c != '-')
            ++matched_cols;
    int min_required = std::max(core_len, (int)std::ceil(L * MIN_MATCH_FRAC));
    if (matched_cols < min_required)
        return {INFINITY, std::string()};

    // Final score using precomputed ln(P): sum(-ln p_j(b)) + gap penalty
    double negloglik = 0.0;
    for (int j = 0; j < L; ++j)
    {
        char c = aligned[j];
        if (c == '-')
            continue;
        int k = base_idx(c);
        const double lnP = (k < 0) ? std::log(0.25) : pwm.position_log_probs[j][k];
        negloglik += (-lnP);
    }
    int gap_total = gaps_left + gaps_right;
    double final_score = negloglik + 4.0 * static_cast<double>(gap_total);

    return {final_score, aligned};
}

// 无gap评分：使用与最终评分一致的代价 Σ(-ln p_j(b))，不添加gap惩罚
static std::pair<double, std::string>
score_without_gaps(const std::string &seq,
                   int p_core,
                   const PositionWeightMatrix &pwm,
                   int core_start, int core_len)
{
    const int L = static_cast<int>(pwm.get_length());
    const int n = static_cast<int>(seq.size());
    // 以 core_start 对齐：序列窗口起点
    const int s0 = p_core - core_start;
    if (s0 < 0 || s0 + L > n)
    {
        return {INFINITY, std::string()}; // 边界不足时无gap模式不可匹配
    }
    double negloglik = 0.0;
    for (int j = 0; j < L; ++j)
    {
        char b = seq[s0 + j];
        // 使用预计算 ln(P)
        negloglik += -logp_col_ln(pwm.position_log_probs[j], b);
    }
    // 无gap惩罚
    std::string aligned = seq.substr(s0, L);
    return {negloglik, aligned};
}

// 前置声明：供 lambda 使用
static double core_cost_contiguous(const std::string &seq, int p_core,
                                   const PositionWeightMatrix &pwm,
                                   int core_start, int core_len);

// Search for motif matches in a single sequence using sliding window
std::vector<MotifSearchResult> search_motif_in_sequence(
    const DNASequence &target_sequence,
    const PositionWeightMatrix &search_pwm,
    const MotifDetectionThresholds &detection_thresholds,
    int core_region_length)
{
    std::vector<MotifSearchResult> sequence_matches;
    sequence_matches.reserve(100);

    const auto &forward_seq = target_sequence.forward_sequence;
    const auto &reverse_seq = target_sequence.reverse_complement;
    const auto &seq_name = target_sequence.sequence_header;

    const int L = static_cast<int>(search_pwm.get_length());
    const int core_start = detection_thresholds.core_start_position;
    const int core_len = core_region_length;
    (void)L; // unused in this function body after refactor

    // Two-pass scan per strand:
    // 1) Effective pass: use GC-effective thresholds (full/core) from thr.txt
    // 2) If any consecutive hits (p_core and p_core+1) exist for this motif on the strand, re-scan using strict thresholds only
    auto pass_scan = [&](const std::string &strand, const std::string &strand_name,
                         bool strict_only) -> std::pair<std::vector<MotifSearchResult>, std::vector<int>>
    {
        std::vector<MotifSearchResult> hits;
        std::vector<int> p_core_list;
        const int n = static_cast<int>(strand.size());
        for (int p_core = 0; p_core + core_len <= n; ++p_core)
        {
            // Core cost first
            double core_cost = core_cost_contiguous(strand, p_core, search_pwm, core_start, core_len);
            if (strict_only)
            {
                if (!(core_cost <= detection_thresholds.strict_core_threshold))
                    continue;
            }
            else
            {
                if (!(core_cost <= detection_thresholds.core_region_threshold))
                    continue;
            }

            // Score full motif at this alignment
            std::pair<double, std::string> scored;
            if (g_allow_gaps)
                scored = score_with_core_anchored_alignment(strand, p_core, search_pwm, core_start, core_len);
            else
                scored = score_without_gaps(strand, p_core, search_pwm, core_start, core_len);
            double final_score = scored.first;

            if (strict_only)
            {
                if (!(final_score <= detection_thresholds.strict_full_threshold))
                    continue;
            }
            else
            {
                if (!(final_score <= detection_thresholds.full_motif_threshold))
                    continue;
            }

            // Emit hit
            int hit_pos = p_core - core_start;
            std::string matched_seq = scored.second.empty()
                                          ? strand.substr(std::max(0, hit_pos), std::min((int)search_pwm.get_length(), n - std::max(0, hit_pos)))
                                          : scored.second;
            hits.emplace_back(seq_name, search_pwm.motif_identifier, hit_pos, strand_name, final_score, matched_seq);
            p_core_list.push_back(p_core);
        }
        return {std::move(hits), std::move(p_core_list)};
    };

    auto process_strand_two_pass = [&](const std::string &strand, const std::string &strand_name)
    {
        // Pass 1: effective thresholds (already GC-aware)
        auto [hits1, pcores1] = pass_scan(strand, strand_name, /*strict_only=*/false);
        bool has_consecutive = false;
        if (!pcores1.empty())
        {
            std::sort(pcores1.begin(), pcores1.end());
            for (size_t i = 1; i < pcores1.size(); ++i)
            {
                if (pcores1[i] == pcores1[i - 1] + 1)
                {
                    has_consecutive = true;
                    break;
                }
            }
        }

        if (has_consecutive)
        {
            // 触发严格阈值条件：无论是否有严格命中，都放弃第一轮结果（保持优先级一致性）
            auto [hits2, _] = pass_scan(strand, strand_name, /*strict_only=*/true);
            if (!hits2.empty())
            {
                sequence_matches.insert(sequence_matches.end(),
                                        std::make_move_iterator(hits2.begin()),
                                        std::make_move_iterator(hits2.end()));
            }
            return; // 不回退
        }
        else if (!hits1.empty())
        {
            sequence_matches.insert(sequence_matches.end(),
                                    std::make_move_iterator(hits1.begin()),
                                    std::make_move_iterator(hits1.end()));
        }
    };

    process_strand_two_pass(forward_seq, "positive");
    process_strand_two_pass(reverse_seq, "negative");

    return sequence_matches;
}

// Process multiple sequences in parallel using thread pool
void process_sequences_in_parallel(
    const std::vector<DNASequence> &sequence_collection,
    const std::vector<PositionWeightMatrix> &pwm_collection,
    const std::map<std::string, MotifDetectionThresholds> &threshold_map,
    int core_region_length,
    SearchResultCollector &result_collector)
{

    const int thread_count = std::thread::hardware_concurrency();
    std::vector<std::future<void>> processing_futures;

    // Process sequences in chunks across threads
    auto process_sequence_chunk = [&](size_t chunk_start, size_t chunk_end)
    {
        for (size_t seq_index = chunk_start; seq_index < chunk_end; ++seq_index)
        {
            const auto &sequence = sequence_collection[seq_index];

            // Search each motif in current sequence
            for (const auto &pwm : pwm_collection)
            {
                auto threshold_iter = threshold_map.find(pwm.motif_identifier);
                if (threshold_iter != threshold_map.end())
                {
                    auto matches = search_motif_in_sequence(sequence, pwm, threshold_iter->second, core_region_length);
                    if (!matches.empty())
                    {
                        result_collector.add_search_results(std::move(matches));
                    }
                }
            }
        }
    };

    // Distribute work across available threads
    size_t sequences_per_thread = (sequence_collection.size() + thread_count - 1) / thread_count;
    for (int thread_id = 0; thread_id < thread_count; ++thread_id)
    {
        size_t chunk_start = thread_id * sequences_per_thread;
        size_t chunk_end = std::min(chunk_start + sequences_per_thread, sequence_collection.size());
        if (chunk_start < chunk_end)
        {
            processing_futures.push_back(std::async(std::launch::async,
                                                    process_sequence_chunk,
                                                    chunk_start, chunk_end));
        }
    }

    // Wait for all threads to complete
    for (auto &future : processing_futures)
    {
        future.wait();
    }
}

// Process motifs in parallel: split pwm_collection across threads, each scans all sequences for its motif chunk
void process_motifs_in_parallel(
    const std::vector<DNASequence> &sequence_collection,
    const std::vector<PositionWeightMatrix> &pwm_collection,
    const std::map<std::string, MotifDetectionThresholds> &threshold_map,
    int core_region_length,
    SearchResultCollector &result_collector)
{

    unsigned int hw = std::thread::hardware_concurrency();
    const int thread_count = hw ? static_cast<int>(hw) : 4;
    std::vector<std::future<void>> futures;

    auto process_motif_chunk = [&](size_t m_begin, size_t m_end)
    {
        for (size_t mi = m_begin; mi < m_end; ++mi)
        {
            const auto &pwm = pwm_collection[mi];
            auto it = threshold_map.find(pwm.motif_identifier);
            if (it == threshold_map.end())
                continue;
            const auto &thr = it->second;

            // scan all sequences for this motif
            for (const auto &sequence : sequence_collection)
            {
                auto matches = search_motif_in_sequence(sequence, pwm, thr, core_region_length);
                if (!matches.empty())
                {
                    result_collector.add_search_results(std::move(matches));
                }
            }
        }
    };

    size_t motifs_per_thread = (pwm_collection.size() + thread_count - 1) / thread_count;
    for (int t = 0; t < thread_count; ++t)
    {
        size_t m_begin = static_cast<size_t>(t) * motifs_per_thread;
        size_t m_end = std::min(m_begin + motifs_per_thread, pwm_collection.size());
        if (m_begin < m_end)
        {
            futures.push_back(std::async(std::launch::async, process_motif_chunk, m_begin, m_end));
        }
    }
    for (auto &f : futures)
        f.wait();
}

// Auto choose parallel strategy: by sequences (existing) or by motifs (new), whichever yields larger chunks
void process_in_parallel_auto(
    const std::vector<DNASequence> &sequence_collection,
    const std::vector<PositionWeightMatrix> &pwm_collection,
    const std::map<std::string, MotifDetectionThresholds> &threshold_map,
    int core_region_length,
    SearchResultCollector &result_collector)
{

    if (pwm_collection.size() >= sequence_collection.size())
    {
        // many motifs or few sequences: parallel by motif
        process_motifs_in_parallel(sequence_collection, pwm_collection, threshold_map,
                                   core_region_length, result_collector);
    }
    else
    {
        // many sequences or few motifs: parallel by sequence (existing)
        process_sequences_in_parallel(sequence_collection, pwm_collection, threshold_map,
                                      core_region_length, result_collector);
    }
}

static std::vector<std::string> get_sequence_files_from_directory(const std::string &dir_path)
{
    if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path))
    {
        throw std::runtime_error("Sequence directory does not exist or is not a directory: " + dir_path);
    }
    std::vector<std::string> files;
    for (auto &entry : std::filesystem::directory_iterator(dir_path))
    {
        if (!entry.is_regular_file())
            continue;
        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".fa" || ext == ".fasta" || ext == ".fna")
        {
            files.push_back(entry.path().string());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

static std::string generate_output_filename(const std::string &input_file, const std::string &output_dir)
{
    std::filesystem::path in(input_file);
    auto stem = in.stem().string();
    std::filesystem::path out_dir(output_dir);
    // output file format: <stem>_results.txt
    return (out_dir / (stem + "_results.txt")).string();
}

// 计算碱基 c 在 PWM 某列的 log 概率（含平滑处理）
// 取值：[0, -∞)，越大表示越可能是正样本；未定义碱基返回 0.25 的均匀先验 log 概率。
static inline double logp_col_base(const std::array<double, 4> &col, char c)
{
    int k = base_idx(c);
    const double eps = 1e-12;
    if (k < 0)
        return std::log(0.25);
    double p = col[k];
    if (p <= 0.0)
        p = eps;
    return std::log(p);
}

// 计算 core 段（连续、不含 gap）在 seq[p_core .. p_core+core_len) 上的 -log 概率和（cost）
static double core_cost_contiguous(const std::string &seq, int p_core,
                                   const PositionWeightMatrix &pwm,
                                   int core_start, int core_len)
{
    const int n = static_cast<int>(seq.size());
    const int L = static_cast<int>(pwm.get_length());
    if (core_len <= 0 || core_start < 0 || core_start + core_len > L)
        return INFINITY;
    if (p_core < 0 || p_core + core_len > n)
        return INFINITY;

    double s = 0.0;
    for (int t = 0; t < core_len; ++t)
    {
        const int j = core_start + t;
        const char c = seq[p_core + t];
        const int k = base_idx(c);
        const double lnP = (k < 0) ? std::log(0.25) : pwm.position_log_probs[j][k];
        s += (-lnP);
    }
    return s;
}

int main(int argc, char *argv[])
{
    // 新增 --cluster-file=path : TSV (motif_id<TAB>cluster_index)
    if (argc < 6)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <seq_path> <pwm.txt> <thr.txt> <out> <core_len> [--batch] [--no-gaps|--allow-gaps=0|1] [--cluster-file=path]\n";
        return 1;
    }

    bool batch_mode = false;
    std::string cluster_file;
    std::map<std::string, int> motif_cluster;

    // 解析必选参数
    std::string seq_path = argv[1];
    std::string pwm_file = argv[2];
    std::string thr_file = argv[3];
    std::string out_path = argv[4];
    int core_len = std::stoi(argv[5]);

    // 解析可选flags
    for (int i = 6; i < argc; ++i)
    {
        std::string flag = argv[i];
        if (flag == "--batch")
        {
            batch_mode = true;
        }
        else if (flag == "--no-gaps")
        {
            g_allow_gaps = false;
        }
        else if (flag.rfind("--allow-gaps=", 0) == 0)
        {
            std::string v = flag.substr(std::string("--allow-gaps=").size());
            g_allow_gaps = !(v == "0" || v == "false" || v == "False");
        }
        else if (flag.rfind("--cluster-file=", 0) == 0)
        {
            cluster_file = flag.substr(std::string("--cluster-file=").size());
        }
        else
        {
            std::cerr << "Warning: unknown flag ignored: " << flag << "\n";
        }
    }

    try
    {
        std::string seq_path = argv[1];
        std::string pwm_file = argv[2];
        std::string thr_file = argv[3];
        std::string out_path = argv[4];
        int core_region_length = std::stoi(argv[5]);
        bool batch_mode = (argc == 7 && std::string(argv[6]) == "--batch") || std::filesystem::is_directory(seq_path);

        // Preload PWM definitions and thresholds; no need to reload for each file in batch mode
        auto pwm_collection = load_position_weight_matrices(pwm_file);
        auto threshold_map = load_detection_thresholds(thr_file);
        // 读取 cluster 文件并调整阈值：仅 cluster==3 使用 relaxed 阈值，其余使用常规（unadjusted）
        if (!cluster_file.empty())
        {
            std::ifstream cf(cluster_file);
            if (!cf)
            {
                std::cerr << "Warning: cannot open cluster file: " << cluster_file << "\n";
            }
            else
            {
                std::string line;
                while (std::getline(cf, line))
                {
                    if (line.empty())
                        continue;
                    std::istringstream iss(line);
                    std::string mid;
                    int cid;
                    if (iss >> mid >> cid)
                        motif_cluster[mid] = cid;
                }
                std::cerr << "Loaded " << motif_cluster.size() << " cluster assignments" << std::endl;
            }
            for (auto &kv : threshold_map)
            {
                auto &thr = kv.second;
                auto itc = motif_cluster.find(kv.first);
                if (itc != motif_cluster.end())
                {
                    if (itc->second == 3)
                    {
                        if (thr.relaxed_full_threshold > 0)
                            thr.full_motif_threshold = thr.relaxed_full_threshold;
                        if (thr.relaxed_core_threshold > 0)
                            thr.core_region_threshold = thr.relaxed_core_threshold;
                    }
                    else
                    {
                        // 其它 cluster 恢复为未调整阈值（full_motif_threshold_unadjusted 已保存）
                        if (thr.full_motif_threshold_unadjusted > 0)
                            thr.full_motif_threshold = thr.full_motif_threshold_unadjusted;
                    }
                }
            }
        }
        if (pwm_collection.empty() || threshold_map.empty())
        {
            throw std::runtime_error("Failed to load PWM definitions or thresholds; they may be empty");
        }

        if (batch_mode)
        {
            // Batch mode: seq_path is a directory, out_path is created if it does not exist
            if (!std::filesystem::exists(out_path))
            {
                std::filesystem::create_directories(out_path);
            }

            auto seq_files = get_sequence_files_from_directory(seq_path);
            if (seq_files.empty())
            {
                std::cerr << "No FASTA files found in directory: " << seq_path << "\n";
                return 1;
            }
            std::cerr << "Batch mode: found " << seq_files.size() << " FASTA files\n";

            int idx = 0;
            for (const auto &seq_file : seq_files)
            {
                ++idx;
                std::string result_file = generate_output_filename(seq_file, out_path);
                std::cerr << "Processing [" << idx << "/" << seq_files.size() << "]: "
                          << std::filesystem::path(seq_file).filename().string()
                          << " -> " << std::filesystem::path(result_file).filename().string()
                          << "\n";

                // Load sequences from the current FASTA file
                auto sequences = load_dna_sequences_from_fasta(seq_file);
                if (sequences.empty())
                {
                    std::cerr << "  ⚠️  No sequences found in " << seq_file << ", skipping\n";
                    continue;
                }

                // Perform parallel motif search
                SearchResultCollector collector;
                auto t0 = std::chrono::high_resolution_clock::now();
                process_in_parallel_auto(sequences, pwm_collection, threshold_map,
                                         core_region_length, collector);
                auto t1 = std::chrono::high_resolution_clock::now();

                // Write results to file
                collector.write_results_to_file(result_file);
                auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
                std::cerr << "  ✓ Done, matched " << collector.get_result_count()
                          << " entries, took " << ms << " ms\n";
            }
        }
        else
        {
            auto sequences = load_dna_sequences_from_fasta(seq_path);
            if (sequences.empty())
            {
                throw std::runtime_error("Can't find sequences in " + seq_path);
            }
            std::cerr << "Loaded " << sequences.size() << " sequences\n";

            SearchResultCollector collector;
            auto t0 = std::chrono::high_resolution_clock::now();
            process_in_parallel_auto(sequences, pwm_collection, threshold_map,
                                     core_region_length, collector);
            auto t1 = std::chrono::high_resolution_clock::now();

            collector.write_results_to_file(out_path);
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
            std::cerr << "Motif search completed successfully!\n";
            std::cerr << "Execution time: " << ms << " ms\n";
            std::cerr << "Total matches found: " << collector.get_result_count() << "\n";
            std::cerr << "Results saved to: " << out_path << "\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}