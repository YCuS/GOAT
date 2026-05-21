#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

struct Pwm {
    std::string id;
    std::vector<std::array<double, 4>> columns; // A, C, G, T
};

struct Threshold {
    double full = 0.0;
    double core = 0.0;
    int core_start = 0;   // 0-based
    int motif_length = 0; // optional metadata
};

struct Sequence {
    std::string id;
    std::string bases;
};

struct Hit {
    std::string sequence_id;
    std::string motif_id;
    char strand = '+';
    int start = 0; // 1-based inclusive on the input sequence
    int end = 0;   // 1-based inclusive on the input sequence
    double full_score = 0.0;
    double core_score = 0.0;
    std::string matched_sequence;
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
    const double sum = row[0] + row[1] + row[2] + row[3];
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

static std::map<std::string, Threshold> load_thresholds(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open threshold file: " + path);

    std::map<std::string, Threshold> thresholds;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        std::istringstream row(line);
        std::string motif_id;
        Threshold threshold;
        if (row >> motif_id >> threshold.full >> threshold.core >> threshold.core_start) {
            row >> threshold.motif_length;
            thresholds[motif_id] = threshold;
        }
    }
    return thresholds;
}

static std::vector<Sequence> load_fasta(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open FASTA file: " + path);

    std::vector<Sequence> sequences;
    std::string line;
    Sequence current;

    auto flush_current = [&]() {
        if (!current.id.empty() && !current.bases.empty()) {
            sequences.push_back(std::move(current));
        }
        current = Sequence{};
    };

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        if (line[0] == '>') {
            flush_current();
            current.id = trim(line.substr(1));
        } else {
            current.bases += line;
        }
    }
    flush_current();
    return sequences;
}

static std::string reverse_complement(const std::string &sequence) {
    std::string rc;
    rc.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
        case 'A': case 'a': rc.push_back('T'); break;
        case 'C': case 'c': rc.push_back('G'); break;
        case 'G': case 'g': rc.push_back('C'); break;
        case 'T': case 't':
        case 'U': case 'u': rc.push_back('A'); break;
        default: rc.push_back('N'); break;
        }
    }
    return rc;
}

static double score_window(const std::string &sequence,
                           int sequence_start,
                           const Pwm &pwm,
                           int pwm_start,
                           int length) {
    static constexpr double eps = 1e-12;
    double score = 0.0;

    for (int i = 0; i < length; ++i) {
        const int idx = base_index(sequence[sequence_start + i]);
        const double probability = (idx >= 0) ? pwm.columns[pwm_start + i][idx] : 0.25;
        score += -std::log(std::max(probability, eps));
    }
    return score;
}

static void scan_strand(const Sequence &sequence,
                        const std::string &strand_sequence,
                        char strand,
                        const Pwm &pwm,
                        const Threshold &threshold,
                        int core_len,
                        std::vector<Hit> &hits) {
    const int motif_len = static_cast<int>(pwm.columns.size());
    const int seq_len = static_cast<int>(strand_sequence.size());
    if (motif_len <= 0 || core_len <= 0 || threshold.core_start < 0) return;
    if (threshold.core_start + core_len > motif_len) return;
    if (motif_len > seq_len) return;

    for (int pos = 0; pos + motif_len <= seq_len; ++pos) {
        const double core_score =
            score_window(strand_sequence, pos + threshold.core_start, pwm, threshold.core_start, core_len);
        if (core_score > threshold.core) continue;

        const double full_score = score_window(strand_sequence, pos, pwm, 0, motif_len);
        if (full_score > threshold.full) continue;

        Hit hit;
        hit.sequence_id = sequence.id;
        hit.motif_id = pwm.id;
        hit.strand = strand;
        hit.full_score = full_score;
        hit.core_score = core_score;
        hit.matched_sequence = strand_sequence.substr(pos, motif_len);

        if (strand == '+') {
            hit.start = pos + 1;
            hit.end = pos + motif_len;
        } else {
            hit.start = static_cast<int>(sequence.bases.size()) - pos - motif_len + 1;
            hit.end = static_cast<int>(sequence.bases.size()) - pos;
        }
        hits.push_back(std::move(hit));
    }
}

static std::vector<Hit> search_sequences(const std::vector<Sequence> &sequences,
                                         const std::vector<Pwm> &motifs,
                                         const std::map<std::string, Threshold> &thresholds,
                                         int core_len) {
    std::vector<Hit> hits;
    for (const auto &sequence : sequences) {
        const std::string rc = reverse_complement(sequence.bases);
        for (const auto &pwm : motifs) {
            const auto threshold_it = thresholds.find(pwm.id);
            if (threshold_it == thresholds.end()) continue;

            scan_strand(sequence, sequence.bases, '+', pwm, threshold_it->second, core_len, hits);
            scan_strand(sequence, rc, '-', pwm, threshold_it->second, core_len, hits);
        }
    }

    std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
        return std::tie(a.sequence_id, a.start, a.end, a.motif_id, a.strand, a.full_score) <
               std::tie(b.sequence_id, b.start, b.end, b.motif_id, b.strand, b.full_score);
    });
    return hits;
}

static void write_hits(const std::string &path, const std::vector<Hit> &hits) {
    const std::filesystem::path out_path(path);
    if (!out_path.parent_path().empty()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot create output file: " + path);

    out << "sequence_id\tmotif_id\tstrand\tstart\tend\tfull_score\tcore_score\tmatched_sequence\n";
    for (const auto &hit : hits) {
        out << hit.sequence_id << '\t'
            << hit.motif_id << '\t'
            << hit.strand << '\t'
            << hit.start << '\t'
            << hit.end << '\t'
            << std::fixed << std::setprecision(6)
            << hit.full_score << '\t'
            << hit.core_score << '\t'
            << hit.matched_sequence << '\n';
    }
}

static bool is_fasta_file(const std::filesystem::path &path) {
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return ext == ".fa" || ext == ".fasta" || ext == ".fna";
}

static std::vector<std::string> list_fasta_files(const std::string &dir_path) {
    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(dir_path)) {
        if (entry.is_regular_file() && is_fasta_file(entry.path())) {
            files.push_back(entry.path().string());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

static std::string batch_output_name(const std::string &input_file, const std::string &output_dir) {
    const std::filesystem::path input_path(input_file);
    return (std::filesystem::path(output_dir) /
            (input_path.stem().string() + "_motif_hits.tsv")).string();
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <seq.fa|seq_dir> <pwm.txt> <thresholds.tsv> <out.tsv|out_dir>"
                  << " <core_len> [--batch]\n";
        return 1;
    }

    try {
        const std::string seq_path = argv[1];
        const std::string pwm_path = argv[2];
        const std::string threshold_path = argv[3];
        const std::string out_path = argv[4];
        const int core_len = std::stoi(argv[5]);

        bool batch_mode = std::filesystem::is_directory(seq_path);
        for (int i = 6; i < argc; ++i) {
            const std::string flag = argv[i];
            if (flag == "--batch") {
                batch_mode = true;
            } else if (flag == "--no-gaps" || flag.rfind("--allow-gaps=", 0) == 0) {
                std::cerr << "Ignoring deprecated gap flag: " << flag << '\n';
            } else {
                std::cerr << "Warning: unknown flag ignored: " << flag << '\n';
            }
        }

        const auto motifs = load_pwms(pwm_path);
        const auto thresholds = load_thresholds(threshold_path);
        if (motifs.empty()) throw std::runtime_error("No motifs found in PWM file");
        if (thresholds.empty()) throw std::runtime_error("No thresholds found");

        std::cerr << "Loaded " << motifs.size() << " motifs and "
                  << thresholds.size() << " threshold rows\n";

        const auto started = std::chrono::steady_clock::now();
        size_t total_hits = 0;

        if (batch_mode) {
            std::filesystem::create_directories(out_path);
            const auto fasta_files = list_fasta_files(seq_path);
            if (fasta_files.empty()) throw std::runtime_error("No FASTA files found in: " + seq_path);

            for (size_t i = 0; i < fasta_files.size(); ++i) {
                const auto sequences = load_fasta(fasta_files[i]);
                const auto hits = search_sequences(sequences, motifs, thresholds, core_len);
                const std::string result_file = batch_output_name(fasta_files[i], out_path);
                write_hits(result_file, hits);
                total_hits += hits.size();

                std::cerr << "[" << (i + 1) << "/" << fasta_files.size() << "] "
                          << std::filesystem::path(fasta_files[i]).filename().string()
                          << ": " << hits.size() << " hits -> "
                          << std::filesystem::path(result_file).filename().string() << '\n';
            }
        } else {
            const auto sequences = load_fasta(seq_path);
            if (sequences.empty()) throw std::runtime_error("No sequences found in: " + seq_path);

            const auto hits = search_sequences(sequences, motifs, thresholds, core_len);
            write_hits(out_path, hits);
            total_hits = hits.size();
            std::cerr << "Results: " << hits.size() << " hits -> " << out_path << '\n';
        }

        const auto finished = std::chrono::steady_clock::now();
        const auto ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(finished - started).count();
        std::cerr << "Search complete. Total hits: " << total_hits
                  << ". Elapsed: " << ms << " ms\n";
    } catch (const std::exception &error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
