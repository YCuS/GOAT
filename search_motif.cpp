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

// Position Weight Matrix with cache-optimized storage
struct PositionWeightMatrix {
    std::string motif_identifier;
    std::vector<std::array<double, 4>> position_probabilities;  // [A,C,G,T] for each position
    
    void reserve_positions(size_t position_count) { position_probabilities.reserve(position_count); }
    size_t get_length() const { return position_probabilities.size(); }
};

// Threshold parameters for motif detection
struct MotifDetectionThresholds {
    std::string motif_identifier;
    double full_motif_threshold;
    double core_region_threshold;
    int core_start_position;
};

// DNA sequence with precomputed reverse complement
struct DNASequence {
    std::string sequence_header;
    std::string forward_sequence;
    std::string reverse_complement;
    
    DNASequence(std::string header, std::string sequence) 
        : sequence_header(std::move(header)), forward_sequence(std::move(sequence)) {
        reverse_complement = compute_reverse_complement(forward_sequence);
    }
    
private:
    // Compute reverse complement using lookup table for performance
    static std::string compute_reverse_complement(const std::string& sequence) {
        std::string complement;
        complement.reserve(sequence.size());
        
        // Initialize complement lookup table once
        static char complement_lookup[256];
        static bool lookup_initialized = false;
        
        if (!lookup_initialized) {
            std::fill(complement_lookup, complement_lookup + 256, 0);
            complement_lookup['A'] = 'T'; complement_lookup['T'] = 'A';
            complement_lookup['C'] = 'G'; complement_lookup['G'] = 'C';
            complement_lookup['a'] = 't'; complement_lookup['t'] = 'a';
            complement_lookup['c'] = 'g'; complement_lookup['g'] = 'c';
            lookup_initialized = true;
        }
        
        // Build reverse complement
        for (auto nucleotide_iter = sequence.rbegin(); nucleotide_iter != sequence.rend(); ++nucleotide_iter) {
            char complement_base = complement_lookup[static_cast<unsigned char>(*nucleotide_iter)];
            complement += (complement_base != 0) ? complement_base : *nucleotide_iter;
        }
        return complement;
    }
};

// Motif hit result: sequence_name, motif_id, position, strand, score, matched_sequence
using MotifSearchResult = std::tuple<std::string, std::string, int, std::string, double, std::string>;

// Thread-safe collector for search results
class SearchResultCollector {
private:
    std::vector<MotifSearchResult> collected_results;
    std::mutex collection_mutex;
    
public:
    void add_search_results(std::vector<MotifSearchResult>&& new_results) {
        std::lock_guard<std::mutex> lock(collection_mutex);
        collected_results.insert(collected_results.end(), 
                                std::make_move_iterator(new_results.begin()),
                                std::make_move_iterator(new_results.end()));
    }
    
    void write_results_to_file(const std::string& output_filename) {
        std::ofstream output_file(output_filename);
        if (!output_file) {
            throw std::runtime_error("Failed to create output file: " + output_filename);
        }
        
        for (const auto& result : collected_results) {
            output_file << std::get<0>(result) << '\t'    // sequence_name
                       << std::get<1>(result) << '\t'    // motif_id
                       << std::get<2>(result) << '\t'    // position
                       << std::get<3>(result) << '\t'    // strand
                       << std::get<4>(result) << '\t'    // score
                       << std::get<5>(result) << '\n';   // matched_sequence
        }
    }
    
    size_t get_result_count() const { return collected_results.size(); }
};

// Load PWM definitions from file with memory preallocation
std::vector<PositionWeightMatrix> load_position_weight_matrices(const std::string& pwm_file) {
    std::ifstream file(pwm_file);
    if (!file) throw std::runtime_error("Failed to open PWM file: " + pwm_file);
    
    std::vector<PositionWeightMatrix> pwm_collection;
    pwm_collection.reserve(100);  // Reserve space for typical motif count
    
    PositionWeightMatrix current_pwm;
    std::string line;
    
    while (std::getline(file, line)) {
        // Clean input line
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;
        
        // Remove leading/trailing whitespace
        auto content_start = line.find_first_not_of(" \t");
        auto content_end = line.find_last_not_of(" \t");
        if (content_start == std::string::npos) continue;
        line = line.substr(content_start, content_end - content_start + 1);
        
        // Process motif header
        if (line.compare(0, 6, "Motif:") == 0) {
            // Save previous motif if complete
            if (!current_pwm.position_probabilities.empty()) {
                pwm_collection.push_back(std::move(current_pwm));
            }
            
            // Initialize new motif
            current_pwm = PositionWeightMatrix();
            current_pwm.motif_identifier = line.substr(6);
            current_pwm.reserve_positions(30);  // Reserve for typical motif length
        }
        // Process probability data
        else {
            std::istringstream stream(line);
            double prob_A, prob_C, prob_G, prob_T;
            if (stream >> prob_A >> prob_C >> prob_G >> prob_T) {
                current_pwm.position_probabilities.push_back({prob_A, prob_C, prob_G, prob_T});
            }
        }
    }
    
    // Add final motif
    if (!current_pwm.position_probabilities.empty()) {
        pwm_collection.push_back(std::move(current_pwm));
    }
    
    std::cerr << "Loaded " << pwm_collection.size() << " PWMs" << std::endl;
    return pwm_collection;
}

// Load detection thresholds from file
std::map<std::string, MotifDetectionThresholds> load_detection_thresholds(const std::string& threshold_file) {
    std::ifstream file(threshold_file);
    if (!file) throw std::runtime_error("Failed to open threshold file: " + threshold_file);
    
    std::map<std::string, MotifDetectionThresholds> threshold_map;
    std::string line;
    
    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;
        
        std::istringstream stream(line);
        MotifDetectionThresholds thresholds;
        if (stream >> thresholds.motif_identifier >> thresholds.full_motif_threshold >> 
                      thresholds.core_region_threshold >> thresholds.core_start_position) {
            threshold_map[thresholds.motif_identifier] = thresholds;
        }
    }
    
    std::cerr << "Loaded " << threshold_map.size() << " threshold sets" << std::endl;
    return threshold_map;
}

// Load DNA sequences from FASTA file
std::vector<DNASequence> load_dna_sequences_from_fasta(const std::string& fasta_file) {
    std::ifstream file(fasta_file);
    if (!file) throw std::runtime_error("Failed to open sequence file: " + fasta_file);
    
    std::vector<DNASequence> sequence_collection;
    std::string line, current_header, current_sequence;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Save previous sequence if complete
            if (!current_header.empty() && !current_sequence.empty()) {
                sequence_collection.emplace_back(std::move(current_header), std::move(current_sequence));
            }
            
            // Start new sequence
            current_header = line.substr(1);
            current_sequence.clear();
        } else {
            current_sequence += line;
        }
    }
    
    // Add final sequence
    if (!current_header.empty() && !current_sequence.empty()) {
        sequence_collection.emplace_back(std::move(current_header), std::move(current_sequence));
    }
    
    std::cerr << "Loaded " << sequence_collection.size() << " DNA sequences" << std::endl;
    return sequence_collection;
}

// Calculate log-likelihood score using optimized lookup table
double calculate_motif_score_optimized(const std::string& sequence, 
                                      int start_position, 
                                      const PositionWeightMatrix& pwm) {
    double total_score = 0.0;
    // Nucleotide to array index mapping: A=0, C=1, G=2, T=3
    static int nucleotide_to_index[256];
    static bool lookup_initialized = false;
    
    if (!lookup_initialized) {
        std::fill(nucleotide_to_index, nucleotide_to_index + 256, -1);
        nucleotide_to_index['A'] = 0; nucleotide_to_index['a'] = 0;
        nucleotide_to_index['C'] = 1; nucleotide_to_index['c'] = 1;
        nucleotide_to_index['G'] = 2; nucleotide_to_index['g'] = 2;
        nucleotide_to_index['T'] = 3; nucleotide_to_index['t'] = 3;
        lookup_initialized = true;
    }
    
    // Calculate log-likelihood score for each position
    for (size_t position = 0; position < pwm.position_probabilities.size(); ++position) {
        int nucleotide_index = nucleotide_to_index[static_cast<unsigned char>(sequence[start_position + position])];
        if (nucleotide_index >= 0) {
            double probability = pwm.position_probabilities[position][nucleotide_index];
            if (probability > 0.0) {
                total_score += -std::log(probability);
            } else {
                std::cerr << "Warning: PWM for motif '"
                          << pwm.motif_identifier
                          << "' contains a zero probability at position "
                          << position
                          << ". Consider using Bayesian mode for more robust scoring."
                          << std::endl;
                total_score += -std::log(0.000001);  // Use a small fallback value
            }
        }
    }
    
    return total_score;
}

// Search for motif matches in a single sequence using sliding window
std::vector<MotifSearchResult> search_motif_in_sequence(
    const DNASequence& target_sequence,
    const PositionWeightMatrix& search_pwm,
    const MotifDetectionThresholds& detection_thresholds,
    int core_region_length) {
    
    std::vector<MotifSearchResult> sequence_matches;
    sequence_matches.reserve(100);  // Reserve space for typical match count
    
    const auto& forward_seq = target_sequence.forward_sequence;
    const auto& reverse_seq = target_sequence.reverse_complement;
    const auto& seq_name = target_sequence.sequence_header;  // Get sequence name
    
    // Create core region PWM for initial filtering
    PositionWeightMatrix core_pwm;
    core_pwm.motif_identifier = search_pwm.motif_identifier;
    int core_start = detection_thresholds.core_start_position;
    int core_end = std::min(core_start + core_region_length, static_cast<int>(search_pwm.get_length()));
    
    for (int position = core_start; position < core_end; ++position) {
        core_pwm.position_probabilities.push_back(search_pwm.position_probabilities[position]);
    }
    
    // Search function for both strands
    auto search_strand = [&](const std::string& strand, const std::string& strand_name) {
        if (strand.size() < search_pwm.get_length()) return;
        
        size_t max_search_position = strand.size() - search_pwm.get_length() + 1;
        for (size_t pos = 0; pos < max_search_position; ++pos) {
            // First check core region for efficiency
            double core_score = calculate_motif_score_optimized(strand, pos + core_start, core_pwm);
            if (core_score < detection_thresholds.core_region_threshold) {
                // Core passed, check full motif
                double full_score = calculate_motif_score_optimized(strand, pos, search_pwm);
                if (full_score < detection_thresholds.full_motif_threshold) {
                    std::string matched_sequence = strand.substr(pos, search_pwm.get_length());
                    sequence_matches.emplace_back(seq_name,                             // sequence_name
                                                 detection_thresholds.motif_identifier,  // motif_id
                                                 static_cast<int>(pos),                  // position
                                                 strand_name,                            // strand
                                                 full_score,                             // score
                                                 std::move(matched_sequence));           // matched_sequence
                }
            }
        }
    };
    
    // Search both strands
    search_strand(forward_seq, "positive");
    search_strand(reverse_seq, "negative");
    
    return sequence_matches;
}

// Process multiple sequences in parallel using thread pool
void process_sequences_in_parallel(
    const std::vector<DNASequence>& sequence_collection,
    const std::vector<PositionWeightMatrix>& pwm_collection,
    const std::map<std::string, MotifDetectionThresholds>& threshold_map,
    int core_region_length,
    SearchResultCollector& result_collector) {
    
    const int thread_count = std::thread::hardware_concurrency();
    std::vector<std::future<void>> processing_futures;
    
    // Process sequences in chunks across threads
    auto process_sequence_chunk = [&](size_t chunk_start, size_t chunk_end) {
        for (size_t seq_index = chunk_start; seq_index < chunk_end; ++seq_index) {
            const auto& sequence = sequence_collection[seq_index];
            
            // Search each motif in current sequence
            for (const auto& pwm : pwm_collection) {
                auto threshold_iter = threshold_map.find(pwm.motif_identifier);
                if (threshold_iter != threshold_map.end()) {
                    auto matches = search_motif_in_sequence(sequence, pwm, threshold_iter->second, core_region_length);
                    if (!matches.empty()) {
                        result_collector.add_search_results(std::move(matches));
                    }
                }
            }
        }
    };
    
    // Distribute work across available threads
    size_t sequences_per_thread = (sequence_collection.size() + thread_count - 1) / thread_count;
    for (int thread_id = 0; thread_id < thread_count; ++thread_id) {
        size_t chunk_start = thread_id * sequences_per_thread;
        size_t chunk_end = std::min(chunk_start + sequences_per_thread, sequence_collection.size());
        if (chunk_start < chunk_end) {
            processing_futures.push_back(std::async(std::launch::async, 
                                                   process_sequence_chunk, 
                                                   chunk_start, chunk_end));
        }
    }
    
    // Wait for all threads to complete
    for (auto& future : processing_futures) {
        future.wait();
    }
}


static std::vector<std::string> get_sequence_files_from_directory(const std::string& dir_path) {
    if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path)) {
        throw std::runtime_error("Sequence directory does not exist or is not a directory: " + dir_path);
    }
    std::vector<std::string> files;
    for (auto& entry : std::filesystem::directory_iterator(dir_path)) {
        if (!entry.is_regular_file()) continue;
        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".fa" || ext == ".fasta" || ext == ".fna") {
            files.push_back(entry.path().string());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

static std::string generate_output_filename(const std::string& input_file, const std::string& output_dir) {
    std::filesystem::path in(input_file);
    auto stem = in.stem().string();
    std::filesystem::path out_dir(output_dir);
    // output file format: <stem>_results.txt
    return (out_dir / (stem + "_results.txt")).string();
}

int main(int argc, char* argv[]) {
    // support single-file mode and batch mode
    if (argc != 6 && argc != 7) {
        std::cerr << "Usage (single file mode):\n  " << argv[0]
                  << " <input_sequences.fasta> <input_pwms.txt> <input_thresholds.txt> "
                  << "<output_results.txt> <core_region_length>\n";
        std::cerr << "Usage (batch mode):\n  " << argv[0]
                  << " <input_seq_directory> <input_pwms.txt> <input_thresholds.txt> "
                  << "<output_results_directory> <core_region_length> --batch\n";
        return 1;
    }

    try {
        std::string seq_path       = argv[1];
        std::string pwm_file       = argv[2];
        std::string thr_file       = argv[3];
        std::string out_path       = argv[4];
        int core_region_length     = std::stoi(argv[5]);
        bool batch_mode = (argc == 7 && std::string(argv[6]) == "--batch")
                       || std::filesystem::is_directory(seq_path);

        // Preload PWM definitions and thresholds; no need to reload for each file in batch mode
        auto pwm_collection = load_position_weight_matrices(pwm_file);
        auto threshold_map  = load_detection_thresholds(thr_file);
        if (pwm_collection.empty() || threshold_map.empty()) {
            throw std::runtime_error("Failed to load PWM definitions or thresholds; they may be empty");
        }

        if (batch_mode) {
            // Batch mode: seq_path is a directory, out_path is created if it does not exist
            if (!std::filesystem::exists(out_path)) {
                std::filesystem::create_directories(out_path);
            }

            auto seq_files = get_sequence_files_from_directory(seq_path);
            if (seq_files.empty()) {
                std::cerr << "No FASTA files found in directory: " << seq_path << "\n";
                return 1;
            }
            std::cerr << "Batch mode: found " << seq_files.size() << " FASTA files\n";

            int idx = 0;
            for (const auto& seq_file : seq_files) {
                ++idx;
                std::string result_file = generate_output_filename(seq_file, out_path);
                std::cerr << "Processing [" << idx << "/" << seq_files.size() << "]: "
                          << std::filesystem::path(seq_file).filename().string()
                          << " -> " << std::filesystem::path(result_file).filename().string()
                          << "\n";

                // Load sequences from the current FASTA file
                auto sequences = load_dna_sequences_from_fasta(seq_file);
                if (sequences.empty()) {
                    std::cerr << "  ⚠️  No sequences found in " << seq_file << ", skipping\n";
                    continue;
                }

                // Perform parallel motif search
                SearchResultCollector collector;
                auto t0 = std::chrono::high_resolution_clock::now();
                process_sequences_in_parallel(sequences, pwm_collection, threshold_map,
                                             core_region_length, collector);
                auto t1 = std::chrono::high_resolution_clock::now();

                // Write results to file
                collector.write_results_to_file(result_file);
                auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
                std::cerr << "  ✓ Done, matched " << collector.get_result_count()
                          << " entries, took " << ms << " ms\n";
            }

        } else {
            auto sequences = load_dna_sequences_from_fasta(seq_path);
            if (sequences.empty()) {
                throw std::runtime_error("Can't find sequences in " + seq_path);
            }
            std::cerr << "Loaded " << sequences.size() << " sequences\n";

            SearchResultCollector collector;
            auto t0 = std::chrono::high_resolution_clock::now();
            process_sequences_in_parallel(sequences, pwm_collection, threshold_map,
                                         core_region_length, collector);
            auto t1 = std::chrono::high_resolution_clock::now();

            collector.write_results_to_file(out_path);
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
            std::cerr << "Motif search completed successfully!\n";
            std::cerr << "Execution time: " << ms << " ms\n";
            std::cerr << "Total matches found: " << collector.get_result_count() << "\n";
            std::cerr << "Results saved to: " << out_path << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}