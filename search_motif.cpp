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

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_sequences.fasta> <input_pwms.txt> <input_thresholds.txt> "
                  << "<output_results.txt> <core_region_length>\n";
        return 1;
    }
    
    try {
        std::string input_sequence_file = argv[1];
        std::string input_pwm_file = argv[2];
        std::string input_threshold_file = argv[3];
        std::string output_results_file = argv[4];
        int core_region_length = std::stoi(argv[5]);
        
        std::cerr << "Loading input data..." << std::endl;
        
        // Load all required data
        auto sequence_collection = load_dna_sequences_from_fasta(input_sequence_file);
        auto pwm_collection = load_position_weight_matrices(input_pwm_file);
        auto threshold_map = load_detection_thresholds(input_threshold_file);
        
        // Validate input data
        if (sequence_collection.empty() || pwm_collection.empty() || threshold_map.empty()) {
            std::cerr << "ERROR: One or more input files contain no data" << std::endl;
            return 1;
        }
        
        std::cerr << "Starting parallel motif search..." << std::endl;
        
        // Execute parallel search
        SearchResultCollector result_collector;
        auto search_start_time = std::chrono::high_resolution_clock::now();
        
        process_sequences_in_parallel(sequence_collection, pwm_collection, threshold_map, 
                                     core_region_length, result_collector);
        
        auto search_end_time = std::chrono::high_resolution_clock::now();
        auto total_search_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            search_end_time - search_start_time);
        
        // Save results to file
        result_collector.write_results_to_file(output_results_file);
        
        std::cerr << "Motif search completed successfully!" << std::endl;
        std::cerr << "Execution time: " << total_search_time.count() << " ms" << std::endl;
        std::cerr << "Total matches found: " << result_collector.get_result_count() << std::endl;
        std::cerr << "Results saved to: " << output_results_file << std::endl;
        
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << std::endl;
        return 1;
    }
    
    return 0;
}