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

// Position Weight Matrix structure for threshold calculation
struct PositionWeightMatrix {
    std::map<char, std::vector<double>> nucleotide_probabilities;
    int get_length() const { return nucleotide_probabilities.at('A').size(); }
};

// Results container for calculated thresholds
struct MotifThresholds {
    std::string motif_identifier;
    double full_motif_threshold;
    double core_region_threshold;
    int core_start_position;
    bool calculation_successful = false;
};

// Parse PWM file and extract probability matrices with memory optimization
std::map<std::string, std::vector<std::vector<double>>> 
load_probability_matrices(const std::string& pwm_file) {
    std::ifstream file(pwm_file);
    if (!file) throw std::runtime_error("Failed to open PWM file: " + pwm_file);
    
    std::map<std::string, std::vector<std::vector<double>>> motif_matrices;
    std::string line, current_motif_id;
    std::vector<std::vector<double>> current_matrix;
    current_matrix.reserve(30);  // Reserve for typical motif length
    
    // Pattern to match probability matrix lines (4 decimal numbers)
    std::regex probability_line_pattern(R"(^\s*[\d.eE+-]+(\s+[\d.eE+-]+){3}\s*$)");
    
    while (std::getline(file, line)) {
        // Remove carriage returns
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        
        // Detect motif header
        if (line.compare(0, 6, "Motif:") == 0) {
            // Store previous motif if complete
            if (!current_motif_id.empty() && !current_matrix.empty()) {
                motif_matrices[current_motif_id] = std::move(current_matrix);
            }
            
            // Initialize new motif
            current_motif_id = line.substr(6);
            current_matrix.clear();
            current_matrix.reserve(30);
        }
        // Parse probability values
        else if (std::regex_match(line, probability_line_pattern)) {
            std::istringstream stream(line);
            std::vector<double> position_probabilities;
            position_probabilities.reserve(4);
            
            double probability;
            while (stream >> probability) {
                position_probabilities.push_back(probability);
            }
            
            if (position_probabilities.size() == 4) {
                current_matrix.push_back(std::move(position_probabilities));
            }
        }
    }
    
    // Store final motif
    if (!current_motif_id.empty() && !current_matrix.empty()) {
        motif_matrices[current_motif_id] = std::move(current_matrix);
    }
    
    return motif_matrices;
}

// Convert matrix format to PWM structure for calculations
PositionWeightMatrix create_pwm_from_matrix(const std::vector<std::vector<double>>& probability_matrix) {
    PositionWeightMatrix pwm;
    
    // Pre-allocate vectors for nucleotides
    for (char nucleotide : {'A', 'C', 'G', 'T'}) {
        pwm.nucleotide_probabilities[nucleotide].reserve(probability_matrix.size());
    }
    
    // Fill PWM structure: rows become columns
    for (const auto& position_probs : probability_matrix) {
        pwm.nucleotide_probabilities['A'].push_back(position_probs[0]);
        pwm.nucleotide_probabilities['C'].push_back(position_probs[1]);
        pwm.nucleotide_probabilities['G'].push_back(position_probs[2]);
        pwm.nucleotide_probabilities['T'].push_back(position_probs[3]);
    }
    return pwm;
}

// Generate random sequence according to PWM probabilities
std::string generate_sequence_from_pwm(const PositionWeightMatrix& pwm, 
                                      int sequence_length, 
                                      std::mt19937& random_generator) {
    const std::string nucleotides = "ACGT";
    std::string generated_sequence;
    generated_sequence.reserve(sequence_length);
    
    for (int position = 0; position < sequence_length; ++position) {
        std::vector<double> position_probabilities;
        position_probabilities.reserve(4);
        double probability_sum = 0;
        
        // Collect probabilities for this position
        for (char nucleotide : nucleotides) {
            position_probabilities.push_back(pwm.nucleotide_probabilities.at(nucleotide)[position]);
            probability_sum += position_probabilities.back();
        }
        
        // Generate nucleotide based on probabilities
        if (probability_sum > 0) {
            // Normalize probabilities
            for (auto& prob : position_probabilities) prob /= probability_sum;
            
            std::discrete_distribution<> nucleotide_distribution(
                position_probabilities.begin(), position_probabilities.end());
            generated_sequence += nucleotides[nucleotide_distribution(random_generator)];
        } else {
            // Fallback for zero probabilities
            generated_sequence += nucleotides[random_generator() % 4];
        }
    }
    return generated_sequence;
}

// Calculate log-likelihood score for sequence against PWM
double calculate_log_likelihood_score(const std::string& sequence, const PositionWeightMatrix& pwm) {
    double total_score = 0;
    const double zero_probability_penalty = -std::log(0.000001);  // Penalty for zero probabilities
    
    for (size_t position = 0; position < sequence.size(); ++position) {
        double nucleotide_probability = pwm.nucleotide_probabilities.at(sequence[position])[position];
        total_score += (nucleotide_probability > 0) ? 
                       -std::log(nucleotide_probability) : zero_probability_penalty;
    }
    return total_score;
}

// Calculate thresholds for a single motif using parallel Monte Carlo simulation
MotifThresholds calculate_motif_thresholds_parallel(const std::string& motif_id,
                                                   const PositionWeightMatrix& full_pwm,
                                                   double full_threshold_percentile, 
                                                   double core_threshold_percentile,
                                                   int simulation_iterations, 
                                                   int core_region_length) {
    MotifThresholds result{motif_id, 0, 0, 0, false};
    int motif_length = full_pwm.get_length();
    
    if (motif_length < core_region_length) return result;
    
    // Use thread-local random generator for thread safety
    thread_local std::mt19937 random_generator(std::random_device{}());
    
    std::vector<double> full_motif_scores, core_region_scores;
    full_motif_scores.reserve(simulation_iterations);
    core_region_scores.reserve(simulation_iterations);
    
    int core_start_position = 0;  // Start core region from beginning
    
    // Create core region PWM
    PositionWeightMatrix core_pwm;
    for (char nucleotide : {'A', 'C', 'G', 'T'}) {
        core_pwm.nucleotide_probabilities[nucleotide].reserve(core_region_length);
        for (int pos = core_start_position; 
             pos < core_start_position + core_region_length && pos < motif_length; ++pos) {
            core_pwm.nucleotide_probabilities[nucleotide].push_back(
                full_pwm.nucleotide_probabilities.at(nucleotide)[pos]);
        }
    }
    
    // Monte Carlo simulation
    for (int iteration = 0; iteration < simulation_iterations; ++iteration) {
        // Generate sequence from full PWM
        auto simulated_sequence = generate_sequence_from_pwm(full_pwm, motif_length, random_generator);
        full_motif_scores.push_back(calculate_log_likelihood_score(simulated_sequence, full_pwm));
        
        // Score core region
        if (simulated_sequence.size() >= core_region_length) {
            auto core_sequence = simulated_sequence.substr(core_start_position, core_region_length);
            core_region_scores.push_back(calculate_log_likelihood_score(core_sequence, core_pwm));
        }
    }
    
    // Calculate threshold values from percentiles
    std::sort(full_motif_scores.begin(), full_motif_scores.end());
    std::sort(core_region_scores.begin(), core_region_scores.end());
    
    int full_threshold_index = std::max(0, std::min(
        static_cast<int>(full_motif_scores.size() * full_threshold_percentile / 100.0) - 1, 
        static_cast<int>(full_motif_scores.size()) - 1));
    int core_threshold_index = std::max(0, std::min(
        static_cast<int>(core_region_scores.size() * core_threshold_percentile / 100.0) - 1, 
        static_cast<int>(core_region_scores.size()) - 1));
    
    result.full_motif_threshold = full_motif_scores[full_threshold_index];
    result.core_region_threshold = core_region_scores[core_threshold_index];
    result.core_start_position = core_start_position;
    result.calculation_successful = true;
    
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_pwm_file> <output_threshold_file> <full_threshold_percentile> "
                  << "<core_threshold_percentile> <simulation_iterations> <core_length> <output_precision>\n";
        return 1;
    }
    
    try {
        std::string input_pwm_file = argv[1];
        std::string output_threshold_file = argv[2];
        double full_threshold_percentile = std::stod(argv[3]);
        double core_threshold_percentile = std::stod(argv[4]);
        int simulation_iterations = std::stoi(argv[5]);
        int core_region_length = std::stoi(argv[6]);
        int output_precision = std::stoi(argv[7]);
        
        // std::cerr << "Loading PWM probability matrices..." << std::endl;
        auto motif_matrices = load_probability_matrices(input_pwm_file);
        std::cerr << "Found " << motif_matrices.size() << " motifs for processing" << std::endl;
        
        // std::cerr << "Starting parallel threshold calculations..." << std::endl;
        auto calculation_start_time = std::chrono::high_resolution_clock::now();
        
        // Launch parallel calculations for all motifs
        std::vector<std::future<MotifThresholds>> calculation_futures;
        for (const auto& [motif_id, probability_matrix] : motif_matrices) {
            auto pwm = create_pwm_from_matrix(probability_matrix);
            calculation_futures.push_back(std::async(std::launch::async, 
                                                    calculate_motif_thresholds_parallel, 
                                                    motif_id, pwm, 
                                                    full_threshold_percentile, core_threshold_percentile, 
                                                    simulation_iterations, core_region_length));
        }
        
        // Collect results and write to output file
        std::ofstream output_file(output_threshold_file);
        if (!output_file) {
            throw std::runtime_error("Failed to create output file: " + output_threshold_file);
        }
        
        int completed_calculations = 0;
        for (auto& future : calculation_futures) {
            auto threshold_result = future.get();
            if (threshold_result.calculation_successful) {
                output_file << threshold_result.motif_identifier << '\t'
                           << std::fixed << std::setprecision(output_precision) 
                           << threshold_result.full_motif_threshold << '\t'
                           << threshold_result.core_region_threshold << '\t'
                           << threshold_result.core_start_position << '\n';
            }
            
            completed_calculations++;
            // if (completed_calculations % 10 == 0) {
            //     std::cerr << "Completed " << completed_calculations << "/" 
            //              << motif_matrices.size() << " motif calculations" << std::endl;
            // }
        }
        
        auto calculation_end_time = std::chrono::high_resolution_clock::now();
        auto total_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            calculation_end_time - calculation_start_time);
        
        std::cerr << "Threshold calculation completed successfully!" << std::endl;
        std::cerr << "Total execution time: " << total_execution_time.count() << " ms" << std::endl;
        
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << std::endl;
        return 1;
    }
    
    return 0;
}