#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <numeric>
#include <cmath>
#include <chrono>
#include <iomanip>

// Remove leading/trailing whitespace and carriage returns
static std::string clean_string(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// Position Weight Matrix structure
struct PositionWeightMatrix {
    std::string motif_id;
    std::vector<std::vector<double>> probability_matrix;
    
    void reserve_rows(size_t size) { probability_matrix.reserve(size); }
};

// Detect input file format
enum class FileFormat {
    MEME_FORMAT,      // Contains "letter-probability matrix:"
    SIMPLE_FORMAT     // Direct motif format with "Motif:" lines
};

FileFormat detect_file_format(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Failed to open file: " + filename);
    
    std::string line;
    while (std::getline(file, line)) {
        line = clean_string(line);
        if (line.find("letter-probability matrix:") != std::string::npos) {
            return FileFormat::MEME_FORMAT;
        }
        // If we find a Motif line without seeing letter-probability matrix first,
        // it's likely the simple format
        if (line.compare(0, 6, "Motif:") == 0) {
            return FileFormat::SIMPLE_FORMAT;
        }
    }
    
    // Default to simple format if no clear indicators found
    return FileFormat::SIMPLE_FORMAT;
}

// Parse MEME format file (original method)
std::vector<PositionWeightMatrix> extract_pwms_from_meme_format(const std::string& meme_file) {
    std::ifstream file(meme_file);
    if (!file) throw std::runtime_error("Failed to open MEME file: " + meme_file);
    
    std::vector<PositionWeightMatrix> pwm_collection;
    pwm_collection.reserve(100);  // Reserve space for typical motif count
    
    PositionWeightMatrix current_pwm;
    bool parsing_matrix = false;
    std::string line;
    
    // Precompiled regex for numerical matrix lines (4 floating point numbers)
    static const std::regex matrix_line_pattern(R"(^\s*[\d.eE+-]+(\s+[\d.eE+-]+){3}\s*$)");
    
    while (std::getline(file, line)) {
        // Clean line: remove carriage returns and trim whitespace
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty()) continue;
        
        line = clean_string(line);
        if (line.empty()) continue;
        
        // Detect new motif definition
        if (line.compare(0, 5, "MOTIF") == 0) {
            // Save previous motif if exists
            if (!current_pwm.motif_id.empty()) {
                pwm_collection.push_back(std::move(current_pwm));
            }
            
            // Initialize new motif
            current_pwm = PositionWeightMatrix();
            current_pwm.reserve_rows(30);  // Reserve for typical motif length
            parsing_matrix = false;
            
            // Extract motif ID (last word in MOTIF line)
            auto space_pos = line.find_last_of(' ');
            current_pwm.motif_id = (space_pos == std::string::npos) ? 
                                   line.substr(6) : line.substr(space_pos + 1);
        }
        // Detect start of probability matrix
        else if (line.find("letter-probability matrix:") != std::string::npos) {
            parsing_matrix = true;
        }
        // Parse matrix data lines
        else if (parsing_matrix && std::regex_match(line, matrix_line_pattern)) {
            std::istringstream stream(line);
            std::vector<double> probability_row;
            probability_row.reserve(4);  // A, C, G, T probabilities
            
            double prob_value;
            while (stream >> prob_value) {
                probability_row.push_back(prob_value);
            }
            
            if (probability_row.size() == 4) {
                current_pwm.probability_matrix.push_back(std::move(probability_row));
            }
        }
        // End of matrix when empty line encountered
        else if (parsing_matrix && line.empty()) {
            parsing_matrix = false;
        }
    }
    
    // Add final motif if exists
    if (!current_pwm.motif_id.empty()) {
        pwm_collection.push_back(std::move(current_pwm));
    }
    
    return pwm_collection;
}

// Parse simple format file (new method for direct motif format)
std::vector<PositionWeightMatrix> extract_pwms_from_simple_format(const std::string& input_file) {
    std::ifstream file(input_file);
    if (!file) throw std::runtime_error("Failed to open input file: " + input_file);
    
    std::vector<PositionWeightMatrix> pwm_collection;
    pwm_collection.reserve(100);
    
    PositionWeightMatrix current_pwm;
    std::string line;
    bool reading_matrix = false;
    
    // Precompiled regex for numerical matrix lines (4 floating point numbers)
    static const std::regex matrix_line_pattern(R"(^\s*[\d.eE+-]+(\s+[\d.eE+-]+){3}\s*$)");
    
    while (std::getline(file, line)) {
        // Clean line: remove carriage returns and trim whitespace
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line = clean_string(line);
        
        if (line.empty()) {
            // Empty line indicates end of current motif
            if (!current_pwm.motif_id.empty() && !current_pwm.probability_matrix.empty()) {
                pwm_collection.push_back(std::move(current_pwm));
                current_pwm = PositionWeightMatrix();
                reading_matrix = false;
            }
            continue;
        }
        
        // Detect motif header line
        if (line.compare(0, 6, "Motif:") == 0) {
            // Save previous motif if exists
            if (!current_pwm.motif_id.empty() && !current_pwm.probability_matrix.empty()) {
                pwm_collection.push_back(std::move(current_pwm));
            }
            
            // Initialize new motif
            current_pwm = PositionWeightMatrix();
            current_pwm.reserve_rows(30);
            
            // Extract motif ID (everything after "Motif:")
            std::string motif_part = line.substr(6);
            motif_part = clean_string(motif_part);
            current_pwm.motif_id = motif_part;
            reading_matrix = true;
        }
        // Parse matrix data lines
        else if (reading_matrix && std::regex_match(line, matrix_line_pattern)) {
            std::istringstream stream(line);
            std::vector<double> probability_row;
            probability_row.reserve(4);  // A, C, G, T probabilities
            
            double prob_value;
            while (stream >> prob_value) {
                probability_row.push_back(prob_value);
            }
            
            if (probability_row.size() == 4) {
                current_pwm.probability_matrix.push_back(std::move(probability_row));
            }
        }
    }
    
    // Add final motif if exists
    if (!current_pwm.motif_id.empty() && !current_pwm.probability_matrix.empty()) {
        pwm_collection.push_back(std::move(current_pwm));
    }
    
    return pwm_collection;
}

// Main extraction function that automatically detects format
std::vector<PositionWeightMatrix> extract_pwms_from_file(const std::string& input_file) {
    FileFormat format = detect_file_format(input_file);
    
    if (format == FileFormat::MEME_FORMAT) {
        std::cerr << "Detected MEME format, using original parsing method..." << std::endl;
        return extract_pwms_from_meme_format(input_file);
    } else {
        std::cerr << "Detected simple format, using direct parsing method..." << std::endl;
        return extract_pwms_from_simple_format(input_file);
    }
}

// Apply Bayesian smoothing or simple normalization to PWM probabilities
void normalize_probability_matrices(std::vector<PositionWeightMatrix>& pwm_list, 
                                   bool apply_bayesian_smoothing, 
                                   double pseudocount) {
    for (auto& pwm : pwm_list) {
        for (auto& position_probs : pwm.probability_matrix) {
            double row_sum = std::accumulate(position_probs.begin(), position_probs.end(), 0.0);
            
            if (row_sum > 0) {  // Prevent division by zero
                if (apply_bayesian_smoothing) {
                    double adjusted_sum = row_sum * pseudocount + 2.0;
                    for (auto& prob : position_probs) {
                        prob = (prob * pseudocount + 0.5) / adjusted_sum;
                    }
                } else {
                    // Simple normalization: divide by row sum
                    for (auto& prob : position_probs) {
                        prob /= row_sum;
                    }
                }
            }
        }
    }
}

// Write PWM data to output file with optimized I/O
void save_pwms_to_file(const std::vector<PositionWeightMatrix>& pwm_list, 
                       const std::string& output_file) {
    std::ofstream file(output_file);
    if (!file) throw std::runtime_error("Failed to create output file: " + output_file);
    
    // Set larger buffer for improved I/O performance
    file.rdbuf()->pubsetbuf(nullptr, 8192);
    
    for (const auto& pwm : pwm_list) {
        // Write motif header
        file << "Motif:" << pwm.motif_id << '\n';
        
        // Write probability matrix (A, C, G, T columns)
        for (const auto& position_probs : pwm.probability_matrix) {
            file << std::fixed << std::setprecision(6)
                 << position_probs[0] << '\t' << position_probs[1] << '\t'
                 << position_probs[2] << '\t' << position_probs[3] << '\n';
        }
        file << '\n';  // Blank line between motifs
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file> <output_pwm_file> <bayesian_smoothing(0|1)> <pseudocount>\n";
        std::cerr << "Input file can be either MEME format or simple motif format\n";
        return 1;
    }
    
    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::string input_file = argv[1];
        std::string output_pwm_file = argv[2];
        bool apply_bayesian_smoothing = std::stoi(argv[3]) != 0;
        double pseudocount_value = std::stod(argv[4]);
        
        // Extract PWMs using automatic format detection
        auto pwm_collection = extract_pwms_from_file(input_file);
        
        std::cerr << "Normalizing probability matrices..." << std::endl;
        normalize_probability_matrices(pwm_collection, apply_bayesian_smoothing, pseudocount_value);
        
        std::cerr << "Writing PWM file..." << std::endl;
        save_pwms_to_file(pwm_collection, output_pwm_file);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        std::cout << "Successfully converted " << pwm_collection.size() << " motifs" << std::endl;
        std::cerr << "Execution time: " << execution_time.count() << " ms" << std::endl;
        
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << std::endl;
        return 1;
    }
    
    return 0;
}