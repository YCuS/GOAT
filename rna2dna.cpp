#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <sequence_file>\n";
        return 1;
    }
    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Cannot open file: " << argv[1] << "\n";
        return 1;
    }

    std::string line;
    while (std::getline(infile, line)) {
        for (char& c : line) {
            if (c == 'U') c = 'T';
        }
        std::cout << line << "\n";
    }

    return 0;
}