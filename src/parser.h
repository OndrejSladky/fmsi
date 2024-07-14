#pragma once

#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

/// Load masked superstring from the fasta file in the mask-cased format.
std::string read_masked_superstring(std::string fn) {
    std::vector<char> masked_superstring;
    std::ifstream fasta(fn);
    if (fasta.is_open()) {
        char c;
        bool readingHeader = false;
        while (fasta >> std::noskipws >> c) {
            if (c == '>') {
                readingHeader = true;
            }
            else if (c == '\n') readingHeader = false;
            if (readingHeader) continue;
            // Disregard white space.
            if (c == '\n' || c == '\r' || c == ' ') continue;
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'a' && c != 'c' && c != 'g' && c != 't') {
                throw std::invalid_argument("invalid character in masked superstring");
            }
            masked_superstring.push_back(c);
        }
        fasta.close();
    } else {
        throw std::invalid_argument("couldn't open file " + fn);
    }
    return std::string(masked_superstring.begin(), masked_superstring.end());
}
