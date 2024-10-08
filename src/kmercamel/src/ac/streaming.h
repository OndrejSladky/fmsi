#pragma once
#include <unordered_set>
#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>

#include "kmers_ac.h"

void Push(std::ostream &of, const std::string &current, int k, int32_t used) {
    if (current.size() == static_cast<size_t>(k) && used != 0) {
        auto output = current.substr(1, k - 1 - __builtin_ctzll(used));
        std::transform(output.begin(), output.end(), output.begin(), ::tolower);
        of << output;
    }
}


void Streaming(std::string &path, std::ostream &of, int k, bool complements) {
    std::ifstream fasta(path);
    std::unordered_set <int64_t> kMers;
    std::string line;
    std::string kMer;
    int32_t used = 0;
    int32_t usedMask = (1 << (k - 1)) - 1;
    if (fasta.is_open()) {
        while (std::getline(fasta, line)) {
            if (!line.empty() && line[0] == '>') {
                Push(of, kMer, k, used);
                used = 0;
                kMer = "";
            } else {
                for (size_t i = 0; i < line.size(); ++i) {
                    if (NucleotideToInt(line[i]) == -1) {
                        Push(of, kMer, k, used);
                        used = 0;
                        kMer = "";
                    } else {
                        kMer += (char)::toupper(line[i]);
                        if (kMer.size() == size_t(k + 1)) kMer=kMer.substr(1);
                        assert(kMer.size() <= size_t(k));;
                        if (kMer.length() == size_t (k) ) {
                            uint64_t encoded = KMerToNumber(KMer{kMer});
                            auto rc = ReverseComplement(encoded, k);
                            auto rc2 = NumberToKMer(rc, k);
                            bool contained = (kMers.count(encoded) > 0) || ((complements) && kMers.count(ReverseComplement(encoded, k)) > 0);
                            if (!contained) kMers.insert(encoded);
                            if (!contained) of << kMer[0];
                            else if (used) {
                                of << (char)::tolower(kMer[0]);
                            }
                            used <<= 1;
                            used |= !contained;
                            used &= usedMask;
                        }
                    }
                }
            }
        }
        Push(of, kMer, k, used);
        fasta.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
}