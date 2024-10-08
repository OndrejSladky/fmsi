#pragma once

#include <string>
#include <iostream>
#include <glpk.h>

#include "parser.h"
#include "khash_utils.h"
#include "kmers.h"

/// For the given masked superstring output the same superstring with mask with minimal/maximal number of ones.
void OptimizeOnes(std::ifstream &in, std::ostream &of, kh_S64_t *kMers, int k, bool complements, bool minimize) {
    char c;
    int beforeKMerEnd = k;
    kmer_t currentKMer = 0;
    kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
    bool readingHeader = false;
    while (in >> std::noskipws >> c) {
        // Start of a header.
        if (c == '>') {
            readingHeader = true;
            currentKMer &= mask >> 1;
            for (char ch : NumberToKMer(currentKMer, k - beforeKMerEnd)) of << UpperToLower(ch);
            currentKMer = 0;
            beforeKMerEnd = k;
        }
        // Reprint the header.
        if (readingHeader) of << c;
        if (c == '\n') readingHeader = false;
        if (readingHeader) continue;
        auto data = NucleotideToInt(c);
        // Disregard white space.
        if (c == '\n' || c == '\r' || c == ' ') continue;
        // Break at an unknown character.
        if (data == -1) {
            currentKMer &= mask >> 1;
            for (char ch : NumberToKMer(currentKMer, k - beforeKMerEnd)) of << UpperToLower(ch);
            currentKMer = 0;
            beforeKMerEnd = k;
            continue;
        }
        currentKMer <<= 2;
        currentKMer &= mask;
        currentKMer |= data;
        --beforeKMerEnd;
        kmer_t reverseComplement = ReverseComplement(currentKMer, k);
        if (beforeKMerEnd == 0) {
            ++beforeKMerEnd;
            char to_print = NucleotideAtIndex(currentKMer, k, 0);
            auto key = kh_get_S64(kMers, currentKMer);
            if (key != kh_end(kMers)) {
                if (minimize) {
                    kh_del_S64(kMers, key);
                }
            } else if (complements && (key = kh_get_S64(kMers, reverseComplement)) != kh_end(kMers)) {
                if (minimize) {
                    kh_del_S64(kMers, key);
                }
            } else {
                to_print = UpperToLower(to_print);
            }
            of << to_print;
        }
    }
    currentKMer &= mask >> 1;
    for (char ch : NumberToKMer(currentKMer, k - beforeKMerEnd)) of << UpperToLower(ch);
    of << std::endl;
}

/// Pre-solve to the ILP in minimizing the number of runs.
/// Set the intervals with single-occurring k-mers to 1 and remove k-mers from those intervals.
std::pair<std::vector<int>, std::vector<int>> HeuristicPreSolve(std::vector<std::list<size_t>> &intervalsForKMer, size_t intervalCount, int &mappedSize, size_t &totalIntervals, int &newIntervals) {
    std::vector<int> mapping(intervalsForKMer.size());
    std::vector<int> intervalMapping(intervalCount);

    // Set undecided intervals with a single-occurring k-mer to 1s.
    for (size_t i = 0; i < intervalsForKMer.size(); ++i) {
        if (intervalsForKMer[i].size() == 1) {
            mapping[i] = -1;
            intervalMapping[intervalsForKMer[i].front()] = -1;
        }
    }

    // Remove k-mers from decided intervals.
    for (size_t i = 0; i < intervalsForKMer.size(); ++i) {
        for (auto j : intervalsForKMer[i]) {
            if (intervalMapping[j] == -1) {
                mapping[i] = -1;
            }
        }
    }

    totalIntervals = 0;
    int newIndex = 0;
    for (size_t i = 0; i < intervalsForKMer.size(); ++i) {
        if (mapping[i] == -1) continue;
        mapping[i] = newIndex++;
        totalIntervals += intervalsForKMer[i].size();
    }
    mappedSize = newIndex;

    int newIntervalIndex = 0;
    for (size_t i = 0; i < intervalCount; ++i) {
        if (intervalMapping[i] == -1) continue;
        intervalMapping[i] = newIntervalIndex++;
    }
    newIntervals = newIntervalIndex;


    return {mapping, intervalMapping};
}


/// For the given masked superstring output the same superstring with mask with minimal number of runs of ones.
void OptimizeRuns(std::string path, kh_S64_t *kMers, std::ostream &of, int k, bool complements, bool approximate) {
    kh_O64_t *intervals = kh_init_O64();
    std::vector<std::list<size_t>> intervalsForKMer;
    auto [size, rows] = ReadIntervals(intervals, kMers, intervalsForKMer, path, k, complements, of, nullptr);
    int mappedSize, newIntervals; size_t totalIntervals;
    auto [mapping, intervalMapping] = HeuristicPreSolve(intervalsForKMer, rows, mappedSize, totalIntervals, newIntervals);
    glp_prob *lp;
    lp = glp_create_prob();
    if (mappedSize != 0 && !approximate) {
        auto *ia = new int[totalIntervals + 1];
        auto *ja = new int[totalIntervals + 1];
        auto *ar = new double[totalIntervals + 1];
        glp_set_obj_dir(lp, GLP_MIN);
        // Add a row per each k-mer.
        glp_add_rows(lp, mappedSize);
        for (int i = 0; i < mappedSize; ++i) {
            glp_set_row_bnds(lp, i + 1, GLP_LO, 1.0, 0.0);
        }
        // Add a column per each undecided interval.
        glp_add_cols(lp, newIntervals);
        for (int i = 0; i < newIntervals; ++i) {
            glp_set_col_bnds(lp, i + 1, GLP_LO, 0.0, 1.0);
            glp_set_col_kind(lp, i + 1, GLP_IV);
            glp_set_obj_coef(lp, i + 1, 1.0);
        }
        int index = 0;
        for (auto key = kh_begin(intervals); key != kh_end(intervals); ++key) {
            if (!kh_exist(intervals, key)) continue;
            size_t i = kh_value(intervals, key);
            if (mapping[i] == -1) continue;
            for (auto j: intervalsForKMer[i]) {
                ja[index + 1] = intervalMapping[j] + 1;
                ia[index + 1] = mapping[i] + 1;
                ar[index + 1] = 1.0;
                ++index;
            }
        }
        // Supress glpk output.
        glp_term_out(GLP_OFF);

        glp_load_matrix(lp, index, ia, ja, ar);
        glp_simplex(lp, nullptr);
    }

    bool *intervalsSet = new bool[rows];

    for (size_t i = 0; i < rows; ++i) {
        if (intervalMapping[i] == -1) intervalsSet[i] = true;
        else if (approximate) intervalsSet[i] = mappedSize != 0;
        else intervalsSet[i] = mappedSize == 0 ? false : (glp_get_col_prim(lp, intervalMapping[i] + 1) > 0.5);
    }

    ReadIntervals(nullptr, kMers, intervalsForKMer, path, k, complements, of, intervalsSet);
    of << std::endl;
}

int Optimize(std::string &algorithm, std::string path, std::ostream &of,  int k, bool complements) {
    kh_S64_t *kMers = kh_init_S64();
    ReadKMers(kMers, path, k, complements, true);

    std::ifstream in(path);
    if (in.is_open()) {
        if (algorithm == "ones") {
            OptimizeOnes(in, of, kMers, k, complements, false);
        } else if (algorithm == "zeros") {
            OptimizeOnes(in, of, kMers, k, complements, true);
        } else if (algorithm == "runs") {
            OptimizeRuns(path, kMers, of, k, complements, false);
        } else if (algorithm == "runsapprox") {
            OptimizeRuns(path, kMers, of, k, complements, true);
        } else {
            std::cerr << "Algorithm '" + algorithm + "' not recognized." << std::endl;
            in.close();
            return 1;
        }
        in.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
    return 0;
}