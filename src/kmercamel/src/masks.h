#pragma once

#include <string>
#include <iostream>
#include <glpk.h>

#include "parser.h"
#include "khash_utils.h"
#include "kmers.h"

/// Print a warning to stderr if the mask after optimization violates the mask convention
/// (i.e., if there are more than k-1 last characters OFF).
void PrintMaskConventionWarning() {
    std::cerr << "Warning: the mask after optimization violates the mask convention" <<
        "as there are more than k-1 trailing zeros (more than k-1 trailing lowercase characters)." << std::endl <<
        "Depending on your case, you might want to consider one of the following manual fixes:" << std::endl <<
        "  - Set the k-th last mask symbol to 1 (corresponds to setting the letter to uppercase)." <<
        "The masked superstring will still represent the same set, but the mask will not be optimal" << std::endl <<
        "  - Truncate the masked superstring so that there are exactly k-1 trailing lowercase characters." <<
        "The masked superstring will still represent the same set, but the superstring will differ from the original" <<
        " (it will be shorter)." << std::endl <<
        "  - Ignore the warning if the previous fixes are not suitable" <<
        "Ensure that for your masked superstring, k is always explicitly provided and not inferred from the mask." << std::endl;
}

/// Return the given character in the correct case corresponding to the mask symbol.
inline char Masked(char c, bool mask) {
    int masked_difference = (c <= 'Z') - (int) mask;
    return c + (char) masked_difference * ('a' - 'A');
}

/// Reprint the sequence header as it was in the original fasta file.
void ReprintSequenceHeader(kseq_t* masked_superstring, std::ostream &of) {
    of << ">";
    if (masked_superstring->name.s) {
        of << masked_superstring->name.s;
    }
    if (masked_superstring->comment.s) {
        of << " " << masked_superstring->comment.s;
    }
    of << std::endl;
}

/// For the given masked superstring output the same superstring with mask with minimal/maximal number of ones.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void OptimizeOnes(kseq_t* masked_superstring, std::ostream &of, kh_S_t *kMers, kh_wrapper_t wrapper,
                  [[maybe_unused]] kmer_t _, int k,
                  bool complements, bool minimize) {
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t mask = ((kmer_t(1)) << (2 * k)) - 1;
    kmer_t shift = 2 * (k - 1);
    ReprintSequenceHeader(masked_superstring, of);
    uint8_t ms_validation = 0;
    for (size_t i = 0; i < masked_superstring->seq.l; ++i) {
        auto data = nucleotideToInt[(uint8_t) masked_superstring->seq.s[i]];
        ms_validation |= data;
        currentKMer = ((currentKMer << 2) | data) & mask;
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        if (i >= (size_t)k - 1) {
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            auto kmer_pointer = wrapper.kh_get_from_set(kMers, canonical);
            bool contained = kmer_pointer != kh_end(kMers);
            of << Masked(masked_superstring->seq.s[i - k + 1], contained);
            // If minimizing, erase the k-mer once set.
            if (minimize && contained) {
                wrapper.kh_del_from_set(kMers, kmer_pointer);
            }
            // Print a warning if the mask convention is violated.
            if (i == masked_superstring->seq.l - 1 && !contained) {
                PrintMaskConventionWarning();
            }
        }
    }
    // Print the remaining k-1 characters.
    for (size_t i = masked_superstring->seq.l - k + 1; i < masked_superstring->seq.l; ++i) {
        of << Masked(masked_superstring->seq.s[i], false);
    }
    of << std::endl;
    // Check that characters were only ACGTacgt.
    if (ms_validation >= 4) {
        throw std::invalid_argument("Masked superstring contains invalid characters.");
    }
}

/// Read or set the intervals.
/// If [setIntervals] is provided reprint the given files with the corresponding intervals set to 1.
/// Otherwise, read the intervals in which each k-mer occurs.
template <typename kmer_t, typename kh_P_t, typename kh_S_t, typename kh_wrapper_t>
std::pair<size_t, size_t> ReadWriteIntervals(kh_P_t *intervals, kh_S_t *kMers, kh_wrapper_t wrapper,
                             [[maybe_unused]] kmer_t _, std::vector<std::list<size_t>> &intervalsForKmer,
                             kseq_t* masked_superstring, int k, bool complements, std::ostream &of,
                             const bool* setIntervals = nullptr) {
    bool reading = setIntervals == nullptr;
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t mask = ((kmer_t(1)) << (2 * k)) - 1;
    kmer_t shift = 2 * (k - 1);
    size_t currentInterval = 0;
    size_t occurrences = 0;
    bool interval_used = false;
    uint8_t ms_validation = 0;
    for (size_t i = 0; i < masked_superstring->seq.l; ++i) {
        auto data = nucleotideToInt[(uint8_t) masked_superstring->seq.s[i]];
        ms_validation |= data;
        currentKMer = ((currentKMer << 2) | data) & mask;
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        if (i >= (size_t)k - 1) {
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            auto kmer_pointer = wrapper.kh_get_from_set(kMers, canonical);
            bool contained = kmer_pointer != kh_end(kMers);
            bool set = false;
            if (contained) {
                interval_used = true;
                if (reading) occurrences += appendInterval(intervals, wrapper, intervalsForKmer, currentKMer, currentInterval, k, complements);
                else set = setIntervals[currentInterval];
            } else {
                currentInterval += interval_used;
                interval_used = false;
            }
            if (!reading) {
                char toPrint = masked_superstring->seq.s[i - k + 1];
                of << Masked(toPrint, set);
                // Print a warning if the mask convention is violated.
                if (i == masked_superstring->seq.l - 1 && !contained) {
                    PrintMaskConventionWarning();
                }
            }
        }
    }
    // Print the remaining k-1 characters.
    if (!reading) {
        for (size_t i = masked_superstring->seq.l - k + 1; i < masked_superstring->seq.l; ++i) {
            of << Masked(masked_superstring->seq.s[i], false);
        }
    }
    if (ms_validation >= 4) {
        throw std::invalid_argument("Masked superstring contains invalid characters.");
    }
    return {occurrences, currentInterval + interval_used};
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
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void OptimizeRuns(kh_wrapper_t wrapper, kmer_t _, kseq_t* masked_superstring, kh_S_t *kMers, std::ostream &of, int k, bool complements, bool approximate) {
    auto *intervals = wrapper.kh_init_map();
    std::vector<std::list<size_t>> intervalsForKMer;
    auto [size, rows] = ReadWriteIntervals(intervals, kMers, wrapper, _, intervalsForKMer, masked_superstring, k, complements, of, nullptr);
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

    ReprintSequenceHeader(masked_superstring, of);
    ReadWriteIntervals(intervals, kMers, wrapper, _, intervalsForKMer, masked_superstring, k, complements, of, intervalsSet);
    of << std::endl;
}

template <typename kmer_t, typename kh_wrapper_t>
int Optimize(kh_wrapper_t wrapper, kmer_t _, std::string &algorithm, std::string path, std::ostream &of,  int k, bool complements) {
    kseq_t* masked_superstring = ReadMaskedSuperstring(path);
    auto *kMers = wrapper.kh_init_set();
    AddKMers(kMers, wrapper, _, masked_superstring->seq.l, masked_superstring->seq.s, k, complements, true);

    if (algorithm == "ones") {
        OptimizeOnes(masked_superstring, of, kMers, wrapper, _, k, complements, false);
    } else if (algorithm == "zeros") {
        OptimizeOnes(masked_superstring, of, kMers,  wrapper, _,k, complements, true);
    } else if (algorithm == "runs") {
        OptimizeRuns(wrapper, _, masked_superstring, kMers, of, k, complements, false);
    } else if (algorithm == "runsapprox") {
        OptimizeRuns(wrapper, _,masked_superstring, kMers, of, k, complements, true);
    } else {
        std::cerr << "Algorithm '" + algorithm + "' not recognized." << std::endl;
        kseq_destroy(masked_superstring);
        return 1;
    }
    AssertEOF(masked_superstring, "Expecting only a single FASTA record -- the masked superstring.");
    kseq_destroy(masked_superstring);
    return 0;
}
