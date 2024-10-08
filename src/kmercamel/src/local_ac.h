#pragma once

#include <string>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <fstream>

#include "models.h"
#include "ac_automaton.h"

/// Find the index of the first extending k-mer from the incidentKMers which is not forbidden.
/// Mark this k-mer forbidden and remove all the k-mers in incidentKMers which are forbidden.
/// The latter is done in order not to increase time complexity by repeatedly iterating over forbidden k-mers.
/// If complements is true, it is expected that index i and i + size/2 in forbidden represent complementary k-mers.
/// Return -1 if only forbidden k-mers were found.
size_t ExtensionAC(std::vector<bool> &forbidden, std::list<size_t> &incidentKMers, bool complements) {
    size_t n = forbidden.size() / (1 + complements);
    while(!incidentKMers.empty()) {
        if(forbidden[incidentKMers.front()]) {
            incidentKMers.pop_front();
        } else {
            forbidden[incidentKMers.front()] = true;
            // Mark complementary k-mer as used (if complements is set to false, this is k-mer itself).
            forbidden[(incidentKMers.front() + n) % forbidden.size()] = true;
            size_t ret = incidentKMers.front();
            incidentKMers.pop_front();
            return ret;
        }
    }
    return -1;
}

/// Get the approximated shortest superstring of the given k-mers using the local greedy algorithm with Aho-Corasick automaton.
///
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, it is expected that kMers do not contain both k-mer and its reverse complement.
void LocalAC(std::vector<KMer> kMers, std::ostream& of, int k, int d_max, bool complements) {
    // Add complementary k-mers.
    size_t n = kMers.size();
    kMers.resize(n * (1 + complements));
    if (complements) for (size_t i = 0; i < n; ++i) {
        kMers[i + n] = ReverseComplement(kMers[i]);
    }

    ACAutomaton a;
    a.Construct(kMers);

    // suffixes[i][j] is the state in the AC automaton given by the suffix of kMers[i] of size j (or -1 is none).
    std::vector<std::vector<int>> suffixes(kMers.size(), std::vector<int> (k + 1, -1));
    // suffixes[i][j] is the state in the AC automaton given by the prefix of kMers[i] of size j.
    std::vector<std::vector<int>> prefixes(kMers.size(), std::vector<int>(k + 1, 0));
    // For each state the list of k-mers which have the given state as a suffix.
    std::vector<std::list<size_t>> incidentKMers(a.states.size());
    // true if the given k-mer has already been used.
    std::vector<bool> forbidden(kMers.size(), false);

    for (size_t i = 0; i < kMers.size(); ++i) {
        for(int j = 0; j < k; ++j) {
            prefixes[i][j + 1] = a.states[prefixes[i][j]].forwardEdges[NucleotideToInt(kMers[i].value[j])];
        }
        for (int s = a.endStateIndices[i]; ; s = a.states[s].backwardEdge) {
            suffixes[i][a.states[s].depth] = s;
            incidentKMers[s].push_back(i);
            if (s == 0) break;
        }
    }

    size_t firstUnused = 0;
    for(;;) {
        // Find the first unused k-mer.
        while(forbidden[firstUnused]) {
            ++firstUnused;
            if (firstUnused == kMers.size()) {
                firstUnused = size_t(-1);
                break;
            }
        }
        if (firstUnused == size_t(-1)) break;
        std::list<char> simplitig = {kMers[firstUnused].value[0]};
        // Maintain the left and right most k-mer of the generalized simplitig.
        size_t firstKMer = firstUnused;
        size_t lastKMer = firstUnused;
        forbidden[firstUnused] = true;
        // Forbid the complementary k-mer.
        forbidden[(firstUnused + n) % forbidden.size()] = true;
        int d_l = 1, d_r = 1;
        while (d_l <= d_max || d_r <= d_max) {
            if (d_r <= d_l) {
                int state = suffixes[lastKMer][k - d_r];
                size_t ext = -1;
                if (state != -1) ext = ExtensionAC(forbidden, a.states[state].supporters, complements);
                if (ext == size_t(-1)) {
                    // No right extension found.
                    ++d_r;
                } else {
                    // Extend the generalized simplitig to the right.
                    for (int i = 1; i < d_r; ++i) simplitig.emplace_back((char)std::tolower(kMers[lastKMer].value[i]));
                    simplitig.emplace_back(kMers[lastKMer].value[d_r]);
                    lastKMer = ext;
                    d_r = 1;
                }
            } else {
                int state = prefixes[firstKMer][k - d_l];
                size_t ext = ExtensionAC(forbidden, incidentKMers[state], complements);
                if (ext == size_t(-1)) {
                    // No left extension found.
                    ++d_l;
                } else {
                    // Extend the simplitig to the left.
                    for (int i = d_l - 1; i > 0; --i) simplitig.emplace_front((char)std::tolower(kMers[ext].value[i]));
                    simplitig.emplace_front(kMers[ext].value[0]);
                    firstKMer = ext;
                    d_l = 1;
                }
            }
        }
        for (int i = 1; i < k; ++i) simplitig.emplace_back((char)std::tolower(kMers[lastKMer].value[i]));
        of << std::string(simplitig.begin(), simplitig.end());
    }
}

