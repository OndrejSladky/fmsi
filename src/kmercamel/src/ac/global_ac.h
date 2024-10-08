#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <list>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "kmers_ac.h"
#include "ac_automaton.h"


/// Represents oriented edge in the overlap graph.
struct OverlapEdge {
    // Index of the first k-mer.
    size_t firstIndex;
    // Index of the second k-mer.
    size_t secondIndex;
    // Length of the overlap of the two k-mers.
    int overlapLength;
};


/// Greedily find the approximate overlapPath path with longest overlaps using the AC automaton.
std::vector<OverlapEdge> OverlapHamiltonianPathAC (const std::vector<KMer> &kMers, bool complements) {
    size_t n = kMers.size() / (1 + complements);
    ACAutomaton automaton;
    automaton.Construct(kMers);
    std::vector<bool> forbidden(kMers.size(), false);
    std::vector<bool> prefixForbidden(kMers.size(), false);
    std::vector<std::list<size_t>> incidentKMers(automaton.states.size());
    std::vector<OverlapEdge> hamiltonianPath;
    std::vector<size_t> first(kMers.size());
    std::vector<size_t> last(kMers.size());
    for (size_t i = 0; i < kMers.size(); ++i) {
        first[i] = last[i] = i;
        incidentKMers[automaton.states[automaton.endStateIndices[i]].backwardEdge].push_back(i);
    }
    for (int s : automaton.reversedOrdering) {
        if (incidentKMers[s].empty()) continue;
        for (size_t j : automaton.states[s].supporters) {
            if (incidentKMers[s].empty()) continue;
            if (forbidden[j]) continue;
            auto i = incidentKMers[s].begin();
            // If the path forms a cycle, or is between k-mer and its reverse complement, skip this path.
            while (i != incidentKMers[s].end() && (first[*i]%n == j%n || first[*i]%n == last[j]%n|| prefixForbidden[*i])) {
                auto new_i = i;
                new_i++;
                if (prefixForbidden[*i]) incidentKMers[s].erase(i);
                i = new_i;
            }
            if (i == incidentKMers[s].end()) {
                continue;
            }
            std::vector<std::pair<size_t,size_t>> new_edges ({{*i, j}});
            if (complements) new_edges.emplace_back((j + n) % kMers.size(), (*i + n) % kMers.size());
            for (auto [x, y] : new_edges) {
                hamiltonianPath.push_back(OverlapEdge{x, y, automaton.states[s].depth});
                forbidden[y] = true;
                first[last[y]] = first[x];
                last[first[x]] = last[y];
                prefixForbidden[x] = true;
            }
            incidentKMers[s].erase(i);
        }
        incidentKMers[automaton.states[s].backwardEdge].splice(incidentKMers[automaton.states[s].backwardEdge].end(), incidentKMers[s]);
    }
    return hamiltonianPath;
}


/// Return the suffix of the given kMer without the first *overlap* chars.
std::string Suffix(const KMer &kMer, const int overlap) {
    return kMer.value.substr(overlap, kMer.length() - overlap);
}


/// Construct the superstring and the path from the given overlapPath path in the overlap graph.
void SuperstringFromPath(const std::vector<OverlapEdge> &hamiltonianPath, const std::vector<KMer> &kMers, std::ostream& of, const int k) {
    std::vector<OverlapEdge> edgeFrom (kMers.size(), OverlapEdge{size_t(-1),size_t(-1), -1});
    std::vector<bool> isStart(kMers.size(), false);
    for (auto edge : hamiltonianPath) {
        isStart[edge.firstIndex] = true;
        edgeFrom[edge.firstIndex] = edge;
    }
    for (auto edge : hamiltonianPath) {
        isStart[edge.secondIndex] = false;
    }

    // Find the vertex in the overlap graph with in-degree 0.
    size_t start = 0;
    for (; start < kMers.size() && !isStart[start]; ++start);
    // Handle the edge case of only one k-mer.
    start %= kMers.size();

    // Print the first character.
    of << kMers[start].value[0];

    // Move from the first k-mer to the last which has no successor.
    while(edgeFrom[start].secondIndex != size_t(-1)) {
        int overlapLength = edgeFrom[start].overlapLength;
        if (overlapLength != k - 1) {
            std::string unmaskedNucleotides = kMers[start].value.substr(1, k - 1 - overlapLength);
            std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
            of << unmaskedNucleotides;
        }
        of << kMers[edgeFrom[start].secondIndex].value[0];
        start = edgeFrom[start].secondIndex;
    }

    // Print the trailing k-1 characters.
    std::string unmaskedNucleotides = kMers[start].value.substr(1);
    std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
    of << unmaskedNucleotides;
}

/// Get the approximated shortest superstring of the given k-mers using the global greedy algorithm with Aho-Corasick automaton.
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, it is expected that kMers do not contain both k-mer and its reverse complement.
void GlobalAC(std::vector<KMer> kMers, std::ostream& of, bool complements) {
	if (kMers.empty()) {
		throw std::invalid_argument("input cannot be empty");
	}
    // Add complementary k-mers.
    size_t n = kMers.size();
    kMers.resize(n * (1 + complements));
    if (complements) for (size_t i = 0; i < n; ++i) {
        kMers[i + n] = ReverseComplement(kMers[i]);
    }

	const int k = (int)kMers[0].length();
    auto hamiltonianPath = OverlapHamiltonianPathAC(kMers, complements);
    SuperstringFromPath(hamiltonianPath, kMers, of, k);
}
