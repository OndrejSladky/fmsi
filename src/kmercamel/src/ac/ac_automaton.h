#pragma once

#include <list>
#include <vector>

#include "kmers_ac.h"

constexpr int INVALID_STATE = -1;
struct ACState {
    // List of k-mers whose prefix this state is.
    std::list<size_t> supporters;
    // Indexes of the states longer by the corresponding nucleotide.
    int forwardEdges[4] =  {INVALID_STATE, INVALID_STATE, INVALID_STATE, INVALID_STATE};
    // Where to go if searching failed.
    int backwardEdge = 0;
    // Length of the corresponding string.
    int depth = 0;
    // Index of the state in the automaton's *states* list.
    int id = 0;
};

struct ACAutomaton {
    std::vector<ACState> states;
    // Indices of the states where i-th k-mer ends.
    std::vector<int> endStateIndices;
    // Reversed BFS order.
    std::vector<int> reversedOrdering;

    /// Append a new state and return its ID.
    int AddState(const int depth) {
        ACState newState;
        newState.depth = depth;
        newState.id = (int)states.size();
        states.push_back(newState);
        return (int)states.size() - 1;
    }

    /// Generate the trie from the given k-mers and set *endStateIndices*.
    void ConstructTrie(const std::vector<KMer> &kMers) {
        endStateIndices = std::vector<int>(kMers.size());
        // Initialize the root.
        AddState(0);
        for (size_t i = 0; i < kMers.size(); ++i) states[0].supporters.push_back(i);
        for (size_t i = 0; i < kMers.size(); ++i) {
            int state = 0;
            for (char c : kMers[i].value) {
                int index = NucleotideToInt(c);
                // If the next state does not yet exist, create it.
                if (states[state].forwardEdges[index] == INVALID_STATE) {
                    int newState =  AddState(states[state].depth + 1);
                    states[state].forwardEdges[index] = newState;
                }
                state = states[state].forwardEdges[index];
                states[state].supporters.push_back(i);
            }
            endStateIndices[i] = state;
        }

        // Create a forward edge from the root to itself so that the AC Step always finds a valid forward edge.
        for (int & forwardEdge : states[0].forwardEdges) {
            if (forwardEdge == INVALID_STATE) {
                forwardEdge = 0;
            }
        }
    }

    /// Do one step of the AC algorithm from the given state with a nucleotide with a given index.
    int Step (int state, const int index) {
        while (states[state].forwardEdges[index] == INVALID_STATE) {
            state = states[state].backwardEdge;
        }
        return states[state].forwardEdges[index];
    }

    /// Construct the fail edges for the trie that already had been created.
    /// Also construct the reversed BFS ordering.
    void ConstructBackwardEdges () {
        reversedOrdering = std::vector<int> (states.size(), 0);
        std::queue<int> q;
        for (int & forwardEdge : states[0].forwardEdges) {
            if (forwardEdge != 0) {
                q.push(forwardEdge);
            }
        }
        size_t orderingIndex = reversedOrdering.size() - size_t(2);
        while (!q.empty()) {
            int state = q.front();
            reversedOrdering[orderingIndex--] = state;
            q.pop();
            for (int i = 0; i < 4; ++i) {
                int nextState = states[state].forwardEdges[i];
                if (nextState != INVALID_STATE && nextState != 0) {
                    states[nextState].backwardEdge = Step(states[state].backwardEdge, i);
                    q.push(nextState);
                }
            }
        }
    }

    /// Construct the Aho-Corasick automaton.
    void Construct (const std::vector<KMer> &kMers) {
        ConstructTrie(kMers);
        ConstructBackwardEdges();
    }
};
