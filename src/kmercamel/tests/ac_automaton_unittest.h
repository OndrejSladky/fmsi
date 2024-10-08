#pragma once
#include "../src/ac/ac_automaton.h"

#include "gtest/gtest.h"


TEST(ACAutomaton, ConstructTrie) {
    ACAutomaton a;
    std::vector<ACState> wantStates = {
            ACState{std::list<size_t>{0, 1}, 1, 0, 0, 0, 0, 0, 0},
            ACState{std::list<size_t>{0, 1}, -1, 2, -1, -1, 0, 1, 1},
            ACState{std::list<size_t>{0, 1}, -1, -1, 4, 3, 0, 2, 2},
            ACState{std::list<size_t>{0}, -1, -1, -1, -1, 0, 3, 3},
            ACState{std::list<size_t>{1}, -1, -1, -1, -1, 0, 3, 4},
    };
    std::vector<int> wantEndStateIndices = {3, 4};

    a.ConstructTrie(std::vector<KMer>{KMer{"ACT"}, KMer{"ACG"}});

    ASSERT_EQ(wantStates.size(), a.states.size());
    for (size_t i = 0; i < wantStates.size(); ++i) {
        EXPECT_EQ(wantStates[i].backwardEdge, a.states[i].backwardEdge);
    for (int j = 0; j < 4; ++j) {
        EXPECT_EQ(wantStates[i].forwardEdges[j], a.states[i].forwardEdges[j]);
    }
    EXPECT_EQ(wantStates[i].depth, a.states[i].depth);
    EXPECT_EQ(wantStates[i].id, a.states[i].id);
    EXPECT_EQ(wantStates[i].supporters, a.states[i].supporters);
}

ASSERT_EQ(wantEndStateIndices, a.endStateIndices);
}

TEST(ACAutomaton, ConstructBackwardEdges) {
    ACAutomaton a;
    // Trie representing ["ACT", "ACA"].
    a.states = {
            ACState{std::list<size_t>{0, 1}, 1, 0, 0, 0, 0, 0, 0},
            ACState{std::list<size_t>{0, 1}, -1, 2, -1, -1, 0, 1, 1},
            ACState{std::list<size_t>{0, 1}, 4, -1, -1, 3, 0, 2, 2},
            ACState{std::list<size_t>{0}, -1, -1, -1, -1, 0, 3, 3},
            ACState{std::list<size_t>{1}, -1, -1, -1, -1, 0, 3, 4},
    };
    std::vector<ACState> wantStates = {
            ACState{std::list<size_t>{0, 1}, 1, 0, 0, 0, 0, 0, 0},
            ACState{std::list<size_t>{0, 1}, -1, 2, -1, -1, 0, 1, 1},
            ACState{std::list<size_t>{0, 1}, 4, -1, -1, 3, 0, 2, 2},
            ACState{std::list<size_t>{0}, -1, -1, -1, -1, 0, 3, 3},
            ACState{std::list<size_t>{1}, -1, -1, -1, -1, 1, 3, 4},
    };
    std::vector<int> wantReversedOrdering = {3,4,2,1,0};

    a.ConstructBackwardEdges();

    for (size_t i = 0; i < wantStates.size(); ++i) {
        EXPECT_EQ(wantStates[i].backwardEdge, a.states[i].backwardEdge);
    }
    EXPECT_EQ(wantReversedOrdering, a.reversedOrdering);
}
