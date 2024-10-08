#pragma once
#include "../src/ac/local_ac.h"

#include "gtest/gtest.h"

namespace {
    TEST(LocalAC, ExtensionAC) {
        struct TestCase {
            std::vector<bool> forbidden;
            std::list<size_t> incidentKMers;
            bool complements;
            int wantResult;
            std::vector<bool> wantForbidden;
            std::list<size_t> wantIncidentKMers;
        };
        std::vector<TestCase> tests = {
                {{0,0,0}, {1, 0}, false, 1, {0, 1, 0}, {0}},
                {{0,1,1}, {1, 0}, false, 0, {1, 1, 1}, {}},
                {{0,0,0, 0}, {}, false, -1, {0,0,0,0}, {}},
                {{1,1,0, 1}, {3, 0, 1}, false, -1, {1, 1, 0, 1}, {}},
                {{0,0,0,0}, {3, 0, 1}, true, 3, {0, 1, 0, 1}, {0,1}},
        };

        for (auto t: tests) {
            size_t gotResult = ExtensionAC(t.forbidden, t.incidentKMers, t.complements);

            EXPECT_EQ(t.wantResult, gotResult);
            EXPECT_EQ(t.wantForbidden, t.forbidden);
            EXPECT_EQ(t.wantIncidentKMers, t.incidentKMers);
        }
    }

    TEST(LocalAC, LocalAC) {
        struct TestCase {
            std::vector<KMer> kMers;
            int k;
            int d_max;
            bool complements;
            std::string wantSuperstring;
        };
        std::vector<TestCase> tests = {
                // As behavior of the unordered_set.begin() is not specified, some tests are commented, as they could fail otherwise.
                // Uncommenting them may add additional check but could also add false positives.
                //{ {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, 4, 3,
                //  "AACAATTTCCCC", {1,1,0,0,1, 0,0,0, 1, 0,0,0}},
                //{ {KMer{"G"}, KMer{"T"}, KMer{"A"}}, 1, 0, false,
                //        "GTA"},
                { {KMer{"GCT"}, KMer{"TAA"}, KMer{"AAA"}}, 3, 2, false,
                        "GcTAaa"},
                { {KMer{"TAA"}, KMer{"AAA"}, KMer{"GCT"}}, 3, 2, false,
                        "GcTAaa"},
                {{KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}, KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"}}, 31, 5, false,
                        "TTtcttttttttttttttttttttttttttga"},
                { {KMer{"TAA"}, KMer{"TTT"}}, 3, 2, true,
                        "TAaa"},
        };

        for (auto t: tests) {
            std::stringstream of;

            LocalAC(t.kMers, of, t.k, t.d_max, t.complements);

            EXPECT_EQ(t.wantSuperstring, of.str());
        }
    }


}
