// An example of the C++ API of FMSI 2 (fmsi_api.h): build an index in memory, query it, write it to a file
// and load it again. Compile with `make example` in src/ and run `../fmsi-example`.

#include <iostream>

#include "fmsi_api.h"

int main() {
    // A masked superstring of the 3-mers ACG, CGG, GGT and TAA: the characters at which the k-mers of the set
    // start are in upper case. Masked superstrings of real data are computed by KmerCamel.
    const std::string masked_superstring = "ACGgTaa";
    const int k = 3;

    fmsi::BuildOptions options;
    options.klcp = true;
    options.access = true;
    fmsi::Index index = fmsi::Index::build(masked_superstring, k, options);
    std::cout << "k = " << index.k() << ", " << index.number_of_hashes() << " k-mers, mask " << index.mask_type()
              << ", " << index.bytes() << " bytes in memory" << std::endl;

    // Membership queries: one k-mer, and all k-mers of a sequence. A k-mer and its reverse complement are the same.
    std::cout << "ACG: " << index.contains("ACG") << ", CGT (its reverse complement): " << index.contains("CGT")
              << ", ACT: " << index.contains("ACT") << std::endl;
    std::cout << "k-mers of ACGGTAAC:";
    for (uint8_t present : index.contains_all("ACGGTAAC")) std::cout << " " << (int)present;
    std::cout << std::endl;

    // Dictionary queries: with a min-one mask, lookup is a minimal perfect hash and access inverts it.
    std::cout << "minimal perfect hash: " << (index.minimal_hashes() ? "yes" : "no") << std::endl;
    for (const char *kmer : {"ACG", "CGG", "GGT", "TAA", "TTA", "ACT"}) {
        const int64_t hash = index.lookup(kmer);
        std::cout << kmer << " -> " << hash;
        if (hash >= 0) std::cout << " -> " << index.access(hash);
        std::cout << std::endl;
    }

    // Many k-mers at once are searched in lockstep, which is several times faster per k-mer on a large index.
    const std::string kmers = "ACGCGGGGTTAAACT";
    std::vector<uint8_t> present(kmers.size() / k);
    index.contains_batch(kmers.data(), present.size(), present.data());
    std::cout << "ACG CGG GGT TAA ACT:";
    for (uint8_t p : present) std::cout << " " << (int)p;
    std::cout << std::endl;

    // The index as one file (and its kLCP array next to it), loaded again in the fastest layout.
    const std::string file = index.save("example");
    fmsi::LoadOptions load;
    load.layout = fmsi::Layout::fast;
    fmsi::Index loaded = fmsi::Index::load(file, load);
    std::cout << "loaded " << file << ": GGT -> " << loaded.lookup("GGT") << std::endl;
    return 0;
}
