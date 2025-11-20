#include "fmsi_api.h"
#include <iostream>

int main() {
    std::string ms = "ACGgTaa";
    int k = 3;
    fms_index index = fmsi_construct<int64_t>(ms, k, true, true);


    std::cout << ms << std::endl;
    std::cout << "sa-t-mask: " << index.sa_transformed_mask << std::endl;
    std::cout << "index.counts: " << index.counts[0] << " ";
    std::cout << index.counts[1] << " ";
    std::cout << index.counts[2] << " ";
    std::cout << index.counts[3] << std::endl;


    std::string query0 = "ACG";
    std::string query1 = "CGG";
    std::string query2 = "GGT";
    std::string query3 = "TAA";
    std::string query4 = "ACT";
    std::string query5 = "GTA";
    
    // Single membership queries taking advantage of mask maximality (faster & recommended).
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<true>(index, (char*)query0.c_str(), k) << std::endl;
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<true>(index, (char*)query1.c_str(), k) << std::endl;
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<true>(index, (char*)query2.c_str(), k) << std::endl;
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<true>(index, (char*)query3.c_str(), k) << std::endl;
    std::cout << "expects: 0; got: " << fmsi_membership_single_query<true>(index, (char*)query4.c_str(), k) << std::endl;
    std::cout << "expects: 0; got: " << fmsi_membership_single_query<true>(index, (char*)query5.c_str(), k) << std::endl;

    // Single membership queries for general masks (only possibility if masks are not optimized).
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<false>(index, (char*)query0.c_str(), k) << std::endl;
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<false>(index, (char*)query1.c_str(), k) << std::endl;
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<false>(index, (char*)query2.c_str(), k) << std::endl;
    std::cout << "expects: 1; got: " << fmsi_membership_single_query<false>(index, (char*)query3.c_str(), k) << std::endl;
    std::cout << "expects: 0; got: " << fmsi_membership_single_query<false>(index, (char*)query4.c_str(), k) << std::endl;
    std::cout << "expects: 0; got: " << fmsi_membership_single_query<false>(index, (char*)query5.c_str(), k) << std::endl;
    
    // Minimal unique k-mer hashes (slower but guaranteedly minimal).
    int64_t hash0 = fmsi_lookup_single_query<true, false>(index, (char*)query0.c_str(), k);
    int64_t hash1 = fmsi_lookup_single_query<true, false>(index, (char*)query1.c_str(), k);
    int64_t hash2 = fmsi_lookup_single_query<true, false>(index, (char*)query2.c_str(), k);
    int64_t hash3 = fmsi_lookup_single_query<true, false>(index, (char*)query3.c_str(), k);
    int64_t hash4 = fmsi_lookup_single_query<true, false>(index, (char*)query4.c_str(), k);
    int64_t hash5 = fmsi_lookup_single_query<true, false>(index, (char*)query5.c_str(), k);
    std::cout << "expects: 0; got: " << hash0 << std::endl;
    std::cout << "expects: 1; got: " << hash1 << std::endl;
    std::cout << "expects: 2; got: " << hash2 << std::endl;
    std::cout << "expects: 3; got: " << hash3 << std::endl;
    std::cout << "expects: -1; got: " << hash4 << std::endl;
    std::cout << "expects: -1; got: " << hash5 << std::endl;

    // Corresponding access to minimal hashes.
    std::string kmer0 = fmsi_access<true>(index, hash0, k);
    std::string kmer1 = fmsi_access<true>(index, hash1, k); 
    std::string kmer2 = fmsi_access<true>(index, hash2, k); 
    std::string kmer3 = fmsi_access<true>(index, hash3, k); 
    std::cout << "expects: " << query0 << "; got: " << kmer0 << std::endl;
    std::cout << "expects: " << query1 << "; got: " << kmer1 << std::endl;
    std::cout << "expects: " << query2 << "; got: " << kmer2 << std::endl;
    std::cout << "expects: " << query3 << "; got: " << kmer3 << std::endl;

    // Non-minimal unique k-mer hashes (faster but the numbers are not minimal possible)
    hash0 = fmsi_lookup_single_query<false, true>(index, (char*)query0.c_str(), k);
    hash1 = fmsi_lookup_single_query<false, true>(index, (char*)query1.c_str(), k);
    hash2 = fmsi_lookup_single_query<false, true>(index, (char*)query2.c_str(), k);
    hash3 = fmsi_lookup_single_query<false, true>(index, (char*)query3.c_str(), k);
    hash4 = fmsi_lookup_single_query<false, true>(index, (char*)query4.c_str(), k);
    hash5 = fmsi_lookup_single_query<false, true>(index, (char*)query5.c_str(), k);
    std::cout << "expects: 3; got: " << hash0 << std::endl;
    std::cout << "expects: 4; got: " << hash1 << std::endl;
    std::cout << "expects: 5; got: " << hash2 << std::endl;
    std::cout << "expects: 7; got: " << hash3 << std::endl;
    std::cout << "expects: -1; got: " << hash4 << std::endl;
    std::cout << "expects: -1; got: " << hash5 << std::endl;

    // Corresponding access to non-minimal hashes.
    kmer0 = fmsi_access<false>(index, hash0, k);
    kmer1 = fmsi_access<false>(index, hash1, k); 
    kmer2 = fmsi_access<false>(index, hash2, k); 
    kmer3 = fmsi_access<false>(index, hash3, k); 
    std::cout << "expects: " << query0 << "; got: " << kmer0 << std::endl;
    std::cout << "expects: " << query1 << "; got: " << kmer1 << std::endl;
    std::cout << "expects: " << query2 << "; got: " << kmer2 << std::endl;
    std::cout << "expects: " << query3 << "; got: " << kmer3 << std::endl;

    std::string stream = "ACGGTACC";
    std::string rc_stream = "GGTACCGA";
    // Streamed membership queries
    auto results = fmsi_membership_streamed_query<true>(index, (char*)stream.c_str(), (char*)rc_stream.c_str(), stream.size(), k);

    std::cout << "expects: 1 1 1 0 0 1; got:";
    for (auto &&x : results) {
        std::cout << " " << x;
    }
    std::cout << std::endl;

    // Streamed membership queries
    auto hashes = fmsi_lookup_streamed_query<true, false>(index, (char*)stream.c_str(), (char*)rc_stream.c_str(), stream.size(), k);

    std::cout << "expects: 0 1 2 -1 -1 2; got:";
    for (auto &&x : hashes) {
        std::cout << " " << x;
    }
    std::cout << std::endl;

}