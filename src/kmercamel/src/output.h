#pragma once

#include <iostream>
#include <string>
#include <vector>

/// Print the fasta file header.
void WriteName(const int k, std::ostream &of) {
    of << ">superstring ";
    of << "k=" << k << std::endl;
}
