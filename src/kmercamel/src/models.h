#pragma once
#include <string>
#include <vector>

struct KMer {
    std::string value;
    size_t length() const{ return value.size(); }
};
