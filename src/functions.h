#pragma once

#include <functional>
#include <string>
#include <stdexcept>

bool f_or(size_t ones, [[maybe_unused]] size_t total) { return ones; }

bool f_and(size_t ones, size_t total) { return ones == total; }

bool f_xor(size_t ones, [[maybe_unused]] size_t total) { return ones % 2; }

bool f_r_to_s(size_t ones, [[maybe_unused]] size_t total, size_t r, size_t s) {
  return ones <= s && ones >= r;
}

std::function<bool(int, int)> mask_function(std::string name) {
  if (name == "or")
    return &f_or;
  if (name == "and")
    return &f_and;
  if (name == "xor")
    return &f_xor;
  if (name == "default")
    return nullptr;

  for (size_t i = 1; i + 1 < name.size(); ++i) {
    if (name[i] == '-') {
      std::string r_string = name.substr(0, i),
                  s_string = name.substr(i + 1, name.size() - i - 1);
      bool valid = true;
      for (size_t j = 0; j < r_string.size(); ++j)
        if (r_string[j] < '0' || r_string[j] > '9')
          valid = false;
      for (size_t j = 0; j < s_string.size(); ++j)
        if (s_string[j] < '0' || s_string[j] > '9')
          valid = false;
      if (!valid)
        break;
      size_t r = std::stoi(r_string), s = std::stoi(s_string);
      return [r, s](size_t ones, size_t total) {
        return f_r_to_s(ones, total, r, s);
      };
    }
  }
  throw std::invalid_argument("unknown function name");
}