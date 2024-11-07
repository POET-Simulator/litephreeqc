#pragma once

#include "Solution.h"
#include "WrapperBase.hpp"
#include <array>
#include <cstddef>
#include <string>
#include <vector>

class SolutionWrapper : public WrapperBase {
public:
  SolutionWrapper(cxxSolution *soln,
                  const std::vector<std::string> &solution_order);

  void get(std::span<LDBLE> &data) const;

  void set(const std::span<LDBLE> &data);

  static std::vector<std::string>
  names(cxxSolution *solution, std::vector<std::string> &solution_order);

  std::vector<std::string> getEssentials() const;

private:
  cxxSolution *solution;
  const std::vector<std::string> solution_order;

  static constexpr std::array<const char *, 5> ESSENTIALS = {"H", "O", "Charge",
                                                             "H(0)", "O(0)"};

  static constexpr std::size_t NUM_ESSENTIALS = ESSENTIALS.size();
};