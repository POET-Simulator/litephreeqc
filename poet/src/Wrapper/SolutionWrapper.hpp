#pragma once

#include "NameDouble.h"
#include "Solution.h"
#include "WrapperBase.hpp"
#include <cstddef>
#include <vector>

class SolutionWrapper : public WrapperBase {
public:
  SolutionWrapper(cxxSolution *soln,
                  const std::vector<std::string> &solution_order);

  void get(std::span<LDBLE> &data) const;

  void set(const std::span<LDBLE> &data);

  static std::vector<std::string>
  names(cxxSolution *solution, std::vector<std::string> &solution_order);

private:
  cxxSolution *solution;
  const std::vector<std::string> solution_order;

  static constexpr std::size_t NUM_ESSENTIALS = 3;
};