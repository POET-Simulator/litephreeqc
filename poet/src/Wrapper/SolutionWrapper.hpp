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
                  const std::vector<std::string> &solution_order,
                  bool with_redox);

  void get(std::span<LDBLE> &data) const;

  void set(const std::span<LDBLE> &data);

  static std::vector<std::string>
  names(cxxSolution *solution, bool include_h0_o0,
        std::vector<std::string> &solution_order, bool with_redox);

  std::vector<std::string> getEssentials() const;

private:
  cxxSolution *solution;
  const std::vector<std::string> solution_order;

  static constexpr std::array<const char *, 6> ESSENTIALS = {"H", "O",
                                                             "Charge",
                                                             "SolVol", "pH","pe"}; // MDL

  static constexpr std::size_t NUM_ESSENTIALS = ESSENTIALS.size();

  const bool _with_redox;
};
