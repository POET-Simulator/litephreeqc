/*
 * This project is subject to the original PHREEQC license. `litephreeqc` is a
 * version of the PHREEQC code that has been modified to be used as a library.
 *
 * It adds a C++ interface on top of the original PHREEQC code, with small
 * changes to the original code base.
 *
 * Authors of Modifications:
 * - Max Luebke (mluebke@uni-potsdam.de) - University of Potsdam
 * - Marco De Lucia (delucia@gfz.de) - GFZ Helmholz Centre for Geosciences
 *
 */

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

  static constexpr std::array<const char *, 8> ESSENTIALS = {
      "H",      "O",  "Charge", "tc", "patm",

      "SolVol", "pH", "pe"}; // MDL; ML: only output

  static constexpr std::size_t NUM_ESSENTIALS = ESSENTIALS.size();

  const bool _with_redox;
};
