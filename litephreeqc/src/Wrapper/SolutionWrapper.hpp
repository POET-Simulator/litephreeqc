/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
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
