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

#include "PhreeqcMatrix.hpp"
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
                  bool with_redox,
                  PhreeqcMatrix::OptionalEssentials optional_essentials =
                      PhreeqcMatrix::OptionalEssentials::All);

  void get(std::span<LDBLE> &data) const override;

  void set(const std::span<LDBLE> &data) override;

  static std::vector<std::string>
  names(cxxSolution *solution, bool include_h0_o0,
        std::vector<std::string> &solution_order, bool with_redox,
        PhreeqcMatrix::OptionalEssentials optional_essentials =
            PhreeqcMatrix::OptionalEssentials::All);

  std::vector<std::string> getEssentials() const;

private:
  cxxSolution *solution;
  const std::vector<std::string> solution_order;

  // Mandatory essentials (always included): H, O, Charge, tc, patm
  static constexpr std::array<const char *, 5> MANDATORY_ESSENTIALS = {
      "H", "O", "Charge", "tc", "patm"};

  // Optional essentials (configurable): SolVol, pH, pe, MassH2O, Viscosity,
  // Density
  static constexpr std::array<const char *, 6> OPTIONAL_ESSENTIAL_NAMES = {
      "SolVol", "pH", "pe", "MassH2O", "Viscosity", "Density"};

  // Corresponding flags for each optional essential
  static constexpr std::array<PhreeqcMatrix::OptionalEssentials, 6>
      OPTIONAL_ESSENTIAL_FLAGS = {PhreeqcMatrix::OptionalEssentials::SolVol,
                                  PhreeqcMatrix::OptionalEssentials::pH,
                                  PhreeqcMatrix::OptionalEssentials::pe,
                                  PhreeqcMatrix::OptionalEssentials::MassH2O,
                                  PhreeqcMatrix::OptionalEssentials::Viscosity,
                                  PhreeqcMatrix::OptionalEssentials::Density};

  static constexpr std::size_t NUM_MANDATORY_ESSENTIALS =
      MANDATORY_ESSENTIALS.size();

  const bool _with_redox;
  const PhreeqcMatrix::OptionalEssentials _optional_essentials;

  // Compute the number of active optional essentials
  std::size_t countActiveOptionals() const;
};
