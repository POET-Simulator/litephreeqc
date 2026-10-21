/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */

#pragma once

#include "IPhreeqc.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <string>

class IPhreeqcReader {
public:
  IPhreeqcReader(const std::string &database, const std::string &script);

  void setOutputID(std::uint32_t cell_id);
  double operator[](const std::string &name) const;

  void run(double dt, std::vector<std::uint32_t> &&cell_ids);

private:
  std::unique_ptr<IPhreeqc> _m_pqc;

  std::string _m_raw_output = "";

  const std::map<std::string, std::string> _m_rename_map = {
      {"H", "-total_h"}, {"O", "-total_o"},     {"Charge", "-cb"},
      {"tc", "-temp"},   {"patm", "-pressure"}, {"SolVol", "-soln_vol"},
      {"pH", "-pH"},     {"pe", "-pe"}};
};
