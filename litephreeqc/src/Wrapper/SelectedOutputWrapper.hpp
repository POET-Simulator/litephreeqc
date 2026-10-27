/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */

#pragma once

#include "CSelectedOutput.hxx"
#include "WrapperBase.hpp"
#include <cstddef>
#include <span>
#include <string>
#include <vector>

class SelectedOutputWrapper : public WrapperBase {
public:
  SelectedOutputWrapper(CSelectedOutput *selected_output,
                        const std::vector<std::string> &headings,
                        std::size_t target_row = 0);

  void get(std::span<LDBLE> &data) const override;

  void set(const std::span<LDBLE> &data) override;

  static std::vector<std::string>
  names(const CSelectedOutput *selected_output,
        std::vector<std::string> &base_names);

  static std::vector<std::string>
  names(const CSelectedOutput *selected_output);

private:
  void resolve_col_indices() const;

  CSelectedOutput *selected_output;
  std::vector<std::string> headings;
  std::size_t target_row;
  mutable std::vector<int> col_indices;
};
