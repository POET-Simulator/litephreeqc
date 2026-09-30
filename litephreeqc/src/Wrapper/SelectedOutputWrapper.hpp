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
