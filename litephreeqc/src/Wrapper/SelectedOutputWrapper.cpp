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

#include "SelectedOutputWrapper.hpp"
#include <limits>
#include <map>

static inline std::string trim(const std::string &str) {
  const std::size_t first = str.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const std::size_t last = str.find_last_not_of(" \t\r\n");
  return str.substr(first, last - first + 1);
}

SelectedOutputWrapper::SelectedOutputWrapper(
    CSelectedOutput *selected_output,
    const std::vector<std::string> &headings,
    std::size_t target_row)
    : selected_output(selected_output), headings(headings),
      target_row(target_row) {
  this->num_elements = headings.size();
  this->resolve_col_indices();
}

void SelectedOutputWrapper::resolve_col_indices() const {
  if (this->selected_output == nullptr ||
      this->selected_output->GetColCount() == 0) {
    return;
  }

  std::map<std::string, int> heading_to_col;
  for (size_t col = 0; col < this->selected_output->GetColCount(); ++col) {
    CVar v = this->selected_output->Get(0, static_cast<int>(col));
    std::string h = (v.type == TT_STRING && v.sVal) ? v.sVal : "";
    heading_to_col[trim(h)] = static_cast<int>(col);
  }

  this->col_indices.clear();
  this->col_indices.reserve(this->headings.size());
  for (const auto &name : this->headings) {
    const std::string trimmed_name = trim(name);
    auto it = heading_to_col.find(trimmed_name);
    if (it != heading_to_col.end()) {
      this->col_indices.push_back(it->second);
    } else {
      // Try stripping "_SO" suffix if present
      std::string stripped = trimmed_name;
      if (stripped.ends_with("_SO")) {
        stripped = stripped.substr(0, stripped.size() - 3);
      }
      auto it2 = heading_to_col.find(stripped);
      if (it2 != heading_to_col.end()) {
        this->col_indices.push_back(it2->second);
      } else {
        this->col_indices.push_back(-1);
      }
    }
  }
}

void SelectedOutputWrapper::get(std::span<LDBLE> &data) const {
  if (this->col_indices.empty() ||
      this->col_indices.size() != this->num_elements) {
    this->resolve_col_indices();
  }

  if (this->selected_output == nullptr ||
      this->selected_output->GetColCount() == 0) {
    for (size_t i = 0; i < data.size(); ++i) {
      data[i] = std::numeric_limits<LDBLE>::quiet_NaN();
    }
    return;
  }

  const size_t row_count = this->selected_output->GetRowCount();
  // row_count includes row 0 (headings). If row_count <= 1, no data rows have been punched.
  if (row_count <= 1) {
    for (size_t i = 0; i < data.size(); ++i) {
      data[i] = std::numeric_limits<LDBLE>::quiet_NaN();
    }
    return;
  }

  // target_row == 0 means latest available data row (row_count - 1)
  size_t row_to_read = this->target_row;
  if (row_to_read == 0) {
    row_to_read = row_count - 1;
  } else if (row_to_read >= row_count) {
    for (size_t i = 0; i < data.size(); ++i) {
      data[i] = std::numeric_limits<LDBLE>::quiet_NaN();
    }
    return;
  }

  for (size_t i = 0; i < data.size() && i < this->col_indices.size(); ++i) {
    const int col = this->col_indices[i];
    if (col < 0 ||
        static_cast<size_t>(col) >= this->selected_output->GetColCount()) {
      data[i] = std::numeric_limits<LDBLE>::quiet_NaN();
      continue;
    }

    const CVar var =
        this->selected_output->Get(static_cast<int>(row_to_read), col);
    switch (var.type) {
    case TT_DOUBLE:
      data[i] = var.dVal;
      break;
    case TT_LONG:
      data[i] = static_cast<LDBLE>(var.lVal);
      break;
    default:
      data[i] = std::numeric_limits<LDBLE>::quiet_NaN();
      break;
    }
  }
}

void SelectedOutputWrapper::set(const std::span<LDBLE> &data) {
  // Selected output is output-only and cannot be written back to Phreeqc state.
  (void)data;
}

std::vector<std::string>
SelectedOutputWrapper::names(const CSelectedOutput *selected_output,
                             std::vector<std::string> &base_names) {
  if (selected_output == nullptr || selected_output->GetColCount() == 0) {
    return {};
  }

  const size_t col_count = selected_output->GetColCount();
  base_names.reserve(col_count);
  std::vector<std::string> names;
  names.reserve(col_count);

  for (size_t col = 0; col < col_count; ++col) {
    const CVar heading_var =
        selected_output->Get(0, static_cast<int>(col));
    const std::string raw =
        (heading_var.type == TT_STRING && heading_var.sVal)
            ? heading_var.sVal
            : "";
    const std::string h = trim(raw);
    base_names.push_back(h);
    names.push_back(h + "_SO");
  }

  return names;
}

std::vector<std::string>
SelectedOutputWrapper::names(const CSelectedOutput *selected_output) {
  std::vector<std::string> placeholder;
  return names(selected_output, placeholder);
}
