/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include "PhreeqcEngine.hpp"
#include "PhreeqcMatrix.hpp"
#include "PhreeqcRunner.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <cstddef>
#include <memory>
#include <set>
#include <vector>

PhreeqcRunner::PhreeqcRunner(const PhreeqcMatrix &matrix) {
  this->_ncols = matrix.get().names.size();
  // first make sure to have enough space in our buffer
  this->_buffer.reserve(this->_ncols);

  // Create a PhreeqcEngine for each id
  for (const auto &id : matrix.getIds()) {
    this->_engineStorage[id] = std::make_unique<PhreeqcEngine>(matrix, id);
    this->_cellColumns[id] = matrix.getCellColumns(id);
  }
}

// The values of the cell are taken by column index, not by skipping NaNs: a
// value of the cell may itself be NaN (e.g. a text SELECTED_OUTPUT column)
void PhreeqcRunner::run_row(std::vector<double> &row, const double time_step) {
  if (row.size() != this->_ncols) {
    throw std::invalid_argument(
        "PhreeqcRunner::run: a row has " + std::to_string(row.size()) +
        " columns, the matrix has " + std::to_string(this->_ncols));
  }
  const auto pqc_id = static_cast<int>(row[0]);
  const auto &columns = this->_cellColumns.at(pqc_id);

  this->_buffer.clear();
  for (const auto j : columns) {
    this->_buffer.push_back(row[j]);
  }

  this->_engineStorage.at(pqc_id)->runCell(this->_buffer, time_step);

  for (std::size_t k = 0; k < columns.size(); k++) {
    row[columns[k]] = this->_buffer[k];
  }
}

void PhreeqcRunner::run(std::vector<std::vector<double>> &simulationInOut,
                        const double time_step) {
  for (auto &row : simulationInOut) {
    this->run_row(row, time_step);
  }
}

void PhreeqcRunner::run(std::vector<std::vector<double>> &simulationInOut,
                        const double time_step,
                        const std::vector<std::size_t> &to_ignore) {
  const std::set<std::size_t> to_ignore_set(to_ignore.begin(), to_ignore.end());

  for (std::size_t i = 0; i < simulationInOut.size(); i++) {
    if (to_ignore_set.find(i) != to_ignore_set.end()) {
      continue;
    }
    this->run_row(simulationInOut[i], time_step);
  }
}