#include "PhreeqcEngine.hpp"
#include "PhreeqcMatrix.hpp"
#include "PhreeqcRunner.hpp"
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

PhreeqcRunner::PhreeqcRunner(const PhreeqcMatrix &matrix) {
  // first make sure to have enough space in our buffer
  this->_buffer.reserve(matrix.get().names.size());

  // Create a PhreeqcEngine for each id
  for (const auto &id : matrix.getIds()) {
    this->_engineStorage[id] = std::make_unique<PhreeqcEngine>(matrix, id);
  }
}

static void copy_to_buffer(std::vector<double> &buffer,
                           const std::vector<double> &current_inout) {
  for (std::size_t j = 1; j < current_inout.size(); j++) {
    if (std::isnan(current_inout[j])) {
      continue;
    }
    buffer.push_back(current_inout[j]);
  }
}

static void copy_from_buffer(const std::vector<double> &buffer,
                             std::vector<double> &current_inout) {
  std::size_t buffer_index = 0;
  for (std::size_t j = 1; j < current_inout.size(); j++) {
    if (std::isnan(current_inout[j])) {
      continue;
    }
    current_inout[j] = buffer[buffer_index++];
  }
}

void PhreeqcRunner::run(std::vector<std::vector<double>> &simulationInOut,
                        const double time_step) {
  for (std::size_t i = 0; i < simulationInOut.size(); i++) {
    const auto &current_inout = simulationInOut[i];
    const auto pqc_id = static_cast<int>(current_inout[0]);

    this->_buffer.clear();

    // Copy the input to the buffer while ignoring the first element and NaNs
    copy_to_buffer(this->_buffer, current_inout);

    this->_engineStorage.at(pqc_id)->runCell(this->_buffer, time_step);

    // Copy the buffer back to the output while ignoring the first element and
    // NaNs
    copy_from_buffer(this->_buffer, simulationInOut[i]);
  }
}