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

#include "PhreeqcMatrix.hpp"

#include <IPhreeqc.hpp>

std::vector<int> PhreeqcMatrix::getIds() const {
  std::vector<int> ids;
  for (const auto &[id, _] : _m_map) {
    ids.push_back(id);
  }

  return ids;
}

// std::array<std::size_t, 5> PhreeqcMatrix::getComponentCount(int cell_id)
// const {
//   std::array<std::size_t, 5> counts = {0, 0, 0, 0, 0};

//   if (!this->_m_internal_names.contains(cell_id)) {
//     return counts;
//   }

//   for (const auto &element : this->_m_internal_names.at(cell_id)) {
//     counts[static_cast<std::size_t>(element.type)]++;
//   }

//   return counts;
// }

std::map<int, std::string> PhreeqcMatrix::getDumpStringsPQI() const {
  std::map<int, std::string> dumps;

  for (const auto &cell_id : this->getIds()) {
    dumps[cell_id] = this->getDumpStringsPQI(cell_id);
  }

  return dumps;
}

std::string PhreeqcMatrix::getDumpStringsPQI(int cell_id) const {
  this->_m_pqc->SetDumpStringOn(true);

  const std::string call_string =
      "DUMP\n -cells " + std::to_string(cell_id) + "\nEND\n";

  this->_m_pqc->RunString(call_string.c_str());

  const std::string dump_string = this->_m_pqc->GetDumpString();

  this->_m_pqc->SetDumpStringOn(false);

  return dump_string;
}

std::string PhreeqcMatrix::getDatabase() const { return _m_database; }

bool PhreeqcMatrix::checkIfExists(int cell_id) const {
  return this->_m_map.find(cell_id) != this->_m_map.end();
}