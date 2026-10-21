/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include <PhreeqcMatrix.hpp>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <vector>

PhreeqcMatrix PhreeqcMatrix::subset(const std::vector<int> &indices) const {
  PhreeqcMatrix result(*this);

  result._m_map.clear();
  result._m_internal_names.clear();

  for (const auto &index : indices) {
    result._m_map[index] = _m_map.at(index);
    result._m_internal_names[index] = _m_internal_names.at(index);
  }

  result.remove_NaNs();

  return result;
}

PhreeqcMatrix PhreeqcMatrix::erase(const std::vector<int> &indices) const {
  PhreeqcMatrix result(*this);

  for (const auto &index : indices) {
    result._m_map.erase(index);
    result._m_internal_names.erase(index);
  }

  result.remove_NaNs();

  return result;
}

// PhreeqcMatrix &PhreeqcMatrix::subset(const std::vector<int> &indices) {
//   std::map<int, std::vector<element>> new_map;
//   std::map<int, std::vector<base_names>> new_internal_names;

//   for (const auto &index : indices) {
//     new_map[index] = _m_map.at(index);
//     new_internal_names[index] = _m_internal_names.at(index);
//   }

//   this->_m_map = new_map;
//   this->_m_internal_names = new_internal_names;

//   this->remove_NaNs();

//   return *this;
// }

// PhreeqcMatrix &PhreeqcMatrix::erase(const std::vector<int> &indices) {
//   for (const auto &index : indices) {
//     _m_map.erase(index);
//     _m_internal_names.erase(index);
//   }

//   this->remove_NaNs();

//   return *this;
// }

PhreeqcMatrix::STLExport PhreeqcMatrix::get(VectorExportType type,
                                            bool include_id) const {
  const auto &first_element = _m_map.begin()->second;

  const std::size_t cols = first_element.size();
  const std::size_t rows = _m_map.size();
  const std::size_t total_elements = cols * rows;

  STLExport result;

  if (include_id) {
    result.names.push_back("ID");
  }

  for (const auto &element : first_element) {
    if (element.type != PhreeqcMatrix::PhreeqcComponent::SOLUTION) {
      break;
    }
    result.names.push_back(element.name);
  }

  for (std::size_t component = 1;
       component <= static_cast<std::size_t>(
                        PhreeqcMatrix::PhreeqcComponent::SELECTED_OUTPUT);
       component++) {
    std::vector<std::string> values;
    for (const auto &[id, elements] : _m_map) {
      std::vector<std::string> names;
      for (const auto &element : elements) {
        if (static_cast<std::size_t>(element.type) == component) {
          names.push_back(element.name);
        }
      }
      std::vector<std::string> union_names;
      std::set_union(values.begin(), values.end(), names.begin(), names.end(),
                     std::back_inserter(union_names));
      values = union_names;
    }
    result.names.insert(result.names.end(), values.begin(), values.end());
  }

  result.roles = this->roles_of(result.names, include_id);

  result.values.reserve(total_elements + (include_id ? rows : 0));

  if (type == VectorExportType::COLUMN_MAJOR) {
    std::size_t column_index = 0;

    if (include_id) {
      for (const auto &[id, _] : _m_map) {
        result.values.push_back(id);
      }
      column_index++;
    }

    for (; column_index < result.names.size(); column_index++) {
      for (const auto &[_, elements] : _m_map) {
        double value_to_add = std::numeric_limits<double>::quiet_NaN();
        // double value_to_add;
        for (const auto &curr_element : elements) {
          const std::string &curr_element_name = curr_element.name;

          if (curr_element_name == result.names[column_index]) {
            value_to_add = curr_element.value;
            break;
          }
        }
        result.values.push_back(value_to_add);
      }
    }

    return result;
  }

  for (const auto &[id, element_vec] : _m_map) {
    std::size_t column_index = 0;
    if (include_id) {
      result.values.push_back(id);
      column_index++;
    }

    for (; column_index < result.names.size(); column_index++) {
      double value_to_add = std::numeric_limits<double>::quiet_NaN();
      for (const auto &curr_element : element_vec) {
        const std::string &curr_element_name = curr_element.name;

        if (curr_element_name == result.names[column_index]) {
          value_to_add = curr_element.value;
          break;
        }
      }
      result.values.push_back(value_to_add);
    }
  }

  return result;
}

std::vector<PhreeqcMatrix::ColumnRole>
PhreeqcMatrix::roles_of(const std::vector<std::string> &names,
                        bool include_id) const {
  std::vector<ColumnRole> roles;
  roles.reserve(names.size());

  for (std::size_t i = 0; i < names.size(); i++) {
    if (include_id && i == 0) {
      roles.push_back(ColumnRole::ID);
      continue;
    }
    // a column has the same role in all cells that have it
    ColumnRole role = ColumnRole::IMMOBILE;
    for (const auto &[_, elements] : _m_map) {
      const auto it = std::find_if(
          elements.begin(), elements.end(),
          [&](const element &el) { return el.name == names[i]; });
      if (it != elements.end()) {
        role = it->role;
        break;
      }
    }
    roles.push_back(role);
  }

  return roles;
}

std::vector<PhreeqcMatrix::ColumnRole>
PhreeqcMatrix::getColumnRoles(bool include_id) const {
  return this->get(VectorExportType::ROW_MAJOR, include_id).roles;
}

std::vector<std::size_t> PhreeqcMatrix::getCellColumns(int id,
                                                       bool include_id) const {
  const auto names = this->get(VectorExportType::ROW_MAJOR, include_id).names;
  const auto &elements = _m_map.at(id);

  std::vector<std::size_t> columns;
  for (std::size_t j = include_id ? 1 : 0; j < names.size(); j++) {
    if (std::any_of(elements.begin(), elements.end(),
                    [&](const element &el) { return el.name == names[j]; })) {
      columns.push_back(j);
    }
  }
  return columns;
}

std::string PhreeqcMatrix::roleName(ColumnRole role) {
  switch (role) {
  case ColumnRole::ID:
    return "id";
  case ColumnRole::TRANSPORTED:
    return "transported";
  case ColumnRole::IMMOBILE:
    return "immobile";
  case ColumnRole::PARAMETER:
    return "parameter";
  case ColumnRole::HINT:
    return "hint";
  case ColumnRole::DIAGNOSTIC:
    return "diagnostic";
  case ColumnRole::SELECTED_OUTPUT:
    return "selected_output";
  }
  return "unknown";
}

std::vector<std::string> PhreeqcMatrix::getSolutionNames() const {
  std::vector<std::string> names;

  const auto &first_element = _m_map.begin()->second;

  for (const auto &element : _m_map.begin()->second) {
    // assuming the element vector always starts with the solution components
    if (element.type != PhreeqcComponent::SOLUTION) {
      break;
    }

    names.push_back(element.name);
  }

  return names;
}

template <PhreeqcMatrix::base_names::Components comp>
static inline std::vector<std::string> get_names_from_internal_vector(
    const std::map<int, std::vector<PhreeqcMatrix::base_names>> &map, int id) {

  const auto &it = map.find(id);

  if (it == map.end()) {
    return {};
  }

  std::vector<std::string> names;

  for (const auto &element : it->second) {
    if (element.type == comp) {
      names.push_back(element.name);
    }
  }

  return names;
}

std::vector<std::string> PhreeqcMatrix::getExchanger(int id) const {
  return get_names_from_internal_vector<base_names::Components::EXCHANGER>(
      this->_m_internal_names, id);
}

std::vector<std::string> PhreeqcMatrix::getKineticsNames(int id) const {
  return get_names_from_internal_vector<base_names::Components::KINETICS>(
      this->_m_internal_names, id);
}

std::vector<std::string> PhreeqcMatrix::getEquilibriumNames(int id) const {
  return get_names_from_internal_vector<base_names::Components::EQUILIBRIUM>(
      this->_m_internal_names, id);
}

std::vector<std::string> PhreeqcMatrix::getSurfaceCompNames(int id) const {
  return get_names_from_internal_vector<base_names::Components::SURACE_COMP>(
      this->_m_internal_names, id);
}

std::vector<std::string> PhreeqcMatrix::getSurfaceChargeNames(int id) const {
  return get_names_from_internal_vector<base_names::Components::SURFACE_CHARGE>(
      this->_m_internal_names, id);
}

std::vector<std::string> PhreeqcMatrix::getSolutionPrimaries() const {
  return std::vector<std::string>(_m_surface_primaries.begin(),
                                  _m_surface_primaries.end());
}

double PhreeqcMatrix::operator()(int cell_id, const std::string &name) const {
  const auto &elements = _m_map.at(cell_id);

  const auto it = std::find_if(
      elements.begin(), elements.end(),
      [&name](const element &element) { return element.name == name; });

  if (it == elements.end()) {
    throw std::runtime_error("Element not found");
  }

  return it->value;
}

// MDL
// Column names (in matrix order) whose component type is `type` in any cell
std::vector<std::string>
PhreeqcMatrix::getMatrixNamesOfType(PhreeqcComponent type) const {
  std::set<std::string> of_type;
  for (const auto &[id, elements] : _m_map) {
    for (const auto &element : elements) {
      if (element.type == type) {
        of_type.insert(element.name);
      }
    }
  }

  std::vector<std::string> names;
  for (const auto &name : this->get().names) {
    if (of_type.contains(name)) {
      names.push_back(name);
    }
  }
  return names;
}

std::vector<std::string> PhreeqcMatrix::getMatrixKinetics() const {
  return getMatrixNamesOfType(PhreeqcComponent::KINETIC);
}

// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixEquilibrium() const {
  return getMatrixNamesOfType(PhreeqcComponent::EQUILIBRIUM);
}

std::vector<std::string> PhreeqcMatrix::getMatrixImmobileExtensive() const {
  // Intensive entries of the immobile components, by column suffix (see the
  // names() functions of the wrappers). Equilibrium phases and kinetics are
  // extensive only for the moles (_eq, _kin); _si and the kinetic
  // parameters _pN are not.
  static const std::vector<std::string> exchange_intensive = {
      "_la", "_phase_proportion", "_formular_z"};
  static const std::vector<std::string> surface_intensive = {"_la", "_area"};

  auto ends_with_any = [](const std::string &name,
                          const std::vector<std::string> &suffixes) {
    return std::any_of(
        suffixes.begin(), suffixes.end(),
        [&name](const std::string &suffix) { return name.ends_with(suffix); });
  };

  std::vector<std::string> names;
  for (const auto &name : getMatrixNamesOfType(PhreeqcComponent::EXCHANGE)) {
    if (!ends_with_any(name, exchange_intensive)) {
      names.push_back(name);
    }
  }
  for (const auto &name : getMatrixKinetics()) {
    if (name.ends_with("_kin")) {
      names.push_back(name);
    }
  }
  for (const auto &name : getMatrixEquilibrium()) {
    if (name.ends_with("_eq")) {
      names.push_back(name);
    }
  }
  for (const auto &name :
       getMatrixNamesOfType(PhreeqcComponent::SURFACE_COMPS)) {
    if (!ends_with_any(name, surface_intensive)) {
      names.push_back(name);
    }
  }
  return names;
}

// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixMinerals() const {

    std::vector<std::string> nam_eq  = this->getMatrixEquilibrium();
    std::vector<std::string> nam_kin = this->getMatrixKinetics();
    std::vector<std::string> all;
    for (auto name : nam_eq){
	if (name.ends_with("_eq")) {
	    all.push_back(name);
	}
    }
    for (auto name : nam_kin){
	if (name.ends_with("_kin")) {
	    all.push_back(name);
	}
    }
    
    return all;
}

// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixTransported() const {
  using OE = PhreeqcMatrix::OptionalEssentials;

  // Build the removal list dynamically based on optional essentials
  std::vector<std::string> to_remove = {"tc", "patm"};

  // Add optional essentials that are enabled to the removal list
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::SolVol))
    to_remove.push_back("SolVol");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::pH))
    to_remove.push_back("pH");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::pe))
    to_remove.push_back("pe");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::MassH2O))
    to_remove.push_back("MassH2O");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::Viscosity))
    to_remove.push_back("Viscosity");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::Density))
    to_remove.push_back("Density");

  std::vector<std::string> names;
  auto sols = this->getSolutionNames();
  for (const auto &name : sols) {
    if (std::find(to_remove.begin(), to_remove.end(), name) == to_remove.end()) {
      names.push_back(name);
    }
  }

  return names;
}

// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixOutOnly() const {
  using OE = PhreeqcMatrix::OptionalEssentials;

  // Build the list dynamically based on optional essentials
  std::vector<std::string> ret = {"tc", "patm"};

  // Add only the enabled optional essentials
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::SolVol))
    ret.push_back("SolVol");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::pH))
    ret.push_back("pH");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::pe))
    ret.push_back("pe");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::MassH2O))
    ret.push_back("MassH2O");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::Viscosity))
    ret.push_back("Viscosity");
  if (PhreeqcMatrix::hasFlag(_m_optional_essentials, OE::Density))
    ret.push_back("Density");

  return ret;
}
