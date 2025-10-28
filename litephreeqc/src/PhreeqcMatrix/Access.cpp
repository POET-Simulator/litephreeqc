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
       component <=
       static_cast<std::size_t>(PhreeqcMatrix::PhreeqcComponent::SURFACE_COMPS);
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
std::vector<std::string> PhreeqcMatrix::getMatrixKinetics() const {
    std::vector<std::string> names;
    
    auto n = this->getIds().size();
    for (auto i = 0; i<n; ++i) {
	auto pqc_kinnames = this->getKineticsNames(i);
	for (auto nam : pqc_kinnames ) {
	    for (auto mat_name : this->get().names){
		if (mat_name.starts_with(nam)) {
		    // check if we already have this mat_name
		    if (std::find(names.begin(), names.end(), mat_name) == names.end()) {
			names.push_back(mat_name);
		    }
		}
	    }
	}
    }
    return names;
}


// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixEquilibrium() const {

    std::vector<std::string> names;
    std::vector<std::string> mat_names = this->get().names;
    auto n = this->getIds().size();
    for (auto i = 0; i<n; ++i) {
	auto pqc_eqnames = this->getEquilibriumNames(i);
	for (auto nam : pqc_eqnames ) {
	    for (auto mat_name : mat_names){
		if (mat_name.starts_with(nam)) {
		    // check if we already have this mat_name
		    if (std::find(names.begin(), names.end(), mat_name) == names.end()) {
			names.push_back(mat_name);
		    }
		}
	    }
	}
    }
    return names;
}

// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixTransported() const {
    std::vector<std::string> names;
    
    const std::vector<std::string> to_remove = {
	"tc", "patm", "SolVol", "pH", "pe"
    };

    // sols contains all solutes; we must remove { tc, patm, SolVol, pH, pe }
    auto sols = this->getSolutionNames();
    for (auto name : sols) {
	if (std::find(to_remove.begin(), to_remove.end(), name) == to_remove.end()) {
	    names.push_back(name);
	}
    }
    
    return names;
}

// MDL
std::vector<std::string> PhreeqcMatrix::getMatrixOutOnly() const {
    // MDL we must append here selected_output / user_punch
    std::vector<std::string> defaultnames = {
	"tc", "patm", "SolVol", "pH", "pe"
    };
    std::vector<std::string> ret;
    for (auto nm : defaultnames) {
	ret.push_back(nm);
    }
    return ret;
}
