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

#include "../Wrapper/EquilibriumWrapper.hpp"
#include "../Wrapper/ExchangeWrapper.hpp"
#include "../Wrapper/KineticWrapper.hpp"
#include "../Wrapper/SolutionWrapper.hpp"
#include "../Wrapper/SurfaceWrapper.hpp"

#include <IPhreeqc.hpp>
#include <Phreeqc.h>
#include <Solution.h>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

bool include_h0_o0 = false;
bool with_redox = false;

static std::vector<std::string> dump_solution_names(cxxSolution *solution,
                                                    Phreeqc *phreeqc) {
  std::vector<std::string> placeholder;

  std::vector<std::string> solnames =
      SolutionWrapper::names(solution, include_h0_o0, placeholder, with_redox);

  if (with_redox) {
    solnames = phreeqc->find_all_valence_states(solnames);
  }

  return solnames;
}

template <enum PhreeqcMatrix::PhreeqcComponent comp, class T>
static void
base_add_to_element_vector(const T &wrapper,
                           const std::vector<std::string> &names,
                           std::vector<PhreeqcMatrix::element> &elements) {
  std::vector<double> values(wrapper.size());
  std::span<double> values_span(values);

  wrapper.get(values_span);

  for (std::size_t i = 0; i < names.size(); ++i) {
    elements.push_back({names[i], comp, values[i]});
  }
}

template <enum PhreeqcMatrix::PhreeqcComponent comp, class Wrapper_t, class T>
static void ex_ki_eq_add_to_element_vector(
    T *phreeqc_component, std::vector<PhreeqcMatrix::element> &elements,
    std::vector<PhreeqcMatrix::base_names> &b_names) {
  if (phreeqc_component == NULL) {
    return;
  }

  std::vector<std::string> ctor_names;
  std::vector<std::string> names =
      Wrapper_t::names(phreeqc_component, ctor_names);

  for (const auto &name : ctor_names) {
    b_names.push_back(
        {static_cast<PhreeqcMatrix::base_names::Components>(comp), name});
  }

  Wrapper_t wrapper(phreeqc_component, ctor_names);

  base_add_to_element_vector<comp>(wrapper, names, elements);
}

static void
surface_add_to_element_vector(Phreeqc *phreeqc, cxxSurface *surface,
                              const std::vector<std::string> &solution_names,
                              std::vector<PhreeqcMatrix::element> &elements,
                              std::vector<PhreeqcMatrix::base_names> &b_names,
                              std::set<std::string> &surface_primaries) {

  if (surface == NULL) {
    return;
  }
  // H and O are fixed surface primaries
  surface_primaries.insert("H");
  surface_primaries.insert("O");

  for (std::size_t i = 3; i < solution_names.size(); i++) {
    const auto master_primary =
        phreeqc->master_bsearch_primary(solution_names[i].c_str());
    if (master_primary != NULL) {
      surface_primaries.insert(master_primary->elt->name);
    }
  }

  std::vector<std::string> comp_formulas;
  std::vector<std::string> charge_names;

  std::vector<std::string> surf_names = SurfaceWrapper::names(
      surface, surface_primaries, comp_formulas, charge_names);

  for (const auto &comp : comp_formulas) {
    b_names.push_back({static_cast<PhreeqcMatrix::base_names::Components>(
                           PhreeqcMatrix::base_names::Components::SURACE_COMP),
                       comp});
  }

  for (const auto &charge : charge_names) {
    b_names.push_back(
        {static_cast<PhreeqcMatrix::base_names::Components>(
             PhreeqcMatrix::base_names::Components::SURFACE_CHARGE),
         charge});
  }

  SurfaceWrapper wrapper(surface, surface_primaries, comp_formulas,
                         charge_names);

  base_add_to_element_vector<PhreeqcMatrix::PhreeqcComponent::SURFACE_COMPS>(
      wrapper, surf_names, elements);
}

static std::pair<std::vector<PhreeqcMatrix::element>,
                 std::vector<PhreeqcMatrix::base_names>>
create_vector_from_phreeqc(Phreeqc *phreeqc, int id,
                           const std::vector<std::string> &solution_names,
                           std::set<std::string> &surface_primaries) {
  std::vector<PhreeqcMatrix::element> elements;
  std::vector<PhreeqcMatrix::base_names> b_names;

  // Solution
  SolutionWrapper sol_wrapper(
      Utilities::Rxn_find(phreeqc->Get_Rxn_solution_map(), id), solution_names,
      with_redox);

  base_add_to_element_vector<PhreeqcMatrix::PhreeqcComponent::SOLUTION>(
      sol_wrapper, solution_names, elements);

  // keep track if exchange is present
  const size_t vec_size_before_exchange = elements.size();

  // Exchange
  ex_ki_eq_add_to_element_vector<PhreeqcMatrix::PhreeqcComponent::EXCHANGE,
                                 ExchangeWrapper>(
      Utilities::Rxn_find(phreeqc->Get_Rxn_exchange_map(), id), elements,
      b_names);

  const bool has_exchange = elements.size() > vec_size_before_exchange;

  // Kinetic
  ex_ki_eq_add_to_element_vector<PhreeqcMatrix::PhreeqcComponent::KINETIC,
                                 KineticWrapper>(
      Utilities::Rxn_find(phreeqc->Get_Rxn_kinetics_map(), id), elements,
      b_names);
  // Equilibrium
  ex_ki_eq_add_to_element_vector<PhreeqcMatrix::PhreeqcComponent::EQUILIBRIUM,
                                 EquilibriumWrapper>(
      Utilities::Rxn_find(phreeqc->Get_Rxn_pp_assemblage_map(), id), elements,
      b_names);

  // Surface
  surface_add_to_element_vector(
      phreeqc, Utilities::Rxn_find(phreeqc->Get_Rxn_surface_map(), id),
      solution_names, elements, b_names, surface_primaries);

  // substitute exchange names
  if (has_exchange) {
    for (auto &currentElement : elements) {
      if (currentElement.type == PhreeqcMatrix::PhreeqcComponent::EXCHANGE) {
        for (const auto &species : phreeqc->Get_species_list()) {
          const std::string &species_name = species.s->name;
          // check if `name` is a prefix of `species_name`
          if (species_name.compare(0, currentElement.name.size(),
                                   currentElement.name) == 0) {
            currentElement.name = species_name;
            break;
          }
        }
      }
    }
  }

  return {elements, b_names};
}

static std::vector<std::string> find_all_solutions(Phreeqc *phreeqc) {
  std::vector<std::vector<std::string>> all_names;

  for (auto &[id, solution] : phreeqc->Get_Rxn_solution_map()) {
    if (id < 0) {
      continue;
    }
    std::vector<std::string> names = dump_solution_names(&solution, phreeqc);
    all_names.push_back(names);
  }

  std::vector<std::string> union_names;

  for (const auto &vec : all_names) {
    std::vector<std::string> result;
    std::set_union(union_names.begin(), union_names.end(), vec.begin(),
                   vec.end(), std::back_inserter(result));

    union_names = result;
  }

  return union_names;
}

void PhreeqcMatrix::initialize() {

  Phreeqc *phreeqc = this->_m_pqc->GetPhreeqcPtr();

  if (!phreeqc->Get_Rxn_surface_map().empty()) {
    this->_m_pqc->RunString("RUN_CELLS\n-cells 1\nEND");
  }

  include_h0_o0 = this->_m_with_h0_o0;
  with_redox = this->_m_with_redox;

  std::vector<std::string> solutions = find_all_solutions(phreeqc);

  for (auto &[id, solution] : phreeqc->Get_Rxn_solution_map()) {
    if (id < 0) {
      continue;
    }
    const auto &[elements, base_names] = create_vector_from_phreeqc(
        phreeqc, id, solutions, this->_m_surface_primaries);

    _m_map[id] = elements;
    _m_internal_names[id] = base_names;
  }
}

void PhreeqcMatrix::remove_NaNs() {
  std::vector<std::vector<std::string>> elements_to_remove;

  for (const auto &[id, elements] : _m_map) {
    elements_to_remove.push_back({});
    std::vector<std::string> &curr_to_remove = elements_to_remove.back();
    for (const auto &element : elements) {
      if (std::isnan(element.value)) {
        curr_to_remove.push_back(element.name);
      }
    }
  }

  // find intersection of all elements to remove
  std::vector<std::string> intersection = elements_to_remove[0];

  for (const auto &vec : elements_to_remove) {
    std::set_intersection(intersection.begin(), intersection.end(), vec.begin(),
                          vec.end(), std::back_inserter(intersection));
  }

  for (auto &[id, elements] : _m_map) {
    for (const auto &to_remove_element : intersection) {
      elements.erase(std::remove_if(elements.begin(), elements.end(),
                                    [to_remove_element](const element &el) {
                                      return el.name == to_remove_element;
                                    }),
                     elements.end());
    }
  }
}