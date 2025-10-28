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

#include "EquilibriumWrapper.hpp"

#include <algorithm>
#include <memory>

EquilibriumWrapper::EquilibriumWrapper(
    cxxPPassemblage *ppassemblage,
    const std::vector<std::string> &ppassemblage_names)
    : ppassemblage(ppassemblage), ppassemblage_order(ppassemblage_names) {
  for (const auto &ppassemblage_name : ppassemblage_names) {
    auto it = std::find_if(
        ppassemblage->Get_pp_assemblage_comps().begin(),
        ppassemblage->Get_pp_assemblage_comps().end(),
        [&](const auto &comp) { return comp.first == ppassemblage_name; });

    if (it == ppassemblage->Get_pp_assemblage_comps().end()) {
      throw std::runtime_error(
          "Equilibrium component not found in Phreeqc variables");
    }

    equilibrium_comps.push_back(
        std::make_unique<EquilibriumCompWrapper>(it->second));

    num_elements += equilibrium_comps.back()->size();
  }
}

void EquilibriumWrapper::get(std::span<LDBLE> &data) const {
  std::size_t offset = 0;

  for (const auto &comp_name : ppassemblage_order) {
    auto it = std::find_if(
        equilibrium_comps.begin(), equilibrium_comps.end(),
        [&](const auto &comp) { return comp->getCompName() == comp_name; });

    std::span<LDBLE> comp_span = data.subspan(offset, (*it)->size());

    (*it)->get(comp_span);

    offset += (*it)->size();
  }

  // for (const auto &comp : equilibrium_comps) {
  //   std::span<LDBLE> comp_span = data.subspan(offset, comp->size());

  //   comp->get(comp_span);

  //   offset += comp->size();
  // }
}

void EquilibriumWrapper::set(const std::span<LDBLE> &data) {
  std::size_t offset = 0;

  for (const auto &comp_name : ppassemblage_order) {
    auto it = std::find_if(
        equilibrium_comps.begin(), equilibrium_comps.end(),
        [&](const auto &comp) { return comp->getCompName() == comp_name; });

    std::span<LDBLE> comp_span = data.subspan(offset, (*it)->size());

    (*it)->set(comp_span);

    offset += (*it)->size();
  }

  // for (const auto &comp : equilibrium_comps) {
  //   std::span<LDBLE> comp_span = data.subspan(offset, comp->size());

  //   comp->set(comp_span);

  //   offset += comp->size();
  // }
}

std::vector<std::string>
EquilibriumWrapper::names(const cxxPPassemblage *ppassemblage,
                          std::vector<std::string> &ppassemblage_names) {
  std::vector<std::string> names;

  for (const auto &comp : ppassemblage->Get_pp_assemblage_comps()) {
    std::vector<std::string> comp_names =
        EquilibriumCompWrapper::names(comp.second);
    names.insert(names.end(), comp_names.begin(), comp_names.end());

    ppassemblage_names.push_back(comp.first);
  }

  return names;
}