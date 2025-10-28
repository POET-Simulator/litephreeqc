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

EquilibriumWrapper::EquilibriumCompWrapper::EquilibriumCompWrapper(
    cxxPPassemblageComp &comp_)
    : comp(comp_) {
  this->num_elements = 2;
}

void EquilibriumWrapper::EquilibriumCompWrapper::get(
    std::span<LDBLE> &data) const {
  data[0] = this->comp.Get_moles();
  data[1] = this->comp.Get_si();
}

void EquilibriumWrapper::EquilibriumCompWrapper::set(
    const std::span<LDBLE> &data) {
  this->comp.Set_moles(data[0]);
  this->comp.Set_si(data[1]);
}

std::vector<std::string> EquilibriumWrapper::EquilibriumCompWrapper::names(
    const cxxPPassemblageComp &comp) {
  std::vector<std::string> names;

  const std::string &comp_name = comp.Get_name();
  names.push_back(comp_name + "_eq");
  names.push_back(comp_name + "_si");

  return names;
}