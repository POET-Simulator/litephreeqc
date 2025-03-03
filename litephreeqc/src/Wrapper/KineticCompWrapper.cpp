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

#include "KineticWrapper.hpp"
#include <cstddef>
#include <string>
#include <vector>

KineticWrapper::KineticCompWrapper::KineticCompWrapper(cxxKineticsComp &comp)
    : kin_comp(comp) {
  this->num_elements = kin_comp.Get_d_params().size() + 1;
}

void KineticWrapper::KineticCompWrapper::get(std::span<LDBLE> &data) const {
  data[0] = this->kin_comp.Get_m();

  std::size_t next_index = 1;
  for (const auto &param : this->kin_comp.Get_d_params()) {
    data[next_index++] = param;
  }
}
void KineticWrapper::KineticCompWrapper::set(const std::span<LDBLE> &data) {
  this->kin_comp.Set_m(data[0]);

  std::size_t next_index = 1;
  for (auto &param : this->kin_comp.Get_d_params()) {
    param = data[next_index++];
  }
}

std::vector<std::string>
KineticWrapper::KineticCompWrapper::names(const cxxKineticsComp &comp) {
  std::vector<std::string> names;

  const std::string &comp_name = comp.Get_rate_name();
  names.push_back(comp_name + "_kin");

  for (std::size_t i = 0; i < comp.Get_d_params().size(); i++) {
    names.push_back(comp_name + "_p" + std::to_string(i + 1));
  }

  return names;
}
