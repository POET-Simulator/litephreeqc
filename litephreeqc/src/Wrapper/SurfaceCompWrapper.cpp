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

#include "SurfaceWrapper.hpp"

SurfaceWrapper::SurfaceCompWrapper::SurfaceCompWrapper(cxxSurfaceComp &comp)
    : surface_comp(comp) {
  num_elements = NUM_NOT_TOTALS + comp.Get_totals().size();

  for (const auto &[name, value] : comp.Get_totals()) {
    total_names.push_back(name);
  }
}

void SurfaceWrapper::SurfaceCompWrapper::get(std::span<LDBLE> &surface) const {
  surface[0] = this->surface_comp.Get_moles();
  surface[1] = this->surface_comp.Get_la();
  surface[2] = this->surface_comp.Get_charge_balance();

  const auto &totals = this->surface_comp.Get_totals();
  for (std::size_t i = 0; i < this->total_names.size(); i++) {
    surface[NUM_NOT_TOTALS + i] = totals.at(this->total_names[i]);
  }
}

void SurfaceWrapper::SurfaceCompWrapper::set(const std::span<LDBLE> &surface) {
  this->surface_comp.Set_moles(surface[0]);
  this->surface_comp.Set_la(surface[1]);
  this->surface_comp.Set_charge_balance(surface[2]);

  auto &totals = this->surface_comp.Get_totals();
  totals.clear();

  for (std::size_t i = 0; i < this->total_names.size(); i++) {
    const std::size_t index = NUM_NOT_TOTALS + i;

    totals[this->total_names[i]] = surface[index];
  }
}

std::vector<std::string>
SurfaceWrapper::SurfaceCompWrapper::names(cxxSurfaceComp &comp) {
  std::vector<std::string> names;

  const std::string &phase_name = comp.Get_formula();
  names.push_back(phase_name + "_moles");
  names.push_back(phase_name + "_la");
  names.push_back(phase_name + "_cb");

  for (const auto &tot : comp.Get_totals()) {
    names.push_back(phase_name + "_" + tot.first);
  }

  return names;
}