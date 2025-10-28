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

SurfaceWrapper::SurfaceChargeWrapper::SurfaceChargeWrapper(
    cxxSurfaceCharge &charge, const std::set<std::string> &solution_primaries)
    : surface_charge(charge), primaries(solution_primaries) {
  num_elements = NUM_NOT_TOTALS + solution_primaries.size();
}

void SurfaceWrapper::SurfaceChargeWrapper::get(
    std::span<LDBLE> &surface) const {
  surface[0] = this->surface_charge.Get_specific_area();
  surface[1] = this->surface_charge.Get_grams();
  surface[2] = this->surface_charge.Get_charge_balance();
  surface[3] = this->surface_charge.Get_mass_water();
  surface[4] = this->surface_charge.Get_la_psi();

  const auto dl_map = this->surface_charge.Get_diffuse_layer_totals();

  for (auto it = this->primaries.begin(); it != this->primaries.end(); it++) {
    const std::size_t index =
        NUM_NOT_TOTALS + std::distance(this->primaries.begin(), it);

    const auto dl_it = dl_map.find(*it);

    if (dl_it == dl_map.end()) {
      surface[index] = 0;
      continue;
    }
    surface[index] = dl_it->second;
  }
}

void SurfaceWrapper::SurfaceChargeWrapper::set(
    const std::span<LDBLE> &surface) {
  this->surface_charge.Set_specific_area(surface[0]);
  this->surface_charge.Set_grams(surface[1]);
  this->surface_charge.Set_charge_balance(surface[2]);
  this->surface_charge.Set_mass_water(surface[3]);
  this->surface_charge.Set_la_psi(surface[4]);

  auto &dl_map = this->surface_charge.Get_diffuse_layer_totals();
  dl_map.clear();

  for (auto it = this->primaries.begin(); it != this->primaries.end(); it++) {
    const std::size_t index =
        NUM_NOT_TOTALS + std::distance(this->primaries.begin(), it);

    if (surface[index] == 0) {
      continue;
    }

    dl_map[*it] = surface[index];
  }
}

std::vector<std::string> SurfaceWrapper::SurfaceChargeWrapper::names(
    cxxSurfaceCharge *charge, const std::set<std::string> &primaries) {
  std::vector<std::string> names;

  const std::string &charge_name = charge->Get_name();
  names.push_back(charge_name + "_area");
  names.push_back(charge_name + "_grams");
  names.push_back(charge_name + "_cb");
  names.push_back(charge_name + "_mw");
  names.push_back(charge_name + "_la");

  for (const auto &name : primaries) {
    names.push_back(charge_name + "_tot_" + name);
  }

  return names;
}