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

// area, grams, cb, [mw], la_psi, [diffuse layer totals]
SurfaceWrapper::SurfaceChargeWrapper::SurfaceChargeWrapper(
    cxxSurfaceCharge &charge, const std::set<std::string> &solution_primaries,
    const Layout &layout)
    : surface_charge(charge), with_water(layout.dl_water),
      with_totals(layout.dl_totals), num_not_totals(layout.dl_water ? 5 : 4),
      primaries(solution_primaries) {
  num_elements =
      num_not_totals + (with_totals ? solution_primaries.size() : 0);
}

void SurfaceWrapper::SurfaceChargeWrapper::get(
    std::span<LDBLE> &surface) const {
  std::size_t index = 0;
  surface[index++] = this->surface_charge.Get_specific_area();
  surface[index++] = this->surface_charge.Get_grams();
  surface[index++] = this->surface_charge.Get_charge_balance();
  if (this->with_water) {
    surface[index++] = this->surface_charge.Get_mass_water();
  }
  surface[index++] = this->surface_charge.Get_la_psi();

  if (!this->with_totals) {
    return;
  }

  const auto dl_map = this->surface_charge.Get_diffuse_layer_totals();

  for (const auto &name : this->primaries) {
    const auto dl_it = dl_map.find(name);

    surface[index++] = dl_it == dl_map.end() ? 0 : dl_it->second;
  }
}

void SurfaceWrapper::SurfaceChargeWrapper::set(
    const std::span<LDBLE> &surface) {
  std::size_t index = 0;
  this->surface_charge.Set_specific_area(surface[index++]);
  this->surface_charge.Set_grams(surface[index++]);
  this->surface_charge.Set_charge_balance(surface[index++]);
  if (this->with_water) {
    this->surface_charge.Set_mass_water(surface[index++]);
  }
  this->surface_charge.Set_la_psi(surface[index++]);

  if (!this->with_totals) {
    return;
  }

  auto &dl_map = this->surface_charge.Get_diffuse_layer_totals();
  dl_map.clear();

  for (const auto &name : this->primaries) {
    const double value = surface[index++];

    if (value == 0) {
      continue;
    }

    dl_map[name] = value;
  }
}

std::vector<std::string> SurfaceWrapper::SurfaceChargeWrapper::names(
    cxxSurfaceCharge *charge, const std::set<std::string> &primaries,
    const Layout &layout) {
  std::vector<std::string> names;

  const std::string &charge_name = charge->Get_name();
  names.push_back(charge_name + "_area");
  names.push_back(charge_name + "_grams");
  names.push_back(charge_name + "_cb");
  if (layout.dl_water) {
    names.push_back(charge_name + "_mw");
  }
  names.push_back(charge_name + "_la");

  if (!layout.dl_totals) {
    return names;
  }

  for (const auto &name : primaries) {
    names.push_back(charge_name + "_tot_" + name);
  }

  return names;
}
