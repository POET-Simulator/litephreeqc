/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include "SurfaceWrapper.hpp"

// la, [cb], totals
SurfaceWrapper::SurfaceCompWrapper::SurfaceCompWrapper(cxxSurfaceComp &comp,
                                                       const Layout &layout)
    : surface_comp(comp), with_cb(layout.comp_cb),
      num_not_totals(layout.comp_cb ? 2 : 1) {
  num_elements = num_not_totals + comp.Get_totals().size();

  for (const auto &[name, value] : comp.Get_totals()) {
    total_names.push_back(name);
  }
}

void SurfaceWrapper::SurfaceCompWrapper::get(std::span<LDBLE> &surface) const {
  surface[0] = this->surface_comp.Get_la();
  if (this->with_cb) {
    surface[1] = this->surface_comp.Get_charge_balance();
  }

  const auto &totals = this->surface_comp.Get_totals();
  for (std::size_t i = 0; i < this->total_names.size(); i++) {
    surface[num_not_totals + i] = totals.at(this->total_names[i]);
  }
}

void SurfaceWrapper::SurfaceCompWrapper::set(const std::span<LDBLE> &surface) {
  this->surface_comp.Set_la(surface[0]);
  if (this->with_cb) {
    this->surface_comp.Set_charge_balance(surface[1]);
  }

  auto &totals = this->surface_comp.Get_totals();
  totals.clear();

  for (std::size_t i = 0; i < this->total_names.size(); i++) {
    const std::size_t index = num_not_totals + i;

    totals[this->total_names[i]] = surface[index];
  }
}

std::vector<std::string>
SurfaceWrapper::SurfaceCompWrapper::names(cxxSurfaceComp &comp,
                                          const Layout &layout) {
  std::vector<std::string> names;

  const std::string &phase_name = comp.Get_formula();
  names.push_back(phase_name + "_la");
  if (layout.comp_cb) {
    names.push_back(phase_name + "_cb");
  }

  for (const auto &tot : comp.Get_totals()) {
    names.push_back(phase_name + "_" + tot.first);
  }

  return names;
}
