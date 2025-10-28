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
#include "SurfaceComp.h"

SurfaceWrapper::SurfaceWrapper(cxxSurface *surf,
                               const std::set<std::string> &solution_primaries,
                               const std::vector<std::string> &comp_formulas,
                               const std::vector<std::string> &charge_names)
    : surface(surf) {

  for (const auto &comp_name : comp_formulas) {
    auto it = std::find_if(surface->Get_surface_comps().begin(),
                           surface->Get_surface_comps().end(),
                           [&](const cxxSurfaceComp &comp) {
                             return comp.Get_formula() == comp_name;
                           });

    if (it == surface->Get_surface_comps().end()) {
      throw std::runtime_error(
          "Surface component not found in Phreeqc variables");
    }

    surface_comps.push_back(std::make_unique<SurfaceCompWrapper>(*it));

    num_elements += surface_comps.back()->size();
  }

  for (const auto &charge_name : charge_names) {
    auto it = std::find_if(surface->Get_surface_charges().begin(),
                           surface->Get_surface_charges().end(),
                           [&](const cxxSurfaceCharge &charge) {
                             return charge.Get_name() == charge_name;
                           });

    if (it == surface->Get_surface_charges().end()) {
      throw std::runtime_error("Surface charge not found in Phreeqc variables");
    }

    surface_charges.push_back(
        std::make_unique<SurfaceChargeWrapper>(*it, solution_primaries));

    num_elements += surface_charges.back()->size();
  }
}

void SurfaceWrapper::get(std::span<LDBLE> &surface) const {
  std::size_t offset = 0;

  for (const auto &comp : surface_comps) {
    std::span<LDBLE> comp_span = surface.subspan(offset, comp->size());
    comp->get(comp_span);
    offset += comp->size();
  }

  for (const auto &charge : surface_charges) {
    std::span<LDBLE> charge_span = surface.subspan(offset, charge->size());
    charge->get(charge_span);
    offset += charge->size();
  }
}

void SurfaceWrapper::set(const std::span<LDBLE> &surface) {
  std::size_t offset = 0;

  for (const auto &comp : surface_comps) {
    std::span<LDBLE> comp_span = surface.subspan(offset, comp->size());
    comp->set(comp_span);
    offset += comp->size();
  }

  for (const auto &charge : surface_charges) {
    std::span<LDBLE> charge_span = surface.subspan(offset, charge->size());
    charge->set(charge_span);
    offset += charge->size();
  }
}

std::vector<std::string>
SurfaceWrapper::names(cxxSurface *surface,
                      const std::set<std::string> &solution_primaries,
                      std::vector<std::string> &comp_formulas,
                      std::vector<std::string> &charge_names) {
  std::vector<std::string> names;

  for (auto &comp : surface->Get_surface_comps()) {
    comp_formulas.push_back(comp.Get_formula());
    const auto charge_names = SurfaceCompWrapper::names(comp);
    names.insert(names.end(), charge_names.begin(), charge_names.end());
  }

  for (auto &charge : surface->Get_surface_charges()) {
    charge_names.push_back(charge.Get_name());
    const auto comp_names =
        SurfaceChargeWrapper::names(&charge, solution_primaries);
    names.insert(names.end(), comp_names.begin(), comp_names.end());
  }

  return names;
}