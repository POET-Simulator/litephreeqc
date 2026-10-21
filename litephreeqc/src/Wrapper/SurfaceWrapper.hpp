/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#pragma once

#include "Surface.h"
#include "SurfaceCharge.h"
#include "SurfaceComp.h"
#include "WrapperBase.hpp"
#include "phrqtype.h"
#include <cstddef>
#include <memory>
#include <set>
#include <span>
#include <string>
#include <vector>

/*
 * Only the values PHREEQC reads back when it equilibrates the surface are
 * columns (step.cpp: add_surface; mainsubs.cpp: xsurface_save):
 *
 * - per component: la (starting guess), cb (only with -no_edl: it is then
 *   added to the charge balance, otherwise ignored), element totals. The
 *   component moles are never read (saved as 0) and are not a column.
 * - per charge (not with -no_edl): area and grams (parameters), cb, la_psi
 *   (starting guess), and only with a diffuse layer (-donnan,
 *   -diffuse_layer): the diffuse layer totals and, for a fixed thickness,
 *   the water of the layer (recomputed when the thickness is given in Debye
 *   lengths).
 */
class SurfaceWrapper : public WrapperBase {
public:
  SurfaceWrapper(cxxSurface *surf,
                 const std::set<std::string> &solution_primaries,
                 const std::vector<std::string> &comp_formulas,
                 const std::vector<std::string> &charge_names);

  void get(std::span<LDBLE> &surface) const override;

  void set(const std::span<LDBLE> &surface) override;

  static std::vector<std::string>
  names(cxxSurface *surface, const std::set<std::string> &solution_primaries,
        std::vector<std::string> &comp_formulas,
        std::vector<std::string> &charge_names);

private:
  cxxSurface *surface;

  // Which optional values of the surface are read by PHREEQC
  struct Layout {
    bool comp_cb;   // component charge balance (-no_edl)
    bool dl_totals; // diffuse layer totals
    bool dl_water;  // water of the diffuse layer (fixed thickness)

    explicit Layout(const cxxSurface *surface);
  };

  class SurfaceCompWrapper : public WrapperBase {
  public:
    SurfaceCompWrapper(cxxSurfaceComp &comp, const Layout &layout);

    void get(std::span<LDBLE> &surface) const;

    void set(const std::span<LDBLE> &surface);

    static std::vector<std::string> names(cxxSurfaceComp &comp,
                                          const Layout &layout);

  private:
    cxxSurfaceComp &surface_comp;
    const bool with_cb;
    const std::size_t num_not_totals;

    std::vector<std::string> total_names;
  };

  class SurfaceChargeWrapper : public WrapperBase {
  public:
    SurfaceChargeWrapper(cxxSurfaceCharge &charge,
                         const std::set<std::string> &solution_primaries,
                         const Layout &layout);
    void get(std::span<LDBLE> &surface) const;

    void set(const std::span<LDBLE> &surface);

    static std::vector<std::string>
    names(cxxSurfaceCharge *charge, const std::set<std::string> &primaries,
          const Layout &layout);

  private:
    cxxSurfaceCharge &surface_charge;
    const bool with_water;
    const bool with_totals;
    const std::size_t num_not_totals;

    std::set<std::string> primaries;
  };

  std::vector<std::unique_ptr<SurfaceCompWrapper>> surface_comps;
  std::vector<std::unique_ptr<SurfaceChargeWrapper>> surface_charges;
};