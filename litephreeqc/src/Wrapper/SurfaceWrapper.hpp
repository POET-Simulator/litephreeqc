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

  class SurfaceCompWrapper : public WrapperBase {
  public:
    SurfaceCompWrapper(cxxSurfaceComp &comp);

    void get(std::span<LDBLE> &surface) const;

    void set(const std::span<LDBLE> &surface);

    static std::vector<std::string> names(cxxSurfaceComp &comp);

  private:
    cxxSurfaceComp &surface_comp;
    static constexpr std::size_t NUM_NOT_TOTALS = 3;

    std::vector<std::string> total_names;
  };

  class SurfaceChargeWrapper : public WrapperBase {
  public:
    SurfaceChargeWrapper(cxxSurfaceCharge &charge,
                         const std::set<std::string> &solution_primaries);
    void get(std::span<LDBLE> &surface) const;

    void set(const std::span<LDBLE> &surface);

    static std::vector<std::string>
    names(cxxSurfaceCharge *charge, const std::set<std::string> &primaries);

  private:
    cxxSurfaceCharge &surface_charge;
    static constexpr std::size_t NUM_NOT_TOTALS = 5;

    std::set<std::string> primaries;
  };

  std::vector<std::unique_ptr<SurfaceCompWrapper>> surface_comps;
  std::vector<std::unique_ptr<SurfaceChargeWrapper>> surface_charges;
};