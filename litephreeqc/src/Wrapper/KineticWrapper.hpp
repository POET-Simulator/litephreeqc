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

#include "KineticsComp.h"
#include "WrapperBase.hpp"
#include "cxxKinetics.h"
#include <memory>
#include <vector>

class KineticWrapper : public WrapperBase {

public:
  KineticWrapper(cxxKinetics *kinetics,
                 const std::vector<std::string> &kin_comps);

  void get(std::span<LDBLE> &data) const override;
  void set(const std::span<LDBLE> &data) override;

  static std::vector<std::string> names(cxxKinetics *kinetics,
                                        std::vector<std::string> &kin_comps);

private:
  cxxKinetics *kinetics;

  class KineticCompWrapper : public WrapperBase {
  public:
    KineticCompWrapper(cxxKineticsComp &comp);

    void get(std::span<LDBLE> &data) const override;
    void set(const std::span<LDBLE> &data) override;

    static std::vector<std::string> names(const cxxKineticsComp &comp);

  private:
    cxxKineticsComp &kin_comp;
  };

  std::vector<std::unique_ptr<KineticCompWrapper>> kinetic_comps;
};