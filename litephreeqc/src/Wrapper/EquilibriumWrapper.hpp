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

#include "PPassemblage.h"
#include "WrapperBase.hpp"
#include <memory>
#include <vector>

class EquilibriumWrapper : public WrapperBase {
public:
  EquilibriumWrapper(cxxPPassemblage *ppassemblage,
                     const std::vector<std::string> &ppassemblage_names);

  void get(std::span<LDBLE> &data) const override;

  void set(const std::span<LDBLE> &data) override;

  static std::vector<std::string>
  names(const cxxPPassemblage *ppassemblage,
        std::vector<std::string> &ppassemblage_names);

private:
  cxxPPassemblage *ppassemblage;

  class EquilibriumCompWrapper : public WrapperBase {
  public:
    EquilibriumCompWrapper(cxxPPassemblageComp &comp);

    void get(std::span<LDBLE> &data) const override;

    void set(const std::span<LDBLE> &data) override;

    static std::vector<std::string> names(const cxxPPassemblageComp &comp);

    std::string getCompName() const { return this->comp.Get_name(); }

  private:
    cxxPPassemblageComp &comp;
  };

  std::vector<std::unique_ptr<EquilibriumCompWrapper>> equilibrium_comps;
  const std::vector<std::string> ppassemblage_order;
};