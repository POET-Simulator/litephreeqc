#pragma once

#include "PPassemblage.h"
#include "WrapperBase.hpp"
#include <memory>

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

  private:
    cxxPPassemblageComp &comp;
  };

  std::vector<std::unique_ptr<EquilibriumCompWrapper>> equilibrium_comps;
};