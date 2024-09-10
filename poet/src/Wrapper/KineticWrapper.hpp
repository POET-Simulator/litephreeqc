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