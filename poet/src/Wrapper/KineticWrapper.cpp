#include "KineticWrapper.hpp"
#include <vector>

KineticWrapper::KineticWrapper(cxxKinetics *kinetics,
                               const std::vector<std::string> &kin_comps)
    : kinetics(kinetics) {

  for (const auto &kin_name : kin_comps) {
    auto it = std::find_if(kinetics->Get_kinetics_comps().begin(),
                           kinetics->Get_kinetics_comps().end(),
                           [&](const cxxKineticsComp &comp) {
                             return comp.Get_rate_name() == kin_name;
                           });

    if (it == kinetics->Get_kinetics_comps().end()) {
      throw std::runtime_error(
          "Kinetic component not found in Phreeqc variables");
    }

    kinetic_comps.push_back(std::make_unique<KineticCompWrapper>(*it));

    num_elements += kinetic_comps.back()->size();
  }
}

void KineticWrapper::get(std::span<LDBLE> &data) const {
  std::size_t offset = 0;

  for (const auto &comp : kinetic_comps) {
    std::span<LDBLE> comp_span = data.subspan(offset, comp->size());

    comp->get(comp_span);

    offset += comp->size();
  }
}

void KineticWrapper::set(const std::span<LDBLE> &data) {
  std::size_t offset = 0;

  for (const auto &comp : kinetic_comps) {
    std::span<LDBLE> comp_span = data.subspan(offset, comp->size());

    comp->set(comp_span);

    offset += comp->size();
  }
}

std::vector<std::string>
KineticWrapper::names(cxxKinetics *kinetics,
                      std::vector<std::string> &kin_comps) {
  std::vector<std::string> names;

  for (const auto &comp : kinetics->Get_kinetics_comps()) {
    std::vector<std::string> comp_names = KineticCompWrapper::names(comp);
    names.insert(names.end(), comp_names.begin(), comp_names.end());

    kin_comps.push_back(comp.Get_rate_name());
  }

  // for (const auto &kin_name : kin_comps) {
  //   auto it = std::find_if(kinetics->Get_kinetics_comps().begin(),
  //                          kinetics->Get_kinetics_comps().end(),
  //                          [&](const cxxKineticsComp &comp) {
  //                            return comp.Get_rate_name() == kin_name;
  //                          });

  //   if (it == kinetics->Get_kinetics_comps().end()) {
  //     throw std::runtime_error(
  //         "Kinetic component not found in Phreeqc variables");
  //   }

  //   names.push_back(it->Get_rate_name());
  // }

  return names;
}