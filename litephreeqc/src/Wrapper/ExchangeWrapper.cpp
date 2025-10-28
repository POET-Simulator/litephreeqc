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

#include "ExchangeWrapper.hpp"

ExchangeWrapper::ExchangeWrapper(cxxExchange *exch,
                                 const std::vector<std::string> &exchanger)
    : exchange(exch) {

  for (const auto &comp_name : exchanger) {
    auto it = std::find_if(exchange->Get_exchange_comps().begin(),
                           exchange->Get_exchange_comps().end(),
                           [&](const cxxExchComp &comp) {
                             return comp.Get_formula() == comp_name;
                           });

    if (it == exchange->Get_exchange_comps().end()) {
      throw std::runtime_error(
          "Exchange component not found in Phreeqc variables");
    }

    exchange_comps.push_back(std::make_unique<ExchangeCompWrapper>(*it));

    num_elements += exchange_comps.back()->size();
  }

  // const std::size_t defined_comps = exchange->Get_exchange_comps().size();

  // auto header_it = remaining_field_header.begin();

  // while (header_it != remaining_field_header.end() &&
  //        exchange_comps.size() < defined_comps) {
  //   const std::string formular = *header_it;

  //   auto it = std::find_if(exchange->Get_exchange_comps().begin(),
  //                          exchange->Get_exchange_comps().end(),
  //                          [&](const cxxExchComp &comp) {
  //                            return comp.Get_formula() == formular;
  //                          });

  //   if (it != exchange->Get_exchange_comps().end()) {
  //     const size_t i = this->exchange_comps.size();

  //     exchange_comps.push_back(std::make_unique<ExchangeCompWrapper>(*it));
  //     header_it += this->exchange_comps[i]->size();
  //     num_elements += this->exchange_comps[i]->size();
  //     continue;
  //   }

  //   header_it++;
  // }

  // if (exchange_comps.size() != defined_comps) {
  //   throw std::runtime_error(
  //       "Not all exchange components found in Phreeqc variables");
  // }
}

void ExchangeWrapper::get(std::span<LDBLE> &exchange) const {
  std::size_t offset = 0;

  for (const auto &comp : exchange_comps) {
    std::span<LDBLE> comp_span = exchange.subspan(offset, comp->size());
    comp->get(comp_span);
    offset += comp->size();
  }
}

void ExchangeWrapper::set(const std::span<LDBLE> &exchange) {
  std::size_t offset = 0;

  for (const auto &comp : exchange_comps) {
    std::span<LDBLE> comp_span = exchange.subspan(offset, comp->size());
    comp->set(comp_span);
    offset += comp->size();
  }
}

std::vector<std::string>
ExchangeWrapper::names(cxxExchange *exchange,
                       std::vector<std::string> &exchange_formulas) {
  std::vector<std::string> e_names;

  for (const auto &comp : exchange->Get_exchange_comps()) {

    const std::string &formular = comp.Get_formula();
    exchange_formulas.push_back(formular);
    e_names.push_back(formular);

    e_names.push_back(formular + "_cb");
    e_names.push_back(formular + "_la");
    e_names.push_back(formular + "_phase_proportion");
    e_names.push_back(formular + "_formular_z");

    for (const auto &total : comp.Get_totals()) {
      if (total.first == formular) {
        continue;
      }
      e_names.push_back(total.first + formular);
    }
  }

  return e_names;
}