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

ExchangeWrapper::ExchangeCompWrapper::ExchangeCompWrapper(cxxExchComp &comp)
    : exch_comp(comp) {
  const std::size_t totals_size = exch_comp.Get_totals().size() - 1;
  this->num_elements = totals_size + NUM_NOT_TOTALS;
}

void ExchangeWrapper::ExchangeCompWrapper::get(
    std::span<LDBLE> &exchange) const {

  std::size_t exch_index = this->NUM_NOT_TOTALS;

  exchange[0] = exch_comp.Get_totals().find(exch_comp.Get_formula())->second;
  exchange[1] = exch_comp.Get_charge_balance();
  exchange[2] = exch_comp.Get_la();
  exchange[3] = exch_comp.Get_phase_proportion();
  exchange[4] = exch_comp.Get_formula_z();

  for (const auto &[name, value] : exch_comp.Get_totals()) {
    if (name != exch_comp.Get_formula()) {
      exchange[exch_index++] = value;
    }
  }
}

void ExchangeWrapper::ExchangeCompWrapper::set(
    const std::span<LDBLE> &exchange) {

  std::size_t exch_index = this->NUM_NOT_TOTALS;

  exch_comp.Get_totals().find(exch_comp.Get_formula())->second = exchange[0];
  exch_comp.Set_charge_balance(exchange[1]);
  exch_comp.Set_la(exchange[2]);
  exch_comp.Set_phase_proportion(exchange[3]);
  exch_comp.Set_formula_z(exchange[4]);

  for (auto &[name, value] : exch_comp.Get_totals()) {
    if (name != exch_comp.Get_formula()) {
      value = exchange[exch_index++];
    }
  }
}
