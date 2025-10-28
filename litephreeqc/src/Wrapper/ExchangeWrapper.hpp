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

#include "ExchComp.h"
#include "Exchange.h"
#include "WrapperBase.hpp"
#include "phrqtype.h"
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <span>

class ExchangeWrapper : public WrapperBase {
public:
  ExchangeWrapper(cxxExchange *exch, const std::vector<std::string> &exchanger);

  void get(std::span<LDBLE> &exchange) const;

  void set(const std::span<LDBLE> &exchange);

  static std::vector<std::string>
  names(cxxExchange *exchange, std::vector<std::string> &exchange_formulas);

private:
  cxxExchange *exchange;

  class ExchangeCompWrapper : public WrapperBase {
  public:
    ExchangeCompWrapper(cxxExchComp &comp);

    void get(std::span<LDBLE> &exchange) const;

    void set(const std::span<LDBLE> &exchange);

  private:
    cxxExchComp &exch_comp;

    static constexpr std::size_t NUM_NOT_TOTALS = 5;
  };

  std::vector<std::unique_ptr<ExchangeCompWrapper>> exchange_comps;
};