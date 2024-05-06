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