/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#pragma once

#include <phrqtype.h>
#include <span>

class WrapperBase {
public:
  virtual ~WrapperBase() = default;

  std::size_t size() const { return this->num_elements; };

  virtual void get(std::span<LDBLE> &data) const = 0;

  virtual void set(const std::span<LDBLE> &data) = 0;

protected:
  std::size_t num_elements = 0;
};