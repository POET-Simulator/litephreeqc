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