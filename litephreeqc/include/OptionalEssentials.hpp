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

// This header is deprecated. OptionalEssentials is now part of PhreeqcMatrix.
// Include PhreeqcMatrix.hpp directly and use PhreeqcMatrix::OptionalEssentials.

#include "PhreeqcMatrix.hpp"

// Backward compatibility alias (deprecated)
namespace litepqc {
using OptionalEssentials = PhreeqcMatrix::OptionalEssentials;

inline constexpr bool hasFlag(OptionalEssentials set, OptionalEssentials flag) {
  return PhreeqcMatrix::hasFlag(set, flag);
}
} // namespace litepqc
