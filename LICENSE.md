# License Information

This repository contains two main components with differing open-source licenses:

1. **`litephreeqc`**: The C++ wrapper library, matrix abstractions, runner, knob management, tests, and utility functions located in the `litephreeqc/` directory and `src/phreeqcpp/litephreeqc_funcs.cpp`.
2. **IPhreeqc**: The underlying geochemical calculation engine developed by the U.S. Geological Survey (USGS) located in `src/`.

---

## 1. litephreeqc (EUPL-1.2)

The `litephreeqc` library and its added extensions are licensed under the **European Union Public Licence (EUPL) Version 1.2**.

```text
Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
                      Marco De Lucia (GFZ German Research Centre for Geosciences)

SPDX-License-Identifier: EUPL-1.2
```

The full text of the EUPL v1.2 licence is available in [LICENSES/EUPL-1.2.txt](LICENSES/EUPL-1.2.txt) or online at https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12.

Under Article 5 of the EUPL-1.2, derivative or combined works may also be distributed under any compatible licence listed in the EUPL Appendix (including GNU General Public License (GPL) v2/v3, AGPL v3, LGPL v2.1/v3, MPL v2, etc.).

---

## 2. IPhreeqc / PHREEQC (USGS User Rights Notice)

The IPhreeqc engine source code is developed by the U.S. Geological Survey (USGS) and provided under the **USGS User Rights Notice**.

The original notice is available in [phreeqc3-doc/NOTICE.TXT](phreeqc3-doc/NOTICE.TXT).

Any modifications made to original IPhreeqc source files to support `litephreeqc` are prominently noted in the respective files, including authors and nature of changes, as required by the USGS notice.

---

## 3. Third-Party Libraries

* **SUNDIALS / CVODE**: Developed by Lawrence Livermore National Laboratory and distributed under a BSD license (see notices in `src/phreeqcpp/cvode.*`).
* **Chipmunk BASIC**: Embedded basic interpreter routines (see notices in `src/phreeqcpp/PBasic.cpp`).
* **GoogleTest**: Used for unit testing, licensed under the BSD 3-Clause License.

---

## 4. Summary of Modifications to USGS IPhreeqc Code

In compliance with the USGS User Rights Notice, the following modifications were made to the original IPhreeqc source code to enable library integration with `litephreeqc`:

* **`src/IPhreeqc.hpp` & `src/IPhreeqc.cpp`**:
  - Added `GetPhreeqcPtr()` method to expose a direct pointer to the internal `Phreeqc` instance.
  - Added `struct SolutionMapping` and `getSolutionMapping()` to track module and solution index associations during simulation execution.
* **`src/phreeqcpp/Phreeqc.h` & `src/phreeqcpp/litephreeqc_funcs.cpp`**:
  - Added valence state identification helper methods (`find_all_valence_states`, `find_valence_states`) for geochemical species.
* **`src/phreeqcpp/Solution.h` & `src/phreeqcpp/Solution.cxx`**:
  - Added `Update()` overload accepting temperature (`tc`) and pressure (`patm`).
* **Authors of Modifications**:
  - Max Luebke (mluebke@uni-potsdam.de) - University of Potsdam
  - Marco De Lucia (delucia@gfz.de) - GFZ German Research Centre for Geosciences
* **Dates of Modifications**: 2024–2026.

### New Components in `litephreeqc/` (EUPL-1.2)
The files in `litephreeqc/` are newly authored components and not modifications of USGS files:
* **`PhreeqcEngine`**: High-level execution engine wrapping individual simulation cells.
* **`PhreeqcMatrix`**: In-memory matrix representation of multi-cell chemical states.
* **`PhreeqcRunner`**: Multi-cell scheduler and batch execution runner.
* **`PhreeqcKnobs`**: Convergence and numerical parameter controller.
* **`Wrapper/`**: Typed C++ wrappers for solutions, equilibrium phases, ion exchangers, surfaces, and kinetics.
