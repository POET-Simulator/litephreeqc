---
name: phreeqc-license-compliance
description: >-
  Use this skill to ensure license compliance in litephreeqc and PHREEQC/IPhreeqc codebases.
  It guides the verification and application of EUPL-1.2 headers for native components,
  USGS User Rights Notice modification notices for altered USGS files, and synchronization with LICENSE.md.
---

# PHREEQC & litephreeqc License Compliance Runbook

This skill provides procedures and tooling to maintain compliance with:
1. **USGS User Rights Notice** (see `phreeqc3-doc/NOTICE.TXT`): All modified USGS files must prominently state the modifications, authors, and dates, and must be summarized in `LICENSE.md`.
2. **European Union Public Licence (EUPL-1.2)**: All native `litephreeqc` components must carry the standard EUPL-1.2 copyright header and SPDX identifier.

---

## Quick Verification

To verify license compliance against the base branch (`origin/main`):

```bash
python3 .agents/skills/phreeqc-license-compliance/scripts/check_license_compliance.py --base origin/main
```

Parameters:
- `--base <ref>`: Git reference to compare against (default: `origin/main`).
- `--no-working-tree`: Check only committed changes, ignore unstaged/staged working tree files.
- `--all-files`: Audit every source file in the repository.

---

## Standard License Headers

### 1. USGS Modified File Header
Place this at the very top of any file under `src/` (or `phreeqc3-doc/`) modified from original USGS distribution:

```cpp
/*
 * Modified for litephreeqc by Max Luebke (University of Potsdam)
 * and Marco De Lucia (GFZ German Research Centre for Geosciences), 2024-2026.
 * Modifications: <Brief summary of modifications>.
 * Original code is subject to the USGS User Rights Notice (see phreeqc3-doc/NOTICE.TXT).
 */
```

### 2. Native litephreeqc Header (EUPL-1.2)
Place this at the top of all new `litephreeqc` source and header files (`litephreeqc/`, `src/phreeqcpp/litephreeqc_funcs.cpp`, `src/IPhreeqc_litephreeqc.cpp`):

```cpp
/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
```

---

## Step-by-Step Workflow

1. **Identify Changed Files**:
   Run `git diff --name-only origin/main...HEAD` (plus `git status --porcelain`).

2. **Categorize Files**:
   - Files in `litephreeqc/` -> Require EUPL-1.2 header.
   - Files in `src/` originally from USGS -> Require USGS modification header.
   - Build configs / docs / tests -> Standard project conventions.

3. **Apply / Update Headers**:
   - If an existing USGS file already has a modification notice, update the `Modifications:` line to include the new changes.
   - If it does not have a notice, prepend the standard USGS modification header.

4. **Update `LICENSE.md`**:
   Add or update the file and a bulleted list of changes in Section `4. Summary of Modifications to USGS IPhreeqc Code`.

5. **Run Verification Script**:
   ```bash
   python3 .agents/skills/phreeqc-license-compliance/scripts/check_license_compliance.py --base origin/main
   ```

6. **Compile & Test**:
   Ensure modifications didn't break builds or tests:
   ```bash
   cmake --build build -j$(nproc)
   ctest --test-dir build --output-on-failure
   ```

7. **Create Git Commit**:
   Commit the license adjustments using Conventional Commits:
   ```bash
   git commit -m "docs(license): update file headers and modifications summary"
   ```
