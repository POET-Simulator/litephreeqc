#!/usr/bin/env python3
"""
License Compliance Checker for litephreeqc and USGS PHREEQC/IPhreeqc code.

Ensures compliance with:
1. USGS User Rights Notice:
   - Modified USGS source files must include a prominent header notice stating:
     * Author(s) of modifications (Max Luebke, Marco De Lucia)
     * Date of modifications
     * Summary of modifications made
     * Reference to USGS User Rights Notice (phreeqc3-doc/NOTICE.TXT)
   - Modified files must be documented in Section 4 of LICENSE.md.

2. EUPL-1.2:
   - Native litephreeqc files must contain the EUPL-1.2 license header and SPDX tag.
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

# Paths considered native litephreeqc (EUPL-1.2)
NATIVE_LITEPHREEQC_PATTERNS = [
    re.compile(r"^litephreeqc/"),
    re.compile(r"^src/phreeqcpp/litephreeqc_funcs\.(cpp|h)"),
    re.compile(r"^src/IPhreeqc_litephreeqc\.cpp"),
]

# Paths considered third-party / external libraries (exempt from USGS / litephreeqc headers)
THIRD_PARTY_PATTERNS = [
    re.compile(r"^src/phreeqcpp/cvode"),
    re.compile(r"^src/phreeqcpp/PBasic"),
    re.compile(r"^gtest/"),
]

# Source file extensions to check
SOURCE_EXTENSIONS = {".c", ".cpp", ".cxx", ".h", ".hpp", ".hxx"}

USGS_NOTICE_MARKERS = [
    "Modified for litephreeqc",
    "USGS User Rights Notice",
    "Modifications:",
]

EUPL_MARKERS = [
    "SPDX-License-Identifier: EUPL-1.2",
    "Max Luebke",
]


def is_native_litephreeqc(rel_path: str) -> bool:
    return any(p.search(rel_path) for p in NATIVE_LITEPHREEQC_PATTERNS)


def is_third_party(rel_path: str) -> bool:
    return any(p.search(rel_path) for p in THIRD_PARTY_PATTERNS)


def is_usgs_source(rel_path: str) -> bool:
    if not (rel_path.startswith("src/") or rel_path.startswith("phreeqc3-doc/")):
        return False
    if is_native_litephreeqc(rel_path) or is_third_party(rel_path):
        return False
    _, ext = os.path.splitext(rel_path)
    return ext in SOURCE_EXTENSIONS


def get_git_changed_files(repo_root: Path, base_ref: str, include_working_tree: bool = True) -> list:
    files = set()
    try:
        # Changes committed on branch compared to base_ref
        res = subprocess.run(
            ["git", "diff", "--name-only", f"{base_ref}...HEAD"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=False,
        )
        if res.returncode == 0:
            for line in res.stdout.splitlines():
                if line.strip():
                    files.add(line.strip())
        else:
            # Fallback to direct diff
            res2 = subprocess.run(
                ["git", "diff", "--name-only", f"{base_ref}..HEAD"],
                cwd=repo_root,
                capture_output=True,
                text=True,
                check=True,
            )
            for line in res2.stdout.splitlines():
                if line.strip():
                    files.add(line.strip())

        if include_working_tree:
            wt_res = subprocess.run(
                ["git", "status", "--porcelain"],
                cwd=repo_root,
                capture_output=True,
                text=True,
                check=True,
            )
            for line in wt_res.stdout.splitlines():
                if len(line) > 3:
                    path = line[3:].strip()
                    # Handle renames: R old -> new
                    if " -> " in path:
                        path = path.split(" -> ")[1].strip()
                    files.add(path)
    except Exception as e:
        print(f"Error querying git: {e}", file=sys.stderr)
        sys.exit(2)

    return sorted(list(files))


def check_file_headers(repo_root: Path, rel_path: str) -> list:
    issues = []
    file_path = repo_root / rel_path
    if not file_path.exists() or not file_path.is_file():
        return issues

    try:
        content = file_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return [f"Could not read {rel_path}: {e}"]

    # We inspect the first 40 lines
    header = "\n".join(content.splitlines()[:40])

    if is_usgs_source(rel_path):
        for marker in USGS_NOTICE_MARKERS:
            if marker not in header:
                issues.append(f"Missing required USGS modification notice marker '{marker}' in header of {rel_path}")
    elif is_native_litephreeqc(rel_path):
        _, ext = os.path.splitext(rel_path)
        if ext in SOURCE_EXTENSIONS:
            for marker in EUPL_MARKERS:
                if marker not in header:
                    issues.append(f"Missing required EUPL-1.2 marker '{marker}' in header of {rel_path}")

    return issues


def check_license_md(repo_root: Path, modified_usgs_files: list) -> list:
    issues = []
    license_file = repo_root / "LICENSE.md"
    if not license_file.exists():
        return ["LICENSE.md file not found at repository root."]

    content = license_file.read_text(encoding="utf-8", errors="replace")

    sec4_match = re.search(r"## 4\. Summary of Modifications to USGS IPhreeqc Code(.*)", content, re.DOTALL)
    if not sec4_match:
        return ["Section '4. Summary of Modifications to USGS IPhreeqc Code' missing in LICENSE.md"]

    sec4_text = sec4_match.group(1)

    for usgs_file in modified_usgs_files:
        # Check if basename or relative path is mentioned in section 4
        basename = os.path.basename(usgs_file)
        if basename not in sec4_text and usgs_file not in sec4_text:
            issues.append(f"Modified USGS file '{usgs_file}' is not documented in LICENSE.md Section 4.")

    return issues


def main():
    parser = argparse.ArgumentParser(description="Check license compliance for litephreeqc repository.")
    parser.add_argument("--repo-root", default=".", help="Root path of the repository (default: current dir)")
    parser.add_argument("--base", default="origin/main", help="Git base reference for comparison (default: origin/main)")
    parser.add_argument("--no-working-tree", action="store_true", help="Exclude unstaged/staged working tree files")
    parser.add_argument("--all-files", action="store_true", help="Check all source files in repository instead of git diff")

    args = parser.parse_args()
    repo_root = Path(args.repo_root).resolve()

    print(f"Checking license compliance in: {repo_root}")
    print(f"Base ref: {args.base}")

    if args.all_files:
        files = []
        for p in repo_root.rglob("*"):
            if p.is_file() and not any(part.startswith(".") for part in p.relative_to(repo_root).parts):
                files.append(str(p.relative_to(repo_root)))
    else:
        files = get_git_changed_files(repo_root, args.base, include_working_tree=not args.no_working_tree)

    print(f"Files to inspect: {len(files)}")

    all_issues = []
    modified_usgs_files = []

    for f in files:
        if is_usgs_source(f):
            modified_usgs_files.append(f)
        issues = check_file_headers(repo_root, f)
        all_issues.extend(issues)

    # Check documentation in LICENSE.md
    if modified_usgs_files:
        doc_issues = check_license_md(repo_root, modified_usgs_files)
        all_issues.extend(doc_issues)

    if all_issues:
        print("\n[FAILED] License compliance issues detected:")
        for issue in all_issues:
            print(f"  - {issue}")
        sys.exit(1)
    else:
        print(f"\n[PASSED] All checked files ({len(files)}) and documentation meet license compliance requirements.")
        if modified_usgs_files:
            print(f"  Verified {len(modified_usgs_files)} modified USGS source file(s):")
            for mf in modified_usgs_files:
                print(f"    * {mf}")
        sys.exit(0)


if __name__ == "__main__":
    main()
