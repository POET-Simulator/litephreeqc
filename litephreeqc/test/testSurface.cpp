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

// SURFACE columns of the PhreeqcMatrix for the surface models -no_edl, DDL,
// Donnan with fixed thickness and Donnan with the thickness in Debye lengths:
// the columns present, their roles, and that the roles hold (a hint does not
// change the result, an immobile amount counts only through the total of the
// cell, a parameter or an amount alone changes it).

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <testInput.hpp>

#include <gtest/gtest.h>

#include "PhreeqcEngine.hpp"
#include "PhreeqcMatrix.hpp"
#include "utils.hpp"

namespace {

const std::string database = readFile(base_test::phreeqc_database);

using Role = PhreeqcMatrix::ColumnRole;

const std::map<std::string, std::string> models = {
    {"no_edl", "-no_edl"},
    {"ddl", ""},
    {"donnan", "-donnan"},
    {"donnan_debye", "-donnan debye_lengths 1"}};

std::string script(const std::string &option) {
  return "SOLUTION 1\n units mol/kgw\n pH 7\n Na 0.01\n Cl 0.01 charge\n"
         " Zn 1e-5\n"
         "SURFACE 1\n -equilibrate 1\n " +
         option +
         "\n Hfo_wOH 2e-4 600 1\n Hfo_sOH 5e-6\n"
         "RUN_CELLS\n -cells 1\nEND\n";
}

struct Cell {
  PhreeqcMatrix mat;
  std::vector<std::string> names;
  std::vector<Role> roles;
  std::vector<double> values;

  explicit Cell(const std::string &model)
      : mat(database, script(models.at(model))) {
    const auto exported = mat.get();
    names.assign(exported.names.begin() + 1, exported.names.end());
    roles.assign(exported.roles.begin() + 1, exported.roles.end());
    values.assign(exported.values.begin() + 1, exported.values.end());
  }

  bool has(const std::string &name) const {
    return std::find(names.begin(), names.end(), name) != names.end();
  }

  std::size_t at(const std::string &name) const {
    const auto it = std::find(names.begin(), names.end(), name);
    if (it == names.end()) {
      throw std::out_of_range("no column " + name);
    }
    return static_cast<std::size_t>(it - names.begin());
  }

  Role role(const std::string &name) const { return roles[at(name)]; }
};

// pH and total Zn of solution and sites after equilibration
struct Result {
  double pH, zn_aq, zn_w;
};

class Runner {
public:
  explicit Runner(const Cell &cell) : cell(cell), engine(cell.mat, 1) {
    // a non-equilibrium input: Zn added to the solution
    x0 = cell.values;
    x0[cell.at("Zn")] *= 2;
    ref = run(x0);
  }

  Result run(std::vector<double> x) {
    engine.runCell(x, 0);
    return {x[cell.at("pH")], x[cell.at("Zn")], x[cell.at("Hfo_wOH_Zn")]};
  }

  // largest change against the reference: pH (absolute), Zn (relative)
  double change(const std::vector<double> &x) {
    const Result r = run(x);
    return std::max({std::abs(r.pH - ref.pH),
                     std::abs(r.zn_aq / ref.zn_aq - 1),
                     std::abs(r.zn_w / ref.zn_w - 1)});
  }

  std::vector<double> x0;

private:
  const Cell &cell;
  PhreeqcEngine engine;
  Result ref;
};

constexpr double NO_CHANGE = 1e-6;
// a starting guess moves the result within the convergence of the solver,
// which is looser with a diffuse layer
constexpr double HINT_CHANGE = 1e-5;
constexpr double CHANGE = 1e-4;

} // namespace

POET_TEST(SurfaceColumns) {
  for (const auto &[model, _] : models) {
    SCOPED_TRACE(model);
    const Cell cell(model);
    const bool edl = model != "no_edl";
    const bool donnan = model.starts_with("donnan");

    for (const auto &name : cell.names) {
      EXPECT_FALSE(name.ends_with("_moles")) << name;
      // pH and tc are no element names
      EXPECT_FALSE(name.ends_with("_tot_P") || name.ends_with("_tot_Tc"))
          << name;
    }
    EXPECT_EQ(cell.has("Hfo_wOH_cb"), !edl);
    EXPECT_EQ(cell.has("Hfo_cb"), edl);
    EXPECT_EQ(cell.has("Hfo_area"), edl);
    EXPECT_EQ(cell.has("Hfo_tot_Na"), donnan);
    EXPECT_EQ(cell.has("Hfo_mw"), model == "donnan");

    EXPECT_EQ(cell.role("Hfo_wOH_la"), Role::HINT);
    EXPECT_EQ(cell.role("Hfo_wOH_Zn"), Role::IMMOBILE);
    EXPECT_EQ(cell.role("Hfo_wOH_Hfo_w"), Role::IMMOBILE);
    EXPECT_EQ(cell.role("pH"), Role::HINT);
    EXPECT_EQ(cell.role("tc"), Role::PARAMETER);
    EXPECT_EQ(cell.role("Zn"), Role::TRANSPORTED);
    if (edl) {
      EXPECT_EQ(cell.role("Hfo_la"), Role::HINT);
      EXPECT_EQ(cell.role("Hfo_cb"), Role::IMMOBILE);
      EXPECT_EQ(cell.role("Hfo_grams"), Role::PARAMETER);
    } else {
      EXPECT_EQ(cell.role("Hfo_wOH_cb"), Role::IMMOBILE);
    }
    if (model == "donnan") {
      EXPECT_EQ(cell.role("Hfo_mw"), Role::PARAMETER);
    }
    if (donnan) {
      EXPECT_EQ(cell.role("Hfo_tot_Zn"), Role::IMMOBILE);
    }
  }
}

POET_TEST(SurfaceRolesHold) {
  for (const auto &[model, _] : models) {
    SCOPED_TRACE(model);
    const Cell cell(model);
    Runner run(cell);
    const bool edl = model != "no_edl";

    EXPECT_LT(run.change(run.x0), NO_CHANGE) << "rerun";

    // hints are starting guesses only
    for (std::size_t i = 0; i < cell.names.size(); i++) {
      if (cell.roles[i] != Role::HINT) {
        continue;
      }
      auto x = run.x0;
      x[i] += 1;
      EXPECT_LT(run.change(x), HINT_CHANGE) << cell.names[i];
    }

    // the charge balance that is read counts through the cell total only
    const std::string cb = edl ? "Hfo_cb" : "Hfo_wOH_cb";
    auto x = run.x0;
    x[cell.at(cb)] += 1e-7;
    EXPECT_GT(run.change(x), CHANGE) << cb << " alone";
    x[cell.at("Charge")] -= 1e-7;
    EXPECT_LT(run.change(x), NO_CHANGE) << cb << " moved to Charge";

    // sorbed amounts count through the cell total only
    const double dzn = 0.1 * run.x0[cell.at("Hfo_wOH_Zn")];
    x = run.x0;
    x[cell.at("Hfo_wOH_Zn")] += dzn;
    EXPECT_GT(run.change(x), CHANGE) << "Hfo_wOH_Zn alone";
    x[cell.at("Zn")] -= dzn;
    EXPECT_LT(run.change(x), NO_CHANGE) << "Zn moved onto Hfo_wOH";

    if (edl) {
      x = run.x0;
      x[cell.at("Hfo_grams")] *= 2;
      EXPECT_GT(run.change(x), CHANGE) << "Hfo_grams";
    }

    if (model.starts_with("donnan")) {
      // diffuse layer amounts count through the cell total only
      const double dna = 1e-4;
      x = run.x0;
      x[cell.at("Hfo_tot_Na")] += dna;
      EXPECT_GT(run.change(x), CHANGE) << "Hfo_tot_Na alone";
      x[cell.at("Na")] -= dna;
      EXPECT_LT(run.change(x), NO_CHANGE) << "Na moved into the layer";
    }

    if (model == "donnan") {
      // fixed thickness: the water of the layer is an input
      x = run.x0;
      x[cell.at("Hfo_mw")] *= 2;
      EXPECT_GT(run.change(x), CHANGE) << "Hfo_mw";
    }
  }
}

POET_TEST(EngineRejectsWrongRowSize) {
  const Cell cell("ddl");
  PhreeqcEngine engine(cell.mat, 1);

  std::vector<double> x = cell.values;
  EXPECT_NO_THROW(engine.runCell(x, 0));
  x.pop_back();
  EXPECT_THROW(engine.runCell(x, 0), std::invalid_argument);
}
