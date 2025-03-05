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

#include <gtest/gtest.h>
#include <testInput.hpp>

#include "IPhreeqc.hpp"
#include "PhreeqcKnobs.hpp"
#include "utils.hpp"

const std::string barite_db = readFile(barite_test::database);
const std::string barite_script = readFile(barite_test::script);

const std::string knob_input = R"(
KNOBS 
    -iterations 120
    -convergence_tolerance 1e-12
    -tolerance 1e-16
    -step_size 200
    -pe_step_size 20
    -diagonal_scale true
END
)";

POET_TEST(PhreeqcKnobsDefaultParams) {
  IPhreeqc pqc;

  pqc.LoadDatabaseString(barite_db.c_str());
  pqc.RunString(barite_script.c_str());

  PhreeqcKnobs knobs(pqc.GetPhreeqcPtr());

  const PhreeqcKnobsParams params = knobs.getParams();

  EXPECT_EQ(params.iterations, 100);
  EXPECT_DOUBLE_EQ(params.convergence_tolerance, 1e-8);
  EXPECT_DOUBLE_EQ(params.tolerance, 1e-15);
  EXPECT_DOUBLE_EQ(params.step_size, 100);
  EXPECT_DOUBLE_EQ(params.pe_step_size, 10);
  EXPECT_FALSE(params.diagonal_scale);
}

inline void compare_params(const PhreeqcKnobsParams &params) {
  EXPECT_EQ(params.iterations, 120);
  EXPECT_DOUBLE_EQ(params.convergence_tolerance, 1e-12);
  EXPECT_DOUBLE_EQ(params.tolerance, 1e-16);
  EXPECT_DOUBLE_EQ(params.step_size, 200);
  EXPECT_DOUBLE_EQ(params.pe_step_size, 20);
  EXPECT_TRUE(params.diagonal_scale);
}

POET_TEST(PhreeqcKnobsSetFromScript) {
  IPhreeqc pqc;

  pqc.LoadDatabaseString(barite_db.c_str());
  pqc.RunString(knob_input.c_str());
  pqc.RunString(barite_script.c_str());

  PhreeqcKnobs knobs(pqc.GetPhreeqcPtr());

  const PhreeqcKnobsParams params = knobs.getParams();

  compare_params(params);
}

POET_TEST(PhreeqcKnobsSetFromClass) {
  IPhreeqc pqc;

  pqc.LoadDatabaseString(barite_db.c_str());
  pqc.RunString(knob_input.c_str());
  pqc.RunString(barite_script.c_str());

  PhreeqcKnobs knobs(pqc.GetPhreeqcPtr());

  IPhreeqc new_instance;

  new_instance.LoadDatabaseString(barite_db.c_str());
  new_instance.RunString(barite_script.c_str());

  knobs.writeKnobs(new_instance.GetPhreeqcPtr());

  PhreeqcKnobs new_knobs(new_instance.GetPhreeqcPtr());

  const PhreeqcKnobsParams params = new_knobs.getParams();

  compare_params(params);
}