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

#include <cstdint>
/**
 * @brief A struct to hold the parameters for PhreeqcKnobs.
 *
 * @note See more information on the parameters in the [Phreeqc
 * manual](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-46.html)
 */
struct PhreeqcKnobsParams {
  //
  std::uint32_t iterations;
  double convergence_tolerance;
  double tolerance;
  double step_size;
  double pe_step_size;
  bool diagonal_scale;
};

class Phreeqc;

/**
 * @class PhreeqcKnobs
 * @brief A class to manage and manipulate knobs (parameters) for a Phreeqc
 * instance.
 *
 * This class provides functionalities to read and write knobs (parameters) for
 * a given Phreeqc instance. It also allows getting and setting the parameters
 * through accessor and mutator methods.
 *
 * @note The class supports default copy and move constructors and assignment
 * operators.
 *
 * @param pqc_instance A pointer to a Phreeqc instance.
 */
class PhreeqcKnobs {
public:
  /**
   * @brief Constructs a PhreeqcKnobs object from an exisiting phreeqc instance
   *
   * @param pqc_instance A pointer to an instance of the Phreeqc class.
   */
  PhreeqcKnobs(Phreeqc *pqc_instance);

  PhreeqcKnobs(const PhreeqcKnobs &) = default;
  PhreeqcKnobs(PhreeqcKnobs &&) = default;

  PhreeqcKnobs &operator=(const PhreeqcKnobs &) = default;
  PhreeqcKnobs &operator=(PhreeqcKnobs &&) = default;

  /**
   * @brief Reads and store the configuration knobs from the given Phreeqc
   * instance in the object.
   *
   * @param pqc_instance A pointer to the Phreeqc instance for which the knobs
   * are to be read.
   */
  void readKnobs(Phreeqc *pqc_instance);

  /**
   * @brief Writes the current state of the knobs to the given Phreeqc
   * instance.
   *
   * @param pqc_instance A pointer to the Phreeqc instance whose knobs are to be
   * written.
   */
  void writeKnobs(Phreeqc *pqc_instance) const;

  /**
   * @brief Retrieves the current PhreeqcKnobs parameters.
   *
   * @return PhreeqcKnobsParams The current parameters of the PhreeqcKnobs.
   */
  PhreeqcKnobsParams getParams() const { return _params; }

  /**
   * @brief Sets the parameters for PhreeqcKnobs.
   *
   * This function assigns the provided PhreeqcKnobsParams object to the
   * internal _params member variable, updating the configuration of the
   * PhreeqcKnobs instance.
   *
   * @param params The PhreeqcKnobsParams object containing the new parameters
   * to be set.
   */
  void setParams(const PhreeqcKnobsParams &params) { _params = params; }

private:
  PhreeqcKnobsParams _params;
};