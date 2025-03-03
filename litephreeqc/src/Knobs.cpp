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

#include "PhreeqcKnobs.hpp"

#include <Phreeqc.h>

PhreeqcKnobs::PhreeqcKnobs(Phreeqc *pqc_instance) {
  this->readKnobs(pqc_instance);
}

void PhreeqcKnobs::readKnobs(Phreeqc *pqc_instance) {
  this->_params.iterations = pqc_instance->itmax;
  this->_params.convergence_tolerance = pqc_instance->convergence_tolerance;
  this->_params.tolerance = pqc_instance->ineq_tol;
  this->_params.step_size = pqc_instance->step_size;
  this->_params.pe_step_size = pqc_instance->pe_step_size;
  this->_params.diagonal_scale =
      static_cast<bool>(pqc_instance->diagonal_scale);
}

void PhreeqcKnobs::writeKnobs(Phreeqc *pqc_instance) const {
  pqc_instance->itmax = this->_params.iterations;
  pqc_instance->convergence_tolerance = this->_params.convergence_tolerance;
  pqc_instance->ineq_tol = this->_params.tolerance;
  pqc_instance->step_size = this->_params.step_size;
  pqc_instance->pe_step_size = this->_params.pe_step_size;
  pqc_instance->diagonal_scale = this->_params.diagonal_scale;
}