/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
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