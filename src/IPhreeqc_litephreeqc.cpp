/*
 * Copyright (c) 2024-2026 Max Luebke (University of Potsdam)
 *                       Marco De Lucia (GFZ German Research Centre for Geosciences)
 *
 * SPDX-License-Identifier: EUPL-1.2
 *
 * This file is part of litephreeqc, a C++ interface library on top of PHREEQC.
 */
#include "IPhreeqc.hpp"
#include "Phreeqc.h"
#include "ErrorReporter.hxx"

#include <string>
#include <exception>

int IPhreeqc::RunCell(int cell_no, double time_step) {
  static const char *sz_routine = "RunCell";
  try {
    this->ClearAccumulatedLines();
    this->ClearAccumulated = false;

    this->open_output_files(sz_routine);

    this->ErrorReporter->Clear();
    this->WarningReporter->Clear();

    this->PhreeqcPtr->input_error = 0;
    this->io_error_count = 0;

    this->PhreeqcPtr->run_cells(cell_no, time_step);
  } catch (const IPhreeqcStop &) {
    this->close_input_files();
  } catch (std::exception &e) {
    std::string errmsg("RunCell: ");
    errmsg += e.what();
    try {
      this->PhreeqcPtr->error_msg(errmsg.c_str(), STOP);
    } catch (const IPhreeqcStop &) {
    }
  } catch (...) {
    const char *errmsg = "RunCell: An unhandled exception occurred.\n";
    try {
      this->PhreeqcPtr->error_msg(errmsg, STOP);
    } catch (const IPhreeqcStop &) {
    }
  }
  this->update_errors();

  return this->PhreeqcPtr->get_input_errors();
}
