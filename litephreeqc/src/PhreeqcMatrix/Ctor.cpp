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

#include "IPhreeqc.hpp"
#include "PhreeqcKnobs.hpp"
#include "PhreeqcMatrix.hpp"

#include <Phreeqc.h>
#include <Solution.h>
#include <memory>
#include <regex>
#include <sstream>
#include <string>

static std::string getBlockByKeyword(const std::string &input_script,
                                     const std::string &keyword) {
  const std::regex keyword_regex("^" + keyword + "(?:\\s+.*)?$");
  const std::regex block_regex(R"(^[A-Z]+.*$)");

  bool block_found = false;
  std::size_t block_start = 0;
  std::size_t block_end = 0;
  std::size_t current_pos = 0;

  std::istringstream input_stream(input_script);

  for (std::string line; std::getline(input_stream, line);
       current_pos += line.length() + 1) {

    std::size_t first_char_pos = line.find_first_not_of(" \t\r");

    if (first_char_pos == std::string::npos) {
      continue;
    }

    std::string trimmed_line = line.substr(first_char_pos);

    if (std::regex_match(trimmed_line, keyword_regex)) {
      block_start = current_pos;
      block_found = true;
      continue;
    }

    if (!block_found) {
      continue;
    }

    if (!std::regex_search(trimmed_line, block_regex)) {
      continue;
    }

    block_end = current_pos - 1;
    break;
  }

  if (!block_found) {
    return "";
  }

  if (block_end == 0) {
    block_end = input_script.length();
  }

  return std::string(input_script, block_start, block_end - block_start + 1);
}

static std::string extractSelectedOutputBlock(const std::string &input_script) {
  std::string selected_output_block =
      getBlockByKeyword(input_script, "SELECTED_OUTPUT");
  std::string user_punch_block = getBlockByKeyword(input_script, "USER_PUNCH");

  if (selected_output_block.empty() && user_punch_block.empty()) {
    return "";
  }
  return selected_output_block + "\n" + user_punch_block + "\n";
}

PhreeqcMatrix::PhreeqcMatrix(const std::string &database,
                             const std::string &input_script, bool with_h0_o0,
                             bool with_redox)
    : _m_database(database), _m_with_h0_o0(with_h0_o0),
      _m_with_redox(with_redox) {
  this->_m_pqc = std::make_shared<IPhreeqc>();

  this->_m_pqc->LoadDatabaseString(database.c_str());

  this->_m_pqc->RunString(input_script.c_str());

  if (this->_m_pqc->GetErrorStringLineCount() > 0) {
    std::cerr << ":: Error in Phreeqc script: "
              << this->_m_pqc->GetErrorString() << "\n";
    throw std::runtime_error("Phreeqc script error");
  }

  this->_m_selected_output_block_string =
      extractSelectedOutputBlock(input_script);

  this->_m_knobs =
      std::make_shared<PhreeqcKnobs>(this->_m_pqc.get()->GetPhreeqcPtr());

  this->initialize();
}
