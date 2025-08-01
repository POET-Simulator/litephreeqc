#include "PhreeqcSelectedOutputParser.hpp"

#include <limits>
#include <regex>
#include <sstream>

PhreeqcSelectedOutputParser::PhreeqcSelectedOutputParser(
    IPhreeqc *pqc_instance, const std::string &input_script)
    : _m_pqc_instance(pqc_instance) {
  if (!_m_pqc_instance->GetSelectedOutputStringOn()) {
    throw std::runtime_error("Selected output string is not enabled.");
  }

  this->parseSelectedOutputBlock(input_script);

  if (this->_m_has_selected_output) {
    parseHeader();
    if (this->_m_headings.size() != this->getValues(1).size()) {
      throw std::runtime_error(
          "Number of headings does not match number of values in selected "
          "output.");
    }
  }
}

PhreeqcSelectedOutputParser::PhreeqcSelectedOutputParser(
    IPhreeqc *pqc_instance, const PhreeqcSelectedOutputParser &input_parser)
    : _m_pqc_instance(pqc_instance),
      _m_has_selected_output(input_parser._m_has_selected_output),
      _m_selected_output_block_string(
          input_parser._m_selected_output_block_string),
      _m_headings(input_parser._m_headings) {
  if (_m_has_selected_output) {
    this->_m_pqc_instance->SetSelectedOutputStringOn(true);
    this->_m_pqc_instance->RunString(
        this->_m_selected_output_block_string.c_str());

    if (this->_m_pqc_instance->GetErrorStringLineCount() > 0) {
      throw std::runtime_error("Error running selected output block!");
    }
  }
}

std::string getBlockByKeyword(const std::string &input_script,
                              const std::string &&keyword) {
  const std::regex keyword_regex("^" + keyword + "*$");
  const std::regex block_regex(R"(^[A-Z]+.*$)");

  bool block_found = false;
  std::size_t block_start = 0;
  std::size_t block_end = 0;
  std::size_t current_pos = 0;

  std::istringstream input_stream(input_script);

  for (std::string line; std::getline(input_stream, line);
       current_pos += line.length() + 1) {

    // remove trailing whitespaces
    std::size_t first_char_pos = line.find_first_not_of(" \t\r");

    if (first_char_pos == std::string::npos) {
      continue; // Skip empty lines
    }

    if (std::regex_search(line, keyword_regex)) {
      block_start = current_pos;
      block_found = true;
      continue;
    }

    if (!block_found) {
      continue;
    }

    if (!std::regex_search(
            line.substr(first_char_pos, line.length() - first_char_pos),
            block_regex)) {
      continue; // Skip lines that do not start with a capitalized keyword
    }

    block_end = current_pos - 1; // End of KEYWORD block
    break;
  }

  if (!block_found) {
    return {""};
  }

  return std::string(input_script, block_start, block_end - block_start + 1);
}

void PhreeqcSelectedOutputParser::parseSelectedOutputBlock(
    const std::string &input_script) {
  std::string selected_output_block =
      getBlockByKeyword(input_script, "SELECTED_OUTPUT");

  std::string user_punch_block = getBlockByKeyword(input_script, "USER_PUNCH");

  if (selected_output_block.empty() && user_punch_block.empty()) {
    return;
  }

  if (selected_output_block.empty() && !user_punch_block.empty()) {
    throw std::runtime_error(
        "USER_PUNCH block found without a SELECTED_OUTPUT block.");
  }

  this->_m_has_selected_output = true;
  this->_m_selected_output_block_string =
      selected_output_block + user_punch_block;
}

void PhreeqcSelectedOutputParser::parseHeader() {
  const std::string selected_output_string =
      _m_pqc_instance->GetSelectedOutputString();

  std::istringstream stream(selected_output_string);
  std::string header_line;

  if (!std::getline(stream, header_line)) {
    throw std::runtime_error("No headings found in selected output string.");
  }

  std::istringstream header_stream(header_line);

  for (std::string heading; std::getline(header_stream, heading, '\t');) {
    std::size_t first_char_pos = heading.find_first_not_of(" ");
    _m_headings.push_back(
        heading.substr(first_char_pos, heading.length() - first_char_pos) +
        "_SO");
  }
}

std::vector<double>
PhreeqcSelectedOutputParser::getValues(std::uint32_t cell_id) const {
  if (!this->hasSelectedOutput()) {
    return {}; // No selected output defined
  }
  this->_m_pqc_instance->SetCurrentSelectedOutputUserNumber(cell_id);

  const std::string selected_output_string =
      _m_pqc_instance->GetSelectedOutputString();

  if (selected_output_string.empty()) {
    throw std::runtime_error("Selected output string is empty.");
  }

  std::istringstream stream(selected_output_string);
  std::string line, tmp_line;
  std::vector<double> values;

  while (std::getline(stream, tmp_line)) {
    line = tmp_line;
  }

  std::istringstream line_stream(line);

  for (std::string value; std::getline(line_stream, value, '\t');) {
    try {
      values.push_back(std::stod(value));
    } catch (const std::invalid_argument &) {
      values.push_back(
          std::numeric_limits<double>::quiet_NaN()); // If conversion fails,
                                                     // store NaN
    }
  }

  return values;
}
