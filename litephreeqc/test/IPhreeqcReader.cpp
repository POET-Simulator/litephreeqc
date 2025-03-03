#include "IPhreeqcReader.hpp"

#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>

IPhreeqcReader::IPhreeqcReader(const std::string &database,
                               const std::string &script) {

  _m_pqc = std::make_unique<IPhreeqc>();

  _m_pqc->LoadDatabaseString(database.c_str());
  _m_pqc->RunString(script.c_str());

  if (_m_pqc->GetErrorStringLineCount() > 0) {
    std::cerr << ":: Error in Phreeqc script: " << _m_pqc->GetErrorString()
              << "\n";
    throw std::runtime_error("Phreeqc script error");
  }
}

void IPhreeqcReader::setOutputID(std::uint32_t cell_id) {

  _m_pqc->SetDumpStringOn(true);

  const std::string call_string =
      "DUMP\n -cells " + std::to_string(cell_id) + "\nEND\n";

  _m_pqc->RunString(call_string.c_str());

  if (_m_pqc->GetErrorStringLineCount() > 0) {
    std::cerr << ":: Error in Phreeqc script: " << _m_pqc->GetErrorString()
              << "\n";
    throw std::runtime_error("Possibly invalid cell ID: " +
                             std::to_string(cell_id));
  }

  _m_raw_output = _m_pqc->GetDumpString();

  _m_pqc->SetDumpStringOn(false);
}

double IPhreeqcReader::operator[](const std::string &name) const {
  if (_m_raw_output.empty()) {
    throw std::runtime_error("No output set. Call setOutputID first.");
  }

  std::size_t pos;

  auto it = _m_rename_map.find(name);

  if (it != _m_rename_map.end()) {
    pos = _m_raw_output.find(it->second);

    if (pos == std::string::npos) {
      throw std::runtime_error("Renamed name not found in output: " +
                               it->second);
    }

  } else if (name.ends_with("_eq") || name.ends_with("_si")) {
    // remove the suffix and check for the base name
    std::string base_name(name, 0, name.size() - 3); // remove "_eq"

    pos = _m_raw_output.find(base_name);

    if (pos == std::string::npos) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    if (name.ends_with("_eq")) {
      // if it ends with '_eq', we need to find '-moles'
      pos = _m_raw_output.find("-moles", pos);
    } else {
      // if it ends with '_si', we need to find '-si'
      pos = _m_raw_output.find("-si", pos);
    }
  } else {
    // if not, we first need to find '-totals'
    pos = _m_raw_output.find("-totals");
    if (pos == std::string::npos) {
      throw std::runtime_error("Name not found in output: " + name);
    }

    std::size_t next_field_pos = _m_raw_output.find(" -", pos);

    // Then, find the name in the output
    pos = _m_raw_output.find(name, pos);

    if (pos >= next_field_pos) {
      return 0;
    }
  }

  std::size_t end_pos = _m_raw_output.find('\n', pos);

  if (end_pos == std::string::npos) {
    throw std::runtime_error("End of line not found in output for name: " +
                             name);
  }

  std::regex fp_regex(
      R"( +([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?)");
  std::smatch match;

  const std::string line_to_parse = _m_raw_output.substr(pos, end_pos - pos);

  if (!std::regex_search(line_to_parse, match, fp_regex)) {
    throw std::runtime_error("No floating point number found in line: " +
                             line_to_parse);
  }

  return std::stod(match[0].str());
}

void IPhreeqcReader::run(double dt, std::vector<std::uint32_t> &&cell_ids) {
  std::string run_string =
      "RUN_CELLS\n -time_step " + std::to_string(dt) + "\n -cells ";

  for (const auto &cell_id : cell_ids) {
    run_string += std::to_string(cell_id) + " ";
  }

  run_string += "\nEND\n";

  _m_pqc->RunString(run_string.c_str());

  if (_m_pqc->GetErrorStringLineCount() > 0) {
    std::cerr << ":: Error in Phreeqc script: " << _m_pqc->GetErrorString()
              << "\n";
    throw std::runtime_error("Phreeqc run error");
  }

  _m_raw_output.clear();
}
