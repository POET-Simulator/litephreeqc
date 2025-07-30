#include "IPhreeqc.hpp"
#include "PhreeqcKnobs.hpp"
#include "PhreeqcMatrix.hpp"

#include <Phreeqc.h>
#include <Solution.h>
#include <memory>
#include <string>

PhreeqcMatrix::PhreeqcMatrix(const std::string &database,
                             const std::string &input_script, bool with_h0_o0,
                             bool with_redox)
    : _m_database(database), _m_with_h0_o0(with_h0_o0),
      _m_with_redox(with_redox) {
  this->_m_pqc = std::make_shared<IPhreeqc>();

  this->_m_pqc->LoadDatabaseString(database.c_str());

  this->_m_pqc->SetSelectedOutputStringOn(true);

  this->_m_pqc->RunString(input_script.c_str());

  if (this->_m_pqc->GetErrorStringLineCount() > 0) {
    std::cerr << ":: Error in Phreeqc script: "
              << this->_m_pqc->GetErrorString() << "\n";
    throw std::runtime_error("Phreeqc script error");
  }

  this->_m_selected_output_parser =
      std::make_shared<PhreeqcSelectedOutputParser>(this->_m_pqc.get(),
                                                    input_script);

  this->_m_knobs =
      std::make_shared<PhreeqcKnobs>(this->_m_pqc.get()->GetPhreeqcPtr());

  this->initialize();
}
