#include "IPhreeqc.hpp"
#include <cstdint>

/**
 * @class PhreeqcSelectedOutputParser
 * @brief Parses and manages the SELECTED_OUTPUT block from a PHREEQC input
 * script.
 *
 * This class is responsible for extracting and handling the SELECTED_OUTPUT
 * block from a PHREEQC input script. It provides methods to check if a
 * SELECTED_OUTPUT block is present, retrieve its string representation, access
 * the parsed header, and obtain output values for a specific cell.
 *
 * @note The class requires an instance of IPhreeqc to function.
 */
class PhreeqcSelectedOutputParser {
public:
  /**
   * @brief Constructs a PhreeqcSelectedOutputParser object.
   *
   * Initializes the parser with a given IPhreeqc instance and an input script.
   * The parser will use the provided PHREEQC instance to process the specified
   * input script and extract selected output data.
   *
   * @param pqc_instance Pointer to an IPhreeqc instance used for PHREEQC
   * calculations.
   * @param input_script The PHREEQC input script to be processed.
   */
  PhreeqcSelectedOutputParser(IPhreeqc *pqc_instance,
                              const std::string &input_script);

  /**
   * @brief Constructs a PhreeqcSelectedOutputParser by associating it with an
   * existing IPhreeqc instance and initializing it using another
   * PhreeqcSelectedOutputParser.
   *
   * This constructor allows creating a new parser that shares the IPhreeqc
   * context and copies relevant state or configuration from an existing parser
   * instance.
   *
   * @param pqc_instance Pointer to an IPhreeqc instance to be used by the
   * parser.
   * @param input_parser Reference to an existing PhreeqcSelectedOutputParser to
   * initialize from.
   */
  PhreeqcSelectedOutputParser(IPhreeqc *pqc_instance,
                              const PhreeqcSelectedOutputParser &input_parser);

  /**
   * @brief Checks if selected output is available.
   *
   * @return true if selected output is present, false otherwise.
   */
  bool hasSelectedOutput() const { return _m_has_selected_output; }

  /**
   * @brief Retrieves the string representation of the selected output block.
   *
   * @return A constant reference to the string containing the selected output
   * block.
   */
  const std::string &getSelectedOutputBlockString() const {
    return _m_selected_output_block_string;
  }

  /**
   * @brief Retrieves the header row of the selected output.
   *
   * @return A vector of strings containing the column headings.
   */
  std::vector<std::string> getHeader() const { return _m_headings; }

  /**
   * @brief Retrieves the selected output values for a specified cell.
   *
   * This function returns a vector containing the output values associated with
   * the given cell ID.
   *
   * @param cell_id The identifier of the Phreeqc cell for which to retrieve the
   * output values.
   * @return std::vector<double> A vector of output values corresponding to the
   * specified cell.
   */
  std::vector<double> getValues(std::uint32_t cell_id) const;

private:
  void parseSelectedOutputBlock(const std::string &input_script);
  void parseHeader();

  bool _m_has_selected_output{false}; // true if selected output was defined
  std::string _m_selected_output_block_string;

  IPhreeqc *_m_pqc_instance;
  std::vector<std::string> _m_headings;
};
