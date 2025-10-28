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

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "PhreeqcKnobs.hpp"

class IPhreeqc;

/**
 * @brief Class for storing information from Phreeqc
 *
 * PhreeqcMatrix is used as a container for storing **essential** information
 * from Phreeqc in a C++ data structure. The usage of Phreeqc's interpreter is
 * minimized. Thus, values are written directly from/to Phreeqc's internal data
 * structures, minimizing the overhead of parsing and eliminate floating point
 * errors due to conversion.
 *
 * The class is also used to initialize the PhreeqcEngine and PhreeqcRunner
 * class.
 */
class PhreeqcMatrix {
public:
  /**
   * @brief Construct a new Phreeqc Matrix object
   *
   * Default constructor. Does nothing. Used only for assignment operator.
   */
  PhreeqcMatrix() = default;

  /**
   * @brief Destroy the Phreeqc Matrix object
   *
   * There is no special cleanup needed for PhreeqcMatrix, thus the default (and
   * implicit) destructor is used.
   */
  ~PhreeqcMatrix() = default;

  /**
   * @brief Constructs a PhreeqcMatrix object with the specified database and
   * input script.
   *
   * This constructor initializes a PhreeqcMatrix object using the provided
   * database and input script. It sets the default values to exclude H(0) and
   * O(0) in the output and to include redox states.
   *
   * @param database The path to the database file.
   * @param input_script The input script to be used.
   */
  PhreeqcMatrix(const std::string &database, const std::string &input_script)
      : PhreeqcMatrix(database, input_script, false, true) {}

  /**
   * @brief Construct a new Phreeqc Matrix object
   *
   * Construct a new Phreeqc Matrix object by reading the database and input
   * script already present as a string.
   *
   * @param database Phreeqc database as a string.
   * @param input_script Phreeqc input script as a string.
   * @param with_h0_o0 Whether to include H(0) and O(0) in the output or not.
   * @param with_redox Whether to include redox states in the output or not.
   */
  PhreeqcMatrix(const std::string &database, const std::string &input_script,
                bool with_h0_o0, bool with_redox);

  /**
   * @brief Construct a new Phreeqc Matrix object
   *
   * Copy constructor. The interal used Phreeqc instance is reused by the new
   * object!
   * @param other
   */
  PhreeqcMatrix(const PhreeqcMatrix &other) = default;

  /**
   * @brief Construct a new Phreeqc Matrix object
   *
   * Copy constructor. The interal used Phreeqc instance is reused by the new
   * object!
   * @param other
   */
  PhreeqcMatrix(PhreeqcMatrix &&other) = default;

  /**
   * @brief Assignment operator
   *
   * The interal used Phreeqc instance is reused by the new object!
   * @param other
   * @return PhreeqcMatrix&
   */
  PhreeqcMatrix &operator=(const PhreeqcMatrix &other) = default;

  /**
   * @brief Assignment operator
   *
   * The interal used Phreeqc instance is reused by the new object!
   * @param other
   * @return PhreeqcMatrix&
   */
  PhreeqcMatrix &operator=(PhreeqcMatrix &&other) = default;

  /**
   * @brief Access the value of a given cell by name.
   *
   * @param cell_id ID of the cell (user id from Phreeqc script) to get the
   * value from.
   * @param name Name of the component to get the value from.
   * @return double The stored value.
   * @throw std::runtime_error if the component is not found.
   */
  double operator()(int cell_id, const std::string &name) const;

  /**
   * @brief Subset the PhreeqcMatrix to given cell IDs.
   *
   * With a given set of cell IDs, a new PhreeqcMatrix is created containing
   * only the given cell IDs. All entries, which values refer to NaN, are
   * removed.
   *
   * @param indices Cell IDs to subset the PhreeqcMatrix to.
   * @return PhreeqcMatrix A new PhreeqcMatrix containing only the given cell
   * IDs.
   */
  PhreeqcMatrix subset(const std::vector<int> &indices) const;

  /**
   * @brief Erase the given cell IDs from the PhreeqcMatrix.
   *
   * With a given set of cell IDs, the PhreeqcMatrix is modified to contain only
   * the cell IDs not in the given set. All entries, which values refer to NaN,
   * are removed.
   *
   * @param indices Cell IDs to erase from the PhreeqcMatrix.
   * @return PhreeqcMatrix A new PhreeqcMatrix containing only the cell IDs not
   * in the given set.
   */
  PhreeqcMatrix erase(const std::vector<int> &indices) const;

  /**
   * @brief Type of vector export
   *
   */
  enum class VectorExportType { COLUMN_MAJOR, ROW_MAJOR };

  /**
   * @brief Struct holding a format of the exported data
   *
   */
  struct STLExport {
    std::vector<std::string> names;
    std::vector<double> values;
  };

  /**
   * @brief Export the internal data to consecutive vectors.
   *
   * @param type Type of the order of the exported data
   * @param include_id Whether to include a column with the cell IDs or not
   * @return STLExport Exported data
   */
  STLExport get(VectorExportType type = VectorExportType::ROW_MAJOR,
                bool include_id = true) const;

  enum class PhreeqcComponent {
    SOLUTION = 0,
    EXCHANGE,
    KINETIC,
    EQUILIBRIUM,
    SURFACE_COMPS
  };

  struct element {
    std::string name;
    PhreeqcComponent type;
    double value;
  };

  struct base_names {
    enum class Components {
      EXCHANGER = static_cast<int>(PhreeqcComponent::EXCHANGE),
      KINETICS,
      EQUILIBRIUM,
      SURACE_COMP,
      SURFACE_CHARGE
    } type;

    std::string name;
  };

  /**
   * @brief Get all found solution names.
   *
   * @return std::vector<std::string> Vector containing all solution names.
   */
  std::vector<std::string> getSolutionNames() const;

  /**
   * @brief Get solution total names of all found solutions (excluding H, O,
   * Charge, H(0), O(0))
   *
   * @return std::vector<std::string> Names of all found solutions (excluding H,
   * O, Charge, H(0), O(0))
   */
  std::vector<std::string> getSolutionPrimaries() const;

  /**
   * @brief Get the exchange names for a given cell
   *
   * @param cell_id ID of the cell to get the exchange names for
   * @return std::vector<std::string> Whole vector of exchange names for the
   * cell. Empty if no exchange is defined.
   */
  std::vector<std::string> getExchanger(int cell_id) const;

  /**
   * @brief Get the kinetics names for a given cell
   *
   * @param cell_id ID of the cell to get the kinetics names for
   * @return std::vector<std::string> Whole vector of kinetics names for the
   * cell. Empty if no kinetics are defined.
   */
  std::vector<std::string> getKineticsNames(int cell_id) const;

  /**
   * @brief Get the equilibrium names for a given cell
   *
   * @param cell_id ID of the cell to get the equilibrium names for
   * @return std::vector<std::string> Whole vector of equilibrium names for the
   * cell. Empty if no equilibrium is defined.
   */
  std::vector<std::string> getEquilibriumNames(int cell_id) const;

  /**
   * @brief Get the surface component names for a given cell
   *
   * @param cell_id ID of the cell to get the surface component names for
   * @return std::vector<std::string> Whole vector of surface component names
   * for the cell. Empty if no surface is defined.
   */
  std::vector<std::string> getSurfaceCompNames(int cell_id) const;

  /**
   * @brief Get the surface charge names for a given cell
   *
   * @param cell_id ID of the cell to get the surface charge names for
   * @return std::vector<std::string> Whole vector of surface charge names for
   * the cell. Empty if no surface is defined.
   */
  std::vector<std::string> getSurfaceChargeNames(int cell_id) const;

  /**
   * @brief Get all cell IDs stored in the PhreeqcMatrix.
   *
   * @return std::vector<int> IDs of all cells stored in the PhreeqcMatrix.
   */
  std::vector<int> getIds() const;

  // std::array<std::size_t, 5> getComponentCount(int cell_id) const;

  /**
   * @brief Dump all cells into a **DUMP** format of Phreeqc.
   *
   * @return std::map<int, std::string> Map containing the cell ID as key and
   * the exported DUMP string as value.
   */
  std::map<int, std::string> getDumpStringsPQI() const;

  /**
   * @brief Get the **DUMP** string for a given cell.
   *
   * @param cell_id Cell ID to get the **DUMP** string for.
   * @return std::string Phreeqc **DUMP** string for the given cell.
   */
  std::string getDumpStringsPQI(int cell_id) const;

  /**
   * @brief Get the Database used to initialize the PhreeqcMatrix.
   *
   * @return std::string Database string.
   */
  std::string getDatabase() const;

  /**
   * @brief Check if a cell with given ID exists in the PhreeqcMatrix.
   *
   * @param cell_id ID of the cell (user id from Phreeqc script) to check for.
   * @return true Entry exists
   * @return false Entry doesn't exist
   */
  bool checkIfExists(int cell_id) const;

  /**
   * @brief Retrieves the current PhreeqcKnobs settings.
   *
   * This function returns a copy of the PhreeqcKnobs object that contains
   * the current configuration settings for the Phreeqc instance.
   *
   * @return PhreeqcKnobs The current configuration settings.
   */
  PhreeqcKnobs getKnobs() const { return *_m_knobs; }

  /**
   * @brief Checks if the redox states are included in solution.
   *
   * This function returns a boolean value indicating whether the redox states
   * of a species are considered in the solution.
   *
   * @return true if redox states are included, false otherwise.
   */
  bool withRedox() const { return _m_with_redox; }

  // MDL
  /**
   * @brief Returns all column names of the Matrix pertaining to KINETICS
   *
   * This function returns a string vector.
   *
   * @return std::vector<std::string> Whole vector of names. Empty if no KINETICS
   * is defined
   */
   std::vector<std::string> getMatrixKinetics() const;
   
  /**
   * @brief Returns all column names of the Matrix pertaining to EQUILIBRIUM
   *
   * This function returns a string vector.
   *
   * @return std::vector<std::string> Whole vector of names. Empty if no EQUILIBRIUM
   * is defined
   */
   std::vector<std::string> getMatrixEquilibrium() const;

      
   /*
   
   @brief Returns all column names of the Matrix pertaining to
   quantities that must be transported
   
   @return std::vector<std::string> vector of names

   */
   std::vector<std::string> getMatrixTransported() const;

   /*

   @brief Returns all column names of the Matrix pertaining to
   quantities that must NOT be transported but have to be included in
   the output
   
   @return std::vector<std::string> vector of names
   */
   std::vector<std::string> getMatrixOutOnly() const;

private:
  std::map<int, std::vector<element>> _m_map;
  std::map<int, std::vector<base_names>> _m_internal_names;

  std::set<std::string> _m_surface_primaries;

  void initialize();

  void remove_NaNs();

  std::shared_ptr<IPhreeqc> _m_pqc;
  std::shared_ptr<PhreeqcKnobs> _m_knobs;

  std::string _m_database;

  bool _m_with_h0_o0;
  bool _m_with_redox;
};
