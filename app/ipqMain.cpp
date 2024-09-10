//  Time-stamp: "Last modified 2024-08-27 09:00:18 delucia"
#include <cstddef>
#include <linux/limits.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>
#include <wordexp.h>

#include "POETInit.hpp"
#include "PhreeqcEngine.hpp"
#include "PhreeqcInit.hpp"

static std::string readFile(const std::string &path) {
  std::string string_rpath(PATH_MAX, '\0');

  if (realpath(path.c_str(), string_rpath.data()) == nullptr) {
    throw std::runtime_error("Failed to resolve the realpath to file " + path);
  }

  std::ifstream file(string_rpath);

  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  std::stringstream buffer;
  buffer << file.rdbuf();

  return buffer.str();
}

int main(int argc, char *argv[]) {

  const std::string &script = argv[1];
  const std::string &database = argv[2];
  // wordexp_t exp_result;

  std::cout << "Reading script from file " << script << " using database "
            << database << std::endl;

  // expand the script and database path
  const std::string read_script = readFile(script);
  const std::string read_database = readFile(database);

  // initialize the module instance "ipqcmod"
  PhreeqcInit ipqcmod(read_database, read_script);

  // allocate tabular data
  PhreeqcInit::PhreeqcMat pqc_mat = ipqcmod.getPhreeqcMat();

  const std::size_t ncols = pqc_mat.names.size();
  const std::size_t nrows = pqc_mat.values.size();

  std::cout << "Script contains " << nrows << " rows and " << ncols
            << " columns" << std::endl;

  // pqc_mat.names.insert(pqc_mat.names.begin(), "ID");

  // Rcpp::NumericMatrix mat(nrows, ncols + 1);

  std::cout << std::endl << "Number of CELLS/SOLUTIONS: " << nrows << std::endl;

  for (std::size_t i = 0; i < nrows; ++i) {
    auto pure_names = ipqcmod.getEquilibriumNames(pqc_mat.ids[i]);
    auto totals_names = ipqcmod.getSolutionNames(pqc_mat.ids[i]);

    std::cout << std::endl << "SOL " << pqc_mat.ids[i] << " - Pure phases : ";
    for (std::size_t j = 0; j < pure_names.size(); j++) {
      std::cout << pure_names[j] << ", ";
    }

    std::cout << std::endl << "SOL " << pqc_mat.ids[i] << " - Total concs: ";
    for (std::size_t j = 0; j < totals_names.size(); j++) {
      std::cout << totals_names[j] << ", ";
    }
    std::cout << std::endl;
  }

  // alternatively, use the pqc_mat object
  std::cout << std::endl << "Total pqc_mat: " << std::endl;
  for (std::size_t i = 0; i < nrows; ++i) {
    std::cout << "ID: " << pqc_mat.ids[i] << " - ";
    for (std::size_t j = 0; j < ncols; ++j) {
      std::cout << pqc_mat.names[j] << ": " << pqc_mat.values[i][j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
