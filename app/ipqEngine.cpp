#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <iomanip>
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

void printVectorsAsTable(const std::vector<std::string> &names,
                         const std::vector<double> &values) {
  if (names.size() != values.size()) {
    std::cerr << "Error: Vectors are not of the same length." << std::endl;
    return;
  }

  std::cout << std::left << std::setw(20) << "Array Index" << std::setw(20)
            << "Name" << std::setw(20) << "Value" << std::endl;
  std::cout << std::string(40, '-') << std::endl;

  for (std::size_t i = 0; i < names.size(); ++i) {
    std::cout << std::left << std::setw(20) << std::to_string(i)
              << std::setw(20) << names[i] << std::setw(20) << values[i]
              << std::endl;
  }
}

template <class T> void printVector(const std::vector<T> &vec) {
  for (const auto &elem : vec) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
}

int main(int argc, char *argv[]) {

  if (argc < 3) {
    std::cout << "Two args needed, PHREEQC script and database" << std::endl;
    return 1;
  }

  const std::string &script = argv[1];
  const std::string &database = argv[2];
  wordexp_t exp_result;

  std::cout << "Reading script from file " << script << " using database "
            << database << std::endl;

  // expand the script and database path
  const std::string read_script = readFile(script);
  const std::string read_database = readFile(database);

  // initialize the module instance "ipqcmod"
  PhreeqcInit ipqcmod(read_database, read_script);

  // allocate tabular data
  PhreeqcInit::PhreeqcMat pqc_mat = ipqcmod.getPhreeqcMat();

  // for (std::size_t i = 0; i < nrows; ++i) {
  auto pure_names = ipqcmod.getEquilibriumNames(pqc_mat.ids[0]);
  auto totals_names = ipqcmod.getSolutionNames(pqc_mat.ids[0]);

  POETConfig MyConfig;
  MyConfig.database = read_database;
  MyConfig.input_script = read_script;
  MyConfig.cell = POETInitCell{
      // H, O and charge need to be removed. I'm not sure why I did it this way.
      std::vector<std::string>(totals_names.begin() + 3, totals_names.end()),
      {},
      {},
      {},
      pure_names,
      {},
      {}};

  std::vector<std::string> full_names = {"ID"};
  full_names.insert(full_names.end(), pqc_mat.names.begin(),
                    pqc_mat.names.end());

  PhreeqcEngine MyEngine(MyConfig);
  std::vector<double> NewInp = pqc_mat.values[0];
  NewInp.insert(NewInp.begin(), 1);

  std::cout << "Before running the cell:\n";
  printVectorsAsTable(full_names, NewInp);

  // std::for_each(NewInp.begin()+3, NewInp.end(),
  //               [](double a) { return a * 0.8;});

  NewInp[5] += 0.2; // Cl
  NewInp[7] += 0.2; // Na
  MyEngine.runCell(NewInp, 200);

  std::cout << "\n\nAfter running the cell:\n";
  printVectorsAsTable(full_names, NewInp);

  return 0;
}
