//  Time-stamp: "Last modified 2025-08-01 10:54:30 delucia"
#include <iostream>
#include <iomanip>
#include <linux/limits.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>
#include <memory>
#include <cmath>

#include "PhreeqcMatrix.hpp"
#include "PhreeqcEngine.hpp"
#include "PhreeqcRunner.hpp"


std::string readFile(const std::string &path) {
  std::string string_rpath(PATH_MAX, '\0');

  if (realpath(path.c_str(), string_rpath.data()) == nullptr) {
    throw std::runtime_error(":: Failed to resolve the realpath to file " + path);
  }

  std::ifstream file(string_rpath);

  if (!file.is_open()) {
    throw std::runtime_error(":: Failed to open file: " + path);
  }

  std::stringstream buffer;
  buffer << file.rdbuf();

  return buffer.str();
}

// pretty print a vector, standard implementation from stackoverflow 
template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os << "{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, ", "));
    os << " }";
    return os;
}


int main(int argc, char *argv[]) {
    
    if (argc < 3) {
	std::cout << "::" << argv[0] << ": two args needed, script and database\n";
	return 1;
    }


    ////// INITIALISATION
    // read Script and Database and put it in a std::string
    auto script = readFile(argv[1]);
    auto db     = readFile(argv[2]);

    // Create the matrix directly from database and init script
    std::cout << ":: Creating a PhreeqcMatrix with valence states and H(0)/O(0) \n";
    PhreeqcMatrix pqc_mat1(db, script, true, true);
    
    // How many different SOLUTIONS ("CELLS") are defined in the script?
    const auto ids = pqc_mat1.getIds();
    
    int n = ids.size();
    
    std::cout << ":: Found " << n << " distinct PHREEQC problems \n";

    std::cout << ":: getSolutionsNames(): the common solutes across all problems: \n";
    const auto solutes1  = pqc_mat1.getSolutionNames();
    std::cout << solutes1 << "\n";

    // auto expmat = pqc_mat1.get();
    auto allvars = pqc_mat1.get().names;
    std::cout << ":: pqc_mat1.get().names (all names in the PhreeqcMatrix): \n";
    std::cout << allvars << "\n\n";

    std::cout << "\n-- Now the new getMatrix*() --\n\n";
    
    auto transported = pqc_mat1.getMatrixTransported(); 
    std::cout << ":: pqc_mat1.getMatrixTransported(): \n";
    std::cout << transported << "\n\n";

    auto MatNamesKin = pqc_mat1.getMatrixKinetics(); 
    std::cout << ":: pqc_mat1.getMatrixKinetics(): \n";
    std::cout << MatNamesKin << "\n\n";
    
    auto MatNamesEqui = pqc_mat1.getMatrixEquilibrium(); 
    std::cout << ":: pqc_mat1.getMatrixEquilibrium(): \n";
    std::cout << MatNamesEqui << "\n\n";

    auto outonly = pqc_mat1.getMatrixOutOnly(); 
    std::cout << ":: pqc_mat1.getMatrixOutOnly(): \n";
    std::cout << outonly << "\n\n";

    auto selout = pqc_mat1.getSelectedOutputNames(); 
    std::cout << ":: pqc_mat1.getSelectedOutputNames(): \n";
    std::cout << selout << "\n\n";
    return 0;
}

// Oneliner for rz-vm278 relative to iphreeqc/poet/test!!

// g++ testGolemRunner.cpp -o testG -Wall -I../../poet/include -I../../src -I../../src/phreeqcpp -I../../src/phreeqcpp/common -I../../src/phreeqcpp/PhreeqcKeywords -lIPhreeqc -lIPhreeqcPOET -L../../bbuild/ -L../../bbuild/poet
