//  Time-stamp: "Last modified 2025-06-11 15:30:00 delucia"
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
	std::cout << "mm ::" << argv[0] << ": two args needed, script and database\n";
	return 1;
    }


    ////// INITIALISATION
    // read Script and Database and put it in a std::string
    auto script = readFile(argv[1]);
    auto db     = readFile(argv[2]);

    // Create the matrix directly from database and init script
    PhreeqcMatrix pqc_mat(db, script, false, false);
    // How many different SOLUTIONS ("CELLS") are defined in the script?
    const auto ids = pqc_mat.getIds();
    
    int n = ids.size();
    
    std::cout << ":: Found " << n << " distinct PHREEQC problems \n";
    std::cout << ids << "\n";

    const auto solutes  = pqc_mat.getSolutionNames();
    std::cout << ":: These are the common solutes across all the " << n << " problems: \n";
    std::cout << solutes << "\n";
    
    // iterate on the ids (THEY start at 1!!)
    for (const auto & i : ids) {
	auto pphases  = pqc_mat.getEquilibriumNames(i);
	if (!pphases.empty()) {
	    std::cout << ":: Equilibrium phases [" << (int) i <<  "]: \n";
	    std::cout << pphases << "\n";
	}
	
	auto kinetics = pqc_mat.getKineticsNames(i);
	if (!kinetics.empty()) {
	    std::cout << ":: Kinetics [" << i <<  "]: \n";
	    std::cout << kinetics << "\n";
	}
    }

    // The exported data type holds the matrix in a "STL format" with
    // a "header" of names and their accompanying values. The values
    // are stored in a row-major order per default.
    auto exported_mat = pqc_mat.get();
    // Get the total number of solutes
    const int len = exported_mat.names.size();
    // Get the values as reference to modify them in place
    std::vector<double> &cell_values = exported_mat.values;


    // optional SUBSET the matrix (i.e., the unique ids defined in
    // golem map as input)
    //

    const auto subsetted_pqc_mat1 = pqc_mat.subset({1});
    const auto subsetted_pqc_mat2 = pqc_mat.subset({2});

    const auto stl_mat1 = subsetted_pqc_mat1.get();
    const auto stl_mat2 = subsetted_pqc_mat2.get();
    
    const auto matrix_values1 = stl_mat1.values;
    const auto num_columns1 = stl_mat1.names.size();
    const auto spec_names1 = stl_mat1.names;

    const auto matrix_values2 = stl_mat2.values;
    const auto num_columns2 = stl_mat2.names.size();
    const auto spec_names2 = stl_mat2.names;

    std::cout <<  "Subset 1: \n";
    for (std::size_t i = 0; i < num_columns1; ++i) {
	std::cout << i << ") " << spec_names1[i] << ", ";
    }

    std::cout <<  "\n";
    std::cout <<  "Subset 2: \n";
    for (std::size_t i = 0; i < num_columns2; ++i) {
	std::cout << i << ") " << spec_names2[i] << ", ";
    }

    std::cout <<  "\n";

    // // container to pass in/out
    // std::vector<std::vector<double>> simulationInOut;

    // // // grid cells
    // // const std::size_t num_cells = 10;
    // // const std::size_t half_cells = 5;

    // // // copy the values to the InOut vector. We replicate cell 1 
    // // for (std::size_t index = 0; index < num_cells; ++index) {
    // // 	if (index < half_cells) {
    // // 	    simulationInOut.push_back(std::vector<double>(
    // //         matrix_values.begin(), matrix_values.begin() + num_columns));
    // // 	} else {
    // // 	    simulationInOut.push_back(std::vector<double>(
    // //         matrix_values.begin() + num_columns, matrix_values.end()));
    // // 	}
    // // }

    // // copy the values from the matrix to the InOut vector
    // for (std::size_t index = 0; index < 2; ++index) {
    // 	simulationInOut.push_back(std::vector<double>(
    //     matrix_values.begin()+ num_columns*index, matrix_values.begin() + num_columns*(index +1)));
    // }


    // std::cout << "\n:: Values in the PhreeqcMatrix after initialisation: \n\n";

    // // output the values from the matrix
    // for (std::size_t cell_index = 0; cell_index < simulationInOut.size(); ++cell_index) {
    // 	std::cout << "Grid element: " << cell_index << " \n";
    // 	for (std::size_t spec = 0; spec < num_columns; ++spec) {
    // 	    std::cout << spec_names[spec] << "=" << simulationInOut[cell_index][spec] << ", ";
    // 	}
    // 	std::cout << " \n";

    // }

    // // now we compute 1 timestep
    // const double timestep = 1000.;

    // // create the runner
    // PhreeqcRunner runner(subsetted_pqc_mat);
    // // run
    // runner.run(simulationInOut, timestep);

    // // output the values returned after simulation
    // std::cout << "\n:: Values in InOut after one time step: \n\n";
    // for (std::size_t cell_index = 0; cell_index < simulationInOut.size(); ++cell_index) {
    // 	std::cout << "Grid element: " << cell_index << " \n";
    // 	for (std::size_t spec = 0; spec < num_columns; ++spec) {
    // 	    std::cout << spec_names[spec] << "=" << simulationInOut[cell_index][spec] << ", ";
    // 	}
    // 	std::cout << " \n";

    // }

    
    return 0;
}

// Oneliner for rz-vm278 relative to iphreeqc/poet/test!!

// g++ testGolemRunner.cpp -o testG -Wall -I../../poet/include -I../../src -I../../src/phreeqcpp -I../../src/phreeqcpp/common -I../../src/phreeqcpp/PhreeqcKeywords -lIPhreeqc -lIPhreeqcPOET -L../../bbuild/ -L../../bbuild/poet
