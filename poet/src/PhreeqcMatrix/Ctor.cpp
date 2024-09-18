#include "IPhreeqc.hpp"
#include "PhreeqcMatrix.hpp"

#include <Phreeqc.h>
#include <Solution.h>
#include <cmath>
#include <string>

PhreeqcMatrix::PhreeqcMatrix(const std::string &database,
                             const std::string &input_script)
    : _m_database(database) {
  this->_m_pqc = std::make_shared<IPhreeqc>();

  this->_m_pqc->LoadDatabaseString(database.c_str());
  this->_m_pqc->RunString(input_script.c_str());

  this->initialize();
}

PhreeqcMatrix::PhreeqcMatrix(const PhreeqcMatrix &other)
    : _m_map(other._m_map), _m_internal_names(other._m_internal_names),
      _m_surface_primaries(other._m_surface_primaries), _m_pqc(other._m_pqc),
      _m_database(other._m_database) {}

PhreeqcMatrix::PhreeqcMatrix(PhreeqcMatrix &&other)
    : _m_map(other._m_map), _m_internal_names(other._m_internal_names),
      _m_surface_primaries(other._m_surface_primaries), _m_pqc(other._m_pqc),
      _m_database(other._m_database) {}

PhreeqcMatrix &PhreeqcMatrix::operator=(const PhreeqcMatrix &other) {
  _m_map = other._m_map;
  _m_internal_names = other._m_internal_names;
  _m_surface_primaries = other._m_surface_primaries;
  _m_pqc = other._m_pqc;
  _m_database = other._m_database;

  return *this;
}

PhreeqcMatrix &PhreeqcMatrix::operator=(PhreeqcMatrix &&other) {
  _m_map = other._m_map;
  _m_internal_names = other._m_internal_names;
  _m_surface_primaries = other._m_surface_primaries;
  _m_pqc = other._m_pqc;
  _m_database = other._m_database;

  return *this;
}