#include <IPhreeqcPOET.hpp>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

static inline std::vector<std::string>
unionStringVectors(const std::vector<std::string> &a,
                   const std::vector<std::string> &b) {
  std::vector<std::string> result;

  std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                 std::back_inserter(result));

  return result;
}

static std::vector<double>
createConcVector(const std::vector<std::string> &conc_names,
                 const std::vector<double>::const_iterator &conc_it,
                 const std::vector<std::string> &new_names) {
  std::vector<double> conc_vec;

  for (const auto &i : new_names) {
    auto it = std::find(conc_names.begin(), conc_names.end(), i);

    auto conc_index = std::distance(conc_names.begin(), it);

    if (it != conc_names.end()) {
      conc_vec.push_back(*(conc_it + conc_index));
    } else {
      conc_vec.push_back(std::numeric_limits<double>::quiet_NaN());
    }
  }

  return conc_vec;
}

void IPhreeqcPOET::parseInit() {

  std::map<int, std::pair<essential_names, std::vector<double>>> init_values;

  for (const auto &[id, val] : this->PhreeqcPtr->Get_Rxn_solution_map()) {
    // A key less than zero indicates an internal solution
    if (id < 0) {
      continue;
    }

    const essential_names curr_conc_names = this->dump_essential_names(id);

    auto pair = std::make_pair(
        curr_conc_names, this->get_essential_values(id, curr_conc_names[0]));
    raw_initials[id] = pair;
  }
}

IPhreeqcPOET::essential_names IPhreeqcPOET::union_raws(const RawMap &raws) {
  IPhreeqcPOET::essential_names names;
  for (const auto &[sol_id, val] : raws) {
    for (std::size_t i = 0; i < names.size(); i++) {
      names[i] = unionStringVectors(names[i], val.first[i]);
    }

    for (auto &specie : names[POET_EXCH]) {
      specie = this->subExchangeName(specie);
    }
  }

  return names;
}

std::vector<std::vector<double>>
IPhreeqcPOET::conc_from_essentials(const RawMap &raws,
                                   const essential_names &names) {
  std::vector<std::vector<double>> values;

  for (const auto &[key, val] : raws) {
    std::vector<double> combined_conc_vec;

    for (std::size_t i = 0, offset = 0; i < names.size(); i++) {
      std::vector<double> union_vec =
          createConcVector(val.first[i], val.second.begin() + offset, names[i]);

      combined_conc_vec.insert(combined_conc_vec.end(), union_vec.begin(),
                               union_vec.end());

      offset += val.first[i].size();
    }

    values.push_back(combined_conc_vec);
  }

  return values;
}

IPhreeqcPOET::PhreeqcMat IPhreeqcPOET::getPhreeqcMat() {
  std::vector<std::vector<double>> values;

  const IPhreeqcPOET::essential_names ess_names =
      this->union_raws(this->raw_initials);

  const std::vector<std::vector<double>> conc_values =
      this->conc_from_essentials(this->raw_initials, ess_names);

  std::vector<std::string> conc_names;

  for (const auto &i : ess_names) {
    conc_names.insert(conc_names.end(), i.begin(), i.end());
  }

  std::vector<int> ids(this->raw_initials.size());

  std::transform(this->raw_initials.begin(), this->raw_initials.end(),
                 ids.begin(), [](const auto &pair) { return pair.first; });

  return {.names = conc_names, .ids = ids, .values = conc_values};
}

IPhreeqcPOET::ModulesArray
IPhreeqcPOET::getModuleSizes(const std::vector<int> &cell_ids) {
  ModulesArray module_sizes;
  RawMap raws;

  for (const auto &id : cell_ids) {
    raws[id] = this->raw_initials[id];
  }

  const IPhreeqcPOET::essential_names ess_names = this->union_raws(raws);

  std::transform(ess_names.begin(), ess_names.end(), module_sizes.begin(),
                 [](const auto &mod) { return mod.size(); });

  return module_sizes;
}

// IPhreeqcPOET::PhreeqcMat
// IPhreeqcPOET::getPhreeqcMat(const std::vector<int> &use_ids) {
//   std::vector<std::vector<double>> values;
//   RawMap raws;

//   for (const auto &id : use_ids) {
//     raws[id] = this->raw_initials[id];
//   }

//   const IPhreeqcPOET::essential_names ess_names = this->union_raws(raws);

//   const std::vector<std::vector<double>> conc_values =
//       this->conc_from_essentials(raws, ess_names);

//   std::vector<std::string> conc_names;

//   for (const auto &i : ess_names) {
//     conc_names.insert(conc_names.end(), i.begin(), i.end());
//   }

//   return {.names = conc_names, .ids = use_ids, .values = conc_values};
// }