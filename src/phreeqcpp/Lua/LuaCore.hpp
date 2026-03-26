#pragma once

#include <lua.hpp>
#include <string>

class Phreeqc;

class LuaCore {
private:
  lua_State *L;

public:
  LuaCore(Phreeqc *phreeqc);
  ~LuaCore();

  int run(const std::string &chunk_name);

  int parse_chunk(const std::string &script, const std::string &chunk_name);

  double pop_value();

  const std::string Phreeqc_Lua_Expose_Name = "Pqc";

private:
  int initialize_runtime();
};
