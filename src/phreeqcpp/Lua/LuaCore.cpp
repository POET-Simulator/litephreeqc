#include "LuaInitializer.hpp"
#include "NameDouble.h"

#include <iostream>
#include <lua.hpp>

static int error_read_only(lua_State *L) {
  return luaL_error(L, "read-only variable");
}

LuaCore::LuaCore(Phreeqc *phreeqc) : L(luaL_newstate()) {
  luaL_openlibs(this->L);

  luaJIT_setmode(this->L, 0, LUAJIT_MODE_ENGINE | LUAJIT_MODE_ON);

  Phreeqc **phreeqc_lua_ptr =
      (Phreeqc **)lua_newuserdata(this->L, sizeof(Phreeqc *));
  *phreeqc_lua_ptr = phreeqc;

  initialize_runtime();

  lua_setglobal(L, Phreeqc_Lua_Expose_Name.c_str());
}

int LuaCore::initialize_runtime() {
  luaL_newmetatable(L, "Meta.Pqc");
  luaL_setfuncs(L, pqc_func_vec, 0);
  lua_pushvalue(L, -1);
  lua_setfield(L, -2, "__index");

  lua_setmetatable(L, -2);

  return 1;
}

LuaCore::~LuaCore() { lua_close(this->L); }

int LuaCore::run(const std::string &chunk_name) {
  lua_getglobal(this->L, chunk_name.c_str());
  int ret_value = lua_pcall(this->L, 0, 1, 0);

  if (ret_value != LUA_OK) {
    std::cerr << lua_tostring(this->L, -1) << std::endl;
  }
  return ret_value;
}

int LuaCore::parse_chunk(const std::string &script,
                         const std::string &chunk_name) {
  std::string in_script = script;
  // remove leading ';' if present
  if (script[0] == ';') {
    in_script = script.substr(1);
  }
  int status = luaL_loadstring(this->L, in_script.c_str());
  if (status != LUA_OK) {
    std::cout << lua_tostring(this->L, -1) << std::endl;
    return status;
  }
  lua_setglobal(this->L, chunk_name.c_str());

  return LUA_OK;
}

double LuaCore::pop_value() { return lua_tonumber(this->L, -1); }
