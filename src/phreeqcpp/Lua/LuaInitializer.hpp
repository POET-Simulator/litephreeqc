#include <Phreeqc.h>
#include <lua.hpp>

template <typename Func, typename... Args>
int push_phreeqc_data(lua_State *L, Func rate_func, Args... args) {
  auto phreeqc = static_cast<Phreeqc **>(lua_touserdata(L, 1));
  lua_pushnumber(L, ((*phreeqc)->*rate_func)(args...));
  return 1;
}

static int err_no_newindex(lua_State *L) {
  return luaL_error(L, "No newindex");
}

static const struct luaL_Reg pqc_func_vec[] = {
    {"M",
     [](lua_State *L) -> int {
       return push_phreeqc_data(L, &Phreeqc::get_rate_m);
     }},
    {"TIME",
     [](lua_State *L) -> int {
       return push_phreeqc_data(L, &Phreeqc::get_rate_time);
     }},
    {"SI",
     [](lua_State *L) -> int {
       const char *name = lua_tostring(L, 2);
       return push_phreeqc_data(L, &Phreeqc::get_saturation_index, name);
     }},
    {"ACT",
     [](lua_State *L) -> int {
       const char *name = lua_tostring(L, 2);
       return push_phreeqc_data(L, &Phreeqc::get_activity, name);
     }},
    {"TK",
     [](lua_State *L) -> int {
       return push_phreeqc_data(L, &Phreeqc::get_temperature_kelvin);
     }},
    {"SR",
     [](lua_State *L) -> int {
       const char *name = lua_tostring(L, 2);
       return push_phreeqc_data(L, &Phreeqc::get_saturation_ratio, name);
     }},
    {"PARM",
     [](lua_State *L) -> int {
       std::size_t index = lua_tointeger(L, 2);
       return push_phreeqc_data(L, &Phreeqc::get_script_parm, index);
     }},
    {"__newindex", err_no_newindex},
    {NULL, NULL}};
