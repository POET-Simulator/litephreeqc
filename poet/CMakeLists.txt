file(GLOB_RECURSE PQC_POET_SRC CONFIGURE_DEPENDS 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(IPhreeqcPOET ${PQC_POET_SRC})
target_link_libraries(IPhreeqcPOET PUBLIC IPhreeqc)
target_include_directories(IPhreeqcPOET PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Set C++20 standard
target_compile_features(IPhreeqcPOET PUBLIC cxx_std_20)

if (BUILD_TESTING AND STANDALONE_BUILD)
    add_subdirectory(test)
endif()