enable_testing()

add_executable(
    poet_test
    testPhreeqcEngine.cpp
    testPhreeqcMatrix.cpp
    utils.cpp
)

target_link_libraries(
    poet_test
    IPhreeqcPOET
    GTest::gtest_main
)

target_include_directories(poet_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# read file and store in variable
file(REAL_PATH "${PROJECT_SOURCE_DIR}/database/phreeqc.dat" POET_PHREEQCDAT_DB)
file(REAL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/barite_db.dat" POET_BARITE_DB)
file(REAL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/barite_het.pqi" POET_BARITE_PQI)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testInput.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/testInput.hpp")

include(GoogleTest)
gtest_discover_tests(poet_test)