enable_testing()

add_executable(
    poet_test
    testPhreeqcMatrix.cpp
)

target_link_libraries(
    poet_test
    IPhreeqcPOET
    GTest::gtest_main
)

target_include_directories(poet_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# read file and store in variable
file(READ "${PROJECT_SOURCE_DIR}/database/phreeqc.dat" POET_PHREEQCDAT_DB)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/barite_db.dat" POET_BARITE_DB)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/barite_het.pqi" POET_BARITE_PQI)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/testPhreeqcMatrix.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/testPhreeqcMatrix.hpp")

include(GoogleTest)
gtest_discover_tests(poet_test)