
set(LPQC_SOURCE_FILES
    src/Engine.cpp
    src/Runner.cpp 
    src/Knobs.cpp
    #Wrappers
    src/Wrapper/EquilibriumWrapper.cpp
    src/Wrapper/EquilibriumCompWrapper.cpp
    src/Wrapper/ExchangeWrapper.cpp
    src/Wrapper/ExchangeCompWrapper.cpp
    src/Wrapper/KineticWrapper.cpp
    src/Wrapper/KineticCompWrapper.cpp
    src/Wrapper/SolutionWrapper.cpp
    src/Wrapper/SurfaceWrapper.cpp
    src/Wrapper/SurfaceCompWrapper.cpp 
    src/Wrapper/SurfaceChargeWrapper.cpp
    #PhreeqcMatrix 
    src/PhreeqcMatrix/Access.cpp
    src/PhreeqcMatrix/Ctor.cpp
    src/PhreeqcMatrix/Init.cpp
    src/PhreeqcMatrix/Misc.cpp
)

add_library(litephreeqc ${LPQC_SOURCE_FILES})
target_link_libraries(litephreeqc PUBLIC IPhreeqc)
target_include_directories(litephreeqc PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Set C++20 standard
target_compile_features(litephreeqc PUBLIC cxx_std_20)

if (BUILD_TESTING AND STANDALONE_BUILD)
    enable_testing()

    set(LPQC_TEST_SOURCE_FILES
        test/testPhreeqcEngine.cpp
        test/testPhreeqcMatrix.cpp
        test/testPhreeqcRunner.cpp
        test/testPhreeqcKnobs.cpp
        test/utils.cpp
        test/IPhreeqcReader.cpp
    )

    add_executable(litephreeqc_test ${LPQC_TEST_SOURCE_FILES})

    target_link_libraries(
        litephreeqc_test
        litephreeqc
        GTest::gtest_main
    )

    target_include_directories(litephreeqc_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

    # read file and store in variable
    file(REAL_PATH "${PROJECT_SOURCE_DIR}/database/phreeqc.dat" POET_PHREEQCDAT_DB)
    file(REAL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/test/barite_db.dat" POET_BARITE_DB)
    file(REAL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/test/barite_het.pqi" POET_BARITE_PQI)

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test/testInput.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/testInput.hpp")

    include(GoogleTest)
    gtest_discover_tests(litephreeqc_test)
endif()

add_executable(golemrunner test/testGolemRunner.cpp)
target_link_libraries(golemrunner litephreeqc)

add_executable(testGetters     test/testGetters.cpp)
target_link_libraries(testGetters litephreeqc)
