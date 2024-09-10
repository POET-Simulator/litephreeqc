## Time-stamp: "Last modified 2024-08-26 17:42:18 delucia"

add_executable(ipqMain 
    ipqMain.cpp
  )
add_executable(ipqEngine 
    ipqEngine.cpp
  )

# add_library(IPhreeqcPOET ${PQC_POET_SRC})
##target_link_libraries(IPhreeqcPOET PUBLIC IPhreeqc)

target_include_directories(IPhreeqcPOET PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src)


target_link_libraries(ipqMain 
    IPhreeqcPOET
)

target_link_libraries(ipqEngine 
    IPhreeqcPOET
)
