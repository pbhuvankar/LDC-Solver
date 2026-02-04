# Source files
set(LDC_SOURCES
    BoundaryConditions.cpp
    Convection.cpp
    Diffusion.cpp
    Initialize.cpp
    Projection.cpp
    Solver.cpp
)

set(LDC_HEADERS
    LDC.hpp
)

# Serial executable
if(BUILD_SERIAL)
    add_executable(ldc_serial main.cpp ${LDC_SOURCES})
    target_include_directories(ldc_serial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    set_target_properties(ldc_serial PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    message(STATUS "Building serial executable: ldc_serial")
endif()

# OpenMP executable
if(BUILD_OPENMP AND OpenMP_CXX_FOUND)
    add_executable(ldc_openmp main.cpp ${LDC_SOURCES})
    target_include_directories(ldc_openmp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(ldc_openmp PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(ldc_openmp PRIVATE _OPENMP)
    set_target_properties(ldc_openmp PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    message(STATUS "Building OpenMP executable: ldc_openmp")
endif()

# Test executable
if(BUILD_TEST)
    add_executable(ldc_test tests.cpp ${LDC_SOURCES})
    target_include_directories(ldc_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    set_target_properties(ldc_test PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    message(STATUS "Building test executable: ldc_test")
endif()

# Install targets
if(BUILD_SERIAL)
    install(TARGETS ldc_serial DESTINATION bin)
endif()

if(BUILD_OPENMP AND OpenMP_CXX_FOUND)
    install(TARGETS ldc_openmp DESTINATION bin)
endif()

if(BUILD_TEST)
    install(TARGETS ldc_test DESTINATION bin)
endif()