# PROJECT SETUP
# ===========================================

cmake_minimum_required(VERSION 3.15)

project(Polaron_MC
    VERSION 1.0.0
    DESCRIPTION " Monte Carlo simulation of polarons"
    LANGUAGES CXX
)

cmake_policy(SET CMP0077 NEW) # Options in target_link_libraries
#cmake_policy(SET CMP0074 NEW)

# ==========================================
# PROJECT OPTIONS & CONFIGURATION
# ==========================================
# User options (can be set via -D flag)
#option(BUILD_SHARED_LIBS "Build shared libraries instead of static" OFF)
#option(BUILD_TESTS "Build unit tests" OFF)
#option(BUILD_EXAMPLES "Build example programs" OFF)
#option(ENABLE_SANITIZERS "Enable address/undefined behavior sanitizers" OFF)

# Configure compilation database for tools like clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# PROJECT STRUCTURE
set(SOURCE_DIR src)
set(INCLUDE_DIR include)
set(THIRDPARTY_DIR thirdparty)

# COMPILER WARNINGS AND FLAGS
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Werror -O3)
endif()

# SOURCE FILES
set(MAIN_SOURCES
    "${SOURCE_DIR}/main.cpp"
    "${SOURCE_DIR}/Diagram.cpp"
    "${SOURCE_DIR}/GreenFuncNph.cpp"
    "${SOURCE_DIR}/GreenFuncNphBands.cpp"
    "${SOURCE_DIR}/utils/MC_Benchmarking.cpp"
    "${SOURCE_DIR}/utils/progressbar.cpp"
)

# HEADERS
set(MAIN_HEADERS
    "${INCLUDE_DIR}/Diagram.hpp"
    "${INCLUDE_DIR}/GreenFuncNph.hpp"
    "${INCLUDE_DIR}/GreenFuncNphBands.hpp"
    "${INCLUDE_DIR}/utils/MC_Benchmarking.hpp"
    "${INCLUDE_DIR}/utils/progressbar.hpp"
    "${INCLUDE_DIR}/utils/computational_methods.hpp"
    "${INCLUDE_DIR}/utils/MC_data_structures.hpp"
)

# Create main executable
add_executable(${PROJECT_NAME}
    ${MAIN_SOURCES}
    ${MAIN_HEADERS}
)

# EXTERNAL LIBRARIES
find_package(Eigen3 3.3.0 REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)


# 4. CREATE LIBRARIES FROM YOUR CODE
# ============================================
# Library 1: Math operations (static library)
add_library(Diagram STATIC 
    "src/Diagram.cpp"
    "include/Diagram.hpp"
)

# Include directories
target_include_directories(${PROJECT_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Compiler definitions
target_compile_definitions(${PROJECT_NAME}
    PRIVATE
        $<$<CONFIG:Debug>:DEBUG_MODE=1>
    PUBLIC
        PROJECT_NAME=\"${PROJECT_NAME}\"
        PROJECT_VERSION=\"${PROJECT_VERSION}\"
)

# Print configuration summary
message(STATUS "=========================================")
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Source files: ${MAIN_SOURCES}")
message(STATUS "=========================================")

