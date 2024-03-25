cmake_minimum_required(VERSION 3.5)
project(FULGOR)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY omp)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

MESSAGE(STATUS "Compiling for processor: " ${CMAKE_HOST_SYSTEM_PROCESSOR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  MESSAGE(STATUS "Compiling with flags: -march=native -mbmi2 -msse4.2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2 -mpopcnt") # for hardware popcount and pdep
endif()

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")

  MESSAGE(STATUS "Compiling with flags: ${CMAKE_CXX_FLAGS}")

  # Flags for PTHash:
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # for multi-threading

  if (FULGOR_USE_SANITIZERS)
    MESSAGE(STATUS "Using sanitizers. Compiling with flags: -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

include_directories(.)
include_directories(external/sshash/external/pthash)
include_directories(external/sshash/external/pthash/external/essentials/include)
include_directories(external/sshash)
include_directories(external/FQFeeder/include)
include_directories(external/sketch/include)
include_directories(external/sketch/include/blaze)

## GGCAT ####

add_custom_target(
        ggcat_cpp_api
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api
)
set(GGCAT_CPP_BINDINGS ${CMAKE_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cpp_bindings.a)
set(GGCAT_CXX_INTEROP ${CMAKE_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cxx_interop.a)
set(GGCAT_API ${CMAKE_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_api.a)


## Fulgor ###

set(Z_LIB_SOURCES
        external/sshash/include/gz/zip_stream.cpp
)

set(SSHASH_SOURCES
        external/sshash/include/dictionary.cpp
        external/sshash/include/info.cpp
        external/sshash/include/builder/build.cpp
)

# Create a static lib for SSHash
add_library(sshash_static STATIC
        ${Z_LIB_SOURCES}
        ${SSHASH_SOURCES}
)

add_executable(fulgor tools/fulgor.cpp)
add_dependencies(fulgor ggcat_cpp_api sshash_static)
target_link_libraries(fulgor
        sshash_static
        z
        ${GGCAT}
        ${GGCAT_API}
        ${GGCAT_CPP_BINDINGS}
        ${GGCAT_CXX_INTEROP}
        ${CMAKE_DL_LIBS}
)

if (UNIX)
  if (APPLE)
    MESSAGE(STATUS "linking with rt should not be necessary on OSX; not adding rt")
  else()
    target_link_libraries(fulgor rt)
  endif()
endif()