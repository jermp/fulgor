cmake_minimum_required(VERSION 3.5...4.0)
project(FULGOR)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2 -mpopcnt")
endif()

if (UNIX)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

  if (FULGOR_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()

endif()

MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Compiling for processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
MESSAGE(STATUS "Compiling with flags:${CMAKE_CXX_FLAGS}")

include_directories(.)
include_directories(external/sshash/external/pthash/include)
include_directories(external/sshash/external/pthash/external/bits/include)
include_directories(external/sshash/external/pthash/external/fastmod)
include_directories(external/sshash/external/pthash/external/xxHash)
include_directories(external/sshash/external/pthash/external/mm_file/include)
include_directories(external/sshash/external/pthash/external/bits/external/essentials/include)
include_directories(external/sshash)
include_directories(external/FQFeeder/include)
include_directories(external/sketch/include)
include_directories(external/sketch/include/blaze)

### GGCAT ####

add_custom_target(
  ggcat_cpp_api
  COMMAND make
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api
)
set(GGCAT_CPP_BINDINGS ${CMAKE_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cpp_bindings.a)
set(GGCAT_CXX_INTEROP ${CMAKE_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cxx_interop.a)
set(GGCAT_API ${CMAKE_SOURCE_DIR}/external/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_api.a)


### Fulgor ###

add_executable(fulgor tools/fulgor.cpp)
add_dependencies(fulgor ggcat_cpp_api)
target_link_libraries(fulgor
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
