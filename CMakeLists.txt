cmake_minimum_required(VERSION 3.16...4.0)
project(FULGOR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mno-avx512f")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2 -mpopcnt")
endif ()

if (UNIX)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    endif ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

    if (FULGOR_USE_SANITIZERS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    endif ()

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
    endif ()
endif ()

MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "Compiling for processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
MESSAGE(STATUS "Compiling with flags:${CMAKE_CXX_FLAGS}")

set(CITYHASH_SOURCES ${CMAKE_SOURCE_DIR}/external/sshash/external/cityhash/cityhash.cpp)

### Fulgor ###
add_subdirectory(external/cdbg-builder cdbg_builder_bin)
add_executable(fulgor tools/fulgor.cpp)

# Cleaned up: Only containing FULGOR's direct internal dependencies
target_include_directories(fulgor PRIVATE
        .
        external/sshash
        external/sshash/external/pthash/include
        external/sshash/external/pthash/external/bits/include
        external/sshash/external/pthash/external/fastmod
        external/sshash/external/pthash/external/xxHash
        external/sshash/external/pthash/external/mm_file/include
        external/sshash/external/pthash/external/bits/external/essentials/include
        external/FQFeeder/include
        external/sketch/include
        external/sketch/include/blaze
)

target_link_libraries(fulgor
        PRIVATE
        z
        ${CMAKE_DL_LIBS}
        ${CITYHASH_SOURCES}
        cdbg_builder_headers
)

if (UNIX)
    if (APPLE)
        MESSAGE(STATUS "linking with rt should not be necessary on OSX; not adding rt")
    else ()
        target_link_libraries(fulgor PRIVATE rt)
    endif ()
endif ()