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

### Subdirectories ###
add_subdirectory(external/cdbg-builder cdbg_builder_bin)

set(SSHASH_BUILD_EXECUTABLES OFF CACHE BOOL "" FORCE)
add_subdirectory(external/sshash)

### Fulgor Target ###
add_executable(fulgor tools/fulgor.cpp)

target_include_directories(fulgor PRIVATE
        .
        external/FQFeeder/include
        external/sketch/include
        external/sketch/include/blaze
)

target_link_libraries(fulgor
        PRIVATE
        cdbg_builder_headers # provides cdbg headers.
        sshash::sshash # provides sshash headers AND the single copy of bits.
        z
        ${CMAKE_DL_LIBS}
)

if (UNIX)
    if (APPLE)
        MESSAGE(STATUS "linking with rt should not be necessary on OSX; not adding rt")
    else ()
        target_link_libraries(fulgor PRIVATE rt)
    endif ()
endif ()