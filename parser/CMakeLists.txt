cmake_minimum_required(VERSION 3.25)
project(parser)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

message(STATUS "===============================================")
message(STATUS "csd_parser ${PROJECT_VERSION}")
message(STATUS "System: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "===============================================")

message(STATUS "Searching for libpq...")


set(LIBPQ_HINTS "")

if(APPLE)
    # macOS
    execute_process(
        COMMAND brew --prefix libpq
        OUTPUT_VARIABLE _brew_libpq
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(_brew_libpq)
        list(APPEND LIBPQ_HINTS "${_brew_libpq}")
        message(STATUS "  Homebrew libpq: ${_brew_libpq}")
    else()
        list(APPEND LIBPQ_HINTS
            "/Applications/Postgres.app/Contents/Versions/latest"
            "/opt/homebrew/opt/libpq"
            "/usr/local/opt/libpq")
    endif()

elseif(WIN32)
    # Windows
    file(GLOB _pg_dirs
        "C:/Program Files/PostgreSQL/*"
        "C:/PostgreSQL/*")
    foreach(_d ${_pg_dirs})
        list(APPEND LIBPQ_HINTS "${_d}")
    endforeach()
    message(STATUS "  Windows PostgreSQL dirs: ${_pg_dirs}")

else()
    # Linux
    list(APPEND LIBPQ_HINTS
        "/usr"
        "/usr/local"
        "/usr/lib/postgresql")
endif()

find_path(LIBPQ_INCLUDE_DIR
    NAMES libpq-fe.h
    HINTS ${LIBPQ_HINTS}
    PATH_SUFFIXES include include/postgresql include/pgsql
    DOC "Path to libpq headers")

find_library(LIBPQ_LIBRARY
    NAMES pq libpq
    HINTS ${LIBPQ_HINTS}
    PATH_SUFFIXES lib lib/x64 lib/x86
    DOC "Path to libpq library")

if(NOT LIBPQ_INCLUDE_DIR)
    message(FATAL_ERROR
    "Didn't find libpq headers (libpq-fe.h)."
    "Please install the PostgreSQL client library"
    "and set CMAKE_PREFIX_PATH to its installation directory.")
endif()

if (NOT LIBPQ_LIBRARY)
    message(FATAL_ERROR
    "Didn't find libpq library (libpq.so, libpq.dylib, or libpq.dll)."
    "Please install the PostgreSQL client library"
    "and set CMAKE_PREFIX_PATH to its installation directory.")
endif()

message(STATUS "Found libpq headers: ${LIBPQ_INCLUDE_DIR}")
message(STATUS "Found libpq library: ${LIBPQ_LIBRARY}")

if(WIN32)
    get_filename_component(_libpq_lib_dir "${LIBPQ_LIBRARY}" DIRECTORY)
    get_filename_component(_libpq_root "${_libpq_lib_dir}/.." ABSOLUTE)
    if(EXISTS "${_libpq_root}/bin/libpq.dll")
        message(STATUS "  DLL: ${_libpq_root}/bin/libpq.dll")
    endif()
endif()

set(SRC
    csd_parser.cpp
    pgutils.cpp
)

set(HEADERS
    pgutils.hpp
    cif.hpp
)

add_executable(csd_parser ${SRC})

target_sources(csd_parser PRIVATE
    ${HEADERS}
)

target_include_directories(csd_parser PRIVATE
    ${LIBPQ_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(csd_parser PRIVATE
    ${LIBPQ_LIBRARY}
)