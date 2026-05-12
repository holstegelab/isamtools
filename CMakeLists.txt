cmake_minimum_required(VERSION 3.16)
project(isamtools C)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-Wall -Wextra -Wno-unused-parameter -Wshadow -fvisibility=hidden)
add_compile_definitions(_POSIX_C_SOURCE=200809L _GNU_SOURCE)

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET htslib)

option(ISAM_LMDB_FETCH "Allow CMake to fetch LMDB source when not provided" OFF)
set(ISAM_LMDB_SOURCE_DIR "" CACHE PATH
    "Path to LMDB libraries/liblmdb source directory containing mdb.c and midl.c")

if(ISAM_LMDB_SOURCE_DIR)
    if(NOT EXISTS "${ISAM_LMDB_SOURCE_DIR}/mdb.c")
        message(FATAL_ERROR "ISAM_LMDB_SOURCE_DIR must contain mdb.c")
    endif()
    set(LMDB_SOURCE_DIR "${ISAM_LMDB_SOURCE_DIR}")
elseif(ISAM_LMDB_FETCH)
    include(FetchContent)
    FetchContent_Declare(
        lmdb_source
        URL https://github.com/LMDB/lmdb/archive/LMDB_0.9.31.tar.gz
        URL_HASH MD5=9d7f059b1624d0a4d4b2f1781d08d600
    )
    FetchContent_MakeAvailable(lmdb_source)
    FetchContent_GetProperties(lmdb_source)
    set(LMDB_SOURCE_DIR "${lmdb_source_SOURCE_DIR}/libraries/liblmdb")
else()
    message(FATAL_ERROR
        "LMDB is required; set ISAM_LMDB_SOURCE_DIR or ISAM_LMDB_FETCH=ON")
endif()

add_library(lmdb_static OBJECT
    "${LMDB_SOURCE_DIR}/mdb.c"
    "${LMDB_SOURCE_DIR}/midl.c"
)
target_include_directories(lmdb_static PUBLIC "${LMDB_SOURCE_DIR}")
target_compile_definitions(lmdb_static PRIVATE MDB_MAXKEYSIZE=0)

set(ISAM_SOURCES
    src/common.c
    src/backend_lmdb.c
    src/hash.c
    src/rebuild.c
    src/store.c
    src/view.c
    src/cmd_index.c
    src/cmd_import.c
    src/cmd_check.c
    src/cmd_sort.c
    src/cmd_view.c
    src/cmd_annotate.c
    src/cmd_update.c
    src/cmd_rehint.c
    src/api.c
    src/htsisam_dispatch.c
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(isam STATIC ${ISAM_SOURCES})
target_include_directories(isam
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(isam PUBLIC PkgConfig::HTSLIB lmdb_static)

add_executable(isamtools src/main.c)
target_link_libraries(isamtools PRIVATE isam)

# Drop-in compat shared library for htslib consumers. Re-compiles the same
# isam sources (including htsisam_dispatch.c) into a self-contained .so so
# downstream packages just need `-lhtsisam -lhts`. Visibility is controlled
# by ISAM_EXPORT / explicit `extern` declarations on htsisam_*.
add_library(htsisam SHARED ${ISAM_SOURCES} $<TARGET_OBJECTS:lmdb_static>)
target_include_directories(htsisam
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${LMDB_SOURCE_DIR}
)
target_link_libraries(htsisam PUBLIC PkgConfig::HTSLIB)
set_target_properties(htsisam PROPERTIES C_VISIBILITY_PRESET default)

install(TARGETS isamtools RUNTIME DESTINATION bin)
install(TARGETS htsisam   LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
