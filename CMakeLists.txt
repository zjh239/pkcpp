# Remember to use Ninja as generator for modules support.

cmake_minimum_required(VERSION 3.29)

project(pkcpp
    VERSION 0.3
    DESCRIPTION "Polikit: Polyhedral analyzing toolkits."
    LANGUAGES CXX)

# Force Ninja as the build system

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_SCAN_FOR_MODULES ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")

file(GLOB SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/basic/*.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/io/*.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/func/*.cxx"
)

# find_package(OpenMP REQUIRED)

add_executable(pkcpp)

target_include_directories(pkcpp PRIVATE
    ext/tomlpp/include
)

# first compile the modules
target_sources(pkcpp
    PRIVATE
        FILE_SET modules TYPE CXX_MODULES FILES
        ${SRCS}
)

# then compile the excutable file
target_sources(pkcpp PRIVATE src/main.cxx)

option(DEBUG "Enable Thread Sanitizer" OFF)

if(DEBUG)
    add_compile_options(-g -O1 -fsanitize=thread -fno-omit-frame-pointer)
    add_link_options(-fsanitize=thread)
    message(STATUS "Building in debug mode, thread sanitizer enabled")
endif()

# if (OpenMP_CXX_FOUND)
#     target_link_libraries(pkcpp PRIVATE OpenMP::OpenMP_CXX)
# endif()
