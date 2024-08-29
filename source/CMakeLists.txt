#git clone --recurse-submodules https://github.com/seqan/seqan3.git
cmake_minimum_required (VERSION 3.8)
# set( CMAKE_VERBOSE_MAKEFILE on )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project (tree CXX)


# Find OpenMP
if(APPLE)
    set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/18.1.8/bin/clang")
    set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/18.1.8/bin/clang++")
    set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/18.1.8/lib")
    set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/18.1.8/include")

    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

# find_package(BF_LIB libbf REQUIRED)
find_package(OpenMP)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)


# Shared interface
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include)
 
# build app with seqan3
add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}_interface")

target_compile_definitions(${PROJECT_NAME} PRIVATE
    "-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")

target_compile_definitions(${PROJECT_NAME} PRIVATE
    "-DGIT_BRANCH=\"${GIT_BRANCH}\"")
include_directories(include)


# Build all subcommands
message(STATUS "build_all_subcommands------------------------------------------{")
file(GLOB files "*_main.cpp")
foreach(file ${files})
    get_filename_component(lib_we ${file} NAME_WE)
    message(STATUS "building ${lib_we}")
    add_library(${lib_we} SHARED ${lib_we}.cpp)
    target_link_libraries (${lib_we} PUBLIC "${PROJECT_NAME}_interface")
    target_link_libraries ("${PROJECT_NAME}" PUBLIC "${lib_we}")
    if(NOT APPLE)
        target_link_libraries ("${PROJECT_NAME}" PUBLIC OpenMP::OpenMP_CXX)
    endif()
    install(TARGETS ${lib_we} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endforeach()
message(STATUS "build_all_subcommands------------------------------------------}")


install(TARGETS tree DESTINATION bin)