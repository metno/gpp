cmake_minimum_required(VERSION 3.0)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(gridpp)
enable_language(CXX)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ENABLE_TESTS "build test suite" OFF)
if (NOT CMAKE_BUILD_TYPE)
   message(STATUS "CMAKE_BUILD_TYPE is not set, switching to RELEASE")
   set(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
      "Choose the type of build, options are: NONE DEBUG RELEASE."
      FORCE)
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

if (NOT GTEST_DIR)
   message(STATUS "GTEST_DIR is not set, assuming ${CMAKE_CURRENT_SOURCE_DIR}/gtest/")
   set(GTEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gtest/" CACHE STRING
      "Set the path to gtest tarball/source directory"
      FORCE)
endif()
if(NOT IS_ABSOLUTE "${GTEST_DIR}")
   if (CMAKE_VERSION VERSION_LESS 3.4.3)
      message(FATAL_ERROR "GTEST_DIR path should be absolute")
   else()
      get_filename_component(GTEST_DIR "${GTEST_DIR}" ABSOLUTE
         BASE_DIR ${CMAKE_BINARY_DIR})
   endif()
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare")
if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
   set(CMAKE_CXX_FLAGS_DEBUG "-O0 -fPIC -g -pg -rdynamic --coverage")
   add_definitions("-DDEBUG")
elseif(CMAKE_BUILD_TYPE MATCHES "RELEASE")
   find_package(OpenMP REQUIRED)
   set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${OpenMP_CXX_FLAGS}")
endif()

find_package(Boost  REQUIRED)
find_package(Armadillo 6.5 REQUIRED)
find_package(NetCDF REQUIRED)
find_package(GSL REQUIRED)

function(list_sources sources)
   file(GLOB CORESRC "src/*.cpp")
   file(GLOB CALSRC  "src/Calibrator/*.cpp")
   file(GLOB DOWNSRC "src/Downscaler/*.cpp")
   file(GLOB FILESRC "src/File/*.cpp")
   file(GLOB PARSRC  "src/ParameterFile/*.cpp")

   set(${sources}
      ${CORESRC}
      ${CALSRC}
      ${FILESRC}
      ${DOWNSRC}
      ${PARSRC} PARENT_SCOPE)
endfunction(list_sources)

function(add_gridpp_library target lib_type)
   add_library(${target} ${lib_type} ${SRC})
   target_include_directories(${target} PUBLIC "${Boost_INCLUDE_DIR}")
   target_include_directories(${target} PUBLIC "${NETCDF_INCLUDE_DIRS}")
   target_include_directories(${target} PUBLIC "${GSL_INCLUDE_DIRS}")
   target_include_directories(${target} PUBLIC "${ARMADILLO_INCLUDE_DIRS}")

   set_target_properties(${target} PROPERTIES OUTPUT_NAME gridpp)
endfunction()

list_sources(SRC)
set(DRVSRC  "src/Driver/Gridpp.cpp")
set(TRAINRC "src/Driver/Train.cpp")

add_gridpp_library(gridpp-static STATIC)
add_executable(gridpp-bin   ${DRVSRC} )
add_executable(gridpp-train ${TRAINRC})
set_target_properties(gridpp-bin   PROPERTIES OUTPUT_NAME gridpp)
set_target_properties(gridpp-train PROPERTIES OUTPUT_NAME gridpp_train)

foreach(item gridpp-bin gridpp-train)
   target_link_libraries(${item} gridpp-static)

   target_link_libraries(${item} "${Boost_LIBRARIES}")
   target_link_libraries(${item} "${NETCDF_LIBRARIES}")
   target_link_libraries(${item} "${GSL_LIBRARIES}")
   target_link_libraries(${item} "${ARMADILLO_LIBRARIES}")
endforeach(item)

if (ENABLE_TESTS)
   include(Testing)
endif()

function(lookup_program result program)
   message(STATUS "Looking for ${program}")
   find_program(${result} "${program}")
   if(${result})
      message(STATUS "Looking for ${program} - found")
   else()
      message(STATUS "Looking for ${program} - not found")
   endif()
endfunction()

lookup_program(WC_TOOL wc)
if(NOT WC_TOOL)
   message(STATUS "Line count statistics disabled")
else()
   add_custom_target(count
      COMMAND wc src/*.h src/*.cpp src/*/*.h src/*/*.cpp -l | tail -1
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
lookup_program(CTAGS_TOOL ctags)
if(NOT CTAGS_TOOL)
   message(STATUS "Tags generation disabled")
else()
   add_custom_target(tags
      COMMAND ctags -R --c++-kinds=+pl --fields=+iaS --extra=+q -f tags src/*.h src/*.cpp src/*/*.h src/*/*.cpp
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
lookup_program(DOXYGEN_TOOL doxygen)
if(NOT DOXYGEN_TOOL)
   message(STATUS "Documentation generation disabled")
else()
   configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doxygen/config.in" "doxygen/config" @ONLY)
   add_custom_target(doxygen
      COMMAND doxygen doxygen/config
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

install(TARGETS gridpp-bin gridpp-train
   RUNTIME DESTINATION bin)
