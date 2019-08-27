function(setup_googletest includes libs)
   include(ExternalProject)

   if(IS_DIRECTORY "${GTEST_DIR}")
      if(EXISTS "${GTEST_DIR}/src" AND
         IS_DIRECTORY "${GTEST_DIR}/src")
         # It seems that we have old gtest without
         # install target defined.
         set(GTEST_HAS_INSTALL_TARGET FALSE)
      else()
         set(GTEST_HAS_INSTALL_TARGET TRUE)
      endif()
   elseif(NOT GTEST_HAS_INSTALL_TARGET)
      message(FATAL_ERROR
         "Unable to guess whether gtest has install target or not "
         "provide -DGTEST_HAS_INSTALL_TARGET=(TRUE|FALSE)")
   endif()

   if(GTEST_HAS_INSTALL_TARGET)
      externalproject_add(gtest
         URL ${GTEST_DIR}
         INSTALL_DIR ${CMAKE_BINARY_DIR}/auxiliary
         CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
            -DCMAKE_INSTALL_LIBDIR:PATH=lib
         CMAKE_CACHE_ARGS
            "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
            "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
         )
      externalproject_get_property(gtest install_dir)

      set(${includes} "${install_dir}/include" PARENT_SCOPE)
      set(${libs} "-L${install_dir}/lib/ -lgtest" PARENT_SCOPE)
   else()
      externalproject_add(gtest
         URL ${GTEST_DIR}
         CMAKE_CACHE_ARGS
            "-DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}"
            "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
         INSTALL_COMMAND ""
         )
      externalproject_get_property(gtest source_dir)
      externalproject_get_property(gtest binary_dir)

      set(${includes} "${source_dir}/include" PARENT_SCOPE)
      set(${libs} "-L${binary_dir} -lgtest" PARENT_SCOPE)
   endif()
endfunction()

enable_testing()

find_package(Threads REQUIRED)
setup_googletest(GTEST_INCLUDE_DIRS GTEST_LIBRARIES)

add_gridpp_library(gridpp-shared SHARED)
add_custom_target(all-tests)
add_custom_target(test-data
   COMMAND ${CMAKE_COMMAND} -E copy_directory
      ${CMAKE_CURRENT_SOURCE_DIR}/testing
      ${CMAKE_BINARY_DIR}/testing)

file(GLOB TESTS  "src/Testing/*.cpp")
foreach(test_src ${TESTS})
   get_filename_component(test ${test_src} NAME_WE)
   set(test "${test}.exe")
   add_executable(${test} ${test_src})
   add_dependencies(${test} gtest)

   target_include_directories(${test} PUBLIC "${GTEST_INCLUDE_DIRS}")
   target_link_libraries(${test} gridpp-shared)

   target_link_libraries(${test} "${Boost_LIBRARIES}")
   target_link_libraries(${test} "${NETCDF_LIBRARIES}")
   target_link_libraries(${test} "${GSL_LIBRARIES}")
   target_link_libraries(${test} "${ARMADILLO_LIBRARIES}")

   target_link_libraries(${test} "${GTEST_LIBRARIES}")
   target_link_libraries(${test} "${CMAKE_THREAD_LIBS_INIT}")

   add_test(NAME ${test} COMMAND ${test})
   add_dependencies(${test} test-data)
   add_dependencies(all-tests ${test})
endforeach(test_src)
