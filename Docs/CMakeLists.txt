cmake_minimum_required (VERSION 3.10)
project(PeleLMeX-Docs NONE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Sphinx REQUIRED)
#find_package(Doxygen)

add_custom_target(sphinx ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/html
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/source/landing
                                             ${CMAKE_BINARY_DIR}/html
  COMMAND ${SPHINX_EXECUTABLE} -M html "${CMAKE_SOURCE_DIR}/source/manual"
                                       "${CMAKE_BINARY_DIR}/html/manual"
                                       -c ${CMAKE_SOURCE_DIR}/source/manual
  )
