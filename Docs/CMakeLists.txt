cmake_minimum_required (VERSION 3.10)
project(PeleLMeX-Docs NONE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Sphinx REQUIRED)

add_custom_target(sphinx ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/landing
                                             ${CMAKE_CURRENT_BINARY_DIR}/html
  COMMAND ${SPHINX_EXECUTABLE} -M html "${CMAKE_CURRENT_SOURCE_DIR}/sphinx/manual"
                                       "${CMAKE_CURRENT_BINARY_DIR}/html/manual"
                                       -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/manual
                                       -W --keep-going -n
  )
