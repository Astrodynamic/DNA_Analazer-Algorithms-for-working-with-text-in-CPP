cmake_minimum_required(VERSION 3.5)

project(DNA_Analyzer VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
  ${CMAKE_SOURCE_DIR}/interface
  ${CMAKE_SOURCE_DIR}/algorithm
)

set(HEADERS
  algorithm/DNA_Analyzer.h
  interface/ainterface.h
  interface/interface.h
)

set(SOURCES
  algorithm/DNA_Analyzer.cc
  interface/ainterface.cc
  interface/interface.cc
)

add_executable(
  ${PROJECT_NAME}
  ${HEADERS}
  ${SOURCES}
  main.cc
)