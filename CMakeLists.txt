cmake_minimum_required(VERSION 2.8)

file(GLOB SOURCE_CUBICALGRID "*.h" "*.cpp")
add_library(libCUBICALGRID STATIC ${SOURCE_CUBICALGRID})