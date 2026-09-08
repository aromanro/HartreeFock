cmake_minimum_required(VERSION 3.24)
project(HartreeFock LANGUAGES CXX)
include(CTest)

# The MFC GUI remains available through HartreeFock.sln. This build exposes the
# numerical library and its boundary regression without a dependency on DMRG.
add_subdirectory(HartreeFock)
