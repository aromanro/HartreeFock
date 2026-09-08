# Build the numerical library with CMake

Requires MSVC with MFC installed, CMake 3.24 or newer, and Eigen 5. Native x64
builds use AVX2 and require an AVX2-capable CPU. The CMake target propagates the
same Eigen headers and AVX2 settings to its clients to preserve allocation ABI.

From the HartreeFock repository directory:

```powershell
cmake -S . -B build-cmake -G "Visual Studio 18 2026" -A x64 -DEIGEN3_INCLUDE_DIR=C:/LIBs/eigen
cmake --build build-cmake --config Release
ctest --test-dir build-cmake -C Release --output-on-failure
```

The library is `build-cmake/HartreeFock/Release/hartreefock_numerics.lib`.
Use `-DBUILD_TESTING=OFF` at configuration time to build only the library, or
build the `hartreefock_numerics` target explicitly. An installed Eigen package
can be selected using `Eigen3_DIR` instead of `EIGEN3_INCLUDE_DIR`.

Other CMake projects can add the `HartreeFock` source subdirectory and link
`HartreeFock::Numerics`. The sibling DMRG project does this directly: its CMake
file no longer duplicates the 28-file numerical source list. The standalone
build has no dependency on the DMRG project.

The MFC GUI remains built through `HartreeFock.sln`; this CMake build covers the
numerical library and its allocation/computation boundary regression.
