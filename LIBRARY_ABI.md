# Numerical library and Eigen ABI

The x64 Release UI crash after the library split was caused by mixed Eigen SIMD
settings. The UI used AVX2 (32-byte alignment and Eigen's custom aligned
allocator); the numerical library used the default SSE configuration (16-byte
alignment and the system allocator). Eigen matrices cross this boundary, so a
destructor could free a buffer using the wrong allocation convention.

The Windows application log recorded heap corruption, `0xC0000374`. A small
AVX2 client allocating `RestrictedHartreeFock::h` and letting the library destroy
it reproduced that error. The same client with matching settings passed.

Native x64 projects now import `HartreeFock/HartreeFockBuild.props` for AVX2,
including the numerical library, GUI, DMRG client, and regression client.
MSVC x64 CMake builds also use AVX2. These builds require an AVX2-capable CPU.
`HartreeFockEigenABI.h`, included by the public mathematical headers, emits
MSVC link-time checks for dynamic/static alignment and allocator compatibility.
An intentionally incompatible SSE client was verified to fail with `LNK2038`.

The user's additional CCSD/CIS/spin-orbital source moves remain in the library;
the CMake source list has been updated to include the same 28 numerical files.

## Regression

The native solution includes `HartreeFockNumericsTests`. From this directory:

```powershell
./tests/x64/Release/HartreeFockNumericsTests.exe ./HartreeFock/sto3g.txt
```

Alternatively, after building the sibling DMRG CMake project:

```powershell
ctest --test-dir ../dmrgqc/build -C Release -R hartreefock_boundary --output-on-failure
```

The regression exercises allocation/deallocation in both directions across the
library boundary, H2 RHF/UHF calculations, MP2, CIS, and CCSD initialization.
It checks the HF energy and finite post-HF outputs; it is a boundary/crash
regression, not a comprehensive validation of the post-HF algorithms.

The native Release GUI and library rebuilt successfully. Both native and CMake
boundary tests passed, as did the DMRG core and standard CTest suites. The UI
computation click itself was not automated; validation used the numerical API
and allocation boundary that reproduced the crash.
