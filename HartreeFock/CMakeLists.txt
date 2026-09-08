# Reusable numerical target, also consumed by dmrgqc via add_subdirectory().
if(NOT MSVC)
  message(FATAL_ERROR "HartreeFock numerical sources currently include MFC. Use MSVC with MFC installed.")
endif()
set(EIGEN3_INCLUDE_DIR "" CACHE PATH "Directory containing Eigen/Core")
if(EIGEN3_INCLUDE_DIR)
  if(NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
    message(FATAL_ERROR "EIGEN3_INCLUDE_DIR must contain Eigen/Core")
  endif()
elseif(NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 5 CONFIG REQUIRED)
endif()

set(CMAKE_MFC_FLAG 2)
add_library(hartreefock_numerics STATIC
  Atom.cpp
  Basis.cpp
  BoysFunction.cpp
  BoysFunctions.cpp
  ChemUtils.cpp
  Coefficients.cpp
  ContractedGaussianOrbital.cpp
  GaussianKinetic.cpp
  GaussianMoment.cpp
  GaussianNuclear.cpp
  GaussianOrbital.cpp
  GaussianOverlap.cpp
  GaussianTwoElectrons.cpp
  HartreeFockAlgorithm.cpp
  IntegralsRepository.cpp
  MathUtils.cpp
  Molecule.cpp
  PrimitiveShell.cpp
  QuantumMatrix.cpp
  QuantumNumbers.cpp
  RestrictedHartreeFock.cpp
  RestrictedCCSD.cpp
  RestrictedCCSDTriples.cpp
  RestrictedCCSD_F.cpp
  RestrictedCCSD_W.cpp
  RestrictedConfigurationIInteractionSingles.cpp
  SpinOrbitalsElectronElectronIntegralsRepository.cpp
  UnrestrictedHartreeFock.cpp
)
add_library(HartreeFock::Numerics ALIAS hartreefock_numerics)
target_compile_features(hartreefock_numerics PUBLIC cxx_std_20)
target_include_directories(hartreefock_numerics PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
if(EIGEN3_INCLUDE_DIR)
  target_include_directories(hartreefock_numerics SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}")
else()
  target_link_libraries(hartreefock_numerics PUBLIC Eigen3::Eigen)
endif()
target_compile_definitions(hartreefock_numerics PRIVATE _AFXDLL UNICODE _UNICODE)
target_compile_options(hartreefock_numerics PUBLIC /bigobj /FC)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  # Propagate AVX2 to clients: Eigen's allocator/alignment ABI must match.
  target_compile_options(hartreefock_numerics PUBLIC /arch:AVX2)
endif()

if(BUILD_TESTING)
  add_executable(hartreefock_boundary_tests NumericsBoundaryTests.cpp)
  target_link_libraries(hartreefock_boundary_tests PRIVATE HartreeFock::Numerics)
  add_test(NAME hartreefock_boundary COMMAND hartreefock_boundary_tests "${CMAKE_CURRENT_SOURCE_DIR}/sto3g.txt")
endif()
