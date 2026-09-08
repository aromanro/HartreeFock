#pragma once

#include <Eigen/Core>

// Eigen allocates differently for SSE and AVX. Catch incompatible library
// clients at link time, before a destructor uses the wrong allocator.
#ifdef _MSC_VER
#define HF_ABI_STRING_IMPL(value) #value
#define HF_ABI_STRING(value) HF_ABI_STRING_IMPL(value)
#pragma detect_mismatch("HartreeFock.Eigen.DynamicAlignment", HF_ABI_STRING(EIGEN_DEFAULT_ALIGN_BYTES))
#pragma detect_mismatch("HartreeFock.Eigen.StaticAlignment", HF_ABI_STRING(EIGEN_MAX_STATIC_ALIGN_BYTES))
#pragma detect_mismatch("HartreeFock.Eigen.SystemAllocator", HF_ABI_STRING(EIGEN_MALLOC_ALREADY_ALIGNED))
#undef HF_ABI_STRING
#undef HF_ABI_STRING_IMPL
#endif
