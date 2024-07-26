#pragma once

#include "Tensor.h"

namespace Tensors {

	template <class T = double> class TensorOrder4 : public Tensor<T, 4>
	{
	public:
		TensorOrder4(size_t dim1 = 1, size_t dim2 = 1, size_t dim3 = 1, size_t dim4 = 1)
			: Tensor(std::array<size_t,4>{ dim1, dim2, dim3, dim4 })
		{
			assert(dim1);
			assert(dim2);
			assert(dim3);
			assert(dim4);
		}

		T& operator()(size_t index1, size_t index2, size_t index3, size_t index4) {
			const std::array<size_t, 4> indices{ index1, index2, index3, index4 };

			return m_values[GetOffset(indices)];
		}
		
		const T& operator()(size_t index1, size_t index2, size_t index3, size_t index4) const {
			const std::array<size_t, 4> indices{ index1, index2, index3, index4 };

			return m_values[GetOffset(indices)];
		}
	};


}