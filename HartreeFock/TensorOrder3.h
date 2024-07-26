#pragma once


#include "Tensor.h"

namespace Tensors {

	template <class T = double> class TensorOrder3 : public Tensor<T, 3>
	{
	public:
		TensorOrder3(size_t dim1 = 1, size_t dim2 = 1, size_t dim3 = 1)
			: Tensor(std::array<size_t, 3>{ dim1, dim2, dim3 })
		{
			assert(dim1);
			assert(dim2);
			assert(dim3);
		}

		T& operator()(size_t index1, size_t index2, size_t index3) {
			const std::array<size_t, 3> indices{ index1, index2, index3 };

			return m_values[GetOffset(indices)];
		}

		const T& operator()(size_t index1, size_t index2, size_t index3) const {
			std::array<size_t, 3> indices{ index1, index2, index3 };

			return m_values[GetOffset(indices)];
		}
	};


}