#pragma once

#include "Tensor.h"

namespace Tensors {

	template <class T> class TensorOrder4 : public Tensor<T, 4>
	{
	public:
		TensorOrder4(unsigned int dim1 = 1, unsigned int dim2 = 1, unsigned int dim3 = 1, unsigned int dim4 = 1)
			: Tensor(std::array<unsigned int,4>{ { dim1, dim2, dim3, dim4 } })
		{
			assert(dim1);
			assert(dim2);
			assert(dim3);
			assert(dim4);
		}

		T& operator()(unsigned int index1, unsigned int index2, unsigned int index3, unsigned int index4) { 
			const std::array<unsigned int, 4> indices{ {index1, index2, index3, index4} };

			return m_values[GetOffset(indices)];
		}
		
		constexpr const T& operator()(unsigned int index1, unsigned int index2, unsigned int index3) const {
			const std::array<unsigned int, 4> indices{ { index1, index2, index3, index4 } };

			return m_values[GetOffset(indices)];
		}
	};


}