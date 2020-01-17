#pragma once


#include "Tensor.h"

namespace Tensors {

	template <class T> class TensorOrder3 : public Tensor<T, 3>
	{
	public:
		TensorOrder3(unsigned int dim1 = 1, unsigned int dim2 = 1, unsigned int dim3 = 1)
			: Tensor(std::array<unsigned int, 3>{ { dim1, dim2, dim3 } })
		{
			assert(dim1);
			assert(dim2);
			assert(dim3);
		}

		T& operator()(unsigned int index1, unsigned int index2, unsigned int index3) { 
			const std::array<unsigned int, 3> indices{ { index1, index2, index3 } };

			return m_values[GetOffset(indices)];
		}

		const T& operator()(unsigned int index1, unsigned int index2, unsigned int index3) const {
			std::array<unsigned int, 3> indices{ { index1, index2, index3 } };

			return m_values[GetOffset(indices)];
		}
	};


}