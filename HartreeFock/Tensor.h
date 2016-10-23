#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

namespace Tensors {

	template<class T, unsigned int O> class Tensor
	{
	protected:
		std::vector<T> m_values;
	    std::array<unsigned int, O> m_dims;

		unsigned int GetOffset(const std::array<unsigned int, O>& indices) const {
			unsigned int result = 0;

			for (unsigned int i = 0; i < m_dims.size(); ++i)
				result = result * m_dims[i] + indices[i];

			return result;
		}
	public:
		Tensor(const std::array<unsigned int, O>& dims) : m_dims(dims)
		{
			m_values.resize(GetSize());
		}

		Tensor(const Tensor& other) : m_values(other.m_values), m_dims(other.m_dims) {}
		
		Tensor(Tensor&& other) noexcept {
			m_dims.swap(other.m_dims);
			m_values.swap(other.m_values);
		}

		Tensor& operator=(const Tensor& other)
		{
			Tensor temp(other);
			*this = std::move(temp);

			return *this;
		}

		Tensor& operator=(Tensor&& other) noexcept
		{
			m_dims.swap(other.m_dims);
			m_values.swap(other.m_values);

			return *this;
		}

		unsigned int GetSize() const { return std::accumulate(m_dims.begin(), m_dims.end(), 1, std::multiplies<unsigned int>()); }

		unsigned int GetDim(unsigned int index) const { assert(index < m_dims.size());  return m_dims[index]; }

		void Clear() 
		{ 
			for (unsigned int i = 0; i < m_dims.size(); ++i) m_dims[i] = 1;

			m_values.resize(GetSize());
		}
	};

}

