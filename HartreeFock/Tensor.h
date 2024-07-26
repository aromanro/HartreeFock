#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

namespace Tensors {

	template<class T, size_t O> class Tensor
	{
	protected:
		std::vector<T> m_values;
		std::array<size_t, O> m_dims;
		size_t m_sz;

		size_t GetOffset(const std::array<size_t, O>& indices) const {
			size_t result = 0;

			for (size_t i = 0; i < m_dims.size(); ++i)
				result = result * m_dims[i] + indices[i];

			return result;
		}

	public:
		Tensor(const std::array<size_t, O>& dims) : m_dims(dims), m_sz(0)
		{
			m_values.resize(GetSize());
		}

		Tensor(const Tensor& other) : m_values(other.m_values), m_dims(other.m_dims), m_sz(other.m_sz) {}

		Tensor(Tensor&& other) noexcept {
			m_dims.swap(other.m_dims);
			m_values.swap(other.m_values);
			m_sz = other.m_sz;
		}

		virtual ~Tensor() = default;

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
			m_sz = other.m_sz;

			return *this;
		}
		
		size_t GetSize()
		{
			if (!m_sz) m_sz = std::accumulate(m_dims.begin(), m_dims.end(), 1, std::multiplies<size_t>());

			return m_sz;
		}

		size_t GetDim(size_t index) const { assert(index < m_dims.size());  return m_dims[index]; }

		void Clear()
		{
			for (size_t i = 0; i < m_dims.size(); ++i) m_dims[i] = 1;
			m_sz = 0;
			m_values.resize(GetSize());
		}
	};

}

