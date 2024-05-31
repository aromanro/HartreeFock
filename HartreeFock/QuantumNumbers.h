#pragma once

#include <cassert>

namespace Orbitals {

	namespace QuantumNumbers {

		class QuantumNumbers
		{
		public:
			unsigned int l, m, n;

			QuantumNumbers(unsigned int L = 0, unsigned int M = 0, unsigned int N = 0) noexcept;

			char AtomicOrbital() const noexcept;

			inline unsigned int AngularMomentum() const noexcept
			{
				return l + m + n;
			}

			inline unsigned int N(unsigned int ind) const noexcept
			{
				if (0 == ind) return l;
				else if (1 == ind) return m;

				assert(2 == ind);

				return n;
			}

			inline unsigned int MaxComponentVal() const noexcept { return max(max(l, m), n); }

			inline unsigned int MinComponentVal() const noexcept { return min(min(l, m), n); }

			inline unsigned int NumOrbitals() const noexcept
			{
				return NumOrbitals(AngularMomentum());
			}

			inline unsigned int GetCanonicalIndex() const noexcept
			{
				return (m + n + 1) * (m + n + 2) / 2 - m - 1;
			}

			inline unsigned int GetTotalCanonicalIndex() const noexcept
			{
				const unsigned int L = AngularMomentum();

				const unsigned int val = L * (L * (L + 3) + 2) / 6;

				return val + GetCanonicalIndex();
			}

			inline void CanonicalIncrement() noexcept
			{
				const unsigned int L = AngularMomentum();

				if (n == L)
				{
					l = L + 1;
					m = n = 0;
				}
				else if (m > 0)
				{
					--m;
					++n;
				}
				else if (l > 0)
				{
					--l;
					m = L - l;
					n = 0;
				}
			}

			// we're using Cartesian Gaussians, there are (L + 1) * (L + 2) of them in a shell with angular momentum L
			// don't confuse it with the usual 2L+1 number
			static inline unsigned int NumOrbitals(unsigned int L) noexcept
			{
				return (L + 1) * (L + 2) / 2;
			}

			inline operator unsigned int() const noexcept { return AngularMomentum(); }

			inline QuantumNumbers& operator++() noexcept
			{
				CanonicalIncrement();

				return *this;
			}

			inline QuantumNumbers operator++(int) noexcept
			{
				QuantumNumbers temp(*this);

				operator++();

				return temp;
			}
		};

		inline bool operator<(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() < rhs.AngularMomentum(); }
		inline bool operator<(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() < rhs; }
		inline bool operator<(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() < static_cast<unsigned int>(rhs); }
		inline bool operator<(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs < rhs.AngularMomentum(); }


		inline bool operator>(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() > rhs.AngularMomentum(); }
		inline bool operator>(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() > rhs; }
		inline bool operator>(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() > static_cast<unsigned int>(rhs); }
		inline bool operator>(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs > rhs.AngularMomentum(); }


		inline bool operator<=(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() <= rhs.AngularMomentum(); }
		inline bool operator<=(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() <= rhs; }
		inline bool operator<=(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() <= static_cast<unsigned int>(rhs); }
		inline bool operator<=(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs <= rhs.AngularMomentum(); }



		inline bool operator>=(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() >= rhs.AngularMomentum(); }
		inline bool operator>=(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() >= rhs; }
		inline bool operator>=(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() >= static_cast<unsigned int>(rhs); }
		inline bool operator>=(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs >= rhs.AngularMomentum(); }






		class SQuantumNumbers : public QuantumNumbers {
		public:
			SQuantumNumbers() : QuantumNumbers(0, 0, 0) {};
		};



		class PxQuantumNumbers : public QuantumNumbers {
		public:
			PxQuantumNumbers() : QuantumNumbers(1, 0, 0) {};
		};

		class PyQuantumNumbers : public QuantumNumbers {
		public:
			PyQuantumNumbers() : QuantumNumbers(0, 1, 0) {};
		};

		class PzQuantumNumbers : public QuantumNumbers {
		public:
			PzQuantumNumbers() : QuantumNumbers(0, 0, 1) {};
		};





		class Dx2QuantumNumbers : public QuantumNumbers {
		public:
			Dx2QuantumNumbers() : QuantumNumbers(2, 0, 0) {};
		};

		class Dy2QuantumNumbers : public QuantumNumbers {
		public:
			Dy2QuantumNumbers() : QuantumNumbers(0, 2, 0) {};
		};

		class Dz2QuantumNumbers : public QuantumNumbers {
		public:
			Dz2QuantumNumbers() : QuantumNumbers(0, 0, 2) {};
		};

		class DxyQuantumNumbers : public QuantumNumbers {
		public:
			DxyQuantumNumbers() : QuantumNumbers(1, 1, 0) {};
		};

		class DxzQuantumNumbers : public QuantumNumbers {
		public:
			DxzQuantumNumbers() : QuantumNumbers(1, 0, 1) {};
		};

		class DyzQuantumNumbers : public QuantumNumbers {
		public:
			DyzQuantumNumbers() : QuantumNumbers(0, 1, 1) {};
		};




		class Fx3QuantumNumbers : public QuantumNumbers {
		public:
			Fx3QuantumNumbers() : QuantumNumbers(3, 0, 0) {};
		};

		class Fy3QuantumNumbers : public QuantumNumbers {
		public:
			Fy3QuantumNumbers() : QuantumNumbers(0, 3, 0) {};
		};

		class Fz3QuantumNumbers : public QuantumNumbers {
		public:
			Fz3QuantumNumbers() : QuantumNumbers(0, 0, 3) {};
		};

		class FxyzQuantumNumbers : public QuantumNumbers {
		public:
			FxyzQuantumNumbers() : QuantumNumbers(1, 1, 1) {};
		};

		class Fx2yQuantumNumbers : public QuantumNumbers {
		public:
			Fx2yQuantumNumbers() : QuantumNumbers(2, 1, 0) {};
		};

		class Fx2zQuantumNumbers : public QuantumNumbers {
		public:
			Fx2zQuantumNumbers() : QuantumNumbers(2, 0, 1) {};
		};

		class Fxy2QuantumNumbers : public QuantumNumbers {
		public:
			Fxy2QuantumNumbers() : QuantumNumbers(1, 2, 0) {};
		};

		class Fy2zQuantumNumbers : public QuantumNumbers {
		public:
			Fy2zQuantumNumbers() : QuantumNumbers(0, 2, 1) {};
		};

		class Fxz2QuantumNumbers : public QuantumNumbers {
		public:
			Fxz2QuantumNumbers() : QuantumNumbers(1, 0, 2) {};
		};

		class Fyz2QuantumNumbers : public QuantumNumbers {
		public:
			Fyz2QuantumNumbers() : QuantumNumbers(0, 1, 2) {};
		};
	}
}

