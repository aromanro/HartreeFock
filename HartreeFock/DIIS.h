#pragma once


#include <unsupported/Eigen/CXX11/Tensor>

// will need specialization for types different than eigen matrices!

template <typename ValueType, int maxRetained = 6, int firstEstimate = 5> class DIISStorage
{
public:
	std::list<ValueType> errors;
	std::list<ValueType> values;

	// returns true if there is enough data to estimate the next value
	bool AddValueAndError(const ValueType& value, ValueType& error)
	{
		errors.emplace_back(error);
		values.push_back(value);
		if (errors.size() > maxRetained)
		{
			errors.pop_front();
			values.pop_front();
		}

		return errors.size() >= firstEstimate;
	}
};


template <typename ValueType, int maxRetained = 6, int firstEstimate = 5> class DIIS : public DIISStorage<ValueType, maxRetained, firstEstimate>
{
public:

	double Estimate(ValueType &value) const
	{
		const size_t nrMatrices = errors.size();
		ValueType B = ValueType::Zero(nrMatrices + 1, nrMatrices + 1);

		double lastErrorEst = 0;

		auto errorIter1 = errors.begin();
		for (size_t i = 0; i < nrMatrices; ++i)
		{
			auto errorIter2 = errors.begin();

			for (size_t j = 0; j < i; ++j)
			{
				B(i, j) = B(j, i) = (*errorIter1).cwiseProduct(*errorIter2).sum();

				++errorIter2;
			}

			B(i, i) = (*errorIter1).cwiseProduct(*errorIter1).sum();

			if (i == nrMatrices - 1 || i == nrMatrices - 2) lastErrorEst += B(i, i);

			B(nrMatrices, i) = B(i, nrMatrices) = 1;

			++errorIter1;
		}

		lastErrorEst = sqrt(lastErrorEst);

		// Solve the system of linear equations

		Eigen::VectorXd CDIIS = Eigen::VectorXd::Zero(nrMatrices + 1);
		CDIIS(nrMatrices) = 1;

		CDIIS = B.colPivHouseholderQr().solve(CDIIS);

		// compute the new value

		value = ValueType::Zero(value.rows(), value.cols());

		auto iter = values.begin();
		for (size_t i = 0; i < nrMatrices; ++i, ++iter)
			value += CDIIS(i) * *iter;
			
		return lastErrorEst;
	}
};

template<int maxRetained, int firstEstimate> class DIIS<Eigen::Tensor<double, 4>, maxRetained, firstEstimate> : public DIISStorage<Eigen::Tensor<double, 4>, maxRetained, firstEstimate>
{
public:
	double Estimate(Eigen::Tensor<double, 4>& value) const
	{
		const size_t nrMatrices = errors.size();
		Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nrMatrices + 1, nrMatrices + 1);

		double lastErrorEst = 0;

		auto errorIter1 = errors.begin();
		for (size_t i = 0; i < nrMatrices; ++i)
		{
			auto errorIter2 = errors.begin();

			for (size_t j = 0; j < i; ++j)
			{
				const Eigen::Tensor<double, 0> val = (*errorIter1 * *errorIter2).sum();
				const double dval = val();
				B(i, j) = B(j, i) = dval;

				++errorIter2;
			}

			const Eigen::Tensor<double, 0> diagVal = (*errorIter1 * *errorIter1).sum();
			const double ddiagVal = diagVal();
			B(i, i) = ddiagVal;

			if (i == nrMatrices - 1 || i == nrMatrices - 2) lastErrorEst += B(i, i);

			B(nrMatrices, i) = B(i, nrMatrices) = 1;

			++errorIter1;
		}

		lastErrorEst = sqrt(lastErrorEst);

		// Solve the system of linear equations

		Eigen::VectorXd CDIIS = Eigen::VectorXd::Zero(nrMatrices + 1);
		CDIIS(nrMatrices) = 1;

		CDIIS = B.colPivHouseholderQr().solve(CDIIS);

		// compute the new value

		value.setZero();

		auto iter = values.begin();
		for (size_t i = 0; i < nrMatrices; ++i, ++iter)
			value += CDIIS(i) * *iter;

		return lastErrorEst;
	}

};

