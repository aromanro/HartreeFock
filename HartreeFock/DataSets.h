#pragma once

#include <list>
#include <vector>

class DataSets {
public:
	class DataSet {
	public:
		DataSet() : color(RGB(0, 0, 0)), lineWidth(0) {}

		std::vector<Gdiplus::PointF> points;
		COLORREF color;
		float lineWidth;

		double getXMin() const;
		double getYMin() const;
		double getXMax() const;
		double getYMax() const;

		void Draw(Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, const Gdiplus::RectF& dataRect, bool spline = true) const;

	protected:
		static double ConvertValue(double val, double valMin, double valMax, double chartMin, double chartMax);
	};

	std::list<DataSet> dataSets;

	double getXMin() const;
	double getYMin() const;
	double getXMax() const;
	double getYMax() const;

	void Draw(Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, const Gdiplus::RectF& dataRect, bool spline = true) const;
};
