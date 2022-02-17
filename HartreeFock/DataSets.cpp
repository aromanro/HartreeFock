#include "stdafx.h"
#include "DataSets.h"


double DataSets::DataSet::getXMin() const
{
	double result = DBL_MAX;

	for (const auto& p : points)
		if (p.X < result) result = p.X;

	return result;
}

double DataSets::DataSet::getYMin() const
{
	double result = DBL_MAX;

	for (const auto& p : points)
		if (p.Y < result) result = p.Y;

	return result;
}

double DataSets::DataSet::getXMax() const
{
	double result = DBL_MIN;

	for (const auto& p : points)
		if (p.X > result) result = p.X;

	return result;
}

double DataSets::DataSet::getYMax() const
{
	double result = DBL_MIN;

	for (const auto& p : points)
		if (p.Y > result) result = p.Y;

	return result;
}

double DataSets::DataSet::ConvertValue(double val, double valMin, double valMax, double chartMin, double chartMax)
{
	return chartMin + (val - valMin) / (valMax - valMin) * (chartMax - chartMin);
}


void DataSets::DataSet::Draw(Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, const Gdiplus::RectF& dataRect, bool spline) const
{
	if (points.size() < static_cast<unsigned int>(spline ? 3 : 2)) return;
	Gdiplus::Pen pen(Gdiplus::Color::MakeARGB(255, GetRValue(color), GetGValue(color), GetBValue(color)));
	pen.SetWidth(lineWidth);

	const double valMinX = dataRect.X;
	const double valMinY = dataRect.Y;
	const double valMaxX = static_cast<double>(dataRect.X) + dataRect.Width;
	const double valMaxY = static_cast<double>(dataRect.Y) + dataRect.Height;

	const double chartMinX = boundRect.X;
	const double chartMinY = boundRect.Y;
	const double chartMaxX = static_cast<double>(boundRect.X) + boundRect.Width;
	const double chartMaxY = static_cast<double>(boundRect.Y) + boundRect.Height;


	Gdiplus::PointF* pnts = new Gdiplus::PointF[points.size()];
	int index = 0;
	for (const auto& pnt : points)
	{
		pnts[index].X = static_cast<float>(ConvertValue(pnt.X, valMinX, valMaxX, chartMinX, chartMaxX));
		pnts[index].Y = static_cast<float>(ConvertValue(pnt.Y, valMinY, valMaxY, chartMinY, chartMaxY));
		++index;
	}

	if (spline) g.DrawCurve(&pen, pnts, index);
	else g.DrawLines(&pen, pnts, index);

	delete[] pnts;
}



double DataSets::getXMin() const
{
	double result = DBL_MAX;

	for (const auto& p : dataSets)
	{
		const double val = p.getXMin();
		if (val < result) result = val;
	}

	return result;
}

double DataSets::getYMin() const
{
	double result = DBL_MAX;

	for (const auto& p : dataSets)
	{
		const double val = p.getYMin();
		if (val < result) result = val;
	}

	return result;
}

double DataSets::getXMax() const
{
	double result = DBL_MIN;

	for (const auto& p : dataSets)
	{
		const double val = p.getXMax();
		if (val > result) result = val;
	}

	return result;
}

double DataSets::getYMax() const
{
	double result = DBL_MIN;

	for (const auto& p : dataSets)
	{
		const double val = p.getYMax();
		if (val > result) result = val;
	}

	return result;
}

void DataSets::Draw(Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, const Gdiplus::RectF& dataRect, bool spline) const
{
	g.SetClip(boundRect);
	for (const auto& dataSet : dataSets)
		dataSet.Draw(g, boundRect, dataRect, spline);
}

