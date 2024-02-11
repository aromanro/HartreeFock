#include "stdafx.h"
#include "Chart.h"

Chart::Chart()
	: X(this, true), Y(this, false)
{
}

// will calculate from data if called
int Chart::GetNumTicks() const
{
	return 10;
}

int Chart::GetNumBigTicks() const
{
	return 5;
}


float Chart::GetLabelHeight(bool XAxis) const
{
	if (XAxis) return static_cast<float>(min(chartRect.Height() / 10., maxLabelHeight));

	const int nrticks = Y.GetNumBigTicks();

	return static_cast<float>(min(static_cast<double>(chartRect.Height()) / nrticks, maxLabelHeight));
}

float Chart::GetLabelWidth(bool XAxis) const
{
	if (XAxis)
	{
		const int nrticks = X.GetNumBigTicks();
		const double ticklen = static_cast<double>(chartRect.Width()) / nrticks;

		return static_cast<float>(ticklen);
	}

	return  static_cast<float>(chartRect.Width() / 10.);
}


std::list<CString> Chart::GetXLabels() const
{
	std::list<CString> l;

	const float xmin = static_cast<float>((XAxisMin == DBL_MIN) ? dataSets.getXMin() : XAxisMin);
	const float xmax = static_cast<float>((XAxisMax == DBL_MAX) ? dataSets.getXMax() : XAxisMax);

	if (isinf(xmin) || isinf(xmax) || isinf(-xmin) || isinf(-xmax)) return l;

	const int ticks = X.GetNumBigTicks();

	float interval = abs (xmax - xmin) / ticks;

	for (int i = (drawStartTickX ? 0 : 1); i <= ticks; ++i)
	{
		CString str;
		str.Format(L"%.2f", static_cast<double>(xmin) + static_cast<double>(interval)*i);
		l.push_back(str);
	}

	return l;
}


std::list<CString> Chart::GetYLabels() const
{
	std::list<CString> l;

	const float ymin = static_cast<float>((YAxisMin == DBL_MIN) ? dataSets.getYMin() : YAxisMin);
	const float ymax = static_cast<float>((YAxisMax == DBL_MAX) ? dataSets.getYMax() : YAxisMax);

	if (isinf(ymin) || isinf(ymax) || isinf(-ymin) || isinf(-ymax)) return l;

	const int ticks = Y.GetNumBigTicks();

	const float interval = abs (ymax - ymin) / ticks;


	for (int i = (drawStartTickY ? 0 : 1); i <= ticks; ++i)
	{
		CString str;
		str.Format(L"%.2f", static_cast<double>(ymin) + static_cast<double>(interval)*i);
		l.push_back(str);
	}

	return l;
}


int Chart::GetNeededFontSize(CString& str, const Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect) const
{
	double maxSize = 128;
	double minSize = 4;
	
	Gdiplus::PointF origin(boundRect.X, boundRect.Y);
	Gdiplus::RectF bounds;

	while (maxSize - minSize >= 1.)
	{
		double mid = (maxSize + minSize) / 2.;
		Gdiplus::Font font(L"Arial", static_cast<float>(mid));

		g.MeasureString(str, str.GetLength(), &font, origin, &bounds);

		if (bounds.Width <= boundRect.Width && bounds.Height <= boundRect.Height)
			minSize = mid;
		else
			maxSize = mid;
	}

	return static_cast<int>(floor(minSize));
}

void Chart::DrawText(CString &str, Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, Gdiplus::StringAlignment align, float fontSize)
{
	if (fontSize == 0) fontSize = static_cast<float>(GetNeededFontSize(str, g, boundRect));
	Gdiplus::Font font(L"Arial", fontSize);

	Gdiplus::SolidBrush brush((Gdiplus::ARGB)Gdiplus::Color::Black);
	Gdiplus::StringFormat format;
	format.SetAlignment(align);

	g.DrawString(str, str.GetLength(), &font, boundRect, &format, &brush);
}


void Chart::Draw(const CDC* pDC, CRect& rect)
{
	// save rectangle position
	int posX = rect.left;
	int posY = rect.top;

	// adjust rectangle
	rect.bottom -= rect.top;
	rect.right -= rect.left;
	rect.left = 0;
	rect.top = 0;

	// draw in a buffer
	Gdiplus::Bitmap buffer(rect.Width(), rect.Height());
	Gdiplus::Graphics g(&buffer);

	if (antialias) g.SetSmoothingMode(Gdiplus::SmoothingMode::SmoothingModeAntiAlias);

	Gdiplus::SolidBrush brush((Gdiplus::ARGB)Gdiplus::Color::White);
	g.FillRectangle(&brush, rect.left, rect.top, rect.Width(), rect.Height());

	int titleHeight = min(rect.Height() / 10, maxTitleHeight);
	Gdiplus::RectF boundRect(static_cast<float>(rect.left + rect.Width() / 10), static_cast<float>(rect.top + rect.Height() / 20), static_cast<float>(rect.Width() - 2 * rect.Width() / 10), static_cast<float>(titleHeight));

	DrawText(title, g, boundRect);

	DrawAxis(g, rect, titleHeight);

	DrawData(g, rect);

	// now draw the buffer on the screen
	Gdiplus::Graphics gdi(pDC->GetSafeHdc());
	gdi.DrawImage(&buffer, posX, posY);
}

void Chart::DrawAxis(Gdiplus::Graphics& g, CRect& rect, int titleHeight)
{
	// make room for title and labels
	int leftSide = min(max(rect.Width() / 20, 70), 100);
	int top = max(titleHeight * 2, rect.Height() / 10);
	rect.DeflateRect(leftSide, top, rect.Width() / 20, rect.Height() / 10);

	chartRect = rect;

	float fontSize = GetLabelFontSize(g);

	// draw horizontal axis
	g.TranslateTransform(static_cast<float>(rect.left), static_cast<float>(rect.bottom));
	Gdiplus::Point zero(0, 0);
	X.Draw(g, zero, rect.Width(), rect.Height(), fontSize);

	// vertical
	g.ScaleTransform(1, -1);
	g.RotateTransform(90.);
	Y.Draw(g, zero, rect.Height(), rect.Width(), fontSize);

	g.ResetTransform();

	DrawXLabel(g, rect);
	DrawYLabel(g, rect, leftSide);

	g.ResetTransform();
	g.TranslateTransform(static_cast<float>(rect.left), static_cast<float>(rect.bottom));
	g.ScaleTransform(1, -1);
}


float Chart::GetLabelFontSize(Gdiplus::Graphics& g)
{
	// try to find out the labels font size
	Gdiplus::RectF labelBound(0, 0, GetLabelWidth(true), GetLabelHeight(true));

	float fontSize = static_cast<float>(maxLabelHeight);
	std::list<CString> labels = X.GetLabels();
	for (auto label : labels)
		fontSize = min(fontSize, static_cast<float>(GetNeededFontSize(label, g, labelBound)));

	labelBound.Width = GetLabelWidth(false);
	labelBound.Height = GetLabelHeight(false);
	labels = Y.GetLabels();
	for (auto label : labels)
		fontSize = min(fontSize, static_cast<float>(GetNeededFontSize(label, g, labelBound)));

	return fontSize;
}


void Chart::DrawXLabel(Gdiplus::Graphics& g, const CRect& rect)
{
	if (XAxisLabel.GetLength())
	{
		Gdiplus::RectF localLabelBound(static_cast<float>(rect.left), static_cast<float>(rect.top + rect.Height() + rect.Height() / 12.), static_cast<float>(rect.Width()), min(static_cast<float>(min(rect.Width(), rect.Height()) / 12.), maxAxisLabelHeight));
		DrawText(XAxisLabel, g, localLabelBound);
	}
}

void Chart::DrawYLabel(Gdiplus::Graphics& g, const CRect& rect, int leftSide)
{
	if (YAxisLabel.GetLength())
	{
		g.TranslateTransform(static_cast<float>(rect.left), static_cast<float>(rect.bottom));
		g.RotateTransform(-90.);
		Gdiplus::RectF localLabelBound(static_cast<float>(rect.Height() / 10.), -static_cast<float>(leftSide), static_cast<float>(rect.Height()), min(static_cast<float>(min(rect.Width(), rect.Height()) / 12.), maxAxisLabelHeight));
		DrawText(YAxisLabel, g, localLabelBound);
	}
}

void Chart::DrawData(Gdiplus::Graphics& g, const CRect& rect)
{
	Gdiplus::RectF boundRect(0, 0, static_cast<float>(rect.Width()), static_cast<float>(rect.Height()));
	Gdiplus::RectF dataRect;

	dataRect.X = static_cast<float>((XAxisMin == DBL_MIN) ? dataSets.getXMin() : XAxisMin);
	dataRect.Y = static_cast<float>((YAxisMin == DBL_MIN) ? dataSets.getYMin() : YAxisMin);

	dataRect.Width = static_cast<float>(((XAxisMax == DBL_MAX) ? dataSets.getXMax() : XAxisMax) - dataRect.X);
	dataRect.Height = static_cast<float>(((YAxisMax == DBL_MAX) ? dataSets.getYMax() : YAxisMax) - dataRect.Y);

	if (dataRect.Width > 0 && dataRect.Height > 0)
		dataSets.Draw(g, boundRect, dataRect, useSpline);
}



void Chart::AddDataSet(const double *dataX, const double *dataY, unsigned int len, float lineWidth, COLORREF color)
{
	if (0 == len) return;

	dataSets.dataSets.emplace_back();
	DataSets::DataSet &data = dataSets.dataSets.back();
	data.points.reserve(len);
	data.color = color;
	data.lineWidth = lineWidth;

	Gdiplus::PointF point;
	for (unsigned int i = 0; i < len; ++i)
	{
		point.X = static_cast<float>(dataX[i]);
		point.Y = static_cast<float>(dataY[i]);
		data.points.push_back(point);
	}
}

void Chart::AddDataSet(const std::vector<std::pair<double, double>> *data, float lineWidth, COLORREF color)
{
	if (!data || !data->size()) return;

	dataSets.dataSets.emplace_back();
	DataSets::DataSet &dataset = dataSets.dataSets.back();
	unsigned int len = static_cast<unsigned int>(data->size());
	dataset.points.reserve(len);
	dataset.color = color;
	dataset.lineWidth = lineWidth;

	Gdiplus::PointF point;
	for (const auto &p : *data)
	{
		point.X = static_cast<float>(p.first);
		point.Y = static_cast<float>(p.second);
		dataset.points.push_back(point);
	}
}


void Chart::AddDataSet(double start, double step, const double *dataY, unsigned int len, float lineWidth, COLORREF color)
{
	if (0 == len) return;

	dataSets.dataSets.emplace_back();
	DataSets::DataSet &data = dataSets.dataSets.back();
	data.points.reserve(len);
	data.color = color;
	data.lineWidth = lineWidth;

	Gdiplus::PointF point;
	for (unsigned int i = 0; i < len; ++i)
	{
		point.X = static_cast<float>(start + step * i);
		point.Y = static_cast<float>(dataY[i]);
		data.points.push_back(point);
	}
}

void Chart::AddDataSlice(double Xval, const double *slice, unsigned int len)
{
	if (0 == len) return;

	Gdiplus::PointF point;

	unsigned int i = 0;
	for (auto &set : dataSets.dataSets)
	{
		point.X = static_cast<float>(Xval);
		point.Y = static_cast<float>(slice[i]);
		set.points.push_back(point);

		++i;
		if (i == len || i > 120) break;
	}

	for (; i < min(len, 120); ++i)
	{
		point.X = static_cast<float>(Xval);
		point.Y = static_cast<float>(slice[i]);
		dataSets.dataSets.emplace_back();
		dataSets.dataSets.back().points.push_back(point);
	}
}


void Chart::clear()
{
	dataSets.dataSets.clear();

	// clear the labels
	std::list<CString> l;

	X.SetLabels(l);
	Y.SetLabels(l);
}
