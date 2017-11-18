#include "stdafx.h"
#include "Chart.h"


Chart::Axis::Axis(Chart *p, bool isx)
	: parent(p), numTicks(-1), numBigTicks(-1), isX(isx)
{
}

void Chart::Axis::Draw(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2, float fontSize)
{
	Gdiplus::Pen *pen = new Gdiplus::Pen((Gdiplus::ARGB)Gdiplus::Color::Black);
	pen->SetWidth(2);

	Gdiplus::Point end(start.X + length + (isX ? length/30 : length/10), start.Y);
	g.DrawLine(pen, start, end);
	delete pen;

	end.X += 4;

	Gdiplus::Point end2(end.X - 15, end.Y - 6);
	Gdiplus::Point end3(end.X - 15, end.Y + 6);

	const Gdiplus::Point points[] = { end, end2, end3 };
	
	Gdiplus::SolidBrush *brush = new Gdiplus::SolidBrush((Gdiplus::ARGB)Gdiplus::Color::Black);
	g.FillPolygon(brush, points, 3);
	delete brush;

	DrawTicks(g, start, length);	

	if ((isX && parent->XAxisGrid) || (!isX && parent->YAxisGrid))
		DrawGrid(g, start, length, length2);

	DrawLabels(g, start, length, fontSize);
}

int Chart::Axis::GetNumTicks() const
{
	if (numTicks > 0) return numTicks;
	else if (nullptr == parent) return 10;

	return parent->GetNumTicks();
}

int Chart::Axis::GetNumBigTicks() const
{
	if (numBigTicks > 0) return numBigTicks;
	else if (nullptr == parent) return 5;

	return parent->GetNumBigTicks();
}

bool Chart::Axis::ShouldDrawFirstTick() const
{
	if (nullptr == parent) return false;

	return (isX ? parent->drawStartTickX : parent->drawStartTickY);
}

void Chart::Axis::DrawTicks(Gdiplus::Graphics& g, Gdiplus::Point& start, int length)
{
	int nrticks = GetNumTicks();
	double ticklen = static_cast<double>(length) / nrticks;

	Gdiplus::Pen pen((Gdiplus::ARGB)Gdiplus::Color::Black);
	pen.SetWidth(2);

	Gdiplus::Point pt(start);
	Gdiplus::Point pte(start);
	pte.Y += 6;
	for (int i = (ShouldDrawFirstTick() ? 0 : 1); i <= nrticks; ++i)
	{
		pt.X = static_cast<int>(start.X + i * ticklen);
		pte.X = pt.X;
		g.DrawLine(&pen, pt, pte);
	}

	nrticks = GetNumBigTicks();
	ticklen = static_cast<double>(length) / nrticks;

	pt = start;
	pte = start;
	pte.Y += 9;
	for (int i = (ShouldDrawFirstTick() ? 0 : 1); i <= nrticks; ++i)
	{
		pt.X = static_cast<int>(start.X + i * ticklen);
		pte.X = pt.X;
		g.DrawLine(&pen, pt, pte);
	}
}

std::list<CString> Chart::Axis::GetLabels() const
{
	if (labels.size()) return labels;
	else if (isX) return parent->GetXLabels();
	
	return parent->GetYLabels();
}

void Chart::Axis::DrawLabels(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, float fontSize)
{
	Gdiplus::GraphicsState state = g.Save();

	if (!isX) g.ScaleTransform(1, -1);

	std::list<CString> l = GetLabels();

	if (0 == l.size()) return;

	const int nrticks = GetNumBigTicks();
	double ticklen = static_cast<double>(length) / nrticks;

	Gdiplus::RectF rect;
	
	rect.Height = parent->GetLabelHeight(isX);
	rect.Width = parent->GetLabelWidth(isX);

	if (isX) {
		rect.X = -static_cast<float>(rect.Width/2.) - static_cast<float>(ShouldDrawFirstTick() ? ticklen : 0);
		rect.Y = static_cast<float>(start.Y + 12);
	}
	else {
		rect.X = static_cast<float>(start.Y - rect.Width - 12);
		rect.Y = -10 + static_cast<float>(ShouldDrawFirstTick() ? ticklen : 0);
	}
	
	auto label = l.begin();
	for (int i = (ShouldDrawFirstTick() ? 0 : 1); i <= nrticks; ++i)
	{
		if (label == l.end()) break;
				
		g.TranslateTransform(static_cast<float>(ticklen),0);
		
		if (!isX) g.RotateTransform(90.);

		// draw the label...
		parent->DrawText(*label, g, rect, isX ? Gdiplus::StringAlignment::StringAlignmentCenter : Gdiplus::StringAlignment::StringAlignmentFar, fontSize);

		if (!isX) g.RotateTransform(-90.);

		++label;
	}

	g.Restore(state);
}

void Chart::Axis::DrawGrid(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2)
{
	const int nrticks = GetNumBigTicks();
	const double ticklen = static_cast<double>(length) / nrticks;

	Gdiplus::Pen pen((Gdiplus::ARGB)Gdiplus::Color::LightGray);
	//pen.SetDashStyle(DashStyle::DashStyleDashDot);

	static const float pattern[] = { 6, 2, 3, 2 };
	pen.SetDashPattern(pattern, 4);

	Gdiplus::Point pt = start;
	pt.Y-=2;
	Gdiplus::Point pte = start;
	pte.Y -= length2;
	for (int i = 1; i < nrticks + 1; ++i)
	{
		pt.X = static_cast<int>(start.X + i * ticklen);
		pte.X = pt.X;
		g.DrawLine(&pen, pt, pte);
	}
}



double Chart::DataSets::DataSet::getXMin() const
{
	double result = DBL_MAX;

	for (auto &&p : points)
		if (p.X < result) result = p.X;

	return result;
}

double Chart::DataSets::DataSet::getYMin() const
{
	double result = DBL_MAX;

	for (const auto &p : points)
		if (p.Y < result) result = p.Y;

	return result;
}

double Chart::DataSets::DataSet::getXMax() const
{
	double result = DBL_MIN;

	for (const auto &p : points)
		if (p.X > result) result = p.X;

	return result;
}

double Chart::DataSets::DataSet::getYMax() const
{
	double result = DBL_MIN;

	for (const auto &p : points)
		if (p.Y > result) result = p.Y;

	return result;
}

double Chart::DataSets::DataSet::ConvertValue(double val, double valMin, double valMax, double chartMin, double chartMax)
{
	return chartMin + (val - valMin) / (valMax - valMin) * (chartMax - chartMin);
}


void Chart::DataSets::DataSet::Draw(Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, const Gdiplus::RectF& dataRect, bool spline) const
{
	if (points.size() < static_cast<unsigned int>(spline ? 3 : 2)) return;
	Gdiplus::Pen pen(Gdiplus::Color::MakeARGB(255, GetRValue(color), GetGValue(color), GetBValue(color)));
	pen.SetWidth(lineWidth);

	const double valMinX = dataRect.X;
	const double valMinY = dataRect.Y;
	const double valMaxX = dataRect.X + dataRect.Width;
	const double valMaxY = dataRect.Y + dataRect.Height;

	const double chartMinX = boundRect.X;
	const double chartMinY = boundRect.Y;
	const double chartMaxX = boundRect.X + boundRect.Width;
	const double chartMaxY = boundRect.Y + boundRect.Height;


	Gdiplus::PointF *pnts = new Gdiplus::PointF[points.size()];
	int index = 0;
	for (const auto &pnt : points)
	{
		pnts[index].X = static_cast<float>(ConvertValue(pnt.X, valMinX, valMaxX, chartMinX, chartMaxX));
		pnts[index].Y = static_cast<float>(ConvertValue(pnt.Y, valMinY, valMaxY, chartMinY, chartMaxY));
		++index;
	}

	if (spline) g.DrawCurve(&pen, pnts, index);
	else g.DrawLines(&pen, pnts, index);

	delete[] pnts;
}



double Chart::DataSets::getXMin() const
{
	double result = DBL_MAX;

	for (const auto &p : dataSets)
	{
		const double val = p.getXMin();
		if (val < result) result = val;
	}

	return result;
}

double Chart::DataSets::getYMin() const
{
	double result = DBL_MAX;

	for (auto &p : dataSets)
	{
		const double val = p.getYMin();
		if (val < result) result = val;
	}

	return result;
}

double Chart::DataSets::getXMax() const
{
	double result = DBL_MIN;

	for (const auto &p : dataSets)
	{
		const double val = p.getXMax();
		if (val > result) result = val;
	}

	return result;
}

double Chart::DataSets::getYMax() const
{
	double result = DBL_MIN;

	for (auto &p : dataSets)
	{
		const double val = p.getYMax();
		if (val > result) result = val;
	}

	return result;
}

void Chart::DataSets::Draw(Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, const Gdiplus::RectF& dataRect, bool spline) const
{
	g.SetClip(boundRect);
	for (const auto &dataSet : dataSets)
		dataSet.Draw(g, boundRect, dataRect, spline);
}



Chart::Chart()
	: X(this, true), Y(this, false), XAxisGrid(true), YAxisGrid(true),
	XAxisMin(DBL_MIN), XAxisMax(DBL_MAX), YAxisMin(DBL_MIN), YAxisMax(DBL_MAX),
	useSpline(true), antialias(false), drawStartTickX(true), drawStartTickY(true),
	maxTitleHeight(36), maxAxisLabelHeight(28), maxLabelHeight(18)
{
}


Chart::~Chart()
{
}


// will calculate from data if called
int Chart::GetNumTicks()
{
	return 10;
}

int Chart::GetNumBigTicks()
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


std::list<CString> Chart::GetXLabels()
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
		str.Format(L"%.2f", xmin + interval*i);
		l.push_back(str);
	}

	return l;
}


std::list<CString> Chart::GetYLabels()
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
		str.Format(L"%.2f", ymin + interval*i);
		l.push_back(str);
	}

	return l;
}


int Chart::GetNeededFontSize(CString& str, const Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect)
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

	// make room for title and labels
	int leftSide = min(max(rect.Width() / 20, 70), 100);
	int top = max(titleHeight * 2, rect.Height() / 10);
	rect.DeflateRect(leftSide, top, rect.Width() / 20, rect.Height() / 10);

	chartRect = rect;

	// try to find out the labels font size
	Gdiplus::RectF labelBound(0,0,GetLabelWidth(true),GetLabelHeight(true));

	float fontSize = static_cast<float>(maxLabelHeight);
	std::list<CString> labels = X.GetLabels();
	for (auto label : labels)
		fontSize = min(fontSize, static_cast<float>(GetNeededFontSize(label, g, labelBound)));
	
	labelBound.Width = GetLabelWidth(false);
	labelBound.Height = GetLabelHeight(false);
	labels = Y.GetLabels();
	for (auto label : labels)
		fontSize = min(fontSize, static_cast<float>(GetNeededFontSize(label, g, labelBound)));


	// draw horizontal axis
	g.TranslateTransform(static_cast<float>(rect.left), static_cast<float>(rect.bottom));
	Gdiplus::Point zero(0, 0);
	X.Draw(g, zero, rect.Width(), rect.Height(), fontSize);

	g.ScaleTransform(1,-1);
	g.RotateTransform(90.);

	Y.Draw(g, zero, rect.Height(), rect.Width(), fontSize);

	g.ResetTransform();

	// draw X label
	if (XAxisLabel.GetLength())
	{
		Gdiplus::RectF localLabelBound(static_cast<float>(rect.left), static_cast<float>(rect.top + rect.Height() + rect.Height() / 12.), static_cast<float>(rect.Width()), min(static_cast<float>(min(rect.Width(),rect.Height()) / 12.), maxAxisLabelHeight));
		DrawText(XAxisLabel, g, localLabelBound);
	}

	// draw Y label
	if (YAxisLabel.GetLength())
	{
		g.TranslateTransform(static_cast<float>(rect.left), static_cast<float>(rect.bottom));
		g.RotateTransform(-90.);
		Gdiplus::RectF localLabelBound(static_cast<float>(rect.Height() / 10.), -static_cast<float>(leftSide) , static_cast<float>(rect.Height()), min(static_cast<float>(min(rect.Width(), rect.Height()) / 12.), maxAxisLabelHeight));
		DrawText(YAxisLabel, g, localLabelBound);
	}

	g.ResetTransform();
	g.TranslateTransform(static_cast<float>(rect.left), static_cast<float>(rect.bottom));
	g.ScaleTransform(1, -1);

	boundRect.X = 0;
	boundRect.Y = 0;
	boundRect.Width = static_cast<float>(rect.Width());
	boundRect.Height = static_cast<float>(rect.Height());

	Gdiplus::RectF dataRect;

	dataRect.X = static_cast<float>((XAxisMin == DBL_MIN) ? dataSets.getXMin() : XAxisMin);
	dataRect.Y = static_cast<float>((YAxisMin == DBL_MIN) ? dataSets.getYMin() : YAxisMin);

	dataRect.Width = static_cast<float>(((XAxisMax == DBL_MAX) ? dataSets.getXMax() : XAxisMax) - dataRect.X);
	dataRect.Height = static_cast<float>(((YAxisMax == DBL_MAX) ? dataSets.getYMax() : YAxisMax) - dataRect.Y);

	if (dataRect.Width > 0 && dataRect.Height > 0) 
		dataSets.Draw(g, boundRect, dataRect, useSpline);

	// now draw the buffer on the screen
	Gdiplus::Graphics gdi(pDC->GetSafeHdc());
	gdi.DrawImage(&buffer, posX, posY);

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
	for (auto &&p : *data)
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
