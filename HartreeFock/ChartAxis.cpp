#include "stdafx.h"
#include "ChartAxis.h"
#include "Chart.h"

Axis::Axis(Chart* p, bool isx)
	: parent(p), numTicks(-1), numBigTicks(-1), isX(isx)
{
}

void Axis::Draw(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2, float fontSize)
{
	Gdiplus::Pen* pen = new Gdiplus::Pen((Gdiplus::ARGB)Gdiplus::Color::Black);
	pen->SetWidth(2);

	Gdiplus::Point end(start.X + length + (isX ? length / 30 : length / 10), start.Y);
	g.DrawLine(pen, start, end);
	delete pen;

	end.X += 4;

	DrawTip(g, end);
	DrawTicks(g, start, length);

	if ((isX && parent->XAxisGrid) || (!isX && parent->YAxisGrid))
		DrawGrid(g, start, length, length2);

	DrawLabels(g, start, length, fontSize);
}

void Axis::DrawTip(Gdiplus::Graphics& g, const Gdiplus::Point& end)
{
	Gdiplus::Point end2(end.X - 15, end.Y - 6);
	Gdiplus::Point end3(end.X - 15, end.Y + 6);

	const Gdiplus::Point points[] = { end, end2, end3 };

	Gdiplus::SolidBrush* brush = new Gdiplus::SolidBrush((Gdiplus::ARGB)Gdiplus::Color::Black);
	g.FillPolygon(brush, points, 3);
	delete brush;
}


int Axis::GetNumTicks() const
{
	if (numTicks > 0) return numTicks;
	else if (nullptr == parent) return 10;

	return parent->GetNumTicks();
}

int Axis::GetNumBigTicks() const
{
	if (numBigTicks > 0) return numBigTicks;
	else if (nullptr == parent) return 5;

	return parent->GetNumBigTicks();
}

bool Axis::ShouldDrawFirstTick() const
{
	if (nullptr == parent) return false;

	return (isX ? parent->drawStartTickX : parent->drawStartTickY);
}

void Axis::DrawTicks(Gdiplus::Graphics& g, Gdiplus::Point& start, int length)
{
	Gdiplus::Pen pen((Gdiplus::ARGB)Gdiplus::Color::Black);
	pen.SetWidth(2);

	DrawTicks(g, pen, start, GetNumTicks(), length, 6);
	DrawTicks(g, pen, start, GetNumBigTicks(), length, 9);
}

void Axis::DrawTicks(Gdiplus::Graphics& g, Gdiplus::Pen& pen, Gdiplus::Point& start, int nrticks, int length, int offset)
{
	double ticklen = static_cast<double>(length) / nrticks;

	Gdiplus::Point pt(start);
	Gdiplus::Point pte(start);
	pte.Y += offset;
	for (int i = (ShouldDrawFirstTick() ? 0 : 1); i <= nrticks; ++i)
	{
		pt.X = static_cast<int>(start.X + i * ticklen);
		pte.X = pt.X;
		g.DrawLine(&pen, pt, pte);
	}
}


std::list<CString> Axis::GetLabels() const
{
	if (labels.size()) return labels;
	else if (isX) return parent->GetXLabels();

	return parent->GetYLabels();
}

void Axis::DrawLabels(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, float fontSize)
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
		rect.X = -static_cast<float>(rect.Width / 2.) - static_cast<float>(ShouldDrawFirstTick() ? ticklen : 0);
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

		g.TranslateTransform(static_cast<float>(ticklen), 0);

		if (!isX) g.RotateTransform(90.);

		// draw the label...
		parent->DrawText(*label, g, rect, isX ? Gdiplus::StringAlignment::StringAlignmentCenter : Gdiplus::StringAlignment::StringAlignmentFar, fontSize);

		if (!isX) g.RotateTransform(-90.);

		++label;
	}

	g.Restore(state);
}

void Axis::DrawGrid(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2)
{
	const int nrticks = GetNumBigTicks();
	const double ticklen = static_cast<double>(length) / nrticks;

	Gdiplus::Pen pen((Gdiplus::ARGB)Gdiplus::Color::LightGray);
	//pen.SetDashStyle(DashStyle::DashStyleDashDot);

	static const float pattern[] = { 6, 2, 3, 2 };
	pen.SetDashPattern(pattern, 4);

	Gdiplus::Point pt = start;
	pt.Y -= 2;
	Gdiplus::Point pte = start;
	pte.Y -= length2;
	for (int i = 1; i < nrticks + 1; ++i)
	{
		pt.X = static_cast<int>(start.X + i * ticklen);
		pte.X = pt.X;
		g.DrawLine(&pen, pt, pte);
	}
}
