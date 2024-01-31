#pragma once

#include <list>

class Chart;

class Axis {
public:
	Axis(Chart* parent, bool isX);

	int GetNumBigTicks() const;

	void Draw(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2, float fontSize = 0);
	void SetNumTicks(int ticks) { numTicks = ticks; }
	void SetNumBigTicks(int ticks) { numBigTicks = ticks; }

	void SetLabels(const std::list<CString>& l) { labels = l; }
	std::list<CString> GetLabels() const;

private:
	Chart* parent;
	bool isX;

	int numTicks;
	int numBigTicks;

	std::list<CString> labels;

	int GetNumTicks() const;
	bool ShouldDrawFirstTick() const;

	void DrawTicks(Gdiplus::Graphics& g, Gdiplus::Point& start, int length);
	void DrawTicks(Gdiplus::Graphics& g, Gdiplus::Pen& pen, Gdiplus::Point& start, int nrticks, int length, int offset);
	void DrawGrid(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2);
	void DrawLabels(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, float fontSize = 0);
	void DrawTip(Gdiplus::Graphics& g, const Gdiplus::Point& end);
};

