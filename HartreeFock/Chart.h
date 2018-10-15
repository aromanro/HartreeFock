#pragma once

#include <list>
#include <vector>

class Chart
{
protected:
	class Axis {
	protected:
		Chart* parent;
		bool isX;

		int numTicks;
		int numBigTicks;

		std::list<CString> labels;

		int GetNumTicks() const;
		bool ShouldDrawFirstTick() const;

		void DrawTicks(Gdiplus::Graphics& g, Gdiplus::Point& start, int length);
		void DrawGrid(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2);
		void DrawLabels(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, float fontSize = 0);
	public:
		Axis(Chart* parent, bool isX);

		int GetNumBigTicks() const;

		void Draw(Gdiplus::Graphics& g, Gdiplus::Point& start, int length, int length2, float fontSize = 0);
		void SetNumTicks(int ticks) { numTicks = ticks; }
		void SetNumBigTicks(int ticks) { numBigTicks = ticks; }

		void SetLabels(const std::list<CString>& l) { labels = l; }
		std::list<CString> GetLabels() const;
	};

	Axis X;
	Axis Y;

	CRect chartRect;

	int GetNumTicks() const;
	int GetNumBigTicks() const;

	std::list<CString> GetXLabels() const;
	std::list<CString> GetYLabels() const;

	float GetLabelHeight(bool XAxis = true) const;
	float GetLabelWidth(bool XAxis = true) const;

	int GetNeededFontSize(CString& str, const Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect) const;

	void DrawText(CString &str, Gdiplus::Graphics& g, const Gdiplus::RectF& boundRect, Gdiplus::StringAlignment align = Gdiplus::StringAlignmentCenter, float fontSize = 0);
public:
	Chart();
	~Chart();

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

	DataSets dataSets;

	void Draw(const CDC* pDC, CRect& rect);

	void SetNumTicksX(int ticks) { X.SetNumTicks(ticks); }
	void SetNumTicksY(int ticks) { Y.SetNumTicks(ticks); }
	void SetNumBigTicksX(int ticks) { X.SetNumBigTicks(ticks); }
	void SetNumBigTicksY(int ticks) { Y.SetNumBigTicks(ticks); }
	void SetXLabels(const std::list<CString>& l) { X.SetLabels(l); }
	void SetYLabels(const std::list<CString>& l) { Y.SetLabels(l); }

	void AddDataSet(const double *dataX, const double *dataY, unsigned int len, float lineWidth = 0, COLORREF color = RGB(0, 0, 0));
	void AddDataSet(const std::vector<std::pair<double, double>> *data, float lineWidth = 0, COLORREF color = RGB(0, 0, 0));
	void AddDataSet(double start, double step, const double *dataY, unsigned int len, float lineWidth = 0, COLORREF color = RGB(0, 0, 0));
	void AddDataSlice(double X, const double *slice, unsigned int len);
	void clear();


	CString title;
	CString XAxisLabel;
	CString YAxisLabel;
	bool XAxisGrid;
	bool YAxisGrid;

	bool useSpline;
	bool antialias;

	double XAxisMin;
	double XAxisMax;
	double YAxisMin;
	double YAxisMax;

	bool drawStartTickX;
	bool drawStartTickY;

	int maxTitleHeight;
	int maxAxisLabelHeight;
	int maxLabelHeight;
};

