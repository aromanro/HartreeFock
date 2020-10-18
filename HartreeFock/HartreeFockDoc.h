
// HartreeFockDoc.h : interface of the CHartreeFockDoc class
//


#pragma once

#include "Basis.h"

#include "Chart.h"



#include <atomic>
#include <list>
#include <memory>

#include "HartreeFockThread.h"

class CHartreeFockView;

class CHartreeFockDoc : public CDocument
{
protected: // create from serialization only
	CHartreeFockDoc();
	DECLARE_DYNCREATE(CHartreeFockDoc)

// Attributes
public:

	Chemistry::Basis basisSTO3G;
	Chemistry::Basis basisSTO6G;

	Chemistry::Basis basis3_21G;
	Chemistry::Basis basis6_21G;
	Chemistry::Basis basis6_31G;

	Chemistry::Basis basis6_31Gstar;
	Chemistry::Basis basis6_31plusGstarstar;

	Chemistry::Basis basis6_31plusG;
	Chemistry::Basis basis6_31plusGstar;
	Chemistry::Basis basis6_31plusplusG;
	Chemistry::Basis basis6_31plusplusGstar;
	Chemistry::Basis basis6_31plusplusGstarstar;

	Chemistry::Basis basis6_311G;
	Chemistry::Basis basis6_311Gstar;
	Chemistry::Basis basis6_311Gstarstar;
	Chemistry::Basis basis6_311plusG;
	Chemistry::Basis basis6_311plusGstar;
	Chemistry::Basis basis6_311plusGstarstar;
	Chemistry::Basis basis6_311plusplusG;
	Chemistry::Basis basis6_311plusplusGstar;
	Chemistry::Basis basis6_311plusplusGstarstar;

	Chart m_Chart;

	std::atomic_int runningThreads;

	Options options; // to save options into, will let the options from the app to be changed while running threads

	std::list<std::unique_ptr<HartreeFockThread>> threadsList;

	std::vector<std::tuple<double, double, double>> results;
	bool convergenceProblem;

	double atomsEnergy;

// Operations
public:

// Overrides
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// Implementation
public:
	virtual ~CHartreeFockDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// Helper function that sets search content for a Search Handler
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
public:
	afx_msg void OnComputationStart();
	bool isFinished();
	void StartThreads();
	void StopThreads(bool cancel = false);
	CHartreeFockView* GetView();
	afx_msg void OnUpdateComputationStart(CCmdUI *pCmdUI);
	void SetChartBoundsAndTicks();
	void ApplyChartOptions();
	void SetChartData();
};
