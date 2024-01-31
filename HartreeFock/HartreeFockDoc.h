
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

	Chemistry::Basis dz;
	Chemistry::Basis dzp;

	Chart m_Chart;

	std::atomic_int runningThreads;

	Options options; // to save options into, will let the options from the app to be changed while running threads

	std::list<std::unique_ptr<HartreeFockThread>> threadsList;

	std::vector<std::tuple<double, double, double>> results;
	bool convergenceProblem;

	double atomsEnergy;

// Operations
// Overrides
	BOOL OnNewDocument() override;
	void Serialize(CArchive& ar) override;
#ifdef SHARED_HANDLERS
	void InitializeSearchContent() override;
	void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds) override;
#endif // SHARED_HANDLERS

// Implementation
	~CHartreeFockDoc() override;
#ifdef _DEBUG
	void AssertValid() const override;
	void Dump(CDumpContext& dc) const override;
#endif

protected:
// Generated message map functions
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
	Chemistry::Basis* GetBasis(const Options& options);

protected:
	void SetChartData();
	void SetChartTitle();
	Chemistry::Basis* GetBasis6_31(const Options& options);
	Chemistry::Basis* GetBasis6_311(const Options& options);
};
