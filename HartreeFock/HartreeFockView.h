
// HartreeFockView.h : interface of the CHartreeFockView class
//

#pragma once


class CHartreeFockView : public CView
{
protected: // create from serialization only
	CHartreeFockView() = default;
	DECLARE_DYNCREATE(CHartreeFockView)

// Attributes
public:
	CHartreeFockDoc* GetDocument() const;

// Operations   
	void StartTimer();
	void StopTimer();


private:
// Implementation
	UINT_PTR timer = 0;
#ifdef _DEBUG
	void AssertValid() const override;
	void Dump(CDumpContext& dc) const override;
#endif

	// Overrides
	void OnDraw(CDC* pDC) override;  // overridden to draw this view
	BOOL PreCreateWindow(CREATESTRUCT& cs) override;
	void OnPrepareDC(CDC* pDC, CPrintInfo* pInfo = nullptr) override;
	BOOL OnPreparePrinting(CPrintInfo* pInfo) override;
	void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo) override;
	void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo) override;

	// Generated message map functions
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()

	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnDestroy();
	void StartWatching();
	afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
};

#ifndef _DEBUG  // debug version in HartreeFockView.cpp
inline CHartreeFockDoc* CHartreeFockView::GetDocument() const
   { return reinterpret_cast<CHartreeFockDoc*>(m_pDocument); }
#endif

