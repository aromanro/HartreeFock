
// HartreeFockView.h : interface of the CHartreeFockView class
//

#pragma once


class CHartreeFockView : public CView
{
protected: // create from serialization only
	CHartreeFockView();
	DECLARE_DYNCREATE(CHartreeFockView)

// Attributes
public:
	CHartreeFockDoc* GetDocument() const;

// Operations
public:

// Overrides
public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// Implementation
protected:
	UINT_PTR timer;
public:
	virtual ~CHartreeFockView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnDestroy();
	void StartWatching();
	afx_msg BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);

	void StartTimer();
	void StopTimer();
	virtual void OnPrepareDC(CDC* pDC, CPrintInfo* pInfo = NULL);
};

#ifndef _DEBUG  // debug version in HartreeFockView.cpp
inline CHartreeFockDoc* CHartreeFockView::GetDocument() const
   { return reinterpret_cast<CHartreeFockDoc*>(m_pDocument); }
#endif

