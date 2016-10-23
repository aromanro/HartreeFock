
// HartreeFockView.cpp : implementation of the CHartreeFockView class
//

#include "stdafx.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "HartreeFock.h"
#endif

#include "HartreeFockDoc.h"
#include "HartreeFockView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CHartreeFockView

IMPLEMENT_DYNCREATE(CHartreeFockView, CView)

BEGIN_MESSAGE_MAP(CHartreeFockView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CHartreeFockView::OnFilePrintPreview)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_TIMER()
	ON_WM_DESTROY()
	ON_WM_SETCURSOR()
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()

// CHartreeFockView construction/destruction

CHartreeFockView::CHartreeFockView()
	: timer(NULL)
{
	// TODO: add construction code here

}

CHartreeFockView::~CHartreeFockView()
{
}

BOOL CHartreeFockView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CHartreeFockView drawing

void CHartreeFockView::OnDraw(CDC* pDC)
{
	CHartreeFockDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	CRect rect;
	GetClientRect(rect);

	pDoc->m_Chart.Draw(pDC, rect);
}


// CHartreeFockView printing


void CHartreeFockView::OnFilePrintPreview()
{
#ifndef SHARED_HANDLERS
	AFXPrintPreview(this);
#endif
}

BOOL CHartreeFockView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CHartreeFockView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CHartreeFockView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CHartreeFockView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CHartreeFockView::OnContextMenu(CWnd* /* pWnd */, CPoint /*point*/)
{
#ifndef SHARED_HANDLERS
#endif
}


// CHartreeFockView diagnostics

#ifdef _DEBUG
void CHartreeFockView::AssertValid() const
{
	CView::AssertValid();
}

void CHartreeFockView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CHartreeFockDoc* CHartreeFockView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CHartreeFockDoc)));
	return (CHartreeFockDoc*)m_pDocument;
}
#endif //_DEBUG


// CHartreeFockView message handlers


void CHartreeFockView::OnTimer(UINT_PTR nIDEvent)
{
	// TODO: Add your message handler code here and/or call default

	CView::OnTimer(nIDEvent);

	CHartreeFockDoc* pDoc = GetDocument();

	if (pDoc->isFinished())
	{
		KillTimer(timer);
		timer = 0;

		EndWaitCursor();
		Invalidate();
	}
}


void CHartreeFockView::OnDestroy()
{
	CView::OnDestroy();

	StopTimer();
}


void CHartreeFockView::StartWatching()
{
	if (!timer) timer = SetTimer(1, 1000, NULL);
	BeginWaitCursor();
}


BOOL CHartreeFockView::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message)
{
	CHartreeFockDoc* pDoc = GetDocument();

	if (pDoc && !pDoc->isFinished())
	{
		RestoreWaitCursor();

		return TRUE;
	}

	return CView::OnSetCursor(pWnd, nHitTest, message);
}


BOOL CHartreeFockView::OnEraseBkgnd(CDC* /*pDC*/)
{
	return TRUE;
}


void CHartreeFockView::StartTimer()
{
	if (!timer) timer = SetTimer(1, 100, NULL);

	CHartreeFockDoc* pDoc = GetDocument();
	if (pDoc && !pDoc->isFinished()) BeginWaitCursor();
}


void CHartreeFockView::StopTimer()
{
	if (timer) {
		KillTimer(timer);
		timer = 0;

		EndWaitCursor();
	}
}

void CHartreeFockView::OnPrepareDC(CDC* pDC, CPrintInfo* pInfo)
{
	CView::OnPrepareDC(pDC, pInfo);

	if (pDC->IsPrinting())
	{
		CRect rect;
		GetClientRect(rect);

		pDC->SetMapMode(MM_ISOTROPIC);

		int cx = pDC->GetDeviceCaps(HORZRES);
		int cy = pDC->GetDeviceCaps(VERTRES);

		pDC->SetWindowExt(rect.Width(), rect.Height());
		pDC->SetViewportExt(cx, cy);
		pDC->SetViewportOrg(0, 0);
	}
}
