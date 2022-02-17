// NumberEdit.cpp : implementation file
//

#include "stdafx.h"

#include "NumberEdit.h"


// CNumberEdit

IMPLEMENT_DYNAMIC(CNumberEdit, CEdit)

CNumberEdit::CNumberEdit()
	: m_GoodSel(0), m_Cancelling(false), allowNegative(true)
{
}

CNumberEdit::~CNumberEdit()
{
}


BEGIN_MESSAGE_MAP(CNumberEdit, CEdit)
	ON_WM_CHAR()
	ON_CONTROL_REFLECT(EN_UPDATE, &CNumberEdit::OnEnUpdate)
	ON_MESSAGE(WM_PASTE, OnPaste)
END_MESSAGE_MAP()



// CNumberEdit message handlers


void CNumberEdit::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	GetWindowText(m_GoodValue);
	m_GoodSel = GetSel();

	CEdit::OnChar(nChar, nRepCnt, nFlags);
}


void CNumberEdit::OnEnUpdate()
{
	if (m_Cancelling) return;

	CString str;
	GetWindowText(str);
	if (allowNegative && str == L"-")
	{
		m_GoodValue = str;
		return;
	}

	wchar_t *end;
	const double val = wcstod(str, &end);
	if (*end || (val < 0 && !allowNegative))
	{
		m_Cancelling = true;
		SetWindowText(m_GoodValue);
		SetSel(m_GoodSel);
		m_Cancelling = false;
	}
	else
		m_GoodValue = str;
}


LRESULT CNumberEdit::OnPaste(WPARAM /*wParam*/, LPARAM /*lParam*/)
{
	// do nothing for now
	return 0;
}
