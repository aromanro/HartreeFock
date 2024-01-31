#pragma once


// CNumberEdit

class CNumberEdit : public CEdit
{
	DECLARE_DYNAMIC(CNumberEdit)

	CString m_GoodValue;
	DWORD m_GoodSel = 0;
	bool m_Cancelling = false; 
	bool allowNegative = true;

	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnEnUpdate();
	afx_msg	LRESULT OnPaste(WPARAM wParam, LPARAM lParam);

	DECLARE_MESSAGE_MAP()
};


