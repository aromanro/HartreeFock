#pragma once


// CNumberEdit

class CNumberEdit : public CEdit
{
	DECLARE_DYNAMIC(CNumberEdit)

public:
	CNumberEdit();
	virtual ~CNumberEdit();

protected:
	DECLARE_MESSAGE_MAP()
private:
	CString m_GoodValue;
	DWORD m_GoodSel;
	bool m_Cancelling;

	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnEnUpdate();
	afx_msg	LRESULT OnPaste(WPARAM wParam, LPARAM lParam);

public:
	bool allowNegative;
};


