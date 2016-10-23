#pragma once
#include "afxwin.h"

#include "NumberEdit.h"

// HartreeFockPropertyPage

class HartreeFockPropertyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(HartreeFockPropertyPage)

public:
	HartreeFockPropertyPage();
	virtual ~HartreeFockPropertyPage();

#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_HARTREEFOCKPROPERTYPAGE };
#endif

protected:
	DECLARE_MESSAGE_MAP()

	virtual BOOL OnApply();
	void ApplyValues();
	virtual BOOL OnInitDialog();
	afx_msg void OnBnClickedRadio1();
	afx_msg void OnBnClickedRadio2();
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnEnChangeEdit2();
	afx_msg void OnEnChangeEdit3();

	virtual void DoDataExchange(CDataExchange* pDX);

	int m_Method;
	CNumberEdit m_AlphaEdit;
	CNumberEdit m_GuessEdit;
	double m_Alpha;
	double m_Guess;
	int iterations;
};


