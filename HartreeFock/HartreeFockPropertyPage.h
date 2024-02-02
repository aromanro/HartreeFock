#pragma once
#include "afxwin.h"

#include "NumberEdit.h"

// HartreeFockPropertyPage

class HartreeFockPropertyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(HartreeFockPropertyPage)

public:
	HartreeFockPropertyPage();

#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_HARTREEFOCKPROPERTYPAGE };
#endif

private:
	DECLARE_MESSAGE_MAP()

	BOOL OnApply() override;
	void ApplyValues();
	BOOL OnInitDialog() override;
	afx_msg void OnBnClickedRadio1();
	afx_msg void OnBnClickedRadio2();
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnEnChangeEdit2();
	afx_msg void OnEnChangeEdit3();
	afx_msg void OnBnClickedCheck2();
	afx_msg void OnEnChangeEdit4();
	afx_msg void OnBnClickedCheck5();
	afx_msg void OnEnChangeEdit5();
	afx_msg void OnEnChangeEdit6();

	void DoDataExchange(CDataExchange* pDX) override;

	CNumberEdit m_AlphaEdit;
	CNumberEdit m_GuessEdit;
	CNumberEdit m_AsymmetryEdit;

	int m_Method;
	double m_Alpha;
	double m_Guess;
	int iterations;
	double asymmetry;
	int addAsymmetry;
	int useDIIS;
	int maxDIISiterations;
	int normalIterAfterDIIS;
};


