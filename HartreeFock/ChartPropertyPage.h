#pragma once
#include "afxwin.h"

#include "NumberEdit.h"

// ChartPropertyPage

class ChartPropertyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(ChartPropertyPage)

public:
	ChartPropertyPage();
	virtual ~ChartPropertyPage();

#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_CHARTPROPERTYPAGE };
#endif

protected:
	DECLARE_MESSAGE_MAP()

	virtual BOOL OnApply();
	void ApplyValues();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);

	afx_msg void OnEnChangeEdit1();
	afx_msg void OnEnChangeEdit2();
	afx_msg void OnEnChangeEdit3();
	afx_msg void OnEnChangeEdit4();
	afx_msg void OnEnChangeEdit5();
	afx_msg void OnEnChangeEdit6();	
	afx_msg void OnBnClickedCheck2();

	CNumberEdit m_minEnergy;
	CNumberEdit m_maxEnergy;

	int YMaxEnergy; // eV
	int YMinEnergy; // eV
	unsigned int YBigTicksEnergy;
	unsigned int YSmallTicksEnergy;

	unsigned int XBigTicksBondLength;
	unsigned int XSmallTicksBondLength;
	int useSplines;

	int DisplayHOMOEnergy;
public:
	afx_msg void OnBnClickedRadio1();
	afx_msg void OnBnClickedRadio2();
	afx_msg void OnBnClickedRadio3();
};


