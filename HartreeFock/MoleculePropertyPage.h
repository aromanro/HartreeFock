#pragma once
#include "afxwin.h"

#include "NumberEdit.h"

// MoleculePropertyPage

class MoleculePropertyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(MoleculePropertyPage)

public:
	MoleculePropertyPage();
	~MoleculePropertyPage() override;


#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MOLECULEPROPERTYPAGE };
#endif

private:
	DECLARE_MESSAGE_MAP()

	virtual BOOL OnApply();
	void ApplyValues();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);


	void FillCombos();
	void AdjustElectrons();

	afx_msg void OnCbnSelchangeCombo1();
	afx_msg void OnCbnSelchangeCombo2();
	afx_msg void OnBnClickedCheck1();
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedRadio1();
	afx_msg void OnBnClickedRadio2();
	afx_msg void OnBnClickedRadio4();
	afx_msg void OnBnClickedRadio5();
	afx_msg void OnBnClickedRadio6();
	afx_msg void OnBnClickedRadio7();
	afx_msg void OnBnClickedRadio8();
	afx_msg void OnEnChangeEdit5();
	afx_msg void OnEnChangeEdit6();
	afx_msg void OnEnChangeEdit7();
	afx_msg void OnEnChangeEdit8();
	afx_msg void OnBnClickedRadio9();
	afx_msg void OnBnClickedRadio10();
	afx_msg void OnBnClickedRadio11();
	afx_msg void OnBnClickedRadio12();
	afx_msg void OnBnClickedRadio13();
	afx_msg void OnBnClickedRadio14();
	afx_msg void OnBnClickedRadio15();
	afx_msg void OnBnClickedRadio16();
	afx_msg void OnBnClickedRadio17();
	afx_msg void OnBnClickedRadio18();
	afx_msg void OnBnClickedRadio19();
	afx_msg void OnBnClickedRadio20();
	afx_msg void OnBnClickedRadio21();
	afx_msg void OnBnClickedRadio22();
	afx_msg void OnBnClickedRadio23();
	afx_msg void OnBnClickedRadio24();

	CString m_atom1;
	CString m_atom2;
	int twoAtom1;

	int basis; // 0 - STO3G, 1 - STO6G, 2 - 3-21G, 3 - 6-21G, 4 - 6-31G, 5 - 6-31G*, 6 - 6-31+G**
	double bondAngle;

	int alphaElectrons;
	int betaElectrons;

	// also used for charting
	double XMaxBondLength;
	double XMinBondLength; // not really for chart, but for calculations

	CNumberEdit bondAngleEdit;
	CNumberEdit minBondEdit;
	CNumberEdit maxBondEdit;
	CComboBox comboBox1;
	CComboBox comboBox2;
};


