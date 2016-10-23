#pragma once


// ComputationPropertyPage

class ComputationPropertyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(ComputationPropertyPage)

public:
	ComputationPropertyPage();
	virtual ~ComputationPropertyPage();


	// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_COMPUTATIONPROPERTYPAGE };
#endif

protected:
	DECLARE_MESSAGE_MAP()

	virtual BOOL OnApply();
	void ApplyValues();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedCheck1();
	afx_msg void OnEnChangeEdit3();

	int m_nrThreads;
	int m_useLotsOfMemory;
	int m_nrPoints;
};


