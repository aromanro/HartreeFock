#pragma once


// ComputationPropertyPage

class ComputationPropertyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(ComputationPropertyPage)

public:
	ComputationPropertyPage();
	~ComputationPropertyPage() override;


	// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_COMPUTATIONPROPERTYPAGE };
#endif

private:
	DECLARE_MESSAGE_MAP()

	BOOL OnApply() override;
	void ApplyValues();
	BOOL OnInitDialog() override;
	void DoDataExchange(CDataExchange* pDX) override;
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedCheck1();
	afx_msg void OnEnChangeEdit3();

	int m_nrThreads;
	int m_useLotsOfMemory;
	int m_nrPoints;
};


