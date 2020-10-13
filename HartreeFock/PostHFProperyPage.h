#pragma once


// PostHFProperyPage

class PostHFProperyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(PostHFProperyPage)

public:
	PostHFProperyPage();
	virtual ~PostHFProperyPage();

	// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_POSTHFPROPERTYPAGE };
#endif
protected:
	DECLARE_MESSAGE_MAP()

	virtual BOOL OnApply();
	void ApplyValues();
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);
	afx_msg void OnBnClickedCheck1();

	int computePostHF;
};


