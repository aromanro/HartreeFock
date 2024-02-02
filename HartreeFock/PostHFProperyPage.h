#pragma once


// PostHFProperyPage

class PostHFProperyPage : public CMFCPropertyPage
{
	DECLARE_DYNAMIC(PostHFProperyPage)

public:
	PostHFProperyPage();

	// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_POSTHFPROPERTYPAGE };
#endif
private:
	DECLARE_MESSAGE_MAP()

	BOOL OnApply() override;
	void ApplyValues();
	BOOL OnInitDialog() override;
	void DoDataExchange(CDataExchange* pDX) override;
	afx_msg void OnBnClickedCheck1();

	int computePostHF;
};


