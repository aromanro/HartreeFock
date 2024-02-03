#pragma once

#include "ChartPropertyPage.h"
#include "HartreeFockPropertyPage.h"
#include "ComputationPropertyPage.h"
#include "MoleculePropertyPage.h"
#include "PostHFProperyPage.h"

// OptionsPropertySheet

class COptionsPropertySheet : public CMFCPropertySheet
{
	DECLARE_DYNAMIC(COptionsPropertySheet)

public:
	COptionsPropertySheet(UINT nIDCaption, CWnd* pParentWnd = nullptr, UINT iSelectPage = 0);
	COptionsPropertySheet(LPCTSTR pszCaption, CWnd* pParentWnd = nullptr, UINT iSelectPage = 0);
	~COptionsPropertySheet() override;

private:
	HICON hIcon;

	MoleculePropertyPage page1;
	HartreeFockPropertyPage page2;
	ComputationPropertyPage page3;
	PostHFProperyPage page4;
	ChartPropertyPage page5;

	DECLARE_MESSAGE_MAP()

	BOOL OnInitDialog() override;
	void AddPages();
};


