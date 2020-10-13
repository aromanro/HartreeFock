// PostHFProperyPage.cpp : implementation file
//

#include "stdafx.h"
#include "HartreeFock.h"
#include "PostHFProperyPage.h"


// PostHFProperyPage

IMPLEMENT_DYNAMIC(PostHFProperyPage, CMFCPropertyPage)

PostHFProperyPage::PostHFProperyPage()
	: CMFCPropertyPage(IDD_POSTHFPROPERTYPAGE)
{
	computePostHF = (theApp.options.computePostHF ? BST_CHECKED : BST_UNCHECKED);
}

PostHFProperyPage::~PostHFProperyPage()
{
}


BEGIN_MESSAGE_MAP(PostHFProperyPage, CMFCPropertyPage)
	ON_BN_CLICKED(IDC_CHECK1, &PostHFProperyPage::OnBnClickedCheck1)
END_MESSAGE_MAP()



// PostHFProperyPage message handlers




BOOL PostHFProperyPage::OnApply()
{
	UpdateData();

	ApplyValues();

	return CMFCPropertyPage::OnApply();
}


void PostHFProperyPage::ApplyValues()
{
	theApp.options.computePostHF = (computePostHF == BST_CHECKED ? true : false);

	theApp.options.Save();
}


BOOL PostHFProperyPage::OnInitDialog()
{
	CMFCPropertyPage::OnInitDialog();

	// TODO:  Add extra initialization here

	return TRUE;  // return TRUE unless you set the focus to a control
				  // EXCEPTION: OCX Property Pages should return FALSE
}


void PostHFProperyPage::DoDataExchange(CDataExchange* pDX)
{
	CMFCPropertyPage::DoDataExchange(pDX);

	DDX_Check(pDX, IDC_CHECK1, computePostHF);
}

void PostHFProperyPage::OnBnClickedCheck1()
{
	SetModified();
}
