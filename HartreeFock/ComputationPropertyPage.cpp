// ComputationPropertyPage.cpp : implementation file
//

#include "stdafx.h"
#include "HartreeFock.h"
#include "ComputationPropertyPage.h"

#include"HartreeFock.h"

// ComputationPropertyPage

IMPLEMENT_DYNAMIC(ComputationPropertyPage, CMFCPropertyPage)

ComputationPropertyPage::ComputationPropertyPage()
	: CMFCPropertyPage(IDD_COMPUTATIONPROPERTYPAGE)
{
	m_nrThreads = theApp.options.nrThreads;
	m_useLotsOfMemory = (theApp.options.useLotsOfMemory ? BST_CHECKED : BST_UNCHECKED);
	m_nrPoints = theApp.options.numberOfPoints;
}

BEGIN_MESSAGE_MAP(ComputationPropertyPage, CMFCPropertyPage)
	ON_EN_CHANGE(IDC_EDIT1, &ComputationPropertyPage::OnEnChangeEdit1)
	ON_BN_CLICKED(IDC_CHECK1, &ComputationPropertyPage::OnBnClickedCheck1)
	ON_EN_CHANGE(IDC_EDIT3, &ComputationPropertyPage::OnEnChangeEdit3)
END_MESSAGE_MAP()



// ComputationPropertyPage message handlers




BOOL ComputationPropertyPage::OnApply()
{
	UpdateData();

	ApplyValues();

	return CMFCPropertyPage::OnApply();
}


void ComputationPropertyPage::ApplyValues()
{
	theApp.options.nrThreads = m_nrThreads;
	theApp.options.useLotsOfMemory = m_useLotsOfMemory == BST_CHECKED;
	theApp.options.numberOfPoints = m_nrPoints;

	theApp.options.Save();
}


BOOL ComputationPropertyPage::OnInitDialog()
{
	CMFCPropertyPage::OnInitDialog();

	// TODO:  Add extra initialization here

	return TRUE;  // return TRUE unless you set the focus to a control
				  // EXCEPTION: OCX Property Pages should return FALSE
}


void ComputationPropertyPage::DoDataExchange(CDataExchange* pDX)
{
	CMFCPropertyPage::DoDataExchange(pDX);

	DDX_Text(pDX, IDC_EDIT1, m_nrThreads);
	DDX_Check(pDX, IDC_CHECK1, m_useLotsOfMemory);
	DDX_Text(pDX, IDC_EDIT3, m_nrPoints);

	DDV_MinMaxUInt(pDX, m_nrThreads, 1, 256);
	DDV_MinMaxUInt(pDX, m_nrPoints, theApp.options.nrThreads > 0 ? theApp.options.nrThreads * 2 : 2, 1000);
}


void ComputationPropertyPage::OnEnChangeEdit1()
{
	SetModified();
}


void ComputationPropertyPage::OnBnClickedCheck1()
{
	SetModified();
}


void ComputationPropertyPage::OnEnChangeEdit3()
{
	SetModified();
}
