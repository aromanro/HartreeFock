// ChartPropertyPage.cpp : implementation file
//

#include "stdafx.h"
#include "HartreeFock.h"
#include "ChartPropertyPage.h"

#include "MainFrm.h"

// ChartPropertyPage

IMPLEMENT_DYNAMIC(ChartPropertyPage, CMFCPropertyPage)

ChartPropertyPage::ChartPropertyPage()
	: CMFCPropertyPage(IDD_CHARTPROPERTYPAGE)
{
	YMaxEnergy = theApp.options.YMaxEnergy; // eV
	YMinEnergy = theApp.options.YMinEnergy; // eV
	YBigTicksEnergy = theApp.options.YBigTicksEnergy;
	YSmallTicksEnergy = theApp.options.YSmallTicksEnergy;

	XBigTicksBondLength = theApp.options.XBigTicksBondLength;
	XSmallTicksBondLength = theApp.options.XSmallTicksBondLength;

	useSplines = (theApp.options.useSplines ? BST_CHECKED : BST_UNCHECKED);

	DisplayHOMOEnergy = theApp.options.DisplayHOMOEnergy;
}

ChartPropertyPage::~ChartPropertyPage()
{
}


BEGIN_MESSAGE_MAP(ChartPropertyPage, CMFCPropertyPage)
	ON_EN_CHANGE(IDC_EDIT1, &ChartPropertyPage::OnEnChangeEdit1)
	ON_EN_CHANGE(IDC_EDIT2, &ChartPropertyPage::OnEnChangeEdit2)
	ON_EN_CHANGE(IDC_EDIT3, &ChartPropertyPage::OnEnChangeEdit3)
	ON_EN_CHANGE(IDC_EDIT4, &ChartPropertyPage::OnEnChangeEdit4)
	ON_EN_CHANGE(IDC_EDIT5, &ChartPropertyPage::OnEnChangeEdit5)
	ON_EN_CHANGE(IDC_EDIT6, &ChartPropertyPage::OnEnChangeEdit6)
	ON_BN_CLICKED(IDC_CHECK2, &ChartPropertyPage::OnBnClickedCheck2)
	ON_BN_CLICKED(IDC_RADIO1, &ChartPropertyPage::OnBnClickedRadio1)
	ON_BN_CLICKED(IDC_RADIO2, &ChartPropertyPage::OnBnClickedRadio2)
END_MESSAGE_MAP()



// ChartPropertyPage message handlers




BOOL ChartPropertyPage::OnApply()
{
	UpdateData();

	ApplyValues();

	return CMFCPropertyPage::OnApply();
}


void ChartPropertyPage::ApplyValues()
{
	theApp.options.YMaxEnergy = YMaxEnergy; // eV
	theApp.options.YMinEnergy = YMinEnergy; // eV
	theApp.options.YBigTicksEnergy = YBigTicksEnergy;
	theApp.options.YSmallTicksEnergy = YSmallTicksEnergy;

	theApp.options.XBigTicksBondLength = XBigTicksBondLength;
	theApp.options.XSmallTicksBondLength = XSmallTicksBondLength;

	theApp.options.useSplines = (useSplines == BST_CHECKED ? true : false);

	theApp.options.DisplayHOMOEnergy = DisplayHOMOEnergy;

	theApp.options.Save();

	dynamic_cast<CMainFrame*>(theApp.m_pMainWnd)->GetDocument()->ApplyChartOptions();
}


BOOL ChartPropertyPage::OnInitDialog()
{
	CMFCPropertyPage::OnInitDialog();

	m_minEnergy.allowNegative = true;
	m_maxEnergy.allowNegative = true;

	return TRUE;  // return TRUE unless you set the focus to a control
				  // EXCEPTION: OCX Property Pages should return FALSE
}


void ChartPropertyPage::DoDataExchange(CDataExchange* pDX)
{
	CMFCPropertyPage::DoDataExchange(pDX);

	DDX_Control(pDX, IDC_EDIT3, m_minEnergy);
	DDX_Control(pDX, IDC_EDIT4, m_maxEnergy);


	DDX_Text(pDX, IDC_EDIT1, XBigTicksBondLength);
	DDX_Text(pDX, IDC_EDIT2, XSmallTicksBondLength);
	DDX_Text(pDX, IDC_EDIT3, YMinEnergy);
	DDX_Text(pDX, IDC_EDIT4, YMaxEnergy);
	DDX_Text(pDX, IDC_EDIT5, YBigTicksEnergy);
	DDX_Text(pDX, IDC_EDIT6, YSmallTicksEnergy);

	DDX_Check(pDX, IDC_CHECK2, useSplines);

	DDX_Radio(pDX, IDC_RADIO1, DisplayHOMOEnergy);
}


void ChartPropertyPage::OnEnChangeEdit1()
{
	SetModified();
}


void ChartPropertyPage::OnEnChangeEdit2()
{
	SetModified();
}


void ChartPropertyPage::OnEnChangeEdit3()
{
	SetModified();
}


void ChartPropertyPage::OnEnChangeEdit4()
{
	SetModified();
}


void ChartPropertyPage::OnEnChangeEdit5()
{
	SetModified();
}


void ChartPropertyPage::OnEnChangeEdit6()
{
	SetModified();
}


void ChartPropertyPage::OnBnClickedCheck2()
{
	SetModified();
}


void ChartPropertyPage::OnBnClickedRadio1()
{
	SetModified();
}


void ChartPropertyPage::OnBnClickedRadio2()
{
	SetModified();
}
