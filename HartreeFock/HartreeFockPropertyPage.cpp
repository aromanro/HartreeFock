// HartreeFockPropertyPage.cpp : implementation file
//

#include "stdafx.h"
#include "HartreeFock.h"
#include "HartreeFockPropertyPage.h"



// HartreeFockPropertyPage

IMPLEMENT_DYNAMIC(HartreeFockPropertyPage, CMFCPropertyPage)

HartreeFockPropertyPage::HartreeFockPropertyPage()
	: CMFCPropertyPage(IDD_HARTREEFOCKPROPERTYPAGE)
{
	m_Method = (theApp.options.restricted ? 0 : 1);
	m_Alpha = theApp.options.alpha;
	m_Guess = theApp.options.initialGuess;
	iterations = theApp.options.iterations;
}

HartreeFockPropertyPage::~HartreeFockPropertyPage()
{
}


BEGIN_MESSAGE_MAP(HartreeFockPropertyPage, CMFCPropertyPage)
	ON_BN_CLICKED(IDC_RADIO1, &HartreeFockPropertyPage::OnBnClickedRadio1)
	ON_BN_CLICKED(IDC_RADIO2, &HartreeFockPropertyPage::OnBnClickedRadio2)
	ON_EN_CHANGE(IDC_EDIT1, &HartreeFockPropertyPage::OnEnChangeEdit1)
	ON_EN_CHANGE(IDC_EDIT2, &HartreeFockPropertyPage::OnEnChangeEdit2)
	ON_EN_CHANGE(IDC_EDIT3, &HartreeFockPropertyPage::OnEnChangeEdit3)
END_MESSAGE_MAP()



// HartreeFockPropertyPage message handlers




BOOL HartreeFockPropertyPage::OnApply()
{
	UpdateData();

	ApplyValues();

	return CMFCPropertyPage::OnApply();
}


void HartreeFockPropertyPage::ApplyValues()
{
	theApp.options.restricted = (m_Method == 0 ? true : false);
	theApp.options.alpha = m_Alpha;
	theApp.options.initialGuess = m_Guess;
	theApp.options.iterations = iterations;

	theApp.options.Save();
}


BOOL HartreeFockPropertyPage::OnInitDialog()
{
	CMFCPropertyPage::OnInitDialog();

	// TODO:  Add extra initialization here
	m_AlphaEdit.allowNegative = false;
	m_GuessEdit.allowNegative = false;

	return TRUE;  // return TRUE unless you set the focus to a control
				  // EXCEPTION: OCX Property Pages should return FALSE
}


void HartreeFockPropertyPage::OnBnClickedRadio1()
{
	SetModified();
}


void HartreeFockPropertyPage::OnBnClickedRadio2()
{
	SetModified();
}


void HartreeFockPropertyPage::OnEnChangeEdit1()
{
	SetModified();
}


void HartreeFockPropertyPage::OnEnChangeEdit2()
{
	SetModified();
}


void HartreeFockPropertyPage::DoDataExchange(CDataExchange* pDX)
{
	// TODO: Add your specialized code here and/or call the base class

	CMFCPropertyPage::DoDataExchange(pDX);
	
	DDX_Control(pDX, IDC_EDIT1, m_AlphaEdit);
	DDX_Control(pDX, IDC_EDIT2, m_GuessEdit);

	DDX_Radio(pDX, IDC_RADIO1, m_Method);
	DDX_Text(pDX, IDC_EDIT1, m_Alpha);
	DDX_Text(pDX, IDC_EDIT2, m_Guess);
	DDX_Text(pDX, IDC_EDIT3, iterations);

	DDV_MinMaxDouble(pDX, m_Alpha, 0.01, 1.);
	DDV_MinMaxDouble(pDX, m_Guess, 0., 10.);
	DDV_MinMaxUInt(pDX, iterations, 50, 50000);
}


void HartreeFockPropertyPage::OnEnChangeEdit3()
{
	SetModified();
}
