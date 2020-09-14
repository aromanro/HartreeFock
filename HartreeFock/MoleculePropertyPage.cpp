// MoleculePropertyPage.cpp : implementation file
//

#include "stdafx.h"

#include "HartreeFock.h"
#include "MoleculePropertyPage.h"
#include "MainFrm.h"
#include "Basis.h"
#include "ChemUtils.h"
#include "HartreeFockDoc.h"


// MoleculePropertyPage

IMPLEMENT_DYNAMIC(MoleculePropertyPage, CMFCPropertyPage)

MoleculePropertyPage::MoleculePropertyPage()
	: CMFCPropertyPage(IDD_MOLECULEPROPERTYPAGE)
{
	m_atom1 = theApp.options.m_atom1;
	m_atom2 = theApp.options.m_atom2;
	twoAtom1 = (theApp.options.twoAtom1 ? BST_CHECKED : BST_UNCHECKED);

	basis = theApp.options.basis; // 0 - STO3G, 1 - STO6G, 2 - 3-21G, 3 - 6-21G, 4 - 6-31G, 5 - 6-31G*, 6 - 6-31+G**
	bondAngle = theApp.options.bondAngle;

	alphaElectrons = theApp.options.alphaElectrons;
	betaElectrons = theApp.options.betaElectrons;

	// also used for charting
	XMaxBondLength = theApp.options.XMaxBondLength;
	XMinBondLength = theApp.options.XMinBondLength; // not really for chart, but for calculations
}

MoleculePropertyPage::~MoleculePropertyPage()
{
}


BEGIN_MESSAGE_MAP(MoleculePropertyPage, CMFCPropertyPage)
	ON_CBN_SELCHANGE(IDC_COMBO1, &MoleculePropertyPage::OnCbnSelchangeCombo1)
	ON_CBN_SELCHANGE(IDC_COMBO2, &MoleculePropertyPage::OnCbnSelchangeCombo2)
	ON_BN_CLICKED(IDC_CHECK1, &MoleculePropertyPage::OnBnClickedCheck1)
	ON_EN_CHANGE(IDC_EDIT1, &MoleculePropertyPage::OnEnChangeEdit1)
	ON_BN_CLICKED(IDC_RADIO1, &MoleculePropertyPage::OnBnClickedRadio1)
	ON_BN_CLICKED(IDC_RADIO2, &MoleculePropertyPage::OnBnClickedRadio2)
	ON_EN_CHANGE(IDC_EDIT5, &MoleculePropertyPage::OnEnChangeEdit5)
	ON_EN_CHANGE(IDC_EDIT6, &MoleculePropertyPage::OnEnChangeEdit6)
	ON_EN_CHANGE(IDC_EDIT7, &MoleculePropertyPage::OnEnChangeEdit7)
	ON_EN_CHANGE(IDC_EDIT8, &MoleculePropertyPage::OnEnChangeEdit8)
	ON_BN_CLICKED(IDC_RADIO4, &MoleculePropertyPage::OnBnClickedRadio4)
	ON_BN_CLICKED(IDC_RADIO5, &MoleculePropertyPage::OnBnClickedRadio5)
	ON_BN_CLICKED(IDC_RADIO6, &MoleculePropertyPage::OnBnClickedRadio6)
	ON_BN_CLICKED(IDC_RADIO7, &MoleculePropertyPage::OnBnClickedRadio7)
	ON_BN_CLICKED(IDC_RADIO8, &MoleculePropertyPage::OnBnClickedRadio8)
END_MESSAGE_MAP()



// MoleculePropertyPage message handlers




BOOL MoleculePropertyPage::OnApply()
{
	UpdateData();

	ApplyValues();

	return CMFCPropertyPage::OnApply();
}


void MoleculePropertyPage::ApplyValues()
{
	theApp.options.m_atom1 = m_atom1;
	theApp.options.m_atom2 = m_atom2;
	theApp.options.twoAtom1 = (twoAtom1 == BST_CHECKED ? true : false);

	theApp.options.basis = basis; // 0 - STO3G, 1 - STO6G, 2 - 3-21G, 3 - 6-21G, 4 - 6-31G, 5 - 6-31G*, 6 - 6-31+G**
	theApp.options.bondAngle = bondAngle;

	theApp.options.alphaElectrons = alphaElectrons;
	theApp.options.betaElectrons = betaElectrons;

	// also used for charting
	theApp.options.XMaxBondLength = XMaxBondLength;
	theApp.options.XMinBondLength = XMinBondLength; // not really for chart, but for calculations

	theApp.options.Save();

	dynamic_cast<CMainFrame*>(theApp.m_pMainWnd)->GetDocument()->ApplyChartOptions();
}


BOOL MoleculePropertyPage::OnInitDialog()
{
	CMFCPropertyPage::OnInitDialog();

	bondAngleEdit.allowNegative = true;
	minBondEdit.allowNegative = false;
	maxBondEdit.allowNegative = false;

	FillCombos();

	bondAngleEdit.EnableWindow(twoAtom1 == BST_CHECKED);

	return TRUE;  // return TRUE unless you set the focus to a control
				  // EXCEPTION: OCX Property Pages should return FALSE
}


void MoleculePropertyPage::DoDataExchange(CDataExchange* pDX)
{
	CMFCPropertyPage::DoDataExchange(pDX);

	DDX_Control(pDX, IDC_EDIT1, bondAngleEdit);
	DDX_Control(pDX, IDC_EDIT7, minBondEdit);
	DDX_Control(pDX, IDC_EDIT8, maxBondEdit);
	DDX_Control(pDX, IDC_COMBO1, comboBox1);
	DDX_Control(pDX, IDC_COMBO2, comboBox2);

	DDX_Text(pDX, IDC_EDIT1, bondAngle);
	DDX_Text(pDX, IDC_EDIT5, alphaElectrons);
	DDX_Text(pDX, IDC_EDIT6, betaElectrons);
	DDX_Text(pDX, IDC_EDIT7, XMinBondLength);
	DDX_Text(pDX, IDC_EDIT8, XMaxBondLength);

	DDX_Check(pDX, IDC_CHECK1, twoAtom1);
	DDX_Radio(pDX, IDC_RADIO1, basis);

	DDX_CBStringExact(pDX, IDC_COMBO1, m_atom1);
	DDX_CBStringExact(pDX, IDC_COMBO2, m_atom2);
}


void MoleculePropertyPage::OnCbnSelchangeCombo1()
{	
	AdjustElectrons();

	SetModified();

	// post a message to adjust electrons
}


void MoleculePropertyPage::OnCbnSelchangeCombo2()
{
	AdjustElectrons();

	SetModified();

	// post a message to adjust electrons
}


void MoleculePropertyPage::OnBnClickedCheck1()
{
	AdjustElectrons();

	SetModified();

	bondAngleEdit.EnableWindow(twoAtom1 == BST_CHECKED);
}


void MoleculePropertyPage::OnEnChangeEdit1()
{
	SetModified();
}


void MoleculePropertyPage::OnBnClickedRadio1()
{
	UpdateData();
	ApplyValues();

	FillCombos();

	SetModified();
}


void MoleculePropertyPage::OnBnClickedRadio2()
{
	OnBnClickedRadio1();
}


void MoleculePropertyPage::OnEnChangeEdit5()
{
	SetModified();
}


void MoleculePropertyPage::OnEnChangeEdit6()
{
	SetModified();
}


void MoleculePropertyPage::OnEnChangeEdit7()
{
	SetModified();
}


void MoleculePropertyPage::OnEnChangeEdit8()
{
	SetModified();
}


void MoleculePropertyPage::FillCombos()
{
	comboBox1.ResetContent();
	comboBox2.ResetContent();

	CHartreeFockDoc* doc = dynamic_cast<CMainFrame*>(theApp.m_pMainWnd)->GetDocument();
	Chemistry::Basis *basisPtr = &doc->basisSTO6G; // some default

	if (0 == theApp.options.basis)
	{
		basisPtr = &doc->basisSTO3G;
	}
	else if (1 == theApp.options.basis)
	{
		basisPtr = &doc->basisSTO6G;
	}
	else if (2 == theApp.options.basis)
	{
		basisPtr = &doc->basis3_21G;
	}
	else if (3 == theApp.options.basis)
	{
		basisPtr = &doc->basis6_21G;
	}
	else if (4 == theApp.options.basis)
	{
		basisPtr = &doc->basis6_31G;
	}	
	else if (5 == theApp.options.basis)
	{
		basisPtr = &doc->basis6_31Gstar;
	}
	else if (6 == theApp.options.basis)
	{
		basisPtr = &doc->basis6_31plusGstarstar;
	}	

	bool atom1Found = false;
	bool atom2Found = false;

	for (const auto &atom : basisPtr->atoms)
	{
		CString atomName(Chemistry::ChemUtils::GetAtomNameForZ(atom.Z).c_str());

		if (atomName == m_atom1) atom1Found = true;
		if (atomName == m_atom2) atom2Found = true;

		comboBox1.AddString(atomName);
		comboBox2.AddString(atomName);
	}

	if (!atom1Found) m_atom1 = L"H";
	if (!atom2Found) m_atom2 = L"O";

	UpdateData(false);
}


void MoleculePropertyPage::AdjustElectrons()
{
	UpdateData();

	CT2CA psz1(m_atom1);
	std::string str1(psz1);
	const unsigned int Z1 = Chemistry::ChemUtils::GetZForAtom(str1);

	CT2CA psz2(m_atom2);
	std::string str2(psz2);
	const unsigned int Z2 = Chemistry::ChemUtils::GetZForAtom(str2);

	const int totalElectrons = (twoAtom1 ? 2 : 1) * Z1 + Z2;

	alphaElectrons = totalElectrons / 2;
	betaElectrons = alphaElectrons;

	if (alphaElectrons + betaElectrons < totalElectrons) ++betaElectrons;

	UpdateData(FALSE);

	SetModified();
}

void MoleculePropertyPage::OnBnClickedRadio4()
{
	OnBnClickedRadio1();
}


void MoleculePropertyPage::OnBnClickedRadio5()
{
	OnBnClickedRadio1();
}


void MoleculePropertyPage::OnBnClickedRadio6()
{
	OnBnClickedRadio1();
}

void MoleculePropertyPage::OnBnClickedRadio7()
{
	OnBnClickedRadio1();
}


void MoleculePropertyPage::OnBnClickedRadio8()
{
	OnBnClickedRadio1();
}

