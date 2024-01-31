
// HartreeFock.h : main header file for the HartreeFock application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols

#include <objidl.h>
#include <gdiplus.h>

#pragma comment (lib,"Gdiplus.lib")

#include "Options.h"


// CHartreeFockApp:
// See HartreeFock.cpp for the implementation of this class
//

class CHartreeFockApp : public CWinAppEx
{
public:
	CHartreeFockApp();


// Overrides
	BOOL InitInstance() override;

// Implementation

	Options options;

	UINT  m_nAppLook;
	BOOL  m_bHiColorIcons;

private:
	ULONG_PTR gdiplusToken;

	void PreLoadState() override;
	void LoadCustomState() override;
	void SaveCustomState() override;
	int ExitInstance() override;

	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CHartreeFockApp theApp;
