// File created: 12/16/2002                          last modified: 11/12/2007
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2008 Dr. Christian Stratowa
 *
 *  Written by: Christian Stratowa, Vienna, Austria <cstrato@aon.at>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy of the GNU General Public
 *  License is available at http://www.gnu.org/copyleft/gpl.html. You
 *  can also obtain it by writing to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA.
 *
 *******************************************************************************
 * Based on: "The ROOT System", http://root.cern.ch/                           *
 * ROOT:     An Object-Oriented Data Analysis Framework                        *
 * Authors:  Rene Brun and Fons Rademakers.                                    *
 * For the licensing terms of "The ROOT System" see $ROOTSYS/LICENSE.          *
 * For the list of contributors to "The ROOT System" see http://root.cern.ch/  *
 *******************************************************************************
 */

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include "XPSAnalyzer.h"

#include "StatUtils.h"

#include "TString.h"


//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XAnalyser);
ClassImp(XUniTester);
ClassImp(XMultiTester);
ClassImp(XClusterizer);
ClassImp(XRegressor);



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniTester                                                           //
//                                                                      //
// Class for analysis of normalized data                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUniTester::XUniTester()
           :XAnalyser()
{
   // Default UniTester constructor
   if(kCS) cout << "---XUniTester::XUniTester(default)------" << endl;

   fUniTest = 0;
}//Constructor

//______________________________________________________________________________
XUniTester::XUniTester(const char *name, const char *type)
           :XAnalyser(name, type)
{
   // Normal UniTester constructor
   if(kCS) cout << "---XUniTester::XUniTester------" << endl;

   fUniTest = 0;
}//Constructor

//______________________________________________________________________________
XUniTester::~XUniTester()
{
   // UniTester destructor
   if(kCS) cout << "---XUniTester::~XUniTester------" << endl;

   SafeDelete(fUniTest);
}//Destructor

//______________________________________________________________________________
Int_t XUniTester::InitType(const char *type, Option_t *options, Int_t npars,
                  Double_t *pars)
{
   // Initialize analyser type
   if(kCS) cout << "------XUniTester::InitType------" << endl;

// Delete default analysis setting
   SafeDelete(fUniTest);

   if (npars != 5) return errInitSetting;

	TString optcpy = options;
   char   *opt    = (char*)optcpy.Data();
   TString opt1   = strtok(opt, ":");
   TString opt2   = strtok(NULL, ":");

// Create and initialize test class
   if (strcmp(type, "normaltest") == 0) {
      if (fHasNA) fUniTest = new TUnivariateTest(type, kExtenUTst[0], fNA);
      else        fUniTest = new TUnivariateTest(type, kExtenUTst[0]);
      // Init(alt,nperm, mu, paired, conflevel)
      fUniTest->Init((Int_t)pars[0],pars[1],(Bool_t)pars[2],pars[3],opt1.Data());
   } else if (strcmp(type, "ttest") == 0) {
      if (fHasNA) fUniTest = new TStudentTest(type, kExtenUTst[1], fNA);
      else        fUniTest = new TStudentTest(type, kExtenUTst[1]);
      // Init(alt,nperm, mu, paired, conflevel, varequ)
      ((TStudentTest*)fUniTest)->Init((Int_t)pars[0],pars[1],(Bool_t)pars[2],
                                      pars[3],(Bool_t)pars[4],opt1);
   } else if (strcmp(type, "wilcoxtest") == 0) {
//to do
      cout << "Note: Wilcox-Test not yet implemented" << endl;
      return errAbort;
   } else if (strcmp(type, "vartest") == 0) {
//to do
      cout << "Note: Variance-Test not yet implemented" << endl;
      return errAbort;
   } else {
      cerr << "Error: Analysis algorithm <" << type << "> not known" << endl;
      return errAbort;
   }//if
   if (fUniTest == 0) return errInitMemory;

//???? ev adjpval as parameter to import???
   Bool_t adjpval = (pars[0] > 0) ? kFALSE : kTRUE;
// Set adjustment option
   fUniTest->SetAdjustment(opt2.Data(), adjpval);

   return errNoErr;
}//InitType

//______________________________________________________________________________
Int_t XUniTester::Analyse(const char *infile, const char *outfile, 
                  const char *varlist, Int_t nrows, Int_t nperm, Double_t mu,
                  const char *sepi, const char *sepo, char delim)
{
   // Analyse data
   if(kCS) cout << "------XUniTester::Analyse------" << endl;

   return fUniTest->Test(infile, outfile, varlist, nrows, nperm, mu,
                    sepi, sepo, delim);
}//Analyse

//______________________________________________________________________________
Int_t XUniTester::Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                  TTree *outtree, const char *varlist, Int_t nperm, Double_t mu)
{
   // Analyse data
   // Note: Since outtree is filled in TUnivariateTest, mn1 and mn2 will be
   //       logarithmic if fLogBase != 0
   if(kCS) cout << "------XUniTester::Analyse------" << endl;

   return fUniTest->Test(n, gid, intree, leafname, outtree, varlist, nperm, mu);
}//Analyse


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiTester                                                         //
//                                                                      //
// Class for multivariate analysis of normalized data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMultiTester::XMultiTester()
             :XAnalyser()
{
   // Default MultiTester constructor
   if(kCS) cout << "---XMultiTester::XMultiTester(default)------" << endl;

//   fMultiTest = 0;
}//Constructor

//______________________________________________________________________________
XMultiTester::XMultiTester(const char *name, const char *type)
             :XAnalyser(name, type)
{
   // Normal MultiTester constructor
   if(kCS) cout << "---XMultiTester::XMultiTester------" << endl;

//   fMultiTest = 0;
}//Constructor

//______________________________________________________________________________
XMultiTester::~XMultiTester()
{
   // MultiTester destructor
   if(kCS) cout << "---XMultiTester::~XMultiTester------" << endl;

//   SafeDelete(fMultiTest);
}//Destructor

//______________________________________________________________________________
Int_t XMultiTester::InitType(const char *type, Option_t *options, Int_t npars,
                    Double_t *pars)
{
   // Initialize analyser type
   if(kCS) cout << "------XMultiTester::InitType------" << endl;

// Delete default analysis setting
//   SafeDelete(fMultiTest);

   return errNoErr;
}//InitType

//______________________________________________________________________________
Int_t XMultiTester::Analyse(const char *infile, const char *outfile, 
                    const char *varlist, Int_t nrows, Int_t nperm, Double_t mu,
                    const char *sepi, const char *sepo, char delim, Int_t linebuf)
{
   // Analyse data
   if(kCS) cout << "------XMultiTester::Analyse------" << endl;

//   return fMultiTest->Test(infile, outfile, varlist, nrows, nperm, mu,
//                      sepi, sepo, delim, linebuf);
   return errNoErr;
}//Analyse

//______________________________________________________________________________
Int_t XMultiTester::Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                    TTree *outtree, const char *varlist, Int_t nperm, Double_t mu)
{
   // Analyse data
   if(kCS) cout << "------XMultiTester::Analyse------" << endl;

//   return fMultiTest->Test(n, gid, intree, leafname, outtree, varlist, nperm, mu);
   return errNoErr;
}//Analyse


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XClusterizer                                                         //
//                                                                      //
// Class for multivariate analysis of normalized data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XClusterizer::XClusterizer()
             :XAnalyser()
{
   // Default Clusterizer constructor
   if(kCS) cout << "---XClusterizer::XClusterizer(default)------" << endl;

//   fCluster = 0;
}//Constructor

//______________________________________________________________________________
XClusterizer::XClusterizer(const char *name, const char *type)
             :XAnalyser(name, type)
{
   // Normal Clusterizer constructor
   if(kCS) cout << "---XClusterizer::XClusterizer------" << endl;

//   fCluster = 0;
}//Constructor

//______________________________________________________________________________
XClusterizer::~XClusterizer()
{
   // Clusterizer destructor
   if(kCS) cout << "---XClusterizer::~XClusterizer------" << endl;

//   SafeDelete(fCluster);
}//Destructor

//______________________________________________________________________________
Int_t XClusterizer::InitType(const char *type, Option_t *options, Int_t npars,
                    Double_t *pars)
{
   // Initialize analyser type
   if(kCS) cout << "------XClusterizer::InitType------" << endl;

// Delete default analysis setting
//   SafeDelete(fCluster);

   return errNoErr;
}//InitType

//______________________________________________________________________________
Int_t XClusterizer::Analyse(const char *infile, const char *outfile, 
                    const char *varlist, Int_t nrows, Int_t nperm, Double_t mu,
                    const char *sepi, const char *sepo, char delim, Int_t linebuf)
{
   // Analyse data
   if(kCS) cout << "------XClusterizer::Analyse------" << endl;

//   return fCluster->Test(infile, outfile, varlist, nrows, nperm, mu,
//                    sepi, sepo, delim, linebuf);
   return errNoErr;
}//Analyse

//______________________________________________________________________________
Int_t XClusterizer::Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                    TTree *outtree, const char *varlist, Int_t nperm, Double_t mu)
{
   // Analyse data
   if(kCS) cout << "------XClusterizer::Analyse------" << endl;

//   return fCluster->Test(n, gid, intree, leafname, outtree, varlist, nperm, mu);
   return errNoErr;
}//Analyse


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRegressor                                                           //
//                                                                      //
// Class for regression analysis of normalized data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRegressor::XRegressor()
           :XAnalyser()
{
   // Default Regressor constructor
   if(kCS) cout << "---XRegressor::XRegressor(default)------" << endl;

//   fRegression = 0;
}//Constructor

//______________________________________________________________________________
XRegressor::XRegressor(const char *name, const char *type)
           :XAnalyser(name, type)
{
   // Normal Regressor constructor
   if(kCS) cout << "---XRegressor::XRegressor------" << endl;

//   fRegression = 0;
}//Constructor

//______________________________________________________________________________
XRegressor::~XRegressor()
{
   // Regressor destructor
   if(kCS) cout << "---XRegressor::~XRegressor------" << endl;

//   SafeDelete(fRegression);
}//Destructor

//______________________________________________________________________________
Int_t XRegressor::InitType(const char *type, Option_t *options, Int_t npars,
                  Double_t *pars)
{
   // Initialize analyser type
   if(kCS) cout << "------XRegressor::InitType------" << endl;

// Delete default analysis setting
//   SafeDelete(fRegression);

   return errNoErr;
}//InitType

//______________________________________________________________________________
Int_t XRegressor::Analyse(const char *infile, const char *outfile, 
                  const char *varlist, Int_t nrows, Int_t nperm, Double_t mu,
                  const char *sepi, const char *sepo, char delim, Int_t linebuf)
{
   // Analyse data
   if(kCS) cout << "------XRegressor::Analyse------" << endl;

//   return fRegression->Test(infile, outfile, varlist, nrows, nperm, mu,
//                       sepi, sepo, delim, linebuf);
   return errNoErr;
}//Analyse

//______________________________________________________________________________
Int_t XRegressor::Analyse(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                  TTree *outtree, const char *varlist, Int_t nperm, Double_t mu)
{
   // Analyse data
   if(kCS) cout << "------XRegressor::Analyse------" << endl;

//   return fRegression->Test(n, gid, intree, leafname, outtree, varlist, nperm, mu);
   return errNoErr;
}//Analyse








