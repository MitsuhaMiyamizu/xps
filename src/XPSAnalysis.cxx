// File created: 12/16/2002                          last modified: 06/20/2010
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2010 Dr. Christian Stratowa
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

/******************************************************************************
* Major Revision History:
* Dec 2002 - Initial versions finished
*          - Completely reorganize source file structure
* Apr 2004 - Add support for filtering: XFilter, XFilterSet, XFilterHeader
* May 2004 - Put XFilter in own source file
* Jun 2004 - Add classes for algorithms XAnalyser, XUniTester
* May 2005 - Replace XContent with XFolder.
* Jul 2005 - Replace TTree::AddFriend() with TList fTrees.
*          - Store trees for every XTreeSet in separate TDirectory.
*
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include "TBranch.h"
#include "THashTable.h"
#include "TLeaf.h"
#include "TFriendElement.h"
#include "TROOT.h"

#include "StatUtils.h"

#include "XPSAnalysis.h"
#include "XPSFilter.h"
#include "XPSAnalyzer.h"

// Tree header type for XUniFilterHeader
enum EHeaderType {
   HEADER_UNITEST     = -1001,
   HEADER_UNIFILTER   = -1002,
   HEADER_MULTITEST   = -1003,
   HEADER_MULTIFILTER = -1004,
};

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XAnalysisManager);
ClassImp(XAnalySetting);
ClassImp(XPreFilterHeader);
ClassImp(XUniFilterHeader);
ClassImp(XMultiFilterHeader);
ClassImp(XUniTestHeader);
ClassImp(XMultiTestHeader);
ClassImp(XClusterHeader);
ClassImp(XRegressionHeader);
ClassImp(XAnalySet);
ClassImp(XPreFilterSet);
ClassImp(XUnivarSet);
ClassImp(XMultivarSet);
ClassImp(XClusterSet);
ClassImp(XRegressionSet);
ClassImp(XScore);
ClassImp(XGrpMn);
ClassImp(XChance);
ClassImp(XAdjust);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalysisManager                                                     //
//                                                                      //
// Class for statistical analysis of microarray data                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAnalysisManager::XAnalysisManager()
                 :XProcessManager()
{
   // Default AnalysisManager constructor
   if(kCS) cout << "---XAnalysisManager::XAnalysisManager(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XAnalysisManager::XAnalysisManager(const char *name, const char *type, Int_t verbose)
                 :XProcessManager(name, type, verbose)
{
   // Normal AnalysisManager constructor
   if(kCS) cout << "---XAnalysisManager::XAnalysisManager------" << endl;

}//Constructor

//______________________________________________________________________________
XAnalysisManager::~XAnalysisManager()
{
   // AnalysisManager destructor
   if(kCS) cout << "---XAnalysisManager::~XAnalysisManager------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XAnalysisManager::InitSetting(Int_t min, const char *logbase, Double_t neglog)
{
   // Re-initialize settings for filters:
   // - clear filters, i.e. reset filter methods to kFALSE
   // - min: minimum number of selected filter methods to satisfy
   //        min = 0: filter tree(s) will be ignored, i.e. all tree entries
   //                 will be used for the analysis
   //        min = 1: equivalent to "OR", i.e. satisfy at least one filter
   //        min = 9: equivalent to "AND", i.e. satisfy all 9 selected filters
   //        e.g. 5 filters selected and min=3: 3 of 5 filters must be satisfied
   // - logbase: convert data to logarithm of base: "0", "log", "log2", "log10"
   // - neglog:  positive value to replace negative values for logarithm
   // Note: if InitSetting() is not called then the default setting is:
   //       min=0:       if no filter tree was added using AddTree() 
   //       min=1:       if at least one filter tree was added using AddTree() 
   //       logbase="0": linear tree entries will be used 
   // Note: after resetting, filter methods must be initialized explicitely
   // Note: use InitSetting() before calling Export() in following case:
   //       - logbase: to convert tree data from logbase before exporting
   if(kCS) cout << "------XAnalysisManager::InitSetting------" << endl;

   if (fAbort) return errAbort;
   if (!fSetting) {HandleError(errInitSetting); return errAbort;}

   return ((XAnalySetting*)fSetting)->InitFilters(min, logbase, neglog);
}//InitSetting

//______________________________________________________________________________
Int_t XAnalysisManager::Analyse(const char *setname, const char *leafname, 
                        const char *outtree, const char *varlist)
{
//BETTER:
// (..., const char *filter)
// e.g. filter = "mn2/mn1 >= 1.3 && pval <= 0.01"
   if(kCS) cout << "------XAnalysisManager::Analyse(tree)------" << endl;

   if (fAbort) return errAbort;

// Create directory with name of treeset
   TDirectory *dir = fFile->GetDirectory(setname);
   if (!dir) fFile->mkdir(setname, fDataType);
   fFile->cd();

// Get treeset "setname" (set created in AddTree())
   XAnalySet *set = 0;
   set = (XAnalySet*)fContent->FindObject(setname);
   if (!set) {
      return HandleError(errGetTreeSet, setname);
   }//if

// Set must be inherited from XAnalySet
   if (!set->InheritsFrom("XAnalySet")) {
      return HandleError(errClassTreeSet, setname, set->ClassName());
   }//if

   Int_t err = errNoErr;
   if (set->GetNumSelections() > 0) {
      TString vars = varlist; vars.ToLower();

      if (!err) err = set->Initialize(fFile, fSetting);
      if (!err) err = set->Analyse(leafname, outtree, vars.Data());

      HandleError(err, "in XAnalysisManager::Analyse");
   } else {
      cerr << "Error: No tree selected." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

//not allowed, deleted in XManager::Close
//   SafeDelete(set);

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XAnalysisManager::Analyse(const char *infile, const char *outfile,
                        const char *varlist, Int_t nrows, const char *sepi,
                        const char *sepo, char delim)
{
   // 1. Prefiltering:
   // xxxx
   //
   // 2. Univariate analysis:
   // Do multiple univariate tests for infile and export result as outfile.
   // infile must contain:  1st colum: row name or identifier
   //    1st row:   header row
   //    2nd row:   "Group" in 1st column and group id 1, 2,...
   //               two-sample test: arbitrary distribution of groups allowed
   //               paired test of n pairs: pair is in col[i] and col[n+i]
   //                                   or: pair is in col[i] and col[i+1]
   //               if 2nd row is missing, a one-sample test is done
   // varlist is the list of variables to calculate and to export:
   //    "stat" - univariate statistic
   //    "mn1"  - mean value of group1
   //    "mn2"  - mean value of group2  (not used for one-sample test)
   //    "se"   - standard error
   //    "df"   - degrees of freedom
   //    "pval" - p-value
   //    "nper" - number of permutations used to determine p-chance
   //    "pcha" - p-chance
   //    "padj" - p-value or p-chance adjusted for multiple comparisons
   //    e.g. for usual test, varlist = "stat:mn1:mn2:se:df:pval"
   //    "*" is equal to "stat:mn1:mn2:se:df:pval:nper:pcha:padj"
   // nrows is number of data rows (minus header and group)
   // Note: No need to create root file and to open data/scheme files
   if(kCS) cout << "------XAnalysisManager::Analyse(table)------" << endl;

   if (fAbort) return errAbort;

   XAnalySet *set = 0;
   set = (XAnalySet*)(this->NewTreeSet(GetTitle()));
   if (!set) {
      return HandleError(errCreateTreeSet, GetTitle());
   }//if

   TString vars = varlist; vars.ToLower();

   Int_t err = errNoErr;
   if (!err) err = set->Initialize(fFile, fSetting, infile);
   if (!err) err = set->Analyse(infile, outfile, vars.Data(), nrows, sepi,
                        sepo, delim);
   else      HandleError(err, "in XAnalysisManager::Analyse");

//not allowed, deleted in XManager::Close
//   SafeDelete(set);

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XAnalysisManager::ImportDefaults(const char *infile)
{
   // Import default values from infile
   if(kCS) cout << "------XAnalysisManager::ImportDefaults------" << endl;

   return errNoErr;
}//ImportDefaults

///////////////
//BETTER: delete InitDefaults() since InitAlgorithm() must be called in macros!?
//______________________________________________________________________________
Int_t XAnalysisManager::InitDefaults()
{
   // Initialize default values
   if(kCS) cout << "------XAnalysisManager::InitDefaults------" << endl;

   Int_t err = errNoErr;

// Initialize algorithms
   if (strcmp(fSetting->GetName(),"PreFilter") == 0) {
      /*ok*/;
//ev delete?? since problem: no SafeDelete(fFilter)!!
//      err = InitAlgorithm("PreFilter","variation","var2mn",0, 3, 0.3, 0.0, 0.01);
//      if (!err) err = InitAlgorithm("prefilter","LowerThreshold","percent",0, 2, 5.0, 100.0);
//      if (!err) err = InitAlgorithm("prefilter","UpperThreshold","percent",0, 2, 14.5, 100.0);
   } else if (strcmp(fSetting->GetName(),"UnivariateAnalysis") == 0) {
//??      err = InitAlgorithm("UniTest","ttest", "twosided:none",0, 5, 0, 0, 0, 0.95, 0);
   } else if (strcmp(fSetting->GetName(),"MultivariateAnalysis") == 0) {
//      err = InitAlgorithm("MultiTest");
   } else if (strcmp(fSetting->GetName(),"ClusterAnalysis") == 0) {
//      err = InitAlgorithm("Cluster");
   } else if (strcmp(fSetting->GetName(),"RegressionAnalysis") == 0) {
//      err = InitAlgorithm("Regression");
   } else {
      err = HandleError(errUnknownType, "Analysis", fSetting->GetName());
   }//if

   return err;
}//InitDefaults

//______________________________________________________________________________
XSetting *XAnalysisManager::NewSetting(const char *type, const char *infile)
{
   // Create new setting for type
   if(kCS) cout << "------XAnalysisManager::NewSetting------" << endl;

   XAnalySetting *setting = new XAnalySetting(type, infile);
   return setting;
}//NewSetting

//______________________________________________________________________________
XTreeSet *XAnalysisManager::NewTreeSet(const char *type)
{
   // Create new tree set of analysis type
   if(kCS) cout << "------XAnalysisManager::NewTreeSet------" << endl;

//NOTE: Suggested class types, could be subject to change!!
// type should/could be analysis type and NOT chip type,
// since analysis is independent of chip type??
   XTreeSet *set = 0;
   if (strcmp(type,"PreFilter") == 0) {
      set = new XPreFilterSet("pfr", type);
   } else if (strcmp(type,"UnivariateAnalysis") == 0) {
      set = new XUnivarSet("uvt", type);
   } else if (strcmp(type,"MultivariateAnalysis") == 0) {
      set = new XMultivarSet("mvt", type);
   } else if (strcmp(type,"ClusterAnalysis") == 0) {
      set = new XClusterSet("cls", type);
   } else if (strcmp(type,"RegressionAnalysis") == 0) {
      set = new XRegressionSet("reg", type);
   } else {
      HandleError(errUnknownType, "Analysis", type);
   }//if

   return set;
}//NewTreeSet

//______________________________________________________________________________
XPlot *XAnalysisManager::NewPlotter(const char *name, const char *title)
{
   // Create new plotter
   if(kCS) cout << "------XAnalysisManager::NewPlotter------" << endl;

//NOTE: Possible class types, could be subject to change!!
   XPlot *plotter = 0;
   if (strcmp(title,"PreFilter") == 0) {
      plotter = new XAnalysisPlot(name, title);
   } else if (strcmp(title,"UnivariateAnalysis") == 0) {
      plotter = new XAnalysisPlot(name, title);
//??      plotter = new XUnivarPlot(name, title);
   } else if (strcmp(title,"MultivariateAnalysis") == 0) {
      plotter = new XAnalysisPlot(name, title);
//??      plotter = new XMultivarPlot(name, title);
   } else if (strcmp(title,"ClusterAnalysis") == 0) {
      plotter = new XAnalysisPlot(name, title);
//??      plotter = new XClusterPlot(name, title);
   } else if (strcmp(title,"RegressionAnalysis") == 0) {
      plotter = new XAnalysisPlot(name, title);
//??      plotter = new XRegressionPlot(name, title);
   } else {
      HandleError(errUnknownType, "Analysis", title);
   }//if
   return plotter;
}//NewPlotter


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalySetting                                                        //
//                                                                      //
// Class for initialization of parameter settings for analysis          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAnalySetting::XAnalySetting()
              :XProcesSetting()
{
   // Default AnalySetting constructor
   if(kCS) cout << "---XAnalySetting::XAnalySetting(default)------" << endl;

   fFilter     = 0;
   fAnalyser   = 0;
   fMinFilters = 1;
   fLogBase    = "0";
}//Constructor

//______________________________________________________________________________
XAnalySetting::XAnalySetting(const char *type, const char *infile)
              :XProcesSetting(type, infile)
{
   // Normal AnalySetting constructor
   if(kCS) cout << "---XAnalySetting::XAnalySetting------" << endl;

   fFilter     = 0;
   fAnalyser   = 0;
   fMinFilters = 1;
   fLogBase    = "0";
}//Constructor

//______________________________________________________________________________
XAnalySetting::~XAnalySetting()
{
   // AnalySetting destructor
   if(kCS) cout << "---XAnalySetting::~XAnalySetting------" << endl;

   SafeDelete(fAnalyser);
   SafeDelete(fFilter);
}//Destructor

//______________________________________________________________________________
Int_t XAnalySetting::InitAlgorithm(const char *name, const char *type,
                     Option_t *options, const char * /*filename*/,
                     Int_t npars, Double_t *pars)
{
   // Initialize algorithm "name" with "type" and "options"
   if(kCS) cout << "------XAnalySetting::InitAlgorithm------" << endl;

   if (strcmp(name, "prefilter") == 0) {
      return this->InitPreFilter(type, options, npars, pars);
   } else if (strcmp(name, "unifilter")   == 0) {
      return this->InitUniFilter(type, options, npars, pars);
   } else if (strcmp(name, "multifilter")   == 0) {
      return this->InitMultiFilter(type, options, npars, pars);
   } else if (strcmp(name, "unitest")   == 0) {
      return this->InitUniTest(type, options, npars, pars);
   } else if (strcmp(name, "multitest") == 0) {
      return this->InitMultiTest(type, options, npars, pars);
   } else if (strcmp(name, "cluster") == 0) {
      return this->InitClusterizer(type, options, npars, pars);
   } else if (strcmp(name, "regression") == 0) {
      return this->InitRegressor(type, options, npars, pars);
   }//if

   return errInitSetting;
}//InitAlgorithm

//______________________________________________________________________________
void XAnalySetting::ResetAlgorithm(const char *name, const char *type)
{
   // Reset algorithm with name and optional with type
   if(kCS) cout << "------XAnalySetting::ResetAlgorithm------" << endl;

   if ((strcmp(name, "prefilter")   == 0)  ||
       (strcmp(name, "unifilter")   == 0)  ||
       (strcmp(name, "multifilter") == 0)) {
      SafeDelete(fFilter);
   } else if ((strcmp(name, "unitest")    == 0)  ||
              (strcmp(name, "multitest")  == 0)  ||
              (strcmp(name, "cluster")    == 0)  ||
              (strcmp(name, "regression") == 0)) {
      SafeDelete(fAnalyser);
   }//if
}//ResetAlgorithm

//______________________________________________________________________________
Int_t XAnalySetting::InitFilters(Int_t min, const char *logbase, Double_t neglog)
{
   // Re-initialize filters and reset filter methods to kFALSE
   if(kCS) cout << "------XAnalySetting::InitFilters------" << endl;

   fNegLog  = neglog;
   fLogBase = logbase;
   fLogBase.ToLower();

   fMinFilters = min;

   return errNoErr;
}//InitFilters

//______________________________________________________________________________
Int_t XAnalySetting::InitPreFilter(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize prefilter algorithm of type with options:
   if(kCS) cout << "------XAnalySetting::InitPreFilter------" << endl;

   Int_t err = errNoErr;

   if (!fFilter) fFilter = new XPreFilter(GetName(), kExtenFltr[0]);
   if (!fFilter) return errInitMemory;

   if (fHasNA) fFilter->InitNA(fNA);

   if (err == errNoErr) err = fFilter->Initialize(fMinFilters);
   if (err == errNoErr) err = fFilter->InitType(type, options, npars, pars);

   return err;
}//InitPreFilter

//______________________________________________________________________________
Int_t XAnalySetting::InitUniFilter(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize unifilter algorithm of type with options:
   if(kCS) cout << "------XAnalySetting::InitUniFilter------" << endl;

   Int_t err = errNoErr;

   if (!fFilter) fFilter = new XUniFilter(GetName(), kExtenFltr[1]);
   if (!fFilter) return errInitMemory;

   if (fHasNA) fFilter->InitNA(fNA);

   if (err == errNoErr) err = fFilter->Initialize(fMinFilters);
   if (err == errNoErr) err = fFilter->InitType(type, options, npars, pars);

   return err;
}//InitUniFilter

//______________________________________________________________________________
Int_t XAnalySetting::InitMultiFilter(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize multifilter algorithm of type with options:
   if(kCS) cout << "------XAnalySetting::InitMultiFilter------" << endl;

   Int_t err = errNoErr;

   if (!fFilter) fFilter = new XMultiFilter(GetName(), kExtenFltr[2]);
   if (!fFilter) return errInitMemory;

   if (fHasNA) fFilter->InitNA(fNA);

   if (err == errNoErr) err = fFilter->Initialize(fMinFilters);
   if (err == errNoErr) err = fFilter->InitType(type, options, npars, pars);

   return err;
}//InitMultiFilter

//______________________________________________________________________________
Int_t XAnalySetting::InitUniTest(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize univariate test algorithm of type with options:
   if(kCS) cout << "------XAnalySetting::InitUniTest------" << endl;

// Delete default analyser setting
   SafeDelete(fAnalyser);

   TString exten = Type2Extension(type, kTypeUTst, kExtenUTst);
   TString stype = Extension2Type(type, kTypeUTst, kExtenUTst);

   fAnalyser = new XUniTester(stype.Data(), exten.Data());
   if (!fAnalyser) return errInitMemory;

   if (fHasNA) fAnalyser->InitNA(fNA);
   return fAnalyser->InitType(type, options, npars, pars);
}//InitUniTest

//______________________________________________________________________________
Int_t XAnalySetting::InitMultiTest(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize multivariate test algorithm of type with options:
   if(kCS) cout << "------XAnalySetting::InitMultiTest------" << endl;

   TString exten = Type2Extension(type, kTypeMTst, kExtenMTst);
   TString stype = Extension2Type(type, kTypeMTst, kExtenMTst);

   fAnalyser = new XMultiTester(stype.Data(), exten.Data());
   if (!fAnalyser) return errInitMemory;

   if (fHasNA) fAnalyser->InitNA(fNA);
   return fAnalyser->InitType(type, options, npars, pars);
}//InitMultiTest

//______________________________________________________________________________
Int_t XAnalySetting::InitClusterizer(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize cluster analysis algorithms of type with options:
   if(kCS) cout << "------XAnalySetting::InitClusterizer------" << endl;

   TString exten = Type2Extension(type, kTypeClus, kExtenClus);
   TString stype = Extension2Type(type, kTypeClus, kExtenClus);

   fAnalyser = new XClusterizer(stype.Data(), exten.Data());
   if (!fAnalyser) return errInitMemory;

   if (fHasNA) fAnalyser->InitNA(fNA);
   return fAnalyser->InitType(type, options, npars, pars);
}//InitClusterizer

//______________________________________________________________________________
Int_t XAnalySetting::InitRegressor(const char *type, Option_t *options,
                     Int_t npars, Double_t *pars)
{
   // Initialize regression analysis algorithms of type with options:
   if(kCS) cout << "------XAnalySetting::InitRegressor------" << endl;

   TString exten = Type2Extension(type, kTypeRgrs, kExtenRgrs);
   TString stype = Extension2Type(type, kTypeRgrs, kExtenRgrs);

   fAnalyser = new XRegressor(stype.Data(), exten.Data());
   if (!fAnalyser) return errInitMemory;

   if (fHasNA) fAnalyser->InitNA(fNA);
   return fAnalyser->InitType(type, options, npars, pars);
}//InitRegressor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreFilterHeader                                                     //
//                                                                      //
// Class for storing nonspecific filter information in tree header      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPreFilterHeader::XPreFilterHeader()
                 :XTreeHeader()
{
   // Default PreFilterHeader constructor
   if(kCS) cout << "---XPreFilterHeader::XPreFilterHeader(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XPreFilterHeader::XPreFilterHeader(const char *str, Int_t treeid)
                 :XTreeHeader(str, treeid)
{
   // Normal PreFilterHeader constructor
   if(kCS) cout << "---XPreFilterHeader::XPreFilterHeader------" << endl;

   fMAD           = 0.0;
   fCov2mn        = 0.0;
   fVar2mn        = 0.0;
   fDif2mn        = 0.0;
   fMax2min       = 0.0;
   fGap2mn        = 0.0;
   fTrim          = 0.0;
   fWindow        = 0.0;
   fLoThreshold   = 0.0;
   fLoCondition   = "percent";
   fLoSamples     = 100.0;
   fUpThreshold   = 0.0;
   fUpCondition   = "percent";
   fUpSamples     = 100.0;
   fQRatio        = 1.0;
   fLoQ           = 0.05;
   fHiQ           = 0.95;
   fEntropy       = 1.0;
   fNQuantiles    = 10;
   fCallPValue    = 1.0;
   fCallCondition = "percent";
   fCallSamples   = 100.0;

   fHasMAD = fHasCov = fHasVar = fHasDif = fHasM2m = fHasGap = 0;
   fHasLoT = fHasUpT = fHasQua = fHasEnt = fHasCal = kFALSE;
}//Constructor

//______________________________________________________________________________
XPreFilterHeader::~XPreFilterHeader()
{
   // PreFilterHeader destructor
   if(kCS) cout << "---XPreFilterHeader::~XPreFilterHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniFilterHeader                                                     //
//                                                                      //
// Class for storing univariate filter information in tree header       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUniFilterHeader::XUniFilterHeader()
                 :XTreeHeader()
{
   // Default UniFilterHeader constructor
   if(kCS) cout << "---XUniFilterHeader::XUniFilterHeader(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XUniFilterHeader::XUniFilterHeader(const char *str, Int_t treeid)
                 :XTreeHeader(str, treeid)
{
   // Normal UniFilterHeader constructor
   if(kCS) cout << "---XUniFilterHeader::XUniFilterHeader------" << endl;

   fFCValue        = 1.0;
   fFCDirection    = 0;
   fStat           = 1.0;
   fPValue         = 0.02;
   fPChance        = 0.05;
   fPAdjust        = 0.05;
   fCallPValue     = 1.0;
   fCallCondition1 = "percent";
   fCallSamples1   = 100.0;
   fCallCondition2 = "percent";
   fCallSamples2   = 100.0;

   fHasStat = fHasPVal = fHasPCha = fHasPAdj = kFALSE;
   fHasFdCh = fHasUniT = fHasCall = kFALSE;
}//Constructor

//______________________________________________________________________________
XUniFilterHeader::~XUniFilterHeader()
{
   // UniFilterHeader destructor
   if(kCS) cout << "---XUniFilterHeader::~XUniFilterHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiFilterHeader                                                   //
//                                                                      //
// Class for storing multivariate filter information in tree header     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMultiFilterHeader::XMultiFilterHeader()
                   :XTreeHeader()
{
   // Default MultiFilterHeader constructor
   if(kCS) cout << "---XMultiFilterHeader::XMultiFilterHeader(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XMultiFilterHeader::XMultiFilterHeader(const char *str, Int_t treeid)
                   :XTreeHeader(str, treeid)
{
   // Normal MultiFilterHeader constructor
   if(kCS) cout << "---XMultiFilterHeader::XMultiFilterHeader------" << endl;

}//Constructor

//______________________________________________________________________________
XMultiFilterHeader::~XMultiFilterHeader()
{
   // MultiFilterHeader destructor
   if(kCS) cout << "---XMultiFilterHeader::~XMultiFilterHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniTestHeader                                                       //
//                                                                      //
// Class for storing univariate-test information in tree header         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUniTestHeader::XUniTestHeader()
               :XTreeHeader()
{
   // Default UniTestHeader constructor
   if(kCS) cout << "---XUniTestHeader::XUniTestHeader(default)------" << endl;

   fConfLevel   = 0.95;
   fMu          = 0.0;
   fNPerm       = 0;
   fAlternative = "twosided";
   fAdjustment  = "none";
   fTwoSample   = kTRUE;
   fPaired      = kFALSE;
   fEqualVar    = kFALSE;
}//Constructor

//______________________________________________________________________________
XUniTestHeader::XUniTestHeader(const char *str, Int_t treeid)
               :XTreeHeader(str, treeid)
{
   // Normal UniTestHeader constructor
   if(kCS) cout << "---XUniTestHeader::XUniTestHeader------" << endl;

   fConfLevel   = 0.95;
   fMu          = 0.0;
   fNPerm       = 0;
   fAlternative = "twosided";
   fAdjustment  = "none";
   fTwoSample   = kTRUE;
   fPaired      = kFALSE;
   fEqualVar    = kFALSE;
}//Constructor

//______________________________________________________________________________
XUniTestHeader::~XUniTestHeader()
{
   // UniTestHeader destructor
   if(kCS) cout << "---XUniTestHeader::~XUniTestHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiTestHeader                                                     //
//                                                                      //
// Class for storing multivariate-test information in tree header       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMultiTestHeader::XMultiTestHeader()
                 :XTreeHeader()
{
   // Default MultiTestHeader constructor
   if(kCS) cout << "---XMultiTestHeader::XMultiTestHeader(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XMultiTestHeader::XMultiTestHeader(const char *str, Int_t treeid)
                 :XTreeHeader(str, treeid)
{
   // Normal MultiTestHeader constructor
   if(kCS) cout << "---XMultiTestHeader::XMultiTestHeader------" << endl;

}//Constructor

//______________________________________________________________________________
XMultiTestHeader::~XMultiTestHeader()
{
   // MultiTestHeader destructor
   if(kCS) cout << "---XMultiTestHeader::~XMultiTestHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XClusterHeader                                                       //
//                                                                      //
// Class for storing cluster analysis information in tree header        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XClusterHeader::XClusterHeader()
               :XTreeHeader()
{
   // Default ClusterHeader constructor
   if(kCS) cout << "---XClusterHeader::XClusterHeader(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XClusterHeader::XClusterHeader(const char *str, Int_t treeid)
               :XTreeHeader(str, treeid)
{
   // Normal ClusterHeader constructor
   if(kCS) cout << "---XClusterHeader::XClusterHeader------" << endl;

}//Constructor

//______________________________________________________________________________
XClusterHeader::~XClusterHeader()
{
   // ClusterHeader destructor
   if(kCS) cout << "---XClusterHeader::~XClusterHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRegressionHeader                                                    //
//                                                                      //
// Class for storing regression-test information in tree header         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRegressionHeader::XRegressionHeader()
                  :XTreeHeader()
{
   // Default RegressionHeader constructor
   if(kCS) cout << "---XRegressionHeader::XRegressionHeader(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XRegressionHeader::XRegressionHeader(const char *str, Int_t treeid)
                  :XTreeHeader(str, treeid)
{
   // Normal RegressionHeader constructor
   if(kCS) cout << "---XRegressionHeader::XRegressionHeader------" << endl;

}//Constructor

//______________________________________________________________________________
XRegressionHeader::~XRegressionHeader()
{
   // XRegressionHeader destructor
   if(kCS) cout << "---XRegressionHeader::~XRegressionHeader------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalySet                                                            //
//                                                                      //
// Class for analysis of microarray data                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAnalySet::XAnalySet()
          :XProcesSet()
{
   // Default AnalySet constructor
   if(kCS) cout << "---XAnalySet::XAnalySet(default)------" << endl;

   fCalls      = 0;
   fFilters    = 0;
   fFilterTree = 0;
   fMinFilters = 0;
   fLogBase    = "0";
}//Constructor

//______________________________________________________________________________
XAnalySet::XAnalySet(const char *name, const char *type)
          :XProcesSet(name, type)
{
   // Normal AnalySet constructor
   if(kCS) cout << "---XAnalySet::XAnalySet------" << endl;

   fCalls      = 0;   //created on demand only
   fFilters    = 0;   //created on demand only
   fFilterTree = 0;
   fMinFilters = 0;   //ignore filters
   fLogBase    = "0"; //linear
}//Constructor

//______________________________________________________________________________
XAnalySet::~XAnalySet()
{
   // AnalySet destructor
   if(kCS) cout << "---XAnalySet::~XAnalySet------" << endl;

   fFilterTree = 0;

   if(fFilters) {fFilters->Delete(); delete fFilters; fFilters = 0;}
   if(fCalls)   {fCalls->Delete();   delete fCalls;   fCalls   = 0;}

   // delete temporary trees
   TFile *file = 0;
   file = (TFile*)gROOT->FindObject("tmp_exprtrees.root");
   if (file) {file->Close(); delete file; file = 0;}
   file = (TFile*)gROOT->FindObject("tmp_calltrees.root");
   if (file) {file->Close(); delete file; file = 0;}
}//Destructor

//______________________________________________________________________________
Int_t XAnalySet::GetCallMask(Int_t ntree, TTree **tree, Int_t n, Int_t *msk)
{
   // Fill array msk with entries from call tree(s)
   if(kCS) cout << "------XAnalySet::GetCallMask------" << endl;

   return errNoErr;
}//GetCallMask

//______________________________________________________________________________
Int_t XAnalySet::GetFilterMask(Int_t ntree, TTree **tree, Int_t n, Int_t *msk)
{
   // Fill array msk of size n with entries from ntree filter tree(s)
   // The msk flag is set to one only if the flag is set for a number of filter
   // trees greater than on equal to the minimum number minftrs defined by the
   // filter tree id.
   // Note: minftrs=1 is equal to OR; AND: minftrs=ntree is equal to AND
   if(kCS) cout << "------XAnalySet::GetFilterMask------" << endl;

   if ((ntree == 0) || (tree == 0)) {
      for (Int_t i=0; i<n; i++) msk[i] = 1;
      return errNoErr;
   }//if

// Get tree branches for leaves "fFlag"
   TBranch **brch = new TBranch*[ntree];
   TLeaf   **leaf = new TLeaf*[ntree];
   for (Int_t k=0; k<ntree; k++) {
      leaf[k] = tree[k]->FindLeaf("fFlag");
      if (leaf[k] == 0) break;
      brch[k] = leaf[k]->GetBranch();
   }//for_k

   Int_t   minftrs  = 1;
   TString treename = tree[0]->GetName();
   XIdxString  *str = 0;

// Get minimal number of filters to satisfy from filtertree id
   TIter next(fSelections);
   while ((str = (XIdxString*)next())) {
      if (strcmp(treename.Data(),str->GetName()) == 0) {
         minftrs = (str->GetIndex() <= ntree) ? str->GetIndex() : ntree;
         break;
      }//if
   }//while

// Set msk if flag >= minimal number of filters (OR: minftrs=1; AND: minftrs=ntree)
   for (Int_t i=0; i<n; i++) {
      Int_t flag = 0;
      for (Int_t k=0; k<ntree; k++) {
         brch[k]->GetEntry(i);
         flag += (Int_t)leaf[k]->GetValue();
      }//for_k

      msk[i] = (Int_t)(flag >= minftrs);
   }//for_i

   delete [] leaf;
   delete [] brch;

   return errNoErr;
}//GetFilterMask

//______________________________________________________________________________
Int_t XAnalySet::CopyExprTrees(Int_t ntree, TTree **fromtree, TTree **totree,
                 Int_t n, Int_t *msk, Int_t base, Bool_t save)
{
   // Copy expression trees from fromtree to trees totree.
//???   // If new file is created for analysis tree then store new tree(s) in the
//???   // same file, otherwise store tree(s) in temporary file "tmp_trees.root"
   // The copied trees will only contain the subset of tree entries defined by
   // mask array msk.
   // If base>0 then convert tree entries to fLogBase
   // If save=kTRUE the new trees will be stored in temporary file "tmp_exprtrees.root"
//??Problem with leafname is not fLevel in Analyse!! 
   if(kCS) cout << "------XAnalySet::CopyExprTrees------" << endl;

   if ((fromtree == 0) || (msk == 0)) return errNoErr;

   TFile *file = 0;
   if (save) {
      file = new TFile("tmp_exprtrees.root", "RECREATE");
      if (!file || file->IsZombie()) {
         cerr << "Error: Could not create temporary file <tmp_exprtrees.root>."
              << endl;
         SafeDelete(file);
         return errCreateFile;
      }//if
   }//if

   Int_t numneg = 0;
   for (Int_t k=0; k<ntree; k++) {
      TTree *tmptree = new TTree(fromtree[k]->GetName(), fromtree[k]->GetTitle());
      if (tmptree == 0) return errCreateTree;
      Int_t split = 99;
      XExpression *expr = 0;
      expr = new XExpression();
      tmptree->Branch("ExprBranch", "XExpression", &expr, 64000, split);

      XExpression *exprk = 0;
      fromtree[k]->SetBranchAddress("ExprBranch",&exprk);

      Double_t  v;
      if (base == 0) {
         for (Int_t i=0; i<n; i++) {
            // filter mask!
            if (msk[i] == 0) continue;

            fromtree[k]->GetEntry(i);
            expr->SetUnitID(exprk->GetUnitID());
            expr->SetLevel(exprk->GetLevel());
            tmptree->Fill();
         }//for_i
      } else if (base == 2) {
         for (Int_t i=0; i<n; i++) {
            // filter mask!
            if (msk[i] == 0) continue;

            fromtree[k]->GetEntry(i);
            expr->SetUnitID(exprk->GetUnitID());
            v = exprk->GetLevel();
            if (v > 0) {expr->SetLevel(TMath::Log2(v));}
            else       {expr->SetLevel(fNegLog); numneg++;}
            tmptree->Fill();
         }//for_i
      } else if (base == 10) {
         for (Int_t i=0; i<n; i++) {
            // filter mask!
            if (msk[i] == 0) continue;

            fromtree[k]->GetEntry(i);
            expr->SetUnitID(exprk->GetUnitID());
            v = exprk->GetLevel();
            if (v > 0) {expr->SetLevel(TMath::Log10(v));}
            else       {expr->SetLevel(fNegLog); numneg++;}
            tmptree->Fill();
         }//for_i
      } else if (base == 1) {
         for (Int_t i=0; i<n; i++) {
            // filter mask!
            if (msk[i] == 0) continue;

            fromtree[k]->GetEntry(i);
            expr->SetUnitID(exprk->GetUnitID());
            v = exprk->GetLevel();
            if (v > 0) {expr->SetLevel(TMath::Log(v));}
            else       {expr->SetLevel(fNegLog); numneg++;}
            tmptree->Fill();
         }//for_i
      }//if

      totree[k] = tmptree;
      if (save) tmptree->Write();

      SafeDelete(expr);
      tmptree->ResetBranchAddress(tmptree->GetBranch("ExprBranch"));
      SafeDelete(exprk);
      fromtree[k]->ResetBranchAddress(fromtree[k]->GetBranch("ExprBranch"));
   }//for_k

   if (numneg > 0) {
      cout << "Warning: <" << numneg << "> data<=0 replaced with <" << fNegLog
           << ">." << endl;
   }//if

   return errNoErr;
}//CopyExprTrees

/////////////
//PROBLEM
////////////
// need to file->Close() for tmp_exprtrees.root and tmp_calltrees.root
//ev in:  void XAnalysisManager::Close()

//______________________________________________________________________________
Int_t XAnalySet::CopyCallTrees(Int_t ntree, TTree **fromtree, TTree **totree,
                 Int_t n, Int_t *msk, Bool_t save)
{
   // Copy call trees from fromtree to trees totree.
//???   // If new file is created for analysis tree then store new tree(s) in the
//???   // same file, otherwise store tree(s) in temporary file "tmp_trees.root"
   // The copied trees will only contain the subset of tree entries defined by
   // mask array msk.
//Why save in tmp_calltrees???
   // If save=kTRUE the new trees will be stored in temporary file "tmp_calltrees.root"
   if(kCS) cout << "------XAnalySet::CopyCallTrees------" << endl;

   if ((fromtree == 0) || (msk == 0)) return errNoErr;

   TFile *file = 0;
   if (save) {
      file = new TFile("tmp_calltrees.root", "RECREATE");
      if (!file || file->IsZombie()) {
         cerr << "Error: Could not create temporary file <tmp_calltrees.root>."
              << endl;
         SafeDelete(file);
         return errCreateFile;
      }//if
   }//if

   for (Int_t k=0; k<ntree; k++) {
      TTree *tmptree = new TTree(fromtree[k]->GetName(), fromtree[k]->GetTitle());
      if (tmptree == 0) return errCreateTree;
      Int_t split = 99;
      XPCall *call = 0;
      call = new XPCall();
      tmptree->Branch("CallBranch", "XPCall", &call, 64000, split);

      XPCall *callk = 0;
      fromtree[k]->SetBranchAddress("CallBranch", &callk);

      for (Int_t i=0; i<n; i++) {
         // filter mask!
         if (msk[i] == 0) continue;

         fromtree[k]->GetEntry(i);
         call->SetUnitID(callk->GetUnitID());
         call->SetCall(callk->GetCall());
         call->SetPValue(callk->GetPValue());

         tmptree->Fill();
      }//for_i

      totree[k] = tmptree;
      if (save) tmptree->Write("", TObject::kOverwrite);

      SafeDelete(call);
      tmptree->ResetBranchAddress(tmptree->GetBranch("CallBranch"));
      SafeDelete(callk);
      fromtree[k]->ResetBranchAddress(fromtree[k]->GetBranch("CallBranch"));
   }//for_k

   return errNoErr;
}//CopyCallTrees

//______________________________________________________________________________
Bool_t XAnalySet::IsFilterTree(TTree *tree)
{
   // Return kTRUE if tree is filter tree, i.e. has extension kExtenFltr
   if(kCS) cout << "------XAnalySet::IsFilterTree------" << endl;

   TString exten = Path2Name(tree->GetName(),".","");
   if (HasExtension(exten.Data(), kExtenFltr)) return kTRUE;
   return kFALSE;
}//IsFilterTree

//______________________________________________________________________________
Int_t XAnalySet::LogBase()
{
   // Get base of logarithm from fLogBase
   if(kCS) cout << "------XAnalySet::LogBase------" << endl;

   Int_t base = 0; 
   if (strcmp(fLogBase.Data(),"log")        == 0) {base = 1;}
   else if (strcmp(fLogBase.Data(),"log2")  == 0) {base = 2;}
   else if (strcmp(fLogBase.Data(),"log10") == 0) {base = 10;}

   return base;
}//LogBase

//______________________________________________________________________________
Int_t XAnalySet::SelectFilter(TTree *tree, Int_t id)
{
   // Add selected tree to list fFilters
   if(kCS) cout << "------XAnalySet::SelectFilter------" << endl;

   if (!fFilters) fFilters = new TList();

   if (tree != 0) {
   // Add name to list
      XIdxString *str = new XIdxString(id, tree->GetName());
      fFilters->Add(str);
//better:
//??      fFilters->Add(tree);
   }//if

   return errNoErr;
}//SelectFilter

//______________________________________________________________________________
Int_t XAnalySet::SelectCall(TTree *tree, Int_t id)
{
   // Add selected call tree to list fCalls
   if(kCS) cout << "------XAnalySet::SelectCall------" << endl;

   if (!fCalls) fCalls = new TList();

   if (tree != 0) {
   // Add name to list
      XIdxString *str = new XIdxString(id, tree->GetName());
      fCalls->Add(str);
//better:
//??      fCalls->Add(tree);
   }//if

   return errNoErr;
}//SelectCall

//______________________________________________________________________________
Int_t XAnalySet::HandleOption(TTree *tree, Option_t *opt)
{
   // Handle option
   if(kCS) cout << "------XAnalySet::HandleOption------" << endl;

   Int_t   err = errNoErr;
   TString str = opt;
   str.ToLower();

   if (strcmp(str.Data(),"filter") == 0) {
      err = this->SelectFilter(tree);
   } else if (strcmp(str.Data(),"call") == 0) {
      err = this->SelectCall(tree);
   }//if

   return err;
}//HandleOption

//______________________________________________________________________________
Int_t XAnalySet::ExportFilterTrees(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in filter tree to file output
   if(kCS) cout << "------XAnalySet::ExportFilterTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Bool_t hasUnit = kFALSE;
   Bool_t hasFlag = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasUnit = kTRUE;
      hasFlag = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit = kTRUE;}
         if (strcmp(name,"fFlag")     == 0) {hasFlag = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree  **tree = new TTree*[n];
   XMask  **mask = new XMask*[n];
   XUnitID *mskunit = 0;  // get unitID from first tree
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         mask[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         if (k == 0) {
            tree[0]->SetBranchAddress("UnitBranch", &mskunit);
         }//if

         tree[k]->SetBranchAddress("MaskBranch", &mask[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         mask[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         if (k == 0) {
            tree[0]->SetBranchAddress("UnitBranch", &mskunit);
         }//if

         tree[k]->SetBranchAddress("MaskBranch", &mask[k]);
      }//for_k
   }//if

// Get treeinfo and its option for selection of unittree and anntree
   XTreeInfo *info   = (XTreeInfo*)tree[0]->GetUserInfo()->At(0);
   Option_t  *option = info->GetOption();

   Int_t type = eTRANSCRIPT;
   if      (strcmp(option, "exon")     == 0) type = eEXONTYPE; 
   else if (strcmp(option, "probeset") == 0) type = ePROBESET; 

// Find fUnitID from tree[0]
   TLeaf   *idleaf = tree[0]->FindLeaf("fUnitID");
   TBranch *idbrch = (idleaf != 0) ? idleaf->GetBranch() : 0;
   if (idbrch == 0) {
      cout << "Error: Missing UnitBranch for tree <" << tree[0]->GetName() << ">.";
      return errAbort;
   }//if
   Int_t nentries = (Int_t)(tree[0]->GetEntries());

// Get scheme name
   fSchemeName = tree[0]->GetTitle();
   if (strcmp(fSchemeName.Data(), "")  == 0) {
      cerr << "Error: No scheme name is present. Please report error." << endl;
      hasUnit  = 0;
   }//if

// Get chip from scheme file
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XGeneChip *chip = (XGeneChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

// Get unit tree for scheme
   XGCUnit *unit     = 0;
   TTree   *unittree = 0; 
   Int_t numunits    = 0;
   if (hasUnit) {
      unittree = GetUnitTree(chip, type);
      if (unittree == 0) return errGetTree;
      unittree->SetBranchAddress("IdxBranch", &unit);

      numunits = (Int_t)(unittree->GetEntries());
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)    output << sep << "UnitName";
   if (n == 1) {
      if (hasFlag) output << sep << "FLAG";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasFlag)  output << sep << (names[k] + "_FLAG").Data();
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and tree branches
   Int_t unitID = 0;
   Int_t index  = 0;
   for (Int_t i=0; i<nentries; i++) {
      idbrch->GetEntry(i);
      unitID = (Int_t)idleaf->GetValue();
      output << unitID;

      // export unitname
      if (hasUnit) {
         unittree->GetEntry(index++);
         while (unitID != unit->GetUnitID()) {
            if (index == numunits) {
              cerr << "Error: UnitID <" << unitID << "> not found." << endl;
              err = errAbort; goto cleanup;
            }//if

            unittree->GetEntry(index++);
         }//while

         output << sep << unit->GetUnitName();
      }//if

      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // output mask for all trees
         if (hasFlag) output << sep << mask[k]->GetFlag();
      }//for_k
      output << endl;
   }//for_i

//Cleanup
cleanup:
   SafeDelete(schemes);

   // remove trees from RAM
   SafeDelete(mskunit);
   tree[0]->ResetBranchAddress(tree[0]->GetBranch("UnitBranch"));
   for (Int_t k=0; k<n; k++) {
      SafeDelete(mask[k]);
      tree[k]->ResetBranchAddress(tree[k]->GetBranch("MaskBranch"));
      SafeDelete(tree[k]);
   }//for_k

   delete [] mask;
   delete [] tree;

   return err;
}//ExportFilterTrees


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreFilterSet                                                        //
//                                                                      //
// Class for nonspecific filtering of microarray data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPreFilterSet::XPreFilterSet()
              :XAnalySet()
{
   // Default FilterSet constructor
   if(kCS) cout << "---XPreFilterSet::XPreFilterSet(default)------" << endl;

   fFilter = 0;
}//Constructor

//______________________________________________________________________________
XPreFilterSet::XPreFilterSet(const char *name, const char *type)
              :XAnalySet(name, type)
{
   // Normal FilterSet constructor
   if(kCS) cout << "---XPreFilterSet::XPreFilterSet------" << endl;

   fFilter = 0;
}//Constructor

//______________________________________________________________________________
XPreFilterSet::~XPreFilterSet()
{
   // FilterSet destructor
   if(kCS) cout << "---XPreFilterSet::~XPreFilterSet------" << endl;

   fFilter = 0;  //deleted in XAnalySetting
}//Destructor

//______________________________________________________________________________
void XPreFilterSet::AddTreeHeader(const char *treename, Int_t treeid)
{
   // Add tree header with treename to list of trees
   // For selection of tree friends set: treeid > 0 
   if(kCS) cout << "------XPreFilterSet::AddTreeHeader------" << endl;

   if (treeid != 0) {
      Select(treename, treeid);
   } else {
      XPreFilterHeader *header = new XPreFilterHeader(treename, treeid);
      header->SetInfile(fInfile);
      header->SetType(fFilter->GetName());
      header->SetMedAbsDev(fFilter->GetMedAbsDev());
      header->SetCov2Mean(fFilter->GetCov2Mean());
      header->SetVar2Mean(fFilter->GetVar2Mean());
      header->SetDif2Mean(fFilter->GetDif2Mean());
      header->SetMax2Min(fFilter->GetMax2Min());
      header->SetGap2Mean(fFilter->GetGap2Mean());
      header->SetTrimValue(fFilter->GetTrimValue());
      header->SetWindowSize(fFilter->GetWindowSize());
      header->SetLowerCondition(fFilter->GetLowerCondition());
      header->SetLowerThreshold(fFilter->GetLowerThreshold());
      header->SetLowerSamples(fFilter->GetLowerSamples());
      header->SetUpperCondition(fFilter->GetUpperCondition());
      header->SetUpperThreshold(fFilter->GetUpperThreshold());
      header->SetUpperSamples(fFilter->GetUpperSamples());
      header->SetLowQuantile(fFilter->GetLowQuantile());
      header->SetHighQuantile(fFilter->GetHighQuantile());
      header->SetQuantileRatio(fFilter->GetQuantileRatio());
      header->SetEntropy(fFilter->GetEntropy());
      header->SetNumberQuantiles(fFilter->GetNumberQuantiles());
      header->SetCallCondition(fFilter->GetCallCondition());
      header->SetCallPValue(fFilter->GetCallPValue());
      header->SetCallSamples(fFilter->GetCallSamples());
      header->SetMADFilter(fFilter->HasMADFilter());
      header->SetCoefOfVarFilter(fFilter->HasCoefOfVarFilter());
      header->SetVarianceFilter(fFilter->HasVarianceFilter());
      header->SetDifferenceFilter(fFilter->HasDifferenceFilter());
      header->SetRatioFilter(fFilter->HasRatioFilter());
      header->SetGapFilter(fFilter->HasGapFilter());
      header->SetLoThresholdFilter(fFilter->HasLoThresholdFilter());
      header->SetUpThresholdFilter(fFilter->HasUpThresholdFilter());
      header->SetQuantileFilter(fFilter->HasQuantileFilter());
      header->SetEntropyFilter(fFilter->HasEntropyFilter());
      header->SetCallFilter(fFilter->HasCallFilter());

      fHeaders->Add(header);
   }//if
}//AddTreeHeader

//______________________________________________________________________________
Int_t XPreFilterSet::Analyse(const char *infile, const char *outfile, 
                     const char *varlist, Int_t nrows, const char *sepi,
                     const char *sepo, char delim)
{
   // Analyse data
   if(kCS) cout << "------XPreFilterSet::Analyse------" << endl;

   return fFilter->Calculate(infile, outfile, varlist, nrows, sepi, sepo, delim);
}//Analyse

//______________________________________________________________________________
Int_t XPreFilterSet::Analyse(const char *leafname, const char *outtree,
                     const char *varlist)
{
   // Calculate data
   if(kCS) cout << "------XPreFilterSet::Analyse------" << endl;

   Int_t err = errNoErr;

// Get logbase
   Int_t base = LogBase(); 

// Get number of trees with leaf leafname and number of calltrees
   XTreeInfo *info = 0;
   TString tname0 = "";
   TString option = "";
   Int_t numcall  = 0;
   Int_t numfltr  = 0;
   Int_t numleaf  = 0;
   Int_t nentries = 0;
   Int_t numtrees = fTrees->GetSize();
   for (Int_t k=0; k<numtrees; k++) {
      TTree* tree = (TTree*)fTrees->At(k);
      Int_t  size = (Int_t)(tree->GetEntries());

      // count trees (check first for CallBranch to avoid leafname=fCall or fPValue)
      if (tree->GetBranch("CallBranch") != 0) numcall++;
      else if (IsFilterTree(tree)   == kTRUE) numfltr++;
      else if (tree->FindLeaf(leafname) != 0) numleaf++;

      // check for equal number of tree entries
      if (k == 0) {
         nentries = size;
         tname0   = tree->GetName();
         info     = (XTreeInfo*)tree->GetUserInfo()->At(0);
         option   = info->GetOption();
      } else if (size != nentries) {
         return fManager->HandleError(errEQTreeEntries, tname0, tree->GetName());
      }//if
   }//for_k

// Initialize calltrees and leaftrees
   TTree **calltree = new TTree*[numcall];
   TTree **fltrtree = new TTree*[numfltr];
   TTree **leaftree = new TTree*[numleaf];
   for (Int_t k=0; k<numcall; k++) calltree[k] = 0;
   for (Int_t k=0; k<numfltr; k++) fltrtree[k] = 0;
   for (Int_t k=0; k<numleaf; k++) leaftree[k] = 0;

   numcall = 0;
   numfltr = 0;
   numleaf = 0;
   for (Int_t k=0; k<numtrees; k++) {
      TTree* tree = (TTree*)fTrees->At(k);

      // check first for CallBranch to avoid leafname=fCall or fPValue
      if (tree->GetBranch("CallBranch") != 0) {
         calltree[numcall++] = tree;
      } else if (IsFilterTree(tree) == kTRUE) {
         fltrtree[numfltr++] = tree;
      } else if (tree->FindLeaf(leafname) != 0) {
         leaftree[numleaf++] = tree;
      }//if
   }//for_k

// Check for presence of trees
   if ((numleaf == 0) && (numcall == 0)) {
      cerr << "Error: no trees have been selected!" << endl;
      return errAbort;;
   }//if

// Filter trees and/or convert trees to logbase
// (need to set ftr=1 if convert trees to logbase only w/o filtering)
   if ((numfltr > 0) || (base > 0)) {
      Int_t *ftr = 0;
      if (!(ftr = new Int_t[nentries])) {return errInitMemory;}

      err = GetFilterMask(numfltr, fltrtree, nentries, ftr);

      if ((err == errNoErr) && (numleaf > 0)) {
         err = CopyExprTrees(numleaf, leaftree, leaftree, nentries, ftr, base);
      }//if

      if ((err == errNoErr) && (numcall > 0)) {
         err = CopyCallTrees(numcall, calltree, calltree, nentries, ftr);
      }//if

      if (err != errNoErr) {
         cerr << "Error: Could not copy trees, aborting analysis." << endl;
         delete [] ftr; ftr = 0;
         return errGetTree;
      }//if

      if (ftr) {delete [] ftr; ftr = 0;}
   }//if

// Change directory
   if (!fFile->cd(fName)) return errGetDir;

// Create filter tree
   TString tname = TString(outtree) + "." + fFilter->GetTitle();
   // get scheme name from tree
   TString title = leaftree[0] ? leaftree[0]->GetTitle() : calltree[0]->GetTitle();
   TTree  *atree = new TTree(tname, title);
   if (atree == 0) return errCreateTree;

// Apply filter
   if (numleaf > 0 && numcall > 0) {
      err = fFilter->CallFlag(numcall, calltree, varlist, 0);
      err = fFilter->Calculate(numleaf, leaftree, leafname, atree, varlist);
   } else if (numleaf > 0 && numcall == 0) {
      err = fFilter->Calculate(numleaf, leaftree, leafname, atree, varlist);
   } else if (numcall > 0 && numleaf == 0) {
      err = fFilter->CallFlag(numcall, calltree, varlist, atree);
   } else {
   }//if
   if (err) return err;

// Add tree info to tree
   AddTreeInfo(atree, atree->GetName(), option);

// Add expression tree header to list
   this->AddTreeHeader(tname, 0);

// Write expression tree to file 
   err = WriteTree(atree, 0); 
//ev kOverwrite to overwrite old atree with same name!! or keep old atree;cycle?
//better??   err = WriteTree(atree, TObject::kOverwrite);

// Delete created filter tree from RAM
   atree->Delete("");
   atree = 0;

   delete [] leaftree;
   delete [] fltrtree;
   delete [] calltree;

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XPreFilterSet::ExportTreeType(const char *exten, Int_t n, TString *names, 
                     const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output
   if(kCS) cout << "------XPreFilterSet::ExportTreeType------" << endl;

   if (HasExtension(exten, kExtenFltr)) {
      return ExportFilterTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XPreFilterSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                     const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output as XML-file
   if(kCS) cout << "------XPreFilterSet::ExportTreeXML------" << endl;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XPreFilterSet::Initialize(TFile *file, XSetting *setting, const char *infile,
                     const char *treename)
{
   // Initialize filter
   if(kCS) cout << "------XPreFilterSet::Initialize------" << endl;

   if ((setting == 0) || ((file == 0) && (strcmp(infile, "") == 0))) return errAbort;

   fFile     = file;
   fInfile   = infile;
   fTreeName = treename;

   fSetting = (XAnalySetting*)setting;
   if (!fSetting) return errAbort;

// Get files from settings
   fDataFile   = ((XAnalySetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XAnalySetting*)fSetting)->GetSchemeName();
   fSchemeType = ((XAnalySetting*)fSetting)->GetSchemeType();

// Get filter from settings
   fFilter = (XPreFilter*)((XAnalySetting*)fSetting)->GetFilter();
//for export   if (!fFilter) return errAbort;

   fLogBase = ((XAnalySetting*)fSetting)->GetLogBase();
   fNegLog  = ((XAnalySetting*)fSetting)->GetNegLog();

   return errNoErr;
}//Initialize


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnivarSet                                                           //
//                                                                      //
// Class for univariate analysis of microarray data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUnivarSet::XUnivarSet()
           :XAnalySet()
{
   // Default UnivarSet constructor
   if(kCS) cout << "---XUnivarSet::XUnivarSet(default)------" << endl;

   fFilter   = 0;
   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XUnivarSet::XUnivarSet(const char *name, const char *type)
           :XAnalySet(name, type)
{
   // Normal UnivarSet constructor
   if(kCS) cout << "---XUnivarSet::XUnivarSet------" << endl;

   fFilter   = 0;
   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XUnivarSet::~XUnivarSet()
{
   // UnivarSet destructor
   if(kCS) cout << "---XUnivarSet::~XUnivarSet------" << endl;

   fFilter   = 0;  //deleted in ~XAnalySetting
   fAnalyser = 0;  //deleted in ~XAnalySetting
}//Destructor

//______________________________________________________________________________
void XUnivarSet::AddTreeHeader(const char *treename, Int_t treeid)
{
   // Add tree header with treename to list of trees
   // For selection of tree friends set: treeid > 0 
   if(kCS) cout << "------XUnivarSet::AddTreeHeader------" << endl;

   if (treeid > 0) {
      Select(treename, treeid);
   } else {
      this->AddHeader(treename, treeid);
   }//if
}//AddTreeHeader

//______________________________________________________________________________
Int_t XUnivarSet::Analyse(const char *infile, const char *outfile, 
                  const char *varlist, Int_t nrows, const char *sepi,
                  const char *sepo, char delim)
{
   // Analyse data
   if(kCS) cout << "------XUnivarSet::Analyse------" << endl;

   Int_t err = errNoErr;

// Calculate univariate test
   if (fAnalyser) {
      err = fAnalyser->Analyse(infile, outfile, varlist, nrows, -1, 0, sepi, sepo, delim);
      if (err != errNoErr) return err;
   }//if

// Calculate filter to apply for univariate test
   if (fFilter) {
      err = fFilter->Calculate(infile, outfile, varlist, nrows, sepi, sepo, delim);
   }//if

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XUnivarSet::Analyse(const char *leafname, const char *outtree,
                  const char *varlist)
{
   // Analyse data
   if(kCS) cout << "------XUnivarSet::Analyse------" << endl;

   Int_t err = errNoErr;

// Decompose varlist to test for presence of analysis and/or filter types
   Int_t numUTst = 0;
   Int_t numUFtr = 0;

   char *vname = new char[strlen(varlist) + 1];
   char *dname = vname;
   vname = strtok(strcpy(vname,varlist),":");
   while(vname) {
      if (strcmp(vname, kTypeUTst[0]) == 0) {numUTst++;}
      if (strcmp(vname, kTypeUTst[1]) == 0) {numUTst++;}
      if (strcmp(vname, kTypeUTst[2]) == 0) {numUTst++;}
      if (strcmp(vname, kTypeUTst[3]) == 0) {numUTst++;}
      if (strcmp(vname, kUniFltr[0])  == 0) {numUFtr++;}
      if (strcmp(vname, kUniFltr[1])  == 0) {numUFtr++;}
      if (strcmp(vname, kUniFltr[2])  == 0) {numUFtr++;}
      if (strcmp(vname, kUniFltr[3])  == 0) {numUFtr++;}
      if (strcmp(vname, kUniFltr[4])  == 0) {numUFtr++;}
      if (strcmp(vname, kUniFltr[5])  == 0) {numUFtr++;}
      vname = strtok(NULL, ":");
      if (vname == 0) break;
   }//while
   delete [] dname;

   // analysis/filter types are set by calling corresponding InitAlgorithm()
//   if (!numUTst && !numUFtr) numUTst = numUFtr = 1;
   if (fAnalyser && !numUTst) numUTst = 1;
   if (fFilter   && !numUFtr) numUFtr = 1;

// Get logbase
   Int_t base = LogBase(); 

// Get number of trees with leaf leafname and number of filtertrees and calltrees
   XTreeInfo *info = 0;
   TString tname0 = "";
   TString option = "";
   Int_t numcall  = 0;
   Int_t numfltr  = 0;
   Int_t numleaf  = 0;
   Int_t nentries = 0;
   Int_t numtrees = fTrees->GetSize();
   for (Int_t k=0; k<numtrees; k++) {
      TTree* tree = (TTree*)fTrees->At(k);
      Int_t  size = (Int_t)(tree->GetEntries());

      // count trees (check first for CallBranch to avoid leafname=fCall or fPValue)
      if (tree->GetBranch("CallBranch") != 0) numcall++;
      else if (IsFilterTree(tree)   == kTRUE) numfltr++;
      else if (tree->FindLeaf(leafname) != 0) numleaf++;

      // check for equal number of tree entries
      if (k == 0) {
         nentries = size;
         tname0   = tree->GetName();
         info     = (XTreeInfo*)tree->GetUserInfo()->At(0);
         option   = info->GetOption();
      } else if (size != nentries) {
         return fManager->HandleError(errEQTreeEntries, tname0, tree->GetName());
      }//if
   }//for_k

// Initialize calltrees and leaftrees
   TTree **calltree = new TTree*[numcall];
   TTree **fltrtree = new TTree*[numfltr];
   TTree **leaftree = new TTree*[numleaf];
   for (Int_t k=0; k<numcall; k++) calltree[k] = 0;
   for (Int_t k=0; k<numfltr; k++) fltrtree[k] = 0;
   for (Int_t k=0; k<numleaf; k++) leaftree[k] = 0;

   numcall = 0;
   numfltr = 0;
   numleaf = 0;
   for (Int_t k=0; k<numtrees; k++) {
      TTree* tree = (TTree*)fTrees->At(k);

      // check first for CallBranch to avoid leafname=fCall or fPValue
      if (tree->GetBranch("CallBranch") != 0) {
         calltree[numcall++] = tree;
      } else if (IsFilterTree(tree) == kTRUE) {
         fltrtree[numfltr++] = tree;
      } else if (tree->FindLeaf(leafname) != 0) {
         leaftree[numleaf++] = tree;
      }//if
   }//for_k

// Check for presence of trees
//   if ((numleaf == 0) && (numcall == 0)) {
   if (numleaf == 0) {
      cerr << "Error: no expression trees have been selected!" << endl;
      return errAbort;;
   }//if

// Filter trees and/or convert trees to logbase
// (need to set ftr[i]=1 if convert trees to logbase only w/o filtering)
   if ((numfltr > 0) || (base > 0)) {
      Int_t *ftr = 0;
      if (!(ftr = new Int_t[nentries])) {return errInitMemory;}

      err = GetFilterMask(numfltr, fltrtree, nentries, ftr);

      if ((err == errNoErr) && (numleaf > 0)) {
         err = CopyExprTrees(numleaf, leaftree, leaftree, nentries, ftr, base);
      }//if

//      if ((err == errNoErr) && (numcall > 0)) {
      if ((err == errNoErr) && (numcall > 0) && (numfltr > 0)) {
         err = CopyCallTrees(numcall, calltree, calltree, nentries, ftr);
      }//if

      if (err != errNoErr) {
         cerr << "Error: Could not copy trees, aborting analysis." << endl;
         delete [] ftr; ftr = 0;
         return errGetTree;
      }//if

      if (ftr) {delete [] ftr; ftr = 0;}
   }//if

// Change directory
   if (!fFile->cd(fName)) return errGetDir;

//EV??? decompose leafname
//?? check leafname for fAnalyser is "fLevel"
//?? check leafname for fFilter is "????"
// Calculate univariate test and/or filter for univariate test
   if (fAnalyser && numUTst && numleaf) {
      err = this->UniTest(numleaf, leaftree, leafname, outtree, varlist, option);
   }//if
//??   if (fFilter && numUFtr && numleaf && (err == errNoErr)) {
   if (fFilter && numUFtr && (err == errNoErr)) {
      err = this->Filter(numleaf, leaftree, leafname, numcall, calltree,
                         outtree, varlist, option, base);
   }//if

// Delete fTree from RAM
   if (fTree) {fTree->Delete(""); fTree = 0;}

   delete [] leaftree;
   delete [] fltrtree;
   delete [] calltree;

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XUnivarSet::ExportTreeType(const char *exten, Int_t n, TString *names, 
                  const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output
   if(kCS) cout << "------XUnivarSet::ExportTreeType------" << endl;

   if (HasExtension(exten, kExtenUTst)) {
      return this->ExportUnivarTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kExtenFltr)) {
      return this->ExportFilterTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XUnivarSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                  const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output as XML-file
   if(kCS) cout << "------XUnivarSet::ExportTreeXML------" << endl;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XUnivarSet::Initialize(TFile *file, XSetting *setting, const char *infile,
                  const char *treename)
{
   // Initialize algorithms
   if(kCS) cout << "------XUnivarSet::Initialize------" << endl;

   if ((file == 0) || (setting == 0)) return errAbort;

   fFile     = file;
   fInfile   = infile;
   fTreeName = treename;

   fSetting = (XAnalySetting*)setting;
   if (!fSetting) return errAbort;

// Get files from settings
   fDataFile   = ((XAnalySetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XAnalySetting*)fSetting)->GetSchemeName();
   fSchemeType = ((XAnalySetting*)fSetting)->GetSchemeType();

// Get algorithms from settings
   fAnalyser = (XUniTester*)((XAnalySetting*)fSetting)->GetAnalyser();

   // optional filter (can be: filter = 0) 
   fFilter = (XUniFilter*)((XAnalySetting*)fSetting)->GetFilter();

// not needed for Export!!!
//   if (!(fAnalyser || fFilter)) return errAbort;

   fLogBase = ((XAnalySetting*)fSetting)->GetLogBase();
   fNegLog  = ((XAnalySetting*)fSetting)->GetNegLog();

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XUnivarSet::ExportUnivarTrees(Int_t n, TString *names, const char *varlist,
                  ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   // if varlist has hasFC then fold-change = Mn2/Mn1 is exported
   if(kCS) cout << "------XUnivarSet::ExportUnivarTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Short_t hasUnit   = 0;  //unit name
   Short_t hasTrans  = 0;  //transcript_id
   Short_t hasName   = 0;  //gene name
   Short_t hasSymbol = 0;  //gene symbol
   Short_t hasAccess = 0;  //mRNA accession
   Short_t hasEntrez = 0;  //entrez ID
   Short_t hasChromo = 0;  //chromosome
   Short_t hasStart  = 0;  //start position
   Short_t hasStop   = 0;  //stop position
   Short_t hasStrand = 0;  //strand
   Short_t hasCyto   = 0;  //cytoband

   Short_t hasAnnot  = 0;  //annotation
   Short_t hasData   = 0;  //data

   Bool_t hasStat = kFALSE;
   Bool_t hasMn1  = kFALSE;
   Bool_t hasMn2  = kFALSE;
   Bool_t hasSE   = kFALSE;
   Bool_t hasDF   = kFALSE;
   Bool_t hasPVal = kFALSE;
   Bool_t hasNPer = kFALSE;
   Bool_t hasPCha = kFALSE;
   Bool_t hasPAdj = kFALSE;
   Bool_t hasFC   = kFALSE;

   Bool_t hasFlag = kFALSE;
   Bool_t hasMask = kFALSE;

   Short_t idx = 0;
   if (strcmp(varlist,"*") == 0) {
      hasUnit   = ++idx;
      hasTrans  = ++idx;
      hasName   = ++idx;
      hasSymbol = ++idx;
      hasAccess = ++idx;
      hasEntrez = ++idx;
      hasChromo = ++idx;
      hasStart  = ++idx;
      hasStop   = ++idx;
      hasStrand = ++idx;
      hasCyto   = ++idx;

      hasStat = kTRUE;
      hasMn1  = kTRUE;
      hasMn2  = kTRUE;
      hasSE   = kTRUE;
      hasDF   = kTRUE;
      hasPVal = kTRUE;
      hasNPer = kTRUE;
      hasPCha = kTRUE;
      hasPAdj = kTRUE;
      hasFC   = kTRUE;
      hasFlag = kTRUE;
   } else {
      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;

      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname,"fUnitName")     == 0) {hasUnit   = ++idx;}
         if (strcmp(vname,"fTranscriptID") == 0) {hasTrans  = ++idx;}
         if (strcmp(vname,"fName")         == 0) {hasName   = ++idx;}
         if (strcmp(vname,"fSymbol")       == 0) {hasSymbol = ++idx;}
         if (strcmp(vname,"fAccession")    == 0) {hasAccess = ++idx;}
         if (strcmp(vname,"fEntrezID")     == 0) {hasEntrez = ++idx;}
         if (strcmp(vname,"fChromosome")   == 0) {hasChromo = ++idx;}
         if (strcmp(vname,"fStart")        == 0) {hasStart  = ++idx;}
         if (strcmp(vname,"fStop")         == 0) {hasStop   = ++idx;}
         if (strcmp(vname,"fStrand")       == 0) {hasStrand = ++idx;}
         if (strcmp(vname,"fCytoBand")     == 0) {hasCyto   = ++idx;}

         if (strcmp(vname,"stat") == 0) {hasStat = kTRUE;}
         if (strcmp(vname,"mn1")  == 0) {hasMn1  = kTRUE;}
         if (strcmp(vname,"mn2")  == 0) {hasMn2  = kTRUE;}
         if (strcmp(vname,"se")   == 0) {hasSE   = kTRUE;}
         if (strcmp(vname,"df")   == 0) {hasDF   = kTRUE;}
         if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE;}
         if (strcmp(vname,"nper") == 0) {hasNPer = kTRUE;}
         if (strcmp(vname,"pcha") == 0) {hasPCha = kTRUE;}
         if (strcmp(vname,"padj") == 0) {hasPAdj = kTRUE;}
         if (strcmp(vname,"fc")   == 0) {hasFC   = kTRUE;}

         if (strcmp(vname,"flag") == 0) {hasFlag = kTRUE;}
         if (strcmp(vname,"mask") == 0) {hasMask = kTRUE;}

         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while

      delete [] dname;
   }//if

   // check for presence of at least one annotation variable
   hasAnnot = (hasTrans  + hasName  + hasSymbol + hasAccess + hasEntrez
            +  hasChromo + hasStart + hasStop   + hasStrand + hasCyto);
   // check for presence of at least one of data
   hasData  = (hasStat || hasMn1  || hasMn2  || hasSE   || hasDF
            || hasPVal || hasNPer || hasPCha || hasPAdj || hasFC);
   hasData  = (hasData > 0) ? ++idx : 0;

///////////
//PROBLEM: hasAnnot but not hasUnit => crash!!!
// TO DO: for annotation necessary to init unittree
// if (hasAnnot & !hasUnit) {hasUnit = ++idx;}
//////////

// Get logbase: necessary to convert mn1 and mn2 and get correct fc = mn2/mn1
   XFolder    *data = 0;
   XUnivarSet *set  = 0;
   data = (XFolder*)(fFile->Get(kContent));
   if (data && (set = (XUnivarSet*)data->FindObject(GetName(), "XUnivarSet"))) {
      fLogBase = set->GetLogBase();
      fNegLog  = set->GetNegLog();
   }//if
   Int_t base = LogBase(); 

// Get trees
   TTree **tree = new TTree*[n];
   if (fTree) {
      tree[0] = fTree;
   } else if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;
      }//for_k
   }//if

// Get treeinfo and its option for selection of unittree and anntree
   XTreeInfo *info   = (XTreeInfo*)tree[0]->GetUserInfo()->At(0);
   Option_t  *option = info->GetOption();

   Int_t type = eTRANSCRIPT;
   if      (strcmp(option, "exon")     == 0) type = eEXONTYPE; 
   else if (strcmp(option, "probeset") == 0) type = ePROBESET; 

// Find fUnitID from tree[0]
   TLeaf   *idleaf = tree[0]->FindLeaf("fUnitID");
   TBranch *idbrch = (idleaf != 0) ? idleaf->GetBranch() : 0;
   if (idbrch == 0) {
      cout << "Error: Missing UnitBranch for tree <" << tree[0]->GetName() << ">.";
      return errAbort;
   }//if

// Init flag
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   Int_t *flag = 0;
   if (!(flag = new Int_t[nentries])) {return errInitMemory;}
   for (Int_t i=0; i<nentries; i++) flag[i] = 1;

// Check for presence of unifilter tree
   if (hasFlag || hasMask) {
      TTree *masktree = 0;
      TString ufrname = Path2Name(tree[0]->GetName(),"",".") + "." + kExtenFltr[1];
      masktree = (TTree*)gDirectory->Get(ufrname.Data());
      if (masktree && (masktree->GetEntries() == nentries)) {
         XMask *mask = new XMask();
         masktree->SetBranchAddress("MaskBranch", &mask);
         for (Int_t i=0; i<nentries; i++) {
            masktree->GetEntry(i);
            flag[i] = mask->GetFlag();
         }//for_i
         delete mask;
      } else {
         cout << "Warning: tree <" << ufrname.Data() << "> does not exist or has not <"
              << nentries << "> entries.";
      }//if
   }//if

// Check for existence of tree leafs
   Bool_t hasLeafStat = kFALSE;
   Bool_t hasLeafMn1  = kFALSE;
   Bool_t hasLeafMn2  = kFALSE;
   Bool_t hasLeafSE   = kFALSE;
   Bool_t hasLeafDF   = kFALSE;
   Bool_t hasLeafPVal = kFALSE;
   Bool_t hasLeafNPer = kFALSE;
   Bool_t hasLeafPCha = kFALSE;
   Bool_t hasLeafPAdj = kFALSE;

   TLeaf *leaf = 0;
   if ((leaf = tree[0]->FindLeaf("stat")) != 0) {hasLeafStat = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("mn1"))  != 0) {hasLeafMn1  = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("mn2"))  != 0) {hasLeafMn2  = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("se"))   != 0) {hasLeafSE   = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("df"))   != 0) {hasLeafDF   = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("pval")) != 0) {hasLeafPVal = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("nper")) != 0) {hasLeafNPer = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("pcha")) != 0) {hasLeafPCha = kTRUE;}
   if ((leaf = tree[0]->FindLeaf("padj")) != 0) {hasLeafPAdj = kTRUE;}

   hasStat = hasStat ? hasLeafStat : kFALSE;
   hasMn1  = hasMn1  ? hasLeafMn1  : kFALSE;
   hasMn2  = hasMn2  ? hasLeafMn2  : kFALSE;
   hasSE   = hasSE   ? hasLeafSE   : kFALSE;
   hasDF   = hasDF   ? hasLeafDF   : kFALSE;
   hasPVal = hasPVal ? hasLeafPVal : kFALSE;
   hasNPer = hasNPer ? hasLeafNPer : kFALSE;
   hasPCha = hasPCha ? hasLeafPCha : kFALSE;
   hasPAdj = hasPAdj ? hasLeafPAdj : kFALSE;
   hasFC   = hasFC   ? (hasLeafMn1 && hasLeafMn2) : kFALSE;

// Get scheme name
   fSchemeName = tree[0]->GetTitle();
   if (strcmp(fSchemeName.Data(), "")  == 0) {
      cerr << "Error: No scheme name is present. Please report error." << endl;
      hasUnit  = 0;
      hasAnnot = 0;
   }//if

// Get chip from scheme file
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XGeneChip *chip = (XGeneChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

// Get unit tree for scheme
   XGCUnit *unit     = 0;
   TTree   *unittree = 0; 
   Int_t    numunits = 0;
   if (hasUnit) {
      unittree = GetUnitTree(chip, type);
      if (unittree == 0) return errGetTree;
      unittree->SetBranchAddress("IdxBranch", &unit);

      numunits = (Int_t)(unittree->GetEntries());
   }//if

// Get annotation tree for scheme
   XTransAnnotation *annot = 0;
   TTree          *anntree = 0; 
   if (hasAnnot) {
      anntree = GetAnnotationTree(chip, type);
      if (anntree) {
         anntree->SetBranchAddress("AnnBranch", &annot);
      } else {
         cout << "Warning: Missing annotation, gene info not exported." << endl;
         hasAnnot = 0;
      }//if
   }//if

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (hasAnnot) {
      Int_t numannot = (Int_t)(anntree->GetEntries());
      if (!(htable = new THashTable(2*numannot))) return errInitMemory;
      htable = FillHashTable(htable, anntree, annot, type);
   }//if

// Output header
   output << "UNIT_ID";
   for (Short_t j=1; j<=idx; j++) {
      if (hasUnit == j)      output << sep << "UnitName";

      if (hasAnnot > 0) {
         if (hasTrans  == j) output << sep << "TranscriptID";
         if (hasName   == j) output << sep << "GeneName";
         if (hasSymbol == j) output << sep << "GeneSymbol";
         if (hasAccess == j) output << sep << "GeneAccession";
         if (hasEntrez == j) output << sep << "EntrezID";
         if (hasChromo == j) output << sep << "Chromosome";
         if (hasStart  == j) output << sep << "Start";
         if (hasStop   == j) output << sep << "Stop";
         if (hasStrand == j) output << sep << "Strand";
         if (hasCyto   == j) output << sep << "Cytoband";
      }//if

      if (hasData == j) {
         if (n > 1) {
            for (Int_t k=0; k<n; k++) {
               if (hasStat) output << sep << (names[k] + "_Statistics").Data();
               if (hasMn1)  output << sep << (names[k] + "_Mean1").Data();
               if (hasMn2)  output << sep << (names[k] + "_Mean2").Data();
               if (hasSE)   output << sep << (names[k] + "_StandardError").Data();
               if (hasDF)   output << sep << (names[k] + "_DegreeOfFreedom").Data();
               if (hasPVal) output << sep << (names[k] + "_P-Value").Data();
               if (hasNPer) output << sep << (names[k] + "_NumberPermutations").Data();
               if (hasPCha) output << sep << (names[k] + "_P-Chance").Data();
               if (hasPAdj) output << sep << (names[k] + "_P-Adjusted").Data();
               if (hasFC)   output << sep << (names[k] + "_FoldChange").Data();
               if (hasFlag) output << sep << (names[k] + "_Flag").Data();
            }//for_i
         } else {
            if (hasStat) output << sep << "Statistics";
            if (hasMn1)  output << sep << "Mean1";
            if (hasMn2)  output << sep << "Mean2";
            if (hasSE)   output << sep << "StandardError";
            if (hasDF)   output << sep << "DegreeOfFreedom";
            if (hasPVal) output << sep << "P-Value";
            if (hasNPer) output << sep << "NumberPermutations";
            if (hasPCha) output << sep << "P-Chance";
            if (hasPAdj) output << sep << "P-Adjusted";
            if (hasFC)   output << sep << "FoldChange";
            if (hasFlag) output << sep << "Flag";
         }//if
      }//if
   }//for_j
   output << endl;

// Loop over tree entries and tree branches
   XIdxString *idxstr = 0;
   TBranch *brch = 0;
   Int_t  unitID = 0;
   Double_t mn1, mn2, fc;
   mn1 = mn2 = fc = 0.0;

   Int_t index = 0;
   for (Int_t i=0; i<nentries; i++) {
      if (hasMask && (flag[i] == 0)) continue;

      idbrch->GetEntry(i);
      unitID = (Int_t)idleaf->GetValue();
      output << unitID;

      if (hasUnit) {
         unittree->GetEntry(index++);
         while (unitID != unit->GetUnitID()) {
            if (index == numunits) {
              cerr << "Error: UnitID <" << unitID << "> not found." << endl;
              err = errAbort; goto cleanup;
            }//if

            unittree->GetEntry(index++);
         }//while
      }//if

      for (Short_t j=1; j<=idx; j++) {
         // export unitname
         if (hasUnit == j) {
            output << sep << unit->GetUnitName();
         }//if

         // export annotation
         if (hasAnnot > 0) {
            idxstr = (XIdxString*)(htable->FindObject(unit->GetUnitName()));
            if (idxstr) {
               anntree->GetEntry(idxstr->GetIndex());
               if (hasTrans  == j) output << sep << GetTranscriptID(unit, annot, type);
               if (hasName   == j) output << sep << annot->GetName();
               if (hasSymbol == j) output << sep << annot->GetSymbol();
               if (hasAccess == j) output << sep << annot->GetAccession();
               if (hasEntrez == j) output << sep << annot->GetEntrezID();
               if (hasChromo == j) output << sep << annot->GetChromosome();
               if (hasStart  == j) output << sep << annot->GetStart();
               if (hasStop   == j) output << sep << annot->GetStop();
               if (hasStrand == j) output << sep << annot->GetStrand();
               if (hasCyto   == j) output << sep << annot->GetCytoBand();
            } else {
               if (hasTrans  == j) output << sep << "NA";
               if (hasName   == j) output << sep << "NA";
               if (hasSymbol == j) output << sep << "NA";
               if (hasAccess == j) output << sep << "NA";
               if (hasEntrez == j) output << sep << "NA";
               if (hasChromo == j) output << sep << "NA";
               if (hasStart  == j) output << sep << "-1";
               if (hasStop   == j) output << sep << "-1";
               if (hasStrand == j) output << sep << "?";
               if (hasCyto   == j) output << sep << "NA";
            }//if
         }//if

         // export data
         if (hasData == j) {
            for (Int_t k=0; k<n; k++) {
               if (hasStat) {
                  leaf = tree[k]->FindLeaf("stat");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasMn1 || hasFC) {
                  leaf = tree[k]->FindLeaf("mn1");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  mn1 = leaf->GetValue();
                  mn1 = (base > 0) ? ((base > 1) ? ((base > 2) ? TMath::Power(10, mn1) : TMath::Power(2, mn1)) : TMath::Power(TMath::E(), mn1)) : mn1;
                  if (hasMn1) output << sep << mn1;
               }//if

               if (hasMn2 || hasFC) {
                  leaf = tree[k]->FindLeaf("mn2");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  mn2 = leaf->GetValue();
                  mn2 = (base > 0) ? ((base > 1) ? ((base > 2) ? TMath::Power(10, mn2) : TMath::Power(2, mn2)) : TMath::Power(TMath::E(), mn2)) : mn2;
                  if (hasMn2) output << sep << mn2;
               }//if

               if (hasSE) {
                  leaf = tree[k]->FindLeaf("se");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasDF) {
                  leaf = tree[k]->FindLeaf("df");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasPVal) {
                  leaf = tree[k]->FindLeaf("pval");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasNPer) {
                  leaf = tree[k]->FindLeaf("nper");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasPCha) {
                  leaf = tree[k]->FindLeaf("pcha");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasPAdj) {
                  leaf = tree[k]->FindLeaf("padj");
                  brch = leaf->GetBranch();
                  brch->GetEntry(i);
                  output << sep << leaf->GetValue();
               }//if

               if (hasFC) {
                  fc = (mn1 != 0) ? mn2/mn1 : 0;
                  output << sep << fc;
               }//if
            }//for_k
         }//if
      }//for_j

      // export flag
      if (hasFlag) {
         output << sep << flag[i];
      }//if
      output << endl;
   }//for_i

//Cleanup
cleanup:
   if (flag)   {delete [] flag;   flag = 0;}
   if (htable) {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   if (anntree) {
      SafeDelete(annot);
      anntree->ResetBranchAddress(anntree->GetBranch("AnnBranch"));
      SafeDelete(anntree);
   }//if

   if (unittree) {
      SafeDelete(unit);
      unittree->ResetBranchAddress(unittree->GetBranch("IdxBranch"));
      SafeDelete(unittree);
   }//if

   delete [] tree;

   return err;
}//ExportUnivarTrees

//______________________________________________________________________________
void XUnivarSet::AddHeader(const char *treename, Int_t treeid)
{
   // Add tree header with treename to list of trees
   if(kCS) cout << "------XUnivarSet::AddHeader------" << endl;

   switch (treeid) {
      case HEADER_UNITEST: {
         XUniTestHeader *header = new XUniTestHeader(treename, treeid);
         header->SetInfile(fInfile);
         header->SetType(fAnalyser->GetUniTest()->GetName());
         header->SetConfidenceLevel(fAnalyser->GetUniTest()->GetConfidenceLevel());
         header->SetMu(fAnalyser->GetUniTest()->GetMu());
         header->SetNumPerm(fAnalyser->GetUniTest()->GetNumPerm());
         header->SetAlternative(fAnalyser->GetUniTest()->GetAlternative());
         header->SetAdjustment(fAnalyser->GetUniTest()->GetAdjustment());
         header->SetIsPaired(fAnalyser->GetUniTest()->GetIsPaired());
         header->SetTwoSample(fAnalyser->GetUniTest()->GetTwoSample());
         if (strcmp(fAnalyser->GetUniTest()->GetName(), "ttest") == 0) {
            header->SetEqualVariance(((TStudentTest*)(fAnalyser->GetUniTest()))->GetEqualVariance());
         } else {
            header->SetEqualVariance(kFALSE);
         }//if

         fHeaders->Add(header);
         break;
      }//case

      case HEADER_UNIFILTER: {
         XUniFilterHeader *header = new XUniFilterHeader(treename, treeid);
         header->SetInfile(fInfile);
         header->SetType(fFilter->GetName());
         header->SetFCValue(fFilter->GetFCValue());
         header->SetFCDirection(fFilter->GetFCDirection());
         header->SetStatistic(fFilter->GetStatistic());
         header->SetPValue(fFilter->GetPValue());
         header->SetPChance(fFilter->GetPChance());
         header->SetPAdjust(fFilter->GetPAdjust());
         header->SetCallPValue(fFilter->GetCallPValue());
         header->SetCallCondition1(fFilter->GetCallCondition1());
         header->SetCallSamples1(fFilter->GetCallSamples1());
         header->SetCallCondition2(fFilter->GetCallCondition2());
         header->SetCallSamples2(fFilter->GetCallSamples2());
         header->SetHasStatistic(fFilter->HasStatistic());
         header->SetHasPValue(fFilter->HasPValue());
         header->SetHasPChance(fFilter->HasPChance());
         header->SetHasPAdjust(fFilter->HasPAdjust());
         header->SetFoldChangeFilter(fFilter->HasFoldChangeFilter());
         header->SetUniTestFilter(fFilter->HasUniTestFilter());
         header->SetCallFilter(fFilter->HasCallFilter());

         fHeaders->Add(header);
         break;
      }//case

      default:
         break;
   }//switch
}//AddTreeHeader

//______________________________________________________________________________
Int_t XUnivarSet::UniTest(Int_t n, TTree **tree, const char *leafname,
                  const char *outtree, const char *varlist, Option_t *option)
{
   // Univariate test
   if(kCS) cout << "------XUnivarSet::UniTest------" << endl;

// Create analysis tree
   Int_t   err   = errNoErr;
   TString tname = TString(outtree) + "." + fAnalyser->GetTitle();

   // get scheme name from leaftree
   TTree *atree = new TTree(tname, tree[0]->GetTitle());
   if (atree == 0) return errCreateTree;

// Get group indices
   Int_t *gid = 0;
   if (!(gid = new Int_t[n])) return errInitMemory;
   for (Int_t k=0; k<n; k++) gid[k] = 0;

   err = this->InitGroups(n, gid, tree, kExtenNorm);
   if (err) goto cleanup;

// Do univariate test
   err = fAnalyser->Analyse(n, gid, tree, leafname, atree, varlist, -1, 0);
   if (err) goto cleanup;

// Add tree info to tree
   AddTreeInfo(atree, atree->GetName(), option);

// Add expression tree header to list
   this->AddTreeHeader(tname, HEADER_UNITEST);

// Write expression tree to file 
   err = WriteTree(atree, 0); 
//ev kOverwrite to overwrite old atree with same name!! or keep old atree;cycle?
//better??   err = WriteTree(atree, TObject::kOverwrite);

// Add unit IDs from fTree as UnitBranch to atree
   CopyUnitBranch(tree[0], atree, TObject::kOverwrite);

// Set atree to current analysis tree fTree, used in this->Filter()
   fTree = atree;

// Clean up
cleanup:
   if (gid) {delete [] gid; gid = 0;}

   return err;
}//UniTest

//______________________________________________________________________________
Int_t XUnivarSet::Filter(Int_t n, TTree **tree, const char *leafname, Int_t nc,
                  TTree **calltree, const char *outtree, const char *varlist,
                  Option_t *option, Int_t base)
{
   // Univariate filter
   if(kCS) cout << "------XUnivarSet::Filter------" << endl;

   Int_t err = errNoErr;

// Apply unifilter to call trees
   if (nc > 0) {
      // get group indices
      Int_t *gid = 0;
      if (!(gid = new Int_t[nc])) return errInitMemory;
      for (Int_t i=0; i<nc; i++) gid[i] = 0;

      // init group indices
      err = this->InitGroups(nc, gid, calltree, kExtenCall);
      if (err != errNoErr) {delete [] gid; return err;}

      // apply unifilter
      err = fFilter->CallFlag(nc, gid, calltree, varlist, 0);
      if (err != errNoErr) {delete [] gid; return err;}

      if (gid) {delete [] gid; gid = 0;}
   }//if

// Get unitest tree: either from this->UniTest() or extern via AddTree()
   if (fTree == 0) {
   // get unitest tree from from list of trees
      TTree* tree = 0;
      for (Int_t k=0; k<fTrees->GetSize(); k++) {
         tree = (TTree*)fTrees->At(k);
         TString exten = Path2Name(tree->GetName(),".","");
         if (HasExtension(exten.Data(), kExtenUTst)) {
            fTree = tree;
            break;
         }//if
      }//for_k

      if (fTree == 0) return errGetTree;
   }//if

// Create unifilter tree
   TString tname = TString(outtree) + "." + fFilter->GetTitle();
   // get scheme name from fTree
   TTree *fltrtree = new TTree(tname, fTree->GetTitle());
   if (fltrtree == 0) return errCreateTree;

// Apply unifilter to unitest tree
   if (err == errNoErr) err = fFilter->Calculate(fTree, leafname, fltrtree, varlist, base);
   if (err != errNoErr) goto cleanup;

// Add tree info to tree
   AddTreeInfo(fltrtree, fltrtree->GetName(), option);

// Add expression tree header to list
   this->AddTreeHeader(tname, HEADER_UNIFILTER);

// Write expression tree to file 
   err = WriteTree(fltrtree, 0);
//ev kOverwrite to overwrite old fltrtree with same name!! or keep old tree;cycle?
//??   err = WriteTree(fltrtree, TObject::kOverwrite);

// Clean up
cleanup:
   // delete created tree from RAM
   if (fltrtree) {fltrtree->Delete(""); fltrtree = 0;}

   return err;
}//Filter


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultivarSet                                                         //
//                                                                      //
// Class for multivariate analysis of microarray data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMultivarSet::XMultivarSet()
             :XAnalySet()
{
   // Default MultivarSet constructor
   if(kCS) cout << "---XMultivarSet::XMultivarSet(default)------" << endl;

   fFilter   = 0;
   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XMultivarSet::XMultivarSet(const char *name, const char *type)
             :XAnalySet(name, type)
{
   // Normal MultivarSet constructor
   if(kCS) cout << "---XMultivarSet::XMultivarSet------" << endl;

   fFilter   = 0;
   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XMultivarSet::~XMultivarSet()
{
   // MultivarSet destructor
   if(kCS) cout << "---XMultivarSet::~XMultivarSet------" << endl;

   fFilter   = 0;  //deleted in ~XAnalySetting
   fAnalyser = 0;  //deleted in ~XAnalySetting
}//Destructor

//______________________________________________________________________________
Int_t XMultivarSet::Analyse(const char *infile, const char *outfile, 
                    const char *varlist, Int_t nrows, const char *sepi,
                    const char *sepo, char delim)
{
   // Analyse data
   if(kCS) cout << "------XMultivarSet::Analyse------" << endl;

   Int_t err = errNoErr;

// Calculate multivariate test
   if (fAnalyser) {
      err = fAnalyser->Analyse(infile, outfile, varlist, nrows, -1, 0, sepi, sepo, delim);
      if (err != errNoErr) return err;
   }//if

// Calculate filter to apply for multivariate test
   if (fFilter) {
      err = fFilter->Calculate(infile, outfile, varlist, nrows, sepi, sepo, delim);
   }//if

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XMultivarSet::Analyse(const char *leafname, const char *outtree,
                    const char *varlist)
{
   // Analyse data
   if(kCS) cout << "------XMultivarSet::Analyse------" << endl;

   Int_t err = errNoErr;

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XMultivarSet::ExportTreeType(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output
   if(kCS) cout << "------XMultivarSet::ExportTreeType------" << endl;

   if (HasExtension(exten, kExtenMTst)) {
      return this->ExportMultivarTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XMultivarSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output as XML-file
   if(kCS) cout << "------XMultivarSet::ExportTreeXML------" << endl;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XMultivarSet::Initialize(TFile *file, XSetting *setting, const char *infile,
                    const char *treename)
{
   // Initialize algorithms
   if(kCS) cout << "------XMultivarSet::Initialize------" << endl;

   if ((file == 0) || (setting == 0)) return errAbort;

   fFile     = file;
   fInfile   = infile;
   fTreeName = treename;

   fSetting = (XAnalySetting*)setting;
   if (!fSetting) return errAbort;

// Get files from settings
   fDataFile   = ((XAnalySetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XAnalySetting*)fSetting)->GetSchemeName();
   fSchemeType = ((XAnalySetting*)fSetting)->GetSchemeType();

// Get algorithms from settings
   fAnalyser = (XMultiTester*)((XAnalySetting*)fSetting)->GetAnalyser();

   // optional filter (can be: filter = 0) 
   fFilter = (XMultiFilter*)((XAnalySetting*)fSetting)->GetFilter();

// not needed for Export!!!
//   if (!(fAnalyser || fFilter)) return errAbort;

//?   fLogBase = ((XAnalySetting*)fSetting)->GetLogBase();
//?   fNegLog  = ((XAnalySetting*)fSetting)->GetNegLog();

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XMultivarSet::ExportMultivarTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XMultivarSet::ExportMultivarTrees------" << endl;

   Int_t err = errNoErr;

   return err;
}//ExportUnivarTrees


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XClusterSet                                                          //
//                                                                      //
// Class for cluster analysis of microarray data                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XClusterSet::XClusterSet()
            :XAnalySet()
{
   // Default ClusterSet constructor
   if(kCS) cout << "---XClusterSet::XClusterSet(default)------" << endl;

   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XClusterSet::XClusterSet(const char *name, const char *type)
            :XAnalySet(name, type)
{
   // Normal ClusterSet constructor
   if(kCS) cout << "---XClusterSet::XClusterSet------" << endl;

   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XClusterSet::~XClusterSet()
{
   // ClusterSet destructor
   if(kCS) cout << "---XClusterSet::~XClusterSet------" << endl;

   fAnalyser = 0;  //deleted in ~XAnalySetting
}//Destructor

//______________________________________________________________________________
Int_t XClusterSet::Analyse(const char *infile, const char *outfile, 
                   const char *varlist, Int_t nrows, const char *sepi,
                   const char *sepo, char delim)
{
   // Analyse data
   if(kCS) cout << "------XClusterSet::Analyse------" << endl;

   Int_t err = errNoErr;

// Calculate cluster analysis
   if (fAnalyser) {
      err = fAnalyser->Analyse(infile, outfile, varlist, nrows, -1, 0, sepi, sepo, delim);
      if (err != errNoErr) return err;
   }//if

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XClusterSet::Analyse(const char *leafname, const char *outtree,
                   const char *varlist)
{
   // Analyse data
   if(kCS) cout << "------XClusterSet::Analyse------" << endl;

   Int_t err = errNoErr;

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XClusterSet::ExportTreeType(const char *exten, Int_t n, TString *names, 
                   const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output
   if(kCS) cout << "------XClusterSet::ExportTreeType------" << endl;

   if (HasExtension(exten, kExtenClus)) {
      return this->ExportClusterTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XClusterSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                   const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output as XML-file
   if(kCS) cout << "------XClusterSet::ExportTreeXML------" << endl;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XClusterSet::Initialize(TFile *file, XSetting *setting, const char *infile,
                   const char *treename)
{
   // Initialize algorithms
   if(kCS) cout << "------XClusterSet::Initialize------" << endl;

   if ((file == 0) || (setting == 0)) return errAbort;

   fFile     = file;
   fInfile   = infile;
   fTreeName = treename;

   fSetting = (XAnalySetting*)setting;
   if (!fSetting) return errAbort;

// Get files from settings
   fDataFile   = ((XAnalySetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XAnalySetting*)fSetting)->GetSchemeName();
   fSchemeType = ((XAnalySetting*)fSetting)->GetSchemeType();

// Get algorithms from settings
   fAnalyser = (XClusterizer*)((XAnalySetting*)fSetting)->GetAnalyser();

//?   fLogBase = ((XAnalySetting*)fSetting)->GetLogBase();
//?   fNegLog  = ((XAnalySetting*)fSetting)->GetNegLog();

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XClusterSet::ExportClusterTrees(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XClusterSet::ExportClusterTrees------" << endl;

   Int_t err = errNoErr;

   return err;
}//ExportClusterTrees


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRegressionSet                                                       //
//                                                                      //
// Class for regression analysis of microarray data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRegressionSet::XRegressionSet()
               :XAnalySet()
{
   // Default RegressionSet constructor
   if(kCS) cout << "---XRegressionSet::XRegressionSet(default)------" << endl;

   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XRegressionSet::XRegressionSet(const char *name, const char *type)
               :XAnalySet(name, type)
{
   // Normal RegressionSet constructor
   if(kCS) cout << "---XRegressionSet::XRegressionSet------" << endl;

   fAnalyser = 0;
}//Constructor

//______________________________________________________________________________
XRegressionSet::~XRegressionSet()
{
   // RegressionSet destructor
   if(kCS) cout << "---XRegressionSet::~XRegressionSet------" << endl;

   fAnalyser = 0;  //deleted in ~XAnalySetting
}//Destructor

//______________________________________________________________________________
Int_t XRegressionSet::Analyse(const char *infile, const char *outfile, 
                      const char *varlist, Int_t nrows, const char *sepi,
                      const char *sepo, char delim)
{
   // Analyse data
   if(kCS) cout << "------XRegressionSet::Analyse------" << endl;

   Int_t err = errNoErr;

// Calculate cluster analysis
   if (fAnalyser) {
      err = fAnalyser->Analyse(infile, outfile, varlist, nrows, -1, 0, sepi, sepo, delim);
      if (err != errNoErr) return err;
   }//if

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XRegressionSet::Analyse(const char *leafname, const char *outtree,
                      const char *varlist)
{
   // Analyse data
   if(kCS) cout << "------XRegressionSet::Analyse------" << endl;

   Int_t err = errNoErr;

   return err;
}//Analyse

//______________________________________________________________________________
Int_t XRegressionSet::ExportTreeType(const char *exten, Int_t n, TString *names, 
                      const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output
   if(kCS) cout << "------XRegressionSet::ExportTreeType------" << endl;

   if (HasExtension(exten, kExtenRgrs)) {
      return this->ExportRegressionTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XRegressionSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                      const char *varlist, ofstream &output, const char *sep)
{
   // Export results stored in tree treename to file output as XML-file
   if(kCS) cout << "------XRegressionSet::ExportTreeXML------" << endl;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XRegressionSet::Initialize(TFile *file, XSetting *setting, const char *infile,
                      const char *treename)
{
   // Initialize algorithms
   if(kCS) cout << "------XRegressionSet::Initialize------" << endl;

   if ((file == 0) || (setting == 0)) return errAbort;

   fFile     = file;
   fInfile   = infile;
   fTreeName = treename;

   fSetting = (XAnalySetting*)setting;
   if (!fSetting) return errAbort;

// Get files from settings
   fDataFile   = ((XAnalySetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XAnalySetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XAnalySetting*)fSetting)->GetSchemeName();
   fSchemeType = ((XAnalySetting*)fSetting)->GetSchemeType();

// Get algorithms from settings
   fAnalyser = (XRegressor*)((XAnalySetting*)fSetting)->GetAnalyser();

//?   fLogBase = ((XAnalySetting*)fSetting)->GetLogBase();
//?   fNegLog  = ((XAnalySetting*)fSetting)->GetNegLog();

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XRegressionSet::ExportRegressionTrees(Int_t n, TString *names, const char *varlist,
                      ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XRegressionSet::ExportRegressionTrees------" << endl;

   Int_t err = errNoErr;

   return err;
}//ExportRegressionTrees





