// File created: 08/05/2002                          last modified: 02/06/2011

// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2011 Dr. Christian Stratowa
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
* May 2002 - Initial versions finished
* Jun 2002 - Remove InitXXX() (replaced with InitAlgorithm() in XManager)
* Dec 2002 - Completely reorganize source file structure
* Sep 2004 - Move XSelector classes and XNormalizer to own source files
* May 2005 - Replace XContent with XFolder.
* Jul 2005 - Replace TTree::AddFriend() with TList fTrees.
*          - Store trees for every XTreeSet in separate TDirectory.
* Jan 2006 - To decrease size of saved trees replace  Int_t with Short_t for 
*            mask/flag arrays
* Oct 2006 - Add support for exon arrays, class XNormedExonSet
* Aug 2007 - Add support for whole genome arrays, class XNormedGenomeSet
*
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <new>  //needed for new (nothrow)

#include "TBranch.h"
#include "THashTable.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TSystem.h"

//for Benchmark test only:
#include <TBenchmark.h>

#include "XPSNormation.h"
#include "XPSSelector.h"
#include "XPSNormalizer.h"
#include "TStat.h"

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XNormationManager);
ClassImp(XNormationSetting);
ClassImp(XNormedSet);
ClassImp(XNormedGCSet);
ClassImp(XNormedGenomeSet);
ClassImp(XNormedExonSet);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormationManager                                                    //
//                                                                      //
// Class for normalization of microarray data                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormationManager::XNormationManager()
                  :XProcessManager()
{
   // Default NormationManager constructor
   if(kCS) cout << "---XNormationManager::XNormationManager(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XNormationManager::XNormationManager(const char *name, const char *arraytype,
                   Int_t verbose)
                  :XProcessManager(name, arraytype, verbose)
{
   // Normal NormationManager constructor
   // If initfile is given, import default values
   if(kCS) cout << "---XNormationManager::XNormationManager------" << endl;

//Test Benchmark
//   if (!gBenchmark) gBenchmark = new TBenchmark();
}//Constructor

//______________________________________________________________________________
XNormationManager::~XNormationManager()
{
   // NormationManager destructor
   if(kCS) cout << "---XNormationManager::~XNormationManager------" << endl;

//Test Benchmark
//   SafeDelete(gBenchmark); 
}//Destructor

//______________________________________________________________________________
Int_t XNormationManager::Normalize(const char *setname, const char *method)
{
   // Normalize all expression trees added to tree set "setname" with AddTree().
   // The following methods can be used to normalize all initialized algorithms:
   // method = "normalize": normalize all initialized algorithms (default)
   if(kCS) cout << "------XNormationManager::Normalize------" << endl;

   if (fAbort) return errAbort;

// Create directory with name of treeset
   TDirectory *dir = fFile->GetDirectory(setname);
   if (!dir) fFile->mkdir(setname, fDataType);
   fFile->cd();

// Get treeset "setname" (set created in AddTree())
   XNormedSet *set = 0;
//??   set = (XNormedSet*)fContent->FindObject(setname, "XNormedSet");
   set = (XNormedSet*)fContent->FindObject(setname);
   if (!set) {
      return HandleError(errGetTreeSet, setname);
   }//if

// Set must be inherited from XNormedSet
   if (!set->InheritsFrom("XNormedSet")) {
      return HandleError(errClassTreeSet, setname, set->ClassName());
   }//if

// to lower
   TString smethod = method; smethod.ToLower();

// Initialize algorithms for methods "???"
   Int_t err = errNoErr;
   if (strcmp(smethod.Data(), "???") == 0) {
//?? place to initialize certain named normalization algorithm
   }//if

// Normalize trees using method smethod
//??   if ((set->GetNumTrees() > 1) && (set->GetNumSelections() > 1)) {
   if (set->GetNumSelections() > 1) {
      if (!err) err = set->Initialize(fFile, fSetting);
      if (!err) err = set->Normalize(smethod.Data());

      HandleError(err, "in XNormationManager::Normalize");
   } else {
      cerr << "Error: At least two trees need to be selected." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

//ev??   SafeDelete(set);  //no???
   return err;
}//Normalize

//______________________________________________________________________________
Int_t XNormationManager::ImportDefaults(const char *infile)
{
   // Import default values from infile
   if(kCS) cout << "------XNormationManager::ImportDefaults------" << endl;

   return errNoErr;
}//ImportDefaults

//______________________________________________________________________________
Int_t XNormationManager::InitDefaults()
{
   // Initialize default values
   if(kCS) cout << "------XNormationManager::InitDefaults------" << endl;

   Int_t err = errNoErr;

// Initialize algorithms
   if (!err) err = InitAlgorithm("selector","rank","separate",0, 4,0,0.3,400,0);
   if (!err) err = InitAlgorithm("normalizer","supsmu","log10",0, 2,0.0,0.0);
   if (!err) err = InitAlgorithm("normalizer","approx","linear:mean",0, 2,0.0,0.0);

   return err;
}//InitDefaults

//______________________________________________________________________________
XSetting *XNormationManager::NewSetting(const char *type, const char *infile)
{
   // Create new setting for type
   if(kCS) cout << "------XNormationManager::NewSetting------" << endl;

   XNormationSetting *setting = new XNormationSetting(type, infile);
   return setting;
}//NewSetting

//______________________________________________________________________________
XTreeSet *XNormationManager::NewTreeSet(const char *type)
{
   // Create new tree set of type
   if(kCS) cout << "------XNormationManager::NewTreeSet------" << endl;

   XTreeSet *set = 0;
   if (strcmp(type,"GeneChip") == 0) {
      set = new XNormedGCSet("NormedGCSet", type);
   } else if (strcmp(type,"GenomeChip") == 0) {
      set = new XNormedGenomeSet("NormedGenomeSet", type);
   } else if (strcmp(type,"ExonChip") == 0) {
      set = new XNormedExonSet("NormedExonSet", type);
   } else if (strcmp(type,"GenePix") == 0) {
//      set = new XNormedGPSet("NormedGPSet", type);
cout << "Need to be done: GenePix analysis" << endl;
   } else {
//      cerr << "Error: Array type <" << type << "> is not known" << endl;
      HandleError(errUnknownType, "Array", type);
   }//if

   return set;
}//NewTreeSet


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormationSetting                                                    //
//                                                                      //
// Class for initialization of normalization parameter settings         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormationSetting::XNormationSetting()
                  :XProcesSetting()
{
   // Default NormationSetting constructor
   if(kCS) cout << "---XNormationSetting::XNormationSetting(default)------" << endl;

   fSelector   = 0;
   fNormalizer = 0;
}//Constructor

//______________________________________________________________________________
XNormationSetting::XNormationSetting(const char *arraytype, const char *infile)
                  :XProcesSetting(arraytype, infile)
{
   // Normal NormationSetting constructor
   if(kCS) cout << "---XNormationSetting::XNormationSetting------" << endl;

   fSelector   = 0;
   fNormalizer = 0;
}//Constructor

//______________________________________________________________________________
XNormationSetting::~XNormationSetting()
{
   // NormationSetting destructor
   if(kCS) cout << "---XNormationSetting::~XNormationSetting------" << endl;

   SafeDelete(fSelector);
   SafeDelete(fNormalizer);
}//Destructor

//______________________________________________________________________________
Int_t XNormationSetting::InitAlgorithm(const char *name, const char *type,
                         Option_t *options, const char * /*filename*/,
                         Int_t npars, Double_t *pars)
{
   // Initialize algorithm "name" with "type" and "options"
   if(kCS) cout << "------XNormationSetting::InitAlgorithm------" << endl;

   if (strcmp(name, "selector") == 0) {
      return this->InitSelector(type, options, npars, pars);
   } else if ((strcmp(name, "normalizer") == 0) &&
              (strcmp(type, "approx")     != 0)) {
      return this->InitNormalizer(type, options, npars, pars);
   } else if ((strcmp(name, "normalizer") == 0) &&
              (strcmp(type, "approx")     == 0)) {
      return this->InitApprox(options, npars, pars);
   } else {
      cerr << "Error: Algorithm <" << name << "> is not known/applicable." << endl;
   }//if

   return errInitSetting;
}//InitAlgorithm

//______________________________________________________________________________
Int_t XNormationSetting::InitApprox(Option_t *options, Int_t npar, Double_t *pars)
{
   // Initialize normalizer algorithm with type = "approx" and parameters
   // and with options of the form "method:ties"
   //    - method: "linear" or "constant" (interpolation method used)
   //    - ties:   "mean", "min", "max", "ordered" (handling of tied values)
   //    parameters are: numpars, rule, f
   //    - numpars: number of other parameters as integer, i.e. numpars = 2:
   //    - rule:    how interpolation should be done at ends
   //    - f:       for method "constant" a value between 0 and 1
   // Note: approx must only be initialized for normalizers "lowess" and "supsmu".
   if(kCS) cout << "------XNormationSetting::InitApprox------" << endl;

   if (fNormalizer == 0) {
      cerr << "Error: Need to initialize Normalizer first" << endl;
      return errGeneral;
   }//if

   TString optcpy = options;
   TString method = "linear";
   TString ties   = "mean";

	char *opt = (char*)optcpy.Data();
   if (NumSeparators(options, ":") == 1) {
      method = strtok(opt, ":");
      ties   = strtok(NULL,":");
   } else {
      cout << "Warning: InitAlgorithm() must have two options for approx." << endl;
      cout << "         Using default options." << endl;
   }//if

   return fNormalizer->InitApprox(method.Data(), ties.Data(), npar, pars);
}//InitApprox

//______________________________________________________________________________
Int_t XNormationSetting::InitNormalizer(const char *type, Option_t *options,
                         Int_t npar, Double_t *pars)
{
   // Initialize normalizer algorithm and parameters
   //
   // type = "mean": scale trimmed mean of data with options "option:logbase"
   //      - "all":     use all data to calculate trimmed mean
   //      - "sel":     use data selected with selector to calculate trimmed mean
   //      - "logbase": "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars, trim, targetinten
   //    - numpars:     number of other parameters as integer, i.e. numpars = 2:
   //    - trim:        trim value for mean, in range [0, 0.5]
   //    - targetinten: if > 0 then scale trimmed mean to value of target intensity
   //
   // type = "median": scale median of data with options "option:logbase"
   //      - "all":     use all data to calculate median
   //      - "sel":     use data selected with selector to calculate median
   //      - "logbase": "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  targetinten
   //    - numpars:     number of other parameters as integer, i.e. numpars = 1:
   //    - targetinten: if > 0 then scale median to value of target intensity
   //
   // type = "ksmooth": apply kernel smoother with options:
   //      - "logbase": "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  bandwidth
   //    - numpars:   number of other parameters as integer, i.e. numpars = 1:
   //    - bandwidth: bandwidth of kernel smoother
   //
   // type = "lowess": apply Lowess smoother with options:
   //      - "logbase": "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  span, iter
   //    - numpars: number of other parameters as integer, i.e. numpars = 2:
   //    - span:    smoother span
   //    - iter:    number of robustifying iterations
   //
   // type = "supsmu": apply super smoother with options:
   //      - "logbase": "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //    parameters are: numpars,  bass, span
   //    - numpars: number of other parameters as integer, i.e. numpars = 2:
   //    - bass:    controls smoothness of fitted curve
   //    - span:    fraction of observations in span
   //
   //    type = "quantile": apply quantile normalization with options:
   //      - "together": quantile will be applied to PM and MM together
   //      - "separate": quantile will be applied to PM and MM separately
   //      - "logbase":  "0" - linear, "log", "log2", "log10" - log with base e, 2, 10
   //      Note: if background should be subtracted first, then options must be
   //            "option:bgrdoption:logbase" where bgrdoption is:
   //      - "subtractbg": subtract bgrd from intensity - result can be negative
   //      - "correctbg":  correct bgrd with noise fraction to avoid negative results
   //    parameters are: numpars, trim, (nfrac)
   //    - numpars: number of other parameters as integer, i.e. numpars = 1-2:
   //    - trim:    trim value for mean, in range [0, 0.5]
   //    - delta:   for robust ties, default 1.0 or 0.4
   if(kCS) cout << "------XNormationSetting::InitNormalizer------" << endl;

// Delete default normalizer setting
   SafeDelete(fNormalizer);

   TString exten = Type2Extension(type, kTypeNorm, kExtenNorm);
   TString stype = Extension2Type(type, kTypeNorm, kExtenNorm);

   if (strcmp(exten.Data(), kExtenNorm[0]) == 0) {
      fNormalizer = new XMeanNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenNorm[1]) == 0) {
      fNormalizer = new XMedianNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenNorm[2]) == 0) {
      fNormalizer = new XKernelNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenNorm[3]) == 0) {
      fNormalizer = new XLowessNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenNorm[4]) == 0) {
      fNormalizer = new XSuperNormalizer(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenNorm[5]) == 0) {
      fNormalizer = new XQuantileNormalizer(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Normalizer <" << type << "> is not known." << endl;
      return errInitSetting;
   }//if
   if (fNormalizer == 0) return errInitMemory;

   fNormalizer->SetOptions(options);

   return fNormalizer->InitParameters(npar, pars);
}//InitNormalizer

//______________________________________________________________________________
Int_t XNormationSetting::InitSelector(const char *type, Option_t *options,
                         Int_t npar, Double_t *pars)
{
   // Initialize selector algorithm of type with options:
   // type = "default": default selector with options 
   //      - "none": does not change the selection mask
   //      - "all":  select all entries, i.e. flag = 1 for all units
   //    parameters are: numpars
   //    - numpars:  number of other parameters as integer, i.e. numpars = 0:
   //
   // type = "rank": selector based on difference of ranks, with options
   //      - "separate": use separate mask for each normalization
   //      - "together": use numflags masks for normalization
   //    parameters are: numpars, cutoff, trim, percent, numflags
   //    - numpars:  number of other parameters as integer, i.e. numpars = 4:
   //    - cutoff:   cutoff value for selection of units, if = 0 then calculated
   //    - trim:     trim value for the fitting line to calculate cutoff
   //    - percent:  minimum percentage of units to be selected
//replace with percent   //    - minunits: minimum number of units to be selected
   //    - numflags: number of flags, e.g. 9999 (=all) for intersect, 1 for union
   //
   // type = "user": user can either select "all" entries or supply a mask file
   //      for the selection of tree entries, with following options:
   //      - "all:separate" or "all:together", all units are used
   //      - "fullname.txt:separate" or "fullname.txt:together", where
   //        "fullname.txt" is the filename of the mask to import, and mask
   //         need only contain unitIDs and flag = 1 of selected units
   //    parameters are: numpars
   //    - numpars:  number of other parameters as integer, i.e. numpars = 0:
   if(kCS) cout << "------XNormationSetting::InitSelector------" << endl;

// Delete default selector setting
   SafeDelete(fSelector);

   TString exten = Type2Extension(type, kTypeSlct, kExtenSlct);
   TString stype = Extension2Type(type, kTypeSlct, kExtenSlct);

   if (strcmp(exten.Data(), kExtenSlct[0]) == 0) {
      fSelector = new XSelector(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenSlct[1]) == 0) {
      fSelector = new XRankSelector(stype.Data(), exten.Data());
   } else if (strcmp(exten.Data(), kExtenSlct[4]) == 0) {
      fSelector = new XUserSelector(stype.Data(), exten.Data());
   } else {
      cerr << "Error: Selector <" << type << "> is not known/applicable." << endl;
      return errInitSetting;
   }//if
   if (fSelector == 0) return errInitMemory;

   fSelector->SetOption(options);

   return fSelector->InitParameters(npar, pars);
}//InitSelector


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedSet                                                           //
//                                                                      //
// Base class for normalization of microarray data                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormedSet::XNormedSet()
           :XProcesSet()
{
   // Default NormedSet constructor
   if(kCS) cout << "---XNormedSet::XNormedSet(default)------" << endl;

   fSelector   = 0;
   fNormalizer = 0;
//   fSF         = 1.0;
}//Constructor

//______________________________________________________________________________
XNormedSet::XNormedSet(const char *name, const char *type)
           :XProcesSet(name, type)
{
   // Normal NormedSet constructor
   if(kCS) cout << "---XNormedSet::XNormedSet------" << endl;

   fSelector   = 0;
   fNormalizer = 0;
//   fSF         = 1.0;
}//Constructor

//______________________________________________________________________________
XNormedSet::~XNormedSet()
{
   // NormedSet destructor
   if(kCS) cout << "---XNormedSet::~XNormedSet------" << endl;

  // deleted in ~XNormationSetting
   fSelector   = 0;
   fNormalizer = 0;
}//Destructor

//______________________________________________________________________________
Int_t XNormedSet::Initialize(TFile *file, XSetting *setting,
                  const char *infile, const char *treename)
{
   // Initialize algorithms
   if(kCS) cout << "------XNormedSet::Initialize------" << endl;

   if ((file == 0) || (setting == 0)) return errAbort;

   fFile     = file;
   fSetting  = (XNormationSetting*)setting;
   fInfile   = infile;
   fTreeName = treename;

// Get files from settings
   fDataFile   = ((XNormationSetting*)fSetting)->GetDataFile();
   fSchemeFile = ((XNormationSetting*)fSetting)->GetSchemeFile();
   fSchemeName = ((XNormationSetting*)fSetting)->GetSchemeName();

// Get algorithms from settings
   fSelector   = ((XNormationSetting*)fSetting)->GetSelector();
   fNormalizer = ((XNormationSetting*)fSetting)->GetNormalizer();

   if (fSetting && fSelector && fNormalizer && fFile) {return errNoErr;}
   return errAbort;
}//Initialize

//______________________________________________________________________________
void XNormedSet::AddMaskTreeInfo(TTree *tree, const char *name, Option_t *option,
                 Int_t nunits, Int_t nflags)
{
   // Add mask tree info to list fUserInfo of tree
   if(kCS) cout << "------XNormedSet::AddMaskTreeInfo------" << endl;

// store name of tree set as title
   XSelectionTreeInfo *info = new XSelectionTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info
   info->AddUserInfo(nunits, nflags);

   tree->GetUserInfo()->Add(info);
}//AddMaskTreeInfo

//______________________________________________________________________________
Int_t XNormedSet::ExportTreeInfo(const char *exten, Int_t n, TString *names,  
                  const char *varlist, ofstream &output, const char *sep)
{
   // Export variables from varlist for trees of type exten
   if(kCS) cout << "------XNormedSet::ExportTreeInfo------" << endl;

// Set scheme file to be able to access scheme data for exporting
   if (fSetting) {
      fSchemeFile = ((XNormationSetting*)fSetting)->GetSchemeFile();
   }//if

   // remove "userinfo" from varlist
   TString infolist = RemoveSubString(varlist, "userinfo:", kFALSE);

   if (HasExtension(exten, kExtenNorm)) {
      return ExportExprTreeInfo(n, names, infolist, output, sep);
   } else if (HasExtension(exten, kExtenSlct)) {
      return this->ExportMaskTreeInfo(n, names, infolist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeInfo

//______________________________________________________________________________
Int_t XNormedSet::ExportTreeType(const char *exten, Int_t n, TString *names,  
                  const char *varlist, ofstream &output, const char *sep)
{
   // Export variables from varlist for trees of type exten
   if(kCS) cout << "------XNormedSet::ExportTreeType------" << endl;

// Set scheme file to be able to access scheme data for exporting
   if (fSetting) {
      fSchemeFile = ((XNormationSetting*)fSetting)->GetSchemeFile();
   }//if

   if (HasExtension(exten, kExtenNorm)) {
      return this->ExportExprTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kExtenSlct)) {
      return this->ExportMaskTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return errNoErr;
}//ExportTreeType

//______________________________________________________________________________
Int_t XNormedSet::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                  const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XNormedSet::ExportTreeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   exten = 0; n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XNormedSet::FillExprArray(TTree *tree, Int_t n, Int_t *idx, Double_t *arr)
{
   // Fill arrays idx and arr with unitID and data from expression tree with name
   if(kCS) cout << "------XNormedSet::FillExprArray------" << endl;

   TBranch *brch = tree->FindBranch("ExprBranch");
   if (!brch) return errGetTree;

   XExpression *expr = 0;
   brch->SetAddress(&expr);
   for (Int_t i=0; i<n; i++) {
      brch->GetEntry(i);
      idx[i] = expr->GetUnitID();
      arr[i] = expr->GetLevel();
   }//for_i

   SafeDelete(expr);
   tree->ResetBranchAddress(tree->GetBranch("ExprBranch"));

   return errNoErr;
}//FillExprArray

//______________________________________________________________________________
Int_t XNormedSet::FillExprTree(const char *name, Int_t n, Int_t *idx, Double_t *arr)
{
   // Fill expression tree with arrays idx and arr and write to file
   if(kCS) cout << "------XNormedSet::FillExprTree------" << endl;

   Int_t err = errNoErr;

   Int_t     nquant = 7;
   Double_t  q[]    = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};
   Double_t *quantL = 0; 
   if (!(quantL = new (nothrow) Double_t[nquant])) return errInitMemory;

   TTree *tree = new TTree(name, fSchemeName.Data());
   if (tree == 0) return errCreateTree;

////////////////////
// PROBLEM: ev XGCExpression with fStdev,fNPairs??
// => must be defined in XGCNormedSet::FillExprTree()
////////////////////
   Int_t split = 99;
   XExpression *expr = 0;
   expr = new XExpression();
   tree->Branch("ExprBranch", "XExpression", &expr, 64000, split);

   for (Int_t i=0; i<n; i++) {
      expr->SetUnitID(idx[i]);
      expr->SetLevel(arr[i]);
      tree->Fill();
   }//for_i

// Quantiles for expression trees
   err = ExpressionQuantiles(tree, expr, nquant, q, quantL);
   if (err != errNoErr) goto cleanup;

// Add tree info to tree
   AddExprTreeInfo(tree, name, fNormalizer->GetOption(),
                   n, quantL[0], quantL[nquant-1], nquant, q, quantL);

// Write expression tree to file 
   if ((err = WriteTree(tree, TObject::kOverwrite)) == errNoErr) {
   // Delete tree header in case of overwrite
      XTreeHeader *header = GetTreeHeader(name);
      if (header) {fHeaders->Remove(header); delete header;}  //????

   // Add tree header to list
      AddTreeHeader(tree->GetName(), fNormalizer->GetName(), 0,
                    fNormalizer->GetNumParameters(), fNormalizer->GetParameters());
   }//if

// Cleanup
cleanup:
   SafeDelete(expr);
   tree->ResetBranchAddress(tree->GetBranch("ExprBranch"));
   SafeDelete(tree);

   if (quantL) {delete [] quantL; quantL = 0;}

   return err;
}//FillExprTree

//______________________________________________________________________________
Int_t XNormedSet::FillMaskArray(const char *name, Int_t n, Int_t *arr)
{
   // Fill array arr with mask from tree with name
   // or calculate mask from intersection/union of mask trees
   if(kCS) cout << "------XNormedSet::FillMaskArray------" << endl;

   Int_t err = errNoErr;
   Int_t idx;

   if (strcmp(name, "together") == 0) {
   // Count number of mask trees with selected extension
      Int_t nummask = 0;
      TKey *key = 0;
      TIter next(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TString xten = Path2Name(key->GetName(),".",";");
         if (strcmp(xten.Data(), fSelector->GetTitle()) == 0) {
            nummask++;
         }//if
      }//while

   // Initialize mask trees
      TTree **tree = new TTree*[nummask];
      XMask **mask = new XMask*[nummask];
      for (Int_t k=0; k<nummask; k++) {
         tree[k] = 0;
         mask[k] = 0;
      }//for_k

   // Init arrays
      TString *arrName = 0;  //to store alias names of friend trees
      Int_t   *arrFlag = 0;  //to store flag values of friend branches
      if (!(arrName = new (nothrow) TString[nummask])) {err = errInitMemory; goto cleanup;}
      if (!(arrFlag = new (nothrow) Int_t[nummask]))   {err = errInitMemory; goto cleanup;} 

   // Get mask trees
      key = 0;
      idx = 0;
      next.Reset();
      while ((key = (TKey*)next())) {
         TString xten = Path2Name(key->GetName(),".",";");
         if (strcmp(xten.Data(), fSelector->GetTitle()) == 0) {
            tree[idx] = (TTree*)gDirectory->Get(key->GetName());
            if (!tree[idx]) {err = errGetTree; goto cleanup;}

            tree[idx]->SetBranchAddress("MaskBranch", &mask[idx]);
            idx++;
         }//if
      }//while

   // Fill mask array
      for (Int_t i=0; i<n; i++) {
         for (Int_t k=0; k<nummask; k++) {
            tree[k]->GetEntry(i);
            arrFlag[k] = mask[k]->GetFlag();
         }//for_k
         arr[i] = fSelector->SetFlag(nummask, arrFlag);
      }//for_i

   // Clean up
   cleanup:
      delete [] arrFlag;
      delete [] arrName;

      for (Int_t k=0; k<nummask; k++) {
         SafeDelete(mask[k]);
         tree[k]->ResetBranchAddress(tree[k]->GetBranch("MaskBranch"));
         SafeDelete(tree[k]);
      }//for_k

      delete [] mask;
      delete [] tree;
   } else {
   // Get tree with name.exten
      TString tname = TString(name) + "." + fSelector->GetTitle();
      TTree  *tree  = (TTree*)gDirectory->Get(tname); 
      if (tree == 0) return errGetTree;

      XMask *mask = 0;
      tree->SetBranchAddress("MaskBranch", &mask);

   // Loop over all entries of tree
      if (n != (Int_t)(tree->GetEntries())) {
         TString str = ""; str += n;
         return fManager->HandleError(errNumTreeEntries, str);
      }//if
      for (Int_t i=0; i<n; i++) {
         tree->GetEntry(i);
         arr[i] = mask->GetFlag();
      }//for_i

      SafeDelete(mask);
      tree->ResetBranchAddress(tree->GetBranch("MaskBranch"));
      SafeDelete(tree);
   }//if

   return err;
}//FillMaskArray

//______________________________________________________________________________
Int_t XNormedSet::FillMaskTree(const char *name, Int_t n, Int_t *idx, Int_t *arr)
{
   // Fill mask tree with arrays idx and arr and write to file
   if(kCS) cout << "------XNormedSet::FillMaskTree------" << endl;

   Int_t err = errNoErr;

   TTree *tree = new TTree(name, fSchemeName.Data());
   if (tree == 0) return errCreateTree;

   Int_t split   = 99;
   XUnitID *unit = 0;
   XMask   *mask = 0;
//?? better: XUnitMask (XUnitID not necessary)
   unit = new XUnitID();
   mask = new XMask();
   tree->Branch("UnitBranch", "XUnitID", &unit, 64000, split);
   tree->Branch("MaskBranch", "XMask", &mask, 64000, split);

// Fill mask tree
   Int_t nflags = 0;
   for (Int_t i=0; i<n; i++) {
      // get minimal/maximal expression levels
      if (arr[i] > 0) nflags++;

      unit->SetUnitID(idx[i]);
      mask->SetFlag((Short_t)arr[i]);
      tree->Fill();
   }//for_i

// Add tree info to tree
   AddMaskTreeInfo(tree, name, fSelector->GetOption(), n, nflags);

// Write mask tree to file 
   if ((err = WriteTree(tree, TObject::kOverwrite)) == errNoErr) {
   // Delete tree header in case of overwrite
      XTreeHeader *header = GetTreeHeader(name);
      if (header) {fHeaders->Remove(header); delete header;}  //????

   // Write tree header to list
      AddTreeHeader(tree->GetName(), fSelector->GetName(), 0,
                    fSelector->GetNumParameters(), fSelector->GetParameters());
   }//if

// Cleanup
   SafeDelete(unit);
   SafeDelete(mask);
   tree->ResetBranchAddress(tree->GetBranch("UnitBranch"));
   tree->ResetBranchAddress(tree->GetBranch("MaskBranch"));
   SafeDelete(tree);

   return err;
}//FillMaskTree

//______________________________________________________________________________
Int_t XNormedSet::MeanReference(Int_t numexpr, TTree **exprtree,
                  Int_t n, Int_t *idx, Double_t *arr)
{
   // Fill array arr with trimmed mean values from numexpr exprtrees
   if(kCS) cout << "------XNormedSet::MeanReference------" << endl;

   TBranch     **brch = new TBranch*[numexpr];
   XExpression **expr = new XExpression*[numexpr];
   for (Int_t k=0; k<numexpr; k++) {
      expr[k] = 0;
      brch[k] = exprtree[k]->GetBranch("ExprBranch");
      brch[k]->SetAddress(&expr[k]);
   }//for_k

   Double_t *arrRef = 0;  //to store entries of tree branches
   if (!(arrRef = new (nothrow) Double_t[numexpr])) return errInitMemory; 

// Fill reference array
   if (numexpr > 1) {
      for (Int_t i=0; i<n; i++) {
         brch[0]->GetEntry(i);
         arrRef[0] = expr[0]->GetLevel();

         for (Int_t k=1; k<numexpr; k++) {
            brch[k]->GetEntry(i);
            arrRef[k] = expr[k]->GetLevel();
         }//for_j

         idx[i] = expr[0]->GetUnitID();
         arr[i] = TStat::Mean(numexpr, arrRef, fRefTrim);
      }//for_i
   } else {
      for (Int_t i=0; i<n; i++) {
         brch[0]->GetEntry(i);

         idx[i] = expr[0]->GetUnitID();
         arr[i] = expr[0]->GetLevel();
      }//for_i
   }//if

// Cleanup
   for (Int_t k=0; k<numexpr; k++) {
      SafeDelete(expr[k]);
      exprtree[k]->ResetBranchAddress(exprtree[k]->GetBranch("ExprBranch"));
   }//for_k

   delete [] arrRef;
   delete [] expr;
   delete [] brch;

   return errNoErr;
}//MeanReference

//______________________________________________________________________________
Int_t XNormedSet::MedianReference(Int_t numexpr, TTree **exprtree,
                  Int_t n, Int_t *idx, Double_t *arr)
{
   // Fill array arr with median values from numexpr exprtrees
   if(kCS) cout << "------XNormedSet::MedianReference------" << endl;

   TBranch     **brch = new TBranch*[numexpr];
   XExpression **expr = new XExpression*[numexpr];
   for (Int_t k=0; k<numexpr; k++) {
      expr[k] = 0;
      brch[k] = exprtree[k]->GetBranch("ExprBranch");
      brch[k]->SetAddress(&expr[k]);
   }//for_k

   Double_t *arrRef = 0;  //to store entries of tree branches
   if (!(arrRef = new (nothrow) Double_t[numexpr])) return errInitMemory; 

// Fill reference array
   if (numexpr > 1) {
      for (Int_t i=0; i<n; i++) {
         brch[0]->GetEntry(i);
         arrRef[0] = expr[0]->GetLevel();

         for (Int_t k=1; k<numexpr; k++) {
            brch[k]->GetEntry(i);
            arrRef[k] = expr[k]->GetLevel();
         }//for_j

         idx[i] = expr[0]->GetUnitID();
         arr[i] = TStat::Median(numexpr, arrRef);
      }//for_i
   } else {
      for (Int_t i=0; i<n; i++) {
         brch[0]->GetEntry(i);

         idx[i] = expr[0]->GetUnitID();
         arr[i] = expr[0]->GetLevel();
      }//for_i
   }//if

// Cleanup
   for (Int_t k=0; k<numexpr; k++) {
      SafeDelete(expr[k]);
      exprtree[k]->ResetBranchAddress(exprtree[k]->GetBranch("ExprBranch"));
   }//for_k

   delete [] arrRef;
   delete [] expr;
   delete [] brch;

   return errNoErr;
}//MedianReference


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedGCSet                                                         //
//                                                                      //
// Class for normalization of GeneChip oligonucleotide arrays           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormedGCSet::XNormedGCSet()
             :XNormedSet()
{
   // Default NormedGeneChipSet constructor
   if(kCS) cout << "---XNormedGCSet::XNormedGCSet(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XNormedGCSet::XNormedGCSet(const char *name, const char *type)
             :XNormedSet(name, type)
{
   // Normal NormedGeneChipSet constructor
   if(kCS) cout << "---XNormedGCSet::XNormedGCSet------" << endl;

}//Constructor

//______________________________________________________________________________
XNormedGCSet::~XNormedGCSet()
{
   // NormedGeneChipSet destructor
   if(kCS) cout << "---XNormedGCSet::~XNormedGCSet------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XNormedGCSet::Normalize(const char *method)
{
   // Normalize GeneChip expression data
   if(kCS) cout << "------XNormedGCSet::Normalize------" << endl;

// Informing user
   if (XManager::fgVerbose) {
      cout << "Normalizing expression data using method <" 
           << fNormalizer->GetName() << ">..." << endl;
   }//if

// Change directory
   TDirectory *savedir = gDirectory;
   if (!fFile->cd(fName)) return errGetDir;

   TString treename, selname, refname;

   Int_t err     = errNoErr;
   Int_t numsels = fTrees->GetSize();
   Int_t numrefs = fReferences->GetSize();

// Initialize selected trees and reference trees
   TTree **seltree = new TTree*[numsels];
   TTree **reftree = new TTree*[numrefs];
   for (Int_t k=0; k<numsels; k++) seltree[k] = 0;
   for (Int_t j=0; j<numrefs; j++) reftree[j] = 0;

// Get selected trees
   for (Int_t k=0; k<numsels; k++) {
      seltree[k] = (TTree*)fTrees->At(k);
      if (seltree[k] == 0) return errGetTree;
      if (seltree[k]->GetBranch("ExprBranch") != 0) {
//????
         treename = Path2Name(seltree[k]->GetName(), dSEP, ".");
         seltree[k]->SetName(treename);
      } else {
         cerr << "Error: Tree <" << seltree[k]->GetName() << "> has no ExprBranch."
              << endl;
         return errGetTree;
      }//if
   }//for_k

// Get reference trees
   for (Int_t j=0; j<numrefs; j++) {
      reftree[j] = (TTree*)fReferences->At(j);
      if (reftree[j] == 0) return errGetTree;
      if (reftree[j]->GetBranch("ExprBranch") != 0) {
//????
         treename = Path2Name(reftree[j]->GetName(), dSEP, ".");
         reftree[j]->SetName(treename);
      } else {
         cerr << "Error: Tree <" << reftree[j]->GetName() << "> has no ExprBranch."
              << endl;
         return errGetTree;
      }//if
   }//for_j

////////////
/*TO DO
// Get chip parameters from scheme file (also alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = datatree[0]->GetTitle();
   } else if (!fSchemeName.Contains(datatree[0]->GetTitle())) {
      return fManager->HandleError(errSchemeDerived, fSchemeName, datatree[0]->GetTitle());
   }//if
*/
/////////////
   fSchemeName = seltree[0]->GetTitle();
   Int_t entries = (Int_t)(seltree[0]->GetEntries());
   TObjString *refstr = 0;

// Initialize local arrays
   Double_t *arrX   = 0;  //to store data of reference tree
   Double_t *arrY   = 0;  //to store data of selected tree
   Double_t *dummy  = 0;  //to prevent compilation error
   Int_t    *arrIdx = 0;  //to store unitIDs
   Int_t    *arrMsk = 0;  //mask, set to 1 if unit is selected
   if (!(arrX   = new (nothrow) Double_t[entries])) {err = errInitMemory; goto cleanup;}
   if (!(arrY   = new (nothrow) Double_t[entries])) {err = errInitMemory; goto cleanup;}
   if (!(arrIdx = new (nothrow) Int_t[entries]))    {err = errInitMemory; goto cleanup;}
   if (!(arrMsk = new (nothrow) Int_t[entries]))    {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<entries; i++) {
      arrX[i]   = 0;
      arrY[i]   = 0;
      arrIdx[i] = 0;
      arrMsk[i] = 1;  //?? if selector type = none
   }//for_i

// Compute reference array(s) arrX
   if (strcmp(fNormalizer->GetName(), "quantile") == 0) {
      for (Int_t k=0; k<numsels; k++) {
         if ((err = FillExprArray(seltree[k], entries, arrIdx, arrX))) goto cleanup;

         treename = seltree[k]->GetName();
         if ((err = fNormalizer->AddArray(entries, arrX, arrMsk, treename))) goto cleanup;
      }//for_k
   } else if (numrefs == 1) {
      refstr = new TObjString(reftree[0]->GetName());

      if ((err = FillExprArray(reftree[0], entries, arrIdx, arrX))) goto cleanup;

      refname  = reftree[0]->GetName();
      treename = refname + "." + fNormalizer->GetTitle();
      if ((err = FillExprTree(treename, entries, arrIdx, arrX))) goto cleanup;
   } else if (numrefs > 1) {
      if (strcmp(fRefOpt.Data(), "mean") == 0) {
         if ((err = MeanReference(numrefs, reftree, entries, arrIdx, arrX))) goto cleanup;
      } else if (strcmp(fRefOpt.Data(), "median") == 0) {
         if ((err = MedianReference(numrefs, reftree, entries, arrIdx, arrX))) goto cleanup;
      }//if

      treename = TString(kReference) + "." + TString(fNormalizer->GetTitle());
      if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrX))) goto cleanup;

      // add reference tree to selections list
      Select(kReference, 1);
      refstr = new TObjString(kReference);
//      refstr = (TObjString*)(fSelections->Last());
   } else {
      cerr << "Error: No reference tree is selected." << endl;
      err = errGeneral;
      goto cleanup;
   }//if

//?? Ev draw graph for ranks, H2all ???
//   TCanvas *c1 = new TCanvas("c1","Normalize data",200,10,700,700);
//   c1->Divide(2,2); //??(3,3) dependent on numtrees (+1) ??
//ev. interactively step thru numtrees and stop execution on demand

// Normalization
   if (strcmp(fNormalizer->GetName(), "quantile") == 0) {
      // informing user
      if (XManager::fgVerbose) {
         cout << "   quantile normalization..." << endl;
      }//if

      // quantile normalization
      if ((err = fNormalizer->Calculate(entries, arrX, arrY, arrMsk))) goto cleanup;

      for (Int_t k=0; k<numsels; k++) {
         TString name = Path2Name(seltree[k]->GetName(),"",".");

         arrY = fNormalizer->GetArray(entries, arrY, arrMsk, name);
         if (arrY == 0) {err = errAbort; goto cleanup;}

         treename = name + "." + fNormalizer->GetTitle();
         if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrY))) break;

         treename = name + "." + fSelector->GetTitle();
         if ((err = FillMaskTree(treename.Data(), entries, arrIdx, arrMsk))) break;
      }//for_k
//////////////////////
//PROBLEM: default selector option = "all" or "none"
//////////////////////
//   } else if (strcmp(fSelector->GetOption(), "separate") == 0) {
   } else if ((strcmp(fSelector->GetName(),   "default")  == 0) ||
             ((strcmp(fSelector->GetName(),   "rank")     == 0) &&
              (strcmp(fSelector->GetOption(), "separate") == 0))) {
   // Separate mask calculated for each pair [arrX,arrY] to be normalized
      for (Int_t k=0; k<numsels; k++) {
         selname  = seltree[k]->GetName();

         // informing user
         if (XManager::fgVerbose) {
//            cout << "   normalizing <" << selname.Data() << ">..." << endl;
            cout << "   normalizing tree " << k+1 << " of " << numsels << ": <"
                 << selname.Data() << ">..." << endl;
         }//if

         // for mean/median the tree used as reference must be normalized, too!
         // e.g. to scale tree to targetinten
         if (strcmp(selname.Data(),refstr->GetName()) == 0) {
            // skip if reference tree is kReference (numrefs > 1)
            if (strcmp(selname.Data(), kReference) == 0) {
               continue;
            }//if

            // skip reference tree if not mean/median 
            if (!((strcmp(fNormalizer->GetName(),"mean")   == 0) ||
                  (strcmp(fNormalizer->GetName(),"median") == 0))) {
               continue;
            }//if
         }//if

         if ((err = FillExprArray(seltree[k], entries, arrIdx, arrY))) break;

         // select non-variant units and fill mask tree
         if ((err = fSelector->Calculate(entries, arrX, arrY, arrMsk))) break;
         treename = selname + "." + fSelector->GetTitle();
         if ((err = FillMaskTree(treename.Data(), entries, arrIdx, arrMsk))) break;

         // normalize expression levels arrY and fill expression tree
         if ((err = fNormalizer->Calculate(entries, arrX, arrY, arrMsk))) break;
         treename = selname + "." + fNormalizer->GetTitle();
         if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrY))) break;
      }//for_k
   } else if ((strcmp(fSelector->GetName(),   "rank")     == 0) &&
              (strcmp(fSelector->GetOption(), "together") == 0)) {
   // One combined mask calculated for all pairs [arrX,arrY] to be normalized
//   Int_t length;        
//   do {
      for (Int_t k=0; k<numsels; k++) {
         selname  = seltree[k]->GetName();

         if (strcmp(selname.Data(),refstr->GetName()) == 0) continue;

         if ((err = FillExprArray(seltree[k], entries, arrIdx, arrY))) break;

         // select non-variant units and fill mask tree
         if ((err = fSelector->Calculate(entries, arrX, arrY, arrMsk))) break;

         treename = selname + "." + fSelector->GetTitle();
         if ((err = FillMaskTree(treename.Data(), entries, arrIdx, arrMsk))) break;
      }//for_k
      if (err != errNoErr) goto cleanup;

      // compute combined mask array arrMsk
      if ((err = FillMaskArray("together", entries, arrMsk))) goto cleanup;

      // get number of selected units
      Int_t length = 0;
      for (Int_t i=0; i<entries; i++) {
         length += arrMsk[i];
      }//for_i
//      Bool_t tooSmall = (length < fSelector->GetMinunits());
//      if (tooSmall) fSelector->IncreaseCutoff();
//      if (tooSmall) fSelector->IncreaseCutoff(0.2);
//   } while (tooSmall);
//   } while (length < fSelector->GetMinunits());

      // normalize data and store as trees
      for (Int_t k=0; k<numsels; k++) {
         selname  = seltree[k]->GetName();

         // informing user
         if (XManager::fgVerbose) {
//            cout << "   normalizing <" << selname.Data() << ">..." << endl;
            cout << "   normalizing tree " << k+1 << " of " << numsels << ": <"
                 << selname.Data() << ">..." << endl;
         }//if

         // if not mean/median fill reference tree only
         if (strcmp(selname.Data(),refstr->GetName()) == 0) {
            if (!((strcmp(fNormalizer->GetName(), "mean")   == 0) ||
                  (strcmp(fNormalizer->GetName(), "median") == 0))) {
               if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrX))) break;
               continue;
            }//if
         }//if

         if ((err = FillExprArray(seltree[k], entries, arrIdx, arrY))) break;

         // normalize expression levels arrY and fill expression tree
         if ((err = fNormalizer->Calculate(entries, arrX, arrY, arrMsk))) break;

         treename = selname + "." + fNormalizer->GetTitle();
         if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrY))) break;
      }//for_k
   } else if (strcmp(fSelector->GetName(), "user") == 0) {
   // User-defined mask imported for all pairs [arrX,arrY] to be normalized
      if ((err = fSelector->Calculate(entries, dummy, dummy, arrMsk))) goto cleanup;
      treename = Path2Name(fSelector->GetOption(), dSEP, ".") + "." + fSelector->GetTitle();
      if ((err = FillMaskTree(treename.Data(), entries, arrIdx, arrMsk))) goto cleanup;

      // normalize data and store as trees
      for (Int_t k=0; k<numsels; k++) {
         selname  = seltree[k]->GetName();

         // informing user
         if (XManager::fgVerbose) {
//            cout << "   normalizing <" << selname.Data() << ">..." << endl;
            cout << "   normalizing tree " << k+1 << " of " << numsels << ": <"
                 << selname.Data() << ">..." << endl;
         }//if

         // if not mean/median fill reference tree only
         if (strcmp(selname.Data(),refstr->GetName()) == 0) {
            if (!((strcmp(fNormalizer->GetName(), "mean")   == 0) ||
                  (strcmp(fNormalizer->GetName(), "median") == 0))) {
//?? FillExprTree second time->overwrite!! need to delete first entry in fHeaders!!
               if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrX))) break;
               continue;
            }//if
         }//if

         if ((err = FillExprArray(seltree[k], entries, arrIdx, arrY))) break;

         // normalize expression levels arrY and fill expression tree
         if ((err = fNormalizer->Calculate(entries, arrX, arrY, arrMsk))) break;

         treename = selname + "." + fNormalizer->GetTitle();
         if ((err = FillExprTree(treename.Data(), entries, arrIdx, arrY))) break;
      }//for_k
//PROBLEM: masktree title must be fSchemeName!! ev. here;
//??      if ((err = this->FillMaskTree(treename.Data(), entries, arrMsk))) goto cleanup;
   } else {
      cerr << "Error: Selector <" << fSelector->GetName()
           << "> or option <" << fSelector->GetOption()
           << "> is not known/applicable." << endl;
      err = errGeneral;
   }//if

// Delete locally created variables
cleanup:
   SafeDelete(refstr);
   delete [] arrMsk;
   delete [] arrIdx;
   delete [] arrY;
   delete [] arrX;

   delete [] reftree;
   delete [] seltree;

   savedir->cd();

// Informing user
   if (XManager::fgVerbose) {
      cout << "   normalization finished." << endl;
   }//if

   return err;
}//Normalize

//______________________________________________________________________________
Int_t XNormedGCSet::ExportExprTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export variables from varlist for all expression trees
   if(kCS) cout << "------XNormedGCSet::ExportExprTrees------" << endl;

   Int_t err = errNoErr;

// Decompose varlist
   Short_t hasUnit   = 0;  //unit name
   Short_t hasTrans  = 0;  //transcript_id
   Short_t hasName   = 0;  //gene name
   Short_t hasSymbol = 0;  //gens symbol
   Short_t hasAccess = 0;  //mRNA accession
   Short_t hasEntrez = 0;  //entrez ID
   Short_t hasChromo = 0;  //chromosome
   Short_t hasStart  = 0;  //start position
   Short_t hasStop   = 0;  //stop position
   Short_t hasStrand = 0;  //strand
   Short_t hasCyto   = 0;  //cytoband
   Short_t hasLevel  = 0;  //expression level

   Short_t hasAnnot  = 0;  //annotation

   Short_t idx = 0;
   if (strcmp(varlist,"*")  == 0) {
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
      hasLevel  = ++idx;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName")     == 0) {hasUnit   = ++idx;}
         if (strcmp(name,"fTranscriptID") == 0) {hasTrans  = ++idx;}
         if (strcmp(name,"fName")         == 0) {hasName   = ++idx;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymbol = ++idx;}
         if (strcmp(name,"fAccession")    == 0) {hasAccess = ++idx;}
         if (strcmp(name,"fEntrezID")     == 0) {hasEntrez = ++idx;}
         if (strcmp(name,"fChromosome")   == 0) {hasChromo = ++idx;}
         if (strcmp(name,"fStart")        == 0) {hasStart  = ++idx;}
         if (strcmp(name,"fStop")         == 0) {hasStop   = ++idx;}
         if (strcmp(name,"fStrand")       == 0) {hasStrand = ++idx;}
         if (strcmp(name,"fCytoBand")     == 0) {hasCyto   = ++idx;}
         if (strcmp(name,"fLevel")        == 0) {hasLevel  = ++idx;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

   // check for presence of at least one annotation variable
   hasAnnot = (hasTrans + hasName  + hasSymbol + hasAccess + hasEntrez
            + hasChromo + hasStart + hasStop   + hasStrand + hasCyto);

// Get trees
   TTree       **tree = new TTree*[n];
   XExpression **expr = new XExpression*[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         expr[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("ExprBranch", &expr[k]);
      }//for_k
   }//if

// Get treeinfo and its option for selection of unittree and anntree
   XTreeInfo *info   = (XTreeInfo*)tree[0]->GetUserInfo()->At(0);
   Option_t  *option = info->GetOption();

   Int_t type = eTRANSCRIPT;
   if      (strcmp(option, "exon")     == 0) type = eEXONTYPE; 
   else if (strcmp(option, "probeset") == 0) type = ePROBESET; 

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName.Data() << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      return errAbort;
   }//if

// Get chip from scheme file
   if (fSchemeFile == 0) return errGetScheme;
   fSchemeFile->cd();
   XFolder *schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      return fManager->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   XGeneChip *chip = (XGeneChip*)schemes->FindObject(fSchemeName, kTRUE);
   if (chip == 0) return errAbort;

// Get unit tree for scheme (needed for annotation)
   XGCUnit *unit   = 0;
//?   XExonUnit *unit = 0;
   TTree *unittree = GetUnitTree(chip, type);
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

   Int_t numunits = (Int_t)(unittree->GetEntries());

// Get annotation tree for scheme
   XTransAnnotation *annot = 0;
   TTree *anntree = this->GetAnnotationTree(chip, type);
   Int_t numannot = 0;

// Create hash table to store unit names from anntree
   THashTable *htable = 0;
   if (anntree) {
      anntree->SetBranchAddress("AnnBranch", &annot);
      numannot = (Int_t)(anntree->GetEntries());

      if (!(htable = new THashTable(2*numannot))) return errInitMemory;
      htable = FillHashTable(htable, anntree, annot, type);
   } else {
      if (hasAnnot) {
         cout << "Warning: Missing annotation, gene info not exported." << endl;
      }//if
      hasAnnot = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   for (Short_t j=1; j<=idx; j++) {
      if (hasUnit == j) output << sep << "UnitName";

      if (hasAnnot > 0) {
//         if (hasTrans  == j) output << sep << "ProbesetID";
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

      if (hasLevel == j) {
         if (n == 1) {
            output << sep << "LEVEL";
         } else {
            for (Int_t k=0; k<n; k++) {
               output << sep << (names[k] + "_LEVEL").Data();
            }//for_k
         }//if
      }//if
   }//for_j
   output << endl;

// Loop over tree entries and trees
   XIdxString *idxstr = 0;
   Int_t index = 0;
   Int_t cnt   = 0;
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      tree[0]->GetEntry(i);

      Int_t unitID = expr[0]->GetUnitID();
      output << unitID;

      unittree->GetEntry(index++);
      while (unitID != unit->GetUnitID()) {
         if (index == numunits) {
           cerr << "Error: UnitID <" << unitID << "> not found." << endl;
           err = errAbort; goto cleanup;
         }//if

         unittree->GetEntry(index++);
      }//while

      for (Short_t j=1; j<=idx; j++) {
         // export unitname
         if (hasUnit == j) {
            output << sep << unit->GetUnitName();
         }//if

         // export annotation
         if (hasAnnot > 0) {
            idxstr = FindUnitID(htable, unit);
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
         if (hasLevel == j) {
            for (Int_t k=0; k<n; k++) {
               tree[k]->GetEntry(i);
               output << sep << expr[k]->GetLevel();
            }//for_k
         }//if
      }//for_j
      output << endl;

      cnt++;
      if (XManager::fgVerbose && cnt%10000 == 0) {
         cout << "<" << cnt << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << cnt << "> of " << "<" << nentries << "> records exported." << endl;
   }//if

//Cleanup
cleanup:
   if (htable) {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(schemes);

   // remove trees from RAM
   for (Int_t k=0; k<n; k++) {
      SafeDelete(expr[k]);
      tree[k]->ResetBranchAddress(tree[k]->GetBranch("ExprBranch"));
      SafeDelete(tree[k]);
   }//for_k

   delete [] expr;
   delete [] tree;

   SafeDelete(annot);
   anntree->ResetBranchAddress(anntree->GetBranch("AnnBranch"));
   SafeDelete(anntree);

   SafeDelete(unit);
   unittree->ResetBranchAddress(unittree->GetBranch("IdxBranch"));
   SafeDelete(unittree);

   return errNoErr;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XNormedGCSet::ExportMaskTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export variables from varlist for all mask trees
   if(kCS) cout << "------XNormedGCSet::ExportMaskTrees------" << endl;

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

// Get scheme name (also for alternative CDFs)
   if (strcmp(fSchemeName.Data(), "") == 0) {
      fSchemeName = tree[0]->GetTitle();
   } else if (!fSchemeName.Contains(tree[0]->GetTitle())) {
      cerr << "Error: Scheme <" << fSchemeName.Data() << "> is not derived from <"
           << tree[0]->GetTitle() << ">." << endl;
      hasUnit  = kFALSE;
   }//if

// If scheme file is present get unit tree
   Int_t    numgenes  = 0;
   Int_t    numctrls  = 0;
   XFolder  *schemes  = 0;
   XGCUnit  *unit     = 0;
   TTree    *unittree = 0; 
   if (hasUnit && fSchemeFile) {
   // Get chip from scheme file
      fSchemeFile->cd();
      schemes = (XFolder*)(fSchemeFile->Get(kContent));
      if (!schemes) {
         return fManager->HandleError(errMissingContent, "Scheme", kContent);
      }//if

      XDNAChip *chip = 0;
      chip = (XDNAChip*)schemes->FindObject(fSchemeName, kTRUE);
      if (chip) {
         numgenes = chip->GetNumGenes();
         numctrls = chip->GetNumControls();

         if (!fSchemeFile->cd(fSchemeName)) return errGetDir;

      // Get unit tree for scheme
         if (hasUnit) {
            unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
            if (unittree == 0) return errGetTree;
            unittree->SetBranchAddress("IdxBranch", &unit);
         }//if
      } else {
         cerr << "Error: Could not find scheme <" << fSchemeName.Data() << ">." << endl;
         hasUnit = kFALSE;
      }//if
   } else if (hasUnit) {
      cerr << "Warning: Missing scheme file, unit names not exported." << endl;
      hasUnit = kFALSE;
   }//if

////////////
//better: see XExonProcesSet!!!
////////////
// Check if tree entries is equal to number of genes in unittree
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   if (hasUnit && (entries != numgenes)) {
      if (XManager::fgVerbose) {
         cout << "Note: Number of tree entries <" << entries  
              << "> is not equal to number of units <" << numgenes << ">" << endl;
      }//if
      hasUnit = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)    output << sep << "UNIT_NAME";
   if (n == 1) {
      if (hasFlag) output << sep << "FLAG";
   } else {
      for (Int_t k=0; k<n; k++) {
         if (hasFlag)  output << sep << (names[k] + "_FLAG").Data();
      }//for_k
   }//if
   output << endl;

// Loop over tree entries and tree branches
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         // output unitID for first tree only
         if (k == 0) {
            Int_t unitID = mskunit->GetUnitID();
            output << unitID;

            if (hasUnit) {
               unittree->GetEntry(unitID + numctrls); //unit names for genes only
               output << sep << unit->GetUnitName();
            }//if
         }//if

         // output mask for all trees
         if (hasFlag) output << sep << mask[k]->GetFlag();
      }//for_k
      output << endl;
   }//for_i

//Cleanup
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

   return errNoErr;
}//ExportMaskTrees


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedGenomeSet                                                     //
//                                                                      //
// Class for normalization of GenomeChip oligonucleotide arrays         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormedGenomeSet::XNormedGenomeSet()
                 :XNormedGCSet()
{
   // Default NormedGenomeSet constructor
   if(kCS) cout << "---XNormedGenomeSet::XNormedGenomeSet(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XNormedGenomeSet::XNormedGenomeSet(const char *name, const char *type)
                 :XNormedGCSet(name, type)
{
   // Normal NormedGenomeSet constructor
   if(kCS) cout << "---XNormedGenomeSet::XNormedGenomeSet------" << endl;

}//Constructor

//______________________________________________________________________________
XNormedGenomeSet::~XNormedGenomeSet()
{
   // NormedGenomeSet destructor
   if(kCS) cout << "---XNormedGenomeSet::~XNormedGenomeSet------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedExonSet                                                       //
//                                                                      //
// Class for normalization of ExonChip oligonucleotide arrays           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormedExonSet::XNormedExonSet()
               :XNormedGenomeSet()
{
   // Default NormedExonChipSet constructor
   if(kCS) cout << "---XNormedExonSet::XNormedExonSet(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XNormedExonSet::XNormedExonSet(const char *name, const char *type)
               :XNormedGenomeSet(name, type)
{
   // Normal NormedExonChipSet constructor
   if(kCS) cout << "---XNormedExonSet::XNormedExonSet------" << endl;

}//Constructor

//______________________________________________________________________________
XNormedExonSet::~XNormedExonSet()
{
   // NormedExonChipSet destructor
   if(kCS) cout << "---XNormedExonSet::~XNormedExonSet------" << endl;

}//Destructor
