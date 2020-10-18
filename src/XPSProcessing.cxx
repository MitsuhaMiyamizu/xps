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
* Oct 2002 - Add method InitApprox()
* Dec 2002 - Inherit classes XPreProcessManager, XNormationManager from XProcessManager
*            and created inherited classes from corresponding supporting classes
*          - Replace InitXXX() with InitAlgorithm()
*          - Completely reorganize source file structure
* Aug 2004 - Reorganize tree extensions.
* Nov 2004 - Add classes XExpressionTreeInfo etc and move XAlgorithm to XPSBase.h.
* May 2005 - Replace XContent with XFolder.
* Jun 2005 - Add methods to support database.
* Mar 2006 - Add support for alternative schemes, i.e. CDF files.
* Nov 2007 - Save fSchemeName of class XProcesSet, increase ClassDef to 2.
*          - XProcessManager now inherits from XManager and XProjectHandler
*          - database methods ProjectInfo() etc are moved to XProjectHandler
* Aug 2008 - Add summarization algorithms FARMS and DFW, classes XFARMS, XDFW
* Dec 2010 - Add quality control algorithms 
*
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include "TBranch.h"
#include "THashTable.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TSystem.h"

#include "XPSProcessing.h"
#include "XPSDataTypes.h"
#include "TStat.h"

// Tree extensions and types for preprocessing:
// background trees
const char *kExtenBgrd[] = { "sbg", "wbg", "rbg", "gbg", ""};
const char *kTypeBgrd[]  = { "sector",
                             "weightedsector",
                             "rma",
                             "gccontent",
                             ""};

// intensity trees after background correction
const char *kExtenIntn[] = { "int", ""};
const char *kTypeIntn[]  = { "intensity",
                             ""};

// residual trees
const char *kExtenResd[] = { "res", ""};
const char *kTypeResd[]  = { "residual",
                             ""};

// border trees
const char *kExtenBord[] = { "brd", ""};
const char *kTypeBord[]  = { "border",
                             ""};
/////////////
//EV BETTER: intensity for each bgrd algorithm
//const char *kExtenIntn[] = { "sin", "win", "rin", ""};
//const char *kTypeBgrd[]  = { "sector intensity",
//                             "weighted sector intensity",
//                             "rma intensity",
//                             ""};
/////////////

// intensity trees after low-level normalization of bg-corrected intensities
const char *kExtenCNrm[] = { "cmn", "cmd", "clw", "css", "cqu", ""};
const char *kTypeCNrm[]  = { "mean",
                             "median",
                             "lowess",
                             "supsmu",
                             "quantile",
                             ""};

// expression trees after summarization
const char *kExtenExpr[] = { "amn", "gmn", "wmn", "gcm", "wdf", "adf", "tbw",
                             "mdp", "frm", "dfw", "fir", ""};
const char *kTypeExpr[]  = { "arithmeticmean",
                             "geometricmean",
                             "weightedmean",
                             "gccorrectedmean",
                             "weighteddiff",
                             "avgdiff",
                             "tukeybiweight",
                             "medianpolish",
                             "farms",
                             "dfw",
                             "firma",
                             ""};

// detection call trees
const char *kExtenCall[] = { "mdf", "dc5", "dc4", "dab", "ini", ""};
const char *kTypeCall[]  = { "meandifferencecall",
                             "detectioncall",
                             "mas4call",
                             "dabgcall",
                             "inicall",
                             ""};

// quality trees
//const char *kExtenQual[] = { "mdp", "plm", ""};
const char *kExtenQual[] = { "rlm", "plm", ""};
const char *kTypeQual[]  = { "rlm",
                             "plm",
                             ""};
/*???
const char *kExtenQual[] = { "rma", "rmar", "rmaa", "rman", "plm", "plmr", "plma", "plmn", ""};
const char *kTypeQual[]  = { "rma",    //all
                             "rmaraw",
                             "rmaadjusted",
                             "rmanormalized",
                             "plm",    //all
                             "plmraw",
                             "plmadjusted",
                             "plmnormalized",
                             ""};
*/

//////////////
//?? Problem?? "mdp" in kExtenExpr AND in kExtenNorm
//////////////
// normalized trees after high-level normalization of expression trees
const char *kExtenNorm[] = { "tmn", "med", "ksm", "low", "sup", "qua", "mdp", ""};
const char *kTypeNorm[]  = { "mean",
                             "median",
                             "ksmooth",
                             "lowess",
                             "supsmu",
                             "quantile",
                             "medianpolish",  //??allowed?? since in kTypeExpr
                             ""};

// selector trees
const char *kExtenSlct[] = { "sel", "rnk", "prb", "idx", "usr", ""};
const char *kTypeSlct[]  = { "default",
                             "rank",
                             "probe",
                             "unit",
                             "user",
                             ""};

// filter trees
const char *kExtenFltr[] = { "pfr", "ufr", "mfr", ""};
const char *kTypeFltr[]  = { "prefilter",
                             "unifilter",
                             "multifilter",
                             ""};

// pre-filter trees
const char *kPreFltr[]   = { "mad",
                             "cv",
                             "variance",
                             "difference",
                             "ratio",
                             "gap",
                             "lowerthreshold",
                             "upperthreshold",
                             "quantile",
                             "entropy",
                             "call",
                             ""};

// univariate filter trees
const char *kUniFltr[]   = { "foldchange",
                             "statistic",
                             "pvalue",
                             "pchance",
                             "padjust",
                             "call",
                             ""};

// multivariate filter trees
const char *kMultiFltr[]   = { "foldchange",
                               "statistic",
                               "pvalue",
                               "call",
                               ""};

// univariate test trees
const char *kExtenUTst[] = { "uvt", "stt", "wil", "var", ""};
const char *kTypeUTst[]  = { "normaltest",
                             "ttest",
                             "wilcoxon",
                             "variance",
                             ""};

// multivariate test trees
const char *kExtenMTst[] = { "aov", "kwt", ""};
const char *kTypeMTst[]  = { "anova",
                             "kruskalwallis",
                             ""};

// cluster analysis trees
const char *kExtenClus[] = { "hcl", "kmn", "som", "pam", "cla", ""};
const char *kTypeClus[]  = { "hierarchicalcluster",
                             "kmeans",
                             "som",
                             "pam",
                             "clara",
                             ""};

// regression trees
const char *kExtenRgrs[] = { "xxx", "yyy", ""};
const char *kTypeRgrs[]  = { "xxxxxx",
                             "yyyyyy",
                             ""};

// processing methods
const char *kProcessMethod[] = { "preprocess",
                                 "rma",
                                 "mas4",
                                 "mas5",
                                 "adjustbgrd",
                                 "normalize",
                                 "express",
                                 "detectcall",
                                 "qualify",
                                 ""};

// options for detection call
const char *kCallOption[] = { "raw",
                              "adjusted",
                              "normalized",
                              ""};

// options for quality control
const char *kQualOption[] = { "raw",
                              "adjusted",
                              "normalized",
                              "all", //??
                              ""};

// name of reference tree
const char*  kReference = "Reference";

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XProcessManager);
ClassImp(XExpressionTreeInfo);
ClassImp(XSelectionTreeInfo);
ClassImp(XProcesSetting);
ClassImp(XProcesSet);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProcessManager                                                      //
//                                                                      //
// Base class for processing of microarray data                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProcessManager::XProcessManager()
                :XManager(), XProjectHandler()
{
   // Default ProcessManager constructor
   if(kCS) cout << "---XProcessManager::XProcessManager(default)------" << endl;

   fSchemeFile = 0;
   fDataFile   = 0;
   fSchemes    = 0;
   fData       = 0;
   fIsSchemeOwner = kFALSE;
   fIsDataOwner   = kFALSE;
}//Constructor

//______________________________________________________________________________
XProcessManager::XProcessManager(const char *name, const char *arraytype, Int_t verbose)
                :XManager(name, arraytype, verbose), XProjectHandler(name, arraytype)
{
   // Normal ProcessManager constructor
   if(kCS) cout << "---XProcessManager::XProcessManager------" << endl;

   fSchemeFile = 0;
   fDataFile   = 0;
   fSchemes    = 0;
   fData       = 0;
   fIsSchemeOwner = kFALSE;
   fIsDataOwner   = kFALSE;
}//Constructor

//______________________________________________________________________________
XProcessManager::~XProcessManager()
{
   // ProcessManager destructor
   if(kCS) cout << "---XProcessManager::~XProcessManager------" << endl;

//   this->Close(); //not here but in base destructor
}//Destructor

//______________________________________________________________________________
Int_t XProcessManager::AddTree(const char *setname, const char *intree,
                       Int_t treeid, Option_t *option)
{
   // Select intree for processing, i.e. add to tree set setname
   // Argument intree can be: /path/filename.root/treename.exten
   // For treename is "*.exten", all trees with exten are added from root file
   // The following  options can be set:
   // option = "reference": intree is used as reference tree for processing
   // option = "baseline":  intree is used as baseline for calculating ratios
   // option = "baseref":   intree is used both as reference and as baseline
   // If more than one intree is set to reference and/or baseline
   // than the mean/median of selected intrees is used as reference 
   // and/or mean/median of selected AND processed trees is used as baseline
   // Note: If argument intree of the first call to AddTree() does not contain
   //       filename.root then OpenData() must be called first.
   // Note: For selection of tree friends set: treeid > 0 in order to prevent
   //       the addition of tree friend headers to list fTrees of XTreeSet!
   if(kCS) cout << "------XProcessManager::AddTree------" << endl;

   if (fAbort) return errAbort;

   Int_t err = errNoErr;

   // try to extract missing fDataFile from intree
   if (!fDataFile) {
      if (strstr(intree,".root")) {
         err = this->OpenData(GetROOTName(intree) + ".root");
      } else {
         err = errAbort;
      }//if
   }//if

   if (err != errNoErr) {
      cerr << "Error: Could not find data file. Need to call OpenData() first."
           << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   return XManager::AddTree(setname, intree, treeid, option);
}//AddTree

//______________________________________________________________________________
void XProcessManager::Close(Option_t *option)
{
   // Close root file
   if(kCS) cout << "------XProcessManager::Close------" << endl;

   // check if fDataFile and fFile are identical
   Int_t hasSamePointer = (fDataFile == fFile) ? 1 : 0;

   // allow Close/Save in XPSApp w/o calling destructor
   XManager::Close(option);

   if (hasSamePointer) fDataFile = 0;

   this->CloseData();
   this->CloseSchemes();
}//Close

//______________________________________________________________________________
Int_t XProcessManager::HandleError(Int_t err, const char *name1, const char *name2)
{
   // Handle error messages
   if(kCS) cout << "------XProcessManager::HandleError------" << endl;

   switch (err) {
      case errGetScheme:
         cerr << "Error: Could not get scheme <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errSchemeDerived:
         cerr << "Error: Scheme <" << name1 << "> is not derived from <"
              << name2 << ">." << endl;
         fAbort = kTRUE;
         break;

      default:
         XDataManager manager;
         err = manager.HandleError(err, name1, name2);
//         err = XManager::HandleError(err, name1, name2);
         break;
   }//switch

   return err;
}//HandleError

//______________________________________________________________________________
Int_t XProcessManager::SetBaseLine(const char *intree, Option_t *option,
                       Double_t trim)
{
   // Set intree to be used as baseline for calculating ratios, i.e. fold-change
   // If more than one intree is set to reference and/or baseline
   // than the mean/median of seleceted intrees is used as reference 
   // For option = "mean", trim can be used to calculate a trimmed mean
   // To change only option, call with intree = "".
   if(kCS) cout << "------XProcessManager::SetBaseLine------" << endl;

   if (fAbort) return errAbort;
   if (!fTreeSet) {HandleError(errInitTreeSet); return errAbort;}

// To change option only
   Int_t err = errNoErr;
   if (strcmp(intree, "") == 0) {
      err = ((XProcesSet*)fTreeSet)->SetBaseLine(0, option, trim);
      return err;
   }//if

// Extract tree name from intree
   TString inname = Path2Name(intree, dSEP, "");
   if (strstr(inname.Data(), ".root")) {
      inname = "";
   }//if

// Extract root filename from intree
   TFile * file = 0;
   TString filename = "";
   Bool_t  isOwner  = kFALSE;
   if (strstr(intree, ".root")) {
      filename = Path2Name(intree, "", ".") + ".root";
      if (strcmp(filename.Data(),gDirectory->GetName()) != 0) {
         // open file
         file = this->OpenFile(filename.Data(), "READ", isOwner);
         if (!file) return perrOpenFile;
         file->cd();
      }//if
   } else {
      filename = gDirectory->GetName();
   }//if

   TDirectory *savedir = gDirectory;

// Get name of treeset and change directory
   TString sname  = "";
   if (strstr(intree, ".root")) {
      TString substr = SubString(intree, '.', sSEP, kFALSE);
      if (substr) sname = Path2Name(substr.Data(), dSEP, "");
      if (sname.Contains("root")) sname = "";
   } else if (strstr(intree, dSEP)) {
      sname = Path2Name(intree, "", dSEP);
   }//if
   if (!gDirectory->cd(sname)) return HandleError(errGetDir, sname);

// Add trees to treeset
   TString name  = Path2Name(intree, dSEP, ".");
   TString exten = Path2Name(intree, ".", "");
   if (strcmp(name.Data(), "*") == 0) {
   // Loop over all trees with extension exten
      Int_t numtrees = 0;
      TKey *key = 0;
      TIter next(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TString xten  = Path2Name(key->GetName(), ".", ";");
         TString kname = Path2Name(key->GetName(), "", ".");
         if (strcmp(xten.Data(), exten) == 0) {
            TTree* tree = (TTree*)gDirectory->Get(key->GetName());
            err = ((XProcesSet*)fTreeSet)->SetBaseLine(tree, option, trim);
            numtrees++;
         }//if
      }//while
   } else {
   // Add tree with name inname
      TTree* tree = (TTree*)gDirectory->Get(inname);
      err = ((XProcesSet*)fTreeSet)->SetBaseLine(tree, option, trim);
   }//if

   savedir->cd();

   if (err != errNoErr) fAbort = kTRUE;
   return err;
}//SetBaseLine

//______________________________________________________________________________
Int_t XProcessManager::SetReference(const char *intree, Option_t *option,
                       Double_t trim)
{
   // Set intree(s) to be used as reference tree for further processing
   // For "*.exten" all trees with extension exten are used as reference
   // If option is set to "mean" or "median" then the mean/median of all
   // tree levels is used as reference for processing
   // For option = "mean", trim can be used to calculate a trimmed mean
   // To change only option, call with intree = "".
   if(kCS) cout << "------XProcessManager::SetReference------" << endl;

   if (fAbort) return errAbort;
   if (!fTreeSet) {HandleError(errInitTreeSet); return errAbort;}

// To change option only
   Int_t err = errNoErr;
   if (strcmp(intree, "") == 0) {
      err = ((XProcesSet*)fTreeSet)->SetReference(0, option, trim);
      return err;
   }//if

// Extract tree name from intree
   TString inname = Path2Name(intree, dSEP, "");
   if (strstr(inname.Data(), ".root")) {
      inname = "";
   }//if

// Extract root filename from intree
   TFile * file = 0;
   TString filename = "";
   Bool_t  isOwner  = kFALSE;
   if (strstr(intree, ".root")) {
      filename = Path2Name(intree, "", ".") + ".root";
      if (strcmp(filename.Data(),gDirectory->GetName()) != 0) {
         // open file
         file = this->OpenFile(filename.Data(), "READ", isOwner);
         if (!file) return perrOpenFile;
         file->cd();
      }//if
   } else {
      filename = gDirectory->GetName();
   }//if

   TDirectory *savedir = gDirectory;

// Get name of treeset and change directory
   TString sname  = "";
   if (strstr(intree, ".root")) {
      TString substr = SubString(intree, '.', sSEP, kFALSE);
      if (substr) sname = Path2Name(substr.Data(), dSEP, "");
      if (sname.Contains("root")) sname = "";
   } else if (strstr(intree, dSEP)) {
      sname = Path2Name(intree, "", dSEP);
   }//if
   if (!gDirectory->cd(sname)) return HandleError(errGetDir, sname);

// Add trees to treeset
   TString name  = Path2Name(intree, dSEP, ".");
   TString exten = Path2Name(intree, ".", "");
   if (strcmp(name.Data(), "*") == 0) {
   // Loop over all trees with extension exten
      Int_t numtrees = 0;
      TKey *key = 0;
      TIter next(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TString xten  = Path2Name(key->GetName(), ".", ";");
         TString kname = Path2Name(key->GetName(), "", ".");
         if (strcmp(xten.Data(), exten) == 0) {
            TTree* tree = (TTree*)gDirectory->Get(key->GetName());
            err = ((XProcesSet*)fTreeSet)->SetReference(tree, option, trim);
            numtrees++;
         }//if
      }//while
   } else {
   // Add tree with name inname
      TTree* tree = (TTree*)gDirectory->Get(inname);
      err = ((XProcesSet*)fTreeSet)->SetReference(tree, option, trim);
   }//if

   savedir->cd();

   if (err != errNoErr) fAbort = kTRUE;
   return err;
}//SetReference

//______________________________________________________________________________
Int_t XProcessManager::InitData(TFile *datafile, Bool_t isOwner)
{
   // Initialize processmanager with external root datafile
   // Set isOwner = kTRUE if datafile should be deleted by XProcessManager
   if(kCS) cout << "------XProcessManager::InitData------" << endl;

   if (fAbort) return errAbort;
   TDirectory *savedir = gDirectory;

   // check if datafile is already open
   if (this->IsOpen(fDataFile, datafile->GetName())) {
      // only one data file can be open
      if (XManager::fgVerbose) {
         cout << "Closing existing data file <" << fDataFile->GetName()
              << ">..." << endl;
      }//if
      this->CloseData();
   }//if

   fDataFile = datafile;
   if (!fDataFile) {fAbort = kTRUE; return errGetFile;}
   fIsDataOwner = isOwner;

   fDataFile->cd();
   fData = (XFolder*)(fDataFile->Get(kContent));
   if (!fData) {
      return HandleError(errMissingContent, "Data", kContent);
   }//if

   if (fSetting) {
      ((XProcesSetting*)fSetting)->SetDataFile(fDataFile);
   }//if

   savedir->cd();
   return errNoErr;
}//InitData

//______________________________________________________________________________
Int_t XProcessManager::InitSchemes(TFile *schemefile, Bool_t isOwner,
                       const char *schemename, const char *schemetype)
{
   // Initialize processmanager with external root schemefile
   // Set isOwner = kTRUE if schemefile should be deleted by XProcessManager
   // Optional select alternative scheme "schemename"
   if(kCS) cout << "------XProcessManager::InitSchemes------" << endl;

   if (fAbort) return errAbort;
   TDirectory *savedir = gDirectory;

   // check if schemefile is already open
   if (this->IsOpen(fSchemeFile, schemefile->GetName())) {
      // only one scheme file can be open
      if (XManager::fgVerbose) {
         cout << "Closing existing scheme file <" << fSchemeFile->GetName()
              << ">..." << endl;
      }//if
      this->CloseSchemes();
   }//if

   fSchemeFile = schemefile;
   if (!fSchemeFile) {fAbort = kTRUE; return errGetFile;}
   fIsSchemeOwner = isOwner;

   fSchemeFile->cd();
   fSchemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!fSchemes) {
      return HandleError(errMissingContent, "Scheme", kContent);
   }//if

   if (fSetting) {
      ((XProcesSetting*)fSetting)->SetSchemeFile(fSchemeFile);
      ((XProcesSetting*)fSetting)->SetSchemeName(schemename);
      ((XProcesSetting*)fSetting)->SetSchemeType(schemetype);
   }//if

   savedir->cd();
   return errNoErr;
}//InitSchemes

//______________________________________________________________________________
Int_t XProcessManager::OpenData(const char *fullname, Option_t *option)
{
   // Open root data file
   if(kCS) cout << "------XProcessManager::OpenData------" << endl;

   if (fAbort) return errAbort;

   // check if datafile is already open
   if (this->IsOpen(fDataFile, fullname)) {
      // only one scheme file can be open
      if (XManager::fgVerbose) {
         cout << "Closing existing data file <" << fDataFile->GetName()
              << ">..." << endl;
      }//if
      this->CloseData(); //deletes also fData!
   }//if

   Bool_t isOwner = kFALSE;
   fDataFile = OpenFile(fullname, option, isOwner);
   if (!fDataFile) {fAbort = kTRUE; return errCreateFile;}
   // assure that manager remains owner if it calls OpenFile multiple times
   if (!fIsDataOwner) fIsDataOwner = isOwner;

   fDataFile->cd();
   fData = (XFolder*)(fDataFile->Get(kContent));
   if (!fData) {
      return HandleError(errMissingContent, "Data", kContent);
   }//if

   if (fSetting) {
      ((XProcesSetting*)fSetting)->SetDataFile(fDataFile);
   }//if

   return errNoErr;
}//OpenData

//______________________________________________________________________________
Int_t XProcessManager::OpenSchemes(const char *fullname, const char *schemename,
                       const char *schemetype)
{
   // Open root scheme file using full path 
   // Optional select alternative scheme "schemename"
   if(kCS) cout << "------XProcessManager::OpenSchemes------" << endl;

   if (fAbort) return errAbort;
   TDirectory *savedir = gDirectory;

   // check if schemefile is already open
   if (this->IsOpen(fSchemeFile, fullname)) {
      // only one scheme file can be open
      if (XManager::fgVerbose) {
         cout << "Closing existing scheme file <" << fSchemeFile->GetName()
              << ">..." << endl;
      }//if
      this->CloseSchemes();
   }//if

   Bool_t isOwner = kFALSE;
   fSchemeFile = OpenFile(fullname, "READ", isOwner);
   if (!fSchemeFile) {fAbort = kTRUE; return errCreateFile;}
   // assure that manager remains owner if it calls OpenFile multiple times
   if (!fIsSchemeOwner) fIsSchemeOwner = isOwner;

   fSchemeFile->cd();
   fSchemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!fSchemes) {
      return HandleError(errMissingContent, "Scheme", kContent);
   }//if

   if (fSetting) {
      ((XProcesSetting*)fSetting)->SetSchemeFile(fSchemeFile);
      ((XProcesSetting*)fSetting)->SetSchemeName(schemename);
      ((XProcesSetting*)fSetting)->SetSchemeType(schemetype);
   }//if

   savedir->cd();
   return errNoErr;
}//OpenSchemes

//______________________________________________________________________________
void XProcessManager::CloseData()
{
   // Close root data file
   if(kCS) cout << "------XProcessManager::CloseData------" << endl;

   if (fDataFile) {
   // Write content to data file only if new file or file is updated
      if (fData && (strcmp(fDataFile->GetOption(), "READ") != 0)) {
         fDataFile->cd();
         fData->Write("", TObject::kWriteDelete);
      }//if

      if (fIsDataOwner) SafeDelete(fDataFile);
      fDataFile = 0;  //if not isOwner!
   }//if

   SafeDelete(fData);
}//CloseData

//______________________________________________________________________________
void XProcessManager::CloseSchemes()
{
   // Close scheme file
   if(kCS) cout << "------XProcessManager::CloseSchemes------" << endl;

   SafeDelete(fSchemes);

   if (fIsSchemeOwner) SafeDelete(fSchemeFile);
   fSchemeFile = 0;  //if not isOwner!
}//CloseSchemes

//______________________________________________________________________________
Int_t XProcessManager::BeginTransaction(const char *name)
{
   // Begin transaction
   if(kCS) cout << "------XProcessManager::BeginTransaction------" << endl;

   if (fAbort) return errAbort;
   
   return errNoErr;
}//BeginTransaction

//______________________________________________________________________________
Int_t XProcessManager::CommitTransaction()
{
   // Commit transaction
   if(kCS) cout << "------XProcessManager::CommitTransaction------" << endl;

   if (fAbort) return errAbort;

   if (fList && fList->GetSize() > 0) {
      for (Int_t i=0; i<fList->GetSize(); i++) {
         XDataTypeInfo *info = (XDataTypeInfo*)(fList->At(i));

         if (strcmp(info->ClassName(), "XDatasetInfo") == 0) {
            // set data type for data contained in dataset
            info->SetDataType(fDataType);

            if (info->Replace() == kTRUE) {
               XDatasetInfo *oldinfo = 0;
               oldinfo = (XDatasetInfo*)fContent->FindObject("Dataset", "XDatasetInfo");
               if (oldinfo == 0) {
                  fContent->Add(info);
                  return errNoErr;
               }//if

               TString oldname = oldinfo->GetDatasetName();
               TString newname = ((XDatasetInfo*)info)->GetDatasetName();
               if (strcmp(oldname.Data(), newname.Data()) != 0) {
                  cout << "Warning: Currently it is not possible to change dataset name <"
                       << oldname.Data() << "> to dataset name <" << newname.Data() <<">."
                       << endl;
                  ((XDatasetInfo*)info)->SetDatasetName(oldname);
               }//if

               fContent->Remove(oldinfo);
            }//if
         } else if (strcmp(info->ClassName(), "XHybridizationList") == 0) {
            //check for existance of info in case of update file and/or replace info
            XHybridizationList *oldlist = 0;
            oldlist = (XHybridizationList*)fContent->FindObject(info->GetName(), info->ClassName());

            if (oldlist) {
               XHybridizationList *newlist = (XHybridizationList*)info;

               XHybInfo *oldhyb = 0;
               XHybInfo *newhyb = 0;
               for (Int_t i=(oldlist->GetSize()-1); i>-1; i--) {
                  oldhyb = (XHybInfo*)(oldlist->At(i));
                  newhyb = (XHybInfo*)(newlist->FindDataTypeInfo(oldhyb->GetHybName()));

                  if (newhyb == 0) { //add only if not replaced by new info
                     newlist->AddAt(oldhyb, 0);
                  }//if
               }//for_i

               fContent->Remove(oldlist);
            }//if
         } else if (strcmp(info->ClassName(), "XTreatmentList") == 0) {
            //check for existance of info in case of update file and/or replace info
            XTreatmentList *oldlist = 0;
            oldlist = (XTreatmentList*)fContent->FindObject(info->GetName(), info->ClassName());

            if (oldlist) {
               XTreatmentList *newlist = (XTreatmentList*)info;

               XTreatmentInfo *oldtreat = 0;
               XTreatmentInfo *newtreat = 0;
               for (Int_t i=(oldlist->GetSize()-1); i>-1; i--) {
                  oldtreat = (XTreatmentInfo*)(oldlist->At(i));
                  newtreat = (XTreatmentInfo*)(newlist->FindDataTypeInfo(oldtreat->GetTreatmentName()));

                  if (newtreat == 0) { //add only if not replaced by new info
                     newlist->AddAt(oldtreat, 0);
                  }//if
               }//for_i

               fContent->Remove(oldlist);
            }//if
         } else if (info->Replace() == kTRUE) {
            fContent->Remove(fContent->FindObject(info->GetName(), info->ClassName()));
         }//if

         fContent->Add(info);
      }//for_i
   } else {
      cerr << "Error: Could not add DataTypes to Content!" << endl;
   }//if
   
   return errNoErr;
}//CommitTransaction

//______________________________________________________________________________
XSetting *XProcessManager::NewSetting(const char *type, const char *infile)
{
   // Create new process setting for type
   if(kCS) cout << "------XProcessManager::NewSetting------" << endl;

   XProcesSetting *setting = new XProcesSetting(type, infile);

   if (setting && fSchemeFile) {
      setting->SetSchemeFile(fSchemeFile);
   }//if

   return setting;
}//NewSetting

//______________________________________________________________________________
XPlot *XProcessManager::NewPlotter(const char *name, const char *title)
{
   // Create new plotter
   if(kCS) cout << "------XProcessManager::NewPlotter------" << endl;

   XPlot *plotter = 0;

   if (strcmp(title,"GeneChip") == 0) {
      plotter = new XGeneChipPlot(name, title);
   } else if (strcmp(title,"GenePix") == 0) {
      plotter = new XGenePixPlot(name, title);
   } else {
      cerr << "Error: Chip type <" << title << "> not known" << endl;
   }//if

   return plotter;
}//NewPlotter

//______________________________________________________________________________
Int_t XProcessManager::DeleteTreeSetInfo(const char *name)
{
   // Delete additional information of treeset "name"
   if(kCS) cout << "------XProcessManager::DeleteTreeSetInfo------" << endl;

   XDatasetInfo *datasetinfo = 0;
   datasetinfo = (XDatasetInfo*)(fContent->FindObject("Dataset", name, "XDatasetInfo"));
   if (datasetinfo) fContent->Remove(datasetinfo);

   return errNoErr;
}//DeleteTreeSetInfo


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExpressionTreeInfo                                                  //
//                                                                      //
// Class containing info about expression tree stored in fUserInfo      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExpressionTreeInfo::XExpressionTreeInfo() 
                    :XTreeInfo()
{
   // Default ExpressionTreeInfo constructor
   if(kCS) cout << "---XExpressionTreeInfo::XExpressionTreeInfo(default)------" << endl;

   fNUnits     = 0;
   fMinLevel   = 0;
   fMaxLevel   = 0;
   fNQuantiles = 0;
   fQuantiles  = 0;
   fLevelQuant = 0;
}//Constructor

//______________________________________________________________________________
XExpressionTreeInfo::XExpressionTreeInfo(const char *name, const char *title) 
                    :XTreeInfo(name, title)
{
   // Normal ExpressionTreeInfo constructor
   if(kCS) cout << "---XExpressionTreeInfo::XExpressionTreeInfo------" << endl;

   fNUnits     = 0;
   fMinLevel   = 0;
   fMaxLevel   = 0;
   fNQuantiles = 7;
   fQuantiles  = new Double_t[fNQuantiles];
   fLevelQuant = new Double_t[fNQuantiles];
}//Constructor

//______________________________________________________________________________
XExpressionTreeInfo::~XExpressionTreeInfo()
{
   // ExpressionTreeInfo destructor
   if(kCS) cout << "---XExpressionTreeInfo::~XExpressionTreeInfo------" << endl;

   if (fLevelQuant) {delete [] fLevelQuant; fLevelQuant = 0;}
   if (fQuantiles)  {delete [] fQuantiles;  fQuantiles  = 0;}
}//Destructor

//______________________________________________________________________________
void XExpressionTreeInfo::AddUserInfo(Int_t nunits, Double_t min, Double_t max)
{
   // Add user info from tree set
   if(kCS) cout << "------XExpressionTreeInfo::AddUserInfo------" << endl;

   fNUnits   = nunits;
   fMinLevel = min;
   fMaxLevel = max;
}//AddUserInfo

//______________________________________________________________________________
void XExpressionTreeInfo::AddUserInfo(Int_t nquant, Double_t *q, Double_t *quant)
{
   // Add user info from tree set
   if(kCS) cout << "------XExpressionTreeInfo::AddUserInfo------" << endl;

   if (nquant > fNQuantiles) {
      if (fLevelQuant) {delete [] fLevelQuant; fLevelQuant = 0;}
      if (fQuantiles)  {delete [] fQuantiles;  fQuantiles  = 0;}

      fQuantiles  = new Double_t[nquant];
      fLevelQuant = new Double_t[nquant];
   }//if

   fNQuantiles = nquant;

   memcpy(fQuantiles,  q,     nquant*sizeof(Double_t));
   memcpy(fLevelQuant, quant, nquant*sizeof(Double_t));
}//AddUserInfo

//______________________________________________________________________________
Double_t XExpressionTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCS) cout << "------XExpressionTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNUnits") == 0) {
      return fNUnits;
   } else if (strcmp(name, "fMinLevel") == 0) {
      return fMinLevel;
   } else if (strcmp(name, "fMaxLevel") == 0) {
      return fMaxLevel;
   } else if (strcmp(name, "fNQuantiles") == 0) {
      return fNQuantiles;
   }//if
   return 0;
}//GetValue

//______________________________________________________________________________
Double_t *XExpressionTreeInfo::GetQuantiles()
{
   // Return quantiles array
   if(kCS) cout << "------XExpressionTreeInfo::GetQuantiles------" << endl;

   if (fQuantiles == 0) return 0;

   return fQuantiles;
}//GetQuantiles

//______________________________________________________________________________
Double_t *XExpressionTreeInfo::GetLevelQuantiles()
{
   // Return quantiles array for levels
   if(kCS) cout << "------XExpressionTreeInfo::GetLevelQuantiles------" << endl;

   if (fLevelQuant == 0) return 0;

   return fLevelQuant;
}//GetLevelQuantiles


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSelectionTreeInfo                                                   //
//                                                                      //
// Class containing info about selection mask tree stored in fUserInfo  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSelectionTreeInfo::XSelectionTreeInfo() 
              :XTreeInfo()
{
   // Default MaskTreeInfo constructor
   if(kCS) cout << "---XSelectionTreeInfo::XSelectionTreeInfo(default)------" << endl;

   fNUnits = 0;
   fNFlags = 0;
}//Constructor

//______________________________________________________________________________
XSelectionTreeInfo::XSelectionTreeInfo(const char *name, const char *title) 
              :XTreeInfo(name, title)
{
   // Normal MaskTreeInfo constructor
   if(kCS) cout << "---XSelectionTreeInfo::XSelectionTreeInfo------" << endl;

   fNUnits = 0;
   fNFlags = 0;
}//Constructor

//______________________________________________________________________________
XSelectionTreeInfo::~XSelectionTreeInfo()
{
   // MaskTreeInfo destructor
   if(kCS) cout << "---XSelectionTreeInfo::~XSelectionTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XSelectionTreeInfo::AddUserInfo(Int_t nunits, Int_t nflags)
{
   // Add user info from tree set
   if(kCS) cout << "------XSelectionTreeInfo::AddUserInfo------" << endl;

   fNUnits = nunits;
   fNFlags = nflags;
}//AddUserInfo

//______________________________________________________________________________
Double_t XSelectionTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCS) cout << "------XSelectionTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNUnits") == 0) {
      return fNUnits;
   } else if (strcmp(name, "fNFlags") == 0) {
      return fNFlags;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProcesSetting                                                       //
//                                                                      //
// Class for initialization of input settings                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProcesSetting::XProcesSetting()
               :XSetting()
{
   // Default ProcesSetting constructor
   if(kCS) cout << "---XProcesSetting::XProcesSetting(default)------" << endl;

   fSchemeFile = 0;
   fSchemeName = "";
}//Constructor

//______________________________________________________________________________
XProcesSetting::XProcesSetting(const char *arraytype, const char *infile)
               :XSetting(arraytype, infile)
{
   // Normal ProcesSetting constructor
   if(kCS) cout << "---XProcesSetting::XProcesSetting------" << endl;

   fSchemeFile = 0;
   fSchemeName = "";
}//Constructor

//______________________________________________________________________________
XProcesSetting::~XProcesSetting()
{
   // ProcesSetting destructor
   if(kCS) cout << "---XProcesSetting::~XProcesSetting------" << endl;

   fSchemeFile = 0;
}//Destructor

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProcesSet                                                           //
//                                                                      //
// Base class containing info about tree sets used for processing of    //
// microarray data stored as trees in root files                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProcesSet::XProcesSet()
           :XTreeSet()
{
   // Default ProcesSet constructor
   if(kCS) cout << "---XProcesSet::XProcesSet(default)------" << endl;

   fSchemeFile = 0;
   fDataFile   = 0;
   fSchemes    = 0;
   fData       = 0;
   fSchemeName = "";
   fBaselines  = new TList();
   fReferences = new TList();
//??   fBaselines  = 0;
//??   fReferences = 0; //not allowed since need to be created!!
   fBaseOpt    = "";
   fRefOpt     = "";
   fBaseTrim   = 0.0;
   fRefTrim    = 0.0;
}//Constructor

//______________________________________________________________________________
XProcesSet::XProcesSet(const char *name, const char *type)
           :XTreeSet(name, type)
{
   // Normal ProcesSet constructor
   if(kCS) cout << "---XProcesSet::XProcesSet------" << endl;

   fSchemeFile = 0;
   fDataFile   = 0;
   fSchemes    = 0;
   fData       = 0;
   fSchemeName = "SchemeName";
   fBaselines  = new TList();
   fReferences = new TList();
   fBaseOpt    = "";
   fRefOpt     = "";
   fBaseTrim   = 0.0;
   fRefTrim    = 0.0;
}//Constructor

//______________________________________________________________________________
XProcesSet::~XProcesSet()
{
   // XrocesSet destructor
   if(kCS) cout << "---XProcesSet::~XProcesSet------" << endl;

   fReferences->Clear("nodelete");  //do not delete trees
   SafeDelete(fReferences);
   fBaselines->Clear("nodelete");  //do not delete trees
   SafeDelete(fBaselines);

   // delete content first, i.e. each chip and hyb
   if (fSchemes)  {fSchemes->Delete(); delete fSchemes; fSchemes = 0;}
   if (fData)     {fData->Delete();    delete fData;    fData    = 0;}

   fSchemeFile = 0;
   fDataFile   = 0;
}//Destructor

//______________________________________________________________________________
Int_t XProcesSet::InitFiles(TFile *schemefile, TFile *datafile)
{
   // Initialize root files
   if(kCS) cout << "------XProcesSet::InitFiles------" << endl;

   fSchemeFile = schemefile;
   fDataFile   = datafile;

   if (fSchemeFile && fDataFile) return errNoErr;

   cerr << "Error: File(s) were not initialized" << endl;
   return errAbort;
}//InitFiles

//______________________________________________________________________________
Int_t XProcesSet::SetBaseLine(TTree *tree, Option_t *option, Double_t trim)
{
   // Add selected tree to baselines list fBaselines
   if(kCS) cout << "------XProcesSet::SetBaseLine------" << endl;

   if (tree != 0) {
   // Add tree to list
      fBaselines->Add(tree);

      // add option only once
      if ((strcmp(fBaseOpt, "") == 0) && (strcmp(option, "") != 0)) {
         fBaseOpt  = option;
         fBaseTrim = trim;
      }//if
   } else {
   // Change option and trim
      fBaseOpt  = option;
      fBaseTrim = trim;
   }//if

   fBaseOpt.ToLower();

   return errNoErr;
}//SetBaseLine

//______________________________________________________________________________
Int_t XProcesSet::SetReference(TTree *tree, Option_t *option, Double_t trim)
{
   // Add selected tree to references list fReferences
//TO DO:
   // Note: Tree(s) with name need not be part of list fSelections, e.g. 
   //       general reference tree(s) stored in different root file(s) possible
   if(kCS) cout << "------XProcesSet::SetReference------" << endl;

   if (tree != 0) {
   // Add tree to list
      fReferences->Add(tree);

      // add option only once
      if ((strcmp(fRefOpt, "") == 0) && (strcmp(option, "") != 0)) {
         fRefOpt  = option;
         fRefTrim = trim;
      }//if
   } else {
   // Change option and trim
      fRefOpt  = option;
      fRefTrim = trim;
   }//if

   fRefOpt.ToLower();

   return errNoErr;
}//SetReference

//______________________________________________________________________________
void XProcesSet::AddExprTreeInfo(TTree *tree, const char *name, Option_t *option,
                 Int_t nunits, Double_t min, Double_t max,
                 Int_t nquant, Double_t *q, Double_t *quant)
{
   // Add expression tree info to list fUserInfo of tree
   if(kCS) cout << "------XProcesSet::AddExprTreeInfo------" << endl;

// store name of tree set as title
   XExpressionTreeInfo *info = new XExpressionTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info
   if (nunits > 0) info->AddUserInfo(nunits, min, max);
   if (nquant > 0) info->AddUserInfo(nquant, q, quant);

   tree->GetUserInfo()->Add(info);
}//AddExprTreeInfo

//______________________________________________________________________________
Int_t XProcesSet::InitGroups(Int_t &n, Int_t *gid, TTree **tree, const char **extens)
{
   // Initialize group IDs for n trees
   // Only group IDs for trees with correct name and extension (from list extens)
   // will be added to array gid.
   // Note: n will return the array size of gid for correct trees only
   if(kCS) cout << "------XProcesSet::InitGroups------" << endl;

// Get group ID for trees with correct name and extension only
   Int_t idx = 0;
   for (Int_t k=0; k<n; k++) {
      TString    name  = tree[k]->GetName();
      TString    exten = Path2Name(tree[k]->GetName(),".","");
      XIdxString *str  = (XIdxString*)(fSelections->FindObject(name));
      if ((strcmp(name.Data(),str->GetName()) == 0) &&
           HasExtension(exten.Data(), extens)) {
         tree[idx] = tree[k];
         gid[idx]  = str->GetIndex();
         idx++;
      }//if
   }//for_k

// Return size of array for correct trees only
   n = idx;

// Check if at least 2 trees are selected
   if (n < 2) {
      cout << "Error: Less than two trees selected" << endl;
      return errGetTree;
   }//if

   return errNoErr;
}//InitGroups

//______________________________________________________________________________
Int_t XProcesSet::ExportExprTreeInfo(Int_t n, TString *names, const char *varlist,
                  ofstream &output, const char *sep)
{
   // Export data stored in userinfo of expression tree to file output
   if(kCS) cout << "------XProcesSet::ExportExprTreeInfo------" << endl;

// Decompose varlist
   Bool_t hasTreeName = kFALSE;
   Bool_t hasSetName  = kFALSE;
   Bool_t hasOption   = kFALSE;
   Bool_t hasNUnits   = kFALSE;
   Bool_t hasMinLevel = kFALSE;
   Bool_t hasMaxLevel = kFALSE;
   Bool_t hasNQuant   = kFALSE;
   Bool_t hasQuant    = kFALSE;
   Bool_t hasLevel    = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasTreeName = kTRUE;
      hasSetName  = kTRUE;
      hasOption   = kTRUE;
      hasNUnits   = kTRUE;
      hasMinLevel = kTRUE;
      hasMaxLevel = kTRUE;
      hasNQuant   = kTRUE;
      hasQuant    = kTRUE;
      hasLevel    = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fName")       == 0) {hasTreeName = kTRUE;}
         if (strcmp(name,"fSetName")    == 0) {hasSetName  = kTRUE;}
         if (strcmp(name,"fOption")     == 0) {hasOption   = kTRUE;}
         if (strcmp(name,"fNUnits")     == 0) {hasNUnits   = kTRUE;}
         if (strcmp(name,"fMinLevel")   == 0) {hasMinLevel = kTRUE;}
         if (strcmp(name,"fMaxLevel")   == 0) {hasMaxLevel = kTRUE;}
         if (strcmp(name,"fNQuantiles") == 0) {hasNQuant   = kTRUE;}
         if (strcmp(name,"fQuantiles")  == 0) {hasQuant    = kTRUE;}
         if (strcmp(name,"fLevelQuant") == 0) {hasLevel    = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree               **tree = new TTree*[n];
   XExpressionTreeInfo **info = new XExpressionTreeInfo*[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         info[k] = (XExpressionTreeInfo*)tree[k]->GetUserInfo()->At(0);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         info[k] = (XExpressionTreeInfo*)tree[k]->GetUserInfo()->At(0);
      }//for_k
   }//if

// Output header
   output << "Parameter";
   for (Int_t k=0; k<n; k++) output << sep << names[k].Data();
   output << endl;

// Output parameters
   if (hasTreeName) {
      output << "TreeName";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetName();
      output << endl;
   }//if

   if (hasSetName) {
      output << "SetName";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetTreeSetName();
      output << endl;
   }//if

   if (hasOption) {
      output << "Option";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetOption();
      output << endl;
   }//if

   if (hasNUnits) {
      output << "NumUnits";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("fNUnits"));
      output << endl;
   }//if

   if (hasMinLevel) {
      output << "MinLevel";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetValue("fMinLevel");
      output << endl;
   }//if

   if (hasMaxLevel) {
      output << "MaxLevel";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetValue("fMaxLevel");
      output << endl;
   }//if

   if (hasNQuant) {
      output << "NumQuantiles";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("fNQuantiles"));
      output << endl;
   }//if

   if (hasQuant) {
      Double_t **quant = new Double_t*[n];
      for (Int_t k=0; k<n; k++) {
         quant[k] = info[k]->GetQuantiles();
      }//for_k

      Int_t nq  = (Int_t)(info[0]->GetValue("fNQuantiles"));
      for (Int_t i=0; i<nq; i++) {
         TString str; str = "Quantile"; str += i;

         output << str.Data();
         for (Int_t k=0; k<n; k++) output << sep << quant[k][i];
         output << endl;
      }//for_i

      delete [] quant;
   }//if

   if (hasLevel) {
      Double_t **quant = new Double_t*[n];
      Double_t **level = new Double_t*[n];
      for (Int_t k=0; k<n; k++) {
         quant[k] = info[k]->GetQuantiles();
         level[k] = info[k]->GetLevelQuantiles();
      }//for_k

      Int_t nq  = (Int_t)(info[0]->GetValue("fNQuantiles"));
      for (Int_t i=0; i<nq; i++) {
         TString str; str.Form("Level_Q%4.2f", quant[0][i]);

         output << str.Data();
         for (Int_t k=0; k<n; k++) output << sep << level[k][i];
         output << endl;
      }//for_i

      delete [] level;
      delete [] quant;
   }//if

// Cleanup
   for (Int_t k=0; k<n; k++) {
//no!      SafeDelete(info[k]);
      SafeDelete(tree[k]);
   }//for_k

   delete [] info;
   delete [] tree;

   return errNoErr;
}//ExportExprTreeInfo

//______________________________________________________________________________
Int_t XProcesSet::ExportMaskTreeInfo(Int_t n, TString *names, const char *varlist,
                  ofstream &output, const char *sep)
{
   // Export data stored in userinfo of mask tree to file output
   if(kCS) cout << "------XProcesSet::ExportMaskTreeInfo------" << endl;

   Int_t err = errNoErr;

   return errNoErr;
}//ExportMaskTreeInfo

//______________________________________________________________________________
const char *XProcesSet::GetTranscriptID(XTransAnnotation *annot)
{
   // Get transcript IDs from annotation tree

   return annot->GetTranscriptID();
}//GetTranscriptID

//______________________________________________________________________________
const char *XProcesSet::GetTranscriptID(XUnit *unit, XTransAnnotation *annot, Int_t type)
{
   // Get transcript IDs from unit tree

   if (strcmp(this->GetTitle(), "GeneChip") == 0) {
      return annot->GetTranscriptID();
   }//if

   Int_t id = (type == eTRANSCRIPT) 
            ? ((XExonUnit*)unit)->GetSubUnitID()
            : ((XProbesetAnnotation*)annot)->GetTranscriptIX();

   return Form("%d", id);
}//GetTranscriptID

//______________________________________________________________________________
Int_t XProcesSet::HandleOption(TTree *tree, Option_t *opt)
{
   // Handle option
   if(kCS) cout << "------XProcesSet::HandleOption------" << endl;

   Int_t err = errNoErr;

   if (strcmp(opt, "") == 0) {
      return err;
   } else if (strcmp(opt, "reference") == 0) {
      err = this->SetReference(tree);
   } else if (strcmp(opt, "baseline")  == 0) {
      err = this->SetBaseLine(tree);
   } else if (strcmp(opt, "baseref")   == 0) {
      err = this->SetBaseLine(tree);
      if (err == errNoErr) err = this->SetReference(tree);
   }//if

   return err;
}//HandleOption

//______________________________________________________________________________
Int_t XProcesSet::CopyUnitBranch(TTree *fromtree, TTree *totree, Int_t writeopt)
{
   // Copy UnitBranch from tree "fromtree" to new UnitBranch from "totree"
   // Write totree to current file:
   // writeopt = -1: do not write tree to file
   // writeopt = 0:  write tree to current file
   // writeopt = TObject::kOverwrite or TObject::kWriteDelete: overwrite tree
   if(kCS) cout << "------XProcesSet::CopyUnitBranch------" << endl;

   if ((fromtree == 0) || (totree == 0)) return errGetTree;

//TO DO: check if totree has already unit branch!!!

   // check if trees have equal entries
   Int_t nentries = (Int_t)(fromtree->GetEntries());
   if ((Int_t)(totree->GetEntries()) != nentries) {
      return fManager->HandleError(errEQTreeEntries, fromtree->GetName(), totree->GetName());
   }//if

   // get fUnitID from fromtree
   TLeaf *leaf = fromtree->FindLeaf("fUnitID");
   if (!leaf) {
      cout << "Warning: Tree <" << fromtree->GetName() << "> has no UnitBranch."
           << endl;
      return errNoErr; //??
   }//if
   TBranch *brch = leaf->GetBranch();

   // create unit branch
   XUnitID *unit = new XUnitID();
   TBranch *idbr = totree->Branch("UnitBranch", "XUnitID", &unit, 64000, 99);

   // fill unit branch with fUnitID
   Int_t id = 0;
   for (Int_t i=0; i<nentries; i++) {
      brch->GetEntry(i);
      id = (Int_t)leaf->GetValue();
      unit->SetUnitID(id);
      idbr->Fill();
   }//for_i

   // write totree with unit branch to file
   if (writeopt != -1) totree->Write("", writeopt);

   SafeDelete(unit);
   totree->ResetBranchAddress(totree->GetBranch("UnitBranch"));

   return errNoErr;
}//CopyUnitBranch

//______________________________________________________________________________
Int_t XProcesSet::ExpressionQuantiles(TTree *exprtree, XExpression *expr,
                  Int_t nquant, Double_t *q, Double_t *quantL)
{
   // Get nquant quantiles for expression levels from expression tree
   if(kCS) cout << "------XProcesSet::ExpressionQuantiles------" << endl;

   Int_t err      = errNoErr;
   Int_t nentries = (Int_t)(exprtree->GetEntries());

   exprtree->SetBranchAddress("ExprBranch", &expr);

// Init arrays
   Double_t *level = 0; 
   Int_t    *index = 0;
   if (!(level = new (nothrow) Double_t[nentries])) {err = errInitMemory; goto cleanup;}
   if (!(index = new (nothrow) Int_t[nentries]))    {err = errInitMemory; goto cleanup;}

// Fill arrays
   for (Int_t i=0; i<nentries; i++) {
      exprtree->GetEntry(i);
      level[i] = expr->GetLevel();
   }//for_i

// Fill quantiles
   quantL = TStat::Quantiles(nentries, level, index, nquant, q, quantL);

// Cleanup
cleanup:
   exprtree->DropBaskets();
   exprtree->ResetBranchAddress(exprtree->GetBranch("ExprBranch"));

   if (index) {delete [] index; index = 0;}
   if (level) {delete [] level; level = 0;}

   return err;
}//ExpressionQuantiles

//______________________________________________________________________________
TTree *XProcesSet::GetUnitTree(XGeneChip *chip, Int_t type)
{
   // Get unit tree of type
   if(kCS) cout << "------XProcesSet::GetUnitTree------" << endl;

   if (!fSchemeFile->cd(fSchemeName)) return 0;

   TTree *unittree = 0; 
   if (type == ePROBESET) { 
      unittree = (TTree*)gDirectory->Get(((XExonChip*)chip)->GetProbesetUnitTree()); 
   } else if (type == eEXONTYPE) {
      unittree = (TTree*)gDirectory->Get(((XExonChip*)chip)->GetExonUnitTree()); 
   } else if (type == eTRANSCRIPT) {
      unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   } else {
      cerr << "Error: Unknown unit tree type" << endl;
   }//if
   if (unittree == 0) return 0;

   return unittree;
}//GetUnitTree

//______________________________________________________________________________
TTree *XProcesSet::GetAnnotationTree(XGeneChip *chip, Int_t type)
{
   // Get annotation tree of type
   if(kCS) cout << "------XProcesSet::GetAnnotationTree------" << endl;

   if (!fSchemeFile->cd(fSchemeName)) return 0;

   TTree *anntree  = 0; 
   if (type == ePROBESET) { 
      anntree = (TTree*)gDirectory->Get(((XExonChip*)chip)->GetProbesetAnnotTree()); 
   } else if (type == eEXONTYPE) {
      anntree = (TTree*)gDirectory->Get(((XExonChip*)chip)->GetExonAnnotTree()); 
   } else if (type == eTRANSCRIPT) {
      anntree = (TTree*)gDirectory->Get(((XExonChip*)chip)->GetAnnotTree()); 
   } else {
      cerr << "Error: Unknown annotation tree type" << endl;
   }//if
   if (anntree == 0) return 0;

   return anntree;
}//GetAnnotationTree

//______________________________________________________________________________
THashTable *XProcesSet::FillHashTable(THashTable *htable, TTree *anntree,
                        XTransAnnotation *annot)
{
   // Fill hash table with IDs from annotation tree
   if(kCS) cout << "------XProcesSet::FillHashTable------" << endl;

   if (XManager::fgVerbose) {
      cout << "Reading entries from <" << anntree->GetName() << "> ...";
   }//if

   XIdxString *idxstr = 0;
   Int_t numannot = (Int_t)(anntree->GetEntries());
   for (Int_t i=0; i<numannot; i++) {
      anntree->GetEntry(i);

      idxstr = new XIdxString(i, annot->GetTranscriptID());
      htable->Add(idxstr);
   }//for_i

   if (XManager::fgVerbose) {
      cout << "Finished" << endl;
   }//if

   return htable;
}//FillHashTable

//______________________________________________________________________________
THashTable *XProcesSet::FillHashTable(THashTable *htable, TTree *anntree,
                        XTransAnnotation *annot, Int_t type)
{
   // Fill hash table with IDs from annotation tree of type
   if(kCS) cout << "------XProcesSet::FillHashTable(type)------" << endl;

   if (XManager::fgVerbose) {
      cout << "Reading entries from <" << anntree->GetName() << "> ...";
   }//if

   TString str;
   XIdxString *idxstr = 0;

   Int_t numannot = (Int_t)(anntree->GetEntries());
   for (Int_t i=0; i<numannot; i++) {
      anntree->GetEntry(i);

      if (type == ePROBESET) { 
         str.Form("%d", ((XProbesetAnnotation*)annot)->GetProbesetID());
         idxstr = new XIdxString(i, str.Data());
      } else if (type == eEXONTYPE) {
         str.Form("%d", ((XExonAnnotation*)annot)->GetExonID());
         idxstr = new XIdxString(i, str.Data());
      } else if (type == eTRANSCRIPT) {
         str = ((XGenomeAnnotation*)annot)->GetTranscriptID();
         idxstr = new XIdxString(i, str.Data());
      } else {
         cerr << "Error: Unknown annotation type" << endl;
      }//if

      htable->Add(idxstr);
   }//for_i

   if (XManager::fgVerbose) {
      cout << "Finished" << endl;
   }//if

   return htable;
}//FillHashTable

//______________________________________________________________________________
XIdxString *XProcesSet::FindUnitID(THashTable *htable, XGCUnit *unit)
{
   // Find unit ID in hash table with IDs from annotation tree of type

   return (XIdxString*)(htable->FindObject(unit->GetUnitName()));
}//FindUnitID

//______________________________________________________________________________
XIdxString *XProcesSet::FindUnitID(THashTable *htable, XExonUnit *unit)
{
   // Find unit ID in hash table with IDs from annotation tree of type

   TString str;
   str.Form("%d", unit->GetSubUnitID());

   return (XIdxString*)(htable->FindObject(str.Data()));
}//FindUnitID
