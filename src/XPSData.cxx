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
* Jan 2003 - Add classes XDataSetting, XGeneChipMetrics, XGeneChipPivot,
*            and plot classes XAnalysisPlot, XGeneChipPlot, XGenePixPlot.
*          - Support to read  *.CEL files in XML format.
* Apr 2005 - Add classes XDataTreeInfo and XMaskTreeInfo.
*          - Support to read binary *.CEL files.
* May 2005 - Replace XContent with XFolder.
* Jun 2005 - Add methods to support database.
* Jul 2005 - Replace TTree::AddFriend() with TList fTrees.
*          - Store trees for every XTreeSet in separate TDirectory.
* Jan 2006 - To decrease size of saved trees replace Double_t with Double32_t
*            in classes XCell, XExpression and their derived classes, and Int_t
*            with Short_t for fMask, fCall and fNPixels
* Mar 2006 - Add support for alternative schemes, i.e. CDF files.
* Aug 2006 - Add support for exon arrays, classes XExonChipHyb, XExonChipMetrics
*            XExonChipPivot.
* Aug 2007 - Add support for whole genome arrays, classes XGenomeChipHyb, 
*            XGenomeChipMetrics, XGenomeChipPivot.
* Nov 2007 - XDataManager now inherits from XManager and XProjectHandler
*          - database methods ProjectInfo() etc are moved to XProjectHandler
* Jan 2008 - Support to read Calvin (AGCC) *.CEL files.
* Feb 2010 - Add support for alternative splicing, class XSpliceExpression.
* Dec 2010 - Add support for quality control, classes XGCQuality, XQCExpression.
*
******************************************************************************/

using namespace std;

#include <cstdlib>
#include <cfloat>
#include <new>  //needed for new (nothrow)

#include "TFriendElement.h"
#include "THashTable.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

//TEST BENCHMARK:
#include <TBenchmark.h>

#include "XPSData.h"
#include "XPSDataTypes.h"
#include "TStat.h"

// Tree extensions and types:
// data trees
const char *kExtenData[] = { "cel", ""};
const char *kTypeData[]  = { "rawdata",
                             ""};

// mask trees
const char *kExtenMask[] = { "msk", ""};
const char *kTypeMask[]  = { "mask",
                             ""};

// Local tree extensions for import of GeneChip Metrics or PivotTable
// expression trees (see kExtenExpr in XPSProcessing)
const char *kPivotExpr[] = { "amn", "gmn", "wmn", "gcm", "wdf", "adf", "tbw", "mdp", ""};

// detection call trees (see kExtenCall in XPSProcessing)
const char *kPivotCall[] = { "mdf", "dc5", ""};


// Allowed synonyms for variables in varlist
const Int_t kNSynExpr = 3;
const char *kSynExpr[3] = { "EXPRESSION",
                            "EXPR", 
                            "SIGNAL" };

const Int_t kNSynCall = 4;
const char *kSynCall[4] = { "CALL",
                            "PRESENT CALL", 
                            "PRESENT_CALL", 
                            "DETECTION" };

const Int_t kNSynPVal = 8;
const char *kSynPVal[8] = { "DETECTION P-VALUE",
                            "DETECTION P_VALUE", 
                            "DETECTION PVALUE", 
                            "DETECTION_PVALUE", 
                            "P-VALUE", 
                            "P_VALUE", 
                            "PVALUE", 
                            "PVAL" };

// struct for binary XDA files
typedef struct CELEntryType
{
	float Intensity;
	float Stdv;
	short Pixels;
} CELEntry;

ClassImp(XDataManager);
ClassImp(XDataTreeInfo);
ClassImp(XMaskTreeInfo);
ClassImp(XDataSetting);
ClassImp(XHybridization);
ClassImp(XGeneChipHyb);
ClassImp(XGeneChipMetrics);
ClassImp(XGeneChipPivot);
ClassImp(XSNPChipHyb);
ClassImp(XSNPChipMetrics);
ClassImp(XSNPChipPivot);
ClassImp(XGenomeChipHyb);
ClassImp(XGenomeChipMetrics);
ClassImp(XGenomeChipPivot);
ClassImp(XExonChipHyb);
ClassImp(XExonChipMetrics);
ClassImp(XExonChipPivot);
ClassImp(XGenePixHyb);
ClassImp(XCellOM);
ClassImp(XMask);
ClassImp(XCellMask);
ClassImp(XUnitMask);
ClassImp(XCell);
ClassImp(XGCCell);
ClassImp(XBgCell);
ClassImp(XRankCell);
ClassImp(XResidual);
ClassImp(XBorder);
ClassImp(XExpression);
ClassImp(XGCExpression);
ClassImp(XQCExpression);
ClassImp(XSpliceExpression);
ClassImp(XCall);
ClassImp(XPCall);
ClassImp(XGPCell);
ClassImp(XFeature635);
ClassImp(XFeature532);
ClassImp(XBg532);
ClassImp(XGPRatio);
ClassImp(XAnalysisPlot);
ClassImp(XGeneChipPlot);
ClassImp(XGenePixPlot);

//?? kBrchBuf=???? better: calculate branchbuffer??
const Int_t  kBrchBuf  = 128000;  //bufsize for tree branch
const Int_t  kLargeBuf = 16635;   //to read lines from PivotData

const Int_t kMagicXDANumber    = 64;
const Int_t kMagicCalvinNumber = 59;

//debug: print function names
const Bool_t kCS  = 0; 
const Bool_t kCSa = 0; //debug: print function names in loops


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataManager                                                         //
//                                                                      //
// Class for managing import of microarray data                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDataManager::XDataManager()
             :XManager(), XProjectHandler()
{
   // Default DataManager constructor
   if(kCS) cout << "---XDataManager::XDataManager(default)------" << endl;

   fSchemeFile    = 0;
   fSchemes       = 0;
   fIsSchemeOwner = kFALSE;
}//Constructor

//______________________________________________________________________________
XDataManager::XDataManager(const char *name, const char *title, Int_t verbose)
             :XManager(name, title, verbose), XProjectHandler(name, title)
{
   // Normal DataManager constructor
   if(kCS) cout << "---XDataManager::XDataManager------" << endl;

   fSchemeFile    = 0;
   fSchemes       = 0;
   fIsSchemeOwner = kFALSE;

//TEST BENCHMARK:
//   gBenchmark = new TBenchmark();
}//Constructor

//______________________________________________________________________________
XDataManager::~XDataManager()
{
   // DataManager destructor
   if(kCS) cout << "---XDataManager::~XDataManager------" << endl;

//   this->Close(); //not here but in base destructor

//TEST BENCHMARK:
//   delete gBenchmark;
}//Destructor

//______________________________________________________________________________
Int_t XDataManager::InitInput(const char *schemename, const char *datatype,
                    const char *varlist, const char *inputtype)
{
   // Initialize settings for input data where:
   // - schemename is the exact name of array type
   // - datatype should describe type of chip data, e.g. cel for rawdata,
   //   tbw, mdp for Metrics or PivotData (but also mas5 or rma)
   // - varlist is list of variables to be stored in tree(s) and should contain
   //   type of variable, e.g. "var1/I:var2/C:var3/D"
   // - inputtype is rawdata (default), Metrics, PivotData
   if(kCS) cout << "------XDataManager::InitInput------" << endl;

   if (fAbort) return errAbort;
   if (!fSetting) {this->HandleError(errInitSetting); return errAbort;}

// Initialize input settings
   Int_t err = errNoErr;
   err = ((XDataSetting*)fSetting)->InitInput(schemename, datatype, varlist, inputtype);
   if (err) fAbort = kTRUE;
   
   return err;
}//InitInput

//______________________________________________________________________________
Int_t XDataManager::InitSchemes(TFile *schemefile, Bool_t isOwner,
                    const char *schemename)
{
   // Initialize data manager with external root schemefile
   // Set isOwner = kTRUE if schemefile should be deleted by datamanager
   // Optional select scheme "schemename", e.g. in GUI application to ensure
   // that user has selected correct scheme.
   // Note: Do not use schemename for alternative scheme!
   if(kCS) cout << "------XDataManager::InitSchemes------" << endl;

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
      return this->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   if (fSetting) {
      ((XDataSetting*)fSetting)->SetSchemeFile(fSchemeFile);
   }//if
   // overwrite InitInput(schemename) only if optional schemename is given
   if (fSetting && (strcmp(schemename, "") != 0)) {
      ((XDataSetting*)fSetting)->SetSchemeName(schemename);
   }//if

   savedir->cd();
   return errNoErr;
}//InitSchemes

//______________________________________________________________________________
Int_t XDataManager::OpenSchemes(const char *fullname, const char *schemename)
{
   // Open root scheme file using full path 
   // Optional select scheme "schemename", e.g. in GUI application to ensure
   // that user has selected correct scheme.
   // Note: Do not use schemename for alternative scheme!
   if(kCS) cout << "------XDataManager::OpenSchemes------" << endl;

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
      return this->HandleError(errMissingContent, "Scheme", kContent);
   }//if

   if (fSetting) {
      ((XDataSetting*)fSetting)->SetSchemeFile(fSchemeFile);
   }//if
   // overwrite InitInput(schemename) only if optional schemename is given
   if (fSetting && (strcmp(schemename, "") != 0)) {
      ((XDataSetting*)fSetting)->SetSchemeName(schemename);
   }//if

   savedir->cd();
   return errNoErr;
}//OpenSchemes

//______________________________________________________________________________
void XDataManager::CloseSchemes()
{
   // Close scheme file
   if(kCS) cout << "------XDataManager::CloseSchemes------" << endl;

   SafeDelete(fSchemes); //?? ev before:

   if (fSchemeFile) {
      if (fIsSchemeOwner) DeleteFile(fSchemeFile);
      fSchemeFile = 0;  //if not isOwner!
   }//if
}//CloseSchemes

//______________________________________________________________________________
Int_t XDataManager::New(const char *name, const char *dir, const char *type,
                    const char *data, Option_t *option)
{
   // Create new root data file with name in directory dir for type
   if(kCS) cout << "------XDataManager::New------" << endl;

// Get extension for file name
   TString datatype = "def";
   if (fSetting) datatype = ((XDataSetting*)fSetting)->GetDataType();
   datatype.ToLower();

// Check array type
   TString filename = TString(name);
   if ((strcmp(type,"GeneChip")   == 0) ||
       (strcmp(type, "SNPChip")   == 0) ||
       (strcmp(type,"GenomeChip") == 0) ||
       (strcmp(type,"ExonChip")   == 0) ||
       (strcmp(type, "GenePix")   == 0) ||
       (strcmp(type,     "GEM")   == 0) ||
       (strcmp(type,  "Custom")   == 0)) {
      filename += "_" + datatype;
   } else {
      cerr << "Error: Array type <" << type << "> not known." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   return XManager::New(filename, dir, type, data, option);
}//New

//______________________________________________________________________________
void XDataManager::Close(Option_t *option)
{
   // Close root files
   if(kCS) cout << "------XDataManager::Close------" << endl;

   this->CloseSchemes();

   // allow Close/Save in XPSApp w/o calling destructor
   XManager::Close(option);
}//Close

//______________________________________________________________________________
Int_t XDataManager::Import(const char *setname, const char *infile,
                    const char *treename, Option_t *option, const char *sep,
                    char delim, Int_t split)
{
   // Import data belonging to dataset setname from file infile. 
   if(kCS) cout << "------XDataManager::Import------" << endl;

   fInterrupt = kFALSE;
   if (fAbort)    {fInterrupt = kTRUE; return errAbort;}
   if (!fSetting) {return this->HandleError(errInitSetting);}

// Add extension for treename of data type to option
   TString datatype = ((XDataSetting*)fSetting)->GetDataType();
   TString opt      = Path2Name(option,"",".");
   opt.ToUpper();
   datatype.ToLower();
   opt = opt + "." + datatype;

   return XManager::Import(setname, infile, treename, opt.Data(), sep, delim, split);
}//Import

//______________________________________________________________________________
Int_t XDataManager::HandleError(Int_t err, const char *name1, const char *name2)
{
   // Handle error messages
   if(kCS) cout << "------XDataManager::HandleError------" << endl;

   switch (err) {
/*      case errExtension:
         cerr << "Error: Tree(s) with extension <" << name1 << "> not known" << endl;
         fInterrupt = kTRUE;
         break;
*/
      case errChipType:
         cerr << "Error: Selected scheme <" << name1 
           << "> is not identical to imported scheme <" << name2 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errCELVersion:
         cerr << "Error: CEL-file with version/magic number <" << name1 
              << "> is not supported." << endl;
         fAbort = kTRUE;
         break;

      case errNameValue:
      cerr << "Error: NameValueType value not found!" << endl;
         fInterrupt = kTRUE;
         break;

      case errNumCells:
         cerr << "Error: number of lines read <" << name1 
              << "> is not equal to to number of cells <" << name2 << ">."
              << endl;
         fInterrupt = kTRUE;
         break;

      default:
         XSchemeManager manager;
         err = manager.HandleError(err, name1, name2);
//         err = XManager::HandleError(err, name1, name2);
         break;
   }//switch

   return err;
}//HandleError

//______________________________________________________________________________
Int_t XDataManager::DrawUnit(const char *canvasname, const char *treename,
                    const char *schemename, const char *unitname, Int_t unitID,
                    const char *varlist, const char *logbase, const char *type,
                    Option_t *opt, Double_t min, Double_t max, const char *var2sort)
{
   // Draw raw data from tree treename for scheme with schemename and unit with
   // unitname or for unitID (if unitname = "").
   // Variables to be drawn are listed in varlist.
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // type: graph, multigraph, grapherr, hist with the corresponding option opt
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Variable var2sort determines the variable to sort for
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XDataManager::DrawUnit------" << endl;

   if (fAbort) return errAbort;

   Int_t err = errNoErr;
   err = ((XAnalysisPlot*)fPlotter)->DrawUnit(canvasname, treename, schemename,
                                     unitname, unitID, varlist, logbase, type,
                                     opt, min, max, var2sort);
   return err;
}//DrawUnit

//______________________________________________________________________________
Int_t XDataManager::BeginTransaction(const char *name)
{
   // Begin transaction
   if(kCS) cout << "------XDataManager::BeginTransaction------" << endl;

   if (fAbort) return errAbort;
   
   return errNoErr;
}//BeginTransaction

//______________________________________________________________________________
Int_t XDataManager::CommitTransaction()
{
   // Commit transaction
   if(kCS) cout << "------XDataManager::CommitTransaction------" << endl;

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
                       << oldname << "> to dataset name <" << newname <<">." << endl;
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
Int_t XDataManager::ImportDefaults(const char *infile)
{
   // Import default values from infile
   // Note that InitDefaults() or ImportDefaults() is called from Initialize()
   if(kCS) cout << "------XDataManager::ImportDefaults------" << endl;

// NOT YET USED - to prevent compiler warnings:
   infile = 0;

   return errNoErr;
}//ImportDefaults

//______________________________________________________________________________
Int_t XDataManager::InitDefaults()
{
   // Initialize default values
   // Note that InitDefaults() or ImportDefaults() is called from Initialize()
   if(kCS) cout << "------XDataManager::InitDefaults------" << endl;

   Int_t err = errNoErr;

// Initialize input data
   err = this->InitInput("", "def", "*", "rawdata");

   return err;
}//InitDefaults

//______________________________________________________________________________
XSetting *XDataManager::NewSetting(const char *type, const char *infile)
{
   // Create new data setting for type
   if(kCS) cout << "------XDataManager::NewSetting------" << endl;

   XDataSetting *setting = new XDataSetting(type, infile);

   if (setting && fSchemeFile) {
      setting->SetSchemeFile(fSchemeFile);
   }//if

   return setting;
}//NewSetting

//______________________________________________________________________________
XTreeSet *XDataManager::NewTreeSet(const char *type)
{
   // Create new Hybridization of selected type 
   if(kCS) cout << "------XDataManager::NewTreeSet------" << endl;

   TString datatype  = ((XDataSetting*)fSetting)->GetDataType();
   TString inputtype = ((XDataSetting*)fSetting)->GetInputType();
   inputtype.ToLower();

   XTreeSet *set = 0;
   if (strcmp(type,"GeneChip") == 0) {
      if (strcmp(inputtype.Data(),"rawdata") == 0) {
         set = new XGeneChipHyb(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"metrics") == 0) {
         set = new XGeneChipMetrics(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"pivotdata") == 0) {
         set = new XGeneChipPivot(datatype.Data(), type);
      } else {
         cerr << "Error: Input type <" << inputtype.Data() << "> not known" << endl;
      }//if
   } else if (strcmp(type,"SNPChip") == 0) {
      if (strcmp(inputtype.Data(),"rawdata") == 0) {
         set = new XSNPChipHyb(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"metrics") == 0) {
         set = new XSNPChipMetrics(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"pivotdata") == 0) {
         set = new XSNPChipPivot(datatype.Data(), type);
      } else {
         cerr << "Error: Input type <" << inputtype.Data() << "> not known" << endl;
      }//if
   } else if (strcmp(type,"GenomeChip") == 0) {
      if (strcmp(inputtype.Data(),"rawdata") == 0) {
         set = new XGenomeChipHyb(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"metrics") == 0) {
         set = new XGenomeChipMetrics(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"pivotdata") == 0) {
         set = new XGenomeChipPivot(datatype.Data(), type);
      } else {
         cerr << "Error: Input type <" << inputtype.Data() << "> not known" << endl;
      }//if
   } else if (strcmp(type,"ExonChip") == 0) {
      if (strcmp(inputtype.Data(),"rawdata") == 0) {
         set = new XExonChipHyb(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"metrics") == 0) {
         set = new XExonChipMetrics(datatype.Data(), type);
      } else if (strcmp(inputtype.Data(),"pivotdata") == 0) {
         set = new XExonChipPivot(datatype.Data(), type);
      } else {
         cerr << "Error: Input type <" << inputtype.Data() << "> not known" << endl;
      }//if
   } else if (strcmp(type,"GenePix") == 0) {
      set = new XGenePixHyb(datatype.Data(), type);
   } else {
      cerr << "Error: Chip type <" << type << "> not known" << endl;
   }//if

//?? not necessary??
   if (((XHybridization*)set)->InitFiles(fSchemeFile) != errNoErr) {
      SafeDelete(set);
      return 0;
   }//if

   return set;
}//NewTreeSet

//______________________________________________________________________________
XPlot *XDataManager::NewPlotter(const char *name, const char *title)
{
   // Create new plotter
   if(kCS) cout << "------XDataManager::NewPlotter------" << endl;

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
Int_t XDataManager::DeleteTreeSetInfo(const char *name)
{
   // Delete additional information of treeset "name"
   if(kCS) cout << "------XDataManager::DeleteTreeSetInfo------" << endl;

   XDatasetInfo *datasetinfo = 0;
   datasetinfo = (XDatasetInfo*)(fContent->FindObject("Dataset", name, "XDatasetInfo"));
   if (datasetinfo) fContent->Remove(datasetinfo);

   return errNoErr;
}//DeleteTreeSetInfo


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTreeInfo                                                        //
//                                                                      //
// Class containing info about data tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDataTreeInfo::XDataTreeInfo() 
              :XTreeInfo()
{
   // Default DataTreeInfo constructor
   if(kCS) cout << "---XDataTreeInfo::XDataTreeInfo(default)------" << endl;

   fNRows      = 0; 
   fNCols      = 0;
   fMinInten   = 0;
   fMaxInten   = 0;
   fNMinInten  = 0;
   fNMaxInten  = 0;
   fMaxNPixels = 0;
   fNQuantiles = 0;
   fQuantiles  = 0;
   fIntenQuant = 0;
}//Constructor

//______________________________________________________________________________
XDataTreeInfo::XDataTreeInfo(const char *name, const char *title) 
              :XTreeInfo(name, title)
{
   // Normal DataTreeInfo constructor
   if(kCS) cout << "---XDataTreeInfo::XDataTreeInfo------" << endl;

   fNRows      = 0; 
   fNCols      = 0;
   fMinInten   = 0;
   fMaxInten   = 0;
   fNMinInten  = 0;
   fNMaxInten  = 0;
   fMaxNPixels = 0;
   fNQuantiles = 7;
   fQuantiles  = new Double_t[fNQuantiles];
   fIntenQuant = new Double_t[fNQuantiles];
}//Constructor

//______________________________________________________________________________
XDataTreeInfo::~XDataTreeInfo()
{
   // DataTreeInfo destructor
   if(kCS) cout << "---XDataTreeInfo::~XDataTreeInfo------" << endl;

   if (fIntenQuant) {delete [] fIntenQuant; fIntenQuant = 0;}
   if (fQuantiles)  {delete [] fQuantiles;  fQuantiles  = 0;}
}//Destructor

//______________________________________________________________________________
void XDataTreeInfo::AddUserInfo(XTreeSet *set)
{
   // Add user info from tree set
   if(kCS) cout << "------XDataTreeInfo::AddUserInfo------" << endl;

   if (strcmp(set->ClassName(), "XGeneChipHyb") == 0) {
      XGeneChipHyb *gcset = (XGeneChipHyb*)set;
      fNRows      = gcset->GetNumRows();
      fNCols      = gcset->GetNumColumns();
      fMinInten   = gcset->GetMinIntensity();
      fMaxInten   = gcset->GetMaxIntensity();
      fNMinInten  = gcset->GetNumMinIntensity();
      fNMaxInten  = gcset->GetNumMaxIntensity();
      fMaxNPixels = gcset->GetMaxNumPixels(); 
   } else if (strcmp(set->ClassName(), "XSNPChipHyb")  == 0) {
      XSNPChipHyb  *scset = (XSNPChipHyb*)set;
      fNRows      = scset->GetNumRows();
      fNCols      = scset->GetNumColumns();
      fMinInten   = scset->GetMinIntensity();
      fMaxInten   = scset->GetMaxIntensity();
      fNMinInten  = scset->GetNumMinIntensity();
      fNMaxInten  = scset->GetNumMaxIntensity();
      fMaxNPixels = scset->GetMaxNumPixels(); 
   } else if (strcmp(set->ClassName(), "XGenomeChipHyb")  == 0) {
      XGenomeChipHyb  *gmset = (XGenomeChipHyb*)set;
      fNRows      = gmset->GetNumRows();
      fNCols      = gmset->GetNumColumns();
      fMinInten   = gmset->GetMinIntensity();
      fMaxInten   = gmset->GetMaxIntensity();
      fNMinInten  = gmset->GetNumMinIntensity();
      fNMaxInten  = gmset->GetNumMaxIntensity();
      fMaxNPixels = gmset->GetMaxNumPixels(); 
   } else if (strcmp(set->ClassName(), "XExonChipHyb")  == 0) {
      XExonChipHyb  *exset = (XExonChipHyb*)set;
      fNRows      = exset->GetNumRows();
      fNCols      = exset->GetNumColumns();
      fMinInten   = exset->GetMinIntensity();
      fMaxInten   = exset->GetMaxIntensity();
      fNMinInten  = exset->GetNumMinIntensity();
      fNMaxInten  = exset->GetNumMaxIntensity();
      fMaxNPixels = exset->GetMaxNumPixels(); 
   } else if (strcmp(set->ClassName(), "XGenePixHyb")  == 0) {
      XGenePixHyb *gpset = (XGenePixHyb*)set;
      fNRows      = gpset->GetNumRows();
      fNCols      = gpset->GetNumColumns();
      fMinInten   = gpset->GetMinIntensity();
      fMaxInten   = gpset->GetMaxIntensity();
      fNMinInten  = gpset->GetNumMinIntensity();
      fNMaxInten  = gpset->GetNumMaxIntensity();
   } else {
//??
   }//if
}//AddUserInfo

//______________________________________________________________________________
void XDataTreeInfo::AddUserInfo(Int_t nrows, Int_t ncols, Int_t nmin, Double_t min,
                    Int_t nmax, Double_t max, Short_t maxnpix)
{
   // Add user info from tree set
   if(kCS) cout << "------XDataTreeInfo::AddUserInfo------" << endl;

   fNRows      = nrows;
   fNCols      = ncols;
   fMinInten   = min;
   fMaxInten   = max;
   fNMinInten  = nmin;
   fNMaxInten  = nmax;
   fMaxNPixels = maxnpix; 
}//AddUserInfo

//______________________________________________________________________________
void XDataTreeInfo::AddUserInfo(Int_t nquant, Double_t *q, Double_t *quant)
{
   // Add user info from tree set
   if(kCS) cout << "------XDataTreeInfo::AddUserInfo------" << endl;

   if (nquant > fNQuantiles) {
      if (fIntenQuant) {delete [] fIntenQuant; fIntenQuant = 0;}
      if (fQuantiles)  {delete [] fQuantiles;  fQuantiles  = 0;}

      fQuantiles  = new Double_t[nquant];
      fIntenQuant = new Double_t[nquant];
   }//if

   fNQuantiles = nquant;

   memcpy(fQuantiles,  q,     nquant*sizeof(Double_t));
   memcpy(fIntenQuant, quant, nquant*sizeof(Double_t));
}//AddUserInfo

//______________________________________________________________________________
Double_t XDataTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCSa) cout << "------XDataTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNRows") == 0) {
      return fNRows;
   } else if (strcmp(name, "fNCols")      == 0) {
      return fNCols;
   } else if (strcmp(name, "fMinInten")   == 0) {
      return fMinInten;
   } else if (strcmp(name, "fMaxInten")   == 0) {
      return fMaxInten;
   } else if (strcmp(name, "fNMinInten")  == 0) {
      return fNMinInten;
   } else if (strcmp(name, "fNMaxInten")  == 0) {
      return fNMaxInten;
   } else if (strcmp(name, "fMaxNPixels") == 0) {
      return fMaxNPixels;
   } else if (strcmp(name, "fNQuantiles") == 0) {
      return fNQuantiles;
   }//if
   return 0;
}//GetValue

//______________________________________________________________________________
Double_t *XDataTreeInfo::GetQuantiles()
{
   // Return quantiles array
   if(kCS) cout << "------XDataTreeInfo::GetQuantiles------" << endl;

   if (fQuantiles == 0) return 0;

   return fQuantiles;
}//GetQuantiles

//______________________________________________________________________________
Double_t *XDataTreeInfo::GetIntenQuantiles()
{
   // Return quantiles array for intensities
   if(kCS) cout << "------XDataTreeInfo::GetIntenQuantiles------" << endl;

   if (fIntenQuant == 0) return 0;

   return fIntenQuant;
}//GetIntenQuantiles


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMaskTreeInfo                                                        //
//                                                                      //
// Class containing info about mask tree stored in fUserInfo            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMaskTreeInfo::XMaskTreeInfo() 
              :XTreeInfo()
{
   // Default MaskTreeInfo constructor
   if(kCS) cout << "---XMaskTreeInfo::XMaskTreeInfo(default)------" << endl;

   fNRows  = 0; 
   fNCols  = 0;
   fNFlags = 0;
}//Constructor

//______________________________________________________________________________
XMaskTreeInfo::XMaskTreeInfo(const char *name, const char *title) 
              :XTreeInfo(name, title)
{
   // Normal MaskTreeInfo constructor
   if(kCS) cout << "---XMaskTreeInfo::XMaskTreeInfo------" << endl;

   fNRows  = 0; 
   fNCols  = 0;
   fNFlags = 0;
}//Constructor

//______________________________________________________________________________
XMaskTreeInfo::~XMaskTreeInfo()
{
   // MaskTreeInfo destructor
   if(kCS) cout << "---XMaskTreeInfo::~XMaskTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XMaskTreeInfo::AddUserInfo(Int_t nrows, Int_t ncols, Int_t nflags)
{
   // Add user info from tree set
   if(kCS) cout << "------XMaskTreeInfo::AddUserInfo------" << endl;

   fNRows  = nrows;
   fNCols  = ncols;
   fNFlags = nflags;
}//AddUserInfo

//______________________________________________________________________________
Double_t XMaskTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCS) cout << "------XMaskTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNRows") == 0) {
      return fNRows;
   } else if (strcmp(name, "fNCols")  == 0) {
      return fNCols;
   } else if (strcmp(name, "fNFlags") == 0) {
      return fNFlags;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataSetting                                                         //
//                                                                      //
// Class for initialization of input dat settings                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDataSetting::XDataSetting()
             :XSetting()
{
   // Default DataSetting constructor
   if(kCS) cout << "---XDataSetting::XDataSetting(default)------" << endl;

   fSchemeFile = 0;
   fSchemeName = "";
   fDataType   = "";
   fInputType  = "";
   fNVar       = 0;
   fVars       = 0;
   fTyps       = 0;
}//Constructor

//______________________________________________________________________________
XDataSetting::XDataSetting(const char *arraytype, const char *infile)
             :XSetting(arraytype, infile)
{
   // Normal DataSetting constructor
   if(kCS) cout << "---XDataSetting::XDataSetting------" << endl;

   fSchemeFile = 0;
   fSchemeName = "";
   fDataType   = "def";
   fInputType  = "RawData";
   fNVar       = 0;
   fVars       = 0;
   fTyps       = 0;
}//Constructor

//______________________________________________________________________________
XDataSetting::~XDataSetting()
{
   // DataSetting destructor
   if(kCS) cout << "---XDataSetting::~XDataSetting------" << endl;

   if (fVars) {delete [] fVars; fVars = 0;}
   if (fTyps) {delete [] fTyps; fTyps = 0;}

   fNVar       = 0;
   fSchemeFile = 0;
}//Destructor

//______________________________________________________________________________
Int_t XDataSetting::InitInput(const char *schemename, const char *datatype,
                    const char *varlist, const char *inputtype)
{
   // Initialize input data with schemename being exact name of array type
   // datatype should describe type of chip data, e.g. cel for rawdata,
   // tbw, mdp for Metrics or PivotData (but also mas5 or rma)
   // varlist is list of variables to be stored in tree(s) and should contain
   // type of variable, e.g. "var1/I:var2/C:var3/D"
   // inputtype is rawdata (default), Metrics, PivotData
   if(kCS) cout << "------XDataSetting::InitInput------" << endl;

   fSchemeName = schemename;
   fDataType   = datatype;
   fInputType  = inputtype;

   fNVar = NumSeparators(varlist, ":") + 1;
   if (fVars) {delete [] fVars; fVars = 0;}
   if (fTyps) {delete [] fTyps; fTyps = 0;}
   if (!(fVars = new (nothrow) TString[fNVar])) return errInitMemory;
   if (!(fTyps = new (nothrow) TString[fNVar])) return errInitMemory;

   for (Int_t i=0; i<fNVar; i++) {
      fVars[i] = SubString(varlist, ":", i);
      fTyps[i] = Path2Name((fVars[i]).Data(),"/","");
      if (fTyps[i] == "") fTyps[i] = "D";
      fVars[i] = Path2Name((fVars[i]).Data(),"","/");
   }//for_i

   return errNoErr;
}//InitInput


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybridization                                                       //
//                                                                      //
// Base class for microarray hybridization                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XHybridization::XHybridization()
               :XTreeSet()
{
   // Default Hybridization constructor
   if(kCS) cout << "---XHybridization::XHybridization(default)------" << endl;

   fSchemeFile   = 0;
   fDataTree     = 0;
   fMaskTree     = 0;
   fDataName     = "";
   fSchemeName   = "";
   fDataTreeName = "";
   fMaskTreeName = "";
   fNRows        = 0; 
   fNCols        = 0;
   fNCells       = 0;
}//Constructor

//______________________________________________________________________________
XHybridization::XHybridization(const char *name, const char *title)
               :XTreeSet(name, title)
{
   // Normal Hybridization constructor
   if(kCS) cout << "---XHybridization::XHybridization------" << endl;

   fSchemeFile   = 0;
   fDataTree     = 0;
   fMaskTree     = 0;
   fDataName     = "";
   fSchemeName   = "";
   fDataTreeName = "NA";
   fMaskTreeName = "NA";
   fNRows        = 0; 
   fNCols        = 0;
   fNCells       = 0;
}//Constructor

//______________________________________________________________________________
XHybridization::~XHybridization()
{
   // DataManager destructor
   if(kCS) cout << "---XHybridization::~XHybridization------" << endl;

   fSchemeFile   = 0;
   fDataTree     = 0;
   fMaskTree     = 0;
}//Destructor

//______________________________________________________________________________
Int_t XHybridization::InitFiles(TFile *schemefile)
{
   // Initialize root scheme file
   if(kCS) cout << "------XHybridization::InitFiles------" << endl;

   fSchemeFile = schemefile;

   if (fSchemeFile) return errNoErr;

   cerr << "Error: Scheme file was not initialized" << endl;
   return errAbort;
}//InitFiles


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipHyb                                                         //
//                                                                      //
// Class for GeneChip hybridization                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGeneChipHyb::XGeneChipHyb()
             :XHybridization()
{
   // Default GeneChipHyb constructor
   if(kCS) cout << "---XGeneChipHyb::XGeneChipHyb(default)------" << endl;

   fMinInten   = 0;
   fMaxInten   = 0;
   fNMinInten  = 0;
   fNMaxInten  = 0;
   fMaxNPixels = 0;
}//Constructor

//______________________________________________________________________________
XGeneChipHyb::XGeneChipHyb(const char *name, const char *title)
             :XHybridization(name, title)
{
   // Normal GeneChipHyb constructor
   if(kCS) cout << "---XGeneChipHyb::XGeneChipHyb------" << endl;

   fMinInten   = 0;
   fMaxInten   = 0;
   fNMinInten  = 0;
   fNMaxInten  = 0;
   fMaxNPixels = 0;
}//Constructor

//______________________________________________________________________________
XGeneChipHyb::~XGeneChipHyb()
{
   // GeneChipHyb destructor
   if(kCS) cout << "---XGeneChipHyb::~XGeneChipHyb------" << endl;

}//Destructor

//______________________________________________________________________________
void XGeneChipHyb::AddDataTreeInfo(TTree *tree, const char *name, Option_t *option,
                   Int_t nquant, Double_t *q, Double_t *quant)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGeneChipHyb::AddDataTreeInfo------" << endl;

   XDataTreeInfo *info = new XDataTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   if (nquant > 0) info->AddUserInfo(nquant, q, quant);

   tree->GetUserInfo()->Add(info);
}//AddDataTreeInfo

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportTreeInfo(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export selected variables given in varlist which are stored in tree
   // treename in output
   if(kCS) cout << "------XGeneChipHyb::ExportTreeInfo------" << endl;

   Int_t err = errNoErr;

   // remove "userinfo" from varlist
   TString infolist = RemoveSubString(varlist, "userinfo:", kFALSE);

   if (strcmp(exten,"cel") == 0) {
      err = this->ExportDataTreeInfo(n, names, infolist, output, sep);
   } else if (strcmp(exten,"msk") == 0) {
      err = this->ExportMaskTreeInfo(n, names, infolist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if

   return err;
}//ExportTreeInfo

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportTreeType(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export selected variables given in varlist which are stored in tree
   // treename in output
   if(kCS) cout << "------XGeneChipHyb::ExportTreeType------" << endl;

   Int_t err = errNoErr;
   if ((n == 1) && (strcmp(varlist,"*") == 0)) {
   // ExportDataTree takes about 20sec to export
      if (strcmp(exten,"cel") == 0) {
         err = this->ExportDataTree(names, output, sep);
      } else if (strcmp(exten,"msk") == 0) {
         err = this->ExportMaskTree(names, output, sep);
      } else {
         return fManager->HandleError(errExtension, exten);
      }//if
   } else {
   // ExportDataTrees(n=1) takes about 40sec to export
      if (strcmp(exten,"cel") == 0) {
         err = this->ExportDataTrees(n, names, varlist, output, sep);
      } else if (strcmp(exten,"msk") == 0) {
         err = this->ExportMaskTrees(n, names, varlist, output, sep);
      } else {
         return fManager->HandleError(errExtension, exten);
      }//if
   }//if

   return err;
}//ExportTreeType

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                    const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XGeneChipHyb::ExportTreeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   exten = 0; n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
Int_t XGeneChipHyb::Import(ifstream &input, Option_t *option, const char *sep,
                    char delim, Int_t split)
{
   // Import data from input. Option is tree extension. 
   if(kCS) cout << "------XGeneChipHyb::Import------" << endl;

   Int_t err = errNoErr;

   if (IsXDAFile(input)) {
      if (!err) err = ReadXDAHeader(input, sep, delim);
      if (!err) err = ReadXDAData(input, option, sep, delim, split);
   } else if (IsCalvinFile(input)) {
      if (!err) err = ReadCalvinGenericFile(input, option, split);
   } else {
      if (!err) err = ReadHeader(input, sep, delim);
      if (!err) err = ReadData(input, option, sep, delim, split);
   }//if

   return err;
}//Import

//______________________________________________________________________________
void XGeneChipHyb::PrintInfo()
{
   // Print GeneChip hybridization information
   if(kCS) cout << "------XGeneChipHyb::PrintInfo------" << endl;

   if (fgPrintHeader) {
      cout << "==============================================================================" << endl;
      cout << setw(14) << "Hybridization" << setw(12) << "DataName"
           << setw(17) << "ChipName"      << setw(17) << "DataTree" 
           << setw(17) << "MaskTree"      << setw(9)  << "Rows" 
           << setw(9)  << "Columns"       << setw(12) << "MinInten"
           << setw(12) << "#MinInten"     << setw(12) << "MaxInten"
           << setw(12) << "#MaxInten"     << setw(12) << "MaxNumPix" << endl;
      cout << "==============================================================================" << endl;
      fgPrintHeader = kFALSE;
   }//if

   cout << setw(14) << this->GetName()      << setw(12) << fDataName.Data()
        << setw(17) << fSchemeName.Data()   << setw(17) << fDataTreeName.Data() 
        << setw(17) << fMaskTreeName.Data() << setw(9)  << fNRows 
        << setw(9)  << fNCols               << setw(12) << fMinInten 
        << setw(12) << fNMinInten           << setw(12) << fMaxInten
        << setw(12) << fNMaxInten           << setw(12) << fMaxNPixels << endl;
   cout << "------------------------------------------------------------------------------" << endl;
}//PrintInfo

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportDataTrees(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGeneChipHyb::ExportDataTrees------" << endl;

// Decompose varlist
   Bool_t hasMean = kFALSE;
   Bool_t hasStdv = kFALSE;
   Bool_t hasNPix = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasMean = kTRUE;
      hasStdv = kTRUE;
      hasNPix = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name, varlist), ":");
      while(name) {
         if (strcmp(name,"fInten")   == 0) {hasMean = kTRUE;}
         if (strcmp(name,"fStdev")   == 0) {hasStdv = kTRUE;}
         if (strcmp(name,"fNPixels") == 0) {hasNPix = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree   **tree = new TTree*[n];
   XGCCell **cell = new XGCCell*[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         cell[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("DataBranch", &cell[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         cell[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("DataBranch", &cell[k]);
      }//for_k
   }//if

// Output header
   output << "X" << sep << "Y";
   if (n > 1) {
      for (Int_t k=0; k<n; k++) {
         if (hasMean) output << sep << (names[k] + "_MEAN").Data();
         if (hasStdv) output << sep << (names[k] + "_STDV").Data();
         if (hasNPix) output << sep << (names[k] + "_NPIXELS").Data();
      }//for_k
   } else {
      if (hasMean) output << sep << "MEAN";
      if (hasStdv) output << sep << "STDV";
      if (hasNPix) output << sep << "NPIXELS";
   }//if
   output << endl;

// Loop over tree entries and tree branches
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);
         if (k == 0)  output << cell[k]->GetX() << sep << cell[k]->GetY();
         if (hasMean) output << sep << cell[k]->GetIntensity();
         if (hasStdv) output << sep << cell[k]->GetStdev();
         if (hasNPix) output << sep << cell[k]->GetNumPixels();
      }//for_j
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported." << endl;
   }//if

   delete [] cell;
   delete [] tree;

   return errNoErr;
}//ExportDataTrees

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportDataTree(TString *name, ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGeneChipHyb::ExportDataTree------" << endl;

// Output header
   output << "X" << sep << "Y" << sep << "MEAN" << sep << "STDV" << sep
          << "NPIXELS" << endl;

// Get tree
   XGCCell *cell = 0;
   fTree = (TTree*)gDirectory->Get((name[0]).Data());
   if (!fTree) return errGetTree;
   fTree->SetBranchAddress("DataBranch", &cell);

// Loop over all entries of tree
   Int_t nentries = (Int_t)(fTree->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << cell->GetX()         << sep
             << cell->GetY()         << sep
             << cell->GetIntensity() << sep 
             << cell->GetStdev()     << sep 
             << cell->GetNumPixels() << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported." << endl;
   }//if

   return errNoErr;
}//ExportDataTree

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportMaskTrees(Int_t n, TString *names, const char * /*varlist*/,
                    ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGeneChipHyb::ExportMaskTrees------" << endl;

// Get trees
   TTree     **tree = new TTree*[n];
   XCellMask **mask = new XCellMask*[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         mask[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("MaskBranch", &mask[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         mask[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("MaskBranch", &mask[k]);
      }//for_k
   }//if

// Output header
   output << "X" << sep << "Y";
   if (n > 1) {
      for (Int_t k=0; k<n; k++) {
         output << sep << (names[k] + "_FLAG").Data();
      }//for_k
   } else {
      output << sep << "FLAG";
   }//if
   output << endl;

// Loop over tree entries and tree branches
   Int_t nentries = (Int_t)(tree[0]->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);
         if (k == 0)  output << mask[k]->GetX() << sep << mask[k]->GetY();
         output << sep << mask[k]->GetFlag();
      }//for_j
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported." << endl;
   }//if

   delete [] mask;
   delete [] tree;

   return errNoErr;
}//ExportMaskTrees

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportMaskTree(TString *name, ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGeneChipHyb::ExportMaskTree------" << endl;

// Output header
   output << "X" << sep << "Y" << sep << "FLAG" << endl;

// Get tree
   XCellMask *mask = 0;
   TTree     *tree = (TTree*)gDirectory->Get((name[0]).Data()); 
   if (tree == 0) return errGetTree;
   tree->SetBranchAddress("MaskBranch", &mask);

// Loop over all entries of tree
   Int_t nentries = (Int_t)(tree->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
      output << mask->GetX()    << sep
             << mask->GetY()    << sep
             << mask->GetFlag() << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported." << endl;
   }//if

   tree->Delete("");
   tree = 0;

   return errNoErr;
}//ExportMaskTree

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportDataTreeInfo(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in userinfo of data tree to file output
   if(kCS) cout << "------XGeneChipHyb::ExportDataTreeInfo------" << endl;

// Decompose varlist
   Bool_t hasTreeName  = kFALSE;
   Bool_t hasSetName   = kFALSE;
   Bool_t hasOption    = kFALSE;
   Bool_t hasNRows     = kFALSE;
   Bool_t hasNCols     = kFALSE;
   Bool_t hasMinInten  = kFALSE;
   Bool_t hasMaxInten  = kFALSE;
   Bool_t hasNMinInten = kFALSE;
   Bool_t hasNMaxInten = kFALSE;
   Bool_t hasMaxNPix   = kFALSE;
   Bool_t hasNQuant    = kFALSE;
   Bool_t hasQuant     = kFALSE;
   Bool_t hasInten     = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasTreeName  = kTRUE;
      hasSetName   = kTRUE;
      hasOption    = kTRUE;
      hasNRows     = kTRUE;
      hasNCols     = kTRUE;
      hasMinInten  = kTRUE;
      hasMaxInten  = kTRUE;
      hasNMinInten = kTRUE;
      hasNMaxInten = kTRUE;
      hasMaxNPix   = kTRUE;
      hasNQuant    = kTRUE;
      hasQuant     = kTRUE;
      hasInten     = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fName")       == 0) {hasTreeName  = kTRUE;}
         if (strcmp(name,"fSetName")    == 0) {hasSetName   = kTRUE;}
         if (strcmp(name,"fOption")     == 0) {hasOption    = kTRUE;}
         if (strcmp(name,"fNRows")      == 0) {hasNRows     = kTRUE;}
         if (strcmp(name,"fNCols")      == 0) {hasNCols     = kTRUE;}
         if (strcmp(name,"fMinInten")   == 0) {hasMinInten  = kTRUE;}
         if (strcmp(name,"fMaxInten")   == 0) {hasMaxInten  = kTRUE;}
         if (strcmp(name,"fNMinInten")  == 0) {hasNMinInten = kTRUE;}
         if (strcmp(name,"fNMaxInten")  == 0) {hasNMaxInten = kTRUE;}
         if (strcmp(name,"fMaxNPixels") == 0) {hasMaxNPix   = kTRUE;}
         if (strcmp(name,"fNQuantiles") == 0) {hasNQuant    = kTRUE;}
         if (strcmp(name,"fQuantiles")  == 0) {hasQuant     = kTRUE;}
         if (strcmp(name,"fIntenQuant") == 0) {hasInten     = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree         **tree = new TTree*[n];
   XDataTreeInfo **info = new XDataTreeInfo*[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         info[k] = (XDataTreeInfo*)tree[k]->GetUserInfo()->At(0);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         info[k] = (XDataTreeInfo*)tree[k]->GetUserInfo()->At(0);
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

   if (hasNRows) {
      output << "Rows";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("fNRows"));
      output << endl;
   }//if

   if (hasNCols) {
      output << "Cols";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("fNCols"));
      output << endl;
   }//if

   if (hasMinInten) {
      output << "MinIntensity";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetValue("fMinInten");
      output << endl;
   }//if

   if (hasMaxInten) {
      output << "MaxIntensity";
      for (Int_t k=0; k<n; k++) output << sep << info[k]->GetValue("fMaxInten");
      output << endl;
   }//if

   if (hasNMinInten) {
      output << "NumMinIntensity";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("fNMinInten"));
      output << endl;
   }//if

   if (hasNMaxInten) {
      output << "NumMaxIntensity";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("fNMaxInten"));
      output << endl;
   }//if

   if (fMaxNPixels) {
      output << "MaxNumPixels";
      for (Int_t k=0; k<n; k++) output << sep << (Int_t)(info[k]->GetValue("hasMaxNPix"));
      output << endl;
   }//if

   if (hasNQuant) {
      output << "NQuantiles";
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
         TString str = "Quantile";

         output << (str+=i).Data();
         for (Int_t k=0; k<n; k++) output << sep << quant[k][i];
         output << endl;
      }//for_i

      delete [] quant;
   }//if

   if (hasInten) {
      Double_t **quant = new Double_t*[n];
      Double_t **inten = new Double_t*[n];
      for (Int_t k=0; k<n; k++) {
         quant[k] = info[k]->GetQuantiles();
         inten[k] = info[k]->GetIntenQuantiles();
      }//for_k

      Int_t nq  = (Int_t)(info[0]->GetValue("fNQuantiles"));
      for (Int_t i=0; i<nq; i++) {
         TString str; str.Form("Intensity_Q%4.2f", quant[0][i]);

         output << str.Data();
         for (Int_t k=0; k<n; k++) output << sep << inten[k][i];
         output << endl;
      }//for_i

      delete [] inten;
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
}//ExportDataTreeInfo

//______________________________________________________________________________
Int_t XGeneChipHyb::ExportMaskTreeInfo(Int_t n, TString *names, const char *varlist,
                    ofstream &output, const char *sep)
{
   // Export data stored in userinfo of mask tree to file output
   if(kCS) cout << "------XGeneChipHyb::ExportMaskTreeInfo------" << endl;

   Int_t err = errNoErr;

   return errNoErr;
}//ExportMaskTreeInfo

//______________________________________________________________________________
Int_t XGeneChipHyb::IsXDAFile(ifstream &input)
{
   if(kCS) cout << "------XGeneChipHyb::IsXDAFile------" << endl;

   // read magic number from input.
   Int_t magic = 0;
   READ_INT(input, magic);

   // need to rewind input to start
   input.seekg(ios::beg);

   return (magic == kMagicXDANumber);
}//IsXDAFile

//______________________________________________________________________________
Int_t XGeneChipHyb::IsCalvinFile(ifstream &input)
{
   if(kCS) cout << "------XGeneChipHyb::IsCalvinFile------" << endl;

   // read magic number from input.
	unsigned char  magic = 0;
	READ_UCHAR(input, magic);

   // need to rewind input to start
   input.seekg(ios::beg);

   return (magic == kMagicCalvinNumber);
}//IsCalvinFile

//______________________________________________________________________________
TString XGeneChipHyb::ChipType(const char *header, Int_t toUpper)
{
   // Determine chip type, i.e. scheme name from header
   // For toUpper = 1 convert first letter to uppercase
   if(kCS) cout << "------XGeneChipHyb::ChipType------" << endl;

   TString name = "";
   if (strstr(header, "DatHeader")) {
      name = strstr(header, "DatHeader");
   } else {
      name = TString(header);
   }//if

   Int_t i;
   for(Int_t j=0;j<2;j++){
      i = name.Index("\x14");
      name = &name[i+1];
   }//for_j

   // remove non-printing chars
   name = RemoveEnds(name);

   // remove extension ".lsq"
   i = name.Index(".");
   name.Remove(i);

   // convert first letter to uppercase
   if (toUpper == 1) {
      char *s = (char*)(name.Data());
      *s = toupper(s[0]);
   }//if

   return name;
}//ChipType

//______________________________________________________________________________
Int_t XGeneChipHyb::CheckChipType(const char *header, const char *name)
{
   // Get scheme name from header, compare with name, and store in fSchemeName
   // Check if first letter is uppercase/lowercase, e.g. test3 or miRNA
   if(kCS) cout << "------XGeneChipHyb::CheckChipType------" << endl;

   TString chipname = this->ChipType(header, 0);
   // check if chipname is identical to name of imported scheme
   if ((strcmp(name, "") != 0) && (strcmp(name, chipname.Data()) != 0)) {
         chipname = this->ChipType(header, 1);
         if ((strcmp(name, "") != 0) && (strcmp(name, chipname.Data()) != 0)) {
            return errChipType;
         }//if
   }//if
   fSchemeName = chipname;

   return errNoErr;
}//CheckChipType

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadHeader(ifstream &input, const char * /*sep*/, char delim)
{
   // Read header from input. 
   if(kCS) cout << "------XGeneChipHyb::ReadHeader------" << endl;

   char  nextline[kBufSize];
   Int_t i;
   Int_t err = errNoErr;

// Check for first line [CEL]
   input.getline(nextline, kBufSize, delim);
   if ( strncmp("[CEL]", nextline, 5) != 0) return errHeaderLine;

// Check version of CEL-file
   input.getline(nextline, kBufSize, delim);
   if (!input.good()) return errReadingInput;
   if ( strncmp("Version=", nextline, 8) == 0) {
      sscanf(&(nextline[8]), "%d", &i);
      if (i != 3) {
         TString str = ""; str += i;
         return fManager->HandleError(errCELVersion, str);
      }//if
   }//if

// Check for header line [HEADER]
   while (strncmp("[HEADER]", nextline, 8) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   // get column number
   input.getline(nextline, kBufSize, delim);
   if ( strncmp("Cols=", nextline, 5) != 0) return errMissingLine;
   sscanf(&(nextline[5]), "%d", &fNCols);

   // get row number
   input.getline(nextline, kBufSize, delim);
   if ( strncmp("Rows=", nextline, 5) != 0) return errMissingLine;
   sscanf(&(nextline[5]), "%d", &fNRows);

// Get chip type from DatHeader
   while (strncmp("DatHeader=", nextline, 10) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   // need to check for correct scheme name
   TString scmname = ((XDataSetting*)fSetting)->GetSchemeName();
   if ((err = this->CheckChipType(&nextline[0], scmname)) != errNoErr) {
      return fManager->HandleError(err, scmname, fSchemeName);
   }//if

//?? Maybe check algorithm?

// Check for intensity line [INTENSITY]
   while (strncmp("[INTENSITY]", nextline, 11) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   // get number of cells
   input.getline(nextline, kBufSize, delim);
   if (strncmp("NumberCells=", nextline, 12) != 0) return errMissingLine;
   sscanf(&(nextline[12]), "%d", &fNCells);

   return err;
}//ReadHeader

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadData(ifstream &input, Option_t *option, const char * /*sep*/,
                    char delim, Int_t split)
{
   // Read data from input and store in data tree. 
   if(kCS) cout << "------XGeneChipHyb::ReadData------" << endl;

   char     nextline[kBufSize];
   Int_t    i, x, y;
   Short_t  numpix;
   Double_t inten, stdev;
   Int_t    err    = 0;
   Int_t    numcel = 0;
   Double_t min    = DBL_MAX;  //defined in float.h
   Double_t max    = 0;
   Int_t    nummin = 0;
   Int_t    nummax = 0;
   Short_t  maxpix = 0;

   Int_t     nquant = 7;
   Double_t  q[]    = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};
   Double_t *quantI = 0; 
   if (!(quantI = new (nothrow) Double_t[nquant])) return errInitMemory;

// Create data tree
   TString exten = Path2Name(option, ".", "");
   fDataTreeName = fTreeName + "." + exten;
   TTree   *datatree = new TTree(fDataTreeName, fSchemeName);
   if (datatree == 0) return errCreateTree;

   XGCCell *cell    = new XGCCell();
   Int_t    bufsize = XManager::GetBufSize();
   datatree->Branch("DataBranch", "XGCCell", &cell, bufsize, split);

// Read header line containing column names
   input.getline(nextline, kBufSize, delim);
   if (strncmp("CellHeader=", nextline, 11) != 0) return errMissingLine;

// Read data and store in data tree
   for (i=0; i<fNCells; i++) {
      input.getline(nextline, kBufSize, delim);
      if (!input.good()) {err = errPrematureEOF; break;}
      sscanf(nextline,"%i %i %lf %lf %hi", &x, &y, &inten, &stdev, &numpix);

      // number of cells with minimal intensity
      if      (inten <  min) {min = inten; nummin = 1;}
      else if (inten == min) {nummin++;}

      // number of cells with maximal intensity
      if      (inten >  max) {max = inten; nummax = 1;}
      else if (inten == max) {nummax++;}

      // maximal pixel number
      if (numpix > maxpix) {maxpix = numpix;}

      // fill data tree
      cell->SetX(x);
      cell->SetY(y);
      cell->SetIntensity(inten);
      cell->SetStdev(stdev);
      cell->SetNumPixels(numpix);
      datatree->Fill(); 

      numcel++;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << fNCells << "> records imported." << endl;
   }//if

   fMinInten   = min;
   fNMinInten  = nummin;
   fMaxInten   = max;
   fNMaxInten  = nummax;
   fMaxNPixels = maxpix;

   if (XManager::fgVerbose) {
      cout << "   hybridization statistics: " << endl;
      cout << "      " << nummin << " cells with minimal intensity " << min << endl;
      cout << "      " << nummax << " cells with maximal intensity " << max << endl;
   }//if

// Check for problems with intensities
   if (max <= min) {
      cout << "Warning: maximal intensity equal or less than minimal intensity!" << endl;
//      cout << "Thus CEL-file will not be imported as <" << fDataTreeName << ">!" << endl;
//      goto cleanup;
   }//if

// Write data tree to file if all data are read
   if (numcel == fNCells) {
      // quantiles for residual trees
      err = this->DataQuantiles(datatree, cell, nquant, q, quantI);
      if (err != errNoErr) goto cleanup;

      // add tree info to tree
//      AddDataTreeInfo(datatree, fDataTreeName);
      AddDataTreeInfo(datatree, datatree->GetName(), "txt", nquant, q, quantI);

      // write data tree to file 
      if ((err = WriteTree(datatree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(datatree->GetName(), 0);
      }//if
   } else {
      fDataTreeName = "NA";
      err = errReadingInput;
      cerr << "Error: number of lines read <" << numcel 
           << "> is not equal to to number of cells <" << fNCells << ">"
           << endl;
   }//if

// Delete data tree from RAM
cleanup:
   datatree->Delete("");
   datatree = 0;
   delete cell;

   if (quantI) {delete [] quantI; quantI = 0;}

   return err;
}//ReadData

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadXDAHeader(ifstream &input, const char * /*sep*/, char /*delim*/)
{
   // Read header from XDA binary input. 
   if(kCS) cout << "------XGeneChipHyb::ReadXDAHeader------" << endl;

   Int_t err = errNoErr;

// Read magic number
   Int_t magic;
   READ_INT(input, magic);
   if (magic != kMagicXDANumber) {
      TString str = ""; str += magic;
      return fManager->HandleError(errCELVersion, str);
   }//if

// Read version
   Int_t version;
   READ_INT(input, version);
   if (version != 4) {
      TString str = ""; str += version;
      return fManager->HandleError(errCELVersion, str);
   }//if

// Read array dimensions
   READ_INT(input, fNRows);
   READ_INT(input, fNCols);
   READ_INT(input, fNCells);

// Read other array parameters
   TString header;
   char   *str = NULL;

   // read DatHeader
   READ_STRING(input, str);
   header = str;
   delete[] str; str = NULL;

   // read Algorithm
   READ_STRING(input, str);
   delete[] str; str = NULL;

   // read AlgorithmParameters
   READ_STRING(input, str);
   delete[] str; str = NULL;

   // read Margin
   int margin;
   READ_INT(input, margin);

   // read Outliers
   unsigned int outliers;
   READ_UINT(input, outliers);

   // read Masked
   unsigned int masked;
   READ_UINT(input, masked);

   // read SubGrids
   int subgrids;
   READ_INT(input, subgrids);

//??   ParseCorners(header);

// Get scheme name from DatHeader
   TString scmname = ((XDataSetting*)fSetting)->GetSchemeName();
   if ((err = this->CheckChipType(header, scmname)) != errNoErr) {
      return fManager->HandleError(err, scmname, fSchemeName);
   }//if

   return err;
}//ReadXDAHeader

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadXDAData(ifstream &input, Option_t *option, const char * /*sep*/,
                    char /*delim*/, Int_t split)
{
   // Read data from XDA binary input and store in data tree. 
   if(kCS) cout << "------XGeneChipHyb::ReadXDAData------" << endl;

   Int_t    x, y;
   Short_t  numpix;
   Double_t inten, stdev;
   Int_t    err    = 0;
   Int_t    numcel = 0;
   Double_t min    = DBL_MAX;  //defined in float.h
   Double_t max    = 0;
   Int_t    nummin = 0;
   Int_t    nummax = 0;
   Short_t  maxpix = 0;

   Int_t     nquant = 7;
   Double_t  q[]    = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};
   Double_t *quantI = 0; 
   if (!(quantI = new (nothrow) Double_t[nquant])) return errInitMemory;

// Create data tree
   TString exten = Path2Name(option, ".", "");
   fDataTreeName = fTreeName + "." + exten;
   TTree *datatree = new TTree(fDataTreeName, fSchemeName);
   if (datatree == 0) return errCreateTree;

   XGCCell *cell    = new XGCCell();
   Int_t    bufsize = XManager::GetBufSize();
   datatree->Branch("DataBranch", "XGCCell", &cell, bufsize, split);

// Read the intensity data
   CELEntry *entries = new CELEntry[fNCells];
   for (Int_t i=0; i<fNCells; i++) {
      CELEntry *entry = entries + i;

      READ_FLOAT(input, entry->Intensity);
      READ_FLOAT(input, entry->Stdv);
      READ_SHORT(input, entry->Pixels);

      x      = Index2X(i);
      y      = Index2Y(i);
      inten  = (Double_t)entry->Intensity;
      stdev  = (Double_t)entry->Stdv;
      numpix = entry->Pixels;

      // number of cells with minimal intensity
      if      (inten <  min) {min = inten; nummin = 1;}
      else if (inten == min) {nummin++;}

      // number of cells with maximal intensity
      if      (inten >  max) {max = inten; nummax = 1;}
      else if (inten == max) {nummax++;}

      // maximal pixel number
      if (numpix > maxpix) {maxpix = numpix;}

      // fill data tree
      cell->SetX(x);
      cell->SetY(y);
      cell->SetIntensity(inten);
      cell->SetStdev(stdev);
      cell->SetNumPixels(numpix);
      datatree->Fill(); 

      numcel++;
      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << fNCells << "> records imported." << endl;
   }//if

	if (entries) {delete[] entries; entries = 0;}

   fMinInten   = min;
   fNMinInten  = nummin;
   fMaxInten   = max;
   fNMaxInten  = nummax;
   fMaxNPixels = maxpix;

   if (XManager::fgVerbose) {
      cout << "   hybridization statistics: " << endl;
      cout << "      " << nummin << " cells with minimal intensity " << min << endl;
      cout << "      " << nummax << " cells with maximal intensity " << max << endl;
   }//if

// Check for problems with intensities
   if (max <= min) {
      cout << "Warning: maximal intensity equal or less than minimal intensity!" << endl;
//      cout << "Thus CEL-file will not be imported as <" << fDataTreeName << ">!" << endl;
//      goto cleanup;
   }//if

// Write data tree to file if all data are read
   if (numcel == fNCells) {
      // quantiles for residual trees
      err = this->DataQuantiles(datatree, cell, nquant, q, quantI);
      if (err != errNoErr) goto cleanup;

      // add tree info to tree
//      AddDataTreeInfo(datatree, fDataTreeName);
      AddDataTreeInfo(datatree, datatree->GetName(), "xda", nquant, q, quantI);

      // write data tree to file 
      if ((err = WriteTree(datatree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(datatree->GetName(), 0);
      }//if
   } else {
      fDataTreeName = "NA";
      TString str1 = ""; str1 += numcel;
      TString str2 = ""; str2 += fNCells;
      err = fManager->HandleError(errNumCells, str1, str2);
   }//if

// Delete data tree from RAM
cleanup:
   datatree->Delete("");
   datatree = 0;
   delete cell;

   if (quantI) {delete [] quantI; quantI = 0;}

   return err;
}//ReadXDAData

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadFileHeader(ifstream &input, Int_t &numgroups, UInt_t &filepos)
{
   // Read Calvin file header 
   if(kCS) cout << "------XGeneChipHyb::ReadFileHeader------" << endl;

// Read magic number
   unsigned char magic;
   READ_UCHAR(input, magic);
   if (magic != kMagicCalvinNumber) {
      TString str = ""; str += magic;
      return fManager->HandleError(errCELVersion, str);
   }//if

// Read version
   unsigned char version;
   READ_UCHAR(input, version);
   if (version != 1) {
      TString str = ""; str += version;
      return fManager->HandleError(errCELVersion, str);
   }//if

// Read number of data groups
   READ_INT(input, numgroups, kTRUE);
   if (numgroups != 1) {
      cerr << "Error: Number of data groups is not 1!" << endl;
      return errCELVersion;
   }//if

// Read file position of the first data group
   READ_UINT(input, filepos, kTRUE);

   return errNoErr;
}//ReadFileHeader

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadGenericDataHeader(ifstream &input, Bool_t isParent)
{
   // Read Calvin generic data header 
   if(kCS) cout << "------XGeneChipHyb::ReadGenericDataHeader------" << endl;

   Int_t    err  = errNoErr;
   char    *str  = 0;
   wchar_t *wstr = 0;

   // data type identifier.
   READ_STRING(input, str, kTRUE);
   delete[] str; str = 0;

   // unique file identifier
   READ_STRING(input, str, kTRUE);
   delete[] str; str = 0;

   // date and time of file creation
   READ_WSTRING(input, wstr, kTRUE);
   delete[] wstr; wstr = 0;

   // locale
   READ_WSTRING(input, str, kTRUE);
   delete[] str; str = 0;

   // number of name/value/type triplets
   Int_t numtriplets = 0;
   READ_INT(input, numtriplets, kTRUE);

   // name/value/type triplets
   AWSTRING *aname  = 0;
   ASTRING  *avalue = 0;
   AWSTRING *atype  = 0;
   for (Int_t i=0; i<numtriplets; i++) {
      aname  = new AWSTRING;
      avalue = new ASTRING;
      atype  = new AWSTRING;

      READ_WSTRING(input, aname, kTRUE);

      READ_STRING(input, avalue, kTRUE);

      // get chip name
      if (wcscmp(aname->value, L"affymetrix-array-type") == 0) {
         str  = new char[avalue->len + 1];
         wstr = DecodeTEXT(avalue);
         wcstombs(str, wstr, avalue->len + 1);
         if (!isParent) fSchemeName = TString(str);
         delete[] wstr; wstr = 0;
         delete[] str;  str  = 0;
      }//if

      // get chip name from dat-header
      // see: https://www.stat.math.ethz.ch/pipermail/bioconductor/2007-October/019665.html
      if ((wcscmp(aname->value, L"affymetrix-dat-header")         == 0) ||
          (wcscmp(aname->value, L"affymetrix-partial-dat-header") == 0)) {
         str  = new char[avalue->len + 1];
         wstr = DecodeTEXT(avalue);
         wcstombs(str, wstr, avalue->len + 1);

         if (strlen(str) > 0) {
            if ((err = this->CheckChipType(str, fSchemeName)) != errNoErr) {
               return fManager->HandleError(err, fSchemeName, str);
            }//if
         }//if

         delete[] wstr; wstr = 0;
         delete[] str;  str  = 0;
      }//if

      // get number of columns
      if (wcscmp(aname->value, L"affymetrix-cel-cols") == 0) {
         fNCols = DecodeINT(avalue);
     }//if

      // get number of rows
      if (wcscmp(aname->value, L"affymetrix-cel-rows") == 0) {
         fNRows = DecodeINT(avalue);
      }//if

      READ_WSTRING(input, atype, kTRUE);

      delete atype;
      delete avalue;
      delete aname;
   }//for_i

// Read number of parents
   Int_t numparents = 0;
   READ_INT(input, numparents, kTRUE);

// Read Generic Data Header recursively
   for (Int_t i=0; i<numparents; i++) {
      err = ReadGenericDataHeader(input, kTRUE);
      if (err != errNoErr) return err;
   }//for_i

// Check if optional scmname is identical to name of imported scheme
   TString scmname = ((XDataSetting*)fSetting)->GetSchemeName();
   if ((strcmp(scmname.Data(), "") != 0) && 
       (strcmp(scmname.Data(), fSchemeName.Data()) != 0)) {
      return fManager->HandleError(errChipType, scmname, fSchemeName);
   }//if

   return err;
}//ReadGenericDataHeader

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadDataGroup(ifstream &input, UInt_t &filepos,
                    Option_t *option, Int_t split)
{
   // Read Calvin data group. 
   if(kCS) cout << "------XGeneChipHyb::ReadDataGroup------" << endl;

   Int_t    x, y;
   Int_t    err    = 0;
   Double_t min    = DBL_MAX;  //defined in float.h
   Double_t max    = 0;
   Int_t    nummin = 0;
   Int_t    nummax = 0;
   Short_t  maxpix = 0;
   char    *str;

   Int_t     nquant = 7;
   Double_t  q[]    = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};

// Create data tree
   TString exten = Path2Name(option, ".", "");
   fDataTreeName = fTreeName + "." + exten;
   TTree *datatree = new TTree(fDataTreeName, fSchemeName);
   if (datatree == 0) return errCreateTree;

   XGCCell *cell    = new XGCCell();
   Int_t    bufsize = XManager::GetBufSize();
   datatree->Branch("DataBranch", "XGCCell", &cell, bufsize, split);

// Data group
   // file position of next data group
   UInt_t nextpos = 0;
   READ_UINT(input, nextpos, kTRUE);

   // file position of first data element in data set
   READ_UINT(input, filepos, kTRUE);

   // number of datasets
   Int_t numdatasets = 0;
   READ_INT(input, numdatasets, kTRUE);

   // data group name
   READ_WSTRING(input, str, kTRUE);
   delete[] str; str = 0;

// File position of first data set within group
   input.clear(); input.seekg(filepos, ios::beg);
   filepos = nextpos;  //pass parameter filepos for next data group

// Data set: ------------ Intensity ------------
   UInt_t datapos = 0;
   // file position of first data element in data set
   READ_UINT(input, datapos, kTRUE);

   // file position of next data set within the data group
   READ_UINT(input, datapos, kTRUE);

   // data set name
   READ_WSTRING(input, str, kTRUE);
   delete[] str; str = 0;

   // number of name/value/type parameters
   Int_t numtriplets = 0;
   READ_INT(input, numtriplets, kTRUE);

   // number of columns in the data set.
   UInt_t numcols = 0;
   READ_UINT(input, numcols, kTRUE);

   // column names, column value types and column type sizes triplets
   char  coltype;
   Int_t colsize = 0;
   for (UInt_t i=0; i<numcols; i++) {
      READ_WSTRING(input, str, kTRUE);
      delete[] str; str = 0;

      READ_CHAR(input, coltype);
      READ_INT(input, colsize, kTRUE);
   }//for_i

   // number of rows in the data set: fNCells
   UInt_t numrows = 0;
   READ_UINT(input, numrows, kTRUE);
   fNCells = (Int_t)numrows;

// Create local arrays to store data from input
   Float_t  *inten  = 0;
   Float_t  *stdev  = 0;
   Short_t  *numpix = 0;
   Double_t *quantI = 0; 

   // initialize memory for local arrays
   if (!(inten  = new (nothrow) Float_t[fNCells])) {err = errInitMemory; goto cleanup;}
   if (!(stdev  = new (nothrow) Float_t[fNCells])) {err = errInitMemory; goto cleanup;}
   if (!(numpix = new (nothrow) Short_t[fNCells])) {err = errInitMemory; goto cleanup;}
   if (!(quantI = new (nothrow) Double_t[nquant])) {err = errInitMemory; goto cleanup;}

   // read intensities
   for (UInt_t i=0; i<numrows; i++) {
      READ_FLOAT(input, inten[i], kTRUE);
   }//for_i

// Data set: ------------ StdDev ------------
   // jump to data set StdDev
   input.clear(); input.seekg(datapos, ios::beg);

   // file position of first data element in data set
   READ_UINT(input, datapos, kTRUE);

   // file position of next data set within the data group
   READ_UINT(input, datapos, kTRUE);

   // data set name
   READ_WSTRING(input, str, kTRUE);
   delete[] str; str = 0;

   // number of name/value/type parameters
   READ_INT(input, numtriplets, kTRUE);

   // number of columns in the data set.
   READ_UINT(input, numcols, kTRUE);

   // column names, column value types and column type sizes triplets
   for (UInt_t i=0; i<numcols; i++) {
      READ_WSTRING(input, str, kTRUE);
      delete[] str; str = 0;

      READ_CHAR(input, coltype);
      READ_INT(input, colsize, kTRUE);
   }//for_i

   // number of rows in the data set
   READ_UINT(input, numrows, kTRUE);
   if (numrows != (UInt_t)fNCells) {
      fDataTreeName = "NA";
      err = errReadingInput;
      cerr << "Error: number of stdev rows <" << numrows 
           << "> is not equal to to number of cells <" << fNCells << ">"
           << endl;
      goto cleanup;
   }//if

   // read stdev
   for (UInt_t i=0; i<numrows; i++) {
      READ_FLOAT(input, stdev[i], kTRUE);
   }//for_i

// Data set: ------------ Pixel ------------
   // jump to data set Pixel
   input.clear(); input.seekg(datapos, ios::beg);

   // file position of first data element in data set
   READ_UINT(input, datapos, kTRUE);

   // file position of next data set within the data group
   READ_UINT(input, datapos, kTRUE);

   // data set name
   READ_WSTRING(input, str, kTRUE);
   delete[] str; str = 0;

   // number of name/value/type parameters
   READ_INT(input, numtriplets, kTRUE);

   // number of columns in the data set.
   READ_UINT(input, numcols, kTRUE);

   // column names, column value types and column type sizes triplets
   for (UInt_t i=0; i<numcols; i++) {
      READ_WSTRING(input, str, kTRUE);
      delete[] str; str = 0;

      READ_CHAR(input, coltype);
      READ_INT(input, colsize, kTRUE);
   }//for_i

   // number of rows in the data set
   READ_UINT(input, numrows, kTRUE);
   if (numrows != (UInt_t)fNCells) {
      fDataTreeName = "NA";
      err = errReadingInput;
      cerr << "Error: number of stdev rows <" << numrows 
           << "> is not equal to to number of cells <" << fNCells << ">"
           << endl;
      goto cleanup;
   }//if

   // numpix
   for (UInt_t i=0; i<numrows; i++) {
      READ_SHORT(input, numpix[i], kTRUE);
   }//for_i

// Data set: ------------ Outlier ------------
   // not used

// Data set: ------------ Mask ---------------
   // not used

// Fill data tree
   for (UInt_t i=0; i<numrows; i++) {
      x = Index2X(i);
      y = Index2Y(i);

      // number of cells with minimal intensity
      if      (inten[i] <  min) {min = inten[i]; nummin = 1;}
      else if (inten[i] == min) {nummin++;}

      // number of cells with maximal intensity
      if      (inten[i] >  max) {max = inten[i]; nummax = 1;}
      else if (inten[i] == max) {nummax++;}

      // fill data tree
      cell->SetX(x);
      cell->SetY(y);
      cell->SetIntensity((Double32_t)inten[i]);
      cell->SetStdev((Double32_t)stdev[i]);
      cell->SetNumPixels(numpix[i]);
      datatree->Fill();
   }//for_i

   fMinInten   = min;
   fNMinInten  = nummin;
   fMaxInten   = max;
   fNMaxInten  = nummax;
   fMaxNPixels = maxpix;

   if (XManager::fgVerbose) {
      cout << "   hybridization statistics: " << endl;
      cout << "      " << nummin << " cells with minimal intensity " << min << endl;
      cout << "      " << nummax << " cells with maximal intensity " << max << endl;
   }//if

// Check for problems with intensities
   if (max <= min) {
      cout << "Warning: maximal intensity equal or less than minimal intensity!" << endl;
//      cout << "Thus CEL-file will not be imported as <" << fDataTreeName << ">!" << endl;
//      goto cleanup;
   }//if

// Write data tree to file
   // quantiles for residual trees
   err = this->DataQuantiles(datatree, cell, nquant, q, quantI);
   if (err != errNoErr) goto cleanup;

   // add tree info to tree
//   AddDataTreeInfo(datatree, fDataTreeName);
   AddDataTreeInfo(datatree, datatree->GetName(), "generic", nquant, q, quantI);

   // write data tree to file 
   if ((err = WriteTree(datatree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(datatree->GetName(), 0);
   }//if

cleanup:
   if (quantI) {delete [] quantI; quantI = 0;}
   if (numpix) {delete [] numpix; numpix = 0;}
   if (stdev)  {delete [] stdev;  stdev  = 0;}
   if (inten)  {delete [] inten;  inten  = 0;}

   // delete data tree from RAM
   datatree->Delete("");
   datatree = 0;
   delete cell;

   return err;
}//ReadDataGroup

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadCalvinGenericFile(ifstream &input, Option_t *option, Int_t split)
{
   // Read Calvin genreric data file 
   if(kCS) cout << "------XGeneChipHyb::ReadCalvinGenericFile------" << endl;

   Int_t  err       = errNoErr;
   Int_t  numgroups = 0;   //number of data groups
   UInt_t filepos   = 0;   //file position of (first) data group

// File Header
   err = ReadFileHeader(input, numgroups, filepos);
   if (err != errNoErr) return err;

// Generic Data Header 
   err = ReadGenericDataHeader(input, kFALSE);
   if (err != errNoErr) return err;

// Data Groups
   for (Int_t i=0; i<numgroups; i++) {
      err = ReadDataGroup(input, filepos, option, split);
      if (err != errNoErr) return err;
   }//for_i

   return err;
}//ReadCalvinGenericFile

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadXMLHeader(ifstream &input, const char * /*sep*/, char delim)
{
   // Read header from input in XML-format. 
   if(kCS) cout << "------XGeneChipHyb::ReadXMLHeader------" << endl;

   char  nextline[kBufSize];
   Int_t err = 0;
   Int_t dif;
   char *pos;
   TString name;

// Check for first line "<?xml"
   input.getline(nextline, kBufSize, delim);
   if ( strncmp("<?xml", nextline, 5) != 0) {
      cerr << "Error: Header line for XML-file is missing!" << endl;
      return errHeaderLine;
   }//if

// Check for line "<MAGE-ML"
   input.getline(nextline, kBufSize, delim);
   input.getline(nextline, kBufSize, delim);
   if ( strncmp("<MAGE-ML", nextline, 8) != 0) {
      cerr << "Error: Header line  for MAGE-ML file is missing!" << endl;
      return errHeaderLine;
   }//if

// Check if <Protocol ... title="CEL Analysis">
   while (strstr(nextline,"<Protocol_assnlist") == 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   input.getline(nextline, kBufSize, delim);
   if ( strstr("CEL Analysis", nextline) != 0) {
      cerr << "Error: File does not describe a CEL Analysis!" << endl;
      return errGeneral;
   }//if

// Get file name from <MeasuredBioAssay ... name="xxx">
   while (strstr(nextline,"<MeasuredBioAssay") == 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   pos = strstr(nextline,"name=");
   if (pos == 0) {
      cerr << "Error: MeasuredBioAssay name not found!" << endl;
      return errGeneral;
   }//if
   dif = pos + 5 - &nextline[0];
   name = &nextline[dif];
   name = RemoveEnds(name);

// Get chip name from <NameValueType ... value="xxx">
   while (strstr(nextline,"<NameValueType") == 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   name = &nextline[dif];
   fSchemeName = RemoveEnds(name);

// Check if optional scmname is identical to name of imported scheme
   TString scmname = ((XDataSetting*)fSetting)->GetSchemeName();
   if ((strcmp(scmname.Data(), "") != 0) && 
       (strcmp(scmname.Data(), fSchemeName.Data()) != 0)) {
      return fManager->HandleError(errChipType, scmname, fSchemeName);
   }//if

// Get values from <SummaryStatistics_assnlist ... value="xxx">
   while (strstr(nextline,"<SummaryStatistics_assnlist") == 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

// Get Number Cells Masked
   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   Int_t fNMasked = atoi(strtok(&nextline[dif], "\""));
//??   fNMasked = atoi(strtok(&nextline[dif], "\""));

// Get Number Outlier Cells
   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   Int_t fNOutlier = atoi(strtok(&nextline[dif], "\""));
//??   fNOutlier = atoi(strtok(&nextline[dif], "\""));

// Get Number Cells Modified
   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   Int_t fNModified = atoi(strtok(&nextline[dif], "\""));
//??   fNModified = atoi(strtok(&nextline[dif], "\""));

// Get Rows
   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   fNRows = atoi(strtok(&nextline[dif], "\""));

// Get Cols
   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   fNCols = atoi(strtok(&nextline[dif], "\""));

// Get Number of Cells
   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"value=");
   if (pos == 0) {
      return fManager->HandleError(errNameValue);
   }//if
   dif = pos + 6 - &nextline[0];
   fNCells = atoi(strtok(&nextline[dif], "\""));

   return err;
}//ReadXMLHeader

//______________________________________________________________________________
Int_t XGeneChipHyb::ReadXMLData(ifstream &input, Option_t *option, const char * /*sep*/,
                    char delim, Int_t split)
{
   // Read data from input. 
   if(kCS) cout << "------XGeneChipHyb::ReadXMLData------" << endl;

   char    nextline[kBufSize];
   char   *pos;
   Int_t   dif;
   TString name;

// Get data filename from <DataExternal ... value="xxx">
   input.getline(nextline, kBufSize, delim);
   while (strstr(nextline,"<DataExternal_assn>") == 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   input.getline(nextline, kBufSize, delim);
   pos = strstr(nextline,"filenameURI=");
   if (pos == 0) {
      cerr << "Error: DataExternal filenameURI not found!" << endl;
      return errGeneral;
   }//if
   dif = pos + 12 - &nextline[0];
   name = &nextline[dif];
   name = RemoveEnds(name);

// Prepend directory containing data file
   name = Name2Path(fInfile, sSEP) + TString(dSEP) + name;

// Open externaldata.txt file containing data
   ifstream data(name.Data(), ios::in);
   if (!data) {
      return fManager->HandleError(errReadingInput, name);
   }//if

   Int_t    i, x, y;
   Short_t  numpix;
   Double_t inten, stdev;
   Bool_t   outlier, masked;
   TString  otl;
   TString  msk;
   Int_t    err    = 0;
   Int_t    numcel = 0;
   Double_t min    = DBL_MAX;
   Double_t max    = 0;
   Int_t    nummin = 0;
   Int_t    nummax = 0;
   Short_t  maxpix = 0;

   Int_t     nquant = 7;
   Double_t  q[]    = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0};
   Double_t *quantI = 0; 
   if (!(quantI = new (nothrow) Double_t[nquant])) return errInitMemory;

// Create data tree
   TString exten = Path2Name(option,".","");
   fDataTreeName = fTreeName + "." + exten;
   TTree *datatree = new TTree(fDataTreeName, fSchemeName);
   if (datatree == 0) return errCreateTree;

   XGCCell *cell    = new XGCCell();
   XCellOM *omsk    = new XCellOM();
   Int_t    bufsize = XManager::GetBufSize();
   datatree->Branch("DataBranch", "XGCCell", &cell, bufsize, split);
   datatree->Branch("OMskBranch", "XCellOM", &omsk, bufsize, split);

// Read data and store in data tree
   for (i=0; i<fNCells; i++) {
      data >> x >> y >> inten >> stdev >> numpix >> otl >> msk;
      if (!data.good()) {err = errPrematureEOF; break;}

      // number of cells with minimal intensity
      if (inten < min) {
         min = inten;
         nummin = 1;
      } else if (inten == min) {
         nummin++;
      }//if

      // number of cells with maximal intensity
      if (inten > max) {
         max = inten;
         nummax = 1;
      } else if (inten == max) {
         nummax++;
      }//if

      // maximal pixel number
      if (numpix > maxpix) {
         maxpix = numpix;
      }//if

      // outlier and mask
      if (strcmp(otl.Data(),"true") == 0) outlier = kTRUE;
      else                                outlier = kFALSE;
      if (strcmp(msk.Data(),"true") == 0) masked  = kTRUE;
      else                                masked  = kFALSE;

      // fill data tree
      cell->SetX(x);
      cell->SetY(y);
      cell->SetIntensity(inten);
      cell->SetStdev(stdev);
      cell->SetNumPixels(numpix);
      omsk->SetOutlier(outlier);
      omsk->SetMasked(masked);
      datatree->Fill(); 

      numcel++;
      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << fNCells << "> records imported." << endl;
   }//if

   fMinInten   = min;
   fNMinInten  = nummin;
   fMaxInten   = max;
   fNMaxInten  = nummax;
   fMaxNPixels = maxpix;

   if (XManager::fgVerbose) {
      cout << "   hybridization statistics: " << endl;
      cout << "      " << nummin << " cells with minimal intensity " << min << endl;
      cout << "      " << nummax << " cells with maximal intensity " << max << endl;
   }//if

// Check for problems with intensities
   if (max <= min) {
      cout << "Warning: maximal intensity equal or less than minimal intensity!" << endl;
//      cout << "Thus CEL-file will not be imported as <" << fDataTreeName << ">!" << endl;
//      goto cleanup;
   }//if

// Write data tree to file if all data are read
   if (numcel == fNCells) {
      // quantiles for residual trees
      err = this->DataQuantiles(datatree, cell, nquant, q, quantI);
      if (err != errNoErr) goto cleanup;

      // add tree info to tree
//      AddDataTreeInfo(datatree, fDataTreeName);
      AddDataTreeInfo(datatree, datatree->GetName(), "xml", nquant, q, quantI);

   // Write data tree to file 
      if ((err = WriteTree(datatree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(datatree->GetName(), 0);
      }//if
   } else {
      fDataTreeName = "NA";
      TString str1 = ""; str1 += numcel;
      TString str2 = ""; str2 += fNCells;
      err = fManager->HandleError(errNumCells, str1, str2);
   }//if

// Cleanup
cleanup:
   // delete data tree from RAM
   datatree->Delete("");
   datatree = 0;
   delete cell;
   delete omsk;

   if (quantI) {delete [] quantI; quantI = 0;}

   data.close();
   return err;
}//ReadXMLData

//______________________________________________________________________________
Int_t XGeneChipHyb::DataQuantiles(TTree *tree, XGCCell *cell,
                    Int_t nquant, Double_t *q, Double_t *quant)
{
   // Get nquant quantiles for intensities from data tree
   if(kCS) cout << "------XGeneChipHyb::DataQuantiles------" << endl;

   Int_t err      = errNoErr;
   Int_t nentries = (Int_t)(tree->GetEntries());

   tree->SetBranchAddress("DataBranch", &cell);

// Init arrays
   Double_t *inten = 0; 
   Int_t    *index = 0;
   if (!(inten = new (nothrow) Double_t[nentries])) {err = errInitMemory; goto cleanup;}
   if (!(index = new (nothrow) Int_t[nentries]))    {err = errInitMemory; goto cleanup;}

// Fill arrays
   for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);
      inten[i] = cell->GetIntensity();
   }//for_i

// Fill quantiles
   quant = TStat::Quantiles(nentries, inten, index, nquant, q, quant);

// Cleanup
cleanup:
   tree->DropBaskets();
   tree->ResetBranchAddress(tree->GetBranch("DataBranch"));

   if (index) {delete [] index; index = 0;}
   if (inten) {delete [] inten; inten = 0;}

   return err;
}//DataQuantiles


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipMetrics                                                     //
//                                                                      //
// Class for processed GeneChip CHP data exported as Metrics.txt        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGeneChipMetrics::XGeneChipMetrics()
                 :XHybridization()
{
   // Default GeneChipMetrics constructor
   if(kCS) cout << "---XGeneChipMetrics::XGeneChipMetrics(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGeneChipMetrics::XGeneChipMetrics(const char *name, const char *title)
                 :XHybridization(name, title)
{
   // Normal GeneChipMetrics constructor
   if(kCS) cout << "---XGeneChipMetrics::XGeneChipMetrics------" << endl;

}//Constructor

//______________________________________________________________________________
XGeneChipMetrics::~XGeneChipMetrics()
{
   // GeneChipMetrics destructor
   if(kCS) cout << "---XGeneChipMetrics::~XGeneChipMetrics------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGeneChipMetrics::ExportTreeType(const char *exten, Int_t n, TString *names, 
                       const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output
   if(kCS) cout << "------XGeneChipMetrics::ExportTreeType------" << endl;

   if (HasExtension(exten, kPivotExpr)) {
      return this->ExportExprTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kPivotCall)) {
      return this->ExportCallTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//ExportTreeType

//______________________________________________________________________________
Int_t XGeneChipMetrics::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                        const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XGeneChipMetrics::ExportTreeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   exten = 0; n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
void XGeneChipMetrics::PrintInfo()
{
   // Print GeneChip Metrics information
   if(kCS) cout << "------XGeneChipMetrics::PrintInfo------" << endl;

   if (fgPrintHeader) {
      cout << "==============================================================================" << endl;
      cout << setw(14) << "ChipName" << setw(12) << "Title"
           << endl;
      cout << "==============================================================================" << endl;
      fgPrintHeader = kFALSE;
   }//if

   cout << setw(14) << this->GetName() << setw(12) << this->GetTitle()

        << endl;

   cout << "------------------------------------------------------------------------------" << endl;
}//PrintInfo

//______________________________________________________________________________
Int_t XGeneChipMetrics::ExportExprTrees(Int_t n, TString *names, const char *varlist,
                        ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGeneChipMetrics::ExportExprTrees------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of GeneChip Metics data is not yet implemented" << endl;
   return 1;
//   return errNoErr;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XGeneChipMetrics::ExportCallTrees(Int_t n, TString *names, const char *varlist,
                        ofstream &output, const char *sep)
{
   // Export mask tree to file output
   if(kCS) cout << "------XGeneChipMetrics::ExportCallTrees------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of GeneChip Metics data is not yet implemented" << endl;
   return 1;
//   return errNoErr;
}//ExportCallTrees

//______________________________________________________________________________
Int_t XGeneChipMetrics::ReadHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input. 
   if(kCS) cout << "------XGeneChipMetrics::ReadHeader------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); sep = 0; delim = 0;

   return 0;
}//ReadHeader

//______________________________________________________________________________
Int_t XGeneChipMetrics::ReadData(ifstream &input, Option_t *option, const char *sep,
                        char delim, Int_t split)
{
   // Read data from input. 
   if(kCS) cout << "------XGeneChipMetrics::ReadData------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of GeneChip Metics data is not yet implemented" << endl;
   return 1;
//   return 0;
}//ReadData

//______________________________________________________________________________
Int_t XGeneChipMetrics::ReadXMLHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input in XML-format. 
   if(kCS) cout << "------XGeneChipMetrics::ReadXMLHeader------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); sep = 0; delim = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ReadXMLHeader

//______________________________________________________________________________
Int_t XGeneChipMetrics::ReadXMLData(ifstream &input, Option_t *option,
                        const char *sep, char delim, Int_t split)
{
   // Read data from input. 
   if(kCS) cout << "------XGeneChipMetrics::ReadXMLData------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return 1;
}//ReadXMLData


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipPivot                                                       //
//                                                                      //
// Class for processed GeneChip CHP data exported as PivotData.txt      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGeneChipPivot::XGeneChipPivot()
               :XHybridization()
{
   // Default GeneChipPivot constructor
   if(kCS) cout << "---XGeneChipPivot::XGeneChipPivot(default)------" << endl;

   fNTrees   = 0;
   fTreeName = 0;
}//Constructor

//______________________________________________________________________________
XGeneChipPivot::XGeneChipPivot(const char *name, const char *title)
               :XHybridization(name, title)
{
   // Normal GeneChipPivot constructor
   if(kCS) cout << "---XGeneChipPivot::XGeneChipPivot------" << endl;

   fNTrees   = 0;
   fTreeName = 0;
}//Constructor

//______________________________________________________________________________
XGeneChipPivot::~XGeneChipPivot()
{
   // GeneChipPivot destructor
   if(kCS) cout << "---XGeneChipPivot::~XGeneChipPivot------" << endl;

   if (fTreeName) {delete [] fTreeName; fTreeName = 0;}
   fNTrees = 0;
}//Destructor

//______________________________________________________________________________
Int_t XGeneChipPivot::ExportTreeType(const char *exten, Int_t n, TString *names, 
                      const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output
   if(kCS) cout << "------XGeneChipPivot::ExportTreeType------" << endl;

// Set scheme file to be able to access scheme data for exporting
   if (fSetting) {
      fSchemeFile = ((XDataSetting*)fSetting)->GetSchemeFile();
   }//if

   if (HasExtension(exten, kPivotExpr)) {
      return this->ExportExprTrees(n, names, varlist, output, sep);
   } else if (HasExtension(exten, kPivotCall)) {
      return this->ExportCallTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//ExportTreeType

//______________________________________________________________________________
Int_t XGeneChipPivot::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                      const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XGeneChipPivot::ExportTreeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   exten = 0; n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTreeXML

//______________________________________________________________________________
void XGeneChipPivot::PrintInfo()
{
   // Print GeneChip PivotData information
   if(kCS) cout << "------XGeneChipPivot::PrintInfo------" << endl;

   if (fgPrintHeader) {
      cout << "==============================================================================" << endl;
      cout << setw(14) << "ChipName" << setw(12) << "Title"
           << endl;
      cout << "==============================================================================" << endl;
      fgPrintHeader = kFALSE;
   }//if

   cout << setw(14) << this->GetName() << setw(12) << this->GetTitle()

        << endl;

   cout << "------------------------------------------------------------------------------" << endl;
}//PrintInfo

//______________________________________________________________________________
Int_t XGeneChipPivot::ExportExprTrees(Int_t n, TString *names, const char *varlist,
                      ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGeneChipPivot::ExportExprTrees------" << endl;

// Decompose varlist
   Bool_t hasUnit   = kFALSE;
   Bool_t hasName   = kFALSE;
   Bool_t hasSymbol = kFALSE;
   Bool_t hasCyto   = kFALSE;
   Bool_t hasAnnot  = kFALSE;
   Bool_t hasLevel  = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasUnit   = kTRUE;
      hasName   = kTRUE;
      hasSymbol = kTRUE;
      hasCyto   = kTRUE;
      hasLevel  = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name, varlist), ":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit   = kTRUE;}
         if (strcmp(name,"fName")     == 0) {hasName   = kTRUE;}
         if (strcmp(name,"fSymbol")   == 0) {hasSymbol = kTRUE;}
         if (strcmp(name,"fCytoBand") == 0) {hasCyto   = kTRUE;}
         if (strcmp(name,"fLevel")    == 0) {hasLevel  = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if
   hasAnnot = (hasName || hasSymbol || hasCyto);

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

// Get scheme name
   TString scmname = tree[0]->GetTitle();
   if (strcmp(scmname.Data(), "")  == 0) {
      cerr << "Error: No scheme name is present. Please report error." << endl;
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

// If scheme file is present get unit tree
   Int_t        numgenes = 0;
   Int_t        numctrls = 0;
   XFolder     *schemes  = 0;
   XGCUnit     *unit     = 0;
   XAnnotation *annot    = 0;
   TTree       *unittree = 0; 
   TTree       *anntree  = 0; 
   TDirectory  *savedir  = gDirectory;

   if ((hasUnit || hasAnnot) && fSchemeFile) {
   // Get chip from scheme file
      fSchemeFile->cd();
      schemes = (XFolder*)(fSchemeFile->Get(kContent));
      if (!schemes) {
         return fManager->HandleError(errMissingContent, "Scheme", kContent);
      }//if

      XDNAChip *chip = 0;
      chip = (XDNAChip*)schemes->FindObject(scmname, kTRUE);
      if (chip) {
         numgenes = chip->GetNumGenes();
         numctrls = chip->GetNumControls();

         if (!fSchemeFile->cd(scmname)) return errGetDir;

      // Get unit tree for scheme
         if (hasUnit) {
            unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
            if (unittree == 0) return errGetTree;
            unittree->SetBranchAddress("IdxBranch", &unit);
         }//if

      // Get annotation tree for scheme
         if (hasAnnot) {
            anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
            if (anntree) {
               anntree->SetBranchAddress("AnnBranch", &annot);
            } else {
               cout << "Warning: Missing annotation, gene info not exported."
                    << endl;
               hasAnnot = kFALSE;
            }//if
         }//if
      } else {
         cerr << "Error: Could not find scheme <" << scmname.Data() << ">." << endl;
         hasUnit  = kFALSE;
         hasAnnot = kFALSE;
      }//if
   } else if (hasUnit || hasAnnot) {
      cout << "Warning: Missing scheme file, unit names not exported." << endl;
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

   savedir->cd();

// Check if tree entries is equal to number of genes in unittree
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   if ((hasUnit || hasAnnot) && (entries != numgenes)) {
      if (XManager::fgVerbose) {
         cout << "Warning: Number of tree entries <" << entries  
              << "> is not equal to number of units <" << numgenes << ">" << endl;
         cout << "         Unit data will not be exported." << endl;
      }//if
      hasUnit  = kFALSE;
      hasAnnot = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit)      output << sep << "UNIT_NAME";
   if (hasAnnot) {
      if (hasName)   output << sep << "GENE_NAME";
      if (hasSymbol) output << sep << "GENE_SYMBOL";
      if (hasCyto)   output << sep << "CYTOBAND";
   }//if
   if (n == 1) {
      if (hasLevel)  output << sep << "LEVEL";
   } else {
      for (Int_t i=0; i<n; i++) {
         if (hasLevel) output << sep << (names[i] + "_LEVEL").Data();
      }//for_i
   }//if
   output << endl;

// Loop over tree entries and tree branches
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         if (k == 0) {
            output << expr[k]->GetUnitID();

            if (hasUnit) {
               unittree->GetEntry(i + numctrls); //unit names for genes only
               output << sep << unit->GetUnitName();
            }//if

            if (hasAnnot) {
               anntree->GetEntry(i + numctrls);
               if (hasName)   output << sep << annot->GetName();
               if (hasSymbol) output << sep << annot->GetSymbol();
               if (hasCyto)   output << sep << annot->GetCytoBand();
            }//if
         }//if

         if (hasLevel) output << sep << expr[k]->GetLevel();
      }//for_j
      output << endl;
   }//for_i

//Cleanup
   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (unittree) {unittree->Delete(""); unittree = 0;}
   SafeDelete(schemes);

   delete [] expr;
   delete [] tree;

   return errNoErr;
}//ExportExprTrees

//______________________________________________________________________________
Int_t XGeneChipPivot::ExportCallTrees(Int_t n, TString *names, const char *varlist,
                      ofstream &output, const char *sep)
{
   // Export mask tree to file output
   if(kCS) cout << "------XGeneChipPivot::ExportCallTrees------" << endl;

// Decompose varlist
   Bool_t hasUnit = kFALSE;
   Bool_t hasCall = kFALSE;
   Bool_t hasPVal = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasUnit = kTRUE;
      hasCall = kTRUE;
      hasPVal = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit = kTRUE;}
         if (strcmp(name,"fCall")     == 0) {hasCall = kTRUE;}
         if (strcmp(name,"fPValue")   == 0) {hasPVal = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get trees
   TTree  **tree = new TTree*[n];
   XPCall **call = new XPCall*[n];
   if (fTrees->GetSize() == 0) {
   // Get trees from names
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)gDirectory->Get((names[k]).Data());
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   } else {
   // Get trees from list fTrees
      for (Int_t k=0; k<n; k++) {
         call[k] = 0;
         tree[k] = (TTree*)fTrees->At(k);
         if (!tree[k]) return errGetTree;

         tree[k]->SetBranchAddress("CallBranch", &call[k]);
      }//for_k
   }//if

// Get scheme name
   TString scmname = tree[0]->GetTitle();
   if (strcmp(scmname.Data(), "")  == 0) {
      cerr << "Error: No scheme name is present. Please report error." << endl;
      hasUnit = kFALSE;
   }//if

// If scheme file is present get unit tree
   Int_t       numgenes  = 0;
   Int_t       numctrls  = 0;
   XFolder     *schemes  = 0;
   XGCUnit     *unit     = 0;
   TTree       *unittree = 0; 
   TDirectory  *savedir  = gDirectory;

   if (hasUnit && fSchemeFile) {
   // Get chip from scheme file
      schemes = (XFolder*)(fSchemeFile->Get(kContent));
      if (!schemes) {
         return fManager->HandleError(errMissingContent, "Scheme", kContent);
      }//if

      XDNAChip *chip = 0;
      chip = (XDNAChip*)schemes->FindObject(scmname, kTRUE);
      if (chip) {
         numgenes = chip->GetNumGenes();
         numctrls = chip->GetNumControls();

         if (!fSchemeFile->cd(scmname)) return errGetDir;

      // Get unit tree for scheme
         unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
         if (unittree == 0) return errGetTree;
         unittree->SetBranchAddress("IdxBranch", &unit);
      } else {
         cerr << "Error: Could not find scheme <" << scmname.Data() << ">." << endl;
         hasUnit = kFALSE;
      }//if
   } else if (hasUnit) {
      cerr << "Warning: Missing scheme file, unit names not exported." << endl;
      hasUnit = kFALSE;
   }//if

   savedir->cd();

// Check if tree entries is equal to number of genes in unittree
   Int_t entries = (Int_t)(tree[0]->GetEntries());
   if (hasUnit && (entries != numgenes)) {
      if (XManager::fgVerbose) {
         cout << "Warning: Number of tree entries <" << entries  
              << "> is not equal to number of units <" << numgenes << ">" << endl;
      }//if
      hasUnit = kFALSE;
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit) output << sep << "UNIT_NAME";
   if (n == 1) {
      if (hasCall) output << sep << "CALL";
      if (hasPVal) output << sep << "PVALUE";
   } else {
      for (Int_t i=0; i<n; i++) {
         if (hasCall) output << sep << (names[i] + "_CALL").Data();
         if (hasPVal) output << sep << (names[i] + "_PVALUE").Data();
      }//for_i
   }//if
   output << endl;

// Loop over tree entries and tree branches
   for (Int_t i=0; i<entries; i++) {
      for (Int_t k=0; k<n; k++) {
         tree[k]->GetEntry(i);

         if (k == 0) output << call[k]->GetUnitID();
         if (k == 0 && hasUnit) {
            unittree->GetEntry(i + numctrls); //unit names for genes only
            output << sep << unit->GetUnitName();
         }//if

         if (hasCall) {
            Short_t cl = call[k]->GetCall();
            const char *ch[1];
            ch[0] = "NA";
            if      (cl == 2) ch[0] = "P";
            else if (cl == 0) ch[0] = "A";
            else if (cl == 1) ch[0] = "M";
            output << sep << ch[0];
         }//if

         if (hasPVal) output << sep << call[k]->GetPValue();
      }//for_j
      output << endl;
   }//for_i

//Cleanup
   // remove trees from RAM
   if (unittree) {unittree->Delete(""); unittree = 0;}
   SafeDelete(schemes);

   delete [] call;
   delete [] tree;

   return errNoErr;
}//ExportCallTrees

//______________________________________________________________________________
Int_t XGeneChipPivot::ReadHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input. 
   if(kCS) cout << "------XGeneChipPivot::ReadHeader------" << endl;

   char nextline[kLargeBuf];

   input.getline(nextline, kLargeBuf, delim);
   if (!input.good()) return errHeaderLine;

   TString *vars  = ((XDataSetting*)fSetting)->GetVariables();
   TString  str   = "";
   Int_t numtrees = 0;
   char *tmpname  = new char[strlen(nextline) + 1];
   char *delname  = tmpname;
   tmpname = strtok(strcpy(tmpname,nextline), sep);
   while(tmpname) {
      str = TString(tmpname);
      if (str.Contains(vars[0], TString::kIgnoreCase)) {numtrees++;}
      tmpname = strtok(NULL, sep);
      if (tmpname == 0) break;
   }//while
   if (numtrees == 0) {delete [] delname; return errHeaderLine;}
   fNTrees = numtrees;

   fTreeName = new (nothrow) TString[numtrees];
   if (!fTreeName) {delete [] delname; return errInitMemory;}
   Int_t idx = 0;
   tmpname = delname;  //reset tmpname
   tmpname = strtok(strcpy(tmpname, nextline), sep);
   while(tmpname) {
      str = TString(tmpname);
      if (str.Contains(vars[0], TString::kIgnoreCase)) {
         str = RemoveSubString(tmpname, (vars[0]).Data(), kFALSE);
         str = SubString(str.Data(), '(', ')', kTRUE); //ev not necessary?
         str = ReplaceNonAlpha(str.Data(), "_");
         fTreeName[idx] = RemoveEnds(str.Data());
         idx++;
      }//if
      tmpname = strtok(NULL, sep);
      if (tmpname == 0) break;
   }//while

// Cleanup
   delete [] delname;

   if (idx != numtrees) return errHeaderLine;
   return errNoErr;
}//ReadHeader

//______________________________________________________________________________
Int_t XGeneChipPivot::ReadData(ifstream &input, Option_t * /*option*/, const char *sep,
                      char delim, Int_t split)
{
   // Read data from input.
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCS) cout << "------XGeneChipPivot::ReadData------" << endl;

   Int_t err = errNoErr;
   Int_t idx = 0;
   Int_t id, k;
   Int_t numgenes = 0;
   Int_t numctrls = 0;
   Int_t numunits = 0;
   Int_t bufsize  = kBrchBuf/fNTrees;  //?? adapt bufsize to #trees in RAM
   Int_t nvar     = 0;

   TString *vars = 0;
   TString *typs = 0;
   TString  xten;
   TString  str;

//ccc char memory problem??
   char *ch;
   char *buf = new char[128];
   char  strg[128];
   char  nextline[kLargeBuf];
   const char *unitname;

   TString     scmname  = "";
   XFolder    *schemes  = 0;
   XDNAChip   *chip     = 0;
   XGCUnit    *unit     = 0;
   TTree      *unittree = 0; 
   TTree      *tmptree  = 0;
   THashTable *htable   = 0;
   XIdxString *idxstr   = 0;

   TDirectory *savedir  = gDirectory;

   Bool_t   hasExpr = kFALSE;
   Bool_t   hasPVal = kFALSE;
   Bool_t   hasCall = kFALSE;
   Double_t *expr   = new Double_t[fNTrees];
   Double_t *pval   = new Double_t[fNTrees];
   Short_t  *call   = new Short_t[fNTrees];
   TBranch *exprBr  = 0;
   TBranch *pvalBr  = 0;
   TBranch *callBr  = 0;
   TTree  *exprtree = 0;
   TTree  *calltree = 0;
   TList  *exprlist = new TList();
   TList  *calllist = new TList();
   XExpression *xpr = new XExpression();
   XPCall      *pcl = new XPCall();

// Get scheme name
   scmname = ((XDataSetting*)fSetting)->GetSchemeName();
   if (strcmp(scmname.Data(), "")  == 0) {
      cerr << "Error: No scheme name given." << endl;
      err = errAbort; goto cleanup;
   }//if

// Test if scheme file is present
   if (fSchemeFile == 0) {
      cerr << "Error: Missing scheme file." << endl;
      err = errAbort; goto cleanup;
   }//if

// Get chip parameters from scheme file
   fSchemeFile->cd();
   schemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!schemes) {
      err = fManager->HandleError(errMissingContent, "Scheme", kContent);
      goto cleanup;
   }//if

   chip = (XDNAChip*)schemes->FindObject(scmname, kTRUE);
   if (!chip) {
      cerr << "Error: Could not find scheme <" << scmname << ">." << endl;
      err = errAbort; goto cleanup;
   }//if
   numgenes = chip->GetNumGenes();
   numctrls = chip->GetNumControls();
   numunits = chip->GetNumUnits();

// Get unit tree for scheme
   if (!fSchemeFile->cd(scmname)) {err = errGetDir; goto cleanup;}

   unittree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (unittree == 0) {err = errGetTree; goto cleanup;}
   unittree->SetBranchAddress("IdxBranch", &unit);

// Create temporary tree
   savedir->cd();
//////////////
//TO DO: save tmptree in tmp.root file!!!!!
//////////////
   tmptree = new TTree("tmptree", "temporary tree");
   if (tmptree == 0) {err = errCreateTree; goto cleanup;}

   // init values to NA (= -1)
   for (Int_t i=0; i<fNTrees; i++) {
      expr[i] = -1;
      pval[i] = -1;
      call[i] = -1;
   }//for_i

   // get variable names and types
   nvar = ((XDataSetting*)fSetting)->GetNumVariables();
   vars = ((XDataSetting*)fSetting)->GetVariables();
   typs = ((XDataSetting*)fSetting)->GetVarTypes();
   xten = ((XDataSetting*)fSetting)->GetDataType();

   xten.ToLower();
   if (xten.Contains("mas4"))      xten = kPivotExpr[5]; // = "adf";
   else if (xten.Contains("mas5")) xten = kPivotExpr[6]; // = "tbw";
   else if (xten.Contains("rma"))  xten = kPivotExpr[7]; // = "mdp";

   // set tree branches
   tmptree->Branch("strg", (void*)strg, "strg[128]/C");
   for (Int_t i=0; i<nvar; i++) {
      ch = (char*)(typs[i]).Data();
      str = vars[i]; str.ToUpper();

      if (StringInList(str.Data(), kSynExpr, kNSynExpr, kFALSE)) {
         sprintf(buf, "expr[%i]/D", fNTrees);
         tmptree->Branch("exprBr", expr, buf);
         hasExpr = kTRUE;
      }//if

      if (StringInList(str.Data(), kSynCall, kNSynCall, kFALSE)) {
         sprintf(buf, "call[%i]/I", fNTrees);
         tmptree->Branch("callBr", call, buf);
         hasCall = kTRUE;
      }//if

      if (StringInList(str.Data(), kSynPVal, kNSynPVal, kFALSE)) {
         sprintf(buf, "pval[%i]/D", fNTrees);
         tmptree->Branch("pvalBr", pval, buf);
         hasPVal = kTRUE;
      }//if
   }//for_i

// Create hash table to store unit names from input table
   if (!(htable = new THashTable(2*numgenes))) {err = errInitMemory; goto cleanup;}

// Fill temporary tree
   while (input.good()) {
      input.getline(nextline, kLargeBuf, delim);
      if (input.fail() || (idx == numgenes)) break;

      // read unitname and store in hash table
      strcpy(strg,strtok(&nextline[0], sep));
      //get ridd of potential exclamation marks from e.g. "1234_at"
      str = RemoveEnds(strg);   //ev? RemoveEnds(strg, "\"")
      idxstr = new XIdxString(idx, str);
      htable->Add(idxstr);

      // read data columns
      for (Int_t i=0; i<fNTrees; i++) {
         if (hasExpr) {expr[i] = atof(strtok(NULL, sep));}

//to do: recognize variable order!!
         if (hasCall) {
            ch = strtok(NULL, sep);
            if      (ch[0] == 'P') call[i] =  2;
            else if (ch[0] == 'A') call[i] =  0;
            else if (ch[0] == 'M') call[i] =  1;
            else                   call[i] = -1;
         }//if

         if (hasPVal) {pval[i] = atof(strtok(NULL, sep));}
      }//for_i

      tmptree->Fill();
      idx++;
   }//while
   tmptree->Write();

// Inform user when input table does not contain all genes
   if (XManager::fgVerbose && idx != numgenes) {
      cout << "Warning: Number of lines read <" << idx  
           << "> is not equal to number of genes <" << numgenes << ">" << endl;
      cout << "         for selected scheme <" << scmname.Data() << ">." << endl;
   }//if

// Get branches from temporary tree
   if (hasExpr) {
      exprBr = tmptree->GetBranch("exprBr");
      exprBr->SetAddress(&expr);
   }//if

   if (hasPVal) {
      pvalBr = tmptree->GetBranch("pvalBr");
      pvalBr->SetAddress(&pval);
   }//if

   if (hasCall) {
      callBr = tmptree->GetBranch("callBr");
      callBr->SetAddress(&call);
   }//if

// Create permanent trees
   for (Int_t i=0; i<fNTrees; i++) {
      if (hasExpr) {
         TString tname = fTreeName[i] + "." + xten;
         exprtree = new TTree(tname, scmname);
         if (exprtree == 0) {err = errCreateTree; goto cleanup;}
         exprtree->Branch("ExprBranch", "XExpression", &xpr, bufsize, split);
         exprlist->Add(exprtree);
      }//if

      if (hasCall || hasPVal) {
         TString tname = fTreeName[i] + ".dc5";
         calltree = new TTree(tname, scmname);
         if (calltree == 0) {err = errCreateTree; goto cleanup;}
         calltree->Branch("CallBranch", "XPCall", &pcl, bufsize, split);
         calllist->Add(calltree);
      }//if
   }//for_i

// Fill trees with data from tmptree in the order stored in unittree
   for (id=0; id<numgenes; id++) {
      unittree->GetEntry(id + numctrls); //get unitnames for genes only
      unitname = unit->GetUnitName();

      idxstr = (XIdxString*)(htable->FindObject(unitname));
      if (idxstr) {
         k = idxstr->GetIndex();
         if (hasExpr) {
            exprBr->GetEntry(k);
            for (Int_t i=0; i<fNTrees; i++) {
               exprtree = (TTree*)exprlist->At(i);
               xpr->SetUnitID(id);
               xpr->SetLevel(expr[i]);
               exprtree->Fill();
            }//for_i
         }//if

         if (hasCall || hasPVal) {
            if (hasCall) callBr->GetEntry(k);
            if (hasPVal) pvalBr->GetEntry(k);
            for (Int_t i=0; i<fNTrees; i++) {
               calltree = (TTree*)calllist->At(i);
               pcl->SetUnitID(id);
               pcl->SetCall((hasCall ? call[i] : -1));
               pcl->SetPValue((hasPVal ? pval[i] : -1));
               calltree->Fill();
            }//for_i
         }//if
      } else {
         if (hasExpr) {
            for (Int_t i=0; i<fNTrees; i++) {
               exprtree = (TTree*)exprlist->At(i);
               xpr->SetUnitID(id);
               xpr->SetLevel(-1);
               exprtree->Fill();
            }//for_i
         }//if

         if (hasCall || hasPVal) {
            for (Int_t i=0; i<fNTrees; i++) {
               calltree = (TTree*)calllist->At(i);
               pcl->SetUnitID(id);
               pcl->SetCall(-1);
               pcl->SetPValue(-1);
               calltree->Fill();
            }//for_i
         }//if
      }//if

      if (XManager::fgVerbose && id%10000 == 0) {
         cout << "   <" << id + 1 << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << id << "> records imported." << endl;
   }//if

// Write trees to file
   for (Int_t i=0; i<fNTrees; i++) {
      if (hasExpr) {
         exprtree = (TTree*)exprlist->At(i);
         // write tree to file
         if ((err = WriteTree(exprtree, TObject::kOverwrite)) == errNoErr) {
            // add tree header to list
            AddTreeHeader(exprtree->GetName(), 0);
         }//if
         delete exprtree;
      }//if

      if (hasCall || hasPVal) {
         calltree = (TTree*)calllist->At(i);
         // write tree to file
         if ((err = WriteTree(calltree, TObject::kOverwrite)) == errNoErr) {
            // add tree header to list
            AddTreeHeader(calltree->GetName(), 0);
         }//if
         delete calltree;
      }//if
   }//for_i

// Clean up
cleanup:
   SafeDelete(pcl);
   SafeDelete(xpr);
   SafeDelete(calllist);
   SafeDelete(exprlist);

   delete [] call;
   delete [] pval;
   delete [] expr;

//wrong:   SafeDelete(idxstr);
   if (htable) {htable->Delete(); delete htable; htable = 0;}
   if (buf)    {delete [] buf; buf = 0;}

   //delete temporary tree
   gDirectory->Delete("tmptree;*");

   // remove trees from RAM
   if (unittree) {unittree->Delete(""); unittree = 0;}
   SafeDelete(schemes);

   return err;
}//ReadData

//______________________________________________________________________________
Int_t XGeneChipPivot::ReadXMLHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input in XML-format. 
   if(kCS) cout << "------XGeneChipPivot::ReadXMLHeader------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); sep = 0; delim = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ReadXMLHeader

//______________________________________________________________________________
Int_t XGeneChipPivot::ReadXMLData(ifstream &input, Option_t *option,
                      const char *sep, char delim, Int_t split)
{
   // Read data from input. 
   if(kCS) cout << "------XGeneChipPivot::ReadXMLData------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return 1;
}//ReadXMLData


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChipHyb                                                          //
//                                                                      //
// Class for Mapping array hybridization                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSNPChipHyb::XSNPChipHyb()
            :XGeneChipHyb()
{
   // Default SNPChipHyb constructor
   if(kCS) cout << "---XSNPChipHyb::XSNPChipHyb(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChipHyb::XSNPChipHyb(const char *name, const char *title)
            :XGeneChipHyb(name, title)
{
   // Normal SNPChipHyb constructor
   if(kCS) cout << "---XSNPChipHyb::XSNPChipHyb------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChipHyb::~XSNPChipHyb()
{
   // SNPChipHyb destructor
   if(kCS) cout << "---XSNPChipHyb::~XSNPChipHyb------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChipMetrics                                                      //
//                                                                      //
// Class for processed Mapping array CHP data exported as Metrics.txt   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSNPChipMetrics::XSNPChipMetrics()
                :XGeneChipMetrics()
{
   // Default SNPChipMetrics constructor
   if(kCS) cout << "---XSNPChipMetrics::XSNPChipMetrics(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChipMetrics::XSNPChipMetrics(const char *name, const char *title)
                :XGeneChipMetrics(name, title)
{
   // Normal SNPChipMetrics constructor
   if(kCS) cout << "---XSNPChipMetrics::XSNPChipMetrics------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChipMetrics::~XSNPChipMetrics()
{
   // SNPChipMetrics destructor
   if(kCS) cout << "---XSNPChipMetrics::~XSNPChipMetrics------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChipPivot                                                        //
//                                                                      //
// Class for processed Mapping array CHP data exported as PivotData.txt //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSNPChipPivot::XSNPChipPivot()
              :XGeneChipPivot()
{
   // Default SNPChipPivot constructor
   if(kCS) cout << "---XSNPChipPivot::XSNPChipPivot(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChipPivot::XSNPChipPivot(const char *name, const char *title)
              :XGeneChipPivot(name, title)
{
   // Normal SNPChipPivot constructor
   if(kCS) cout << "---XSNPChipPivot::XSNPChipPivot------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChipPivot::~XSNPChipPivot()
{
   // SNPChipPivot destructor
   if(kCS) cout << "---XSNPChipPivot::~XSNPChipPivot------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChipHyb                                                       //
//                                                                      //
// Class for whole genome array hybridization                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenomeChipHyb::XGenomeChipHyb()
               :XGeneChipHyb()
{
   // Default GenomeChipHyb constructor
   if(kCS) cout << "---XGenomeChipHyb::XGenomeChipHyb(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeChipHyb::XGenomeChipHyb(const char *name, const char *title)
               :XGeneChipHyb(name, title)
{
   // Normal GenomeChipHyb constructor
   if(kCS) cout << "---XGenomeChipHyb::XGenomeChipHyb------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeChipHyb::~XGenomeChipHyb()
{
   // GenomeChipHyb destructor
   if(kCS) cout << "---XGenomeChipHyb::~XGenomeChipHyb------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChipMetrics                                                   //
//                                                                      //
// Class for processed genome array data exported as Metrics.txt        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenomeChipMetrics::XGenomeChipMetrics()
                   :XGeneChipMetrics()
{
   // Default GenomeChipMetrics constructor
   if(kCS) cout << "---XGenomeChipMetrics::XGenomeChipMetrics(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeChipMetrics::XGenomeChipMetrics(const char *name, const char *title)
                   :XGeneChipMetrics(name, title)
{
   // Normal GenomeChipMetrics constructor
   if(kCS) cout << "---XGenomeChipMetrics::XGenomeChipMetrics------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeChipMetrics::~XGenomeChipMetrics()
{
   // GenomeChipMetrics destructor
   if(kCS) cout << "---XGenomeChipMetrics::~XGenomeChipMetrics------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChipPivot                                                     //
//                                                                      //
// Class for processed genome array data exported as PivotData.txt      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenomeChipPivot::XGenomeChipPivot()
                 :XGeneChipPivot()
{
   // Default GenomeChipPivot constructor
   if(kCS) cout << "---XGenomeChipPivot::XGenomeChipPivot(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeChipPivot::XGenomeChipPivot(const char *name, const char *title)
                 :XGeneChipPivot(name, title)
{
   // Normal GenomeChipPivot constructor
   if(kCS) cout << "---XGenomeChipPivot::XGenomeChipPivot------" << endl;

}//Constructor

//______________________________________________________________________________
XGenomeChipPivot::~XGenomeChipPivot()
{
   // GenomeChipPivot destructor
   if(kCS) cout << "---XGenomeChipPivot::~XGenomeChipPivot------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChipHyb                                                         //
//                                                                      //
// Class for Exon array hybridization                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExonChipHyb::XExonChipHyb()
             :XGeneChipHyb()
{
   // Default ExonChipHyb constructor
   if(kCS) cout << "---XExonChipHyb::XExonChipHyb(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XExonChipHyb::XExonChipHyb(const char *name, const char *title)
             :XGeneChipHyb(name, title)
{
   // Normal ExonChipHyb constructor
   if(kCS) cout << "---XExonChipHyb::XExonChipHyb------" << endl;

}//Constructor

//______________________________________________________________________________
XExonChipHyb::~XExonChipHyb()
{
   // ExonChipHyb destructor
   if(kCS) cout << "---XExonChipHyb::~XExonChipHyb------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChipMetrics                                                     //
//                                                                      //
// Class for processed exon array data exported as Metrics.txt          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExonChipMetrics::XExonChipMetrics()
                 :XGeneChipMetrics()
{
   // Default ExonChipMetrics constructor
   if(kCS) cout << "---XExonChipMetrics::XExonChipMetrics(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XExonChipMetrics::XExonChipMetrics(const char *name, const char *title)
                 :XGeneChipMetrics(name, title)
{
   // Normal ExonChipMetrics constructor
   if(kCS) cout << "---XExonChipMetrics::XExonChipMetrics------" << endl;

}//Constructor

//______________________________________________________________________________
XExonChipMetrics::~XExonChipMetrics()
{
   // ExonChipMetrics destructor
   if(kCS) cout << "---XExonChipMetrics::~XExonChipMetrics------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChipPivot                                                       //
//                                                                      //
// Class for processed exon array data exported as PivotData.txt        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExonChipPivot::XExonChipPivot()
               :XGeneChipPivot()
{
   // Default ExonChipPivot constructor
   if(kCS) cout << "---XExonChipPivot::XExonChipPivot(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XExonChipPivot::XExonChipPivot(const char *name, const char *title)
               :XGeneChipPivot(name, title)
{
   // Normal ExonChipPivot constructor
   if(kCS) cout << "---XExonChipPivot::XExonChipPivot------" << endl;

}//Constructor

//______________________________________________________________________________
XExonChipPivot::~XExonChipPivot()
{
   // ExonChipPivot destructor
   if(kCS) cout << "---XExonChipPivot::~XExonChipPivot------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenePixHyb                                                          //
//                                                                      //
// Class for GenePix microarray hybridization                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenePixHyb::XGenePixHyb()
            :XHybridization()
{
   // Default GenePixHyb constructor
   if(kCS) cout << "---XGenePixHyb::XGenePixHyb(default)------" << endl;

   fNRows = 0; 
   fNCols = 0;
}//Constructor

//______________________________________________________________________________
XGenePixHyb::XGenePixHyb(const char *name, const char *title)
            :XHybridization(name, title)
{
   // Normal GenePixHyb constructor
   if(kCS) cout << "---XGenePixHyb::XGenePixHyb------" << endl;

   fNRows = 0; 
   fNCols = 0;
}//Constructor

//______________________________________________________________________________
XGenePixHyb::~XGenePixHyb()
{
   // GenePixHyb destructor
   if(kCS) cout << "---XGenePixHyb::~XGenePixHyb------" << endl;

}//Destructor

//______________________________________________________________________________
void XGenePixHyb::AddDataTreeInfo(TTree *tree, const char *name, Option_t *option,
                  Int_t nquant, Double_t *q, Double_t *quant)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGenePixHyb::AddDataTreeInfo------" << endl;

// store name of tree set as title
   XDataTreeInfo *info = new XDataTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   if (nquant > 0) info->AddUserInfo(nquant, q, quant);

   tree->GetUserInfo()->Add(info);
}//AddDataTreeInfo

//______________________________________________________________________________
Int_t XGenePixHyb::ExportTreeType(const char *exten, Int_t n, TString *names, 
                   const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output
   if(kCS) cout << "------XGenePixHyb::ExportTreeType------" << endl;

   if (strcmp(exten, "cel") == 0) {
      return this->ExportDataTrees(n, names, varlist, output, sep);
   } else if (strcmp(exten, "msk") == 0) {
      return this->ExportMaskTrees(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//ExportTreeType

//______________________________________________________________________________
Int_t XGenePixHyb::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                   const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XGenePixHyb::ExportTreeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   exten = 0; n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return errNoErr;
}//ExportTreeXML

//______________________________________________________________________________
void XGenePixHyb::PrintInfo()
{
   // Print GenePix hybridization information
   if(kCS) cout << "------XGenePixHyb::PrintInfo------" << endl;

   if (fgPrintHeader) {
      cout << "==============================================================================" << endl;
      cout << setw(14) << "ChipName" << setw(12) << "Title"
           << endl;
      cout << "==============================================================================" << endl;
      fgPrintHeader = kFALSE;
   }//if

   cout << setw(14) << this->GetName() << setw(12) << this->GetTitle()

        << endl;

   cout << "------------------------------------------------------------------------------" << endl;
}//PrintInfo

//______________________________________________________________________________
Int_t XGenePixHyb::ExportDataTrees(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGenePixHyb::ExportDataTrees------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportDataTrees

//______________________________________________________________________________
Int_t XGenePixHyb::ExportDataTree(ofstream &output, const char *sep)
{
   // Export data stored in tree to file output
   if(kCS) cout << "------XGenePixHyb::ExportDataTree------" << endl;

// Output header
//TO DO
   output << "X" << sep << "Y" << endl;

// Get tree
   XFeature635 *f635 = 0;
   XFeature532 *f532 = 0;
   TTree *tree = (TTree*)gDirectory->Get((fName + ".gpr")); 
   if (tree == 0) return errGetTree;
   tree->SetBranchAddress("DataBranch", &f635);
   tree->SetBranchAddress("DataBranch", &f532);

   Int_t entries = (Int_t)(tree->GetEntries());
// Loop over all entries of tree
   for (Int_t i=0; i<entries; i++) {
      tree->GetEntry(i);
//      output << f635->Get() << sep
//             << f532->Get() << endl;
   }//for_i

   tree->Delete("");
   tree = 0;

   return errNoErr;
}//ExportDataTree

//______________________________________________________________________________
Int_t XGenePixHyb::ExportMaskTrees(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export mask tree to file output
   if(kCS) cout << "------XGenePixHyb::ExportMaskTrees------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportMaskTrees

//______________________________________________________________________________
Int_t XGenePixHyb::ExportMaskTree(ofstream &output, const char *sep)
{
   // Export mask tree to file output
   if(kCS) cout << "------XGenePixHyb::ExportMaskTree------" << endl;

// Output header
   output << "X" << sep << "Y" << sep << "Flag" << endl;

// Get tree
   XMask *mask = 0;
   TTree *tree = (TTree*)gDirectory->Get((fName + ".msk")); 
   if (tree == 0) return errGetTree;
   tree->SetBranchAddress("MaskBranch", &mask);

   Int_t entries = (Int_t)(tree->GetEntries());
// Loop over all entries of tree
   for (Int_t i=0; i<entries; i++) {
      tree->GetEntry(i);
//      output << mask->GetFlag() << sep
//             << mask->GetFlag() << endl;
   }//for_i

   tree->Delete("");
   tree = 0;

   return errNoErr;
}//ExportMaskTree

//______________________________________________________________________________
Int_t XGenePixHyb::ReadHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input. 
   if(kCS) cout << "------XGenePixHyb::ReadHeader------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); sep = 0; delim = 0;

   return 0;
}//ReadHeader

//______________________________________________________________________________
Int_t XGenePixHyb::ReadData(ifstream &input, Option_t *option, const char *sep,
                   char delim, Int_t split)
{
   // Read data from input. 
   if(kCS) cout << "------XGenePixHyb::ReadData------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return 0;
}//ReadData

//______________________________________________________________________________
Int_t XGenePixHyb::ReadXMLHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input in XML-format. 
   if(kCS) cout << "------XGenePixHyb::ReadXMLHeader------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); sep = 0; delim = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ReadXMLHeader

//______________________________________________________________________________
Int_t XGenePixHyb::ReadXMLData(ifstream &input, Option_t *option, const char *sep,
                   char delim, Int_t split)
{
   // Read data from input. 
   if(kCS) cout << "------XGenePixHyb::ReadXMLData------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return 1;
}//ReadXMLData


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalysisPlot                                                        //
//                                                                      //
// Class for drawing microarray results                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAnalysisPlot::XAnalysisPlot()
              :XPlot()
{
   // Default AnalysisPlot constructor
   if(kCS) cout << "---XAnalysisPlot::XAnalysisPlot(default)------" << endl;

   fSchemeFile = 0;
   fDataFile   = 0;
   fSchemes    = 0;
   fData       = 0;
   fIsSchemeOwner = kFALSE;
   fIsDataOwner   = kFALSE;
}//Constructor

//______________________________________________________________________________
XAnalysisPlot::XAnalysisPlot(const char *name, const char *title)
              :XPlot(name, title)
{
   // Normal AnalysisPlot constructor
   if(kCS) cout << "---XAnalysisPlot::XAnalysisPlot------" << endl;

   fSchemeFile = 0;
   fDataFile   = 0;
   fSchemes    = 0;
   fData       = 0;
   fIsSchemeOwner = kFALSE;
   fIsDataOwner   = kFALSE;
}//Constructor

//______________________________________________________________________________
XAnalysisPlot::~XAnalysisPlot()
{
   // AnalysisPlot destructor
   if(kCS) cout << "---XAnalysisPlot::~XAnalysisPlot------" << endl;

   this->CloseData();
   this->CloseSchemes();
}//Destructor

//______________________________________________________________________________
Int_t XAnalysisPlot::InitData(TFile *datafile, Bool_t isOwner)
{
   // Initialize analysis plot with external root datafile
   // Set isOwner = kTRUE if datafile should be deleted by processmanager
   if(kCS) cout << "------XAnalysisPlot::InitData------" << endl;

   if (fAbort) return errAbort;
   TDirectory *savedir = gDirectory;

   // check if datafile is already open
   if (this->IsOpen(fDataFile, datafile->GetName())) {
      // only one data file can be open
      cout << "Closing existing data file <" << fDataFile->GetName()
           << ">..." << endl;
      this->CloseData();
   }//if

   fDataFile = datafile;
   if (!fDataFile) {fAbort = kTRUE; return errGetFile;}
   fIsDataOwner = isOwner;

   fDataFile->cd();
   fData = (XFolder*)(fDataFile->Get(kContent));
   if (!fData) {
      cerr << "Error: Data index <" << kContent << "> is missing" << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   savedir->cd();
   return errNoErr;
}//InitData

//______________________________________________________________________________
Int_t XAnalysisPlot::InitSchemes(TFile *schemefile, Bool_t isOwner)
{
   // Initialize analysis plot with external root schemefile
   // Set isOwner = kTRUE if schemefile should be deleted by processmanager
   if(kCS) cout << "------XAnalysisPlot::InitSchemes------" << endl;

   if (fAbort) return errAbort;
   TDirectory *savedir = gDirectory;

   // check if schemefile is already open
   if (this->IsOpen(fSchemeFile, schemefile->GetName())) {
      // only one scheme file can be open
      cout << "Closing existing scheme file <" << fSchemeFile->GetName()
           << ">..." << endl;
      this->CloseSchemes();
   }//if

   fSchemeFile = schemefile;
   if (!fSchemeFile) {fAbort = kTRUE; return errGetFile;}
   fIsSchemeOwner = isOwner;

   fSchemeFile->cd();
   fSchemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!fSchemes) {
      cerr << "Error: Schemes index <" << kContent << "> is missing" << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   savedir->cd();
   return errNoErr;
}//InitSchemes

//______________________________________________________________________________
Int_t XAnalysisPlot::OpenData(const char *fullname)
{
   // Open root data file
   if(kCS) cout << "------XAnalysisPlot::OpenData------" << endl;

   if (fAbort) return errAbort;

   // check if datafile is already open
   if (this->IsOpen(fDataFile, fullname)) {
      // only one scheme file can be open
      cout << "Closing existing data file <" << fDataFile->GetName()
           << ">..." << endl;
      this->CloseData(); //deletes also fData!
   }//if

   Bool_t isOwner = kFALSE;
   fDataFile = OpenFile(fullname, "READ", isOwner);
   if (!fDataFile) {fAbort = kTRUE; return errAbort;}
   // assure that manager remains owner if it calls OpenFile multiple times
   if (!fIsDataOwner) fIsDataOwner = isOwner;

   fDataFile->cd();
   fData = (XFolder*)(fDataFile->Get(kContent));
   if (!fData) {
      cerr << "Error: Data index <" << kContent << "> is missing" << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   return errNoErr;
}//OpenData

//______________________________________________________________________________
Int_t XAnalysisPlot::OpenSchemes(const char *fullname)
{
   // Open root scheme file
   if(kCS) cout << "------XAnalysisPlot::OpenSchemes------" << endl;

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
   if (!fSchemeFile) {fAbort = kTRUE; return errAbort;}
   // assure that manager remains owner if it calls OpenFile multiple times
   if (!fIsSchemeOwner) fIsSchemeOwner = isOwner;

   fSchemeFile->cd();
   fSchemes = (XFolder*)(fSchemeFile->Get(kContent));
   if (!fSchemes) {
      cerr << "Error: Schemes index <" << kContent << "> is missing" << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   savedir->cd();
   return errNoErr;
}//OpenSchemes

//______________________________________________________________________________
void XAnalysisPlot::CloseData()
{
   // Close root data file
   if(kCS) cout << "------XAnalysisPlot::CloseData------" << endl;

   SafeDelete(fData);

   if (fIsDataOwner) SafeDelete(fDataFile);
   fDataFile = 0;  //if not isOwner!
}//CloseData

//______________________________________________________________________________
void XAnalysisPlot::CloseSchemes()
{
   // Close scheme file
   if(kCS) cout << "------XAnalysisPlot::CloseSchemes------" << endl;

   SafeDelete(fSchemes);

   if (fIsSchemeOwner) SafeDelete(fSchemeFile);
   fSchemeFile = 0;  //if not isOwner!
}//CloseSchemes


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipPlot                                                        //
//                                                                      //
// Class for drawing GeneChip data                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGeneChipPlot::XGeneChipPlot()
              :XAnalysisPlot()
{
   // Default GeneChipPlot constructor
   if(kCS) cout << "---XGeneChipPlot::XGeneChipPlot(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGeneChipPlot::XGeneChipPlot(const char *name, const char *title)
              :XAnalysisPlot(name, title)
{
   // Normal GeneChipPlot constructor
   if(kCS) cout << "---XGeneChipPlot::XGeneChipPlot------" << endl;

}//Constructor

//______________________________________________________________________________
XGeneChipPlot::~XGeneChipPlot()
{
   // GeneChipPlot destructor
   if(kCS) cout << "---XGeneChipPlot::~XGeneChipPlot------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGeneChipPlot::DrawUnit(const char *canvasname, const char *treename,
                     const char *schemename, const char *unitname, Int_t unitID,
                     const char *varlist, const char *logbase, const char *type,
                     Option_t *opt, Double_t min, Double_t max,
                     const char *var2sort, Bool_t down)
{
   // Draw raw data from tree treename for scheme with schemename and unit with
   // unitname or for unitID (if unitname = "").
   // Variables to be drawn in one graph are listed in
   // varlist. The following variables are allowed:
   //   inPM, inMM - intensity of PerfectMatch / MisMatch oligonucleotide
   //   sdPM, sdMM - standard deviation of PM / MM oligonucleotide
   //   spPM, spMM - standard deviation in [%] of PM / MM oligonucleotide
   //   pxPM, pxMM - number of pixels of PM / MM oligonucleotide
   //   gcPM, gcMM - number of GCs in PM / MM oligonucleotide
   //   bgPM, bgMM - background value of PM / MM oligonucleotide
   //   ibPM, ibMM - background-subtracted intensity of PM / MM oligonucleotide
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   //   note: number of pixels and GCs not converted to log
   // type: graph, multigraph, grapherr, hist with the corresponding option opt
   //   e.g. for "inPM:sdPM:inMM:sdMM":
   //   tpye = multigraph: four curves are drawn for intensities and stddevs
   //   tpye = hist: four histograms are drawn for intensities and stddevs
   //   tpye = grapherr: varlist = "inPM:inMM"  (sdPM, sdMM added automatically)
   //                    two curves are drawn for intensities +/- stddevs
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Variable var2sort determines the variable to sort for
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XGeneChipPlot::DrawUnit(1)------" << endl;

   if (fAbort) return perrAbort;
   TDirectory *savedir = gDirectory;

   Int_t perr = perrNoErr;

// Create canvas for drawing
   if (strcmp(canvasname, "") != 0) {
      this->NewCanvas(canvasname, "");
   }//if

// Set pad number
   if (fNPads > 1) fPadNr++;
   if (fPadNr > fNPads) {
      cerr << "Error: Number of pads <" << fPadNr << "> is larger than <" 
           << fNPads << ">." << endl;
      return perrNumPads;
   }//if

// Extract tree name
   TString tname = Path2Name(treename, dSEP, "");
   if (strstr(tname.Data(), ".root")) {
      tname = "";
   }//if
   if (strcmp(tname.Data(), "") == 0) {
      cerr << "Error: Treename for tree is missing." << endl;
      return perrTreeName;
   }//if

// Test if tree extension is "cel"
   if (strstr(tname.Data(), ".cel") == 0) {
      cerr << "Error: Treename must be of form <treename.cel>." << endl;
      return perrTreeName;
   }//if

// Extract root filename and get datafile
   TString filename = "";
   if (strstr(treename, ".root")) {
      filename = Path2Name(treename, "", ".root") + ".root";
      this->OpenData(filename);
   } else if (fDataFile) {
//?? not necessary?
      filename = fDataFile->GetName();
   } else {
      cerr << "Error: Data file is not open." << endl;
      return perrOpenFile;
   }//if

// Get name of treeset and change directory
   TString setname  = "";
   if (strstr(treename, ".root")) {
      TString substr = SubString(treename, '.', sSEP, kFALSE);
      if (substr) setname = Path2Name(substr.Data(), dSEP, "");
      if (setname.Contains("root")) setname = "";
   } else if (strstr(treename, dSEP)) {
      setname = Path2Name(treename, "", dSEP);
   }//if

   if (!fDataFile->cd(setname)) return perrGetDir;

// Get datatree
   XGCCell *cell = 0;
   TTree *datatree = (TTree*)gDirectory->Get(tname); 
   if (datatree == 0) return perrGetTree;
   datatree->SetBranchAddress("DataBranch", &cell);

// Get scheme name from datatree and change dir
   TString scmname = TString(schemename);
   if (strcmp(schemename, "") == 0) {
      scmname = datatree->GetTitle();
   } else if (!scmname.Contains(datatree->GetTitle())) {
      return errChipType;
   }//if
   if (!fSchemeFile->cd(scmname)) return errGetDir;

// Get chip parameters from scheme file
   XDNAChip *chip = 0;
   chip = (XDNAChip*)fSchemes->FindObject(scmname, kTRUE);
   if (!chip) {
      cerr << "Error: Could not find scheme <" << scmname.Data() << ">." << endl;
      return perrGeneral;
   }//if
   Int_t numunits = chip->GetNumUnits();
   Int_t numctrls = chip->GetNumControls();
   Int_t numcols  = chip->GetNumColumns();

// Get scheme tree for scheme
   XGCScheme *scheme = 0;
   TTree *scmtree = (TTree*)gDirectory->Get(chip->GetSchemeTree()); 
   if (scmtree == 0) return perrGetTree;
   scmtree->SetBranchAddress("ScmBranch", &scheme);

// Get unit tree for scheme
   XGCUnit *unit = 0;
   TTree *idxtree = (TTree*)gDirectory->Get(chip->GetUnitTree()); 
   if (idxtree == 0) return perrGetTree;
   idxtree->SetBranchAddress("IdxBranch", &unit);

// Get probe tree for scheme
   Bool_t hasProbe = kTRUE;
   XGCProbe *probe = 0;
   TTree *prbtree = (TTree*)gDirectory->Get(chip->GetProbeTree()); 
   if (prbtree == 0) {
      cout << "Warning: Scheme <" << scmname.Data() << "> has no probe information."
           << endl;
      hasProbe = kFALSE;
   } else {
      prbtree->SetBranchAddress("PrbBranch", &probe);
   }//if

// Get annotation tree for scheme
   Bool_t hasAnnot = kTRUE;
   XTransAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get(chip->GetAnnotTree()); 
   if (anntree == 0) {
      cout << "Warning: Scheme <" << scmname.Data() << "> has no annotation."
           << endl;
      hasAnnot = kFALSE;
   } else {
      anntree->SetBranchAddress("AnnBranch", &annot);
   }//if

// Get unitID, unitname and numpairs
   Int_t numcells = 0;
   for (Int_t id=0; id<numunits; id++) { 
      idxtree->GetEntry(id);
      if (strcmp(unitname, "") == 0) {
         if (id == (unitID + numctrls)) {
            unitname = unit->GetUnitName();
            numcells = unit->GetNumCells();
            break;
         }//if
      } else if (strcmp(unitname,unit->GetUnitName()) == 0) {
         unitID = id - numctrls;  //first entries are ctrls with negative IDs 
         numcells = unit->GetNumCells();
         break;
      }//if
   }//for_id
//TO DO
//need to update to work also for Gene/Exon arrays: use THashTable!!

//TO TEST!!!
// Return if unitID/unitname are not correct
   if (numcells == 0) return 1;
//to do   if (numcells == 0) return perrGetUnit;
   Int_t numpairs = (Int_t)(numcells / 2);

// Get annotation for unitID
   TString genename = "";
   TString symbol   = "";
   if (hasAnnot) {
      anntree->GetEntry(unitID + numctrls);
      genename = annot->GetName();
      symbol   = annot->GetSymbol();
   }//if
   TString title = tname + ": " + TString(unitname) + "  -  " + symbol;

   TString *str = new TString[16];
   for (Int_t i=0; i<16; i++) str[i] = "";

   Int_t *ord = new Int_t[16];
   for (Int_t i=0; i<16; i++) ord[i] = kMaxInt;

// Decompose varlist: store order of variables
   Int_t  numvar = 0;
   Int_t  sort   = 0;
//ccc char memory problem??
   char  *var    = new char[strlen(varlist) + 1];
   char  *del    = var;
   var = strtok(strcpy(var,varlist),":");
   while (var) {
      if (strcmp(var,"inPM") == 0) {ord[0]  = numvar; str[numvar] = var;}
      if (strcmp(var,"sdPM") == 0) {ord[1]  = numvar; str[numvar] = var;}
      if (strcmp(var,"spPM") == 0) {ord[2]  = numvar; str[numvar] = var;}
      if (strcmp(var,"bgPM") == 0) {ord[3]  = numvar; str[numvar] = var;}
      if (strcmp(var,"ibPM") == 0) {ord[4]  = numvar; str[numvar] = var;}
      if (strcmp(var,"pxPM") == 0) {ord[5]  = numvar; str[numvar] = var;}
      if (strcmp(var,"gcPM") == 0) {ord[6]  = numvar; str[numvar] = var;}
      if (strcmp(var,"inMM") == 0) {ord[7]  = numvar; str[numvar] = var;}
      if (strcmp(var,"sdMM") == 0) {ord[8]  = numvar; str[numvar] = var;}
      if (strcmp(var,"spMM") == 0) {ord[9]  = numvar; str[numvar] = var;}
      if (strcmp(var,"bgMM") == 0) {ord[10] = numvar; str[numvar] = var;}
      if (strcmp(var,"ibMM") == 0) {ord[11] = numvar; str[numvar] = var;}
      if (strcmp(var,"pxMM") == 0) {ord[12] = numvar; str[numvar] = var;}
      if (strcmp(var,"gcMM") == 0) {ord[13] = numvar; str[numvar] = var;}
      if (strcmp(var2sort,(str[numvar]).Data()) == 0) {sort = numvar + 1;}
      numvar++;
      var = strtok(NULL, ":");
      if (var == 0) break;
   }//while
   delete [] del;

// Get logbase
   Int_t base = -1; 
   if (strcmp(logbase,"0")          == 0) {base = 0;}
   else if (strcmp(logbase,"log")   == 0) {base = 1;}
   else if (strcmp(logbase,"log2")  == 0) {base = 2;}
   else if (strcmp(logbase,"log10") == 0) {base = 10;}
   else {
      cout << "Warning: Logbase not known, using default." << endl;
      base = 0;
   }//if

   Int_t x, y, ij;
   Int_t numpix;
   Short_t flag;
   Int_t  entries = 0;
   Double_t inten = 0;
   Double_t stdev = 0;

// Get GC from probe tree!!!
   Int_t numgc = 0;

   Int_t p = 0;
   Int_t m = 0;
   Double_t bgd      = 0;
   Double_t maxinten = 0;
   Double_t maxstdev = 0;
   Double_t maxbgd   = 0;
   Double_t mininten = DBL_MAX;
   Double_t minbgd   = DBL_MAX;

// Init arrays for graph and hist
   Double_t *index = new Double_t[numpairs];
   Double_t *arrX  = new Double_t[numpairs];
   Double_t *arrY  = new Double_t[numpairs];
   Double_t *arrZ  = new Double_t[numpairs];

// Init table for data
   Double_t **table = 0;
   if (!(table = new (nothrow) Double_t*[numvar])) {
      fAbort = kTRUE; perr = perrInitMemory; goto cleanup;
   }//if 
   for (Int_t i=0; i<numvar; i++) {
      table[i] = 0;
      if (!(table[i] = new (nothrow) Double_t[numpairs])) {
         fAbort = kTRUE; perr = perrInitMemory; goto cleanup;
      }//if 
   }//for_i
   for (Int_t i=0;i<numvar;i++) {
      for (Int_t j=0;j<numpairs;j++) {
         table[i][j] = j;
      }//for_j
   }//for_i

// Get number of entries
   entries = (Int_t)(scmtree->GetEntries());
   if ((Int_t)(datatree->GetEntries() == entries)){
      cout << "Number of entries = " << entries << endl;
   } else {
//PROBLEM???????
/*Alternative CDFs????????????
      cerr << "Error: number of entries in scheme tree not equal to datatree"
           << endl;
      perr = perrNumEntries;
      goto cleanup;
*/
   }//if

// Fill table
   for (Int_t i=0; i<entries; i++) {
      scmtree->GetEntry(i);
      if (scheme->GetUnitID()!= unitID) continue;
      flag = scheme->GetMask();

      // get number of GC
      if (hasProbe) {
         prbtree->GetEntry(i);
         numgc = probe->GetNumberGC();
      }//if

     // get data for entry (x,y)
      x = scheme->GetX();
      y = scheme->GetY();
      ij = y*numcols + x;   //*.cel file loops first for x then for y
//xy   ev:   ij = XY2Index(x, y, numcols);
      datatree->GetEntry(ij);
      if (base == 0) {
         inten  = cell->GetIntensity();
         stdev  = cell->GetStdev();
      } else if (base == 2) {
         inten  = TMath::Log2(cell->GetIntensity());
         stdev  = TMath::Log2(cell->GetStdev());
      } else if (base == 10) {
         inten  = TMath::Log10(cell->GetIntensity());
         stdev  = TMath::Log10(cell->GetStdev());
      } else if (base == 1) {
         inten  = TMath::Log(cell->GetIntensity());
         stdev  = TMath::Log(cell->GetStdev());
      }//if
      numpix = cell->GetNumPixels();

//??      tree->GetEntry(i);
//??      bgd  = ???->GetBackground();  //also log!

      if (flag == 1) {
         if (ord[0] < kMaxInt) table[ord[0]][p] = inten;
         if (ord[1] < kMaxInt) table[ord[1]][p] = stdev;
         if (ord[2] < kMaxInt) table[ord[2]][p] = 100.0*stdev/inten;
         if (ord[3] < kMaxInt) table[ord[3]][p] = bgd;
         if (ord[4] < kMaxInt) table[ord[4]][p] = inten - bgd;
         if (ord[5] < kMaxInt) table[ord[5]][p] = (Double_t)numpix;
         if (ord[6] < kMaxInt) table[ord[6]][p] = (Double_t)numgc;
         p++;
      } else {
         if (ord[7]  < kMaxInt) table[ord[7]][m]  = inten;
         if (ord[8]  < kMaxInt) table[ord[8]][m]  = stdev;
         if (ord[9]  < kMaxInt) table[ord[9]][m]  = 100.0*stdev/inten;
         if (ord[10] < kMaxInt) table[ord[10]][m] = bgd;
         if (ord[11] < kMaxInt) table[ord[11]][m] = inten - bgd;
         if (ord[12] < kMaxInt) table[ord[12]][m] = (Double_t)numpix;
         if (ord[13] < kMaxInt) table[ord[13]][m] = (Double_t)numgc;
         m++;
      }//if

      maxinten = (maxinten < inten) ? inten : maxinten;
      mininten = (mininten > inten) ? inten : mininten;
      maxstdev = (maxstdev < stdev) ? stdev : maxstdev;
      maxbgd   = (maxbgd < bgd)     ? bgd   : maxbgd;
      minbgd   = (minbgd > bgd)     ? bgd   : minbgd;
   }//for_i

   fMin = fMinX = fMinY = fMinZ = DBL_MAX;
   fMax = fMaxX = fMaxY = fMaxZ = 0;
   fNBinsX = fNBinsY = fNBinsZ = p; //??

//to do: needs optimization?
// Fill arrays in the order of varlist
   for (Int_t i=0; i<numpairs; i++) {
      index[i] = i;
      if (numvar > 0) {
         arrX[i] = table[0][i];
         fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
         fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         fMin  = (fMin  < fMinX)   ? fMin  : fMinX;
         fMax  = (fMax  > fMaxX)   ? fMax  : fMaxX;
      }//if
      if (numvar > 1) {
         arrY[i] = table[1][i];
         fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
         fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         fMin  = (fMin  < fMinY)   ? fMin  : fMinY;
         fMax  = (fMax  > fMaxY)   ? fMax  : fMaxY;
      }//if
      if (numvar > 2) {
         arrZ[i] = table[2][i];
         fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
         fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
         fMin  = (fMin  < fMinZ)   ? fMin  : fMinZ;
         fMax  = (fMax  > fMaxZ)   ? fMax  : fMaxZ;
      }//if
   }//for_i

// Set axes range
   if ((min == 0) && (max == 0)) {
      fEqualAxes = kTRUE;
   } else if ((min == -1111) && (max == -1111)){
      fEqualAxes = kFALSE;
   } else if (min < max){
      fMin = min;
      fMax = max;
      fEqualAxes = kTRUE;
   } else {
      cout << "Warning: min > max, thus using computed values" << endl;
      fEqualAxes = kFALSE;
   }//if
//?? in GUI: display fMinX,fMinY,etc

//?? ev title: AffyID: gene name (anntree)
//?? ev X-axis gives index in sorted order, e.g. 2 4 3 1 5 6 9 ...
//?? hist for bar option!! (nbins = npairs!!)
//?? numgc need to scale axis !!

// Draw type
   SetTitleMain(title.Data(), fSetTitle);
   if (strcmp(type, "graph") == 0) {
      switch (numvar) {
         case 1:
            SetTitleX("Index", fSetTitleX, 0);
            SetTitleY(str[0], fSetTitleY, base);
            fMinX = 0;     fMaxX = numpairs - 1;
            fMinY = fMinX; fMaxY = fMaxX;
            if (sort != 0) {sort = down ? -1 : 1;}
            DrawGraph1D(numpairs, index, arrX, opt, sort);
            break;

         case 2:
            SetTitleX(str[0], fSetTitleX, base);
            SetTitleY(str[1], fSetTitleY, base);
            DrawGraph2D(numpairs, arrX, arrY, opt);
            break;

         case 3:
            SetTitleX(str[0], fSetTitleX, base);
            SetTitleY(str[1], fSetTitleY, base);
            SetTitleZ(str[2], fSetTitleZ, base);
            DrawGraph3D(numpairs, arrX, arrY, arrZ, opt);
            break;

         default:
            printf("Dimension <%d> not allowed\n", numvar);
            perr = perrGeneral;
      }//switch
   } else if (strcmp(type, "multigraph") == 0) {
      switch (numvar) {
         case 1:
            SetTitleX("Index", fSetTitleX, base);
            SetTitleY(str[0], fSetTitleY, base);
            fMinX = 0;     fMaxX = numpairs - 1;
            fMinY = fMinX; fMaxY = fMaxX;
            if (sort != 0) {sort = down ? -1 : 1;}
            DrawGraph1D(numpairs, index, arrX, opt, sort);
            break;

         case 2:
            SetTitleX("Index", fSetTitleX, 0);
            SetTitleY(str[0], fSetTitleY, base);
            SetTitleZ(str[1], fSetTitleZ, base);
            DrawMultiGraph(numpairs, arrX, arrY, opt, sort, down);
            break;

         case 3:
            SetTitleX("Index", fSetTitleX, 0);
            SetTitleY(str[0], fSetTitleY, base);
            SetTitleZ(str[1], fSetTitleZ, base);
            DrawMultiGraph(numpairs, arrX, arrY, arrZ, opt, sort, down);
            break;

         default:
            SetTitleX("Index", fSetTitleX, 0);
            SetTitleY(str[0], fSetTitleY, base);
            SetTitleZ(str[1], fSetTitleZ, base);
            DrawMultiGraph(numpairs, numvar, table, opt, sort, down);
      }//switch
   } else if (strcmp(type, "hist") == 0) {
      switch (numvar) {
         case 1:
            SetTitleX(str[0], fSetTitleX, base);
            SetTitleY("", fSetTitleY, base);
            DrawHist1D(numpairs, index, arrX, opt);
            break;

         case 2:
            SetTitleX(str[0], fSetTitleX, base);
            SetTitleY(str[1], fSetTitleY, base);
            DrawHist2D(numpairs, arrX, arrY, opt);
            break;

         case 3:
            SetTitleX(str[0], fSetTitleX, base);
            SetTitleY(str[1], fSetTitleY, base);
            SetTitleZ(str[2], fSetTitleZ, base);
            DrawHist3D(numpairs, arrX, arrY, arrZ, opt);
            break;

         default:
            printf("Dimension <%d> not allowed\n", numvar);
            perr = perrGeneral;
      }//switch
   } else {
      cout << "Error: Drawing type <" << type << "> not known" << endl;
      perr = perrPlotType;
   }//if

// Cleanup
cleanup:
   delete [] ord;
   delete [] str;
   delete [] arrZ;
   delete [] arrY;
   delete [] arrX;
   delete [] index;

   for (Int_t i=0; i<numvar; i++) {
      if (table[i]) {delete [] table[i]; table[i] = 0;}
   }//for_i
   delete [] table;

   // remove trees from RAM
   if (anntree)  {anntree->Delete("");  anntree  = 0;}
   if (prbtree)  {prbtree->Delete("");  prbtree  = 0;}
   if (idxtree)  {idxtree->Delete("");  idxtree  = 0;}
   if (scmtree)  {scmtree->Delete("");  scmtree  = 0;}
   if (datatree) {datatree->Delete(""); datatree = 0;}

// chip and hyb deleted in CloseSchemes()
//   SafeDelete(chip); //not possible if DrawUnit called more than once
//   SafeDelete(hyb); //not possible if DrawUnit called more than once

   savedir->cd();
   return perr;
}//DrawUnit


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenePixPlot                                                         //
//                                                                      //
// Class for drawing GenePix microarray data                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenePixPlot::XGenePixPlot()
             :XAnalysisPlot()
{
   // Default GenePixPlot constructor
   if(kCS) cout << "---XGenePixPlot::XGenePixPlot(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGenePixPlot::XGenePixPlot(const char *name, const char *title)
             :XAnalysisPlot(name, title)
{
   // Normal GenePixPlot constructor
   if(kCS) cout << "---XGenePixPlot::XGenePixPlot------" << endl;

}//Constructor

//______________________________________________________________________________
XGenePixPlot::~XGenePixPlot()
{
   // GenePixPlot destructor
   if(kCS) cout << "---XGenePixPlot::~XGenePixPlot------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGenePixPlot::DrawUnit(const char *canvasname, const char *treename,
                    const char *schemename, const char *unitname, Int_t unitID,
                    const char *varlist, const char *logbase, const char *type,
                    Option_t *opt, Double_t min, Double_t max, const char *var2sort,
                    Bool_t down)
{
   // Draw raw data from tree treename for unit with unitname or for unitID
   // (if unitname = ""). Variables to be drawn are listed in varlist.
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // type: graph, multigraph, grapherr, hist with the corresponding option opt
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Variable var2sort determines the variable to sort for
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XGenePixPlot::DrawUnit(1)------" << endl;

// NOT YET USED - to prevent compiler warnings:
   canvasname = 0; treename = 0; schemename = 0; unitname = 0; unitID = 0;
   varlist = 0; logbase = 0; type = 0; opt = 0; min = 0; max = 0; var2sort = 0; down = 0;

   if (fAbort) return perrAbort;
   TDirectory *savedir = gDirectory;

   Int_t perr = perrNoErr;

   savedir->cd();
   return perr;
}//DrawUnit


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCell                                                                //
//                                                                      //
// Base class describing a single microarray cell                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XCellMask::XCellMask()
          :XPosition()
{
   // Default CellMask constructor
   fFlag = -1;
}//Constructor

//______________________________________________________________________________
XCellMask::~XCellMask()
{
   // CellMask destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCell                                                                //
//                                                                      //
// Base class describing a single microarray cell                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XCell::XCell()
{
   // Default Cell constructor
   fInten   = -1;
   fStdev   = -1;
   fNPixels = -1;
}//Constructor

//______________________________________________________________________________
XCell::~XCell()
{
   // Cell destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCCell                                                              //
//                                                                      //
// Class describing a single GeneChip cell                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGCCell::XGCCell()
        :XPosition()
{
   // Default GCCell constructor
   fInten   = -1;
   fStdev   = -1;
   fNPixels = -1;
}//Constructor

//______________________________________________________________________________
XGCCell::~XGCCell()
{
   // GCCell destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBgCell                                                              //
//                                                                      //
// Class containing background data for microarray cell                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBgCell::XBgCell()
        :XPosition()
{
   // Default BgCell constructor
   fBg    = 0;
   fStdev = 0;
}//Constructor

//______________________________________________________________________________
XBgCell::~XBgCell()
{
   // BgCell destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRankCell                                                            //
//                                                                      //
// Class containing intensity and its rank for microarray cell          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRankCell::XRankCell()
          :XPosition()
{
   // Default BgCell constructor
   fInten = 0;
   fRank  = 0;
}//Constructor

//______________________________________________________________________________
XRankCell::~XRankCell()
{
   // BgCell destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XResidual                                                            //
//                                                                      //
// Class describing microarray residuals an weights                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XResidual::XResidual()
          :XPosition()
{
   // Default Residual constructor
   fResidual = 0.0;
   fWeight   = 0.0;
}//Constructor

//______________________________________________________________________________
XResidual::~XResidual()
{
   // Residual destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBorder                                                              //
//                                                                      //
// Class describing microarrayborder elements                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBorder::XBorder()
        :XPosition()
{
   // Default Border constructor
   fInten = 0.0;
   fFlag  = 0.0;
}//Constructor

//______________________________________________________________________________
XBorder::~XBorder()
{
   // Border destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGPCell                                                              //
//                                                                      //
// GenePix base class describing a single microarray cell               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGPCell::XGPCell()
{
   // Default GenePix Cell constructor
   fMedian = -1;
   fMean   = -1;
   fStdev  = -1;
}//Constructor

//______________________________________________________________________________
XGPCell::~XGPCell()
{
   // GenePix Cell destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFeature635                                                          //
//                                                                      //
// GenePix Feature635 class for a single microarray cell                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XFeature635::XFeature635()
            :XGPCell()
{
   // Default Feature635 constructor
   fPix1SD     = -1;
   fPix2SD     = -1;
   fSaturation = -1;
}//Constructor

//______________________________________________________________________________
XFeature635::~XFeature635()
{
   // Feature635 destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFeature532                                                          //
//                                                                      //
// GenePix Feature532 class for a single microarray cell                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XFeature532::XFeature532()
            :XFeature635()
{
   // Default Feature532 constructor
   fX        = -1;
   fY        = -1;
   fDiameter = -1;
   fNPixels  = -1;
}//Constructor

//______________________________________________________________________________
XFeature532::~XFeature532()
{
   // Feature532 destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBg532                                                               //
//                                                                      //
// GenePix Background532 class for a single microarray cell             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBg532::XBg532()
       :XGPCell()
{
   // Default Background532 constructor
   fNPixels = -1;
}//Constructor

//______________________________________________________________________________
XBg532::~XBg532()
{
   // Background532 destructor 
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGPRatio                                                             //
//                                                                      //
// GenePix class describing ratios of intensities a single cell         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGPRatio::XGPRatio()
{
   // Default GenePix Ratio constructor
   fRatioMd  = 0;
   fRatioMn  = 0;
   fMdRatio  = 0;
   fMnRatio  = 0;
   fRatioSD  = 0;
   fRgnRatio = 0;
   fRgnR2    = 0;
   fLogRatio = 0;
}//Constructor

//______________________________________________________________________________
XGPRatio::~XGPRatio()
{
   // GenePix Ratio destructor 
}//Destructor







