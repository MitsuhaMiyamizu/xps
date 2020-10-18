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

#ifndef __XPSProcessing__
#define __XPSProcessing__

#ifndef __XPSSchemes__
#include "XPSSchemes.h"
#endif

#ifndef __XPSData__
#include "XPSData.h"
#endif

#ifndef __XProjectHandler__
#include "XPSProjectHandler.h"
#endif

// tree preprocessing extensions and types
extern const char *kExtenBgrd[];
extern const char *kTypeBgrd[];
extern const char *kExtenIntn[];
extern const char *kTypeIntn[];
extern const char *kExtenResd[];
extern const char *kTypeResd[];
extern const char *kExtenBord[];
extern const char *kTypeBord[];
extern const char *kExtenCNrm[];
extern const char *kTypeCNrm[];
extern const char *kExtenExpr[];
extern const char *kTypeExpr[];
extern const char *kExtenCall[];
extern const char *kTypeCall[];
extern const char *kExtenQual[];
extern const char *kTypeQual[];

// tree normation extensions and types
extern const char *kExtenNorm[];
extern const char *kTypeNorm[];
extern const char *kExtenSlct[];
extern const char *kTypeSlct[];

// tree analysis extensions and types
extern const char *kExtenFltr[];
extern const char *kTypeFltr[];
extern const char *kPreFltr[];
extern const char *kUniFltr[];
extern const char *kMultiFltr[];
extern const char *kExtenUTst[];
extern const char *kTypeUTst[];
extern const char *kExtenMTst[];
extern const char *kTypeMTst[];
extern const char *kExtenClus[];
extern const char *kTypeClus[];
extern const char *kExtenRgrs[];
extern const char *kTypeRgrs[];

// processing methods
extern const char *kProcessMethod[];

// options
extern const char *kCallOption[];
extern const char *kQualOption[];

// name of reference tree
extern const char *kReference;

// Probeset option
enum EProbesetOption {
   ePROBESET   = 2,
   eEXONTYPE   = 1,
   eTRANSCRIPT = 0, 
};

// Weight
enum EWeight {
   eINITWEIGHT = -16384, 
};

// Error messages for XPSProcessing
enum EErrorProcessing {
   errGetScheme     = -301,
   errSchemeDerived = -302,
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProcessManager                                                      //
//                                                                      //
// Base class for processing of microarray data                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProcessManager: public XManager, public XProjectHandler {

   protected:
      TFile          *fSchemeFile;    //! file containing chip definitions
      XFolder        *fSchemes;       //! content of scheme file
      TFile          *fDataFile;      //! currently open data file
      XFolder        *fData;          //! content of current data file
      Bool_t          fIsSchemeOwner; //! TRUE if fSchemeFile opened by this
      Bool_t          fIsDataOwner;   //! TRUE if fDataFile opened by this

   protected:
      virtual XSetting *NewSetting(const char *type, const char *infile);
      virtual XPlot    *NewPlotter(const char *name, const char *title = "");
      virtual Int_t     DeleteTreeSetInfo(const char *name);

   public:
      XProcessManager();
      XProcessManager(const char *name, const char *arraytype = "",
                      Int_t verbose = kTRUE);
      virtual ~XProcessManager();

      virtual Int_t AddTree(const char *setname, const char *intree,
                       Int_t treeid = 1, Option_t *option = "");
      virtual void  Close(Option_t *option = "");
      virtual Int_t HandleError(Int_t err, const char *name1 = "",
                       const char *name2 = "");

      Int_t SetBaseLine(const char *intree, Option_t *option = "mean",
               Double_t trim = 0.0);
      Int_t SetReference(const char *intree, Option_t *option = "mean",
               Double_t trim = 0.0);

      Int_t InitData(TFile *datafile, Bool_t isOwner = kFALSE);
      Int_t InitSchemes(TFile *schemefile, Bool_t isOwner = kFALSE,
               const char *schemename = "", const char *schemetype = "");
      Int_t OpenData(const char *fullname, Option_t *option = "READ");
      Int_t OpenSchemes(const char *fullname, const char *schemename = "",
               const char *schemetype = "");
      void  CloseData();
      void  CloseSchemes();

      // methods to support database
      virtual Int_t BeginTransaction(const char *name = "");
      virtual Int_t CommitTransaction();

      TFile    *GetDataFile()   const {return fDataFile;}
      TFile    *GetSchemeFile() const {return fSchemeFile;}
      XFolder  *GetData()       const {return fData;}
      XFolder  *GetSchemes()    const {return fSchemes;}
      
      ClassDef(XProcessManager,2) //ProcessManager
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExpressionTreeInfo                                                  //
//                                                                      //
// Class containing info about expression tree stored in fUserInfo      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExpressionTreeInfo: public XTreeInfo {

   protected:
      Int_t      fNUnits;       //number of units
      Double_t   fMinLevel;     //minimal expression level
      Double_t   fMaxLevel;     //maximal expression level
      Int_t      fNQuantiles;   //number of quantiles
      Double_t  *fQuantiles;    //[fNQuantiles] Array of quantile values
      Double_t  *fLevelQuant;   //[fNQuantiles] Array of level quantiles

   public :
      XExpressionTreeInfo();
      XExpressionTreeInfo(const char *name, const char *title);
      virtual ~XExpressionTreeInfo();

      using XTreeInfo::AddUserInfo;
      virtual void     AddUserInfo(Int_t nunits, Double_t min, Double_t max);
      virtual void     AddUserInfo(Int_t nquant, Double_t *q, Double_t *quant);
      virtual Double_t GetValue(const char *name);

      Double_t *GetQuantiles();
      Double_t *GetLevelQuantiles();

      ClassDef(XExpressionTreeInfo,2) //ExpressionTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSelectionTreeInfo                                                   //
//                                                                      //
// Class containing info about selection mask tree stored in fUserInfo  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSelectionTreeInfo: public XTreeInfo {

   protected:
      Int_t      fNUnits;       //number of units
      Int_t      fNFlags;       //number of positive flags

   public :
      XSelectionTreeInfo();
      XSelectionTreeInfo(const char *name, const char *title);
      virtual ~XSelectionTreeInfo();

      using XTreeInfo::AddUserInfo;
      virtual void     AddUserInfo(Int_t nunits, Int_t nflags);
      virtual Double_t GetValue(const char *name);

      ClassDef(XSelectionTreeInfo,1) //SelectionTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProcesSetting                                                       //
//                                                                      //
// Class for initialization of input settings                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProcesSetting: public XSetting {

   protected:
      TFile      *fDataFile;     //! file containing chip data
      TFile      *fSchemeFile;   //! file containing chip definitions
      TString     fSchemeName;   //! name of chip 
      TString     fSchemeType;   //! chip type

   public:
      XProcesSetting();
      XProcesSetting(const char *arraytype, const char *infile);
      virtual ~XProcesSetting();

      void    SetDataFile(TFile *file)    {fDataFile   = file;}
      void    SetSchemeFile(TFile *file)  {fSchemeFile = file;}
      void    SetSchemeName(TString chip) {fSchemeName = chip;}
      void    SetSchemeType(TString type) {fSchemeType = type;}
      TFile  *GetDataFile()         const {return fDataFile;}
      TFile  *GetSchemeFile()       const {return fSchemeFile;}
      TString GetSchemeName()       const {return fSchemeName;}
      TString GetSchemeType()       const {return fSchemeType;}

      ClassDef(XProcesSetting,1) //ProcesSetting
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProcesSet                                                           //
//                                                                      //
// Base class containing info about tree sets used for processing of    //
// microarray data stored as trees in root files                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProcesSet: public XTreeSet {

   protected:
      TFile        *fSchemeFile;  //! file containing chip definitions
      XFolder      *fSchemes;     //! content of scheme file
      TFile        *fDataFile;    //! currently open data file
      XFolder      *fData;        //! content of current data file
//old      TString       fSchemeName;  //! name of chip
      TString       fSchemeName;  //name of chip
      TList        *fBaselines;   //! list of selected baseline trees
      TList        *fReferences;  //! list of selected reference trees
//??      TList        *fBaselines;   //list of selected baseline trees
//??      TList        *fReferences;  //list of selected reference trees
      TString       fBaseOpt;     //baseline option: none, mean or median
      TString       fRefOpt;      //reference option: none, mean or median
      Double_t      fBaseTrim;    //trim value if fBaseOpt is mean
      Double_t      fRefTrim;     //trim value if fRefOpt is mean

   protected:
      virtual void  AddExprTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nunits, Double_t min, Double_t max, 
                       Int_t nquant, Double_t *q, Double_t *quant);
      virtual Int_t InitGroups(Int_t &n, Int_t *gid, TTree **tree, const char **extens);

      virtual Int_t ExportExprTreeInfo(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportMaskTreeInfo(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);

      virtual const char *GetTranscriptID(XTransAnnotation *anno);
      virtual const char *GetTranscriptID(XUnit *unit, XTransAnnotation *annot, Int_t type);

      Int_t  CopyUnitBranch(TTree *fromtree, TTree *totree, Int_t writeopt = -1);
      Int_t  ExpressionQuantiles(TTree *exprtree, XExpression *expr,
                Int_t nquant, Double_t *q, Double_t *quantL);
      TTree *GetUnitTree(XGeneChip *chip, Int_t type);
      TTree *GetAnnotationTree(XGeneChip *chip, Int_t type);

      THashTable *FillHashTable(THashTable *htable, TTree *anntree,
                                XTransAnnotation *annot);
      THashTable *FillHashTable(THashTable *htable, TTree *anntree,
                                XTransAnnotation *annot, Int_t type);
      XIdxString *FindUnitID(THashTable *htable, XGCUnit *unit);
      XIdxString *FindUnitID(THashTable *htable, XExonUnit *unit);

   public:
      XProcesSet();
      XProcesSet(const char *name, const char *type);
      virtual ~XProcesSet();

      virtual Int_t InitFiles(TFile *schemefile, TFile *datafile);

      using XTreeSet::HandleOption;
      virtual Int_t HandleOption(TTree *tree, Option_t *opt);

      Int_t SetBaseLine(TTree *tree, Option_t *option = "", Double_t trim = 0.0);
      Int_t SetReference(TTree *tree, Option_t *option = "", Double_t trim = 0.0);

      void    SetSchemeName(const char *name) {fSchemeName = name;}
      TString GetSchemeName()           const {return fSchemeName;}

      ClassDef(XProcesSet,2) //ProcesSet
};

#endif

