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

#ifndef __XPSData__
#define __XPSData__


#ifndef __XPSSchemes__
#include "XPSSchemes.h"
#endif

#ifndef __XProjectHandler__
#include "XPSProjectHandler.h"
#endif

// tree data extensions
extern const char *kExtenData[];
extern const char *kTypeData[];
extern const char *kExtenMask[];
extern const char *kTypeMask[];

// synonyms for variables
extern const Int_t kNSynExpr;
extern const char *kSynExpr[];
extern const Int_t kNSynCall;
extern const char *kSynCall[];
extern const Int_t kNSynPVal;
extern const char *kSynPVal[];

// Error messages for XPSData
enum EErrorData {
   errChipType   = -201,
   errCELVersion = -202,
   errNameValue  = -203,
   errNumCells   = -204,
};

class XHybridization;
class XGCCell;
class XPosition;
class XPlot;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataManager                                                         //
//                                                                      //
// Class for managing import of microarray data                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDataManager: public XManager, public XProjectHandler {

   private:
      TFile          *fSchemeFile;    //! file containing chip definitions
      XFolder        *fSchemes;       //! content of scheme file
      Bool_t          fIsSchemeOwner; //! delete fSchemeFile if TRUE

   protected:
      virtual Int_t     ImportDefaults(const char *infile);
      virtual Int_t     InitDefaults();
      virtual XSetting *NewSetting(const char *type, const char *infile);
      virtual XTreeSet *NewTreeSet(const char *type); 
      virtual XPlot    *NewPlotter(const char *name, const char *title = "");
      virtual Int_t     DeleteTreeSetInfo(const char *name);

   public:
      XDataManager();
      XDataManager(const char *name, const char *title = "", Int_t verbose = kTRUE);
      virtual ~XDataManager();

      Int_t InitInput(const char *schemename, const char *datatype,
               const char *varlist, const char *inputtype = "rawdata");

      Int_t InitSchemes(TFile *schemefile, Bool_t isOwner = kFALSE,
               const char *schemename = "");
      Int_t OpenSchemes(const char *fullname, const char *schemename = "");
      void  CloseSchemes();

      virtual Int_t New(const char *name, const char *dir, const char *type,
                        const char *data = "Data", Option_t *option = "");
      virtual void  Close(Option_t *option = "");

      virtual Int_t Import(const char *setname, const char *infile, 
                       const char *treename = "", Option_t *option = "NEW",
                       const char *sep = "\t", char delim = '\n', Int_t split = 99);
      virtual Int_t HandleError(Int_t err, const char *name1 = "",
                       const char *name2 = "");

      virtual Int_t DrawUnit(const char *canvasname, const char *treename,
                       const char *schemename, const char *unitname, Int_t unitID,
                       const char *varlist, const char *logbase, const char *type,
                       Option_t *opt, Double_t min = -1111, Double_t max = -1111,
                       const char *var2sort = "");

      // methods to support database
      virtual Int_t BeginTransaction(const char *name = "");
      virtual Int_t CommitTransaction();

      TFile    *GetSchemeFile() const {return fSchemeFile;}
      XFolder  *GetSchemes()    const {return fSchemes;}
      
      ClassDef(XDataManager,2) //DataManager
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTreeInfo                                                        //
//                                                                      //
// Class containing info about data tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDataTreeInfo: public XTreeInfo {

   protected:
      Int_t      fNRows;        //number of rows
      Int_t      fNCols;        //number of columns
      Double_t   fMinInten;     //minimal intensity of hybridized chip
      Double_t   fMaxInten;     //maximal intensity of hybridized chip
      Int_t      fNMinInten;    //number of cells with minimal intensity
      Int_t      fNMaxInten;    //number of cells with minimal intensity
      Short_t    fMaxNPixels;   //maximal pixel number
      Int_t      fNQuantiles;   //number of quantiles
      Double_t  *fQuantiles;    //[fNQuantiles] Array of quantile values
      Double_t  *fIntenQuant;   //[fNQuantiles] Array of intensity quantiles

   public :
      XDataTreeInfo();
      XDataTreeInfo(const char *name, const char *title);
      virtual ~XDataTreeInfo();

      virtual void     AddUserInfo(XTreeSet *set);
      virtual void     AddUserInfo(Int_t nrows, Int_t ncols, Int_t nmin,
                          Double_t min, Int_t nmax, Double_t max, Short_t maxnpix);
      virtual void     AddUserInfo(Int_t nquant, Double_t *q, Double_t *quant);
      virtual Double_t GetValue(const char *name);

      Double_t *GetQuantiles();
      Double_t *GetIntenQuantiles();

      ClassDef(XDataTreeInfo,2) //DataTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMaskTreeInfo                                                        //
//                                                                      //
// Class containing info about mask tree stored in fUserInfo            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMaskTreeInfo: public XTreeInfo {

   protected:
      Int_t      fNRows;        //number of rows
      Int_t      fNCols;        //number of columns
      Int_t      fNFlags;       //number of positive flags

   public :
      XMaskTreeInfo();
      XMaskTreeInfo(const char *name, const char *title);
      virtual ~XMaskTreeInfo();

      using XTreeInfo::AddUserInfo;
      virtual void     AddUserInfo(Int_t nrows, Int_t ncols, Int_t nflags);
      virtual Double_t GetValue(const char *name);

      ClassDef(XMaskTreeInfo,1) //MaskTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataSetting                                                         //
//                                                                      //
// Class for initialization of input data settings                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDataSetting: public XSetting {

   protected:
      TFile     *fSchemeFile;   //! file containing chip definitions
      TString    fSchemeName;   //name of chip type, i.e. scheme
      TString    fDataType;     //type of input data, e.g. cel,adf,tbw
      TString    fInputType;    //type of input file, e.g. PivotData,Metrics
      Int_t      fNVar;         //number of variables to read from input
      TString   *fVars;         //[fNVar] array of fNVar variables
      TString   *fTyps;         //[fNVar] variable type of fVars: C,I,F,D

   public:
      XDataSetting();
      XDataSetting(const char *arraytype, const char *infile);
      virtual ~XDataSetting();

      virtual Int_t InitInput(const char *schemename, const char *datatype,
                       const char *varlist, const char *inputtype);

      void SetSchemeFile(TFile *file)  {fSchemeFile = file;}
      void SetSchemeName(TString chip) {fSchemeName = chip;}

      TFile   *GetSchemeFile()   const {return fSchemeFile;}
      TString  GetSchemeName()   const {return fSchemeName;}
      TString  GetSchemeType()   const {return this->GetName();}
      TString  GetDataType()     const {return fDataType;}
      TString  GetInputType()    const {return fInputType;}
      Int_t    GetNumVariables() const {return fNVar;}
      TString *GetVariables()    const {return fVars;}
      TString *GetVarTypes()     const {return fTyps;}

      ClassDef(XDataSetting,1) //DataSetting
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybridization                                                       //
//                                                                      //
// Base class for microarray hybridization                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XHybridization: public XTreeSet {

   protected:
      TFile     *fSchemeFile;   //! file containing chip definitions
      TTree     *fDataTree;     //! data tree ??del??
      TTree     *fMaskTree;     //! mask tree for data tree ??del??
      TString    fDataName;     //name of input data file ??del??
      TString    fSchemeName;   //name of chip type
      TString    fDataTreeName; //! name of tree containing raw data  ??del??
      TString    fMaskTreeName; //! name of tree containing mask ??del??
      Int_t      fNRows;        //number of rows of hybridization file
      Int_t      fNCols;        //number of columns of hybridization file
      Int_t      fNCells;       //number of cells of hybridization file

   protected:
      Int_t Index2X(Int_t index) {return (index % fNCols);}
      Int_t Index2Y(Int_t index) {return (index / fNCols);}

   public:
      XHybridization();
      XHybridization(const char *name, const char *title);
      virtual ~XHybridization();

      virtual Int_t InitFiles(TFile *schemefile);
      virtual void  PrintInfo() = 0;

      void SetDataName(TString name)     {fDataName = Path2Name(name,"/","");}
      void SetSchemeName(TString name)   {fSchemeName = name;}
      void SetDataTreeName(TString name) {fDataTreeName = name;}
      void SetMaskTreeName(TString name) {fMaskTreeName = name;}
      TString GetDataName()        const {return fDataName;}
      TString GetSchemeName()      const {return fSchemeName;}
      TString GetDataTreeName()    const {return fDataTreeName;}
      TString GetMaskTreeName()    const {return fMaskTreeName;}
      Int_t   GetNumRows()         const {return fNRows;}
      Int_t   GetNumColumns()      const {return fNCols;}
      Int_t   GetNumCells()        const {return fNCells;}

      ClassDef(XHybridization,1) //Hybridization
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipHyb                                                         //
//                                                                      //
// Class for GeneChip hybridization                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGeneChipHyb: public XHybridization {

   protected:
      Double_t   fMinInten;     //! minimal intensity of hybridized chip
      Double_t   fMaxInten;     //! maximal intensity of hybridized chip
      Int_t      fNMinInten;    //! number of cells with minimal intensity
      Int_t      fNMaxInten;    //! number of cells with minimal intensity
      Short_t    fMaxNPixels;   //! maximal pixel number

   protected:
      virtual void    AddDataTreeInfo(TTree *tree, const char *name, Option_t *option,
                         Int_t nquant, Double_t *q, Double_t *quant);

      virtual Int_t   ExportDataTrees(Int_t n, TString *names, const char *varlist,
                         ofstream &output, const char *sep);
      virtual Int_t   ExportDataTree(TString *name, ofstream &output, const char *sep);
      virtual Int_t   ExportMaskTrees(Int_t n, TString *names, const char *varlist,
                         ofstream &output, const char *sep);
      virtual Int_t   ExportMaskTree(TString *name, ofstream &output, const char *sep);

      virtual Int_t   ExportDataTreeInfo(Int_t n, TString *names, const char *varlist,
                         ofstream &output, const char *sep);
      virtual Int_t   ExportMaskTreeInfo(Int_t n, TString *names, const char *varlist,
                         ofstream &output, const char *sep);

      virtual Int_t   IsXDAFile(ifstream &input);
      virtual Int_t   IsCalvinFile(ifstream &input);
      virtual TString ChipType(const char *header, Int_t toUpper = 0);
      virtual Int_t   CheckChipType(const char *header, const char *name);

      virtual Int_t   ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t   ReadData(ifstream &input, Option_t *option, const char *sep,
                         char delim, Int_t split);
      virtual Int_t   ReadXDAHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t   ReadXDAData(ifstream &input, Option_t *option,
                         const char *sep, char delim, Int_t split);
      virtual Int_t   ReadFileHeader(ifstream &input, Int_t &numgroups, UInt_t &filepos);
      virtual Int_t   ReadGenericDataHeader(ifstream &input, Bool_t isParent);
      virtual Int_t   ReadDataGroup(ifstream &input, UInt_t &filepos, Option_t *option,
                         Int_t split);
      virtual Int_t   ReadCalvinGenericFile(ifstream &input, Option_t *option, Int_t split);
      virtual Int_t   ReadXMLHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t   ReadXMLData(ifstream &input, Option_t *option,
                         const char *sep, char delim, Int_t split);

      Int_t  DataQuantiles(TTree *tree, XGCCell *cell, Int_t nquant, Double_t *q, Double_t *quant);

   public:
      XGeneChipHyb();
      XGeneChipHyb(const char *name, const char *title);
      virtual ~XGeneChipHyb();

      virtual Int_t ExportTreeInfo(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t Import(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual void  PrintInfo();

      Double_t GetMinIntensity()    const {return fMinInten;}
      Double_t GetMaxIntensity()    const {return fMaxInten;}
      Int_t    GetNumMinIntensity() const {return fNMinInten;}
      Int_t    GetNumMaxIntensity() const {return fNMaxInten;}
      Short_t  GetMaxNumPixels()    const {return fMaxNPixels;}

      ClassDef(XGeneChipHyb,1) //GeneChip Hybridization
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipMetrics                                                     //
//                                                                      //
// Class for processed GeneChip CHP data exported as Metrics.txt        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGeneChipMetrics: public XHybridization {

   protected:

   protected:
      virtual Int_t ExportExprTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportCallTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      virtual Int_t ReadXMLHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadXMLData(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);

   public:
      XGeneChipMetrics();
      XGeneChipMetrics(const char *name, const char *title);
      virtual ~XGeneChipMetrics();

      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual void  PrintInfo();

      ClassDef(XGeneChipMetrics,1) //GeneChip Metrics
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipPivot                                                       //
//                                                                      //
// Class for processed GeneChip CHP data exported as PivotData.txt      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGeneChipPivot: public XHybridization {

   protected:
      Int_t      fNTrees;       //number of trees to read from input
      TString   *fTreeName;     //[fNTrees] array of fNTrees tree names

   protected:
      virtual Int_t ExportExprTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportCallTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      virtual Int_t ReadXMLHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadXMLData(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);

   public:
      XGeneChipPivot();
      XGeneChipPivot(const char *name, const char *title);
      virtual ~XGeneChipPivot();

      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual void  PrintInfo();

      ClassDef(XGeneChipPivot,1) //GeneChip Pivot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChipHyb                                                          //
//                                                                      //
// Class for Mapping array hybridization                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSNPChipHyb: public XGeneChipHyb {

   protected:

   public:
      XSNPChipHyb();
      XSNPChipHyb(const char *name, const char *title);
      virtual ~XSNPChipHyb();

      ClassDef(XSNPChipHyb,1) //SNPChip Hybridization
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChipMetrics                                                      //
//                                                                      //
// Class for processed Mapping array CHP data exported as Metrics.txt   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSNPChipMetrics: public XGeneChipMetrics {

   protected:

   public:
      XSNPChipMetrics();
      XSNPChipMetrics(const char *name, const char *title);
      virtual ~XSNPChipMetrics();

      ClassDef(XSNPChipMetrics,1) //SNPChip Metrics
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChipPivot                                                        //
//                                                                      //
// Class for processed Mapping array CHP data exported as PivotData.txt //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSNPChipPivot: public XGeneChipPivot {

   protected:

   public:
      XSNPChipPivot();
      XSNPChipPivot(const char *name, const char *title);
      virtual ~XSNPChipPivot();

      ClassDef(XSNPChipPivot,1) //SNPChip Pivot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChipHyb                                                       //
//                                                                      //
// Class for whole genome array hybridization                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeChipHyb: public XGeneChipHyb {

   protected:

   public:
      XGenomeChipHyb();
      XGenomeChipHyb(const char *name, const char *title);
      virtual ~XGenomeChipHyb();

      ClassDef(XGenomeChipHyb,1) //GenomeChip Hybridization
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChipMetrics                                                   //
//                                                                      //
// Class for processed genome array data exported as Metrics.txt        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeChipMetrics: public XGeneChipMetrics {

   protected:

   public:
      XGenomeChipMetrics();
      XGenomeChipMetrics(const char *name, const char *title);
      virtual ~XGenomeChipMetrics();

      ClassDef(XGenomeChipMetrics,1) //GenomeChip Metrics
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChipPivot                                                     //
//                                                                      //
// Class for processed genome array data exported as PivotData.txt      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeChipPivot: public XGeneChipPivot {

   protected:

   public:
      XGenomeChipPivot();
      XGenomeChipPivot(const char *name, const char *title);
      virtual ~XGenomeChipPivot();

      ClassDef(XGenomeChipPivot,1) //GenomeChip Pivot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChipHyb                                                         //
//                                                                      //
// Class for Exon array hybridization                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonChipHyb: public XGeneChipHyb {

   protected:

   public:
      XExonChipHyb();
      XExonChipHyb(const char *name, const char *title);
      virtual ~XExonChipHyb();

      ClassDef(XExonChipHyb,1) //ExonChipHyb Hybridization
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChipMetrics                                                     //
//                                                                      //
// Class for processed exon array data exported as Metrics.txt          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonChipMetrics: public XGeneChipMetrics {

   protected:

   public:
      XExonChipMetrics();
      XExonChipMetrics(const char *name, const char *title);
      virtual ~XExonChipMetrics();

      ClassDef(XExonChipMetrics,1) //ExonChip Metrics
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChipPivot                                                       //
//                                                                      //
// Class for processed exon array data exported as PivotData.txt        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonChipPivot: public XGeneChipPivot {

   protected:

   public:
      XExonChipPivot();
      XExonChipPivot(const char *name, const char *title);
      virtual ~XExonChipPivot();

      ClassDef(XExonChipPivot,1) //ExonChip Pivot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenePixHyb                                                          //
//                                                                      //
// Class for GenePix microarray hybridization                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenePixHyb: public XHybridization {

   protected:
      Int_t      fNRows;        //number of rows  ??
      Int_t      fNCols;        //number of columns ??
      Double_t   fMinInten;     //minimal intensity of hybridized chip
      Double_t   fMaxInten;     //maximal intensity of hybridized chip
      Int_t      fNMinInten;    //number of cells with minimal intensity
      Int_t      fNMaxInten;    //number of cells with minimal intensity

   protected:
      virtual void  AddDataTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nquant, Double_t *q, Double_t *quant);
      virtual Int_t ExportDataTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportDataTree(ofstream &output, const char *sep);
      virtual Int_t ExportMaskTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportMaskTree(ofstream &output, const char *sep);
      virtual Int_t ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      virtual Int_t ReadXMLHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadXMLData(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);

   public:
      XGenePixHyb();
      XGenePixHyb(const char *name, const char *title);
      virtual ~XGenePixHyb();

      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual void  PrintInfo();

      Int_t    GetNumRows()         const {return fNRows;}
      Int_t    GetNumColumns()      const {return fNCols;}
      Double_t GetMinIntensity()    const {return fMinInten;}
      Double_t GetMaxIntensity()    const {return fMaxInten;}
      Int_t    GetNumMinIntensity() const {return fNMinInten;}
      Int_t    GetNumMaxIntensity() const {return fNMaxInten;}

      ClassDef(XGenePixHyb,1) //GenePix Hybridization
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalysisPlot                                                        //
//                                                                      //
// Class for drawing microarray results                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAnalysisPlot: public XPlot {

   protected:
      TFile       *fSchemeFile;     //! file containing chip definitions
      XFolder     *fSchemes;        //! content of scheme file
      TFile       *fDataFile;       //! currently open data file
      XFolder     *fData;           //! content of current data file
      Bool_t       fIsSchemeOwner;  //! TRUE if owner of fSchemeFile
      Bool_t       fIsDataOwner;    //! TRUE if owner of fDataFile

   protected:
//??      Int_t FindUnitID(const char *unitname);

   public :
      XAnalysisPlot();
      XAnalysisPlot(const char *name, const char *title);
      virtual ~XAnalysisPlot();

      virtual Int_t DrawUnit(const char * /*canvasname*/, const char * /*treename*/,
                       const char * /*schemename*/, const char * /*unitname*/, Int_t /*unitID*/,
                       const char * /*varlist*/, const char * /*logbase*/, const char * /*type*/,
                       Option_t * /*opt*/, Double_t min = -1111, Double_t max = -1111,
                       const char *var2sort = "", Bool_t down = kFALSE) {return 0;}

      Int_t InitData(TFile *datafile, Bool_t isOwner = kFALSE);
      Int_t InitSchemes(TFile *schemefile, Bool_t isOwner = kFALSE);
      Int_t OpenData(const char *fullname);
      Int_t OpenSchemes(const char *fullname);
      void  CloseData();
      void  CloseSchemes();

      TFile *GetDataFile()   const {return fDataFile;}
      TFile *GetSchemeFile() const {return fSchemeFile;}
      void   SetDataFile(TFile *file, XFolder *data)
                {fDataFile = file; fData = data; fIsDataOwner = kFALSE;}
      void   SetSchemeFile(TFile *file, XFolder *schemes)
                {fSchemeFile = file; fSchemes = schemes; fIsSchemeOwner = kFALSE;}

      ClassDef(XAnalysisPlot,1) //AnalysisPlot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChipPlot                                                        //
//                                                                      //
// Class for drawing GeneChip data                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGeneChipPlot: public XAnalysisPlot {

   protected:

   protected:

   public:
      XGeneChipPlot();
      XGeneChipPlot(const char *name, const char *title = "");
      virtual ~XGeneChipPlot();

      virtual Int_t DrawUnit(const char *canvasname, const char *treename,
                       const char *schemename, const char *unitname, Int_t unitID,
                       const char *varlist, const char *logbase, const char *type,
                       Option_t *opt, Double_t min = -1111, Double_t max = -1111,
                       const char *var2sort = "", Bool_t down = kFALSE);

      ClassDef(XGeneChipPlot,1) //GeneChipPlot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenePixPlot                                                         //
//                                                                      //
// Class for drawing GenePix microarray data                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenePixPlot: public XAnalysisPlot {

   protected:

   protected:

   public :
      XGenePixPlot();
      XGenePixPlot(const char *name, const char *title);
      virtual ~XGenePixPlot();

      virtual Int_t DrawUnit(const char *canvasname, const char *treename,
                       const char *schemename, const char *unitname, Int_t unitID,
                       const char *varlist, const char *logbase, const char *type,
                       Option_t *opt, Double_t min = -1111, Double_t max = -1111,
                       const char *var2sort = "", Bool_t down = kFALSE);

      ClassDef(XGenePixPlot,1) //GenePixPlot
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCellOM                                                              //
//                                                                      //
// Class containing outlier flag and mask flag of a single  cell        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCellOM {

   private:
      Bool_t     fIsOutlier;    //flagged if cell is outlier
      Bool_t     fIsMasked;     //flagged if cell is masked

   public:
      XCellOM() {}
      virtual ~XCellOM() {}

      void   SetOutlier(Bool_t outlier) {fIsOutlier = outlier;}
      void   SetMasked(Bool_t masked)   {fIsMasked  = masked;}

      Bool_t GetOutlier()         const {return fIsOutlier;}
      Bool_t GetMasked()          const {return fIsMasked;}

      ClassDef(XCellOM,1) //CellOM
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMask                                                                //
//                                                                      //
// Base class for the mask of a single microarray cell or unit          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMask {

   protected:
      Short_t   fFlag;          //type of flag

   public:
      XMask() {}
      virtual ~XMask() {}

      void    SetFlag(Short_t flag) {fFlag = flag;}
      Short_t GetFlag()       const {return fFlag;}

      ClassDef(XMask,1) //Mask
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCellMask                                                            //
//                                                                      //
// Class for the mask of a single microarray cell                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCellMask: public XPosition {

   protected:
      Short_t   fFlag;          //type of flag

   public:
      XCellMask();
      virtual ~XCellMask();

      void    SetFlag(Short_t flag) {fFlag = flag;}
      Short_t GetFlag()       const {return fFlag;}

      ClassDef(XCellMask,1) //CellMask
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnitMask                                                            //
//                                                                      //
// Class for the mask of a single microarray unit                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUnitMask: public XMask {

   protected:
      Int_t      fUnitID;      //unit-ID

   public:
      XUnitMask() {}
      virtual ~XUnitMask() {}

      void  SetUnitID(Int_t unitID) {fUnitID = unitID;}
      Int_t GetUnitID()       const {return fUnitID;}

      ClassDef(XUnitMask,1) //UnitMask
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCell                                                                //
//                                                                      //
// Base class describing a single microarray cell                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCell {

   protected:
      Double32_t fInten;        //cell intensity
      Double32_t fStdev;        //standard deviation
      Short_t    fNPixels;      //number of pixels used for calculation

   public:
      XCell();
      virtual ~XCell();

      void     SetIntensity(Double_t inten) {fInten   = inten;}
      void     SetStdev(Double_t stdev)     {fStdev   = stdev;}
      void     SetNumPixels(Short_t numpix) {fNPixels = numpix;}

      Double_t GetIntensity()         const {return fInten;}
      Double_t GetStdev()             const {return fStdev;}
      Short_t  GetNumPixels()         const {return fNPixels;}

      ClassDef(XCell,1) //Cell
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCCell                                                              //
//                                                                      //
// Class describing a single GeneChip cell                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCCell: public XPosition {

   protected:
      Double32_t fInten;        //cell intensity
      Double32_t fStdev;        //standard deviation
      Short_t    fNPixels;      //number of pixels used for calculation

   public:
      XGCCell();
      virtual ~XGCCell();

      void     SetIntensity(Double_t inten) {fInten   = inten;}
      void     SetStdev(Double_t stdev)     {fStdev   = stdev;}
      void     SetNumPixels(Short_t numpix) {fNPixels = numpix;}

      Double_t GetIntensity()         const {return fInten;}
      Double_t GetStdev()             const {return fStdev;}
      Short_t  GetNumPixels()         const {return fNPixels;}

      ClassDef(XGCCell,1) //GCCell
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBgCell                                                              //
//                                                                      //
// Class containing background data for microarray cell                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XBgCell: public XPosition {

   protected:
      Double32_t fBg;           //cell background
      Double32_t fStdev;        //standard deviation

   public:
      XBgCell();
      virtual ~XBgCell();

      void     SetBackground(Double_t bg) {fBg    = bg;}
      void     SetStdev(Double_t stdev)   {fStdev = stdev;}

      Double_t GetBackground()      const {return fBg;}
      Double_t GetStdev()           const {return fStdev;}

      ClassDef(XBgCell,1) //BgCell
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRankCell                                                            //
//                                                                      //
// Class containing intensity and its rank for microarray cell          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XRankCell: public XPosition {

   protected:
      Double32_t fInten;        //cell intensity
      Int_t      fRank;         //rank of intensity

   public:
      XRankCell();
      virtual ~XRankCell();

      void     SetIntensity(Double_t inten) {fInten = inten;}
      void     SetRank(Int_t rank)          {fRank  = rank;}

      Double_t GetIntensity()         const {return fInten;}
      Int_t    GetRank()              const {return fRank;}

      ClassDef(XRankCell,1) //RankCell
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XResidual                                                            //
//                                                                      //
// Class describing microarray residuals an weights                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XResidual: public XPosition {

   protected:
      Double32_t fResidual;      //residual
      Double32_t fWeight;        //weight

   public:
      XResidual();
      virtual ~XResidual();

      void     SetResidual(Double_t res) {fResidual = res;}
      void     SetWeight(Double_t w)     {fWeight   = w;}

      Double_t GetResidual()       const {return fResidual;}
      Double_t GetWeight()         const {return fWeight;}

      ClassDef(XResidual,1) //Residual
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBorder                                                              //
//                                                                      //
// Class describing microarrayborder elements                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XBorder: public XPosition {

   protected:
      Double32_t fInten;        //cell intensity
      Short_t    fFlag;         //type of flag

   public:
      XBorder();
      virtual ~XBorder();

      void    SetIntensity(Double_t inten) {fInten = inten;}
      void    SetFlag(Short_t flag)        {fFlag = flag;}

      Double_t GetIntensity()   const {return fInten;}
      Short_t GetFlag()         const {return fFlag;}

      ClassDef(XBorder,1) //Border
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExpression                                                          //
//                                                                      //
// Class containing microarray expression data                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExpression {

   protected:
      Int_t      fUnitID;      //unit-ID
      Double32_t fLevel;       //expression level

   public:
      XExpression() {}
      virtual ~XExpression() {}

      void SetUnitID(Int_t unitID)  {fUnitID = unitID;}
      void SetLevel(Double_t level) {fLevel  = level;}

      Int_t    GetUnitID()    const {return fUnitID;}
      Double_t GetLevel()     const {return fLevel;}

      virtual Double_t GetValue(Int_t id = 1) {Int_t idx = id; id = idx; return fLevel;}

      ClassDef(XExpression,1) //Expression
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCExpression                                                        //
//                                                                      //
// Class containing GeneChip expression data                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCExpression: public XExpression {

   private:
      Double32_t fStdev;       //standard deviation
      Int_t      fNPairs;      //number of probe pairs used

   public:
      XGCExpression() {}
      virtual ~XGCExpression() {}

      void SetStdev(Double_t stdev)  {fStdev  = stdev;}
      void SetNumPairs(Int_t npairs) {fNPairs = npairs;}

      Double_t GetStdev()      const {return fStdev;}
      Int_t    GetNumPairs()   const {return fNPairs;}

      virtual Double_t GetValue(Int_t id = 1) 
         {return (id == 2) ? fStdev : ((id == 3) ? fNPairs : fLevel);}

      ClassDef(XGCExpression,1) //GCExpression
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XQCExpression                                                        //
//                                                                      //
// Class containing quality control expression data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XQCExpression: public XExpression {

   private:
      Double32_t fStderr;      //standard error
      Double32_t fNUSE;        //normalized unscaled standard error
      Double32_t fRLE;         //relative log expression

   public:
      XQCExpression() {}
      virtual ~XQCExpression() {}

      void SetStdErr(Double_t se) {fStderr = se;}
      void SetNUSE(Double_t nuse) {fNUSE   = nuse;}
      void SetRLE(Double_t rle)   {fRLE    = rle;}

      Double_t GetStdErr()  const {return fStderr;}
      Double_t GetNUSE()    const {return fNUSE;}
      Double_t GetRLE()     const {return fRLE;}

      virtual Double_t GetValue(Int_t id = 1) 
         {return (id == 2) ? fStderr : ((id == 3) ?  fNUSE : (id == 4) ? fRLE : fLevel);}

      ClassDef(XQCExpression,1) //QCExpression
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSpliceExpression                                                    //
//                                                                      //
// Class containing exon GeneChip expression and splice data            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSpliceExpression: public XExpression {

   private:
      Double32_t fStdev;   //standard deviation of single exon
      Double32_t fScore;   //splice score

   public:
      XSpliceExpression() {}
      virtual ~XSpliceExpression() {}

      void SetStdev(Double_t stdev)  {fStdev = stdev;}
      void SetScore(Double_t score)  {fScore = score;}

      Double_t GetStdev()      const {return fStdev;}
      Double_t GetScore()      const {return fScore;}

      ClassDef(XSpliceExpression,1) //SpliceExpression
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCall                                                                //
//                                                                      //
// Class containing microarray present call data                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCall {

   protected:
      Int_t      fUnitID;      //Unit-ID
      Short_t    fCall;        //Present call: 'P'=2, 'M'=1, 'A'=0, NA = -1

   public:
      XCall() {}
      virtual ~XCall() {}

      void SetUnitID(Int_t unitID)   {fUnitID = unitID;}
      void SetCall(Short_t call)     {fCall   = call;}

      Int_t    GetUnitID()     const {return fUnitID;}
      Short_t  GetCall()       const {return fCall;}

      ClassDef(XCall,1) //Call
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPCall                                                               //
//                                                                      //
// Class containing microarray present call data with statistics        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPCall: public XCall {

   private:
      Double_t   fPValue;      //p-value (for accuracy must not be Double32_t!)

   public:
      XPCall() {}
      virtual ~XPCall() {}

      void SetPValue(Double_t pval)  {fPValue = pval;}
      Double_t GetPValue()     const {return fPValue;}

      ClassDef(XPCall,1) //PCall
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGPCell                                                              //
//                                                                      //
// GenePix base class describing a single microarray cell               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGPCell {

   private:
      Int_t      fMedian;       //median cell intensity
      Int_t      fMean;         //mean cell intensity
      Int_t      fStdev;        //standard deviation of cell intensity

   public:
      XGPCell();
      virtual ~XGPCell();

      void     SetMedian(Int_t inten) {fMedian = inten;}
      void     SetMean(Int_t inten)   {fMean   = inten;}
      void     SetStdev(Int_t stdev)  {fStdev  = stdev;}
      Int_t    GetMedian()      const {return fMedian;}
      Int_t    GetMean()        const {return fMean;}
      Int_t    GetStdev()       const {return fStdev;}

      ClassDef(XGPCell,1) //GPCell
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFeature635                                                          //
//                                                                      //
// GenePix Feature635 class for a single microarray cell                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XFeature635: public XGPCell {

   private:
      Int_t      fPix1SD;        //percent of pixels with Inten > Bg + 1SD
      Int_t      fPix2SD;        //percent of pixels with Inten > Bg + 2SD
      Int_t      fSaturation;    //percentage of saturated feature pixels

   public:
      XFeature635();
      virtual ~XFeature635();

      void     SetPix1SD(Int_t pcpix)   {fPix1SD     = pcpix;}
      void     SetPix2SD(Int_t pcpix)   {fPix2SD     = pcpix;}
      void     SetSaturation(Int_t sat) {fSaturation = sat;}
      Int_t    GetPix1SD()        const {return fPix1SD;}
      Int_t    GetPix2SD()        const {return fPix2SD;}
      Int_t    GetSaturation()    const {return fSaturation;}

      ClassDef(XFeature635,1) //Feature635
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFeature532                                                          //
//                                                                      //
// GenePix Feature532 class for a single microarray cell                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XFeature532: public XFeature635 {

   private:
      Int_t      fX;            //x-coordinate of feature cell
      Int_t      fY;            //y-coordinate of feature cell
      Int_t      fDiameter;     //diameter of feature indicator
      Int_t      fNPixels;      //number of pixels used for calculation

   public:
      XFeature532();
      virtual ~XFeature532();

      void     SetX(Int_t x)              {fX        = x;}
      void     SetY(Int_t y)              {fY        = y;}
      void     SetDiameter(Int_t dia)     {fDiameter = dia;}
      void     SetNumPixels(Int_t numpix) {fNPixels  = numpix;}
      Int_t    GetX()               const {return fX;}
      Int_t    GetY()               const {return fY;}
      Int_t    GetDiameter()        const {return fDiameter;}
      Int_t    GetNumPixels()       const {return fNPixels;}

      ClassDef(XFeature532,1) //Feature532
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBg532                                                               //
//                                                                      //
// GenePix Background532 class for a single microarray cell             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XBg532: public XGPCell {

   private:
      Short_t    fNPixels;     //number of pixels used for calculation

   public:
      XBg532();
      virtual ~XBg532();

      void     SetNumPixels(Short_t numpix) {fNPixels = numpix;}
      Short_t  GetNumPixels()         const {return fNPixels;}

      ClassDef(XBg532,1) //Background532
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGPRatio                                                             //
//                                                                      //
// GenePix class describing ratios of intensities a single cell         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGPRatio {

   private:
      Double32_t fRatioMd;      //ratio of median intensities
      Double32_t fRatioMn;      //ratio of mean intensities
      Double32_t fMdRatio;      //median of pix ratios of pix intensities
      Double32_t fMnRatio;      //mean of pix ratios of pix intensities
      Double32_t fRatioSD;      //standard deviation of pix intensity ratios
      Double32_t fRgnRatio;     //regression ratio
      Double32_t fRgnR2;        //coefficient for regression value
      Double32_t fLogRatio;     //log2 of ratio of medians

   public:
      XGPRatio();
      virtual ~XGPRatio();

      void     SetRatioMedian(Int_t ratio) {fRatioMd  = ratio;}
      void     SetRatioMean(Int_t ratio)   {fRatioMn  = ratio;}
      void     SetMedianRatio(Int_t ratio) {fMdRatio  = ratio;}
      void     SetMeanRatio(Int_t ratio)   {fMnRatio  = ratio;}
      void     SetRatioStdev(Int_t stdev)  {fRatioSD  = stdev;}
      void     SetRgnRatio(Int_t ratio)    {fRgnRatio = ratio;}
      void     SetRgnR2(Int_t ratio)       {fRgnR2    = ratio;}
      void     SetLogRatio(Int_t ratio)    {fLogRatio = ratio;}
      Double_t GetRatioMedian()      const {return fRatioMd;}
      Double_t GetRatioMean()        const {return fRatioMn;}
      Double_t GetMedianRatio()      const {return fMdRatio;}
      Double_t GetMeanRatio()        const {return fMnRatio;}
      Double_t GetRatioStdev()       const {return fRatioSD;}
      Double_t GetRgnRatio()         const {return fRgnRatio;}
      Double_t GetRgnR2()            const {return fRgnR2;}
      Double_t GetLogRatio()         const {return fLogRatio;}

      ClassDef(XGPRatio,1) //GPRatio
};

#endif
