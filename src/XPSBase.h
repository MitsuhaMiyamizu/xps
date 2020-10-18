// File created: 05/18/2002                          last modified: 01/24/2011
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

#ifndef __XPSBase__
#define __XPSBase__

#include "TDatime.h"
#include "TFolder.h"
#include "TMath.h"
#include "TObjString.h"

#include "XPSUtils.h"

const Int_t kBufSize = 1024;
const Int_t kBuf4096 = 4096;

extern const char* kContent;

// Error messages: errNoErr must be zero!!
enum EErrorMessage {
   errNoErr          =  0,
   errFatal          = -1,
   errAbort          = -2,
   errGeneral        = -3,
   errInitMemory     = -4,
   errCreateFile     = -5,
   errCreateDir      = -6,
   errCreateTree     = -7,
   errCreateTreeSet  = -8,
   errGetFile        = -9,
   errGetDir         = -10,
   errGetTree        = -11,
   errOpenFile       = -12,
   errWriteObject    = -13,
   errGetTreeSet     = -14,
   errGetTreeInfo    = -15,
   errPrematureEOF   = -16,
   errHeaderLine     = -17,
   errMissingColumn  = -18,
   errMissingLine    = -19,
   errReadingInput   = -20,
   errOpenOutput     = -21,
   errInitContent    = -22,
   errMissingContent = -23,
   errInitSetting    = -24,
   errInitTreeSet    = -25,
   errInitParameters = -26,
   errInitPlotter    = -27,
   errNumTreeEntries = -28,
   errEQTreeEntries  = -29,
   errClassTreeSet   = -30,
   errAlgorithm      = -31,
   errUnknownType    = -32,
};

// Verbose messages
enum EVERBOSE {
   vSTATS = 1,
   vFILES = 2,
   vREADS = 3,
};

class XManager;
class XTreeSet;

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XIdxString                                                           //
//                                                                      //
// Class derived from TObjString but including integer index            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XIdxString: public TObjString {

   protected:
      Int_t    fIndex;      //index where string is stored originally
      
   public :
      XIdxString() {}
      XIdxString(Int_t id, const char *str): TObjString(str), fIndex(id) {}
      virtual ~XIdxString() {}

      Int_t GetIndex() const {return fIndex;}

      ClassDef(XIdxString,1) //IdxString
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLdxString                                                           //
//                                                                      //
// Class derived from TObjString but including long index               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XLdxString: public TObjString {

   protected:
      Long_t   fIndex;      //index where string is stored originally
      
   public :
      XLdxString() {}
      XLdxString(Long_t id, const char *str): TObjString(str), fIndex(id) {}
      virtual ~XLdxString() {}

      Long_t GetIndex() const {return fIndex;}

      ClassDef(XLdxString,1) //LdxString
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFolder                                                              //
//                                                                      //
// Class containing info for database about objects stored in TFile     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XFolder: public TFolder {

   protected:
      TString    fType;      //type of objects in folder, e.g. GeneChip
      Bool_t     fIsPublic;  //true if folder is made public

   public:
      XFolder();
      XFolder(const char *name, const char *title, const char *type = "",
              Bool_t isPublic = kTRUE, Bool_t newlist = kFALSE);
      virtual ~XFolder();

      XFolder *AddFolder(const char *name, const char *title, const char *type = "",
                  TCollection *collection = 0);

      using TFolder::FindFullPathName;
      virtual  const char *FindFullPathName(const char *name) const;

      using TFolder::FindObject;
      virtual  TObject *FindObject(const char *name) const;
      virtual  TObject *FindObject(const char *name, Bool_t ignoreCase) const;
      virtual  TObject *FindObject(const char *name, const char *classname) const;
      virtual  TObject *FindObject(const char *name, const char *title,
                           const char *classname) const;
      virtual  TObject *FindObjectAny(const char *name) const;

      void     MakePublic()              {fIsPublic = kTRUE;}
      void     MakePrivate()             {fIsPublic = kFALSE;}
      void     SetType(const char *type) {fType  = type;}

      TString  GetType()  const {return fType;}
      Bool_t   IsFolder() const {return kTRUE;}  //not virtual!
      Bool_t   IsPublic() const {return fIsPublic;}

      ClassDef(XFolder,1) //Folder
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreeInfo                                                            //
//                                                                      //
// Base class containing info about tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTreeInfo: public TNamed {

   protected:
      TString    fSetName;      //Name of tree set
      TString    fSetClass;     //Class name of tree set
      TString    fOption;       //Optional option
//      Option_t  *fOption;       //Optional option

   public :
      XTreeInfo();
      XTreeInfo(const char *name, const char *title);
      virtual ~XTreeInfo();

      virtual Option_t *GetOption()                      {return fOption.Data();}
      virtual void      AddUserInfo(XTreeSet * /*set*/)  {}
      virtual Double_t  GetValue(const char * /*name*/)  {return 0;}
      virtual TString   GetString(const char * /*name*/) {return "";}

      void    SetTreeSetName(const char *name)  {fSetName  = name;}
      void    SetTreeSetClass(const char *name) {fSetClass = name;}
      void    SetOption(Option_t *option)       {fOption   = option;}
      TString GetTreeSetName()            const {return fSetName;}
      TString GetTreeSetClass()           const {return fSetClass;}
      TString GetClassName()              const {return GetTitle();}

      ClassDef(XTreeInfo,1) //TreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreeHeader                                                          //
//                                                                      //
// Class storing tree information in header                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTreeHeader: public TObjString {

   protected:
      TString    fInfile;       //Name of infile if data were imported
      TString    fType;         //Type of tree
      TDatime    fDatime;       //Date and time when tree is created
      Int_t      fID;           //ID of tree
      Int_t      fNPar;         //number of parameters
      Double_t  *fPars;         //[fNPar] array of fNPar parameters
      
   public :
      XTreeHeader();
      XTreeHeader(const char *str, Int_t treeid = 0);
      virtual ~XTreeHeader();

      void     SetParameters(Int_t npar, Double_t *pars);

      void     SetInfile(const char *infile) {fInfile = infile;}
      void     SetType(const char *type)     {fType   = type;}
      void     SetID(const Int_t id)         {fID     = id;}

      TString   GetInfile()        const {return fInfile;}
      TString   GetType()          const {return fType;}
      Int_t     GetID()            const {return fID;}
      TDatime  &GetCreationDate()        {return fDatime;}
      Int_t     GetNumParameters() const {return fNPar;}
      Double_t *GetParameters()    const {return fPars;}

      ClassDef(XTreeHeader,1) //TreeHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSetting                                                             //
//                                                                      //
// Class for initialization of parameter settings                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSetting: public TNamed {

   protected:
      Double_t    fNA;           //value to be used as NA
      Bool_t      fHasNA;        //TRUE if data have missing values

   public:
      XSetting();
      XSetting(const char *type, const char *infile);
      virtual ~XSetting();

      virtual Int_t InitAlgorithm(const char * /*name*/, const char * /*type*/,
                       Option_t * /*options*/, const char * /*filename*/,
                       Int_t /*npars*/, Double_t * /*pars*/) {return 0;}
      virtual void  ResetAlgorithm(const char * /*name*/, const char * /*type*/);

      void   InitNA(Double_t na)        {fNA = na; fHasNA = kTRUE;}

      ClassDef(XSetting,1) //Setting
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreeSet                                                             //
//                                                                      //
// Base class containing info about tree sets stored in TFile           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTreeSet: public TNamed {

   protected:
      TList        *fHeaders;      //headers of trees in this set
      TList        *fTrees;        //! list of selected trees in this set
      TList        *fSelections;   //! list of names of selected trees
      TList        *fTrash;        //! objects to be deleted from memory
      TFile        *fFile;         //! currently open file to store tree(s)
      XSetting     *fSetting;      //! parameter settings
      XManager     *fManager;      //! manager owing treeset
      TTree        *fTree;         //! selected tree
      TString       fInfile;       //! Name of infile if data were imported
      TString       fTreeName;     //! Name of tree to be added to this set
      Bool_t        fAsXML;        //! import/export data of type XML
      static Bool_t fgPrintHeader; //! print header for content only once

   protected:
      virtual Int_t IsBinaryFile(std::ifstream &/*input*/)            {return 0;}
      virtual Int_t ReadHeader(std::ifstream &/*input*/, const char * /*sep*/, char /*delim*/)
                                                                 {return 0;}
      virtual Int_t ReadData(std::ifstream &/*input*/, Option_t * /*option*/, const char * /*sep*/,
                       char /*delim*/, Int_t /*split*/)          {return 0;}
      virtual Int_t ReadBinaryHeader(std::ifstream &/*input*/, const char * /*sep*/, char /*delim*/)
                                                                 {return 0;}
      virtual Int_t ReadBinaryData(std::ifstream &/*input*/, Option_t * /*option*/, const char * /*sep*/,
                       char /*delim*/, Int_t /*split*/)          {return 0;}
      virtual Int_t ReadXMLHeader(std::ifstream &/*input*/, const char * /*sep*/, char /*delim*/)
                                                                 {return 0;}
      virtual Int_t ReadXMLData(std::ifstream &/*input*/, Option_t * /*option*/,
                       const char * /*sep*/, char /*delim*/, Int_t /*split*/)
                                                                 {return 0;}
      
   public :
      XTreeSet();
      XTreeSet(const char *name, const char *title = "");
      XTreeSet(const XTreeSet &treeset);
      XTreeSet& operator=(const XTreeSet& rhs);
      virtual ~XTreeSet();
		
      virtual void  AddTree(TTree *tree);
      virtual void  AddTreeInfo(TTree *tree, XTreeInfo *info);
      virtual void  AddTreeInfo(TTree *tree, const char *name, Option_t *option = "");
      virtual void  AddTreeHeader(const char *treename, Int_t treeid);
      virtual void  AddTreeHeader(const char *treename, const char *treetype,
                       Int_t treeid, Int_t npar, Double_t *pars);
      virtual void  RemoveTreeHeader(XTreeHeader *header);
      virtual void  RemoveTreeHeaders();
      virtual Int_t Export(const char *exten, const char *varlist,
                       std::ofstream &output, const char *sep);
      virtual Int_t ExportTrees(const char *exten, const char *varlist, 
                       std::ofstream &output, const char *sep);
      virtual Int_t ExportTree(const char *exten, Int_t n, TString *names,  
                       const char *varlist, std::ofstream &output, const char *sep);
      virtual Int_t ExportTreeInfo(const char * /*exten*/, Int_t /*n*/, TString * /*names*/,  
                       const char * /*varlist*/, std::ofstream &/*output*/, const char * /*sep*/)
                                                                            {return 0;}
      virtual Int_t ExportTreeType(const char * /*exten*/, Int_t /*n*/, TString * /*names*/,  
                       const char * /*varlist*/, std::ofstream &/*output*/, const char * /*sep*/)
                                                                            {return 0;}
      virtual Int_t ExportTreeXML(const char * /*exten*/, Int_t /*n*/, TString * /*names*/,  
                       const char * /*varlist*/, std::ofstream &/*output*/, const char * /*sep*/)
                                                                            {return 0;}
      virtual Int_t HandleOption(const char * /*name*/, Option_t * /*opt*/) {return 0;}
      virtual Int_t HandleOption(TTree * /*tree*/, Option_t * /*opt*/)      {return 0;}
      virtual Int_t Import(std::ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportXML(std::ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");
      void          Select(const char *name, Int_t id = 0);
      virtual void  PrintInfo() {}

      Double_t    **CreateTable(Int_t nrow, Int_t ncol);
      void          DeleteTable(Double_t **table, Int_t nrow);

      Int_t         WriteTree(TTree *tree, Int_t option = 0, Int_t bufsize = 0);
      Int_t         DeleteTree(const char *name, const char *exten,
                       const char *cycle);
      TString       FindTree(const char *name);

      XTreeInfo    *GetTreeInfo(const char *name, TTree *tree);
      XTreeHeader  *GetTreeHeader(const char *name);

      XManager     *GetManager()       const {return fManager;}
      TList        *GetTreeHeaders()   const {return fHeaders;}
      TList        *GetListOfTrees()   const {return fTrees;}
      TList        *GetSelections()    const {return fSelections;}
      Int_t         GetNumHeaders()    const {return fHeaders->GetSize();}
      Int_t         GetNumTrees()      const {return fTrees->GetSize();}
      Int_t         GetNumSelections() const {return fSelections->GetSize();}

      void   SetManager(XManager *manager)   {fManager = manager;}
      void   SetTree(TTree *tree)            {fTree = tree;}
      Bool_t AsXML()                   const {return fAsXML;}
      void   AsXML(Bool_t asXML)             {fAsXML = asXML;}

      static  void   SetPrintHeader(Bool_t print) {fgPrintHeader = print;}
      static  Bool_t GetPrintHeader()             {return fgPrintHeader;}

      ClassDef(XTreeSet,1) //TreeSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAlgorithm                                                           //
//                                                                      //
// Base class for algorithms                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAlgorithm: public TNamed {

   protected:
      TFile     *fFile;        //! optional private file
      XTreeInfo *fTreeInfo;    //! user info for current tree
      TString    fOption;      //option
      Int_t      fNPar;        //number of parameters
      Double_t  *fPars;        //[fNPar] Array of fNpar parameters
      Double_t   fNA;          //value to be used as NA
      Bool_t     fHasNA;       //TRUE if data have missing values
      Bool_t     fIsFileOwner; //! delete fFile only if TRUE

   protected:
      Double_t *Array2Log(Int_t n, Double_t *x, Double_t neglog, const char *base);
      Double_t *Array2Pow(Int_t n, Double_t *x, const char *base);
      Int_t     TestNumParameters(Int_t npar);

      Int_t     XY2Index(Int_t x, Int_t y, Int_t ncol) {return (x + y*ncol);  }
      Int_t     Index2X(Int_t index, Int_t ncol)       {return (index % ncol);}
      Int_t     Index2Y(Int_t index, Int_t ncol)       {return (index / ncol);}

   public:
      XAlgorithm();
      XAlgorithm(const char *name, const char *type);
      XAlgorithm(const XAlgorithm &algorithm);
      XAlgorithm& operator=(const XAlgorithm& rhs);
      virtual ~XAlgorithm();

      virtual Int_t InitType(const char * /*type*/, Option_t * /*options*/,
                       Int_t /*npars*/, Double_t * /*pars*/) {return 0;}

      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);
      virtual Int_t Calculate(Int_t n, Int_t *x, Int_t *msk);
      virtual Int_t Calculate(Int_t n, Double_t *x, Int_t *msk);
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z, Int_t *msk);
      virtual Int_t Calculate(Int_t n, Double_t *x, Int_t *idx, Int_t *msk);
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *idx, Int_t *msk);
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z, Int_t *idx, Int_t *msk);
      virtual Int_t Calculate(Int_t nrow, Int_t ncol, Double_t **table);

      virtual void      SetOptions(Option_t *opt) {fOption = opt; fOption.ToLower();}
      virtual Option_t *GetOptions(const char *sep);
      virtual Option_t *GetOptions(TString &suboption, const char *sep);

      Double_t **CreateTable(Int_t nrow, Int_t ncol);
      void       DeleteTable(Double_t **table, Int_t nrow);

      TFile    *NewFile(const char *name, const char *exten);
      Int_t     InitParameters(Int_t npar, Double_t *pars);
      void      DeleteParameters();

      void      InitNA(Double_t na)           {fNA = na; fHasNA = kTRUE;}
      void      SetFile(TFile *file, Bool_t isOwner = kFALSE)
                                  {fFile = file; fIsFileOwner = isOwner;}
      void      SetTreeInfo(XTreeInfo *info)  {fTreeInfo = info;}
      void      SetOption(Option_t *opt)      {fOption  = opt; fOption.ToLower();}
      Option_t *GetOption()             const {return fOption.Data();}
      TFile    *GetFile()               const {return fFile;}
      Int_t     GetNumParameters()      const {return fNPar;}
      Double_t *GetParameters()         const {return fPars;}

      ClassDef(XAlgorithm,1) //Algorithm
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XManager                                                             //
//                                                                      //
// Base class for managing files                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XManager: public TNamed {

   protected:
      XFolder  *fContent;        //content of file
      TString   fDataType;       //data type
      TString   fOption;         //option
      TFile    *fFile;           //! ROOT file 
      XSetting *fSetting;        //! parameter settings
      XTreeSet *fTreeSet;        //! current tree set
      XPlot    *fPlotter;        //! plotter for drawing
      TList    *fTrash;          //! objects to be deleted from memory
      Bool_t    fIsFileOwner;    //! delete fFile only if TRUE
      Bool_t    fAbort;          //! abort further actions
      Bool_t    fInterrupt;      //! interrupt further input
      Bool_t    fInitFlag;       //! flag set to TRUE in New/Open/Update

   private:
      static Int_t fgBufSize;    //! buffer size of tree baskets

   public:  //better private and use GetVerbose()
      static Int_t fgVerbose;    //! print information if verbose=kTRUE

   protected:
      virtual Int_t     ImportDefaults(const char * /*infile*/)  {return 0;}
      virtual Int_t     InitDefaults()                           {return 0;}
      virtual Bool_t    IsOpen(TFile *file, const char *filename);
      virtual XSetting *NewSetting(const char *type, const char *infile);
      virtual XTreeSet *NewTreeSet(const char * /*type*/)        {return 0;}
      virtual XPlot    *NewPlotter(const char *name, const char *title = "");
      virtual Int_t     ExportTrees(const char *exten, const char *varlist, 
                           std::ofstream &output, Bool_t asXML, const char *sep);
      virtual Int_t     DeleteTreeSetInfo(const char * /*name*/) {return 0;}
      virtual Int_t     HandleOption(const char *setname, const char *treename,
                           Option_t *option);

      TFile   *NewFile(const char *name, const char *title);
      TFile   *OpenFile(const char *name, Option_t *option, Bool_t &isOwner);
      void     CloseFile(TFile *file, Option_t *option = "");
      void     DeleteFile(TFile *file);
      void     DeleteDirectory(const char *name, const char *cycle);
      XFolder *NewContent(const char *name, const char *title, const char *type);

   public:
      XManager();
      XManager(const char *name, const char *title = "", Int_t verbose = kTRUE);
      virtual ~XManager();

      virtual Int_t  Initialize(const char *type = "", const char *data = "",
                        Option_t *option = "", const char *infile = "",
                        Bool_t hasPlotter = kFALSE);
      virtual Int_t  InitNA(Double_t na);
      virtual Int_t  InitAlgorithm(const char *name, const char *type,
                        Option_t *options, const char *filename = 0,
                        Int_t npars = 0, Double_t p1 = 0, Double_t p2 = 0, Double_t p3 = 0,
                        Double_t p4 = 0, Double_t p5 = 0, Double_t p6 = 0, Double_t p7 = 0,
                        Double_t p8 = 0, Double_t p9 = 0, Double_t p10 = 0);
      virtual void   ResetAlgorithm(const char *name, const char *type = "");

      virtual Int_t  New(const char *name, const char *dir = "",
                         const char *type = "", const char *data = "",
                         Option_t *option = "");
      virtual Int_t  Open(const char *fullname, const char *data = "",
                        Option_t *option = "",  Option_t *fileopt = "READ");
      virtual Int_t  Update(const char *fullname, const char *data = "",
                        Option_t *option = "",  const char *userID = "",
                        const char *password = "", Option_t *fileopt = "UPDATE");
      virtual void   Close(Option_t *option = "");
      virtual Bool_t Save();

      virtual Int_t  AddTree(const char *setname, const char *intree,
                        Int_t treeid = 1, Option_t *option = "");
      virtual Int_t  ExportSet(const char *setname, const char *exten, 
                        const char *varlist = "*", const char *outfile = "",
                        const char *sep = "\t");
      virtual Int_t  Export(const char *treename, const char *varlist = "*", 
                        const char *outfile = "", const char *sep = "\t");
      virtual Int_t  Import(const char *setname, const char *infile,
                        const char *treename = "", Option_t *option = "NEW",
                        const char *sep = "\t",  char delim = '\n', Int_t split = 99);
      virtual void   Delete(const char *name);
      virtual Int_t  DeleteTree(const char *namecycle);
      virtual Int_t  DeleteTreeSet(const char *setname);
      virtual Int_t  HandleError(Int_t err, const char *name1 = "",
                        const char *name2 = "");

      virtual void   NewCanvas(const char *name, const char *title,
                        Int_t wtopx = 20, Int_t wtopy = 20, Int_t ww = 500,
                        Int_t wh = 500, Int_t nx = 1, Int_t ny = 1);
      virtual void   CloseCanvas(Option_t *opt);

      using TObject::Draw;
      virtual Int_t  Draw(const char *canvasname, const char *treename,
                        const char *varlist, const char *logbases,
                        const char *type, Option_t *opt,
                        Double_t minX = -1111, Double_t maxX = -1111,
                        Double_t minY = -1111, Double_t maxY = -1111,
                        Double_t minZ = -1111, Double_t maxZ = -1111,
                        const char *var2sort = "", Bool_t down = kFALSE);
      virtual Int_t  Draw(const char *canvasname, const char *treename1,
                        const char *treename2,  const char *varlist,
                        const char *logbases, const char *type, Option_t *opt,
                        Double_t minX = -1111, Double_t maxX = -1111,
                        Double_t minY = -1111, Double_t maxY = -1111,
                        Int_t sort = 0, Bool_t down = kFALSE);
      virtual Int_t  Draw(const char *canvasname, const char *treename1,
                        const char *treename2, const char *treename3,
                        const char *varlist, const char *logbases,
                        const char *type, Option_t *opt, Double_t minX = -1111,
                        Double_t maxX = -1111, Double_t minY = -1111,
                        Double_t maxY = -1111, Double_t minZ = -1111,
                        Double_t maxZ = -1111, Int_t sort = 0, Bool_t down = kFALSE);
      virtual Int_t  DrawImage(const char *canvasname, const char *treename,
                        const char *varlist, const char *logbase, Option_t *opt,
                        Double_t min = 0, Double_t max = 0, Option_t *orientation = "U");
      virtual Int_t  DrawTree(const char *canvasname, const char *treename,
                        const char *varexp, const char *selection, Option_t *opt);
      virtual Int_t  DrawEntries(const char *canvasname, const char *leafname,
                        Int_t n, Int_t *entrylist, const char *logbase,
                        const char *type, Option_t *opt, Double_t min = -1111,
                        Double_t max = -1111, Int_t sort = 0, Bool_t down = kFALSE);
      virtual Int_t  DrawLeaves(const char *canvasname, const char *leafname,
                        const char *logbase, const char *type, Option_t *opt,
                        Double_t min = -1111, Double_t max = -1111,
                        Int_t sort = 0, Bool_t down = kFALSE);

      virtual TTree       *GetTree(const char *treename);
      virtual XTreeHeader *GetTreeHeader(const char *treename);
      virtual XTreeInfo   *GetTreeInfo(const char *treename);

      void PrintContents(const char *setname);
      void PrintContents();
      void SetMaxFileSize(Long64_t maxsize = 1900000000);

      void ResetAbort()                  {fAbort = kFALSE;}  //use with caution!
      void SetFileOwner(Bool_t isOwner)  {fIsFileOwner = isOwner;}
      void SetDataType(const char *type) {fDataType = type;}

      Bool_t   GetFileOwner() const {return fIsFileOwner;}
      TString  GetDataType()  const {return fDataType;}
      TString  GetFileName()  const {return fFile->GetName();}
      TFile   *GetFile()      const {return fFile;}
      XFolder *GetContent()   const {return fContent;}
      XPlot   *GetPlotter()   const {return fPlotter;}

      static void  SetBufSize(Int_t bufsize = 32000);
      static void  SetVerbose(Int_t verbose = 1);
      static Int_t GetBufSize(Int_t n = 1, Int_t m = 10000);
      static Int_t GetVerbose();


      ClassDef(XManager,1) //Manager
};

#endif
