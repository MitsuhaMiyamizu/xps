// File created: 11/02/2002                          last modified: 04/26/2009
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2009 Dr. Christian Stratowa
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

#ifndef __XPSUtils__
#define __XPSUtils__

#include <Riostream.h>
#include <cmath>

#include "TCanvas.h"
#include "TNamed.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TStyle.h"
#include "TArrayI.h"

// define path separators for e.g. Path2Name() and Name2Path()
#define dSEP "/"
#define sSEP '/'
// Note: both ROOT and R allow "/" also on Windows XP and Vista
//       thus "\\" is not necessary but could be activated here:
//#ifdef WIN32
//#   define dSEP "\\"
//#   define sSEP '\\'
//#else
//#   define dSEP "/"
//#   define sSEP '/'
//#endif


// Error messages: plotNoErr must be zero!!
enum EPlotErrors {
   perrNoErr         =  0,
   perrFatal         = -1,
   perrAbort         = -2,
   perrGeneral       = -3,
   perrInitMemory    = -4,
   perrCreateFile    = -5,
   perrCreateDir     = -6,
   perrCreateTree    = -7,
   perrCreateTreeSet = -8,
   perrGetFile       = -9,
   perrGetDir        = -10,
   perrGetTree       = -11,
   perrGetTreeSet    = -12,
   perrOpenFile      = -100,
   perrTreeName      = -101,
   perrNumPads       = -102,
   perrGetLeaf       = -103,
   perrGetBrnch      = -104,
   perrPlotType      = -105,
   perrNumEntries    = -106,
   perrNoImplement   = -107,
   perrNoMsg         = -999,
};

// an affy generic STRING
typedef struct
{
   Int_t len;
   char *value;
} ASTRING;

// an affy generic WSTRING
typedef struct
{
   Int_t    len;
   wchar_t *value;
} AWSTRING;

// String functions
extern Int_t   CheckHeader(const char *header, const char **kHeader,
                  const Int_t ncols, Int_t *index, const char *sep);
extern Int_t   CheckHeaderOrder(const char *header, const char **kHeader,
                  const Int_t ncols, Int_t *index, const char *sep);
extern Int_t   GetHeaderOrder(const char *header, const char **kHeader, 
                  const Int_t ncols, Int_t *index, const char *sep);
extern TString Extension2Type(const char *exten, const char **types,
                  const char **extens);
extern char   *FirstPath(const char *name);
extern TString FirstPath(const char *fullname, char sep);
extern TString FullName(const char *dir, const char *name, const char *sep);
extern TString GetROOTName(const char *fullname);
extern Bool_t  HasExtension(const char *exten, Int_t n, const char **extens);
extern Bool_t  HasExtension(const char *exten, const char **extens);
extern TString Name2Path(const char *fullname, char sep);
extern Bool_t  NameInArray(const char *name, TString *names, Int_t n);
extern Int_t   NumSeparators(const char *name, const char *sep);
extern Int_t   NumSeparators(const char *name, const char sep);
extern TString Path2Name(const char *name, const char *sep, const char *exten);
extern TString RemoveEnds(const char *name);
extern TString RemoveEnds(const char *name, Int_t &begin, Int_t &end);
extern TString RemoveLeadingSpace(const char *name, Int_t &begin);
extern TString RemoveSubString(const char *str, const char *substr, Bool_t exact);
extern TString RemoveNonAlpha(const char *name);
extern TString ReplaceNonAlpha(const char *name, const char *sep);
extern Int_t   StringInList(const char *str, const char **kStrList,
                  const Int_t n, Bool_t exact);
extern TString SubString(const char *str, const char *sep, Int_t n);
extern TString SubString(const char *str, char sep1, char sep2, Bool_t source);
extern Int_t   TokenizeString(const char *cstr, Int_t &n, TString *names,
                  const char *sep);
extern Int_t   TokenizeString(const char *cstr, Int_t &n, TString *names,
                  Int_t len, const char *sep);
extern TString Type2Extension(const char *type, const char **types,
                  const char **extens);

// Functions for file parsing
extern void SwapBytes(char *src, char *dest, Int_t size);
extern void READ_INT(std::ifstream &input, Int_t &value, Bool_t isBE = kFALSE);
extern void READ_UINT(std::ifstream &input, UInt_t &value, Bool_t isBE = kFALSE);
extern void READ_SHORT(std::ifstream &input, Short_t &value, Bool_t isBE = kFALSE);
extern void READ_USHORT(std::ifstream &input, UShort_t &value, Bool_t isBE = kFALSE);
extern void READ_CHAR(std::ifstream &input, char &value);
extern void READ_UCHAR(std::ifstream &input, unsigned char &value);
extern void READ_BOOL(std::ifstream &input, char &value);
extern void READ_FLOAT(std::ifstream &input, Float_t &value, Bool_t isBE = kFALSE);
extern void READ_STRING(std::ifstream &input, char * &value, Bool_t isBE = kFALSE);
extern void READ_STRING(std::ifstream &input, ASTRING *value, Bool_t isBE = kFALSE);
extern void READ_WSTRING(std::ifstream &input, char * &value, Bool_t isBE = kFALSE);
extern void READ_WSTRING(std::ifstream &input, wchar_t * &value, Bool_t isBE = kFALSE);
extern void READ_WSTRING(std::ifstream &input, AWSTRING *value, Bool_t isBE = kFALSE);
extern void READ_FIXED_STRING(std::ifstream &input, char *value, Int_t len);

// Functions to decode Affymetrix MIME types
extern wchar_t *DecodeTEXT(ASTRING *value);
extern Int_t    DecodeINT(ASTRING *value);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBitSet                                                              //
//                                                                      //
// BitSet Class extracted from TObject                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XBitSet {

   private:
      UInt_t fBits;

   public:
      enum {
         kBitClear = 0x00000000,
         kBitMask  = 0x00ffffff
      };

      XBitSet();
      virtual ~XBitSet();

      void    SetBit(UInt_t f)         { fBits |= f & kBitMask; }
      void    ResetBit(UInt_t f)       { fBits &= ~(f & kBitMask); }
      Bool_t  TestBit(UInt_t f)  const { return (Bool_t) ((fBits & f) != 0); }
      Int_t   TestBits(UInt_t f) const { return (Int_t) (fBits & f); }
      void    InvertBit(UInt_t f)      { fBits ^= f & kBitMask; }

      ClassDef(XBitSet,1) //BitSet
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPlot                                                                //
//                                                                      //
// Class for drawing graphs, histograms and images                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPlot: public TNamed {

   protected:
      TFile    *fFile;         //! current ROOT file 
      TList    *fTrees;        //list of trees
//??      TTree    *fTree;         //! selected tree
      TCanvas  *fCanvas;       //! current canvas
      Int_t     fNPads;        //! number of pads in current canvas
      Int_t     fPadNr;        //! number of current pad
      Double_t  fMin;          //minimum value for all axes
      Double_t  fMax;          //maximum value for all axes
      Double_t  fMinX;         //minimum value for X-axis
      Double_t  fMaxX;         //maximum value for X-axis
      Double_t  fMinY;         //minimum value for Y-axis
      Double_t  fMaxY;         //maximum value for Y-axis
      Double_t  fMinZ;         //minimum value for Z-axis
      Double_t  fMaxZ;         //maximum value for Z-axis
      Double_t  fNegLog;       //replacement for negative data in log
      Int_t     fNNegX;        //number of values <=0 on X-axis
      Int_t     fNNegY;        //number of values <=0 on Y-axis
      Int_t     fNNegZ;        //number of values <=0 on Z-axis
      Int_t     fNBinsX;       //number of X-bins for histogram (default = 50)
      Int_t     fNBinsY;       //number of Y-bins for histogram
      Int_t     fNBinsZ;       //number of Z-bins for histogram
      TArrayI   fLineColors;   //array of line colors
      TArrayI   fLineStyles;   //array of line styles
      TArrayI   fMarkerColors; //array of marker colors
      TArrayI   fMarkerStyles; //array of marker styles (see TAttMarker.cxx)
      Int_t     fPriorityLC;   //priority of line colors
      Int_t     fPriorityLS;   //priority of line styles
      Int_t     fPriorityMC;   //priority of marker colors
      Int_t     fPriorityMS;   //priority of marker styles
      TString   fTitle;        //title of frame
      TString   fTitleX;       //title of X-axis
      TString   fTitleY;       //title of Y-axis
      TString   fTitleZ;       //title of Z-axis
      Int_t     fSetTitle;     //set title of frame: auto <0, no title = 0
      Int_t     fSetTitleX;    //set title of X-axis: auto <0, no title = 0
      Int_t     fSetTitleY;    //set title of Y-axis: auto <0, no title = 0
      Int_t     fSetTitleZ;    //set title of Z-axis: auto <0, no title = 0
      Bool_t    fEqualAxes;    //! if TRUE all axes have equal range
      Bool_t    fIsFileOwner;  //! TRUE if owner of fFile
      Bool_t    fAbort;        //! abort further actions

   protected:
     TFile *OpenFile(const char *name, Option_t *option, Bool_t &isOwner);
     Bool_t IsOpen(TFile *file, const char *filename);
     TTree *GetTree(const char *fullname);
     void   FillArrays(Int_t n, TBranch *brch1, TLeaf *leaf1, Double_t *index,
               Double_t *arrX, Int_t *base);
     void   FillArrays(Int_t n, TBranch *brch1, TLeaf *leaf1, TBranch *brch2,
               TLeaf *leaf2, Double_t *arrX, Double_t *arrY, Int_t *base);
     void   FillArrays(Int_t n, TBranch *brch1, TLeaf *leaf1, TBranch *brch2,
               TLeaf *leaf2, TBranch *brch3, TLeaf *leaf3, Double_t *arrX,
               Double_t *arrY, Double_t *arrZ, Int_t *base);
     Int_t  FillEntry(Int_t entry, const char *leafname, Int_t n, Double_t *arr,
               Int_t base);
     void   FillEntrylist(Int_t n, TBranch *brch, TLeaf *leaf, Int_t *entrylist,
               Double_t *arr, Int_t base);

   public:
      XPlot();
      XPlot(const char *name, const char *title = "");
      virtual ~XPlot();

      virtual void  NewCanvas(const char *name, const char *title,
                       Int_t wtopx = 20, Int_t wtopy = 20, Int_t ww = 500,
                       Int_t wh = 500, Int_t nx = 1, Int_t ny = 1);
      virtual void  CloseCanvas(Option_t *opt);

      using TObject::SaveAs;
      virtual Int_t SaveAs(const char *name);

      using TObject::Draw;
      virtual Int_t Draw(const char *canvasname, const char *treename,
                       const char *varlist, const char *logbases, const char *type,
                       Option_t *opt, Double_t minX = -1111, Double_t maxX = -1111,
                       Double_t minY = -1111, Double_t maxY = -1111,
                       Double_t minZ = -1111, Double_t maxZ = -1111,
                       const char *var2sort = "", Bool_t down = kFALSE);
      virtual Int_t Draw(const char *canvasname, const char *treename1,
                       const char *treename2,  const char *varlist,
                       const char *logbases, const char *type, Option_t *opt,
                       Double_t minX = -1111, Double_t maxX = -1111,
                       Double_t minY = -1111, Double_t maxY = -1111,
                       Int_t sort = 0, Bool_t down = kFALSE);
      virtual Int_t Draw(const char *canvasname, const char *treename1,
                       const char *treename2, const char *treename3,
                       const char *varlist, const char *logbases, const char *type,
                       Option_t *opt, Double_t minX = -1111, Double_t maxX = -1111,
                       Double_t minY = -1111, Double_t maxY = -1111,
                       Double_t minZ = -1111, Double_t maxZ = -1111,
                       Int_t sort = 0, Bool_t down = kFALSE);

      virtual Int_t DrawImage(const char *canvasname, const char *treename,
                       const char *varlist, const char *logbase, Option_t *opt,
                       Double_t min = 0, Double_t max = 0, Option_t *orientation = "U");
      virtual Int_t DrawTree(const char *canvasname, const char *treename,
                       const char *varexp, const char *selection, Option_t *opt);

      virtual Int_t DrawEntries(const char *canvasname, const char *leafname,
                       Int_t n, Int_t *entrylist, const char *logbase,
                       const char *type, Option_t *opt, Double_t min = -1111,
                       Double_t max = -1111, Int_t sort = 0, Bool_t down = kFALSE);
      virtual Int_t DrawLeaves(const char *canvasname, const char *leafname,
                       const char *logbase, const char *type, Option_t *opt,
                       Double_t min = -1111, Double_t max = -1111,
                       Int_t sort = 0, Bool_t down = kFALSE);
      virtual Int_t DrawDensity(const char *canvasname, const char *leafname,
                       const char *logbase, const char *kernel, Option_t *opt,
                       Int_t npts, Double_t min, Double_t max);
      virtual Int_t DrawParallelCoord(const char *canvasname, const char *varlist,
                       Double_t min = 0, Double_t max = 0,
                       Bool_t aslog = kTRUE, Bool_t gl = kTRUE, Bool_t can = kFALSE);

      void  DrawGraph1D(Int_t n, Double_t *index, Double_t *x, Option_t *opt,
               Int_t sort = 0);
      void  DrawGraph2D(Int_t n, Double_t *x, Double_t *y, Option_t *opt);
      void  DrawGraph3D(Int_t n, Double_t *x, Double_t *y, Double_t *z,
               Option_t *opt);
      void  DrawMultiGraph(Int_t n, Double_t *x, Double_t *y, Option_t *opt,
               Int_t sort = 0, Bool_t down = kFALSE);
      void  DrawMultiGraph(Int_t n, Double_t *x, Double_t *y, Double_t *z,
               Option_t *opt, Int_t sort = 0, Bool_t down = kFALSE);
      void  DrawMultiGraph(Int_t n, Int_t m, Double_t **table, Option_t *opt,
               Int_t sort = 0, Bool_t down = kFALSE);
      void  DrawDensity(Int_t n, Double_t *index, Double_t *x, Int_t npts,
               Option_t *opt, const char *kernel = "epanechnikov");
      void  DrawDensity(Int_t n, Int_t m, Double_t **table, Int_t npts,
               Option_t *opt, const char *kernel = "epanechnikov");
      void  DrawHist1D(Int_t n, Double_t *index, Double_t *x, Option_t *opt);
      void  DrawHist2D(Int_t n, Double_t *x, Double_t *y, Option_t *opt);
      void  DrawHist3D(Int_t n, Double_t *x, Double_t *y, Double_t *z,
               Option_t *opt);
      void  DrawMVA(Int_t n, Double_t *x, Double_t *y, Int_t base, Option_t *opt);

      Int_t Open(const char *filename, Option_t *option = "READ");
      Int_t AddTree(const char *treename);
      void  ClearTrees();
      void  SetDefault();

      void  SetFillColor(Int_t n = 0, Int_t *colors = 0, Int_t priority = 0);
      void  SetLineColor(Int_t n = 0, Int_t *colors = 0, Int_t priority = 0);
      void  SetLineStyle(Int_t n = 0, Int_t *styles = 0, Int_t priority = 0);
      void  SetMarkerColor(Int_t n = 0, Int_t *colors = 0, Int_t priority = 0);
      void  SetMarkerStyle(Int_t n = 0, Int_t *styles = 0, Int_t priority = 0);
      void  SetTitleMain(const char *title, Int_t setTitle = 1);
      void  SetTitleX(const char *title, Int_t setTitle = 1, Int_t base = 0);
      void  SetTitleY(const char *title, Int_t setTitle = 1, Int_t base = 0);
      void  SetTitleZ(const char *title, Int_t setTitle = 1, Int_t base = 0);

      TString LogTitle(const char *title, Int_t base = -1);

      void  SetNegLog(Double_t neglog) {fNegLog = neglog;}
      void  SetNumBins(Int_t binsX, Int_t binsY = 0, Int_t binsZ = 0);
      void  SetPalette(Int_t ncolors = 1, Int_t *colors = 0)
                                       {gStyle->SetPalette(ncolors,colors);}
      void  SetPadCount(Int_t n = 0)   {fPadNr = n;}
      Int_t IncreasePadCount()         {return ++fPadNr;}
      Int_t DecreasePadCount()         {return --fPadNr;}
      void  ClearPad(Int_t padnr = 0);

      void  SetCanvas(TCanvas *canvas, Int_t npads = 1, Int_t padnr = 0)
                          {fCanvas = canvas; fNPads = npads; fPadNr = padnr;}
      void  SetFile(TFile *file)       {fFile = file; fIsFileOwner = kFALSE;}

      ClassDef(XPlot,1) //Plot
};

//______________________________________________________________________________
inline void XPlot::SetNumBins(Int_t binsX, Int_t binsY, Int_t binsZ)
{
   fNBinsX = (binsX > 0) ? binsX : 50;
   fNBinsY = (binsY > 0) ? binsY : fNBinsX;
   fNBinsZ = (binsZ > 0) ? binsZ : fNBinsY;
}//SetNumBins

#endif
