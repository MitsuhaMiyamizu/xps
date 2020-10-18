// File created: 11/02/2002                          last modified: 12/26/2010
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
* Nov 2002 - Initial versions finished
* Dec 2002 - Move class XStat to own library.
* Aug 2003 - Add enum EPlotErrors.
* Apr 2005 - Add functions for file parsing from Affymetrix GPL code.
* May 2005 - Fuse libraries libXPSUtils and libXPSBase => library libXPSBase.
* Dec 2005 - Replace TTree::AddFriend() with TList fTrees.
*          - Get trees from XTreeSet in separate TDirectory.
* Jan 2008 - Add/update functions for file parsing.
*          - Add functions to decode Affymetrix MIME types*
* Apr 2009 - Add functions GetHeaderOrder() and TokenizeString()
******************************************************************************/


using namespace std;

// IMPORTANT: must be defined for PowerPC and undefined for Intel PCs
#include <RConfig.h>
#if !defined(R__BYTESWAP)
   #define IS_BIG_ENDIAN 1
#endif

#include <cstdlib>
#include <new>  //needed for new (nothrow)

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TFrame.h"
#include "TKey.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TObjString.h"
#include "TParallelCoord.h"
#include "TParallelCoordVar.h"
#include "TParallelCoordRange.h"
#include "TROOT.h"
#include "TSystem.h"

#include "TStat.h"

#include "XPSUtils.h"

#define INT_SIZE sizeof(Int_t)
#define UINT_SIZE sizeof(UInt_t)
#define SHORT_SIZE sizeof(Short_t)
#define USHORT_SIZE sizeof(UShort_t)
#define CHAR_SIZE sizeof(char)
#define UCHAR_SIZE sizeof(unsigned char)
#define FLOAT_SIZE sizeof(Float_t)

//debug: print class method names
const Bool_t kCS  = 0;
const Bool_t kCSa = 0;

ClassImp(XBitSet);
ClassImp(XPlot);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBitSet                                                              //
//                                                                      //
// BitSet Class extracted from TObject                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBitSet::XBitSet(): fBits(kBitClear)
{
}//Constructor

//______________________________________________________________________________
XBitSet::~XBitSet()
{
}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPlot                                                                //
//                                                                      //
// Class for drawing graphs, histograms and images                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPlot::XPlot()
      :TNamed()
{
   // Default Plot constructor
   if(kCS) cout << "---XPlot::XPlot(default)------" << endl;

   fFile        = 0;
   fTrees       = 0;
   fCanvas      = 0;
   fNPads       = 1;
   fPadNr       = 0;
   fMin         = fMinX = fMinY = fMinZ = DBL_MAX;  //defined in float.h
   fMax         = fMaxX = fMaxY = fMaxZ = -DBL_MAX;
   fNBinsX      = fNBinsY = fNBinsZ = 50;
   fNegLog      = 1;
   fNNegX       = fNNegY = fNNegZ = 0;
   fPriorityLC  = 9999999;
   fPriorityLS  = 9999999;
   fPriorityMC  = 9999999;
   fPriorityMS  = 9999999;
   fTitle       = "";
   fTitleX      = "";
   fTitleY      = "";
   fTitleZ      = "";
   fSetTitle    = fSetTitleX = fSetTitleY = fSetTitleZ = -1;
   fEqualAxes   = kTRUE;
   fIsFileOwner = kFALSE;
   fAbort       = kFALSE;
}//Constructor

//______________________________________________________________________________
XPlot::XPlot(const char *name, const char *title)
      :TNamed(name, title)
{
   // Normal Plot constructor
   if(kCS) cout << "---XPlot::XPlot------" << endl;

   fFile        = 0;
   fTrees       = 0;
   fCanvas      = 0;
   fNPads       = 1;
   fPadNr       = 0;
   fMin         = fMinX = fMinY = fMinZ = DBL_MAX;
   fMax         = fMaxX = fMaxY = fMaxZ = -DBL_MAX;
   fNBinsX      = fNBinsY = fNBinsZ = 50;
   fNegLog      = 1;
   fNNegX       = fNNegY = fNNegZ = 0;
   fPriorityLC  = 9999999;
   fPriorityLS  = 9999999;
   fPriorityMC  = 9999999;
   fPriorityMS  = 9999999;
   fTitle       = "";
   fTitleX      = "";
   fTitleY      = "";
   fTitleZ      = "";
   fSetTitle    = fSetTitleX = fSetTitleY = fSetTitleZ = 1;
   fEqualAxes   = kTRUE;
   fIsFileOwner = kFALSE;
   fAbort       = kFALSE;

   this->SetLineColor();
   this->SetLineStyle();
   this->SetMarkerColor();
   this->SetMarkerStyle();
}//Constructor

//______________________________________________________________________________
XPlot::~XPlot()
{
   // Plot destructor
   if(kCS) cout << "---XPlot::~XPlot------" << endl;

   if(fTrees) {fTrees->Delete(); delete fTrees; fTrees = 0;}
}//Destructor

//______________________________________________________________________________
void XPlot::NewCanvas(const char *name, const char *title, Int_t wtopx,
            Int_t wtopy, Int_t ww , Int_t wh, Int_t nx, Int_t ny)
{
   // Create new canvas
   if(kCS) cout << "------XPlot::NewCanvas------" << endl;

   fCanvas = new TCanvas(name, title, wtopx, wtopy, ww, wh);

   fPadNr = 0;
   if ((nx > 1) || (ny > 1)) {
      fCanvas->Divide(nx, ny);
      fNPads = nx * ny;
   }//if

   this->SetPalette();
}//NewCanvas

//______________________________________________________________________________
void XPlot::CloseCanvas(Option_t *opt)
{
   // Close canvas
   if(kCS) cout << "------XPlot::CloseCanvas------" << endl;

   fCanvas->Close(opt);
}//CloseCanvas

//______________________________________________________________________________
Int_t XPlot::SaveAs(const char *name)
{
   // Save content of canvas
   if(kCS) cout << "------XPlot::SaveAs------" << endl;

   fCanvas->cd();
   fCanvas->Print(name);

   return 0;
}//SaveAs

//______________________________________________________________________________
void XPlot::ClearPad(Int_t padnr)
{
   // Clear pad with pad number padnr
   if(kCS) cout << "------XPlot::ClearPad------" << endl;

   TPad *pad = 0;
   pad = (TPad*)(fCanvas->GetPad(padnr));

   if (pad) {
      pad->Clear();
      pad->Modified();
      pad->Update();
   }//if
}//ClearPad

//______________________________________________________________________________
Int_t XPlot::Draw(const char *canvasname, const char *treename, const char *varlist,
             const char *logbases, const char *type, Option_t *opt,
             Double_t minX, Double_t maxX, Double_t minY, Double_t maxY,
             Double_t minZ, Double_t maxZ, const char *var2sort, Bool_t down)
{
   // Draw at most three variables from varlist for tree treename
   // Each variable is converted to logbase given in list logbases
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // The type can be: graph, multigraph, hist, density, profile
   // with the corresponding option opt (TGraph::PaintGraph, THistPainter::Paint)
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Variable var2sort determines the variable to sort for
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::Draw(1)------" << endl;

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

// Get tree
   TTree *tree = this->GetTree(treename);
   if (!tree) return perrGetTree;

// Get logbase for each dimension
   Int_t *base = new Int_t[3]; 
   base[0] = -1;
   base[1] = -1;
   base[2] = -1;

   Int_t dim  = 0;
//ccc char memory problem??
   char *log  = new char[strlen(logbases) + 1];
   char *dlog = log;
   log = strtok(strcpy(log,logbases),":");
//??   log = strtok((char*)logbases,":");
   while(log) {
      if (strcmp(log,"0")     == 0) {base[dim] = 0;}
      if (strcmp(log,"log")   == 0) {base[dim] = 1;}
      if (strcmp(log,"log2")  == 0) {base[dim] = 2;}
      if (strcmp(log,"log10") == 0) {base[dim] = 10;}
      log = strtok(NULL, ":");
      if ((log == 0) || (dim > 2)) break;
      dim++;
   }//while
   delete [] dlog;
   // if logbases contain fewer entries than dimen
   if (base[0] == -1) base[0] = 0;
   if (base[1] == -1) base[1] = base[0];
   if (base[2] == -1) base[2] = base[1];

// Get tree and leaf(s) for varlist
   TLeaf    *leaf1 = 0;
   TLeaf    *leaf2 = 0;
   TLeaf    *leaf3 = 0;
   TBranch  *brch1 = 0;
   TBranch  *brch2 = 0;
   TBranch  *brch3 = 0;
   Double_t *index = 0;
   Double_t *arrX  = 0;
   Double_t *arrY  = 0;
   Double_t *arrZ  = 0;
   TString   var1  = "";
   TString   var2  = "";
   TString   var3  = "";

   Int_t  entries  = (Int_t)(tree->GetEntries());
   Int_t  dimen    = 0;
   Int_t  arr2sort = 0;
   Int_t  sort     = 0;

   // get leaf1 for first variable
//ccc char memory problem??
   char  *varname  = new char[strlen(varlist) + 1];
   char  *dvarname = varname;
   varname = strtok(strcpy(varname,varlist),":");
//??   varname = strtok((char*)varlist,":");
   if (varname != 0) {
      if (!(leaf1 = tree->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch1 = leaf1->GetBranch())) {perr = perrGetBrnch;   goto cleanup;}
      if (!(index = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 
      if (!(arrX  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 
      if (strcmp(varname, var2sort) == 0)   arr2sort = 1;
      var1  = TString(varname);
      dimen = 1;
   } else {
      cerr << "Error: Variable name(s) for tree are missing." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if

   // get leaf2 for second variable
   varname = strtok(NULL, ":");
   if (varname != 0) {
      if (!(leaf2 = tree->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch2 = leaf2->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      if (!(arrY  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;}
      if (strcmp(varname, var2sort) == 0)   arr2sort = 2;
      var2  = TString(varname);
      dimen = 2;
   }//if

   // get leaf3 for third variable
   varname = strtok(NULL, ":");
   if (varname != 0) {
      if (!(leaf3 = tree->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch3 = leaf3->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      if (!(arrZ  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 
      if (strcmp(varname, var2sort) == 0)   arr2sort = 3;
      var3  = TString(varname);
      dimen = 3;
   }//if

// Fill arrays
   if (dimen == 1) {
      this->FillArrays(entries, brch1, leaf1, index, arrX, base);
   } else if (dimen == 2) {
      this->FillArrays(entries, brch1, leaf1, brch2, leaf2, arrX, arrY, base);
   } else if (dimen == 3) {
      this->FillArrays(entries, brch1, leaf1, brch2, leaf2, brch3, leaf3,
                       arrX, arrY, arrZ, base);
   }//if

   if (fNNegX > 0) {
      cout << "Warning: <" << fNNegX << "> data<=0 on X-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if
   if (fNNegY > 0) {
      cout << "Warning: <" << fNNegY << "> data<=0 on Y-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if
   if (fNNegZ > 0) {
      cout << "Warning: <" << fNNegZ << "> data<=0 on Z-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if

// Set axes range
   fMin = (fMinX < fMinY) ? fMinX : fMinY;
   fMin = (fMin  < fMinZ) ? fMin  : fMinZ;
   fMax = (fMaxX > fMaxY) ? fMaxX : fMaxY;
   fMax = (fMax  > fMaxZ) ? fMax  : fMaxZ;
   fEqualAxes = kFALSE;

   if (minX < maxX){
      fMinX = minX;
      fMaxX = maxX;
   } else {
;//      cout << "Warning: minX >= maxX, using computed values" << endl;
   }//if
   if (minY < maxY){
      fMinY = minY;
      fMaxY = maxY;
   } else {
;//      cout << "Warning: minY >= maxY, using computed values" << endl;
   }//if
   if (minZ < maxZ){
      fMinZ = minZ;
      fMaxZ = maxZ;
   } else {
;//      cout << "Warning: minZ >= maxZ, using computed values" << endl;
   }//if

   if ((minX == 0) && (maxX == 0)  &&
       (minY == 0) && (maxY == 0)  &&
       (minZ == 0) && (maxZ == 0)) {
      fMax = fMax + fMax/20;
      fEqualAxes = kTRUE;
   } else {
      if ((minX == -1111) && (maxX == -1111)) fMaxX = fMaxX + fMaxX/20;
      if ((minY == -1111) && (maxY == -1111)) fMaxY = fMaxY + fMaxY/20;
      if ((minZ == -1111) && (maxZ == -1111)) fMaxZ = fMaxZ + fMaxZ/20;
      fEqualAxes = kFALSE;
   }//if
//?? in GUI: display fMinX,fMinY,etc

// Draw type
   this->SetTitleMain(tree->GetName(), fSetTitle);
   if (strcmp(type, "graph") == 0) {
      switch (dimen) {
         case 1:
            this->SetTitleX("Index", fSetTitleX, base[0]);
            this->SetTitleY(var1, fSetTitleY, base[1]);
            if (arr2sort != 0) {sort = down ? -1 : 1;}
            this->DrawGraph1D(entries, index, arrX, opt, sort);
            break;

         case 2:
            this->SetTitleX(var1, fSetTitleX, base[0]);
            this->SetTitleY(var2, fSetTitleY, base[1]);
            this->DrawGraph2D(entries, arrX, arrY, opt);
            break;

         case 3:
            this->SetTitleX(var1, fSetTitleX, base[0]);
            this->SetTitleY(var2, fSetTitleY, base[1]);
            this->SetTitleZ(var3, fSetTitleZ, base[2]);
            this->DrawGraph3D(entries, arrX, arrY, arrZ, opt);
            break;

         default:
            printf("Dimension <%d> not allowed\n", dimen);
            perr = perrGeneral;
      }//switch
   } else if (strcmp(type, "multigraph") == 0) {
      switch (dimen) {
         case 1:
            this->SetTitleX("Index", fSetTitleX, base[0]);
            this->SetTitleY(var1, fSetTitleY, base[1]);
            if (arr2sort != 0) {sort = down ? -1 : 1;}
            this->DrawGraph1D(entries, index, arrX, opt, sort);
            break;

         case 2:
            this->SetTitleX("Index", fSetTitleX, base[0]);
            this->SetTitleY(var1, fSetTitleY, base[1]);
            this->SetTitleZ(var2, fSetTitleZ, base[1]);
            this->DrawMultiGraph(entries, arrX, arrY, opt, arr2sort, down);
            break;

         case 3:
            this->SetTitleX("Index", fSetTitleX, base[0]);
            this->SetTitleY(var1, fSetTitleY, base[1]);
            this->SetTitleZ(var2, fSetTitleZ, base[1]);
            this->DrawMultiGraph(entries, arrX, arrY, arrZ, opt, arr2sort, down);
            break;

         default:
            printf("Dimension <%d> not allowed\n", dimen);
            perr = perrGeneral;
      }//switch
   } else if (strcmp(type, "hist") == 0) {
      switch (dimen) {
         case 1:
            this->SetTitleX(var1, fSetTitleX, base[0]);
            this->SetTitleY("", fSetTitleY, base[1]);
            this->DrawHist1D(entries, index, arrX, opt);
            break;

         case 2:
            this->SetTitleX(var1, fSetTitleX, base[0]);
            this->SetTitleY(var2, fSetTitleY, base[1]);
            this->DrawHist2D(entries, arrX, arrY, opt);
            break;

         case 3:
            this->SetTitleX(var1, fSetTitleX, base[0]);
            this->SetTitleY(var2, fSetTitleY, base[1]);
            this->SetTitleZ(var3, fSetTitleZ, base[2]);
            this->DrawHist3D(entries, arrX, arrY, arrZ, opt);
            break;

         default:
            printf("Dimension <%d> not allowed\n", dimen);
            perr = perrGeneral;
      }//switch
   } else if (strcmp(type, "density") == 0) {
      this->SetTitleX(var1, fSetTitleX, base[0]);
      this->SetTitleY("", fSetTitleY, base[1]);
//TO DO: kernel as variable!!
      this->DrawDensity(entries, index, arrX, 512, opt, "epanechnikov");
   } else if (strcmp(type, "mvaplot") == 0) {
      cout << "Error: Type <" << type << "> not possible for 1D." << endl;
      perr = perrPlotType;
   } else if (strcmp(type, "profile") == 0) {
//TO DO
//      if (dimen == 1)      this->DrawGraph1D(entries, index, arrX, opt);
//      else if (dimen == 2) this->DrawGraph2D(entries, arrX, arrY, opt);
//      else if (dimen == 3) this->DrawGraph3D(entries, arrX, arrY, arrZ, opt);
   } else {
      cout << "Error: Drawing type <" << type << "> not known" << endl;
      perr = perrPlotType;
   }//if

// Cleanup
cleanup:
   if (base)  delete [] base;
   if (index) delete [] index;
   if (arrX)  delete [] arrX;
   if (arrY)  delete [] arrY;
   if (arrZ)  delete [] arrZ;

   delete [] dvarname;  //??

// prevent user from editing content of canvas
// not possible, cannot create more than one pad!!
//   fCanvas->SetEditable(kFALSE); //TEST!!!
   savedir->cd();
   return perr;
}//Draw

//______________________________________________________________________________
Int_t XPlot::Draw(const char *canvasname, const char *treename1,
             const char *treename2,  const char *varlist, const char *logbases,
             const char *type, Option_t *opt, Double_t minX, Double_t maxX,
             Double_t minY, Double_t maxY, Int_t sort, Bool_t down)
{
   // Draw variable from varlist for corresponding tree treename1 or treename2
   // Each variable is converted to logbase given in list logbases
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // type: graph, hist, mvaplot, profile - with corresponding option opt
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::Draw(2)------" << endl;

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

// Get trees
   TTree *tree1 = this->GetTree(treename1);
   if (!tree1) return perrGetTree;

   TTree *tree2 = this->GetTree(treename2);
   if (!tree2) return perrGetTree;

// Get logbase for each dimension
   Int_t *base = new Int_t[2]; 
   base[0] = -1;
   base[1] = -1;

   Int_t dim = 0;
//ccc char memory problem??
   char *log  = new char[strlen(logbases) + 1];
   char *dlog = log;
   log = strtok(strcpy(log,logbases),":");
//??   log = strtok((char*)logbases,":");
   while(log) {
      if (strcmp(log,"0")     == 0) {base[dim] = 0;}
      if (strcmp(log,"log")   == 0) {base[dim] = 1;}
      if (strcmp(log,"log2")  == 0) {base[dim] = 2;}
      if (strcmp(log,"log10") == 0) {base[dim] = 10;}
      log = strtok(NULL, ":");
      if ((log == 0) || (dim > 1)) break;
      dim++;
   }//while
   delete [] dlog;
   // if logbases contains fewer entries than dimen
   if (base[0] == -1) base[0] = 0;
   if (base[1] == -1) base[1] = base[0];

// Get leafs and variables
   TLeaf    *leaf1 = 0;
   TLeaf    *leaf2 = 0;
   TBranch  *brch1 = 0;
   TBranch  *brch2 = 0;
   Double_t *arrX  = 0;
   Double_t *arrY  = 0;
   TString   var1  = "";
   TString   var2  = "";
   TString   title = "";

   Int_t entries = (Int_t)(tree1->GetEntries());

   // get leaf1 for tree1
//ccc char memory problem??
   char *varname  = new char[strlen(varlist) + 1];
   char *dvarname = varname;
   varname = strtok(strcpy(varname,varlist),":");
//??   varname = strtok((char*)varlist,":");
   if (varname != 0) {
      if (!(leaf1 = tree1->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch1 = leaf1->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      if (!(arrX  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 
      var1 = TString(varname);
   } else {
      cerr << "Error: Variable name(s) for trees are missing." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if

   // get leaf2 for tree2
   varname = strtok(NULL, ":");
   if (varname != 0) {var2 = TString(varname);}
   else              {var2 = var1;}
   if (!(leaf2 = tree2->FindLeaf(var2))) {
      cerr << "Error: Leaf <" << var2 << "> not found." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if
   if (!(brch2 = leaf2->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
   if (!(arrY  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 

// Fill arrays
   this->FillArrays(entries, brch1, leaf1, brch2, leaf2, arrX, arrY, base);
   if (fNNegX > 0) {
      cout << "Warning: <" << fNNegX << "> data<=0 on X-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if
   if (fNNegY > 0) {
      cout << "Warning: <" << fNNegY << "> data<=0 on Y-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if

// Set axes range
   fMin = (fMinX < fMinY) ? fMinX : fMinY;
   fMax = (fMaxX > fMaxY) ? fMaxX : fMaxY;
   fEqualAxes = kFALSE;

   if (minX < maxX){
      fMinX = minX;
      fMaxX = maxX;
   } else {
;//      cout << "Warning: minX >= maxX, using computed values" << endl;
   }//if
   if (minY < maxY){
      fMinY = minY;
      fMaxY = maxY;
   } else {
;//      cout << "Warning: minY >= maxY, using computed values" << endl;
   }//if

   if ((minX == 0) && (maxX == 0) &&
       (minY == 0) && (maxY == 0)) {
      fMax = fMax + fMax/20;
      fEqualAxes = kTRUE;
   } else {
      if ((minX == -1111) && (maxX == -1111)) fMaxX = fMaxX + fMaxX/20;
      if ((minY == -1111) && (maxY == -1111)) fMaxY = fMaxY + fMaxY/20;
      fEqualAxes = kFALSE;
   }//if

//TEST
   // prevent x-axis range fMinX == fMaxX
   if ((fMinX >= fMaxX) && (fMaxX == 0)) {
      fMaxX = fMinX + 1.0/10;
      fEqualAxes = kFALSE;
   }//if

//?? in GUI: display fMinX,fMinY,etc

// Draw type
   title = TString(tree1->GetName()) + " vs "
         + TString(tree2->GetName());
   this->SetTitleMain(title.Data(), fSetTitle);
   if (strcmp(type, "graph") == 0) {
      this->SetTitleX(var1, fSetTitleX, base[0]);
      this->SetTitleY(var2, fSetTitleY, base[1]);
      this->DrawGraph2D(entries, arrX, arrY, opt);
   } else if (strcmp(type, "multigraph") == 0) {
      this->SetTitleX("Index", fSetTitleX, base[0]);
      this->SetTitleY(var1, fSetTitleY, base[1]);
      this->SetTitleZ(var2, fSetTitleZ, base[1]);
      this->DrawMultiGraph(entries, arrX, arrY, opt, sort, down);
   } else if (strcmp(type, "hist") == 0) {
      this->SetTitleX(var1, fSetTitleX, base[0]);
      this->SetTitleY(var2, fSetTitleY, base[1]);
      this->DrawHist2D(entries, arrX, arrY, opt);
   } else if (strcmp(type, "mvaplot") == 0) {
      if (base[1] == base[0]) {
         // symmetric M-axis if axes set automatically
         if (minY >= maxY){
            fMaxY = (fMaxY > TMath::Abs(fMinY)) ?  fMaxY : TMath::Abs(fMinY);
            fMinY = (- fMaxY);
         }//if
         this->SetTitleX("A", fSetTitleX, base[0]);
         this->SetTitleY("M", fSetTitleY, base[1]);
         this->DrawMVA(entries, arrX, arrY, base[0], opt);
      } else {
         cout << "Error: Logbase for tree1 and tree2 not identical." << endl;
         perr = perrGeneral;
      }//if
   } else if (strcmp(type, "profile") == 0) {
      this->SetTitleX(var1, fSetTitleX, base[0]);
      this->SetTitleY(var2, fSetTitleY, base[1]);
//to do      this->DrawGraph2D(entries, arrX, arrY, opt);
   } else {
      cout << "Error: Drawing type <" << type << "> not known" << endl;
      perr = perrPlotType;
   }//if

// Cleanup
cleanup:
   if (base) delete [] base;
   if (arrX) delete [] arrX;
   if (arrY) delete [] arrY;

   delete [] dvarname;  //??

   savedir->cd();
   return perr;
}//Draw

//______________________________________________________________________________
Int_t XPlot::Draw(const char *canvasname, const char *treename1,
             const char *treename2, const char *treename3, const char *varlist,
             const char *logbases, const char *type, Option_t *opt,
             Double_t minX, Double_t maxX, Double_t minY, Double_t maxY,
             Double_t minZ, Double_t maxZ, Int_t sort, Bool_t down)
{
   // Draw variable from varlist for corresponding tree treename1, 2 or 3
   // Each variable is converted to logbase given in list logbases
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // type can be: graph, multigraph, hist, profile with corresponding option opt
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::Draw(3)------" << endl;

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

// Get trees
   TTree *tree1 = this->GetTree(treename1);
   if (!tree1) return perrGetTree;

   TTree *tree2 = this->GetTree(treename2);
   if (!tree2) return perrGetTree;

   TTree *tree3 = this->GetTree(treename3);
   if (!tree3) return perrGetTree;

// Get logbase for each dimension
   Int_t *base = new Int_t[3]; 
   base[0] = -1;
   base[1] = -1;
   base[2] = -1;

   Int_t dim  = 0;
//ccc char memory problem??
   char *log  = new char[strlen(logbases) + 1];
   char *dlog = log;
   log = strtok(strcpy(log,logbases),":");
//??   log = strtok((char*)logbases,":");
   while(log) {
      if (strcmp(log,"0")     == 0) {base[dim] = 0;}
      if (strcmp(log,"log")   == 0) {base[dim] = 1;}
      if (strcmp(log,"log2")  == 0) {base[dim] = 2;}
      if (strcmp(log,"log10") == 0) {base[dim] = 10;}
      log = strtok(NULL, ":");
      if ((log == 0) || (dim > 2)) break;
      dim++;
   }//while
   delete [] dlog;
   // if logbases contains fewer entries than dimen
   if (base[0] == -1) base[0] = 0;
   if (base[1] == -1) base[1] = base[0];
   if (base[2] == -1) base[2] = base[1];

// Get leafs and variables
   TLeaf    *leaf1 = 0;
   TLeaf    *leaf2 = 0;
   TLeaf    *leaf3 = 0;
   TBranch  *brch1 = 0;
   TBranch  *brch2 = 0;
   TBranch  *brch3 = 0;
   Double_t *arrX  = 0;
   Double_t *arrY  = 0;
   Double_t *arrZ  = 0;
   TString   var1  = "";
   TString   var2  = "";
   TString   var3  = "";
   TString   title = "";

   Int_t entries  = (Int_t)(tree1->GetEntries());

   // get leaf1 for tree1
//ccc char memory problem??
   char *varname  = new char[strlen(varlist) + 1];
   char *dvarname = varname;
   varname = strtok(strcpy(varname,varlist),":");
//??   varname = strtok((char*)varlist,":");
   if (varname != 0) {
      if (!(leaf1 = tree1->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch1 = leaf1->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      if (!(arrX  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 
      var1  = TString(varname);
   } else {
      cerr << "Error: Variable name(s) for trees are missing." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if

   // get leaf2 for tree2
   varname = strtok(NULL, ":");
   if (varname != 0) {var2 = TString(varname);}
   else              {var2 = var1;}
   if (!(leaf2 = tree2->FindLeaf(var2))) {
      cerr << "Error: Leaf <" << var2 << "> not found." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if
   if (!(brch2 = leaf2->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
   if (!(arrY  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 

   // get leaf3 for tree3
   varname = strtok(NULL, ":");
   if (varname != 0) {var3 = TString(varname);}
   else              {var3 = var2;}
   if (!(leaf3 = tree3->FindLeaf(var3))) {
      cerr << "Error: Leaf <" << var3 << "> not found." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if
   if (!(brch3 = leaf3->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
   if (!(arrZ  = new (nothrow) Double_t[entries])) {perr = perrInitMemory; goto cleanup;} 

// Fill arrays
   this->FillArrays(entries, brch1, leaf1, brch2, leaf2, brch3, leaf3,
                    arrX, arrY, arrZ, base);
   if (fNNegX > 0) {
      cout << "Warning: <" << fNNegX << "> data<=0 on X-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if
   if (fNNegY > 0) {
      cout << "Warning: <" << fNNegY << "> data<=0 on Y-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if
   if (fNNegZ > 0) {
      cout << "Warning: <" << fNNegZ << "> data<=0 on Z-axis replaced with <"
           << fNegLog << ">." << endl;
   }//if

// Set axes range
   fMin = (fMinX < fMinY) ? fMinX : fMinY;
   fMin = (fMin  < fMinZ) ? fMin  : fMinZ;
   fMax = (fMaxX > fMaxY) ? fMaxX : fMaxY;
   fMax = (fMax  > fMaxZ) ? fMax  : fMaxZ;
   fEqualAxes = kFALSE;

   if (minX < maxX){
      fMinX = minX;
      fMaxX = maxX;
   } else {
;//      cout << "Warning: minX >= maxX, using computed values" << endl;
   }//if
   if (minY < maxY){
      fMinY = minY;
      fMaxY = maxY;
   } else {
;//      cout << "Warning: minY >= maxY, using computed values" << endl;
   }//if
   if (minZ < maxZ){
      fMinZ = minZ;
      fMaxZ = maxZ;
   } else {
;//      cout << "Warning: minZ >= maxZ, using computed values" << endl;
   }//if

   if ((minX == 0) && (maxX == 0)  &&
       (minY == 0) && (maxY == 0)  &&
       (minZ == 0) && (maxZ == 0)) {
      fMax = fMax + fMax/20;
      fEqualAxes = kTRUE;
   } else {
      if ((minX == -1111) && (maxX == -1111)) fMaxX = fMaxX + fMaxX/20;
      if ((minY == -1111) && (maxY == -1111)) fMaxY = fMaxY + fMaxY/20;
      if ((minZ == -1111) && (maxZ == -1111)) fMaxZ = fMaxZ + fMaxZ/20;
      fEqualAxes = kFALSE;
   }//if

//?? in GUI: display fMinX,fMinY,etc

// Draw type
   title = TString(tree1->GetName()) + " vs "
         + TString(tree2->GetName()) + " vs "
         + TString(tree3->GetName());
   this->SetTitleMain(title.Data(), fSetTitle);
   if (strcmp(type, "graph") == 0) {
      this->SetTitleX(var1, fSetTitleX, base[0]);
      this->SetTitleY(var2, fSetTitleY, base[1]);
      this->SetTitleZ(var3, fSetTitleZ, base[2]);
      this->DrawGraph3D(entries, arrX, arrY, arrZ, opt);
   } else if (strcmp(type, "multigraph") == 0) {
      this->SetTitleX("Index", fSetTitleX, base[0]);
      this->SetTitleY(var1, fSetTitleY, base[1]);
      this->SetTitleZ(var2, fSetTitleZ, base[1]);
      this->DrawMultiGraph(entries, arrX, arrY, arrZ, opt, sort, down);
   } else if (strcmp(type, "hist") == 0) {
      this->SetTitleX(var1, fSetTitleX, base[0]);
      this->SetTitleY(var2, fSetTitleY, base[1]);
      this->SetTitleZ(var3, fSetTitleZ, base[2]);
      this->DrawHist3D(entries, arrX, arrY, arrZ, opt);
   } else if (strcmp(type, "mvaplot") == 0) {
      cout << "Error: Type <" << type << "> not possible for 3D." << endl;
      perr = perrPlotType;
   } else if (strcmp(type, "profile") == 0) {
//to do      this->DrawGraph2D(entries, arrX, arrY, opt);
   } else {
      cout << "Error: Drawing type <" << type << "> not known" << endl;
      perr = perrPlotType;
   }//if

// Cleanup
cleanup:
   if (base) delete [] base;
   if (arrX) delete [] arrX;
   if (arrY) delete [] arrY;
   if (arrZ) delete [] arrZ;

   delete [] dvarname;  //??

   savedir->cd();
   return perr;
}//Draw

//______________________________________________________________________________
Int_t XPlot::DrawImage(const char *canvasname, const char *treename,
             const char *varlist, const char *logbase, Option_t *opt,
             Double_t min, Double_t max, Option_t *orientation)
{
   // Draw image for tree containing data with (x,y)-coordinates
   // varlist must be of the form "fX:fY:fData" (or equivalent leaf names)
   // logbase: 0(linear), log(ln), log10, log2
   // orientation: U(up), D(down), L(rotate left), R(rotate right)
   //    +M: mirror image, e.g. RM (mirror of rotate right)
   if(kCS) cout << "------XPlot::DrawImage------" << endl;

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

// Get trees
   TTree *tree = this->GetTree(treename);
   if (!tree) return perrGetTree;

// Get logbase
   Int_t base = -1; 
   if (strcmp(logbase,"0")          == 0) {base = 0;}
   else if (strcmp(logbase,"log")   == 0) {base = 1;}
   else if (strcmp(logbase,"log2")  == 0) {base = 2;}
   else if (strcmp(logbase,"log10") == 0) {base = 10;}
   else {
      cerr << "Error: Logbase not known." << endl;
      return perrGeneral;
   }//if

// Get tree and leaf(s) for varlist
   TLeaf    *leaf1 = 0;
   TLeaf    *leaf2 = 0;
   TLeaf    *leaf3 = 0;
   TBranch  *brch1 = 0;
   TBranch  *brch2 = 0;
   TBranch  *brch3 = 0;
   Double_t *img   = 0;
   TH2D     *h2    = 0;
   TString   varX  = "";
   TString   varY  = "";
   TString   varZ  = "";
   TString   str   = ""; 
   TString   title = "";
   Int_t     nrows = 0;
   Int_t     ncols = 0;
   Int_t     size  = 0;
   Int_t     x, y, z;
   Double_t  v;

   Int_t entries = (Int_t)(tree->GetEntries());

   // get leaf1 for X-coordinate
//ccc char memory problem??
   char *varname  = new char[strlen(varlist) + 1];
   char *dvarname = varname;
   varname = strtok(strcpy(varname,varlist),":");
//??   varname = strtok((char*)varlist,":");
   if (varname != 0) {
      if (!(leaf1 = tree->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch1 = leaf1->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      varX = TString(varname);
   } else {
      cerr << "Error: Variable names for tree are missing." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if

   // get leaf2 for Y-coordinate
   varname = strtok(NULL, ":");
   if (varname != 0) {
      if (!(leaf2 = tree->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch2 = leaf2->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      varY = TString(varname);
   } else {
      cerr << "Error: Variable for Y-coordinate is missing." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if

   // get leaf3 for data variable
   varname = strtok(NULL, ":");
   if (varname != 0) {
      if (!(leaf3 = tree->FindLeaf(varname))) {
         cerr << "Error: Leaf <" << varname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch3 = leaf3->GetBranch())) {perr = perrGetBrnch; goto cleanup;}
      varZ = TString(varname);
   } else {
      cerr << "Error: Variable for data is missing." << endl;
      perr = perrGetLeaf;
      goto cleanup;
   }//if

// Get number of rows and columns
   for (Int_t i=0; i<entries; i++) { 
      brch1->GetEntry(i);
      x = (Int_t)(leaf1->GetValue());
      ncols = (ncols > x) ? ncols : x;
      brch2->GetEntry(i);
      y = (Int_t)(leaf2->GetValue());
      nrows = (nrows > y) ? nrows : y;
   }//for_i
   nrows = nrows + 1;
   ncols = ncols + 1;
   size  = nrows * ncols;

   if (size > entries) {
      cout << "Warning: Number of tree entries <" << entries << "> is less than "
           << "number of rows x columns <" << nrows * ncols << ">." << endl;
   } else if (size < entries) {
      cerr << "Error: Number of rows x columns <" << nrows * ncols
           << "> is less than number of tree entries <" << entries << ">." << endl;
      perr = perrNumEntries;
      goto cleanup;
   }//if

// Initialize image array
   if (!(img = new (nothrow) Double_t[size])) {perr = perrInitMemory; goto cleanup;} 
   for (Int_t i=0; i<size; i++) img[i] = 0; 

// Fill arrays: for logarithms, replace data <= 0 with fNegLog
   fMinX = 0;
   fMinY = 0;
   fMaxX = ncols;
   fMaxY = nrows;
   fMinZ  = DBL_MAX;
   fMaxZ  = -DBL_MAX;
   fNNegZ = 0;
   if (base == 0) {
      for (Int_t i=0; i<entries; i++) { 
         brch1->GetEntry(i);
         x = (Int_t)(leaf1->GetValue());
         brch2->GetEntry(i);
         y = (Int_t)(leaf2->GetValue());
         brch3->GetEntry(i);
         z = x * nrows + y;
         img[z] = leaf3->GetValue();

         fMinZ = (fMinZ < img[z]) ? fMinZ : img[z];
         fMaxZ = (fMaxZ > img[z]) ? fMaxZ : img[z];
      }//for
   } else if (base == 2) {
      for (Int_t i=0; i<entries; i++) { 
         brch1->GetEntry(i);
         x = (Int_t)(leaf1->GetValue());
         brch2->GetEntry(i);
         y = (Int_t)(leaf2->GetValue());
         brch3->GetEntry(i);
         z = x * nrows + y;
         v = leaf3->GetValue();
         if (v > 0) {img[z] = TMath::Log2(v);}
         else       {img[z] = fNegLog; fNNegZ++; continue;}

         fMinZ = (fMinZ < img[z]) ? fMinZ : img[z];
         fMaxZ = (fMaxZ > img[z]) ? fMaxZ : img[z];
      }//for
   } else if (base == 10) {
      for (Int_t i=0; i<entries; i++) { 
         brch1->GetEntry(i);
         x = (Int_t)(leaf1->GetValue());
         brch2->GetEntry(i);
         y = (Int_t)(leaf2->GetValue());
         brch3->GetEntry(i);
         z = x * nrows + y;
         v = leaf3->GetValue();
         if (v > 0) {img[z] = TMath::Log10(v);}
         else       {img[z] = fNegLog; fNNegZ++; continue;}

         fMinZ = (fMinZ < img[z]) ? fMinZ : img[z];
         fMaxZ = (fMaxZ > img[z]) ? fMaxZ : img[z];
      }//for
   } else if (base == 1) {
      for (Int_t i=0; i<entries; i++) { 
         brch1->GetEntry(i);
         x = (Int_t)(leaf1->GetValue());
         brch2->GetEntry(i);
         y = (Int_t)(leaf2->GetValue());
         brch3->GetEntry(i);
         z = x * nrows + y;
         v = leaf3->GetValue();
         if (v > 0) {img[z] = TMath::Log(v);}
         else       {img[z] = fNegLog; fNNegZ++; continue;}

         fMinZ = (fMinZ < img[z]) ? fMinZ : img[z];
         fMaxZ = (fMaxZ > img[z]) ? fMaxZ : img[z];
      }//for
   }//if

   if (fNNegZ > 0) {
      cout << "Warning: <" << fNNegZ << "> data<=0 replaced with <"
           << fNegLog << ">." << endl;
   }//if

// Set image range
   if (min < max){
      fMin = min;
      fMax = max;
   } else {
      fMin = fMinZ;
      fMax = fMaxZ;
   }//if

//?? here? problem with colors, see fNPixels!!
//??  change fPadNr to fPadCount!
//better? if(fPadNr > 0) {fCanvas->cd(fPadNr);} else {fCanvas->cd(fPadCount);}
   fCanvas->cd(fPadNr);

   gStyle->SetOptStat(0000000);

// Draw image
   title = "Image: " + TString(tree->GetName());
   this->SetTitleMain(title.Data(), fSetTitle);
   this->SetTitleX(fTitleX, fSetTitleX, base);
   this->SetTitleY(fTitleY, fSetTitleY, base);
   str = "H2_" + TString(fCanvas->GetPad(fPadNr)->GetName()); 
   h2 = new TH2D(str.Data(),fTitle,ncols,0,ncols,nrows,0,nrows);
   h2->SetMinimum(fMin);
   h2->SetMaximum(fMax);
   if (strcmp(orientation, "U") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(i+1,nrows-j,img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "D") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(ncols-i,j+1,img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "L") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(j+1, i+1, img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "R") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(nrows-j,ncols-i,img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "UM") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(ncols-i,nrows-j,img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "DM") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(i+1, j+1, img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "LM") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(j+1,ncols-i,img[i*nrows+j]);
         }//for_j
      }//for_i
   } else if (strcmp(orientation, "RM") == 0) {
      for (Int_t i=0; i<ncols; i++) {
         for (Int_t j=0; j<nrows; j++) {
            h2->SetCellContent(nrows-j,i+1,img[i*nrows+j]);
         }//for_j
      }//for_i
   }//if

   h2->SetXTitle(fTitleX);
   h2->SetYTitle(fTitleY);
   h2->GetXaxis()->CenterTitle(kTRUE);
   h2->GetYaxis()->CenterTitle(kTRUE);
   h2->Draw(opt); 

// Cleanup
cleanup:
   if (img) delete [] img;

   delete [] dvarname;  //??

   savedir->cd();
   return perr;
}//DrawImage

//______________________________________________________________________________
Int_t XPlot::DrawTree(const char *canvasname, const char *treename,
             const char *varexp, const char *selection, Option_t *opt)
{
   // Draw variable expressions for tree(s)
   if(kCS) cout << "------XPlot::DrawTree------" << endl;

   if (fAbort) return perrAbort;
   TDirectory *savedir = gDirectory;

   Int_t perr = perrNoErr;

//to do
//use: tree->Draw(xxx);
//take advantage of tree cut???
cout << "Note: Not yet implemented." << endl;

   savedir->cd();
   return perr;
}//DrawTree

//______________________________________________________________________________
Int_t XPlot::DrawEntries(const char *canvasname, const char *leafname, Int_t n,
             Int_t *entrylist, const char *logbase, const char *type,
             Option_t *opt, Double_t min, Double_t max, Int_t sort, Bool_t down)
{
   // Draw leaf with leafname for all trees listed in fTrees and for n entries 
   // stored in array entrylist. 
   // logbase can be: 0(linear), log(ln), log10, log2
   // type can be: multigraph, hist(box), boxplot with corresponding option opt
   // Axis range is given by min and  max:
   //   min < max:  range is given by user 
   //   min = max:  all axes have same range, min and max are calculated
   // sort options: sort = 0   entries not sorted before drawing
   //               sort = -1  each entry is sorted individually before drawing
   //               sort = k+1 entries are sorted for entry k (first entry: k=0)
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::DrawEntries------" << endl;

   if (fAbort) return perrAbort;
   if (!fTrees || (fTrees->GetSize() == 0)) {
      cerr << "Error: Need to add first trees to tree list." << endl;
      return perrGetTree;
   }//if
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

   TTree     *tree  = 0;
   TLeaf     *leaf  = 0;
   TBranch   *brch  = 0;
   Double_t  *index = 0;
   Double_t  *arr   = 0;
   Double_t **table = 0;
   TString    title = "";

// Draw leaf entries
   Int_t numtrees = fTrees->GetSize();
   if (numtrees == 1) {
   // Draw n leaf entries for one tree as one graph/hist (X-axis: n entries)
      // init arrays
      if (!(arr   = new (nothrow) Double_t[n])) {fAbort = kTRUE; perr = perrInitMemory; goto cleanarr;}  
      if (!(index = new (nothrow) Double_t[n])) {fAbort = kTRUE; perr = perrInitMemory; goto cleanarr;} 
      for (Int_t i=0; i<n; i++) index[i] = i + 1;

      // get leaf and branch
      tree = (TTree*)(fTrees->At(0));
      if (!(leaf = tree->FindLeaf(leafname))) {
         cerr << "Error: Leaf <" << leafname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanarr;
      }//if
      if (!(brch = leaf->GetBranch())) {perr = perrGetBrnch; goto cleanarr;}

      // fill arrays
      this->FillEntrylist(n, brch, leaf, entrylist, arr, base);

      // set axes range
      if (min < max) {
         fMin = min;
         fMax = max;
      } else {
         fMin = (fMinX < fMin) ? fMinX : fMin;
         fMax = (fMaxX > fMax) ? fMaxX : fMax;
      }//if

      // draw type
      title = TString(tree->GetName());
      this->SetTitleMain(title.Data(), fSetTitle);
      if (strcmp(type, "multigraph") == 0) {
//?? X-axis should display entry numbers
         this->SetTitleX("Index", fSetTitleX, 0);
         this->SetTitleY(leafname, fSetTitleY, base);
         this->DrawGraph1D(n, index, arr, opt, sort);
      } else if (strcmp(type, "hist") == 0) {
         cout << "Hist is not yet implemented." << endl;
         perr = perrNoImplement;
      } else if (strcmp(type, "boxplot") == 0) {
         cout << "Boxplot is not yet implemented." << endl;
         perr = perrNoImplement;
      } else {
         cout << "Error: Drawing type <" << type << "> not known" << endl;
         perr = perrPlotType;
      }//if

      // delete arrays
   cleanarr:
      delete [] arr;
      delete [] index;
   } else {
   // Draw n leaf entries for numtrees trees as n graph/hist (X-axis: numtrees trees)
      // initialize memory for table
      if (!(table = new (nothrow) Double_t*[n])) {
         fAbort = kTRUE; perr = perrInitMemory; goto cleantab;
      }//if 
      for (Int_t i=0; i<n; i++) {
         table[i] = 0;
         if (!(table[i] = new (nothrow) Double_t[numtrees])) {
            fAbort = kTRUE; perr = perrInitMemory; goto cleantab;
         }//if 
      }//for_i

      // loop over entrylist
      for (Int_t i=0; i<n; i++) {
         perr = this->FillEntry(entrylist[i], leafname, numtrees, table[i], base);
         if (perr != perrNoErr) {fAbort = kTRUE; goto cleantab;}

         fMin = (fMinX < fMin) ? fMinX : fMin;
         fMax = (fMaxX > fMax) ? fMaxX : fMax;
      }//for_i

      // set axes range
      fMinX = fMin;
      fMaxX = fMax;
      if (min < max) {
         fMin = min;
         fMax = max;
      }//if

      // draw type
      title = "Leaf <" + TString(leafname) + "> for <" ;
      title += numtrees;
      title += "> trees";
      this->SetTitleMain(title.Data(), fSetTitle);
      if (strcmp(type, "multigraph") == 0) {
         this->SetTitleX("Index", fSetTitleX, 0);
         this->SetTitleY(leafname, fSetTitleY, base);
         this->DrawMultiGraph(numtrees, n, table, opt, sort, down);
      } else if (strcmp(type, "hist") == 0) {
         cout << "Hist is not yet implemented." << endl;
         perr = perrNoImplement;
      } else if (strcmp(type, "boxplot") == 0) {
         cout << "Boxplot is not yet implemented." << endl;
         perr = perrNoImplement;
      } else {
         cout << "Error: Drawing type <" << type << "> not known" << endl;
         perr = perrPlotType;
      }//if

      // delete table
   cleantab:
      for (Int_t i=0; i<n; i++) {
         if (table[i]) {delete [] table[i]; table[i] = 0;}
      }//for_i
      delete [] table;
   }//if

// prevent user from editing content of canvas
//not possible   fCanvas->SetEditable(kFALSE); //TEST!!!
   savedir->cd();
   return perr;
}//DrawEntries

//______________________________________________________________________________
Int_t XPlot::DrawLeaves(const char *canvasname, const char *leafname,
             const char *logbase, const char *type, Option_t *opt, Double_t min,
             Double_t max, Int_t sort, Bool_t down)
{
   // Draw leaf with leafname for all trees listed in fTrees. 
   // logbase can be: 0(linear), log(ln), log10, log2
   // type can be: multigraph, density, boxplot with corresponding option opt
   // Axis range is given by min and  max:
   //   min < max:  range is given by user 
   //   min = max:  all axes have same range, min and max are calculated
   // sort options: sort = 0   leaves not sorted before drawing
   //               sort = -1  each leaf is sorted individually before drawing
   //               sort = k   leaves are sorted for leaf of tree k
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::DrawLeaves------" << endl;

   if (fAbort) return perrAbort;
   if (!fTrees || (fTrees->GetSize() == 0)) {
      cerr << "Error: Need to add first trees to tree list." << endl;
      return perrGetTree;
   }//if
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
   Int_t bases[] = {0,base};

   TTree     *tree  = (TTree*)(fTrees->At(0));
   TLeaf     *leaf  = 0;
   TBranch   *brch  = 0;
   Double_t  *index = 0;
   Double_t  *arr   = 0;
   Double_t **table = 0;
   TString    title = "";

   Int_t numtrees = fTrees->GetSize();
   Int_t entries  = (Int_t)(tree->GetEntries());
   if (!(index = new (nothrow) Double_t[entries])) {
      fAbort = kTRUE; perr = perrInitMemory; goto cleanup;
   }//if 

// Draw leaf(s)
   if (numtrees == 1) {
   // Draw leaf for one tree
      // init array
      if (!(arr = new (nothrow) Double_t[entries])) {
         fAbort = kTRUE; perr = perrInitMemory; goto cleanup;
      }//if  

      // get leaf and branch
      if (!(leaf = tree->FindLeaf(leafname))) {
         cerr << "Error: Leaf <" << leafname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch = leaf->GetBranch())) {perr = perrGetBrnch; goto cleanup;}

      // fill array
      this->FillArrays(entries, brch, leaf, index, arr, bases);

      // set axes range
      if (min < max) {
         fMin = min;
         fMax = max;
      } else {
         fMin = (fMinY < fMin) ? fMinY : fMin;
         fMax = (fMaxY > fMax) ? fMaxY : fMax;
      }//if

      // draw type
      title = TString(tree->GetName());
      this->SetTitleMain(title.Data(), fSetTitle);
      if (strcmp(type, "multigraph") == 0) {
         this->SetTitleX("Index", fSetTitleX, bases[0]);
         this->SetTitleY(leafname, fSetTitleY, bases[1]);
         if (sort != 0) {sort = down ? -1 : 1;}
         this->DrawGraph1D(entries, index, arr, opt, sort);
      } else if (strcmp(type, "density") == 0) {
         this->SetTitleX("Index", fSetTitleX, bases[0]);
         this->SetTitleY(leafname, fSetTitleY, bases[1]);
//TO DO: kernel as variable!!
         this->DrawDensity(entries, index, arr, 512, opt, "epanechnikov");
      } else if (strcmp(type, "boxplot") == 0) {
//TO DO
         cout << "Boxplot is not yet implemented." << endl;
         perr = perrNoImplement;
      } else {
         cout << "Error: Drawing type <" << type << "> not known" << endl;
         perr = perrPlotType;
      }//if

      // delete array
      delete [] arr;
   } else {
   // Draw leafs for multiple trees
      // initialize memory for table
      if (!(table = new (nothrow) Double_t*[numtrees])) {
         fAbort = kTRUE; perr = perrInitMemory; goto cleanup;
      }//if 
      for (Int_t i=0; i<numtrees; i++) {
         table[i] = 0;
         if (!(table[i] = new (nothrow) Double_t[entries])) {
            fAbort = kTRUE; perr = perrInitMemory; goto cleanup;
         }//if 
      }//for_i

      // loop over trees
      for (Int_t i=0; i<numtrees; i++) {
         tree = (TTree*)(fTrees->At(i));
         if (!(leaf = tree->FindLeaf(leafname))) {
            cerr << "Error: Leaf <" << leafname << "> not found." << endl;
            perr = perrGetLeaf;
            goto cleanup;
         }//if
         if (!(brch = leaf->GetBranch())) {perr = perrGetBrnch; goto cleanup;}

         // fill table
         this->FillArrays(entries, brch, leaf, index, table[i], bases);

         // set axes range
         if (min < max) {
            fMin = min;
            fMax = max;
         } else {
            fMin = (fMinY < fMin) ? fMinY : fMin;
            fMax = (fMaxY > fMax) ? fMaxY : fMax;
         }//if
      }//for_i

      // draw type
      title = "Leaf <" + TString(leafname) + "> for <" ;
      title += numtrees;
      title += "> trees";
      this->SetTitleMain(title.Data(), fSetTitle);
      if (strcmp(type, "multigraph") == 0) {
         this->SetTitleX("Index", fSetTitleX, bases[0]);
         this->SetTitleY(leafname, fSetTitleY, bases[1]);
         this->DrawMultiGraph(entries, numtrees, table, opt, sort, down);
      } else if (strcmp(type, "density") == 0) {
         this->SetTitleX("Index", fSetTitleX, bases[0]);
         this->SetTitleY(leafname, fSetTitleY, bases[1]);
//TO DO: kernel as variable!!
         this->DrawDensity(entries, numtrees, table, 512, opt, "epanechnikov");
      } else if (strcmp(type, "boxplot") == 0) {
//to do
         cout << "Boxplot is not yet implemented." << endl;
         perr = perrNoImplement;
      } else {
         cout << "Error: Drawing type <" << type << "> not known" << endl;
         perr = perrPlotType;
      }//if

      // delete table
      for (Int_t i=0; i<numtrees; i++) {
         if (table[i]) {delete [] table[i]; table[i] = 0;}
      }//for_i
      delete [] table;
   }//if

// Cleanup
cleanup:
   if (index) delete [] index;

   savedir->cd();
   return perr;
}//DrawLeaves

//______________________________________________________________________________
Int_t XPlot::DrawDensity(const char *canvasname, const char *leafname,
             const char *logbase, const char *kernel, Option_t *opt, Int_t npts,
             Double_t min, Double_t max)
{
   // Draw density histogram of leafname for all trees listed in fTrees. 
   // logbase can be: 0(linear), log(ln), log10, log2
   // kernel: epanechnikov, gaussian, rectangular, triangular, biweight, cosine
   // opt: must not include "A"
   // npts: size npts should be a power of two, preferably 512 or above
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::DrawDensity------" << endl;

   if (fAbort) return perrAbort;
   if (!fTrees || (fTrees->GetSize() == 0)) {
      cerr << "Error: Need to add first trees to tree list." << endl;
      return perrGetTree;
   }//if
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
   Int_t bases[] = {0,base};

   fCanvas->cd(fPadNr);

   TMultiGraph *mgraph = new TMultiGraph();

   TH1F   *frame = 0;
   TGraph *graph = 0;
   Int_t   incrl = -1;
   Int_t   incrm = -1;
   Int_t   ls    = 0;
   Int_t   lc    = 0;
   Int_t   ms    = 0;
   Int_t   mc    = 0;
   Int_t   modls = fLineStyles.GetSize();
   Int_t   modlc = fLineColors.GetSize();
   Int_t   modms = fMarkerStyles.GetSize();
   Int_t   modmc = fMarkerColors.GetSize();

   Double_t minX  =  DBL_MAX;
   Double_t minY  =  DBL_MAX;
   Double_t maxX  = -DBL_MAX;
   Double_t maxY  = -DBL_MAX;
   Double_t value =  0;

   TTree   *tree  = (TTree*)(fTrees->At(0));
   TLeaf   *leaf  = 0;
   TBranch *brch  = 0;
   TString  title = "";

   Int_t numtrees = fTrees->GetSize();
   Int_t entries  = (Int_t)(tree->GetEntries());

// Init local arrays
   Double_t *index = 0;
   Double_t *arr   = 0;
   Double_t *wght  = 0;
   Double_t *xden  = 0;
   Double_t *yden  = 0;
   if (!(index = new (nothrow) Double_t[entries])) {fAbort = kTRUE; perr = perrInitMemory; goto cleanup;} 
   if (!(arr   = new (nothrow) Double_t[entries])) {fAbort = kTRUE; perr = perrInitMemory; goto cleanup;}  
   if (!(wght  = new (nothrow) Double_t[entries])) {fAbort = kTRUE; perr = perrInitMemory; goto cleanup;} 
   if (!(xden  = new (nothrow) Double_t[npts]))    {fAbort = kTRUE; perr = perrInitMemory; goto cleanup;} 
   if (!(yden  = new (nothrow) Double_t[npts]))    {fAbort = kTRUE; perr = perrInitMemory; goto cleanup;} 

   for (Int_t i=0; i<entries; i++) wght[i] = 1.0; 
   for (Int_t i=0; i<npts;    i++) xden[i] = yden[i] = 0;

// Loop over trees and fill graphs
   for (Int_t i=0; i<numtrees; i++) {
      tree = (TTree*)(fTrees->At(i));
      if (!(leaf = tree->FindLeaf(leafname))) {
         cerr << "Error: Leaf <" << leafname << "> not found." << endl;
         perr = perrGetLeaf;
         goto cleanup;
      }//if
      if (!(brch = leaf->GetBranch())) {perr = perrGetBrnch; goto cleanup;}

      // fill array
      this->FillArrays(entries, brch, leaf, index, arr, bases);

      TStat::Density(entries, arr, wght, npts, xden, yden, kernel);

      // fill graphs
      graph = new TGraph(npts, xden, yden);
      graph->SetBit(kCanDelete);    //delete from TPad

      if (fPriorityMS == fPriorityMC) {
         ms = (modms > 0 ? (i % modms) : 0);
         mc = (modmc > 0 ? (i % modmc) : 0);
      } else if (fPriorityMS < fPriorityMC) {
         if (i % modms == 0) incrm++;
         ms = (modms > 0 ? (i % modms) : 0);
         mc = (modmc > 0 ? (incrm % modmc) : 0);
      } else {
         if (i % modmc == 0) incrm++;
         ms = (modms > 0 ? (incrm % modms) : 0);
         mc = (modmc > 0 ? (i % modmc) : 0);
      }//if
      graph->SetMarkerStyle(fMarkerStyles.At(ms));
      graph->SetMarkerColor(fMarkerColors.At(mc));

      if (fPriorityLS == fPriorityLC) {
         ls = (modls > 0 ? (i % modls) : 0);
         lc = (modlc > 0 ? (i % modlc) : 0);
      } else if (fPriorityLS < fPriorityLC) {
         if (i % modls == 0) incrl++;
         ls = (modls > 0 ? (i % modls) : 0);
         lc = (modlc > 0 ? (incrl % modlc) : 0);
      } else {
         if (i % modlc == 0) incrl++;
         ls = (modls > 0 ? (incrl % modls) : 0);
         lc = (modlc > 0 ? (i % modlc) : 0);
      }//if
      graph->SetLineStyle(fLineStyles.At(ls));
      graph->SetLineColor(fLineColors.At(lc));

      value = TMath::MinElement(npts, xden);
      minX = (minX > value) ? value : minX;
      value = TMath::MinElement(npts, yden);
      minY = (minY > value) ? value : minY;
      value = TMath::MaxElement(npts, xden);
      maxX = (maxX < value) ? value : maxX;
      value = TMath::MaxElement(npts, yden);
      maxY = (maxY < value) ? value : maxY;

//      mgraph->Add(graph);
      mgraph->Add(graph, opt);
   }//for_i

// Set titles
   title = "Leaf <" + TString(leafname) + "> for <" ;
   title += numtrees;
   title += "> trees";
   this->SetTitleMain(title.Data(), fSetTitle);
   this->SetTitleX(leafname, fSetTitleX, bases[1]);
   this->SetTitleY("Density", fSetTitleY, bases[0]);

// Draw frame
   frame = gPad->DrawFrame(minX - 0.2*minX, minY - 0.2*minY, maxX + 0.2*maxX, maxY + 0.2*maxY);

   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

// Draw multigraph
   mgraph->Draw(opt); 

// Cleanup
cleanup:
   if (index) {delete [] index; index = 0;}
   if (arr)   {delete [] arr;   arr   = 0;}
   if (yden)  {delete [] yden;  yden  = 0;}
   if (xden)  {delete [] xden;  xden  = 0;}
   if (wght)  {delete [] wght;  wght  = 0;}

   savedir->cd();
   return perr;
}//DrawDensity

//______________________________________________________________________________
Int_t XPlot::DrawParallelCoord(const char *canvasname, const char *varlist,
             Double_t min, Double_t max, Bool_t aslog, Bool_t gl, Bool_t can)
{
   // Draw parallel coordinate chart of varlist for all trees listed in fTrees,
   // where varlist must be one leafname for a list of trees, 
   // but should contain more than one leafname for only one tree.
   // varlist: leafnames of a tree separated by colons, e.g. "fInten:fStdDev"
   // gl: indicates if all axes should be drawn at the same (global) scale
   // can: indicates if a candle chart should be drawn
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XPlot::DrawParallelCoord------" << endl;

   if (fAbort) return perrAbort;
   if (!fTrees || (fTrees->GetSize() == 0)) {
      cerr << "Error: Need to add first trees to tree list." << endl;
      return perrGetTree;
   }//if
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

   fCanvas->cd(fPadNr);

   Int_t  ntrees  = fTrees->GetSize();
   TTree *tree    = (TTree*)(fTrees->At(0));
   Int_t  entries = (Int_t)(tree->GetEntries());

// Add tree friends
   TTree  *treek = 0;
   TString kname, alias;
   if (ntrees > 1) {
      kname = Path2Name(tree->GetName(), "", ".");
      tree->SetName(kname);   //name witout extension
      for (Int_t k=1; k<ntrees; k++) { 
         treek = (TTree*)(fTrees->At(k));
         kname = Path2Name(treek->GetName(), "", ".");
         alias = kname + "=" + TString(treek->GetName());

         tree->AddFriend(treek, alias.Data());
      }//for_i
   }//if

   fCanvas->cd(fPadNr);

// Create parallel coord
   TParallelCoord *para = new TParallelCoord(tree, entries);

   // add varlist
   TString varname;
   if (ntrees > 1) {
      for (Int_t k=0; k<ntrees; k++) {
         treek    = (TTree*)(fTrees->At(k));

         // see documentation for TTree::Draw() and O.Couet at RootTalk
         treek->SetEstimate(treek->GetEntries());

         if (aslog) {varname = "log("; varname += treek->GetName();}
         else       {varname = treek->GetName();}
         varname += ".";
         varname += varlist;
         if (aslog) {varname += ")";}

         para->AddVariable(varname);
      }//for_k
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;

      // see documentation for TTree::Draw() and O.Couet at RootTalk
      tree->SetEstimate(tree->GetEntries());

      name = strtok(strcpy(name, varlist), ":");
      while(name) {
         if (aslog) {varname = "log("; varname += name;}
         else       {varname = name;}
         if (aslog) {varname += ")";}

         para->AddVariable(varname);

         name = strtok(0, ":");
         if (name == 0) break;
      }//while

      delete [] dname;
   }//if

   // set default parameters
   para->SetGlobalScale(gl);
   para->SetDotsSpacing(5);
//no   para->SetCandleChart(can);
   if (can) para->SetCandleChart(can);

// Set axis range
   if (min < max){
      para->SetGlobalMin(min);
      para->SetGlobalMax(max);
   }//if

   fMin = para->GetGlobalMin();
   fMax = para->GetGlobalMax();

// Draw parallel coord
   para->Draw();

   savedir->cd();
   return perr;
}//DrawParallelCoord

//______________________________________________________________________________
void XPlot::DrawGraph1D(Int_t n, Double_t *index, Double_t *x, Option_t *opt,
            Int_t sort)
{
   // Draw 1-dimensional graph
   // sort < 0 decreasing; sort > 0 increasing
   if(kCS) cout << "------XPlot::DrawGraph1D------" << endl;

//??  change fPadNr to fPadCount!
//?? better? if(fPadNr > 0 && fPadNr <= fNPads) {fCanvas->cd(fPadNr);} else {fCanvas->cd(fPadCount);}
   fCanvas->cd(fPadNr);

   TH1F *frame = 0;
   frame = gPad->DrawFrame(fMinX, fMinY, fMaxX, fMaxY);
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

   TGraph *graph = 0;
   if (sort == 0) {
      graph = new TGraph(n, index, x);
   } else {
   // sort:
      Int_t    *idx = 0;
      Double_t *arr = 0;
      if (!(idx = new (nothrow) Int_t[n])) return;
      if (!(arr = new (nothrow) Double_t[n])) {delete [] idx; return;}

      Bool_t down = (sort > 0) ? kFALSE : kTRUE;
      TMath::Sort(n, x, idx, down);
      for (Int_t i=0; i<n; i++) arr[i] = x[idx[i]];
      graph = new TGraph(n, index, arr);

      delete [] idx;
      delete [] arr;
   }//if

   graph->SetBit(kCanDelete);    //delete from TPad
   graph->SetMarkerStyle(fMarkerStyles.At(0));
   graph->SetMarkerColor(fMarkerColors.At(0));
   graph->SetLineStyle(fLineStyles.At(0));
   graph->SetLineColor(fLineColors.At(0));
   graph->Draw(opt);  //option "A" (axis) not allowed if axes set with DrawFrame
}//DrawGraph1D

//______________________________________________________________________________
void XPlot::DrawGraph2D(Int_t n, Double_t *x, Double_t *y, Option_t *opt)
{
   // Draw 2-dimensional graph
   if(kCS) cout << "------XPlot::DrawGraph2D------" << endl;

//??  change fPadNr to fPadCount!
//?? better? if(fPadNr > 0 && fPadNr <= fNPads) {fCanvas->cd(fPadNr);} else {fCanvas->cd(fPadCount);}
   fCanvas->cd(fPadNr);

   Double_t minX, minY, maxX, maxY;
   if (fEqualAxes) {
      minX = minY = fMin;
      maxX = maxY = fMax;
   } else {
      minX = fMinX;
      minY = fMinY;
      maxX = fMaxX;
      maxY = fMaxY;
   }//if

   TH1F *frame = 0;
   frame = gPad->DrawFrame(minX, minY, maxX, maxY);
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

   TGraph *graph = 0 ;
   graph = new TGraph(n, x, y);
   graph->SetBit(kCanDelete);    //delete from TPad

   graph->SetMarkerStyle(fMarkerStyles.At(0));
   graph->SetMarkerColor(fMarkerColors.At(0));
   graph->SetLineStyle(fLineStyles.At(0));
   graph->SetLineColor(fLineColors.At(0));
   graph->Draw(opt);
}//DrawGraph2D

//______________________________________________________________________________
void XPlot::DrawGraph3D(Int_t n, Double_t *x, Double_t *y, Double_t *z,
            Option_t *opt)
{
   // Draw 3-dimensional graph
   if(kCS) cout << "------XPlot::DrawGraph3D------" << endl;

//to do
cout << "Note: 3D graph not yet implemented" << endl;
}//DrawGraph3D

//______________________________________________________________________________
void XPlot::DrawMultiGraph(Int_t n, Double_t *x, Double_t *y, Option_t *opt,
            Int_t sort, Bool_t down)
{
   // Draw two variables x and y in 1-dimensional graph
   // sort options: sort = 0   columns not sorted before drawing
   //               sort = -1  each column is sorted individually before drawing
   //               sort = k   columns are sorted for column k (k=1...m)
   if(kCS) cout << "------XPlot::DrawMultiGraph(2)------" << endl;

   fCanvas->cd(fPadNr);

   TH1F *frame = 0;
   Double_t max = fMax;
   Double_t min = fMin;
//??   Double_t max = fMax + fMax/10;
//??   Double_t min = fMin - fMin/10;
   frame = gPad->DrawFrame(0, min, n-1, max);
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

// Graph for x
   TGraph *graphX = new TGraph();
   graphX->SetMarkerStyle(fMarkerStyles.At(0));
   graphX->SetMarkerColor(fMarkerColors.At(0));
   graphX->SetLineStyle(fLineStyles.At(0));
   graphX->SetLineColor(fLineColors.At(0));

// Graph for y
   TGraph *graphY = new TGraph();
   graphY->SetMarkerStyle(fMarkerStyles.At(1));
   graphY->SetMarkerColor(fMarkerColors.At(1));
   graphY->SetLineStyle(fLineStyles.At(1));
   graphY->SetLineColor(fLineColors.At(1));

   Double_t *index = 0;
   if (!(index = new (nothrow) Double_t[n])) return; 
   for (Int_t i=0; i<n; i++) index[i] = i;

   if (sort == 0) {
      graphX->DrawGraph(n, index, x, opt); 
      graphY->DrawGraph(n, index, y, opt); 
   } else {
      Int_t    *idx  = 0;
      Double_t *arrx = 0;
      Double_t *arry = 0;
      if (!(idx  = new (nothrow) Int_t[n])) goto cleanup;
      if (!(arrx = new (nothrow) Double_t[n])) {delete [] idx; goto cleanup;}
      if (!(arry = new (nothrow) Double_t[n])) {delete [] idx; delete [] arrx; goto cleanup;}

      // sort:
      if (sort > 0) {
         if (sort == 1)      TMath::Sort(n, x, idx, down);
         else if (sort == 2) TMath::Sort(n, y, idx, down);
         for (Int_t i=0; i<n; i++) {
            arrx[i] = x[idx[i]];
            arry[i] = y[idx[i]];
         }//for_i
      } else {
         TMath::Sort(n, x, idx);
         for (Int_t i=0; i<n; i++) {
            arrx[i] = x[idx[i]];
         }//for_i
         TMath::Sort(n, y, idx);
         for (Int_t i=0; i<n; i++) {
            arry[i] = y[idx[i]];
         }//for_i
      }//if

      graphX->DrawGraph(n, index, arrx, opt); 
      graphY->DrawGraph(n, index, arry, opt); 

      delete [] idx;
      delete [] arrx;
      delete [] arry;
   }//if

// Cleanup
cleanup:
   if (index) delete [] index;
   delete graphX;  //possible since DrawGraph() creates newgraph!
   delete graphY;
}//DrawMultiGraph

//______________________________________________________________________________
void XPlot::DrawMultiGraph(Int_t n, Double_t *x, Double_t *y,  Double_t *z,
            Option_t *opt, Int_t sort, Bool_t down)
{
   // Draw three variables x, y and z in 1-dimensional graph
   // sort options: sort = 0   columns not sorted before drawing
   //               sort = -1  each column is sorted individually before drawing
   //               sort = k   columns are sorted for column k (k=1...m)
   if(kCS) cout << "------XPlot::DrawMultiGraph(3)------" << endl;

   fCanvas->cd(fPadNr);

   TH1F *frame = 0;
   Double_t max = fMax;
   Double_t min = fMin;
   frame = gPad->DrawFrame(0, min, n-1, max);
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

// Graph for x
   TGraph *graphX = new TGraph();
   graphX->SetMarkerStyle(fMarkerStyles.At(0));
   graphX->SetMarkerColor(fMarkerColors.At(0));
   graphX->SetLineStyle(fLineStyles.At(0));
   graphX->SetLineColor(fLineColors.At(0));

// Graph for y
   TGraph *graphY = new TGraph();
   graphY->SetMarkerStyle(fMarkerStyles.At(1));
   graphY->SetMarkerColor(fMarkerColors.At(1));
   graphY->SetLineStyle(fLineStyles.At(1));
   graphY->SetLineColor(fLineColors.At(1));

// Graph for z
   TGraph *graphZ = new TGraph();
   graphZ->SetMarkerStyle(fMarkerStyles.At(2));
   graphZ->SetMarkerColor(fMarkerColors.At(2));
   graphZ->SetLineStyle(fLineStyles.At(2));
   graphZ->SetLineColor(fLineColors.At(2));

   Double_t *index = 0;
   if (!(index = new (nothrow) Double_t[n])) return; 
   for (Int_t i=0; i<n; i++) index[i] = i;

   if (sort == 0) {
      graphX->DrawGraph(n, index, x, opt); 
      graphY->DrawGraph(n, index, y, opt); 
      graphZ->DrawGraph(n, index, z, opt); 
   } else {
      Int_t    *idx  = 0;
      Double_t *arrx = 0;
      Double_t *arry = 0;
      Double_t *arrz = 0;
      if (!(idx  = new (nothrow) Int_t[n])) goto cleanup;
      if (!(arrx = new (nothrow) Double_t[n]))
         {delete [] idx; goto cleanup;}
      if (!(arry = new (nothrow) Double_t[n]))
         {delete [] idx; delete [] arrx; goto cleanup;}
      if (!(arrz = new (nothrow) Double_t[n]))
         {delete [] idx; delete [] arrx; delete [] arry; goto cleanup;}

      // sort:
      if (sort > 0) {
         if (sort == 1)      TMath::Sort(n, x, idx, down);
         else if (sort == 2) TMath::Sort(n, y, idx, down);
         else if (sort == 3) TMath::Sort(n, z, idx, down);
         for (Int_t i=0; i<n; i++) {
            arrx[i] = x[idx[i]];
            arry[i] = y[idx[i]];
            arrz[i] = z[idx[i]];
         }//for_i
      } else {
         TMath::Sort(n, x, idx, down);
         for (Int_t i=0; i<n; i++) {
            arrx[i] = x[idx[i]];
         }//for_i
         TMath::Sort(n, y, idx, down);
         for (Int_t i=0; i<n; i++) {
            arry[i] = y[idx[i]];
         }//for_i
         TMath::Sort(n, z, idx, down);
         for (Int_t i=0; i<n; i++) {
            arrz[i] = z[idx[i]];
         }//for_i
      }//if

      graphX->DrawGraph(n, index, arrx, opt); 
      graphY->DrawGraph(n, index, arry, opt); 
      graphZ->DrawGraph(n, index, arrz, opt); 

      delete [] idx;
      delete [] arrx;
      delete [] arry;
      delete [] arrz;
   }//if

// Cleanup
cleanup:
   if (index) delete [] index;
   delete graphX;
   delete graphY;
   delete graphZ;
}//DrawMultiGraph

//______________________________________________________________________________
void XPlot::DrawMultiGraph(Int_t n, Int_t m, Double_t **table, Option_t *opt,
            Int_t sort, Bool_t down)
{
   // Draw m variables with n entries, stored in table[n,m]
   // opt: PaintGraph options, but do not use "A"
   // sort options: sort = 0   columns not sorted before drawing
   //               sort = -1  each column is sorted individually before drawing
   //               sort = k   columns are sorted for column k (k=1...m)
   if(kCS) cout << "------XPlot::DrawMultiGraph(m)------" << endl;

   fCanvas->cd(fPadNr);

   TH1F *frame = 0;
   Double_t max = fMax;
   Double_t min = fMin;
   frame = gPad->DrawFrame(0, min, n-1, max);
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

   Double_t *index = 0;
   if (!(index = new (nothrow) Double_t[n])) return; 
   for (Int_t i=0; i<n; i++) index[i] = i;
//?? better?   for (Int_t i=0; i<n; i++) index[i] = i+1;

   TGraph *graph = 0;
   Int_t   incrl = -1;
   Int_t   incrm = -1;
   Int_t   ls    = 0;
   Int_t   lc    = 0;
   Int_t   ms    = 0;
   Int_t   mc    = 0;
   Int_t   modls = fLineStyles.GetSize();
   Int_t   modlc = fLineColors.GetSize();
   Int_t   modms = fMarkerStyles.GetSize();
   Int_t   modmc = fMarkerColors.GetSize();

   if (sort == 0) {
      for (Int_t i=0; i<m; i++) { 
         graph = new TGraph(n, index, table[i]);
         graph->SetBit(kCanDelete);    //delete from TPad

         if (fPriorityMS == fPriorityMC) {
            ms = (modms > 0 ? (i % modms) : 0);
            mc = (modmc > 0 ? (i % modmc) : 0);
         } else if (fPriorityMS < fPriorityMC) {
            if (i % modms == 0) incrm++;
            ms = (modms > 0 ? (i % modms) : 0);
            mc = (modmc > 0 ? (incrm % modmc) : 0);
         } else {
            if (i % modmc == 0) incrm++;
            ms = (modms > 0 ? (incrm % modms) : 0);
            mc = (modmc > 0 ? (i % modmc) : 0);
         }//if
         graph->SetMarkerStyle(fMarkerStyles.At(ms));
         graph->SetMarkerColor(fMarkerColors.At(mc));

         if (fPriorityLS == fPriorityLC) {
            ls = (modls > 0 ? (i % modls) : 0);
            lc = (modlc > 0 ? (i % modlc) : 0);
         } else if (fPriorityLS < fPriorityLC) {
            if (i % modls == 0) incrl++;
            ls = (modls > 0 ? (i % modls) : 0);
            lc = (modlc > 0 ? (incrl % modlc) : 0);
         } else {
            if (i % modlc == 0) incrl++;
            ls = (modls > 0 ? (incrl % modls) : 0);
            lc = (modlc > 0 ? (i % modlc) : 0);
         }//if
         graph->SetLineStyle(fLineStyles.At(ls));
         graph->SetLineColor(fLineColors.At(lc));

         graph->Draw(opt);
      }//for_i
   } else {
      Int_t    *idx = 0;
      Double_t *arr = 0;
      if (!(idx = new (nothrow) Int_t[n])) goto cleanup;
      if (!(arr = new (nothrow) Double_t[n])) {delete [] idx; goto cleanup;}

      if ((sort > 0) && (sort <= m)) TMath::Sort(n, table[sort-1], idx, down);
      for (Int_t i=0; i<m; i++) { 
         if (sort < 0) TMath::Sort(n, table[i], idx, down);
         for (Int_t j=0; j<n; j++) arr[j] = table[i][idx[j]];

         graph = new TGraph(n, index, arr);
         graph->SetBit(kCanDelete);

         if (fPriorityMS == fPriorityMC) {
            ms = (modms > 0 ? (i % modms) : 0);
            mc = (modmc > 0 ? (i % modmc) : 0);
         } else if (fPriorityMS < fPriorityMC) {
            if (i % modms == 0) incrm++;
            ms = (modms > 0 ? (i % modms) : 0);
            mc = (modmc > 0 ? (incrm % modmc) : 0);
         } else {
            if (i % modmc == 0) incrm++;
            ms = (modms > 0 ? (incrm % modms) : 0);
            mc = (modmc > 0 ? (i % modmc) : 0);
         }//if
         graph->SetMarkerStyle(fMarkerStyles.At(ms));
         graph->SetMarkerColor(fMarkerColors.At(mc));

         if (fPriorityLS == fPriorityLC) {
            ls = (modls > 0 ? (i % modls) : 0);
            lc = (modlc > 0 ? (i % modlc) : 0);
         } else if (fPriorityLS < fPriorityLC) {
            if (i % modls == 0) incrl++;
            ls = (modls > 0 ? (i % modls) : 0);
            lc = (modlc > 0 ? (incrl % modlc) : 0);
         } else {
            if (i % modlc == 0) incrl++;
            ls = (modls > 0 ? (incrl % modls) : 0);
            lc = (modlc > 0 ? (i % modlc) : 0);
         }//if
         graph->SetLineStyle(fLineStyles.At(ls));
         graph->SetLineColor(fLineColors.At(lc));

         graph->Draw(opt);
      }//for_i

      delete [] idx;
      delete [] arr;
   }//if

// Cleanup
cleanup:
   if (index) delete [] index;
}//DrawMultiGraph

//______________________________________________________________________________
void XPlot::DrawDensity(Int_t n, Double_t *index, Double_t *x, Int_t npts,
            Option_t *opt, const char *kernel)
{
   // Draw density histogram
   // Note: for PlotLeaves() opt must include "A", e.g. "AL"
   if(kCS) cout << "------XPlot::DrawDensity------" << endl;

   fCanvas->cd(fPadNr);

   TH1F   *frame = 0;
   TGraph *graph = 0;

// Init local arrays
   Double_t value = 0;
   Double_t *arr  = 0;
   Double_t *xden = 0;
   Double_t *yden = 0;
   if (!(arr  = new (nothrow) Double_t[n]))    goto cleanup;
   if (!(xden = new (nothrow) Double_t[npts])) goto cleanup;
   if (!(yden = new (nothrow) Double_t[npts])) goto cleanup;

   for (Int_t i=0; i<n;    i++) arr[i]   = x[i];
   for (Int_t i=0; i<n;    i++) index[i] = 1.0; // use index as weight = 1
   for (Int_t i=0; i<npts; i++) xden[i]  = yden[i] = 0;
  
   TStat::Density(n, arr, index, npts, xden, yden, kernel);

   fMinX = fMinY = DBL_MAX;
   fMaxX = fMaxY = -DBL_MAX;
   value = TMath::MinElement(npts, xden);
   fMinX = (fMinX > value) ? value : fMinX;
   value = TMath::MinElement(npts, yden);
   fMinY = (fMinY > value) ? value : fMinY;
   value = TMath::MaxElement(npts, xden);
   fMaxX = (fMaxX < value) ? value : fMaxX;
   value = TMath::MaxElement(npts, yden);
   fMaxY = (fMaxY < value) ? value : fMaxY;

   frame = gPad->DrawFrame(fMinX - 0.2*fMinX, fMinY - 0.2*fMinY, fMaxX + 0.2*fMaxX, fMaxY + 0.2*fMaxY);
//cout << "fMinX= " << fMinX << "   fMinY= " << fMinY << "   fMaxX= " << fMaxX << "   fMaxY= " << fMaxY << endl;
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

   graph = new TGraph(npts, xden, yden);

   graph->SetBit(kCanDelete);    //delete from TPad
   graph->SetMarkerStyle(fMarkerStyles.At(0));
   graph->SetMarkerColor(fMarkerColors.At(0));
   graph->SetLineStyle(fLineStyles.At(0));
   graph->SetLineColor(fLineColors.At(0));
   graph->Draw(opt); 

// Cleanup
cleanup:
   if (yden) {delete [] yden; yden = 0;}
   if (xden) {delete [] xden; xden = 0;}
   if (arr)  {delete [] arr;  arr  = 0;}
}//DrawDensity

//______________________________________________________________________________
void XPlot::DrawDensity(Int_t n, Int_t m, Double_t **table, Int_t npts,
            Option_t *opt, const char *kernel)
{
   // Draw density histograms for multiple leaves
   // Note: for PlotLeaves() opt must include "A", e.g. "AL"
   if(kCS) cout << "------XPlot::DrawDensity(table)------" << endl;

   fCanvas->cd(fPadNr);

   TMultiGraph *mgraph = new TMultiGraph();

   TH1F   *frame = 0;
   TGraph *graph = 0;
   Int_t   incrl = -1;
   Int_t   incrm = -1;
   Int_t   ls    = 0;
   Int_t   lc    = 0;
   Int_t   ms    = 0;
   Int_t   mc    = 0;
   Int_t   modls = fLineStyles.GetSize();
   Int_t   modlc = fLineColors.GetSize();
   Int_t   modms = fMarkerStyles.GetSize();
   Int_t   modmc = fMarkerColors.GetSize();

// Init local arrays
   Double_t value = 0;
   Double_t *wght = 0;
   Double_t *xden = 0;
   Double_t *yden = 0;
   if (!(wght = new (nothrow) Double_t[n]))    goto cleanup;
   if (!(xden = new (nothrow) Double_t[npts])) goto cleanup;
   if (!(yden = new (nothrow) Double_t[npts])) goto cleanup;

   for (Int_t i=0; i<n;    i++) wght[i] = 1.0; 
   for (Int_t i=0; i<npts; i++) xden[i] = yden[i] = 0;
  
   fMinX = fMinY = DBL_MAX;
   fMaxX = fMaxY = -DBL_MAX;
   for (Int_t i=0; i<m; i++) { 
      TStat::Density(n, table[i], wght, npts, xden, yden, kernel);

      graph = new TGraph(npts, xden, yden);
//      graph->SetBit(kCanDelete);    //delete from TPad

      if (fPriorityMS == fPriorityMC) {
         ms = (modms > 0 ? (i % modms) : 0);
         mc = (modmc > 0 ? (i % modmc) : 0);
      } else if (fPriorityMS < fPriorityMC) {
         if (i % modms == 0) incrm++;
         ms = (modms > 0 ? (i % modms) : 0);
         mc = (modmc > 0 ? (incrm % modmc) : 0);
      } else {
         if (i % modmc == 0) incrm++;
         ms = (modms > 0 ? (incrm % modms) : 0);
         mc = (modmc > 0 ? (i % modmc) : 0);
      }//if
      graph->SetMarkerStyle(fMarkerStyles.At(ms));
      graph->SetMarkerColor(fMarkerColors.At(mc));

      if (fPriorityLS == fPriorityLC) {
         ls = (modls > 0 ? (i % modls) : 0);
         lc = (modlc > 0 ? (i % modlc) : 0);
      } else if (fPriorityLS < fPriorityLC) {
         if (i % modls == 0) incrl++;
         ls = (modls > 0 ? (i % modls) : 0);
         lc = (modlc > 0 ? (incrl % modlc) : 0);
      } else {
         if (i % modlc == 0) incrl++;
         ls = (modls > 0 ? (incrl % modls) : 0);
         lc = (modlc > 0 ? (i % modlc) : 0);
      }//if
      graph->SetLineStyle(fLineStyles.At(ls));
      graph->SetLineColor(fLineColors.At(lc));

      value = TMath::MinElement(npts, xden);
      fMinX = (fMinX > value) ? value : fMinX;
      value = TMath::MinElement(npts, yden);
      fMinY = (fMinY > value) ? value : fMinY;
      value = TMath::MaxElement(npts, xden);
      fMaxX = (fMaxX < value) ? value : fMaxX;
      value = TMath::MaxElement(npts, yden);
      fMaxY = (fMaxY < value) ? value : fMaxY;

      mgraph->Add(graph);
//      mgraph->Add(graph, opt);
   }//for_i

   this->ClearPad(fPadNr);

   frame = gPad->DrawFrame(fMinX - 0.2*fMinX, fMinY - 0.2*fMinY, fMaxX + 0.2*fMaxX, fMaxY + 0.2*fMaxY);
//cout << "fMinX= " << fMinX << "   fMinY= " << fMinY << "   fMaxX= " << fMaxX << "   fMaxY= " << fMaxY << endl;
   frame->SetTitle(fTitle);
   frame->SetXTitle(fTitleX);
   frame->SetYTitle(fTitleY);
   frame->GetXaxis()->CenterTitle(kTRUE);
   frame->GetYaxis()->CenterTitle(kTRUE);

   mgraph->Draw(opt); 
//   mgraph->Draw("A"); 

// Cleanup
cleanup:
   if (yden) {delete [] yden; yden = 0;}
   if (xden) {delete [] xden; xden = 0;}
   if (wght) {delete [] wght; wght = 0;}
}//DrawDensity

//______________________________________________________________________________
void XPlot::DrawHist1D(Int_t n, Double_t *index, Double_t *x, Option_t *opt)
{
   // Draw 1-dimensional histogram
   if(kCS) cout << "------XPlot::DrawHist1D------" << endl;

//?? ev. do KologorowTest() between two histograms, if opt = ???
//need to pass arrX->index and arrY->arrX??  or: DrawMultiHist()

   fCanvas->cd(fPadNr);

   TH1D *h1 = 0;
// different names for histograms!! to prevent Warning in <TH1::Build>
   TString str = "H1_" + TString(fCanvas->GetPad(fPadNr)->GetName()); 
   // h1 range on x-axis is (minY,maxY)!! (not (minX,maxX))
   h1 = new TH1D(str.Data(), fTitle, fNBinsX, fMinY, fMaxY);
   for (Int_t i=0; i<n; i++) h1->Fill(x[i]);

//TEST!!
//   gStyle->SetOptStat(0000000);

   h1->SetXTitle(fTitleX);
   h1->SetYTitle(fTitleY);
   h1->GetXaxis()->CenterTitle(kTRUE);
   h1->GetYaxis()->CenterTitle(kTRUE);
   h1->SetMarkerStyle(fMarkerStyles.At(0));
   h1->SetMarkerColor(fMarkerColors.At(0));
   h1->SetLineStyle(fLineStyles.At(0));
   h1->SetLineColor(fLineColors.At(0));
   h1->Draw(opt); 

// prevent user from editing content of canvas
// results in opening of default Canvas if it has multiple pads!!
//NOT ALLOWED!!   gPad->SetEditable(kFALSE); //TEST!!!

//no   delete h1;  //not allowed, h1 deleted by fFile
}//DrawHist1D

//______________________________________________________________________________
void XPlot::DrawHist2D(Int_t n, Double_t *x, Double_t *y, Option_t *opt)
{
   // Draw 2-dimensional histogram
   if(kCS) cout << "------XPlot::DrawHist2D------" << endl;

   fCanvas->cd(fPadNr);

   Double_t minX, minY, maxX, maxY;
   if (fEqualAxes) {
      minX = minY = fMin;
      maxX = maxY = fMax;
   } else {
      minX = fMinX;
      minY = fMinY;
      maxX = fMaxX;
      maxY = fMaxY;
   }//if

   TH2D *h2 = 0;
   TString str = "H2_" + TString(fCanvas->GetPad(fPadNr)->GetName()); 
   h2 = new TH2D(str.Data(), fTitle, fNBinsX, minX, maxX, fNBinsY, minY, maxY);
   for (Int_t i=0; i<n; i++) h2->Fill(x[i], y[i]);

   gStyle->SetOptStat(0000000);

   h2->SetXTitle(fTitleX);
   h2->SetYTitle(fTitleY);
   h2->GetXaxis()->CenterTitle(kTRUE);
   h2->GetYaxis()->CenterTitle(kTRUE);
   h2->SetMarkerStyle(fMarkerStyles.At(0));
   h2->SetMarkerColor(fMarkerColors.At(0));
   h2->SetLineStyle(fLineStyles.At(0));
   h2->SetLineColor(fLineColors.At(0));
   h2->Draw(opt); 
}//DrawHist2D

//______________________________________________________________________________
void XPlot::DrawHist3D(Int_t n, Double_t *x, Double_t *y, Double_t *z,
            Option_t *opt)
{
   // Draw 3-dimensional histogram
   if(kCS) cout << "------XPlot::DrawHist3D------" << endl;

   fCanvas->cd(fPadNr);

   Double_t minX, minY, minZ, maxX, maxY, maxZ;
   if (fEqualAxes) {
      minX = minY = minZ = fMin;
      maxX = maxY = maxZ = fMax;
   } else {
      minX = fMinX;
      minY = fMinY;
      minZ = fMinZ;
      maxX = fMaxX;
      maxY = fMaxY;
      maxZ = fMaxZ;
   }//if

   TH3D *h3 = 0;
   TString str = "H3_" + TString(fCanvas->GetPad(fPadNr)->GetName()); 
   h3 = new TH3D(str.Data(), fTitle, fNBinsX, minX, maxX, fNBinsY, minY, maxY,
                 fNBinsZ, minZ, maxZ);
   for (Int_t i=0; i<n; i++) h3->Fill(x[i], y[i], z[i]);

   h3->SetXTitle(fTitleX);
   h3->SetYTitle(fTitleY);
   h3->SetZTitle(fTitleZ);
   h3->GetXaxis()->CenterTitle(kTRUE);
   h3->GetYaxis()->CenterTitle(kTRUE);
   h3->GetZaxis()->CenterTitle(kTRUE);
   h3->SetMarkerStyle(fMarkerStyles.At(0));
   h3->SetMarkerColor(fMarkerColors.At(0));
   h3->SetLineStyle(fLineStyles.At(0));
   h3->SetLineColor(fLineColors.At(0));
   h3->Draw(opt); 
}//DrawHist3D

//______________________________________________________________________________
void XPlot::DrawMVA(Int_t n, Double_t *x, Double_t *y, Int_t base,
            Option_t *opt)
{
   // Draw M vs A-plot: M = logbase(Y/X) vs A = logbase(X*Y)/2
   if(kCS) cout << "------XPlot::DrawMVA------" << endl;

   Double_t min    = DBL_MAX;
   Double_t max    = -DBL_MAX;
   Int_t    fNNegX = 0;

   Double_t *arrA = 0;
   Double_t *arrM = 0;
   if (!(arrA = new (nothrow) Double_t[n])) goto cleanup;
   if (!(arrM = new (nothrow) Double_t[n])) goto cleanup;

   if (base == 0) {
      for (Int_t i=0; i<n; i++) { 
         arrA[i] = (y[i] * x[i]) / 2;
         if (arrA[i] == 0) {arrM[i] = 0; fNNegX++; continue;}

         arrM[i] = (y[i] / x[i] >= 1) ? (y[i] / x[i]) : (-x[i] / y[i]);
         min = (min < arrM[i]) ? min : arrM[i];
         max = (max > arrM[i]) ? max : arrM[i];
      }//for
   } else {
      for (Int_t i=0; i<n; i++) { 
         arrA[i] = (y[i] + x[i]) / 2;
         arrM[i] = (y[i] - x[i]);
         min = (min < arrM[i]) ? min : arrM[i];
         max = (max > arrM[i]) ? max : arrM[i];
      }//for
   }//if

   if (fNNegX > 0) {
      cout << "Warning: For <" << fNNegX 
           << "> data A=0 (x or y=0) M is replaced with <0>." << endl;
   }//if

   fEqualAxes = kFALSE;  // since x: x*y and y: y/x

   if (fMinX >= fMaxX){
      if (base == 0) {
//         fMinX = fMinX;
         fMaxX = (fMaxX * fMaxX) / 2;
      } else {
//         fMinX = fMinX;
         fMaxX = (fMaxX + fMaxX) / 2;
      }//if
   }//if
   if (fMinY >= fMaxY){
      fMaxY = (max > TMath::Abs(min)) ?  max : TMath::Abs(min);
      fMinY = (- fMaxY);
   }//if

   this->DrawGraph2D(n, arrA, arrM, opt);

cleanup:
   if (arrA) delete [] arrA;
   if (arrM) delete [] arrM;
}//DrawMVA

//______________________________________________________________________________
TFile *XPlot::OpenFile(const char *name, Option_t *option, Bool_t &isOwner)
{
   // Open root file
   if(kCS) cout << "------XPlot::OpenFile------" << endl;

   isOwner = kFALSE;

   // convert option toupper
   TString opt = TString(option);
   opt.ToUpper();

// Abort if option is RECREATE file
   if (strcmp(opt.Data(), "RECREATE") == 0) {
      cerr << "Error: Trying to recreate file <" << name << ">" << endl;
      return 0;
   }//if

   TFile *file = 0;
//ccc char memory problem??
   char *fname;
   if ((fname = gSystem->ExpandPathName(name))) {
      file = gROOT->GetFile(fname);
      if (!file) {
         file = TFile::Open(name, opt.Data());
         isOwner = kTRUE;
      }//if

      delete [] fname;
   }//if

   if (!file || file->IsZombie()) {
      fAbort = kTRUE;
   } else if (file->IsOpen()) {
//TO DO: if (verbose) {
//      cout << "Opening file <" << name << "> in <" << option << "> mode..."
//           << endl;
//TO DO: }//if
      return file;
   }//if

   cerr << "Error: Could not open file <" << name << ">" << endl;
   SafeDelete(file);
   fAbort = kTRUE;
   return 0;
}//OpenFile

//______________________________________________________________________________
Bool_t XPlot::IsOpen(TFile *file, const char *filename)
{
   // Check if root file is already open
   // Note: If new file is identical to old file, a warning will be given that the
   //       file is already open but return value is kFALSE to prevent closing file
   if(kCS) cout << "------XPlot::IsOpen------" << endl;

   if (file) {
      TString oldname = file->GetName();
      TString newname = Path2Name(filename, dSEP, ".") + ".root";

      TString xpaname;
//ccc char memory problem??
      const char *fname;
      if ((fname = gSystem->ExpandPathName(filename))) {
         xpaname = TString(fname);
         delete [] (char*)fname;
      }//if

      if ((strcmp(oldname.Data(), newname.Data()) == 0) ||
          (strcmp(oldname.Data(), xpaname.Data()) == 0)) {
         cout << "Warning: File <" << oldname.Data() << "> is already open." << endl;
         return kFALSE;
      }//if

      return kTRUE;
   }//if

   return kFALSE;
}//IsOpen

//______________________________________________________________________________
TTree *XPlot::GetTree(const char *fullname)
{
   // Get tree from tree fullname 
   // Argument fullname can be: /path/filename.root/treeset/treename.exten
   if(kCS) cout << "------XPlot::GetTree------" << endl;

// Extract tree name
   TString treename = Path2Name(fullname, dSEP, "");
   if (strstr(treename.Data(),".root")) {
      treename = "";
   }//if
   if (strcmp(treename.Data(), "") == 0) {
      cerr << "Error: Treename for tree is missing." << endl;
      return 0;
   }//if

// Extract root filename and get file
   TString filename = "";
   Bool_t isOwner = kFALSE;
   if (strstr(fullname,".root")) {
      filename = Path2Name(fullname,"",".root") + ".root";
      fFile = this->OpenFile(filename, "READ", isOwner);
      if (!fFile) return 0;
      fFile->cd();
   } else if (fFile) {
//?? not necessary?
      filename = fFile->GetName();
   } else {
      cerr << "Error: No open file exists." << endl;
      return 0;
   }//if

// Get name of treeset and change directory
   TString setname  = "";
   if (strstr(fullname,".root")) {
      TString substr = SubString(fullname,'.', sSEP, kFALSE);
      if (substr) setname = Path2Name(substr.Data(), dSEP, "");
      if (setname.Contains("root")) setname = "";
   } else if (strstr(fullname, dSEP)) {
      setname = Path2Name(fullname, "", dSEP);
   }//if

   if (!fFile->cd(setname)) return 0;

   TTree *tree = (TTree*)gDirectory->Get(treename);
   if (!tree) {
      cerr << "Error: Tree <" << fullname << "> not found." << endl;
      return 0;
   }//if

   return tree;
}//GetTree

//______________________________________________________________________________
Int_t XPlot::Open(const char *filename, Option_t *option)
{
   // Open file filename
   if(kCS) cout << "------XPlot::Open------" << endl;

   if (fAbort) return perrAbort;

// Open file
   Bool_t isOwner = kFALSE;
   fFile = this->OpenFile(filename, option, isOwner);
   if (!fFile) {fAbort = kTRUE; return perrOpenFile;}
   // assure that manager remains owner if it calls OpenFile multiple times
   if (!fIsFileOwner) fIsFileOwner = isOwner;

// Change dir to file
   fFile->cd();

   return perrNoErr;
}//Open

//______________________________________________________________________________
Int_t XPlot::AddTree(const char *treename)
{
   // Add tree treename to list of trees
   // Argument intree can be: /path/filename.root/treeset/treename.exten
   // Note: If no root filename is given, the current fFile is used as file
   // For treename is "*.exten", all trees with exten are added from root file
   if(kCS) cout << "------XPlot::AddTree------" << endl;

   if (fAbort) return perrAbort;

// Create new tree list
   if (!fTrees) {
      fTrees = new TList();
      if (fTrees == 0) {fAbort = kTRUE; return perrInitMemory;}
   }//if

// Extract tree name
   TString tname = Path2Name(treename, dSEP, "");
   if (strstr(tname.Data(),".root")) {
      tname = "";
   }//if
   if (strcmp(tname.Data(), "") == 0) {
      cerr << "Error: Tree name is missing." << endl;
      return perrTreeName;
   }//if

// Extract root filename and get file
   TString filename = "";
   Bool_t isOwner = kFALSE;
   if (strstr(treename, ".root")) {
      filename = Path2Name(treename, "", ".root") + ".root";
      fFile = this->OpenFile(filename, "READ", isOwner);
      if (!fFile) return perrOpenFile;
      fFile->cd();
   } else if (fFile) {
//?? not necessary?
      filename = fFile->GetName();
   } else {
      cerr << "Error: No open file exists." << endl;
      return perrOpenFile;
   }//if

// Get name of treeset and change directory
   TString setname  = "";
   if (strstr(treename,".root")) {
//old      TString substr = SubString(treename, '.', '/', kFALSE);
      TString substr = SubString(treename, '.', sSEP, kFALSE);
      if (substr) setname = Path2Name(substr.Data(), dSEP, "");
      if (setname.Contains("root")) setname = "";
   } else if (strstr(treename, dSEP)) {
      setname = Path2Name(treename, "", dSEP);
   }//if

   if (!fFile->cd(setname)) return perrGetDir;

// Add trees to fTrees
   TTree  *tree  = 0;
   TString name  = Path2Name(treename, dSEP, ".");
   TString exten = Path2Name(treename, ".", "");
   if (strcmp(name.Data(),"*") == 0) {
   // Loop over all trees with extension exten
      TKey *key = 0;
      TIter next(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TString xten  = Path2Name(key->GetName(),".",";");
         TString kname = Path2Name(key->GetName(),"",".");
         if (strcmp(xten.Data(), exten) == 0) {
            tree = (TTree*)gDirectory->Get(key->GetName());
            fTrees->Add(tree);
         }//if
      }//while
   } else {
   // Add tree with name tname
      tree = (TTree*)gDirectory->Get(tname);
      fTrees->Add(tree);
   }//if

   return perrNoErr;
}//AddTree

//______________________________________________________________________________
void XPlot::ClearTrees()
{
   // Clear list of trees
   if(kCS) cout << "------XPlot::ClearTrees------" << endl;

   if (fTrees) {
      fTrees->Clear();
      delete fTrees;
      fTrees = 0;
   }//if

   this->SetDefault();
}//ClearTrees

//______________________________________________________________________________
void XPlot::SetDefault()
{
   // Set values to default values
   if(kCS) cout << "------XPlot::SetDefault------" << endl;

//TEST: reset also fAbort:
   fAbort = kFALSE;

   fMin         = fMinX = fMinY = fMinZ = DBL_MAX;
   fMax         = fMaxX = fMaxY = fMaxZ = -DBL_MAX;
   fNBinsX      = fNBinsY = fNBinsZ = 50;
   fNegLog      = 1;
   fNNegX       = fNNegY = fNNegZ = 0;
   fPriorityLC  = 9999999;
   fPriorityLS  = 9999999;
   fPriorityMC  = 9999999;
   fPriorityMS  = 9999999;
   fTitle       = "";
   fTitleX      = "";
   fTitleY      = "";
   fTitleZ      = "";
   fSetTitle    = fSetTitleX = fSetTitleY = fSetTitleZ = 1;
}//SetDefault

//______________________________________________________________________________
void XPlot::SetFillColor(Int_t n, Int_t *colors, Int_t priority)
{
   // Set fill color
   if(kCS) cout << "------XPlot::SetFillColor------" << endl;

//to do
cout << "Note: SetFillColor to do!!" << endl;
}//SetFillColor

//______________________________________________________________________________
void XPlot::SetLineColor(Int_t n, Int_t *colors, Int_t priority)
{
   // Set line color
   // If number of lines/curves in pad is larger than number n of colors
   // then cycle thru colors
   // If priority = 0 then do not cycle thru colors
   // If priority is higher then priority for line style (highest priority = 1)
   // then cycle first thru colors then thru styles
   if(kCS) cout << "------XPlot::SetLineColor------" << endl;

   fPriorityLC = priority;

// Set default line color (at least 3 colors for x,y,z)
   if (n <= 0 && priority == 0) {
      fLineColors.Set(3);
      for (Int_t i=0; i<3; i++) fLineColors.fArray[i] = kBlack;
      return;
   }//if

// Set user-defined line colors
   if (n <= 1) {
      fLineColors.Set(3);
      for (Int_t i=0; i<3; i++) fLineColors.fArray[i] = colors[0];
   } else if (n == 2) {
      fLineColors.Set(3);
      fLineColors.fArray[0] = colors[0];
      fLineColors.fArray[1] = colors[1];
      fLineColors.fArray[2] = colors[0];
   } else if (n > 2) {
      fLineColors.Set(n);
      for (Int_t i=0; i<n; i++) fLineColors.fArray[i] = colors[i];
   }//if
}//SetLineColor

//______________________________________________________________________________
void XPlot::SetLineStyle(Int_t n, Int_t *styles, Int_t priority)
{
   // Set line style
   // If number of lines/curves in pad is larger than number n of styles
   // then cycle thru styles
   // If priority = 0 then do not cycle thru styles
   // If priority is higher then priority for line color (highest priority = 1)
   // then cycle first thru styles then thru colorss
   if(kCS) cout << "------XPlot::SetLineStyle------" << endl;

   fPriorityLS = priority;

// Set default line style (at least 3 styles for x,y,z)
   if (n <= 0 && priority == 0) {
      fLineStyles.Set(3);
      for (Int_t i=0; i<3; i++) fLineStyles.fArray[i] = kSolid;
      return;
   }//if

// Set user-defined line styles
   if (n <= 1) {
      fLineStyles.Set(3);
      for (Int_t i=0; i<3; i++) fLineStyles.fArray[i] = styles[0];
   } else if (n == 2) {
      fLineStyles.Set(3);
      fLineStyles.fArray[0] = styles[0];
      fLineStyles.fArray[1] = styles[1];
      fLineStyles.fArray[2] = styles[0];
   } else if (n > 2) {
      fLineStyles.Set(n);
      for (Int_t i=0; i<n; i++) fLineStyles.fArray[i] = styles[i];
   }//if
}//SetLineStyle

//______________________________________________________________________________
void XPlot::SetMarkerColor(Int_t n, Int_t *colors, Int_t priority)
{
   // Set marker color
   // If number of graphs in pad is larger than number n of colors
   // then cycle thru colors
   // If priority = 0 then do not cycle thru colors
   // If priority is higher then priority for marker style (highest priority = 1)
   // then cycle first thru colors then thru styles
   if(kCS) cout << "------XPlot::SetMarkerColor------" << endl;

   fPriorityMC = priority;

// Set default marker color (at least 3 colors for x,y,z)
   if (n <= 0 && priority == 0) {
      fMarkerColors.Set(3);
      for (Int_t i=0; i<3; i++) fMarkerColors.fArray[i] = kBlack;
      return;
   }//if

// Set user-defined marker colors
   if (n <= 1) {
      fMarkerColors.Set(3);
      for (Int_t i=0; i<3; i++) fMarkerColors.fArray[i] = colors[0];
   } else if (n == 2) {
      fMarkerColors.Set(3);
      fMarkerColors.fArray[0] = colors[0];
      fMarkerColors.fArray[1] = colors[1];
      fMarkerColors.fArray[2] = colors[0];
   } else if (n > 2) {
      fMarkerColors.Set(n);
      for (Int_t i=0; i<n; i++) fMarkerColors.fArray[i] = colors[i];
   }//if
}//SetMarkerColor

//______________________________________________________________________________
void XPlot::SetMarkerStyle(Int_t n, Int_t *styles, Int_t priority)
{
   // Set marker style
   // If number of graphs in pad is larger than number n of styles
   // then cycle thru styles
   // If priority = 0 then do not cycle thru styles
   // If priority is higher then priority for marker color (highest priority = 1)
   // then cycle first thru styles then thru colorss
   if(kCS) cout << "------XPlot::SetMarkerStyle------" << endl;

   fPriorityMS = priority;

// Set default marker style (at least 3 styles for x,y,z)
   if (n <= 0 && priority == 0) {
      fMarkerStyles.Set(3);
      for (Int_t i=0; i<3; i++) fMarkerStyles.fArray[i] = kDot;
      return;
   }//if

// Set user-defined marker styles
   if (n <= 1) {
      fMarkerStyles.Set(3);
      for (Int_t i=0; i<3; i++) fMarkerStyles.fArray[i] = styles[0];
   } else if (n == 2) {
      fMarkerStyles.Set(3);
      fMarkerStyles.fArray[0] = styles[0];
      fMarkerStyles.fArray[1] = styles[1];
      fMarkerStyles.fArray[2] = styles[0];
   } else if (n > 2) {
      fMarkerStyles.Set(n);
      for (Int_t i=0; i<n; i++) fMarkerStyles.fArray[i] = styles[i];
   }//if
}//SetMarkerStyle

//______________________________________________________________________________
void XPlot::SetTitleMain(const char *title, Int_t setTitle)
{
   // Set title of frame: setTitle < 0 auto title, = 0 no title, >0 use "title"
   if(kCS) cout << "------XPlot::SetTitleMain------" << endl;

   if (setTitle == 1) {
      fTitle    = title;
      fSetTitle = -1;
   } else if (setTitle == 0) {
      fTitle    = "";
      fSetTitle = 0;
   } else if (setTitle < 0) {
//      fTitle    = title;
      fSetTitle = 1;
   }//if
}//SetTitleMain

//______________________________________________________________________________
void XPlot::SetTitleX(const char *title, Int_t setTitle, Int_t base)
{
   // Set title of frame: setTitle < 0 auto title, = 0 no title, >0 use "title"
   // 
   if(kCS) cout << "------XPlot::SetTitleX------" << endl;

   if (setTitle == 1) {
      fTitleX    = this->LogTitle(title, base);
      fSetTitleX = -1;
   } else if (setTitle == 0) {
      fTitleX    = "";
      fSetTitleX = 0;
   } else if (setTitle < 0) {
      fSetTitleX = 1;
   }//if

//??   fSetTitleX = setTitle; //??
}//SetTitleX

//______________________________________________________________________________
void XPlot::SetTitleY(const char *title, Int_t setTitle, Int_t base)
{
   // Set title of frame: setTitle < 0 auto title, = 0 no title, >0 use "title"
   // 
   if(kCS) cout << "------XPlot::SetTitleY------" << endl;

   if (setTitle == 1) {
      fTitleY    = this->LogTitle(title, base);
      fSetTitleY = -1;
   } else if (setTitle == 0) {
      fTitleY    = "";
      fSetTitleY = 0;
   } else if (setTitle < 0) {
      fSetTitleY = 1;
   }//if
}//SetTitleY

//______________________________________________________________________________
void XPlot::SetTitleZ(const char *title, Int_t setTitle, Int_t base)
{
   // Set title of frame: setTitle < 0 auto title, = 0 no title, >0 use "title"
   // 
   if(kCS) cout << "------XPlot::SetTitleZ------" << endl;

   if (setTitle == 1) {
      fTitleZ    = this->LogTitle(title, base);
      fSetTitleZ = -1;
   } else if (setTitle == 0) {
      fTitleZ    = "";
      fSetTitleZ = 0;
   } else if (setTitle < 0) {
      fSetTitleZ = 1;
   }//if
}//SetTitleZ

//______________________________________________________________________________
TString XPlot::LogTitle(const char *title, Int_t base)
{
   // Prepend logbase to title
   // 
   if(kCS) cout << "------XPlot::LogTitle------" << endl;

   TString logtitle = "";
//ccc char memory problem?? delete??
   char   *tmpname  = new char[256];
   if (base == 0) {
      sprintf(tmpname, "%s", title);
      logtitle = TString(tmpname);
   } else if (base == 10) {
      sprintf(tmpname, "Log10(%s)", title);
      logtitle = TString(tmpname);
   } else if (base == 2) {
      sprintf(tmpname, "Log2(%s)", title);
      logtitle = TString(tmpname);
   } else if (base == 1) {
      sprintf(tmpname, "Ln(%s)", title);
      logtitle = TString(tmpname);
   }//if

   return logtitle;
}//LogTitle

//______________________________________________________________________________
void XPlot::FillArrays(Int_t n, TBranch *brch1, TLeaf *leaf1, Double_t *index,
            Double_t *arrX, Int_t *base)
{
   // Fill arrays index and arrX: for logarithms, replace data <= 0 with fNegLog
   if(kCS) cout << "------XPlot::FillArrays(1)------" << endl;

// index: x-axis
   if (base[0] == 0) {
      for (Int_t i=0; i<n; i++) index[i] = i + 1;
   } else if (base[0] == 2) {
      for (Int_t i=0; i<n; i++) index[i] = TMath::Log2(i + 1);
   } else if (base[0] == 10) {
      for (Int_t i=0; i<n; i++) index[i] = TMath::Log10(i + 1);
   } else if (base[0] == 1) {
      for (Int_t i=0; i<n; i++) index[i] = TMath::Log(i + 1);
   }//if

   // min/max for index, i.e. x-axis
   fMinX  = index[0];
   fMaxX  = index[n-1];

// arrX: data for y-axis
   Double_t v;
   // min/max for y-axis
   fMinY  = DBL_MAX;
   fMaxY  = -DBL_MAX; //??
   fNNegX = 0;
   if (base[1] == 0) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         arrX[i]  = leaf1->GetValue();
         fMinY = (fMinY < arrX[i]) ? fMinY : arrX[i];
         fMaxY = (fMaxY > arrX[i]) ? fMaxY : arrX[i];
      }//for
   } else if (base[1] == 2) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {arrX[i] = TMath::Log2(v);}
         else       {arrX[i] = fNegLog; fNNegX++; continue;}
         fMinY = (fMinY < arrX[i]) ? fMinY : arrX[i];
         fMaxY = (fMaxY > arrX[i]) ? fMaxY : arrX[i];
      }//for
   } else if (base[1] == 10) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {arrX[i] = TMath::Log10(v);}
         else       {arrX[i] = fNegLog; fNNegX++; continue;}
         fMinY = (fMinY < arrX[i]) ? fMinY : arrX[i];
         fMaxY = (fMaxY > arrX[i]) ? fMaxY : arrX[i];
      }//for
   } else if (base[1] == 1) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {arrX[i] = TMath::Log(v);}
         else       {arrX[i] = fNegLog; fNNegX++; continue;}
         fMinY = (fMinY < arrX[i]) ? fMinY : arrX[i];
         fMaxY = (fMaxY > arrX[i]) ? fMaxY : arrX[i];
      }//for
   }//if
///////////////////////
//TEST for speed:
//BETTER HERE:
//   fMinY = TMath::MinElement(n, arrX);
//   fMaxY = TMath::MaxElement(n, arrX);
//also in other methods and other classes!!!!
////////////////////////
}//FillArrays

//______________________________________________________________________________
void XPlot::FillArrays(Int_t n, TBranch *brch1, TLeaf *leaf1, TBranch *brch2,
            TLeaf *leaf2, Double_t *arrX, Double_t *arrY, Int_t *base)
{
   // Fill arrays arrX and arrY: for logarithms, replace data <= 0 with fNegLog
   if(kCS) cout << "------XPlot::FillArrays(2)------" << endl;

   Double_t v;
   fMinX  = fMinY  = DBL_MAX;
   fMaxX  = fMaxY  = -DBL_MAX; //??
   fNNegX = fNNegY = 0;
   if ((base[0] == 0) && (base[1] == 0)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         arrX[i] = leaf1->GetValue();
         brch2->GetEntry(i);
         arrY[i] = leaf2->GetValue();

         fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
         fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
         fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
      }//for
   } else if ((base[0] == 2) && (base[1] == 2)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {
            arrX[i] = TMath::Log2(v);
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (v > 0) {
            arrY[i] = TMath::Log2(v);
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if
      }//for
   } else if ((base[0] == 10) && (base[1] == 10)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {
            arrX[i] = TMath::Log10(v);
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (v > 0) {
            arrY[i] = TMath::Log10(v);
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if
      }//for
   } else if ((base[0] == 1) && (base[1] == 1)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {
            arrX[i] = TMath::Log(v);
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (v > 0) {
            arrY[i] = TMath::Log(v);
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if
      }//for
   } else {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (base[0] == 0) {
            arrX[i] = v;
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else if (v > 0) {
            if (base[0] == 2)       {arrX[i] = TMath::Log2(v);}
            else if (base[0] == 10) {arrX[i] = TMath::Log10(v);}
            else if (base[0] == 1)  {arrX[i] = TMath::Log(v);}
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (base[1] == 0) {
            arrY[i] = v;
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else if (v > 0) {
            if (base[1] == 2)       {arrY[i] = TMath::Log2(v);}
            else if (base[1] == 10) {arrY[i] = TMath::Log10(v);}
            else if (base[1] == 1)  {arrY[i] = TMath::Log(v);}
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if
      }//for
   }//if
}//FillArrays

//______________________________________________________________________________
void XPlot::FillArrays(Int_t n, TBranch *brch1, TLeaf *leaf1, TBranch *brch2,
            TLeaf *leaf2, TBranch *brch3, TLeaf *leaf3, Double_t *arrX,
            Double_t *arrY, Double_t *arrZ, Int_t *base)
{
   // Fill arrays arrX, arrY, arrZ: for logs, replace data <= 0 with fNegLog
   if(kCS) cout << "------XPlot::FillArrays(3)------" << endl;

   Double_t v;
   fMinX = fMinY = fMinZ = DBL_MAX;
   fMaxX = fMaxY = fMaxZ = -DBL_MAX;  //??
   fNNegX = fNNegY = fNNegZ = 0;
   if ((base[0] == 0) && (base[1] == 0) && (base[2] == 0)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         arrX[i] = leaf1->GetValue();
         brch2->GetEntry(i);
         arrY[i] = leaf2->GetValue();
         brch3->GetEntry(i);
         arrZ[i] = leaf3->GetValue();

         fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
         fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
         fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
         fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
      }//for
   } else if ((base[0] == 2) && (base[1] == 2) && (base[2] == 2)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {
            arrX[i] = TMath::Log2(v);
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (v > 0) {
            arrY[i] = TMath::Log2(v);
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if

         brch3->GetEntry(i);
         v = leaf3->GetValue();
         if (v > 0) {
            arrZ[i] = TMath::Log2(v);
            fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
            fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
         } else {
            arrZ[i] = fNegLog;
            fNNegZ++;
         }//if
      }//for
   } else if ((base[0] == 10) && (base[1] == 10) && (base[2] == 10)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {
            arrX[i] = TMath::Log10(v);
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (v > 0) {
            arrY[i] = TMath::Log10(v);
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if

         brch3->GetEntry(i);
         v = leaf3->GetValue();
         if (v > 0) {
            arrZ[i] = TMath::Log10(v);
            fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
            fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
         } else {
            arrZ[i] = fNegLog;
            fNNegZ++;
         }//if
      }//for
   } else if ((base[0] == 1) && (base[1] == 1) && (base[2] == 1)) {
      for (Int_t i=0; i<n; i++) { 
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (v > 0) {
            arrX[i] = TMath::Log(v);
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (v > 0) {
            arrY[i] = TMath::Log(v);
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if

         brch3->GetEntry(i);
         v = leaf3->GetValue();
         if (v > 0) {
            arrZ[i] = TMath::Log(v);
            fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
            fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
         } else {
            arrZ[i] = fNegLog;
            fNNegZ++;
         }//if
      }//for
   } else {
      for (Int_t i=0; i<n; i++) {
         brch1->GetEntry(i);
         v = leaf1->GetValue();
         if (base[0] == 0) {
            arrX[i] = v;
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else if (v > 0) {
            if (base[0] == 2)       {arrX[i] = TMath::Log2(v);}
            else if (base[0] == 10) {arrX[i] = TMath::Log10(v);}
            else if (base[0] == 1)  {arrX[i] = TMath::Log(v);}
            fMinX = (fMinX < arrX[i]) ? fMinX : arrX[i];
            fMaxX = (fMaxX > arrX[i]) ? fMaxX : arrX[i];
         } else {
            arrX[i] = fNegLog;
            fNNegX++;
         }//if

         brch2->GetEntry(i);
         v = leaf2->GetValue();
         if (base[1] == 0) {
            arrY[i] = v;
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else if (v > 0) {
            if (base[1] == 2)       {arrY[i] = TMath::Log2(v);}
            else if (base[1] == 10) {arrY[i] = TMath::Log10(v);}
            else if (base[1] == 1)  {arrY[i] = TMath::Log(v);}
            fMinY = (fMinY < arrY[i]) ? fMinY : arrY[i];
            fMaxY = (fMaxY > arrY[i]) ? fMaxY : arrY[i];
         } else {
            arrY[i] = fNegLog;
            fNNegY++;
         }//if

         brch3->GetEntry(i);
         v = leaf3->GetValue();
         if (base[1] == 0) {
            arrZ[i] = v;
            fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
            fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
         } else if (v > 0) {
            if (base[1] == 2)       {arrZ[i] = TMath::Log2(v);}
            else if (base[1] == 10) {arrZ[i] = TMath::Log10(v);}
            else if (base[1] == 1)  {arrZ[i] = TMath::Log(v);}
            fMinZ = (fMinZ < arrZ[i]) ? fMinZ : arrZ[i];
            fMaxZ = (fMaxZ > arrZ[i]) ? fMaxZ : arrZ[i];
         } else {
            arrZ[i] = fNegLog;
            fNNegZ++;
         }//if
      }//for
   }//if
}//FillArrays

//______________________________________________________________________________
Int_t XPlot::FillEntry(Int_t entry, const char *leafname, Int_t n, Double_t *arr,
             Int_t base)
{
   // Fill array arr with tree entry from n trees from list fTrees:
   // for logarithms, replace data <= 0 with fNegLog
   if(kCS) cout << "------XPlot::FillEntry------" << endl;

   TTree   *tree = 0;
   TLeaf   *leaf = 0;
   TBranch *brch = 0;

   Double_t v;
   fMinX  = DBL_MAX;
   fMaxX  = -DBL_MAX;
   fNNegX = 0;
   for (Int_t i=0; i<n; i++) {
      tree = (TTree*)(fTrees->At(i));
      if (!(leaf = tree->FindLeaf(leafname))) {
         cerr << "Error: Leaf <" << leafname << "> not found." << endl;
         return perrGetLeaf;
      }//if
      if (!(brch = leaf->GetBranch())) return perrGetLeaf;

      if (base == 0) {
         brch->GetEntry(entry);
         arr[i]  = leaf->GetValue();
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      } else if (base == 2) {
         brch->GetEntry(entry);
         v = leaf->GetValue();
         if (v > 0) {arr[i] = TMath::Log2(v);}
         else       {arr[i] = fNegLog; fNNegX++; continue;}
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      } else if (base == 10) {
         brch->GetEntry(entry);
         v = leaf->GetValue();
         if (v > 0) {arr[i] = TMath::Log10(v);}
         else       {arr[i] = fNegLog; fNNegX++; continue;}
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      } else if (base == 1) {
         brch->GetEntry(entry);
         v = leaf->GetValue();
         if (v > 0) {arr[i] = TMath::Log(v);}
         else       {arr[i] = fNegLog; fNNegX++; continue;}
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      }//if
   }//for_i

   return perrNoErr;
}//FillEntry

//______________________________________________________________________________
void XPlot::FillEntrylist(Int_t n, TBranch *brch, TLeaf *leaf, Int_t *entrylist,
            Double_t *arr, Int_t base)
{
   // Fill array arr with n entries listed in entrylist:
   // for logarithms, replace data <= 0 with fNegLog
   if(kCS) cout << "------XPlot::FillEntrylist------" << endl;

   Double_t v;
   fMinX  = DBL_MAX;
   fMaxX  = -DBL_MAX;
   fNNegX = 0;
   if (base == 0) {
      for (Int_t i=0; i<n; i++) { 
         brch->GetEntry(entrylist[i]);
         arr[i]  = leaf->GetValue();
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      }//for
   } else if (base == 2) {
      for (Int_t i=0; i<n; i++) { 
         brch->GetEntry(entrylist[i]);
         v = leaf->GetValue();
         if (v > 0) {arr[i] = TMath::Log2(v);}
         else       {arr[i] = fNegLog; fNNegX++; continue;}
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      }//for
   } else if (base == 10) {
      for (Int_t i=0; i<n; i++) { 
         brch->GetEntry(entrylist[i]);
         v = leaf->GetValue();
         if (v > 0) {arr[i] = TMath::Log10(v);}
         else       {arr[i] = fNegLog; fNNegX++; continue;}
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      }//for
   } else if (base == 1) {
      for (Int_t i=0; i<n; i++) { 
         brch->GetEntry(entrylist[i]);
         v = leaf->GetValue();
         if (v > 0) {arr[i] = TMath::Log(v);}
         else       {arr[i] = fNegLog; fNNegX++; continue;}
         fMinX = (fMinX < arr[i]) ? fMinX : arr[i];
         fMaxX = (fMaxX > arr[i]) ? fMaxX : arr[i];
      }//for
   }//if
}//FillEntrylist


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// String utility functions                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
Int_t CheckHeader(const char *header, const char **kHeader, const Int_t ncols,
                  Int_t *index, const char *sep)
{
   // Get header line from input as "header" and compare to header columns
   // stored in string array "kHeader". Export "index" array containing
   // 1 for present column and 0 for absent column.

   Int_t i, k;
   Int_t size = 0;
   char *str1;

   char *str = new char[strlen(header) + 1];
   strcpy(str, header);

   index[0] = 1;
   if (strcmp(strtok(str, sep), kHeader[0]) != 0) index[0] = 0;
   for (i=1; i<ncols; i++) {
      str1  = strtok(NULL, sep);
      size += 1 << (ncols - i);
      for (k=1; k<ncols; k++) {
         if (str1 && strcmp(str1, kHeader[k]) == 0) {
            index[k] = 1;
            size -= 1 << (ncols - k);
            break;
         }//if
      }//for_k
   }//for_i

   delete [] str;

   return size;
}//CheckHeader

//______________________________________________________________________________
Int_t CheckHeaderOrder(const char *header, const char **kHeader, const Int_t ncols,
                       Int_t *index, const char *sep)
{
   // Get header line from input as "header" and compare to header columns
   // stored in string array "kHeader".
   // Export "index" array containing column number or -1 for absent column.
   // Returns number of wrong ordered or missing columns

   Int_t i, k;
   Int_t size = ncols;
   char *str1;

   char *str = new char[strlen(header) + 1];
   strcpy(str, header);

// Inititialize index
   for (i=1; i<ncols; i++) index[i] = -1;

   if (strcmp(strtok(str, sep), kHeader[0]) == 0) {
      index[0] = 0;
      size -= 1;
   }//if

   for (i=1; i<ncols; i++) {
      str1  = strtok(NULL, sep);
      for (k=1; k<ncols; k++) {
         if (str1 && strcmp(str1, kHeader[k]) == 0) {
            index[i] = k;
            size -= (i == k) ? 1 : 0;
            break;
         }//if
      }//for_k
   }//for_i

   delete [] str;

   return size;
}//CheckHeaderOrder

//______________________________________________________________________________
Int_t GetHeaderOrder(const char *header, const char **kHeader, const Int_t ncols,
                       Int_t *index, const char *sep)
{
   // Get header line from input as "header" and compare to header columns
   // stored in string array "kHeader".
   // Export "index" array containing column number or -1 for absent column.
   // Returns number of present columns

   char *str = new char[strlen(header) + 1];
   strcpy(str, header);

// Inititialize index
   for (Int_t i=1; i<ncols; i++) index[i] = -1;

// Set index to header order
   Int_t size = 0;
   for (Int_t i=0; i<ncols; i++) {
      char *str1  = (i == 0) ? strtok(str, sep) : strtok(0, sep);
      if (!str1) break;
      for (Int_t k=0; k<ncols; k++) {
         if (strcmp(str1, kHeader[k]) == 0) {
            index[k] = i;
            size++;
            break;
         }//if
      }//for_k
   }//for_i

   delete [] str;

   return size;
}//GetHeaderOrder

//______________________________________________________________________________
TString Extension2Type(const char *exten, const char **types, const char **extens)
{
   // Return type contained in array types for exten of array extens
   // Note: exten can be either one of types or one of extens
   // Note: last entry of arrays types and extens must be the empty string ""

   TString type(exten);

   Int_t i = 0;
   while (strcmp(types[i], "") != 0) {
      if ((strcmp(type, types[i]) == 0) || (strcmp(type, extens[i]) == 0)) {
         type = types[i];
         return type;
      }//if
      i++;
   }//while

//x   return 0;
   return TString(0);
}//Extension2Type

//______________________________________________________________________________
char *FirstPath(const char *name)
{
   Int_t i = 0;
   Int_t k = 0;
   if (name == 0) return 0;
   if (name[0] == sSEP) {
      k = (name[1] == sSEP) ? 2 : 1;
      if (!(i = strcspn(name + k, dSEP))) return 0;
   }//if

/*   char *path = "";
   path = strncat(path,name + k,i);
   return path;
*/
   const char *path = "";
   return strncat((char*)path, name + k, i);
}//FirstPath

//______________________________________________________________________________
TString FirstPath(const char *fullname, char sep)
{
   TString outname(fullname);

   Int_t pos = outname.Last(sep);
   while (pos > 0) {
      outname.Resize(pos);
      pos = outname.Last(sep);
   }//if

   return outname.Strip(TString::kBoth, sep);
}//FirstPath

//______________________________________________________________________________
TString FullName(const char *dir, const char *name, const char *sep)
{
   // Return full name, dir-sep-name, e.g. /path/name
   // Note: name is first stripped of potential path and extension

   TString fullname(dir);
   if (fullname.EndsWith(sep)) {
      fullname = TString(dir) + Path2Name(name, sep, ".");
   } else {
      fullname = TString(dir) + TString(dSEP) + Path2Name(name, sep, ".");
   }//if

//ccc char memory problem??
   char *fname;
   if ((fname = gSystem->ExpandPathName(fullname.Data()))) {
      fullname = TString(fname);
      delete [] fname;
   }//if

   return fullname;
}//FullName

//______________________________________________________________________________
TString GetROOTName(const char *fullname)
{
   // Extract path and root name from fullname

   TString outname(fullname);

   Int_t pos = outname.Index(".root", 5);
   if (pos > 0) {
      outname.Resize(pos);
   } else {
      outname = "";
   }//if

   return outname;
}//GetROOTName

//______________________________________________________________________________
Bool_t HasExtension(const char *exten, Int_t n, const char **extens)
{
   // Return kTRUE if extension exten is contained in array extens

   Bool_t hasExten = kFALSE;

   for (Int_t i = 0; i < n; i++) {
      hasExten = (strcmp(exten, extens[i]) == 0);
      if (hasExten) break;
   }//for_i

   return hasExten;
}//HasExtension

//______________________________________________________________________________
Bool_t HasExtension(const char *exten, const char **extens)
{
   // Return kTRUE if extension exten is contained in array extens
   // Note: last entry of array extens must be the empty string ""

   Bool_t hasExten = kFALSE;

   Int_t i = 0;
   while (strcmp(extens[i], "") != 0) {
      hasExten = (strcmp(exten, extens[i]) == 0);
      if (hasExten) break;
      i++;
   }//while

   return hasExten;
}//HasExtension

//______________________________________________________________________________
TString Name2Path(const char *fullname, char sep)
{
   // Extract path from fullname; sep is path separator

   TString outname(fullname);

   Int_t pos = outname.Last(sep);
   if (pos > 0) {
      outname.Resize(pos);
   } else {
      outname = "";
   }//if

   return outname;
}//Name2Path

//______________________________________________________________________________
Bool_t NameInArray(const char *name, TString *names, Int_t n)
{
   // Return TRUE if name is found in array names.

   for (Int_t i=0; i<n; i++) {
      if (strcmp(name, (names[i]).Data()) == 0) return kTRUE;
   }//for_i

   return kFALSE;
}//NameInArray

//______________________________________________________________________________
Int_t NumSeparators(const char *name, const char *sep)
{
   // Return number of separators in name

   char *tmpname = new char[strlen(name) + 1];
   char *delname = tmpname;

   Int_t idx = 0;
   tmpname = strtok(strcpy(tmpname,name),sep);
   while(tmpname) {
      tmpname = strtok(NULL,sep);
      idx++;
   }//while
   idx--;

   delete [] delname;

   return idx;
}//NumSeparators

//______________________________________________________________________________
Int_t NumSeparators(const char *name, const char sep)
{
   // Return number of separators in name

   Int_t idx = 0;
   char *pos = strchr((char*)name,sep);

   while(pos) {
      pos = strchr(pos+1,sep);
      idx++;
   }//while

   return idx;
}//NumSeparators

//______________________________________________________________________________
TString Path2Name(const char *name, const char *sep, const char *exten)
{
   // Extract name from full path
   // sep is path separator and exten is name extension

   TString outname(name);
   char   *tmpname = new char[strlen(name) + 1];
   char   *delname = tmpname;

   tmpname = strtok(strcpy(tmpname,name),sep);
   while(tmpname) {
      outname = tmpname;
      tmpname = strtok(NULL,sep);
   }//while
      
   if (strcmp(exten,"") != 0) {
      Int_t i = outname.Index(exten);
      if (i > 0) {outname = outname.Remove(i);}
   }//if

   delete [] delname;

   return outname;
}//Path2Name

//______________________________________________________________________________
TString RemoveEnds(const char *name)
{
// Test if name exists
   if (!name || (strlen(name) <= 1)) return name;

   TString outname = name;

// Remove non-alpha chars from start
   Int_t i   = 0;
   Int_t len = strlen(name);
   while((i < len) && (isalpha(outname[i]) == 0) && (isdigit(outname[i]) == 0)) {
      i++;
   }

// Return emtpy string if name has no alphanumeric characters
   if (i == len) return "";

   outname = &outname[i];

// Remove non-alpha chars from end
   i = outname.Length() - 1;
   while((isalpha(outname[i]) == 0) && (isdigit(outname[i]) == 0)) {
      i--;
   }
   outname.Resize(i + 1);

   return outname;
}//RemoveEnds

//______________________________________________________________________________
TString RemoveEnds(const char *name, Int_t &begin, Int_t &end)
{
// Test if name exists
   if (!name || (strlen(name) <= 1)) return name;

   TString outname = name;

// Remove non-alpha chars from start
   Int_t i   = 0;
   Int_t len = strlen(name);
   while((i < len) && (isalpha(outname[i]) == 0) && (isdigit(outname[i]) == 0)) {
      i++;
   }
   // number of leading chars removed
   begin = i;

// Return emtpy string if name has no alphanumeric characters
   if (i == len) return "";

   outname = &outname[i];

// Remove non-alpha chars from end
   len = outname.Length();
   i = len - 1;
   while((isalpha(outname[i]) == 0) && (isdigit(outname[i]) == 0)) {
      i--;
   }
   // number of trailing chars removed
   end = len - i - 1;
   outname.Resize(i + 1);

   return outname;
}//RemoveEnds

//______________________________________________________________________________
TString RemoveLeadingSpace(const char *name, Int_t &begin)
{
// Test if name exists
   if (!name || (strlen(name) <= 1)) return name;

   TString outname = name;

// Remove space chars from start
   Int_t i   = 0;
   Int_t len = strlen(name);
   while((i < len) && (isspace(outname[i]) != 0)) {
      i++;
   }
   // number of leading spaces (space, tab, formfeed)
   begin = i;

// Return emtpy string if name has no alphanumeric characters
   if (i == len) return "";

   outname = &outname[i];
   return outname;
}//RemoveLeadingSpace

//______________________________________________________________________________
TString RemoveSubString(const char *name, const char *substr, Bool_t exact = kTRUE)
{
   // Remove substring from name: if exact=false then ignore case of characters

   char   *tmpname = new char[strlen(name) + 1];
   char   *delname = tmpname;
   Int_t   len     = 0;
   char   *pos;

   tmpname = strcpy(tmpname,name);
   if (exact) {
      pos = strstr(tmpname,substr);
      len = tmpname - pos;
   } else {
      TString upname = TString(name);   upname.ToUpper();
      TString upsubs = TString(substr); upsubs.ToUpper();

      pos = strstr((char*)(upname.Data()),upsubs.Data());
      len = upname.Data() - pos;
   }//if

   TString outname = TString(name);
   if (len < 0) outname.Resize(TMath::Abs(len));
   else         outname = tmpname + strlen(substr);

   delete [] delname;

   return outname;
}//RemoveSubString

//______________________________________________________________________________
TString RemoveNonAlpha(const char *name)
{
   TString outname(name);
   for (Int_t i=0, len=outname.Length(); i<len;) {
      if (!isalnum(outname[i])) {
         outname.Remove(i, 1);
         --len;
      } else {
         ++i;
      }//if
   }//for_i
  
   return outname;
}//RemoveNonAlpha

//______________________________________________________________________________
TString ReplaceNonAlpha(const char *name, const char *sep)
{
   Int_t len = strlen(name);

// Test if name exists
   if (!name || (len <= 1)) return name;

   TString outname(name);
// Replace non-alpha chars with sep
   if (strcmp(sep, "") != 0) {
      for (Int_t i=0; i<len; i++) {
         if(isalnum(outname[i]) == 0) outname.Replace(i, 1, sep, 1);
      }//for_i
   } else {
      outname = RemoveNonAlpha(name);
   }//if

   return outname;
}//ReplaceNonAlpha

//______________________________________________________________________________
Int_t StringInList(const char *str, const char **kStrList, const Int_t n,
                   Bool_t exact = kTRUE)
{
   // Test if str is stored in string array "kStrList" of size n.

   if (str == 0) return 0;

   TString name(str);
   if (exact) {
      for (Int_t i=0; i<n; i++)
         if (name.CompareTo(kStrList[i],TString::kExact) == 0) return 1;
   } else {
      for (Int_t i=0; i<n; i++)
         if (name.CompareTo(kStrList[i],TString::kIgnoreCase) == 0) return 1;
   }//if

   return 0;
}//StringInList

//______________________________________________________________________________
TString SubString(const char *str, const char *sep, Int_t n)
{
   // Extract substring from str at separator number n

   TString outname(str);
   char   *tmpname = new char[strlen(str) + 1];
   char   *delname = tmpname;

   Int_t idx = 0;
   tmpname = strtok(strcpy(tmpname,str),sep);
   while(tmpname) {
      outname = tmpname;
      if (idx == n) break;
      tmpname = strtok(NULL,sep);
      idx++;
   }//while

   delete [] delname;

   return outname;
}//SubString

//______________________________________________________________________________
TString SubString(const char *str, char sep1, char sep2, Bool_t source)
{
   // Extract substring between first sep1 and last sep2
   // For source=true return original string if no sep1 and sep2 in string

   TString outname(str);

   Int_t first = outname.First(sep1) + 1;
   Int_t last  = outname.Last(sep2);
   Int_t len   = last - first;

//x   if (len < 0) return source ? outname : 0;
   if (len < 0) return source ? outname : TString(0);

   outname = &outname[first];
   if (last > 0) outname.Resize(len);

   return outname;
}//SubString

//______________________________________________________________________________
Int_t TokenizeString(const char *cstr, Int_t &n, TString *names, const char *sep)
{
   // Tokenize string into n substrings, separated by delimiter sep 

   Int_t idx = 0;

   names[idx++] = RemoveEnds(strtok((char*)cstr, sep));

   for (Int_t i=1; i<n; i++) {
      char *str = strtok(0, sep);
      if (str == 0) break;
      names[idx++] = str;
   }//for_i

   n = idx;

   return idx;
}//TokenizeString

//______________________________________________________________________________
Int_t TokenizeString(const char *cstr, Int_t &n, TString *names,
      Int_t len, const char *sep)
{
   // Tokenize string into n substrings, separated by delimiter sep of length len 

   Ssiz_t  index = 0;
   Ssiz_t  start = 0;
   Int_t   idx   = 0;
   TString str   = TString(cstr);

   do {
      if (idx == n) break;

      start = index + 1;
      index = str.Index(sep, len, index, TString::kExact);

      if (start == 1) {
         names[idx++] = str(start - 1, (index - start + 1));

         index += 1;
         continue;
      }//if

      if (index < 0) {
         names[idx++] = str(start + 2, str.Length() - start);
         break;
      }//if

      names[idx++] = str(start + 2, (index - start - 2));

      index += 1;
   } while (index > 0);

   n = idx;

   return index;
}//TokenizeString

//______________________________________________________________________________
TString Type2Extension(const char *type, const char **types, const char **extens)
{
   // Return extension contained in array extens for type of array types
   // Note: type can be either one of types or one of extens
   // Note: last entry of arrays types and extens must be the empty string ""

   TString exten(type);

   Int_t i = 0;
   while (strcmp(extens[i], "") != 0) {
      if ((strcmp(type, types[i]) == 0) || (strcmp(type, extens[i]) == 0)) {
         exten = extens[i];
         return exten;
      }//if
      i++;
   }//while

//x   return 0;
   return TString(0);
}//Type2Extension


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Utility functions for file parsing                                   //
// Taken from Affymetrix FILE_PARSERS_SDK                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
inline UShort_t Swap16Bit(UShort_t x) {
  return ((((x) >> 8) & 0xff) | (((x) & 0xff) << 8));
}//Swap16Bit

//______________________________________________________________________________
inline UInt_t Swap32Bit(UInt_t x) {
   return ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) |
           (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24));
}//Swap32Bit

//______________________________________________________________________________
void SwapBytes(char *src, char *dest, Int_t size)
{
	for (Int_t i=0; i<size; i++)
		dest[size - i - 1] = src[i];
}//SwapBytes

//______________________________________________________________________________
void READ_INT(std::ifstream &input, Int_t &value, Bool_t isBE)
{
   READ_UINT(input, (UInt_t&)value, isBE);
}//READ_INT

//______________________________________________________________________________
void READ_UINT(std::ifstream &input, UInt_t &value, Bool_t isBE)
{
   UInt_t val = 0;
   input.read((char*)&val, sizeof(val));

   Bool_t swap = kFALSE;
#ifdef IS_BIG_ENDIAN
   swap = kTRUE;
   if (isBE == kTRUE) swap = kFALSE;
#else
   if (isBE == kTRUE) swap = kTRUE;
#endif
   if (swap) value = Swap32Bit(val);
   else      value = val;
/*   if (swap) {
      char str2[UINT_SIZE];
      SwapBytes((char*)&val, str2, UINT_SIZE);
      memcpy((char*)&val, str2, UINT_SIZE);
   }//if
	memcpy(&value, (char*)&val, UINT_SIZE);
*/
}//READ_UINT

//______________________________________________________________________________
void READ_SHORT(std::ifstream &input, Short_t &value, Bool_t isBE)
{
   READ_USHORT(input, (UShort_t&)value, isBE);
}//READ_SHORT

//______________________________________________________________________________
void READ_USHORT(std::ifstream &input, UShort_t &value, Bool_t isBE)
{
   UShort_t val = 0;
   input.read((char*)&val, sizeof(val));

   Bool_t swap = kFALSE;
#ifdef IS_BIG_ENDIAN
   swap = kTRUE;
   if (isBE == kTRUE) swap = kFALSE;
#else
   if (isBE == kTRUE) swap = kTRUE;
#endif
   if (swap) value = Swap16Bit(val);
   else      value = val;
}//READ_USHORT

//______________________________________________________________________________
void READ_CHAR(std::ifstream &input, char &value)
{
   input.read(&value, CHAR_SIZE);
}//READ_CHAR

//______________________________________________________________________________
void READ_UCHAR(std::ifstream &input, unsigned char &value)
{
   input.read((char *) &value, UCHAR_SIZE);
}//READ_UCHAR

//______________________________________________________________________________
void READ_BOOL(std::ifstream &input, char &value)
{
   input.read(&value, CHAR_SIZE);
}//READ_BOOL

//______________________________________________________________________________
void READ_FLOAT(std::ifstream &input, Float_t &value, Bool_t isBE)
{
   READ_UINT(input, (UInt_t&)value, isBE);
}//READ_FLOAT

//______________________________________________________________________________
void READ_STRING(std::ifstream &input, char * &value, Bool_t isBE)
{
   UInt_t len;
   READ_UINT(input, len, isBE);
   value = new char[len+1];
   if (len > 0)
      input.read(value, len);
   value[len] = '\0';
}//READ_STRING

//______________________________________________________________________________
void READ_STRING(std::ifstream &input, ASTRING *value, Bool_t isBE)
{
   READ_INT(input, value->len, isBE);
   value->value = new char[value->len + 1];
   if (value->len > 0)
      input.read(value->value, value->len);
   value->value[value->len] = '\0';
}//READ_STRING

//______________________________________________________________________________
void READ_WSTRING(std::ifstream &input, char * &value, Bool_t isBE)
{
   UInt_t   len = 0;
   UShort_t val = 0;
   READ_UINT(input, len, isBE);
   value = new char[len+1];
   wchar_t *wstr = new wchar_t[len+1];
   if (len > 0) {
      for (UInt_t i=0; i<len; i++){
         READ_USHORT(input, val, isBE);
         wstr[i] = (wchar_t)val;
      }//for_i
   }//if
   wstr[len] = '\00';
   wcstombs(value, wstr, len+1);
   delete[] wstr; wstr = 0;
}//READ_WSTRING

//______________________________________________________________________________
void READ_WSTRING(std::ifstream &input, wchar_t * &value, Bool_t isBE)
{
   UInt_t   len = 0;
   UShort_t val = 0;
   READ_UINT(input, len, isBE);
   value = new wchar_t[len+1];
   if (len > 0) {
      for (UInt_t i=0; i<len; i++){
         READ_USHORT(input, val, isBE);
         value[i] = (wchar_t)val;
      }//for_i
   }//if
   value[len] = '\00';
}//READ_WSTRING

//______________________________________________________________________________
void READ_WSTRING(std::ifstream &input, AWSTRING *value, Bool_t isBE)
{
   UShort_t val = 0;
   READ_INT(input, value->len, isBE);
   value->value = new wchar_t[value->len + 1];
   if (value->len > 0) {
      for (Int_t i=0; i<value->len; i++){
         input.read((char *)&val, sizeof(UShort_t));
#ifdef IS_BIG_ENDIAN
         value->value[i] = val;
#else
         value->value[i] = Swap16Bit(val);
#endif
      }//for_i
   } else {
      value->value = 0;
   }//if
   value->value[value->len] = '\00';
}//READ_WSTRING

//______________________________________________________________________________
void READ_FIXED_STRING(std::ifstream &input, char *value, Int_t len)
{
   input.read(value, len);
}//READ_FIXED_STRING


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Functions to decode Affymetrix MIME types                            //
// Adapted from Ben Bolstad: package "affyio"                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
wchar_t *DecodeTEXT(ASTRING *value)
{
   Int_t    len    = value->len / sizeof(Short_t);
   wchar_t* result = new wchar_t[len+1];

   ASTRING str;
   str.len   = value->len;
   str.value = new char[value->len+1];
   memcpy(str.value, value->value, value->len);
  
   Short_t *content = (Short_t*)str.value;
   for (Int_t i=0; i<len; i++) {
#ifndef IS_BIG_ENDIAN 
      content[i]=(((content[i]>>8)&0xff) | ((content[i]&0xff)<<8));
#endif
      result[i] = content[i];
   }//for_i
   result[len] = '\00';

   delete[] str.value; str.value = 0;

   return result;
}//DecodeTEXT

//______________________________________________________________________________
Int_t DecodeINT(ASTRING *value)
{
   Int_t result = 0;
   memcpy(&result, value->value, sizeof(Int_t));

#ifndef IS_BIG_ENDIAN 
  result = (((result>>24)&0xff)  | ((result&0xff)<<24) |
		      ((result>>8)&0xff00) | ((result&0xff00)<<8));  
#endif 

   return result;
}//DecodeINT



