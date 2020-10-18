// File created: 08/05/2002                          last modified: 08/14/2009
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

/******************************************************************************
* Major Revision History:
* Nov 2002 - Initial versions finished
* ??? ???? - Add XQuantileNormalizer.
* Jul 2009 - Update XQuantileNormalizer to allow improved ties handling like
+            preprocessCore as option, and to reduce memory consumption
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <new>  //needed for new (nothrow)

#include "TBranch.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphSmooth.h"
#include "TH2.h"
#include "TKey.h"
#include "TLeaf.h"

//for Benchmark test only:
#include <TBenchmark.h>

#include "XPSNormalizer.h"
#include "TStat.h"

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XNormalizer);
ClassImp(XMeanNormalizer);
ClassImp(XMedianNormalizer);
ClassImp(XKernelNormalizer);
ClassImp(XLowessNormalizer);
ClassImp(XSuperNormalizer);
ClassImp(XQuantileNormalizer);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormalizer                                                          //
//                                                                      //
// Class containing algorithms used for normalization                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XNormalizer::XNormalizer()
            :XAlgorithm()
{
   // Default Normalizer constructor
   if(kCS) cout << "---XNormalizer::XNormalizer(default)------" << endl;

   fDataOpt = "";
   fBgrdOpt = "";
   fLogBase = "0";
   fMethod  = "";
   fTies    = "";
   fNApar   = 0;
   fApars   = 0;
   fInitApprox = kFALSE;
}//Constructor

//______________________________________________________________________________
XNormalizer::XNormalizer(const char *name, const char *type)
            :XAlgorithm(name, type)
{
   // Normal Normalizer constructor
   if(kCS) cout << "---XNormalizer::XNormalizer------" << endl;

   fOption  = "transcript"; //default for normalizer
   fDataOpt = "";
   fBgrdOpt = "none";
   fLogBase = "0";
   fMethod  = "";
   fTies    = "";
   fNApar   = 0;
   fApars   = 0;
   fInitApprox = kFALSE;
}//Constructor

//______________________________________________________________________________
XNormalizer::~XNormalizer()
{
   // Normalizer destructor
   if(kCS) cout << "---XNormalizer::~XNormalizer------" << endl;

   if (fApars) {delete [] fApars; fApars = 0;}
   fNApar = 0;
   fInitApprox = kFALSE;
}//Destructor

//______________________________________________________________________________
void XNormalizer::SetOptions(Option_t *opt)
{
   // Set options, e.g. "myoption" or "myoption:log2" or "myoption:correctbg:log2"
   if(kCS) cout << "------XNormalizer::SetOptions------" << endl;

	TString optcpy = opt;

   char *options = (char*)optcpy.Data();
   if (NumSeparators(opt, ":") == 0) {
      fOption  = "transcript"; 
      fDataOpt = "all"; 
      fBgrdOpt = "none";
      fLogBase = strtok(options, ":");
   } else if (NumSeparators(opt, ":") == 1) {
      fOption  = "transcript"; 
      fDataOpt = strtok(options, ":");
      fBgrdOpt = "none";
      fLogBase = strtok(NULL, ":");
   } else if (NumSeparators(opt, ":") == 2) {
      fOption  = strtok(options, ":");
      fDataOpt = strtok(NULL, ":");
      fBgrdOpt = "none";
      fLogBase = strtok(NULL, ":");
   } else {
      fOption  = strtok(options, ":");
      fDataOpt = strtok(NULL, ":");
      fBgrdOpt = strtok(NULL, ":");
      fLogBase = strtok(NULL, ":");
   }//if
}//SetOptions

//______________________________________________________________________________
Int_t XNormalizer::InitApprox(const char *method, const char *ties,
                   Int_t npar, Double_t *pars)
{
   // Initialize normalizer algorithm with approx parameters
   if(kCS) cout << "------XNormalizer::InitApprox------" << endl;

   fMethod = method;
   fTies   = ties;

// since fNApar=0 is allowed
   fInitApprox = kTRUE;

   fNApar = npar;
   if (npar == 0) return errNoErr;
   if (pars == 0) return errNoErr;
   if (!fApars) fApars = new (nothrow) Double_t[fNApar];
   if (!fApars) return errInitMemory;
   memcpy(fApars, pars, npar*sizeof(Double_t));

   return errNoErr;
}//InitApprox

//______________________________________________________________________________
Int_t XNormalizer::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Normalize expression data
   // Input is:
   //    x:   reference array of length n,
   //    y:   current array of length n to be normalized,
   //    msk: mask to be used to select non-varying units
   // Output is:
   //    y:   current array after normalization
   if(kCS) cout << "------XNormalizer::Calculate------" << endl;

   Int_t err = errNoErr;

   TH2F *his2all = 0;
   TH2F *his2sel = 0;

// Get number of masked genes
   Int_t length = 0;        
   for (Int_t i=0; i<n; i++) {
      if (msk[i] == 1) length++;
   }//for_i

   if (length == 0) {
      cerr << "Error: Length of non-varying units is zero." << endl;
      return errAbort;
   }//if

// Initialize local arrays
   Double_t *selX = 0;
   Double_t *selY = 0;
   Double_t *arrX = 0;
   Double_t *arrY = 0;

// Create local arrays
   if (!(selX = new (nothrow) Double_t[length])) {err = errInitMemory; goto cleanup;} 
   if (!(selY = new (nothrow) Double_t[length])) {err = errInitMemory; goto cleanup;} 
   if (!(arrX = new (nothrow) Double_t[n]))      {err = errInitMemory; goto cleanup;} 
   if (!(arrY = new (nothrow) Double_t[n]))      {err = errInitMemory; goto cleanup;} 
   for (Int_t i=0; i<length; i++) {selX[i] = selY[i] = 0.0;}

// Select non-varying genes
   length = 0;        
   for (Int_t i=0; i<n; i++) {
      if (msk[i] == 1) {
         selX[length] = x[i];
         selY[length] = y[i];
         length++;
      }//if
   }//for_i

// Convert arrays to logbase
//?? Fill histograms with log-values
//?? delete histograms?? why so they exist?? drawing only??
   his2all = new TH2F("his2all","H2all",200,1,5,200,1,5);
   his2sel = new TH2F("his2sel","H2sel",200,1,5,200,1,5);
///////////
//PROBLEM!! negative values!!
//but: x, y should no longer contain negative values!
//////////
   if (strcmp(fLogBase.Data(), "0") == 0) {
      for (Int_t i=0; i<n; i++) { 
         arrX[i] = x[i];
         arrY[i] = y[i];
         his2all->Fill(arrX[i], arrY[i]);
      }//for_i
      for (Int_t i=0; i<length; i++) { 
         his2sel->Fill(selX[i], selY[i]);
      }//for_i
   } else if (strcmp(fLogBase.Data(), "log10") == 0) {
      for (Int_t i=0; i<n; i++) { 
         arrX[i] = TMath::Log10(x[i]);
         arrY[i] = TMath::Log10(y[i]);
         his2all->Fill(arrX[i], arrY[i]);
      }//for_i
      for (Int_t i=0; i<length; i++) { 
         selX[i] = TMath::Log10(selX[i]);
         selY[i] = TMath::Log10(selY[i]);
         his2sel->Fill(selX[i], selY[i]);
      }//for_i
   } else if (strcmp(fLogBase.Data(), "log2") == 0) {
      for (Int_t i=0; i<n; i++) { 
         arrX[i] = TMath::Log2(x[i]);
         arrY[i] = TMath::Log2(y[i]);
         his2all->Fill(arrX[i], arrY[i]);
      }//for_i
      for (Int_t i=0; i<length; i++) { 
         selX[i] = TMath::Log2(selX[i]);
         selY[i] = TMath::Log2(selY[i]);
         his2sel->Fill(selX[i], selY[i]);
      }//for_i
   } else if (strcmp(fLogBase.Data(), "log") == 0) {
      for (Int_t i=0; i<n; i++) { 
         arrX[i] = TMath::Log(x[i]);
         arrY[i] = TMath::Log(y[i]);
         his2all->Fill(arrX[i], arrY[i]);
      }//for_i
      for (Int_t i=0; i<length; i++) { 
         selX[i] = TMath::Log(selX[i]);
         selY[i] = TMath::Log(selY[i]);
         his2sel->Fill(selX[i], selY[i]);
      }//for_i
   }//if

// Normalize selected data
   err = this->DoNormalize(length, selX, selY, n, arrX, arrY);

// Fill array y with normalized data (converted from logbase)
   if (strcmp(fLogBase.Data(), "0") == 0) {
      for (Int_t i=0; i<n; i++) y[i] = arrY[i];
   } else if (strcmp(fLogBase.Data(), "log10") == 0) {
      for (Int_t i=0; i<n; i++) y[i] = TMath::Power(10, arrY[i]);
   } else if (strcmp(fLogBase.Data(), "log2") == 0) {
      for (Int_t i=0; i<n; i++) y[i] = TMath::Power(2, arrY[i]);
   } else if (strcmp(fLogBase.Data(), "log") == 0) {
      for (Int_t i=0; i<n; i++) y[i] = TMath::Power(TMath::E(), arrY[i]);
   }//if

// Delete local variables
cleanup:
   if (arrY) {delete [] arrY; arrY = 0;}
   if (arrX) {delete [] arrX; arrX = 0;}
   if (selY) {delete [] selY; selY = 0;}
   if (selX) {delete [] selX; selX = 0;}

//?????
   SafeDelete(his2all);
   SafeDelete(his2sel);

   return err;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMeanNormalizer                                                      //
//                                                                      //
// Class using trimmed mean algorithm for normalization                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMeanNormalizer::XMeanNormalizer()
                :XNormalizer()
{
   // Default MeanNormalizer constructor
   if(kCS) cout << "---XMeanNormalizer::XMeanNormalizer(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XMeanNormalizer::XMeanNormalizer(const char *name, const char *type)
                :XNormalizer(name, type)
{
   // Normal MeanNormalizer constructor
   if(kCS) cout << "---XMeanNormalizer::XMeanNormalizer------" << endl;

}//Constructor

//______________________________________________________________________________
XMeanNormalizer::~XMeanNormalizer()
{
   // MeanNormalizer destructor
   if(kCS) cout << "---XMeanNormalizer::~XMeanNormalizer------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XMeanNormalizer::DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout)
{
   // Normalize array yout by scaling trimmed mean(yin) to trimmed mean(xin)
   // If targetinten > 0, scale mean(yin) to targetinten
   if(kCS) cout << "------XMeanNormalizer::DoNormalize------" << endl;

// Test minimum number of parameters
   if (TestNumParameters(2) != errNoErr) return errInitParameters;

// Parameters
   Double_t trim        = fPars[0];
   Double_t targetinten = fPars[1];

   if (strcmp(fLogBase.Data(), "log10") == 0) {
      targetinten = TMath::Log10(targetinten);
   } else if (strcmp(fLogBase.Data(), "log2") == 0) {
      targetinten = TMath::Log2(targetinten);
   } else if (strcmp(fLogBase.Data(), "log") == 0) {
      targetinten = TMath::Log(targetinten);
   }//if

// Calculate scaling factor
   Double_t sf = 1.0;
   if (strcmp(fDataOpt.Data(), "sel") == 0) {
      if (targetinten > 0) {
         sf = targetinten / TStat::Mean(nin, yin, trim);
      } else {
         sf = TStat::Mean(nin, xin, trim) / TStat::Mean(nin, yin, trim);
      }//if
   } else if (strcmp(fDataOpt.Data(), "all") == 0) {
      if (targetinten > 0) {
         sf = targetinten / TStat::Mean(nout, yout, trim);
      } else {
         sf = TStat::Mean(nout, xout, trim) / TStat::Mean(nout, yout, trim);
      }//if
   } else {
      cerr << "Error: Normalization option <" << fDataOpt.Data() << "> is not known."
           << endl;
      return errAbort;
   }//if

   if (XManager::fgVerbose) {
      cout << "      normalization <Mean>: Scaling factor SF is <" << sf << ">" << endl;
   }//if

   for (Int_t i=0; i<nout; i++) {
      yout[i] = sf * yout[i];
   }//for_i

   return errNoErr;
}//DoNormalize


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMedianNormalizer                                                    //
//                                                                      //
// Class using median as algorithm for normalization                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMedianNormalizer::XMedianNormalizer()
                  :XNormalizer()
{
   // Default MedianNormalizer constructor
   if(kCS) cout << "---XMedianNormalizer::XMedianNormalizer(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XMedianNormalizer::XMedianNormalizer(const char *name, const char *type)
                  :XNormalizer(name, type)
{
   // Normal MedianNormalizer constructor
   if(kCS) cout << "---XMedianNormalizer::XMedianNormalizer------" << endl;

}//Constructor

//______________________________________________________________________________
XMedianNormalizer::~XMedianNormalizer()
{
   // MedianNormalizer destructor
   if(kCS) cout << "---XMedianNormalizer::~XMedianNormalizer------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XMedianNormalizer::DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                         Int_t nout, Double_t *xout, Double_t *yout)
{
   // Normalize array yout by scaling median(yin) to median(xin)
   // If targetinten > 0, scale median(yin) to targetinten
   if(kCS) cout << "------XMedianNormalizer::DoNormalize------" << endl;

// Test minimum number of parameters
   if (TestNumParameters(1) != errNoErr) return errInitParameters;

// Parameters
   Double_t targetinten = fPars[0];

   if (strcmp(fLogBase.Data(), "log10") == 0) {
      targetinten = TMath::Log10(targetinten);
   } else if (strcmp(fLogBase.Data(), "log2") == 0) {
      targetinten = TMath::Log2(targetinten);
   } else if (strcmp(fLogBase.Data(), "log") == 0) {
      targetinten = TMath::Log(targetinten);
   }//if

// Calculate scaling factor
   Double_t sf = 1.0;
   if (strcmp(fDataOpt.Data(), "sel") == 0) {
      if (targetinten > 0) {
         sf = targetinten / TStat::Median(nin, yin);
      } else {
         sf = TStat::Median(nin, xin) / TStat::Median(nin, yin);
      }//if
   } else if (strcmp(fDataOpt.Data(), "all") == 0) {
      if (targetinten > 0) {
         sf = targetinten / TStat::Median(nout, yout);
      } else {
         sf = TStat::Median(nout, xout) / TStat::Median(nout, yout);
      }//if
   } else {
      cerr << "Error: Normalization option <" << fDataOpt.Data() << "> is not known."
           << endl;
      return errAbort;
   }//if

   if (XManager::fgVerbose) {
      cout << "      normalization <Median>: Scaling factor SF is <" << sf << ">" << endl;
   }//if

   for (Int_t i=0; i<nout; i++) {
      yout[i] = sf * yout[i];
   }//for_i

   return errNoErr;
}//DoNormalize


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XKernelNormalizer                                                    //
//                                                                      //
// Class using kernel smoother as algorithm for normalization           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XKernelNormalizer::XKernelNormalizer()
                  :XNormalizer()
{
   // Default KernelNormalizer constructor
   if(kCS) cout << "---XKernelNormalizer::XKernelNormalizer(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XKernelNormalizer::XKernelNormalizer(const char *name, const char *type)
                  :XNormalizer(name, type)
{
   // Normal KernelNormalizer constructor
   if(kCS) cout << "---XKernelNormalizer::XKernelNormalizer------" << endl;

}//Constructor

//______________________________________________________________________________
XKernelNormalizer::~XKernelNormalizer()
{
   // KernelNormalizer destructor
   if(kCS) cout << "---XKernelNormalizer::~XKernelNormalizer------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XKernelNormalizer::DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                         Int_t nout, Double_t *xout, Double_t *yout)
{
   // Normalize array yout by applying kernel smoother on (xin,yin)
   if(kCS) cout << "------XKernelNormalizer::DoNormalize------" << endl;

// Check if Approx is initialized
   if (!fInitApprox) {
      cerr << "Error: InitApprox() was not called! Aborting program."  << endl;
      return errAbort;
   }//if

// Test minimum number of parameters
   if (TestNumParameters(1) != errNoErr) return errInitParameters;

// Parameters
   Double_t bandwidth = fPars[0];
   Int_t    rule      = (Int_t)fApars[0];
   Double_t f         = fApars[1];

// Normalize selected data with kernel smoother
   TGraph *grin, *grout;
   TGraphSmooth *gs = new TGraphSmooth("ksmooth");
   grin = new TGraph(nin, yin, xin); //reference xin as y array!
   grout = gs->SmoothKern(grin, fDataOpt.Data() ,bandwidth);

// Approximate y values for corresponding x values
   TGraph *grnor, *grapp;
   grnor = new TGraph(grout->GetN(), grout->GetX(), grout->GetY());
   grapp = gs->Approx(grnor, fMethod, nout, yout, 0, 0, rule, f, fTies);

   memcpy(yout, grapp->GetY(), nout*sizeof(Double_t));

// Cleanup
   SafeDelete(grnor);
   SafeDelete(grin);
   SafeDelete(gs);

   return errNoErr;
}//DoNormalize


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLowessNormalizer                                                    //
//                                                                      //
// Class using lowess as algorithm for normalization                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XLowessNormalizer::XLowessNormalizer()
                  :XNormalizer()
{
   // Default LowessNormalizer constructor
   if(kCS) cout << "---XLowessNormalizer::XLowessNormalizer(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XLowessNormalizer::XLowessNormalizer(const char *name, const char *type)
                  :XNormalizer(name, type)
{
   // Normal LowessNormalizer constructor
   if(kCS) cout << "---XLowessNormalizer::XLowessNormalizer------" << endl;

}//Constructor

//______________________________________________________________________________
XLowessNormalizer::~XLowessNormalizer()
{
   // LowessNormalizer destructor
   if(kCS) cout << "---XLowessNormalizer::~XLowessNormalizer------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XLowessNormalizer::DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                         Int_t nout, Double_t *xout, Double_t *yout)
{
   // Normalize array yout by applying Lowess smoother on (xin,yin)
   if(kCS) cout << "------XLowessNormalizer::DoNormalize------" << endl;

// Check if Approx is initialized
   if (!fInitApprox) {
      cerr << "Error: InitApprox() was not called! Aborting program."  << endl;
      return errAbort;
   }//if

// Test minimum number of parameters
   if (TestNumParameters(2) != errNoErr) return errInitParameters;

// Parameters
   Double_t span = fPars[0];
   Int_t    iter = (Int_t)fPars[1];
   Int_t    rule = (Int_t)fApars[0];
   Double_t f    = fApars[1];

// Normalize selected data with Lowess smoother
   TGraph *grin, *grout;
   TGraphSmooth *gs = new TGraphSmooth("lowess");
   grin = new TGraph(nin, yin, xin); //reference xin as y array!
   grout = gs->SmoothLowess(grin, fDataOpt.Data() ,span, iter);

// Approximate y values for corresponding x values
   TGraph *grnor, *grapp;
   grnor = new TGraph(grout->GetN(), grout->GetX(), grout->GetY());
   grapp = gs->Approx(grnor, fMethod, nout, yout, 0, 0, rule, f, fTies);

//????????? where is xout?????
   memcpy(yout, grapp->GetY(), nout*sizeof(Double_t));

// Cleanup
   SafeDelete(grnor);
   SafeDelete(grin);
   SafeDelete(gs);

   return errNoErr;
}//DoNormalize


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSuperNormalizer                                                     //
//                                                                      //
// Class using super smoother as algorithm for normalization            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSuperNormalizer::XSuperNormalizer()
                 :XNormalizer()
{
   // Default SuperNormalizer constructor
   if(kCS) cout << "---XSuperNormalizer::XSuperNormalizer(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSuperNormalizer::XSuperNormalizer(const char *name, const char *type)
                 :XNormalizer(name, type)
{
   // Normal SuperNormalizer constructor
   if(kCS) cout << "---XSuperNormalizer::XSuperNormalizer------" << endl;

}//Constructor

//______________________________________________________________________________
XSuperNormalizer::~XSuperNormalizer()
{
   // SuperNormalizer destructor
   if(kCS) cout << "---XSuperNormalizer::~XSuperNormalizer------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XSuperNormalizer::DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                        Int_t nout, Double_t *xout, Double_t *yout)
{
   // Normalize array yout by applying super smoother on (xin,yin)
   if(kCS) cout << "------XSuperNormalizer::DoNormalize------" << endl;

// Check if Approx is initialized
   if (!fInitApprox) {
      cerr << "Error: InitApprox() was not called! Aborting program."  << endl;
      return errAbort;
   }//if

// Test minimum number of parameters
   if (TestNumParameters(2) != errNoErr) return errInitParameters;

// Parameters
   Double_t bass = fPars[0];
   Double_t span = fPars[1];
   Int_t    rule = (Int_t)fApars[0];
   Double_t f    = fApars[1];

// Normalize selected data with super smoother
   TGraph *grin, *grout;
   TGraphSmooth *gs = new TGraphSmooth("supsmu");
   grin = new TGraph(nin, yin, xin); //reference xin as y array!
   grout = gs->SmoothSuper(grin, "" ,bass, span);

// Approximate y values for corresponding x values
   TGraph *grnor, *grapp;
   grnor = new TGraph(grout->GetN(), grout->GetX(), grout->GetY());
   grapp = gs->Approx(grnor, fMethod, nout, yout, 0, 0, rule, f, fTies);

//????????? where is xout?????
   memcpy(yout, grapp->GetY(), nout*sizeof(Double_t));

// Cleanup
   SafeDelete(grnor);
   SafeDelete(grin);
   SafeDelete(gs);

   return errNoErr;
}//DoNormalize


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XQuantileNormalizer                                                  //
//                                                                      //
// Class using quantile algorithm for normalization                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XQuantileNormalizer::XQuantileNormalizer()
                    :XNormalizer()
{
   // Default QuantileNormalizer constructor
   if(kCS) cout << "---XQuantileNormalizer::XQuantileNormalizer(default)------" << endl;

   fNData   = 0;
   fCount   = 0;
   fNMean1  = 0;
   fNMean2  = 0;
   fMean1   = 0;
   fMean2   = 0;
   fTmpFile = 0;
}//Constructor

//______________________________________________________________________________
XQuantileNormalizer::XQuantileNormalizer(const char *name, const char *type)
                    :XNormalizer(name, type)
{
   // Normal QuantileNormalizer constructor
   if(kCS) cout << "---XQuantileNormalizer::XQuantileNormalizer------" << endl;

   fNData   = 0;
   fCount   = 0;
   fNMean1  = 0;
   fNMean2  = 0;
   fMean1   = 0;
   fMean2   = 0;
   fTmpFile = 0;
}//Constructor

//______________________________________________________________________________
XQuantileNormalizer::~XQuantileNormalizer()
{
   // QuantileNormalizer destructor
   if(kCS) cout << "---XQuantileNormalizer::~XQuantileNormalizer------" << endl;

   if (fMean2) {delete [] fMean2; fMean2 = 0;}
   if (fMean1) {delete [] fMean1; fMean1 = 0;}
   SafeDelete(fTmpFile);  //??
}//Destructor

//______________________________________________________________________________
Int_t XQuantileNormalizer::AddArray(Int_t n, Double_t *x, Int_t *msk,
                           const char *name)
{
   // Add content of array x as tree with name to temporary root file
   // Array msk contains flag for PM/MM data to distinguish options:
   // option = "together": store all data in same tree
   // option = "separate": store PM and MM data in separate trees
   // Note: fTmpFile <tmp_rkq.root> will be created on demand
   if(kCS) cout << "---XQuantileNormalizer::AddArray------" << endl;

// Parameter trim (default=0)
   Double_t trim = (fNPar > 0) ? fPars[0] : 0.0;

   TDirectory *savedir = gDirectory;
   Int_t       err     = errNoErr;

// Create temporary root file to store trees
   if (!fTmpFile) {
      // save tmp_rkq.root in same system directory as main file
      TString fname = Path2Name(savedir->GetPath(), "", ".root");
      fname  = Name2Path(fname.Data(), sSEP);
      fname += TString(dSEP) + "tmp_rkq.root";

      fTmpFile = new TFile(fname, "RECREATE");
      if (!fTmpFile || fTmpFile->IsZombie()) {
         cerr << "Error: Could not create temporary file for quantile normalizer."
              << endl;
         SafeDelete(fTmpFile);
         return errCreateFile;
      }//if
   }//if

// Fill tree(s) with sorted data array and store in fTmpFile
   fTmpFile->cd();
   if (strcmp(fDataOpt.Data(), "together") == 0) {
      TString   tname = TString(name) + ".rkt";  //rkt - rank together
      Double_t  sort  = 0;
      Float_t   rank  = 0;
      Int_t     id    = 0;
      Int_t     i     = 0;

   // Init arrays
      Double_t *arr   = 0;
      Int_t    *idx   = 0;
      Double_t *ord   = 0;
      Int_t    *rnki  = 0;
      Double_t *rnkd  = 0;

   // Create tree with data and rank branches
      TTree *tree = new TTree(tname, "temporary tree");
      if (tree == 0) {err = errCreateTree; goto clean1;}
      tree->Branch("sortBr", &sort, "sort/D");
      tree->Branch("rankBr", &rank, "rank/F");

   // Get size of arrays for PM and MM
      for (i=0; i<n; i++) {if (msk[i] == 1) id++;}

   // Create mean array and initialize to 0
      if (!fMean1) {
         fNMean1 = id;
         if (!(fMean1 = new (nothrow) Double_t[id])) {err = errInitMemory; goto clean1;}
         for (Int_t i=0; i<id; i++) fMean1[i] = 0;
      }//if

   // Create arrays
      if (!(arr  = new (nothrow) Double_t[id])) {err = errInitMemory; goto clean1;}
      if (!(idx  = new (nothrow) Int_t[id]))    {err = errInitMemory; goto clean1;}
      if (!(ord  = new (nothrow) Double_t[id])) {err = errInitMemory; goto clean1;}
      if (!(rnki = new (nothrow) Int_t[id]))    {err = errInitMemory; goto clean1;}
      if (!(rnkd = new (nothrow) Double_t[id])) {err = errInitMemory; goto clean1;}

   // Fill arr
      id = 0;
      for (i=0; i<n; i++) {
         if (msk[i] == 1) arr[id++] = x[i];
      }//for_i

   // Sort array arr
      TStat::Rank(id, arr, idx, rnki, kFALSE);
      // get rank rnkd analogoulsy to rank used for RMA in Bioconductor
      for (Int_t i=0; i<id; i++) ord[i] = arr[idx[i]];
      TStat::Rank(id, ord, rnkd);

   // Fill tree with ranks 
      for (i=0; i<id; i++) {
         if (trim == 0.0) fMean1[i] += arr[idx[i]];
         else             sort       = arr[idx[i]];
         rank = rnkd[rnki[i]];
         tree->Fill();
      }//for_i

      tree->Write();
      fNData++;

   clean1:
      // need to delete tree from RAM (otherwise it is used by Calculate()!!)
      if (tree) {tree->Delete(""); tree = 0;}
      if (rnkd) {delete [] rnkd;   rnkd = 0;}
      if (rnki) {delete [] rnki;   rnki = 0;}
      if (ord)  {delete [] ord;    ord  = 0;}
      if (idx)  {delete [] idx;    idx  = 0;}
      if (arr)  {delete [] arr;    arr  = 0;}
   } else if (strcmp(fDataOpt.Data(), "separate") == 0) {
      TString   pname  = TString(name) + ".rkp";  //rks - pm rank separate
      TString   mname  = TString(name) + ".rkm";  //rks - mm rank separate
      Double_t  pmsort = 0;
      Double_t  mmsort = 0;
      Float_t   pmrank = 0;
      Float_t   mmrank = 0;
      Int_t     i      = 0;
      Int_t     p      = 0;
      Int_t     m      = 0;

   // Init arrays
      Double_t *arrPM  = 0;
      Double_t *arrMM  = 0;
      Int_t    *pmidx  = 0;
      Int_t    *mmidx  = 0;
      Double_t *pmord  = 0;
      Double_t *mmord  = 0;
      Int_t    *prnki  = 0;
      Int_t    *mrnki  = 0;
      Double_t *prnkd  = 0;
      Double_t *mrnkd  = 0;

   // Create trees with data and rank branches for PM and MM
      TTree *pmtree = 0;
      TTree *mmtree = 0;

      pmtree = new TTree(pname, "temporary PM tree");
      if (pmtree == 0) {err = errCreateTree; goto clean2;}
      pmtree->Branch("pmsortBr", &pmsort, "pmsort/D");
      pmtree->Branch("pmrankBr", &pmrank, "pmrank/F");

      mmtree = new TTree(mname, "temporary MM tree");
      if (mmtree == 0) {err = errCreateTree; goto clean2;}
      mmtree->Branch("mmsortBr", &mmsort, "mmsort/D");
      mmtree->Branch("mmrankBr", &mmrank, "mmrank/F");

   // Get size of arrays for PM and MM
      for (i=0; i<n; i++) {
         if      (msk[i] == 1) p++;
         else if (msk[i] == 0) m++;
      }//for_i

   // Create mean arrays and initialize to 0
      if (!fMean1) {
         fNMean1 = p;
         if (!(fMean1 = new (nothrow) Double_t[p])) {err = errInitMemory; goto clean2;}
         for (Int_t i=0; i<p; i++) fMean1[i] = 0;
      }//if
      if (!fMean2) {
         fNMean2 = m;
         if (!(fMean2 = new (nothrow) Double_t[m])) {err = errInitMemory; goto clean2;}
         for (Int_t i=0; i<m; i++) fMean2[i] = 0;
      }//if

   // Create arrays
      if (!(arrPM = new (nothrow) Double_t[p])) {err = errInitMemory; goto clean2;}
      if (!(arrMM = new (nothrow) Double_t[m])) {err = errInitMemory; goto clean2;}
      if (!(pmidx = new (nothrow) Int_t[p]))    {err = errInitMemory; goto clean2;}
      if (!(mmidx = new (nothrow) Int_t[m]))    {err = errInitMemory; goto clean2;}
      if (!(pmord = new (nothrow) Double_t[p])) {err = errInitMemory; goto clean2;}
      if (!(mmord = new (nothrow) Double_t[m])) {err = errInitMemory; goto clean2;}
      if (!(prnki = new (nothrow) Int_t[p]))    {err = errInitMemory; goto clean2;}
      if (!(mrnki = new (nothrow) Int_t[m]))    {err = errInitMemory; goto clean2;}
      if (!(prnkd = new (nothrow) Double_t[p])) {err = errInitMemory; goto clean2;}
      if (!(mrnkd = new (nothrow) Double_t[m])) {err = errInitMemory; goto clean2;}

   // Fill arrays
      p = 0;
      m = 0;
      for (i=0; i<n; i++) {
         if      (msk[i] == 1) arrPM[p++] = x[i];
         else if (msk[i] == 0) arrMM[m++] = x[i];
      }//for_i

   // Sort array arrPM and fill pmtree
      TStat::Rank(p, arrPM, pmidx, prnki, kFALSE);
      // get rank rnkd analogoulsy to rank used for RMA in Bioconductor
      for (Int_t i=0; i<p; i++) pmord[i] = arrPM[pmidx[i]];
      TStat::Rank(p, pmord, prnkd);

      for (i=0; i<p; i++) {
         if (trim == 0.0) fMean1[i] += arrPM[pmidx[i]];
         else             pmsort     = arrPM[pmidx[i]];
         pmrank = prnkd[prnki[i]];
         pmtree->Fill();
      }//for_i

      pmtree->Write();

   // Sort array arrMM and fill mmtree
      TStat::Rank(m, arrMM, mmidx, mrnki, kFALSE);
      // get rank rnkd analogoulsy to rank used for RMA in Bioconductor
      for (Int_t i=0; i<m; i++) mmord[i] = arrMM[mmidx[i]];
      TStat::Rank(m, mmord, mrnkd);

      for (i=0; i<m; i++) {
         if (trim == 0.0) fMean2[i] += arrMM[mmidx[i]];
         else             mmsort     = arrMM[mmidx[i]];
         mmrank = mrnkd[mrnki[i]];
         mmtree->Fill();
      }//for_i

      mmtree->Write();

      fNData++;

   clean2:
      if (mrnkd)  {delete [] mrnkd; mrnkd = 0;}
      if (prnkd)  {delete [] prnkd; prnkd = 0;}
      if (mrnki)  {delete [] mrnki; mrnki = 0;}
      if (prnki)  {delete [] prnki; prnki = 0;}
      if (mmord)  {delete [] mmord; mmord = 0;}
      if (pmord)  {delete [] pmord; pmord = 0;}
      if (mmidx)  {delete [] mmidx; mmidx = 0;}
      if (pmidx)  {delete [] pmidx; pmidx = 0;}
      if (arrMM)  {delete [] arrMM; arrMM = 0;}
      if (arrPM)  {delete [] arrPM; arrPM = 0;}
      // need to delete trees from RAM (otherwise it is used by Calculate()!!)
      if (mmtree) {mmtree->Delete(""); mmtree = 0;}
      if (pmtree) {pmtree->Delete(""); pmtree = 0;}
   } else { 
      cerr << "Error: Option <" << fDataOpt.Data()
           << "> cannot be used with quantile normalizer." << endl;
      err = errGeneral;
   }//if

   savedir->cd();

   return err;
}//AddArray

//______________________________________________________________________________
Double_t *XQuantileNormalizer::GetArray(Int_t n, Double_t *x, Int_t *msk, 
                               const char *name)
{
   // Return array x of length n with name
   if(kCS) cout << "---XQuantileNormalizer::GetArray------" << endl;

// Parameter delta (default=1.0)
   Float_t delta = (fNPar > 1) ? fPars[1] : 1.0;

   TDirectory *savedir = gDirectory;
   fTmpFile->cd();

   if (strcmp(fDataOpt.Data(), "together") == 0) {
   // Get tree with name
      TString tname = TString(name) + ".rkt";
      Float_t rank  = 0;
      Int_t   frank = 0;
      TTree  *tree  = (TTree*)fTmpFile->Get(tname.Data()); 
      if (tree == 0) return 0;
      tree->SetBranchAddress("rankBr", &rank);

   // Check for presence of array fMean
      if (!fMean1) {
         cerr << "Error: No mean values for quantile normalization." << endl;
         return 0;
      }//if

   // Sort mean in order of original tree data 
      Int_t id = 0;
      for (Int_t i=0; i<n; i++) {
         if (msk[i] == 1) {
            tree->GetEntry(id++);
            frank = (Int_t)floor(rank);
            if (rank - frank > delta) {
               x[i] = 0.5*(fMean1[frank - 1] + fMean1[frank]);
            } else {
               x[i] = fMean1[frank - 1];
            }//if
         } else {
            x[i] = 0;  //not necessary?
         }//if
      }//for_i

   // Cleanup
      // need to delete tree from RAM
      tree->Delete(""); tree = 0;
   } else if (strcmp(fDataOpt.Data(), "separate") == 0) {
   // Get PM tree with name
      TString pmname = TString(name) + ".rkp";
      Float_t pmrank = 0;
      Int_t   pfrank = 0;
      TTree  *pmtree = (TTree*)fTmpFile->Get(pmname.Data()); 
      if (pmtree == 0) return 0;
      pmtree->SetBranchAddress("pmrankBr", &pmrank);

   // Get MM tree with name
      TString mmname = TString(name) + ".rkm";
      Float_t mmrank = 0;
      Int_t   mfrank = 0;
      TTree  *mmtree = (TTree*)fTmpFile->Get(mmname.Data()); 
      if (mmtree == 0) return 0;
      mmtree->SetBranchAddress("mmrankBr", &mmrank);

   // Check for presence of arrays fMean1 and fMean2
      if (!(fMean1 && fMean2)) {
         cerr << "Error: No mean values for quantile normalization." << endl;
         return 0;
      }//if

   // Sort fMean1 and fMean2 in order of original tree data 
      Int_t p = 0;
      Int_t m = 0;
      for (Int_t i=0; i<n; i++) {
         if (msk[i] == 1) {
            pmtree->GetEntry(p++);
            pfrank = (Int_t)floor(pmrank);
            if (pmrank - pfrank > delta) {
               x[i] = 0.5*(fMean1[pfrank - 1] + fMean1[pfrank]);
            } else {
               x[i] = fMean1[pfrank - 1];
            }//if
         } else if (msk[i] == 0) {
            mmtree->GetEntry(m++);
            mfrank = (Int_t)floor(mmrank);
            if (mmrank - mfrank > delta) {
               x[i] = 0.5*(fMean2[mfrank - 1] + fMean2[mfrank]);
            } else {
               x[i] = fMean2[mfrank - 1];
            }//if
         } else {
            x[i] = 0;  //not necessary?
         }//if
      }//for_i

   // Cleanup
      // need to delete trees from RAM
      mmtree->Delete(""); mmtree = 0;
      pmtree->Delete(""); pmtree = 0;
   } else {
      cerr << "Error: Option <" << fDataOpt.Data()
           << "> cannot be used with quantile normalizer." << endl;
      return 0;
   }//if

   savedir->cd();

   return x;
}//GetArray

//______________________________________________________________________________
Int_t XQuantileNormalizer::DoMean(Int_t n, Double_t *x)
{
   // Get mean array 
   // dramatically reduced memory consumption, can be used for many thousand trees
   if(kCS) cout << "------XQuantileNormalizer::DoMean------" << endl;

   if (strcmp(fDataOpt.Data(), "together") == 0) {
      for (Int_t i=0; i<fNMean1; i++) {
         fMean1[i] = fMean1[i]/fNData;
         x[i]      = fMean1[i];
      }//for_i
   } else if (strcmp(fDataOpt.Data(),"separate") == 0) {
      for (Int_t i=0; i<fNMean1; i++) {
         fMean1[i] = fMean1[i]/fNData;
//no         x[i]      = fMean1[i];
      }//for_i

      for (Int_t i=0; i<fNMean2; i++) {
         fMean2[i] = fMean2[i]/fNData;
//no         x[i]      = fMean2[i];
      }//for_i
   } else {
      cerr << "Error: Option <" << fDataOpt.Data()
           << "> cannot be used with quantile normalizer." << endl;
      return errGeneral;
   }//if

   return errNoErr;
}//DoMean

//______________________________________________________________________________
Int_t XQuantileNormalizer::DoTrimmedMean(Int_t n, Double_t *x, Double_t trim)
{
   // Calculate trimmed mean array  
   if(kCS) cout << "------XQuantileNormalizer::DoTrimmedMean------" << endl;

   if (strcmp(fDataOpt.Data(), "together") == 0) {
      TTree   **treek = new TTree*[fNData];
      Double_t *sortk = new Double_t[fNData];
      for (Int_t k=0; k<fNData; k++) sortk[k] = 0;

   // Init trees
      Int_t idx = 0;
      TKey *key = 0;
      TIter next(fTmpFile->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TTree *tmptree = (TTree*)fTmpFile->Get(key->GetName());
         treek[idx] = tmptree;
         treek[idx]->SetBranchAddress("sortBr", &sortk[idx]);
         idx++;
      }//while
      if (idx != fNData) {
         cerr << "Error: Number of trees for quantile normalization is not <"
              << fNData << ">." << endl;
         return errAbort;
      }//if

   // Create data array
      Double_t *arr = 0;
      if (!(arr = new (nothrow) Double_t[fNData])) return errInitMemory;

   // Create mean array
      fNMean1 = treek[0]->GetEntries();
      if (!fMean1) if (!(fMean1 = new (nothrow) Double_t[fNMean1])) return errInitMemory;

   // Calculate trimmed mean values of all tree entries
      for (Int_t i=0; i<fNMean1; i++) {
         // read entry i from trees and fill array
         for (Int_t k=0; k<fNData; k++) {
            treek[k]->GetEntry(i);
            arr[k] = sortk[k];
         }//for_k

         // calculate trimmed mean of array
         fMean1[i] = TStat::Mean(fNData, arr, trim);
         x[i]      = fMean1[i];
      }//for_i

   // Cleanup
      if (arr) {delete [] arr; arr = 0;}
      delete [] sortk;
      // need to delete tree from RAM
      for (Int_t k=0; k<fNData; k++) SafeDelete(treek[k]);
      delete [] treek;
   } else if (strcmp(fDataOpt.Data(),"separate") == 0) {
      TTree   **pmtreek = new TTree*[fNData];
      TTree   **mmtreek = new TTree*[fNData];
      Double_t *pmsortk = new Double_t[fNData];
      Double_t *mmsortk = new Double_t[fNData];
      for (Int_t k=0; k<fNData; k++) pmsortk[k] = mmsortk[k] = 0;

   // Init trees
      Int_t p   = 0;
      Int_t m   = 0;
      TKey *key = 0;
      TIter next(fTmpFile->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TTree *tmptree = (TTree*)fTmpFile->Get(key->GetName());
         if (tmptree->GetBranch("pmsortBr")) {
            pmtreek[p] = tmptree;
            pmtreek[p]->SetBranchAddress("pmsortBr", &pmsortk[p]);
            p++;
            if (p > fNData) return errGetTree;
         } else if (tmptree->GetBranch("mmsortBr")) {
            mmtreek[m] = tmptree;
            mmtreek[m]->SetBranchAddress("mmsortBr", &mmsortk[m]);
            m++;
            if (m > fNData) return errGetTree;
         }//if
      }//while

   // Create data array
      Double_t *arr = 0;
      if (!(arr = new (nothrow) Double_t[fNData])) return errInitMemory;

   // Create mean arrays
      fNMean1 = pmtreek[0]->GetEntries();
      fNMean2 = mmtreek[0]->GetEntries();
      if (!(fMean1 = new (nothrow) Double_t[fNMean1])) return errInitMemory;
      if (!(fMean2 = new (nothrow) Double_t[fNMean2])) return errInitMemory;

   // Calculate trimmed mean values of PM tree entries
      for (Int_t i=0; i<fNMean1; i++) {
         // read entry i from trees and fill array
         for (Int_t k=0; k<fNData; k++) {
            pmtreek[k]->GetEntry(i);
            arr[k] = pmsortk[k];
         }//for_k

         // calculate trimmed mean of array
         fMean1[i] = TStat::Mean(fNData, arr, trim);
//no         x[i]      = fMean1[i];
      }//for_i

   // Calculate trimmed mean values of MM tree entries
      for (Int_t i=0; i<fNMean2; i++) {
         // read entry i from trees and fill array
         for (Int_t k=0; k<fNData; k++) {
            mmtreek[k]->GetEntry(i);
            arr[k] = mmsortk[k];
         }//for_k

         // calculate trimmed mean of array
         fMean2[i] = TStat::Mean(fNData, arr, trim);
//no         x[i]      = fMean2[i];
      }//for_i

   // Cleanup
      if (arr) {delete [] arr; arr = 0;}
      delete [] mmsortk;
      delete [] pmsortk;
      // need to delete trees from RAM
      for (Int_t k=0; k<fNData; k++) {
         SafeDelete(mmtreek[k]);
         SafeDelete(pmtreek[k]);
      }//for_k
      delete [] mmtreek;
      delete [] pmtreek;
   } else {
      cerr << "Error: Option <" << fDataOpt.Data()
           << "> cannot be used with quantile normalizer." << endl;
      return errGeneral;
   }//if

   return errNoErr;
}//DoTrimmedMean

//______________________________________________________________________________
Int_t XQuantileNormalizer::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Normalize data 
   // Input is:
   //    x:   optional array of length n, not used
   //    msk: not used
   // Output is:
   //    y:   array containing mean values of quantile normalization
   if(kCS) cout << "------XQuantileNormalizer::Calculate------" << endl;

   Int_t err = errNoErr;

// Parameter trim (default=0)
   Double_t trim = (fNPar > 0) ? fPars[0] : 0.0;

   TDirectory *savedir = gDirectory;
   fTmpFile->cd();

// Calculate mean values
   if (trim == 0.0) err = this->DoMean(n, y);
   else             err = this->DoTrimmedMean(n, y, trim);

   savedir->cd();

   return err;
}//Calculate


