// File created: 08/05/2002                          last modified: 12/30/2010
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
* Aug 2008 - Add summarization algorithms FARMS and DFW, classes XFARMS, XDFW
* Oct 2008 - Add I/NI-call algorithm class XINICall
* Jan 2010 - Move methods from XHybridizer to XAlgorithm.
* Feb 2010 - Add exon splice/summarization algorithms, class XSpliceExpressor
* Dec 2010 - Add quality control algorithms, classes XQualifier and derived
*
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <new>        //needed for new (nothrow)
#include <algorithm>  //needed for sort (STL vector) for VC++ 

#include "TBranch.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TSystem.h"

//for Benchmark test only:
#include <TBenchmark.h>

#include "TStat.h"
#include "XPSHybridizer.h"

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XHybridizer);
ClassImp(XBackgrounder);
ClassImp(XCallDetector);
ClassImp(XExpressor);
ClassImp(XMultichipExpressor);
ClassImp(XSpliceExpressor);
ClassImp(XQualifier);
ClassImp(XSectorBackground);
ClassImp(XWeightedBackground);
ClassImp(XRMABackground);
ClassImp(XGCBackground);
ClassImp(XMeanDifferenceCall);
ClassImp(XDetectionCall);
ClassImp(XMAS4Call);
ClassImp(XDABGCall);
ClassImp(XINICall);
ClassImp(XArithmeticMean);
ClassImp(XGeometricMean);
ClassImp(XWeightedMean);
ClassImp(XGCCorrectedMean);
ClassImp(XWeightedDiff);
ClassImp(XAvgDif);
ClassImp(XTukeyBiweight);
ClassImp(XMedianPolish);
ClassImp(XFARMS);
ClassImp(XDFW);
ClassImp(XFIRMA);
ClassImp(XRMAQualifier);
ClassImp(XPLMQualifier);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybridizer                                                          //
//                                                                      //
// Base class for microarray hybridization algorithms                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XHybridizer::XHybridizer()
            :XAlgorithm()
{
   // Default Hybridizer constructor
   if(kCS) cout << "---XHybridizer::XHybridizer(default)------" << endl;

   fLength    = 0;
   fInten1    = 0;
   fStdev1    = 0;
   fNPix1     = 0;
   fInten2    = 0;
   fStdev2    = 0;
   fNPix2     = 0;
   fArray     = 0;
   fNDefPar   = 0;
   fMultichip = kFALSE;
}//Constructor

//______________________________________________________________________________
XHybridizer::XHybridizer(const char *name, const char *type)
            :XAlgorithm(name, type)
{
   // Normal Hybridizer constructor
   if(kCS) cout << "---XHybridizer::XHybridizer------" << endl;

   fLength    = 0;
   fInten1    = 0;
   fStdev1    = 0;
   fNPix1     = 0;
   fInten2    = 0;
   fStdev2    = 0;
   fNPix2     = 0;
   fArray     = 0;
   fNDefPar   = 0;
   fMultichip = kFALSE;
}//Constructor

//______________________________________________________________________________
XHybridizer::~XHybridizer()
{
   // Hybridizer destructor
   if(kCS) cout << "---XHybridizer::~XHybridizer------" << endl;

   this->DeleteArray();

   fInten1   = 0;
   fStdev1   = 0;
   fNPix1    = 0;
   fInten2   = 0;
   fStdev2   = 0;
   fNPix2    = 0;
}//Destructor

//______________________________________________________________________________
void XHybridizer::DeleteArray()
{
   // Delete array
   if(kCSa) cout << "------XHybridizer::DeleteArray------" << endl;

   if (fArray) {delete [] fArray; fArray = 0;}
}//DeleteArray

//______________________________________________________________________________
void XHybridizer::InitArrays(Int_t length, Double_t *inten1, Double_t *stdev1,
                  Int_t *npix1, Double_t *inten2, Double_t *stdev2, Int_t *npix2)
{
   // Initialize data arrays for hybridization algorithm
   if(kCSa) cout << "------XAlgorithm::InitArrays------" << endl;

   fLength = length;
   fInten1 = inten1;
   fStdev1 = stdev1;
   fNPix1  = npix1;
   fInten2 = inten2;
   fStdev2 = stdev2;
   fNPix2  = npix2;
}//InitArrays

//______________________________________________________________________________
Int_t XHybridizer::SetArray(Int_t length, Double_t *array)
{
   // Set size of array for calculation to length doubles and set contents
   if(kCSa) cout << "------XHybridizer::SetArray------" << endl;

   if (length == 0) return 1;
   if (array  == 0) return 1;

   if (!fArray || (fLength != length)) {
      this->DeleteArray();

      if (!(fArray = new (nothrow) Double_t[length])) return errInitMemory;
      fLength = length;
   }//if

   memcpy(fArray, array, length*sizeof(Double_t));

   return errNoErr;
}//SetArray


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBackgrounder                                                        //
//                                                                      //
// Base class for microarray background analysis algorithm              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBackgrounder::XBackgrounder()
              :XHybridizer()
{
   // Default Backgrounder constructor
   if(kCS) cout << "---XBackgrounder::XBackgrounder(default)------" << endl;

   fNRows = 0;
   fNCols = 0;
}//Constructor

//______________________________________________________________________________
XBackgrounder::XBackgrounder(const char *name, const char *type)
              :XHybridizer(name, type)
{
   // Normal Backgrounder constructor
   if(kCS) cout << "---XBackgrounder::XBackgrounder------" << endl;

   fNRows = 0;
   fNCols = 0;
}//Constructor

//______________________________________________________________________________
XBackgrounder::~XBackgrounder()
{
   // Backgrounder destructor
   if(kCS) cout << "---XBackgrounder::~XBackgrounder------" << endl;

   fNRows  = 0;
   fNCols  = 0;
   fOption = "subtractbg";  //default for backgrounder
}//Destructor

//______________________________________________________________________________
Double_t *XBackgrounder::AdjustIntensity(Int_t n, Double_t *inten, Double_t *bgrd,
                         Double_t *stdv)
{
   // Adjust intensity based on fOption
   if(kCS) cout << "------XBackgrounder::AdjustIntensity------" << endl;

   if (n     == 0 || inten == 0) return 0;
   if (inten != 0 && bgrd  == 0) return inten;

   if (strcmp(fOption.Data(), "subtractbg") == 0) {
      // (inten - bg)
      for (Int_t i=0; i<n; i++) {
         inten[i] = inten[i] - bgrd[i];
      }//for_i
   } else if (strcmp(fOption.Data(), "correctbg") == 0) {
      Double_t noisefrac = (fPars[fNPar-1] > 0) ? fPars[fNPar-1] : 0.5;

      // max(inten - bg, noisefrac*noise)
      for (Int_t i=0; i<n; i++) {
         inten[i] = TMath::Max(inten[i] - bgrd[i], noisefrac*stdv[i]);
      }//for_i
   } else if (strcmp(fOption.Data(), "attenuatebg") == 0) {
      Double_t h = -1.0;
      Double_t l = 0.005;
      if (fNPar > 1) {
         h = fPars[fNPar-1];
         l = fPars[fNPar-2];
      }//if

      // (X+sqrt(X^2+H))/2 where H = 4*Inten*Bgrd*L and X = Inten-Bgrd
      for (Int_t i=0; i<n; i++) {
         Double_t hh = (h < 0) ? 4*inten[i]*bgrd[i]*l : h;
         Double_t xx = inten[i] - bgrd[i];

         inten[i] = (xx + TMath::Sqrt(TMath::Power(xx, 2) + hh))/2.0;
      }//for_i
   }//if

   return inten;
}//AdjustIntensity

//______________________________________________________________________________
Double_t *XBackgrounder::AdjustError(Int_t n, Double_t *sd1, Double_t *sd2)
{
   // Adjust error
   if(kCSa) cout << "------XBackgrounder::AdjustError------" << endl;

   if (n   == 0 || sd1 == 0) return 0;
   if (sd1 != 0 && sd2 == 0) return sd1;

   for (Int_t i=0; i<n; i++) {
   // error propagation: sum or sum of squares???
      sd1[i] = TMath::Sqrt(sd1[i]*sd1[i] + sd2[i]*sd2[i]);
//??      sd1[i] = sd1[i] + sd2[i];
   }//for_i

   return sd1;
}//AdjustError


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCallDetector                                                        //
//                                                                      //
// Base class for microarray analysis present call algorithm            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XCallDetector::XCallDetector()
              :XHybridizer()
{
   // Default CallDetector constructor
   if(kCS) cout << "---XCallDetector::XCallDetector(default)------" << endl;

   fOption  = "";
   fDataOpt = "";
   fBgrdOpt = "";
}//Constructor

//______________________________________________________________________________
XCallDetector::XCallDetector(const char *name, const char *type)
              :XHybridizer(name, type)
{
   // Normal CallDetector constructor
   if(kCS) cout << "---XCallDetector::XCallDetector------" << endl;

   fOption  = "transcript"; //default for caller
   fDataOpt = "raw";        //default for caller
   fBgrdOpt = "none";
}//Constructor

//______________________________________________________________________________
XCallDetector::~XCallDetector()
{
   // CallDetector destructor
   if(kCS) cout << "---XCallDetector::~XCallDetector------" << endl;

}//Destructor

//______________________________________________________________________________
void XCallDetector::SetOptions(Option_t *opt)
{
   // Set options, e.g. "raw", "transcript:raw", "transcript:raw:correctbg"
   if(kCS) cout << "------XCallDetector::SetOptions------" << endl;

//	TString optcpy = opt;
//   char *options = (char*)optcpy.Data();

   char *options = new char[strlen(opt) + 1];
   char *doption = options;
   options = strcpy(options, opt);

   if (NumSeparators(opt, ":") == 0) {
      fOption  = "transcript"; 
      fDataOpt = strtok(options, ":");
   } else if (NumSeparators(opt, ":") == 1) {
      fOption  = strtok(options, ":");
      fDataOpt = strtok(NULL, ":");
   } else {
      fOption  = strtok(options, ":");
      fDataOpt = strtok(NULL, ":");
      fBgrdOpt = strtok(NULL, ":");
   }//if

   delete [] doption;
}//SetOptions


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExpressor                                                           //
//                                                                      //
// Base class for microarray expression analysis algorithm              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExpressor::XExpressor()
           :XHybridizer()
{
   // Default Expressor constructor
   if(kCS) cout << "---XExpressor::XExpressor(default)------" << endl;

   fOption    = "";
   fBgrdOpt   = "";
   fLogBase   = "0";
   fEstimator = 0;
   fSplicexpr = kFALSE;
}//Constructor

//______________________________________________________________________________
XExpressor::XExpressor(const char *name, const char *type)
           :XHybridizer(name, type)
{
   // Normal Expressor constructor
   if(kCS) cout << "---XExpressor::XExpressor------" << endl;

   fOption    = "transcript"; //default for expressor
   fBgrdOpt   = "none";       //default for expressor
   fLogBase   = "0";
   fEstimator = 0;
   fSplicexpr = kFALSE;
}//Constructor

//______________________________________________________________________________
XExpressor::~XExpressor()
{
   // Expressor destructor
   if(kCS) cout << "---XExpressor::~XExpressor------" << endl;

   SafeDelete(fEstimator);
}//Destructor

//______________________________________________________________________________
Int_t XExpressor::SetArray(Int_t length, Double_t *array)
{
   // Set size of array for calculation to length doubles and set contents
   // and optionally convert array to logarithm of base fLogBase
   if(kCSa) cout << "------XExpressor::SetArray------" << endl;

   Int_t err = XHybridizer::SetArray(length, array);

// Get neglog (default = 1.0)
   Double_t neglog = (fNPar > fNDefPar) ? fPars[fNDefPar] : 1.0;

   fArray = Array2Log(fLength, fArray, neglog, fLogBase);

   return err;
}//SetArray

//______________________________________________________________________________
void XExpressor::SetOptions(Option_t *opt)
{
   // Set options, e.g. "transcript" or "transcript:log2" or "transcript:correctbg:log2"
   if(kCS) cout << "------XExpressor::SetOptions------" << endl;

//	  TString optcpy = opt;
//   char *options = (char*)optcpy.Data();

   char *options = new char[strlen(opt) + 1];
   char *doption = options;
   options = strcpy(options,opt);

   if (NumSeparators(opt, ":") == 0) {
      fOption    = "transcript"; 
      fLogBase   = strtok(options, ":");
   } else if (NumSeparators(opt, ":") == 1) {
      fOption    = strtok(options, ":");
      fLogBase   = strtok(NULL, ":");
   } else if (NumSeparators(opt, ":") == 2) {
      fOption    = strtok(options, ":");
      fBgrdOpt   = strtok(NULL, ":");
      fLogBase   = strtok(NULL, ":");
   } else {
      fOption    = strtok(options, ":");
      fEstimator = TMEstimator::Estimator((const char*)strtok(NULL, ":"));
      fBgrdOpt   = strtok(NULL, ":");
      fLogBase   = strtok(NULL, ":");
   }//if

   delete [] doption;
}//SetOptions


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultichipExpressor                                                  //
//                                                                      //
// Base class for multichip expression analysis algorithm               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMultichipExpressor::XMultichipExpressor()
                    :XExpressor()
{
   // Default XMultichipExpressor constructor
   if(kCS) cout << "---XMultichipExpressor::XMultichipExpressor(default)------" << endl;

   fResiduals = 0;
   fMultichip = kTRUE;
}//Constructor

//______________________________________________________________________________
XMultichipExpressor::XMultichipExpressor(const char *name, const char *type)
                    :XExpressor(name, type)
{
   // Normal MultichipExpressor constructor
   if(kCS) cout << "---XMultichipExpressor::XMultichipExpressor------" << endl;

   fResiduals = 0;
   fMultichip = kTRUE;
}//Constructor

//______________________________________________________________________________
XMultichipExpressor::~XMultichipExpressor()
{
   // MultichipExpressor destructor
   if(kCS) cout << "---XMultichipExpressor::~XMultichipExpressor------" << endl;

   if (fResiduals) {delete [] fResiduals; fResiduals = 0;}
}//Destructor

//______________________________________________________________________________
Int_t XMultichipExpressor::DoMedianPolish(Int_t nrow, Int_t ncol, Double_t *inten, 
                           Double_t *level, Double_t *rowmed, Double_t *colmed,
                           Double_t *residu, Double_t *weight)
{
   // Calculate median polish
   if(kCSa) cout << "------XMultichipExpressor::DoMedianPolish------" << endl;

   Int_t err  = errNoErr;

// Get parameters or use default values
   Int_t    iter   = (fNPar > 0) ? (Int_t)fPars[0] : 10;
   Double_t eps    = (fNPar > 1) ? fPars[1] : 0.01;
   Int_t    medpol = (fNPar > 2) ? (Int_t)fPars[2] : 1;
   Double_t totmed = 0;

   if (iter <= 0 || iter >= 100) {
      cout << "Warning: Number of iterations <" << iter
           << "> is not in range [1,99], setting iter to default, iter = 10."
           << endl;
      iter = 10;
   }//if

// Median polish
   if (medpol == 1) {
      totmed = TStat::MedianPolish(nrow, ncol, inten, rowmed, colmed, residu, iter, eps);
   } else if (medpol == 2) {
      totmed = TStat::MedianPolishTranspose(nrow, ncol, inten, rowmed, colmed, residu, iter, eps);
   } else {
      cout << "Warning: <medpol = " << medpol
           << "> is not valid, setting medpol to default <medpol = 1>."
           << endl;
   }//if

// Convert expression levels
   for (Int_t j=0; j<ncol; j++) level[j] = totmed + colmed[j];
   level = Array2Pow(ncol, level, fLogBase);

// Pseudo error and weight
   if (fEstimator) {
      if (colmed) {
         colmed = this->PseudoError(nrow, ncol, residu, colmed);
         if (colmed == 0) return errInitMemory;
      }//if

      if (weight) {
         weight = this->PseudoWeight(nrow, ncol, residu, weight);
         if (weight == 0) return errInitMemory;
      }//if
   }//if

// Convert standard deviation/error
   colmed = Array2Pow(ncol, colmed, fLogBase);
 
   return err;
}//DoMedianPolish

//______________________________________________________________________________
Double_t *XMultichipExpressor::PseudoError(Int_t nrow, Int_t ncol,  
                               Double_t *residu, Double_t *se)
{
   // Calculate pseudo standard error using a specified M-estimator
   // adapted from Bioconductor package "affyPLM" function "compute_pseudoSE()"
   // created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
   if(kCSa) cout << "------XMultichipExpressor::PseudoError------" << endl;

   Int_t n = nrow*ncol;

   Double_t *arr = 0;
   if (!(arr = new (nothrow) Double_t[n])) return 0;

   for (Int_t i=0; i<n; i++) arr[i] = TMath::Abs(residu[i]);

   // median of absolute residuals (why 1/0.6745?)
   Double_t md = TStat::Median(n, arr)/0.6745;

   // protect against residuals = 0
   if (md == 0) {
      for (Int_t j=0; j<ncol; j++) se[j] = 0.0;
      delete [] arr;
      return se;
   }//if

   Double_t residSE = 0.0;
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t j=0; j<ncol; j++) {
         Double_t r = residu[j + i*ncol];
         residSE += fEstimator->Weight(r/md)*r*r;
      }//for_j
   }//for_i
   residSE = TMath::Sqrt(residSE/(Double_t)(n - (nrow + ncol - 1)));

   for (Int_t j=0; j<ncol; j++) {
      Double_t weight = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         weight += fEstimator->Weight(residu[j + i*ncol]/md);
      }//for_i
      se[j] = residSE/TMath::Sqrt(weight);
   }//for_j

   delete [] arr;

   return se;
}//PseudoError

//______________________________________________________________________________
Double_t *XMultichipExpressor::PseudoWeight(Int_t nrow, Int_t ncol,  
                               Double_t *residu, Double_t *weight)
{
   // Calculate pseudo weights using a specified M-estimator
   // adapted from Bioconductor package "affyPLM" function "compute_pseudoweights()"
   // created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
   if(kCSa) cout << "------XMultichipExpressor::PseudoWeight------" << endl;

   Int_t n = nrow*ncol;

   Double_t *arr = 0;
   if (!(arr = new (nothrow) Double_t[n])) return 0;

   for (Int_t i=0; i<n; i++) arr[i] = TMath::Abs(residu[i]);

   // median of absolute residuals (why 1/0.6745?)
   Double_t md = TStat::Median(n, arr)/0.6745;

   // protect against residuals = 0
   if (md == 0) {
      for (Int_t j=0; j<ncol; j++) {
         for (Int_t i=0; i<nrow; i++) {
            weight[j + i*ncol] = 0.0;
         }//for_i
      }//for_j

      delete [] arr;

      return weight;
   }//if

   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         weight[j + i*ncol] = fEstimator->Weight(residu[j + i*ncol]/md);
      }//for_i
   }//for_j

   delete [] arr;

   return weight;
}//PseudoWeight


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSpliceExpressor                                                     //
//                                                                      //
// Algorithms used for summarization and splice detection               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSpliceExpressor::XSpliceExpressor()
                 :XMultichipExpressor()
{
   // Default SpliceExpressor constructor
   if(kCS) cout << "---XSpliceExpressor::XSpliceExpressor(default)------" << endl;

   fSplicexpr = kTRUE;
}//Constructor

//______________________________________________________________________________
XSpliceExpressor::XSpliceExpressor(const char *name, const char *type)
                 :XMultichipExpressor(name, type)
{
   // Normal SpliceExpressor constructor
   if(kCS) cout << "---XSpliceExpressor::XSpliceExpressor------" << endl;

   fSplicexpr = kTRUE;
}//Constructor

//______________________________________________________________________________
XSpliceExpressor::~XSpliceExpressor()
{
   // SpliceExpressor destructor
   if(kCS) cout << "---XSpliceExpressor::~XSpliceExpressor------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XQualifier                                                           //
//                                                                      //
// Base class for microarray quality control algorithms                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XQualifier::XQualifier()
           :XAlgorithm()
{
   // Default Qualifier constructor
   if(kCS) cout << "---XQualifier::XQualifier(default)------" << endl;

   fCoefHi    = 1.2;
   fCoefLo    = 0.8;
   fQualOpt   = "";
   fExpressor = 0;
}//Constructor

//______________________________________________________________________________
XQualifier::XQualifier(const char *name, const char *type)
           :XAlgorithm(name, type)
{
   // Normal Qualifier constructor
   if(kCS) cout << "---XQualifier::XQualifier------" << endl;

   fCoefHi    = 1.2;
   fCoefLo    = 0.8;
   fQualOpt   = "raw";
   fExpressor = new XMultichipExpressor(name, type);
}//Constructor

//______________________________________________________________________________
XQualifier::~XQualifier()
{
   // Qualifier destructor
   if(kCS) cout << "---XQualifier::~XQualifier------" << endl;

   SafeDelete(fExpressor);
}//Destructor

//______________________________________________________________________________
Int_t XQualifier::Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                  Double_t *w, Int_t *msk)
{
   // Calculate median polish
   if(kCSa) cout << "------XQualifier::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t nrow = (Int_t)(fExpressor->GetLength() / n);
   Int_t ncol = n;

//   if (fResiduals) {delete [] fResiduals; fResiduals = 0;}
//   if (!(fResiduals = new (nothrow) Double_t[fLength])) return errInitMemory;

   Double_t *rowmed = 0;
   if (!(rowmed = new (nothrow) Double_t[nrow])) return errInitMemory;

   err = fExpressor->DoMedianPolish(nrow, ncol, fExpressor->GetArray(),
                                    x, rowmed, y, z, w);

   if (rowmed) {delete [] rowmed; rowmed = 0;}
 
   return err;
}//Calculate

//______________________________________________________________________________
Double_t XQualifier::MeanBorder(Int_t begin, Int_t end, Double_t *arr)
{
   // Compute mean for border arr
   if(kCS) cout << "------XQualifier::MeanBorder------" << endl;

   if (begin == end) return arr[begin];

   Double_t mean = 0.0;
   for (Int_t i=begin; i<end; i++) mean += arr[i];

   return mean/(end - begin);
}//MeanBorder

//______________________________________________________________________________
void XQualifier::HiLoBorder(Int_t begin, Int_t end, Double_t *arr, Short_t *msk,
                            Double_t mean, Double_t &meanhi, Double_t &meanlo)
{
   // Compute high and low mean for border arr
   if(kCS) cout << "------XQualifier::HiLoBorder------" << endl;

   Int_t idxhi = 0;
   Int_t idxlo = 0;
   for (Int_t i=begin; i<end; i++) {
      if (arr[i] > fCoefHi * mean) {msk[i] =  1; meanhi += arr[i]; idxhi++;} else
      if (arr[i] < fCoefLo * mean) {msk[i] = -1; meanlo += arr[i]; idxlo++;}
   }//for_i

   if (idxhi > 0) meanhi = meanhi/idxhi;
   if (idxlo > 0) meanlo = meanlo/idxlo;
}//HiLoBorder

//______________________________________________________________________________
Int_t XQualifier::SetArray(Int_t length, Double_t *array)
{
   // Set size of array for calculation to length doubles and set contents
   // and optionally convert array to logarithm of base fLogBase
   if(kCSa) cout << "------XQualifier::SetArray------" << endl;

   return fExpressor->SetArray(length, array);
}//SetArray

//______________________________________________________________________________
void XQualifier::SetOptions(Option_t *opt)
{
   // Set fDataOpt then pass other options to fExpressor
   if(kCS) cout << "------XQualifier::SetOptions------" << endl;

   fQualOpt = SubString(opt, ":", 0);

   char *options = new char[strlen(opt) + 1];
   char *doption = options;
   options = (char*)(strchr(opt,':') + 1);

   fExpressor->SetOptions(options);

   delete [] doption;
}//SetOptions


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSectorBackground                                                    //
//                                                                      //
// Class for sector background algorithm                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSectorBackground::XSectorBackground()
                  :XBackgrounder()
{
   // Default SectorBackground constructor
   if(kCS) cout << "---XSectorBackground::XSectorBackground(default)------" << endl;

   fNDefPar = 4;
}//Constructor

//______________________________________________________________________________
XSectorBackground::XSectorBackground(const char *name, const char *type)
                  :XBackgrounder(name, type)
{
   // Normal SectorBackground constructor
   if(kCS) cout << "---XSectorBackground::XSectorBackground------" << endl;

   fNDefPar = 4;
}//Constructor

//______________________________________________________________________________
XSectorBackground::~XSectorBackground()
{
   // SectorBackground destructor
   if(kCS) cout << "---XSectorBackground::~XSectorBackground------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XSectorBackground::Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                         Int_t *msk)
{
   // Calculate sector background
   if(kCS) cout << "------XSectorBackground::Calculate------" << endl;

   Int_t err = errNoErr;

   Int_t i, j, ix, jy, ij, xy;
   Int_t num2, maxX, maxY, max;
   Int_t startX, endX, startY, endY;

// Test number of parameters
   if (TestNumParameters(4) != errNoErr) return errInitParameters;

// Set parameter variables
   Double_t percent = fPars[0];
   Int_t numRowSec = (Int_t)fPars[1]; 
   Int_t numColSec = (Int_t)fPars[2];
   Int_t numIter   = (Int_t)fPars[3];
   Int_t numSec    = numRowSec * numColSec;  // number of sectors
   Int_t numCells  = 0;              // number of cells in sector

   Int_t sizeX = (Int_t)TMath::Floor((Double_t)fNRows/numRowSec);
   Int_t sizeY = (Int_t)TMath::Floor((Double_t)fNCols/numColSec);
   Int_t modX  = fNRows % numRowSec;
   Int_t modY  = fNCols % numColSec;

// Init local arrays
   Int_t *lenX = 0; 
   Int_t *lenY = 0; 
   Double_t *sector    = 0;
   Double_t *secSorted = 0;
   Double_t *arrMean   = 0;
   Double_t *arrStdv   = 0;
   Double_t *arrSmooth = 0;


// Calculate initial size lenX x lenY of sectors
   if (!(lenX = new (nothrow) Int_t[numRowSec])) {err = errInitMemory; goto cleanup;}
   if (!(lenY = new (nothrow) Int_t[numColSec])) {err = errInitMemory; goto cleanup;} 

   for (i=0; i<numRowSec; i++) lenX[i] = sizeX;
   for (i=0; i<numColSec; i++) lenY[i] = sizeY;

// Adjust sector sizes if modulo > 0
   num2 = (Int_t)TMath::Floor(numRowSec/2.0);
   if (modX > 0) {
      for (i=0; i<num2; i++) {
         lenX[i] = sizeX + 1;
         modX--;
         if (modX == 0) break;
         lenX[numRowSec - 1 - i] = sizeX + 1;
         modX--;
         if (modX == 0) break;
      }//for_i
   }//if

   num2 = (Int_t)TMath::Floor(numColSec/2.0);
   if (modY > 0) {
      for (i=0; i<num2; i++) {
         lenY[i] = sizeY + 1;
         modY--;
         if (modY == 0) break;
         lenY[numColSec - 1 - i] = sizeY + 1;
         modY--;
         if (modY == 0) break;
      }//for_i
   }//if

// Get size max of largest sector
   maxX = maxY = 0;
   for (i=0; i<numRowSec; i++) {
      maxX = (maxX >= lenX[i] ? maxX : lenX[i]);
   }//for_i
   for (i=0; i<numColSec; i++) {
      maxY = (maxY >= lenY[i] ? maxY : lenY[i]);
   }//for_i
   max = maxX * maxY; // size of largest sector

// Create arrays
   if (!(sector    = new (nothrow) Double_t[max]))    {err = errInitMemory; goto cleanup;} 
   if (!(secSorted = new (nothrow) Double_t[max]))    {err = errInitMemory; goto cleanup;} 
   if (!(arrMean   = new (nothrow) Double_t[numSec])) {err = errInitMemory; goto cleanup;} 
   if (!(arrStdv   = new (nothrow) Double_t[numSec])) {err = errInitMemory; goto cleanup;} 
   if (!(arrSmooth = new (nothrow) Double_t[numSec])) {err = errInitMemory; goto cleanup;} 

// Calculate mean value for each sector
   startX = 0;
   endY   = 0;
   for (i=0; i<numRowSec; i++) {
      startY = 0;
      endX   = startX + (lenX[i] - 1);

      for (j=0; j<numColSec; j++) {
         endY = startY + (lenY[j] - 1);

         numCells = 0;
         for (ix=startX; ix<=endX; ix++) {
            for (jy=startY; jy<=endY; jy++) {
            // sector array of all or only MM intensities
               xy = XY2Index(ix, jy, fNCols);
               if (msk[xy] == 1) sector[numCells++] = x[xy];
            }//for_jy
         }//for_ix

         // sort 
         Int_t *index = 0;
         if (!(index = new (nothrow) Int_t[numCells])) {err = errInitMemory; goto cleanup;}  
         TMath::Sort(numCells, sector, index, kFALSE);
         for (Int_t m=0; m<numCells; m++) {
            secSorted[m] = sector[index[m]];
         }//for_m
         delete [] index;

         // mean and stdev of percent lowest data
         Double_t mean = 0.0;
         Double_t var  = 0.0;
         Int_t numCells4Bg = (Int_t)(numCells * percent);
         if (numCells4Bg > 0) {
            mean = TStat::Mean(numCells4Bg, secSorted);
            var  = TStat::Var(numCells4Bg, secSorted, mean);
         }//if

         ij = XY2Index(j, i, numRowSec);
//or??         ij = XY2Index(i, j, numColSec);
         arrMean[ij] = mean;
         arrStdv[ij] = TMath::Sqrt(var);
         
         startY = endY + 1;
      }//for_j

      startX = endX + 1;
   }//for_i

// Smooth background sectors
   for (Int_t iter=0; iter<numIter; iter++) {
      this->Smooth(arrMean, arrSmooth, numRowSec, numColSec);
      for (i=0; i<numRowSec; i++) {
         for (j=0; j<numColSec; j++) {
            ij = XY2Index(j, i, numRowSec);
//or??            ij = XY2Index(i, j, numColSec);
            arrMean[ij] = arrSmooth[ij];
         }//for_j
      }//for_i
   }//for_iter

// Return background array y and return stdev in array x
   startX = 0;
   endY   = 0;
   for (i=0; i<numRowSec; i++) {
      startY = 0;
      endX   = startX + (lenX[i] - 1);

      for (j=0; j<numColSec; j++) {
         endY = startY + (lenY[j] - 1);

         numCells = 0;
         for (ix=startX; ix<=endX; ix++) {
            for (jy=startY; jy<=endY; jy++) {
               xy = XY2Index(ix, jy, fNCols);
               ij = XY2Index(j, i, numRowSec);
//or??               ij = XY2Index(i, j, numColSec);
               y[xy] = arrMean[ij];
               z[xy] = arrStdv[ij];
            }//for_jy
         }//for_ix

         startY = endY + 1;
      }//for_j

      startX = endX + 1;
   }//for_i

// Return background-corrected intensity in array x
   x = AdjustIntensity(n, x, y, z);

cleanup:
// delete locally created variables
   if (arrSmooth) {delete [] arrSmooth; arrSmooth = 0;}
   if (arrStdv)   {delete [] arrStdv;   arrStdv   = 0;}
   if (arrMean)   {delete [] arrMean;   arrMean   = 0;}
   if (secSorted) {delete [] secSorted; secSorted = 0;}
   if (sector)    {delete [] sector;    sector    = 0;}
   if (lenY)      {delete [] lenY;      lenY      = 0;}
   if (lenX)      {delete [] lenX;      lenX      = 0;}

   return err;
}//Calculate

//______________________________________________________________________________
void XSectorBackground::Smooth(const Double_t *arrIn, Double_t *arrOut, 
                        Int_t numrows, Int_t numcols)
{
   // Smooth sector background
   if(kCS) cout << "------XSectorBackground::Smooth------" << endl;

// Smooth background data: mean(i-1 to i+1, j-1 to j+1)
   Int_t startR, endR, startC, endC;
   Double_t sum;
   for (Int_t i=0; i<numrows; i++) {
      for(Int_t j=0; j<numcols; j++){
         startR = (i == 0) ? 0 : i-1;
         endR   = (i == numrows-1) ? numrows-1 : i+1;
         startC = (j == 0) ? 0 : j-1;
         endC   = (j == numcols-1) ? numcols-1 : j+1;
         sum    = 0;

         for (Int_t r=startR; r<=endR; r++) {
            for(Int_t c=startC; c<=endC; c++){
               sum += arrIn[r*numrows+c];
            }//for_c
         }//for_r

         arrOut[i*numrows+j] = sum / ((endR - startR + 1)*(endC - startC + 1));
      }//for_j
   }//for_i
}//Smooth


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XWeightedBackground                                                  //
//                                                                      //
// Class for weighted sector background algorithm from MAS5             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XWeightedBackground::XWeightedBackground()
                    :XSectorBackground()
{
   // Default WeightedBackground constructor
   if(kCS) cout << "---XWeightedBackground::XWeightedBackground(default)------" << endl;

   fNDefPar = 5;
}//Constructor

//______________________________________________________________________________
XWeightedBackground::XWeightedBackground(const char *name, const char *type)
                    :XSectorBackground(name, type)
{
   // Normal WeightedBackground constructor
   if(kCS) cout << "---XWeightedBackground::XWeightedBackground------" << endl;

   fNDefPar = 5;
}//Constructor

//______________________________________________________________________________
XWeightedBackground::~XWeightedBackground()
{
   // WeightedBackground destructor
   if(kCS) cout << "---XWeightedBackground::~XWeightedBackground------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XWeightedBackground::Calculate(Int_t n, Double_t *x, Double_t *y,
                           Double_t *z, Int_t *msk)
{
   // Calculate weighted sector background as described in MAS5
   if(kCS) cout << "------XWeightedBackground::Calculate------" << endl;

   Int_t err = errNoErr;

   Int_t i, j, k, l, ix, jy, ij, kl;
   Int_t num2, maxX, maxY, max;
   Int_t startX, endX, startY, endY;
   Double_t distX, distY;

// Test number of parameters
   if (TestNumParameters(5) != errNoErr) return errInitParameters;

// Set  parameter variables
   Double_t percent = fPars[0];
   Int_t numRowSec  = (Int_t)fPars[1]; 
   Int_t numColSec  = (Int_t)fPars[2];
   Int_t numIter    = (Int_t)fPars[3];
   Int_t smooth     = (Int_t)fPars[4];
   Int_t numSec     = numRowSec * numColSec;  // number of sectors
   Int_t numCells   = 0;              // number of cells in sector

   Int_t sizeX = (Int_t)TMath::Floor((Double_t)fNRows/numRowSec);
   Int_t sizeY = (Int_t)TMath::Floor((Double_t)fNCols/numColSec);
   Int_t modX  = fNRows % numRowSec;
   Int_t modY  = fNCols % numColSec;

// Init local arrays
   Int_t *lenX = 0; 
   Int_t *lenY = 0; 
   Int_t *cenX = 0; 
   Int_t *cenY = 0; 
   Double_t *sector    = 0;
   Double_t *secSorted = 0;
   Double_t *arrMean   = 0;
   Double_t *arrStdv   = 0;
   Double_t *arrSmooth = 0;
   Float_t **weight    = 0; //use Float_t to save memory

// Calculate initial size lenX x lenY of sectors
   if (!(lenX = new (nothrow) Int_t[numRowSec])) {err = errInitMemory; goto cleanup;}
   if (!(lenY = new (nothrow) Int_t[numColSec])) {err = errInitMemory; goto cleanup;} 

   for (i=0; i<numRowSec; i++) lenX[i] = sizeX;
   for (i=0; i<numColSec; i++) lenY[i] = sizeY;

// Adjust sector sizes if modulo > 0
   num2 = (Int_t)TMath::Floor(numRowSec/2.0);
   if (modX > 0) {
      for (i=0; i<num2; i++) {
         lenX[i] = sizeX + 1;
         modX--;
         if (modX == 0) break;
         lenX[numRowSec - 1 - i] = sizeX + 1;
         modX--;
         if (modX == 0) break;
      }//for_i
   }//if

   num2 = (Int_t)TMath::Floor(numColSec/2.0);
   if (modY > 0) {
      for (i=0; i<num2; i++) {
         lenY[i] = sizeY + 1;
         modY--;
         if (modY == 0) break;
         lenY[numColSec - 1 - i] = sizeY + 1;
         modY--;
         if (modY == 0) break;
      }//for_i
   }//if

// Get centroids for each sector
   if (!(cenX = new (nothrow) Int_t[numRowSec])) {err = errInitMemory; goto cleanup;}
   if (!(cenY = new (nothrow) Int_t[numColSec])) {err = errInitMemory; goto cleanup;} 

   cenX[0] = (Int_t)TMath::Ceil(lenX[0]/2.0);
   cenY[0] = (Int_t)TMath::Ceil(lenY[0]/2.0);
   for (i=1; i<numRowSec; i++) {
      cenX[i] = cenX[i-1] + lenX[i];
   }//for_i
   for (i=1; i<numColSec; i++) {
      cenY[i] = cenY[i-1] + lenY[i];
   }//for_i

// Get size max of largest sector (to init arrays)
   maxX = maxY = 0;
   for (i=0; i<numRowSec; i++) {
      maxX = (maxX >= lenX[i] ? maxX : lenX[i]);
   }//for_i
   for (i=0; i<numColSec; i++) {
      maxY = (maxY >= lenY[i] ? maxY : lenY[i]);
   }//for_i
   max = maxX * maxY; // size of largest sector

// Create arrays
   if (!(sector    = new (nothrow) Double_t[max]))    {err = errInitMemory; goto cleanup;} 
   if (!(secSorted = new (nothrow) Double_t[max]))    {err = errInitMemory; goto cleanup;} 
   if (!(arrMean   = new (nothrow) Double_t[numSec])) {err = errInitMemory; goto cleanup;} 
   if (!(arrStdv   = new (nothrow) Double_t[numSec])) {err = errInitMemory; goto cleanup;} 
   if (!(arrSmooth = new (nothrow) Double_t[numSec])) {err = errInitMemory; goto cleanup;} 

// Create weight table to store weights (use Float_t to save memory)
//PROBLEM: weight: exon array: 6.5 Millions x 16 x 8 byte >= 800 MB RAM!!!!!
//   if (!(weight    = new Double_t[n*numSec])) {err = errInitMemory; goto cleanup;}  
   if (!(weight = new (nothrow) Float_t*[numSec])) {err = errInitMemory; goto cleanup;} 
   for (Int_t i=0; i<numSec; i++) {
      weight[i] = 0;
      if (!(weight[i] = new (nothrow) Float_t[n])) {err = errInitMemory; goto cleanup;} 
   }//for_i

// Calculate mean value for each sector
   startX = 0;
   for (i=0; i<numRowSec; i++) {
      endX = startX + (lenX[i] - 1);

      startY = 0;
      for (j=0; j<numColSec; j++) {
         endY = startY + (lenY[j] - 1);

         // get intensities and weights for current sector
         numCells = 0;
         for (ix=startX; ix<=endX; ix++) {
            for (jy=startY; jy<=endY; jy++) {
               ij = XY2Index(ix, jy, fNCols);

               // weight as reciprocal distance to each centroid
               for (k=0; k<numRowSec; k++) {
                  distX = ix - cenX[k];
                  for (l=0; l<numColSec; l++) {
                     distY = jy - cenY[l];
                     kl = XY2Index(l, k, numRowSec);
//or?                     kl = XY2Index(k, l, numColSec);

                     weight[kl][ij] = 1.0/(distX*distX + distY*distY + smooth);
                  }//for_l
               }//for_k

               // sector array of both, PM only, or only MM intensities
               if (msk[ij] == 1) sector[numCells++] = x[ij];
            }//for_jy
         }//for_ix

         // sort intensities
         Int_t *index = 0;
         if (!(index = new (nothrow) Int_t[numCells])) {err = errInitMemory; goto cleanup;}
         TMath::Sort(numCells, sector, index, kFALSE);
         for (Int_t m=0; m<numCells; m++) {
            secSorted[m] = sector[index[m]];
         }//for_m
         delete [] index;

         // mean and stdev of percent lowest data
         Double_t mean = 0.0;
         Double_t var  = 0.0;
         Int_t numCells4Bg = (Int_t)(numCells * percent);
         if (numCells4Bg > 0) {
            mean = TStat::Mean(numCells4Bg, secSorted);
            var  = TStat::Var(numCells4Bg, secSorted, mean);
         }//if

         ij = XY2Index(j, i, numRowSec);
//or?         ij = XY2Index(i, j, numColSec);
         arrMean[ij] = mean;
         arrStdv[ij] = TMath::Sqrt(var);
         
         startY = endY + 1;
      }//for_j

      startX = endX + 1;
   }//for_i

// Smooth background sectors
   for (Int_t iter=0; iter<numIter; iter++) {
      this->Smooth(arrMean, arrSmooth, numRowSec, numColSec);
      for (i=0; i<numRowSec; i++) {
         for (j=0; j<numColSec; j++) {
            ij = XY2Index(j, i, numRowSec);
//or?            ij = XY2Index(i, j, numColSec);
            arrMean[ij] = arrSmooth[ij];
         }//for_j
      }//for_i
   }//for_iter

// Return background array y and return stdev in array z
   startX = 0;
   for (i=0; i<numRowSec; i++) {
      endX = startX + (lenX[i] - 1);

      startY = 0;
      for (j=0; j<numColSec; j++) {
         endY = startY + (lenY[j] - 1);

         for (ix=startX; ix<=endX; ix++) {
            for (jy=startY; jy<=endY; jy++) {
               ij = XY2Index(ix, jy, fNCols);

               // background and noise corrected by weights
               Double_t summ = 0;
               Double_t sums = 0;
               Double_t sumw = 0;
               for (k=0; k<numRowSec; k++) {
                  for (l=0; l<numColSec; l++) {
                     kl = XY2Index(l, k, numRowSec);
//or?                     kl = XY2Index(k, l, numColSec);

                     summ += weight[kl][ij]*arrMean[kl];
                     sums += weight[kl][ij]*arrStdv[kl];
                     sumw += weight[kl][ij];
                  }//for_l
               }//for_k

               y[ij] = summ/sumw;
               z[ij] = sums/sumw;
            }//for_jy
         }//for_ix

         startY = endY + 1;
      }//for_j

      startX = endX + 1;
   }//for_i

// Return background-corrected intensity in array x
   x = AdjustIntensity(n, x, y, z);

cleanup:
   // delete weight table
   for (Int_t i=0; i<numSec; i++) {
      if (weight[i]) {delete [] weight[i]; weight[i] = 0;}
   }//for_i
   if (weight) delete [] weight;

// delete locally created variables
   if (arrSmooth) {delete [] arrSmooth; arrSmooth = 0;}
   if (arrStdv)   {delete [] arrStdv;   arrStdv   = 0;}
   if (arrMean)   {delete [] arrMean;   arrMean   = 0;}
   if (secSorted) {delete [] secSorted; secSorted = 0;}
   if (sector)    {delete [] sector;    sector    = 0;}
   if (cenY)      {delete [] cenY;      cenY      = 0;}
   if (cenX)      {delete [] cenX;      cenX      = 0;}
   if (lenY)      {delete [] lenY;      lenY      = 0;}
   if (lenX)      {delete [] lenX;      lenX      = 0;}

   return err;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRMABackground                                                       //
//                                                                      //
// Class for global RMA background algorithm                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRMABackground::XRMABackground()
               :XBackgrounder()
{
   // Default RMABackground constructor
   if(kCS) cout << "---XRMABackground::XRMABackground(default)------" << endl;

   fNDefPar = 0;
   fKernel  = "";
}//Constructor

//______________________________________________________________________________
XRMABackground::XRMABackground(const char *name, const char *type)
               :XBackgrounder(name, type)
{
   // Normal RMABackground constructor
   if(kCS) cout << "---XRMABackground::XRMABackground------" << endl;

   fNDefPar = 0;
   fKernel  = "epanechnikov";
}//Constructor

//______________________________________________________________________________
XRMABackground::~XRMABackground()
{
   // RMABackground destructor
   if(kCS) cout << "---XRMABackground::~XRMABackground------" << endl;

}//Destructor

//______________________________________________________________________________
void XRMABackground::SetOptions(Option_t *opt)
{
   // Set options, e.g. "myoption" or "myoption:kernel"
   if(kCS) cout << "------XRMABackground::SetOptions------" << endl;

//   TString optcpy = opt;
//   char *options = (char*)optcpy.Data();

   char *options = new char[strlen(opt) + 1];
   char *doption = options;
   options = strcpy(options, opt);

   if (NumSeparators(opt, ":") == 0) {
      fOption = strtok(options, ":");
   } else {
      fOption = strtok(options, ":");
      fKernel = strtok(NULL, ":");
   }//if

   delete [] doption;
}//SetOptions

//______________________________________________________________________________
Int_t XRMABackground::Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                      Int_t *msk)
{
   // Calculate RMA background using a global model for the distribution of
   // probe intensities
   if(kCS) cout << "------XRMABackground::Calculate------" << endl;

   Int_t err = errNoErr;

   Double_t pars[3];

// Init local arrays
   Int_t     p     = 0;
   Int_t     m     = 0;
   Double_t *arrPM = 0;
   Double_t *arrMM = 0;
   Double_t *weiPM = 0; //weight for PM, currently 1.0
   Double_t *weiMM = 0; //weigth for MM, currently 1.0

// Get size of arrays for PM and MM
   for (Int_t i=0; i<n; i++) {
      if      (msk[i] == 1) p++;
      else if (msk[i] == 0) m++;
   }//for_i

// Create local arrays
   if (!(arrPM = new (nothrow) Double_t[p])) {err = errInitMemory; goto cleanup;}
   if (!(weiPM = new (nothrow) Double_t[p])) {err = errInitMemory; goto cleanup;}
   if (!(arrMM = new (nothrow) Double_t[m])) {err = errInitMemory; goto cleanup;}
   if (!(weiMM = new (nothrow) Double_t[m])) {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<p; i++) arrPM[i] = weiPM[i] = 0.0;
   for (Int_t i=0; i<m; i++) arrMM[i] = weiMM[i] = 0.0;

// Fill PM and MM arrays and set weights to one (weights are ignored)
   p = 0;
   m = 0;
   for (Int_t i=0; i<n; i++) {
      if (msk[i] == 1) {
         if (x[i] == 0) continue;  // protect against broken chips
         arrPM[p] = x[i];
         weiPM[p] = 1.0;
         p++;
      } else if (msk[i] == 0) {
         if (x[i] == 0) continue;
         arrMM[m] = x[i];
         weiMM[m] = 1.0;
         m++;
      }//if
   }//for_i

   if ((strcmp(fOption.Data(), "pmonly") == 0) && (p != 0)) {
      err = this->ComputeParameters(p, arrPM, weiPM, pars);
      if (err == errNoErr) this->Adjust(p, arrPM, pars);
   } else if ((strcmp(fOption.Data(), "mmonly") == 0) && (m != 0)) {
      err = this->ComputeParameters(m, arrMM, weiMM, pars);
      if (err == errNoErr) this->Adjust(m, arrMM, pars);
   } else if ((strcmp(fOption.Data(), "both") == 0)  && (p != 0) && (m != 0)) {
      if (p != m) {
         cout << "Warning: Number of PMs <" << p
              << "> is not equal to number of MMs <" << m << ">." << endl;
      }//if
      err = this->ComputeParameters(p, arrPM, weiPM, m, arrMM, weiMM, pars);
      if (err == errNoErr) this->Adjust(p, arrPM, pars);
   } else if ((p == 0) || (m == 0)) {
      cerr << "Error: Number of PMs or MMs is zero." << endl;
      err = errAbort;
      goto cleanup;
   } else {
      cerr << "Error: Option <" << fOption.Data() << "> is not applicable." << endl;
      err = errAbort;
      goto cleanup;
   }//if

// Fill output arrays: x (bg-corrected intensity), y (background), z (stdev)
   p = 0;
   m = 0;
   for (Int_t i=0; i<n; i++) {
      if (msk[i] == 1) {
         y[i] = x[i] - arrPM[p];
         x[i] = arrPM[p];
         z[i] = pars[2];
         p++;
      } else if (msk[i] == 0) {
         y[i] = x[i] - arrMM[m];
         x[i] = arrMM[m];
         z[i] = pars[2];
         m++;
      }//if
   }//for_i

// Cleanup
cleanup:
   if (weiMM) {delete [] weiMM; weiMM = 0;}
   if (arrMM) {delete [] arrMM; arrMM = 0;}
   if (weiPM) {delete [] weiPM; weiPM = 0;}
   if (arrPM) {delete [] arrPM; arrPM = 0;}

   return err;
}//Calculate

//______________________________________________________________________________
Int_t XRMABackground::ComputeParameters(Int_t n, Double_t *x, Double_t *w,
                      Double_t *pars)
{
   // Estimate the parameters for the background, they will be returned in *pars
   // with: pars[0] is alpha, pars[1] is mu, pars[2] is sigma.
   // Note: parameter estimates are same as those given by Bioconductor affy_1.1
   if(kCS) cout << "------XRMABackground::ComputeParameters(x)------" << endl;

   Int_t err  = errNoErr;
   Int_t npts = (fNPar > 0) ? (Int_t)fPars[0] : 16384;

   Double_t alpha, max, sigma, sum;
   Int_t    nlo, nhi, ntop;

// Init local array
   Double_t *arr = 0;
   if (!(arr = new (nothrow) Double_t[n])) return errInitMemory;
   for (Int_t i=0; i<n; i++) arr[i] = 0.0;
   
   max = TStat::MaxDensity(n, x, w, npts, fKernel);
   
   nlo = 0;
   for (Int_t i=0; i<n; i++) {
      if (x[i] < max) arr[nlo++] = x[i];
   }//for_i
 
   max = TStat::MaxDensity(nlo, arr, w, npts, fKernel);

// Estimate sigma 
   sum  = 0.0;
   ntop = 0;
   for (Int_t i=0; i<n; i++) {
      if (x[i] < max) {
         sum = sum + (x[i] - max)*(x[i] - max);
         ntop++;
      }//if
   }//for_i
   sigma = TMath::Sqrt(sum/(ntop - 1))*TMath::Sqrt(2.0);
 
   nhi = 0;
   for (Int_t i=0; i<n; i++) {
      if (x[i] > max) arr[nhi++] = x[i];
   }//for_i
 
// Estimate alpha 
   for (Int_t i=0; i<nhi; i++) arr[i] = arr[i] - max;
   alpha = (nhi > 0) ? 1.0 / TStat::MaxDensity(nhi, arr, w, npts, fKernel) : 0.0;
 
   pars[0] = alpha;
   pars[1] = max;
   pars[2] = sigma; 

// Cleanup
   if (arr) {delete [] arr; arr = 0;}

   return err;
}//ComputeParameters

//______________________________________________________________________________
Int_t XRMABackground::ComputeParameters(Int_t n1, Double_t *x1, Double_t *w1,
                      Int_t n2, Double_t *x2, Double_t *w2, Double_t *pars)
{
   // Estimate the parameters for the background, they will be returned in *pars
   // with: pars[0] is alpha, pars[1] is mu, pars[2] is sigma.
   // Note: parameter estimates are same as those given by Bioconductor affy_1.1.1
   if(kCS) cout << "------XRMABackground::ComputeParameters(x1,x2)------" << endl;

// Number of points for density
   Int_t npts = (fNPar > 0) ? (Int_t)fPars[0] : 16384;
   
   Double_t x1max = TStat::MaxDensity(n1, x1, w1, npts, fKernel);
   Double_t x2max = TStat::MaxDensity(n2, x2, w2, npts, fKernel);

// Estimate alpha 
   Double_t sum  = 0.0;
   Int_t    ntop = 0;
   for (Int_t i=0; i<n1; i++){
      if (x1[i] > x1max){
         sum = sum + (x1[i] - x1max);
         ntop++;
      }
   }//for_i
   pars[0] = ntop/sum;

// Estimate mu
   pars[1] = x2max;

// Estimate sigma 
   sum  = 0.0;
   ntop = 0;
   for (Int_t i=0; i<n2; i++){
      if (x2[i] < x2max){
         sum = sum + (x2[i] - x2max)*(x2[i] - x2max);
         ntop++;
      }//if
   }//for_i
   pars[2] = TMath::Sqrt(sum/(ntop - 1))*TMath::Sqrt(2.0) / 0.85;

   return errNoErr;
}//ComputeParameters

//______________________________________________________________________________
void XRMABackground::Adjust(Int_t n, Double_t *x, Double_t *pars)
{
   // Adjustment of array x where parameters are assumed to be:
   // pars[0] is alpha, pars[1] is mu, pars[2] is sigma.
   if(kCS) cout << "------XRMABackground::Adjust------" << endl;

   for (Int_t i=0; i<n; i++) {
      Double_t a = x[i] - pars[1] - pars[0]*pars[2]*pars[2];
      x[i] = (pars[2] != 0) 
           ? a + pars[2]*TMath::Gaus(a/pars[2],0,1,kTRUE)/TMath::Freq(a/pars[2])
           : a;   // protect against broken chips
//old      x[i] = a + pars[2]*TMath::Gaus(a/pars[2],0,1,kTRUE)/TMath::Freq(a/pars[2]);
//R      x[i] = a + pars[2]*TMLMath::DNorm(a/pars[2],0,1,kFALSE)/TMLMath::PNorm(a/pars[2],0,1,kTRUE,kFALSE);
   }//for_i
}//Adjust


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCBackground                                                        //
//                                                                      //
// Class for background algorithm based on similar GC content           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGCBackground::XGCBackground()
              :XBackgrounder()
{
   // Default GCBackground constructor
   if(kCS) cout << "---XGCBackground::XGCBackground(default)------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XGCBackground::XGCBackground(const char *name, const char *type)
              :XBackgrounder(name, type)
{
   // Normal GCBackground constructor
   if(kCS) cout << "---XGCBackground::XGCBackground------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XGCBackground::~XGCBackground()
{
   // GCBackground destructor
   if(kCS) cout << "---XGCBackground::~XGCBackground------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGCBackground::Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                     Int_t *msk)
{
   // Calculate background using background probes with GC-content identical
   // to PM probes intensities
   // Note: array msk contains the GC-content for the probes with:
   //       0 <= GC-content <= probe-length for PM probes 
   //       eINITMASK - 1 >= GC-content >= eINITMASK - probe-length - 1 for MM probes 
   //       msk = -1 for unused probes 
   if(kCS) cout << "------XGCBackground::Calculate------" << endl;

   Int_t err    = errNoErr;
   Int_t size   = kProbeLength+1;  //0 <= GC-content <= kProbeLength
   Int_t maxbin = 0;               //probesize for GC-content with largest number of probes
   Int_t mingc  = 0;               //minimum GC-content of probes
   Int_t maxgc  = kProbeLength;    //maximum GC-content of probes

   Int_t    len = 0;
   Double_t var = 0;

// Get trim parameter (default = 0.5  i.e. median)
   Double_t trim = (fNPar > 0) ? fPars[0] : 0.5;  //trim vlaue

// Init local arrays
   Int_t     *nbins = 0;  //number of bgrd probes with certain GC-content
   Double_t  *mean  = 0;  //mean value of bgrd probes with certain GC-content
   Double_t  *stdv  = 0;  //standard deviation of bgrd probes with certain GC-content
   Double_t **table = 0;  //table containing bgrd intensities for each GC-content

// Init memory local arrays
   if (!(nbins = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
   if (!(mean  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(stdv  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}

   for (Int_t k=0; k<size; k++) {
      nbins[k] = 0;
      mean[k] = stdv[k] = 0.0;
   }//for_k

// Fill bins with number of bgrd probes for each GC-content 
   for (Int_t i=0; i<n; i++) {
      if (msk[i] < eINITMASK) {
//no      if (msk[i] <= eINITMASK) {
//no         Int_t gcbin = eINITMASK - msk[i];
         //need to set msk+1, since negative msk was set to -(numGC+1)
         Int_t gcbin = eINITMASK - (msk[i] + 1);
         nbins[gcbin]++;
      }//if
   }//for_i

// Get bin with minimum/maximum GC-content
   for (Int_t k=0; k<size; k++) {
      if (nbins[k] != 0) {
         mingc = k;
         break;
      }//if
   }//for_k
   for (Int_t k=mingc; k<size; k++) {
      if (nbins[k] == 0) {
         maxgc = k - 1;
//?         break;
         continue;
      }//if
      maxgc = k;
   }//for_k

   if (XManager::fgVerbose) {
      cout << "      range of background GC-content: " << endl;
      cout << "         minimum GC-content is " << mingc << endl;
      cout << "         maximum GC-content is " << maxgc << endl;
   }//if

// Init memory for table containing bgrd intensities for each GC-content
   // get maximum number of bgrd probes to reserve memory for table
   maxbin = TMath::MaxElement(size, nbins);
   if (!(table = CreateTable(size, maxbin))) {err = errInitMemory; goto cleanup;}

// Fill table with intensity for each bin consecutively
   for (Int_t k=0; k<size; k++) nbins[k] = 0;  //reset nbins
   for (Int_t i=0; i<n; i++) {
      if (msk[i] < eINITMASK) {
         Int_t gcbin = eINITMASK - (msk[i] + 1);
         table[gcbin][nbins[gcbin]++] = x[i];
      }//if
   }//for_i

////////////////
//TO DO: need to check if ok when e.g. only GC=4 with nbins[4]=0
///////////////

// Calculate trimmed mean as bgrd estimate for probes with same GC-content
   mean[mingc] = TStat::Mean(nbins[mingc], &table[mingc][0], trim, var, len);
   stdv[mingc] = TMath::Sqrt(var);
   for (Int_t k=0; k<mingc; k++) {
      mean[k] = mean[mingc];
      stdv[k] = stdv[mingc];
   }//for_k
   for (Int_t k=mingc+1; k<maxgc; k++) {
      mean[k] = TStat::Mean(nbins[k], &table[k][0], trim, var, len);
      stdv[k] = TMath::Sqrt(var);
   }//for_k
   mean[maxgc] = TStat::Mean(nbins[maxgc], &table[maxgc][0], trim, var, len);
   stdv[maxgc] = TMath::Sqrt(var);
   for (Int_t k=maxgc+1; k<size; k++) {
      mean[k] = mean[maxgc];
      stdv[k] = stdv[maxgc];
   }//for_k

// Return background array y and return stdev in array z
   for (Int_t i=0; i<n; i++) {
      if (msk[i] < eINITMASK) {
         Int_t gcbin = eINITMASK - (msk[i] + 1);
         y[i] = mean[gcbin];
         z[i] = stdv[gcbin];
      } else if (msk[i] >= 0) {
         Int_t gcbin = msk[i];
         y[i] = mean[gcbin];
         z[i] = stdv[gcbin];
      }//if
   }//for_i

// Return background-corrected intensity in array x
   x = AdjustIntensity(n, x, y, z);

cleanup:
   // delete table
   DeleteTable(table, size);

   // delete arrays
   if (stdv)  {delete [] stdv;  stdv  = 0;}
   if (mean)  {delete [] mean;  mean  = 0;}
   if (nbins) {delete [] nbins; nbins = 0;}

   return err;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMeanDifferenceCall                                                  //
//                                                                      //
// Present call algorithm based on difference of means                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMeanDifferenceCall::XMeanDifferenceCall()
                    :XCallDetector()
{
   // Default PresentCall constructor
   if(kCS) cout << "---XMeanDifferenceCall::XMeanDifferenceCall(default)------" << endl;

   fNDefPar = 4;
   fMean    = 0;
}//Constructor

//______________________________________________________________________________
XMeanDifferenceCall::XMeanDifferenceCall(const char *name, const char *type)
                    :XCallDetector(name, type)
{
   // Normal PresentCall constructor
   if(kCS) cout << "---XMeanDifferenceCall::XMeanDifferenceCall------" << endl;

   fNDefPar = 4;
   fMean    = new XArithmeticMean("ArithmeticMean", type);
}//Constructor

//______________________________________________________________________________
XMeanDifferenceCall::~XMeanDifferenceCall()
{
   // PresentCall destructor
   if(kCS) cout << "---XMeanDifferenceCall::~XMeanDifferenceCall------" << endl;

   SafeDelete(fMean);
}//Destructor

//______________________________________________________________________________
Int_t XMeanDifferenceCall::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate present call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XMeanDifferenceCall::Calculate------" << endl;

// Test number of parameters
   if (TestNumParameters(4) != errNoErr) return errInitParameters;

// Get parameters
   Double_t trim    = fPars[0];
   Double_t trim1   = fPars[1];
   Double_t trim2   = fPars[2];
   Double_t percent = fPars[3];
   Double_t max     = fTreeInfo->GetValue("fMaxInten");

   Int_t numpar = 1;
   Double_t *ptr;

   ptr = &trim;
   Double_t meanPM = this->GetMean(numpar, ptr, fLength, fInten1);
   if (meanPM <= 0.0) {
      value1 = 0;
      value2 = 0;  //pvalue???
      return errNoErr;
   }//if

   ptr = &trim1;
   Double_t meanPM1 = this->GetMean(numpar, ptr, fLength, fInten1);
   Double_t meanMM1 = this->GetMean(numpar, ptr, fLength, fInten2);

   ptr = &trim2;
   Double_t meanPM2 = this->GetMean(numpar, ptr, fLength, fInten1);
   Double_t meanMM2 = this->GetMean(numpar, ptr, fLength, fInten2);

   Bool_t conditionA = ((meanPM1 - meanMM1) < TMath::Abs(meanPM1 * percent));
   Bool_t conditionB = ((meanPM2 - meanMM2) < TMath::Abs(meanPM2 * percent));
   Bool_t conditionC = ((meanPM1 >= max) && (meanMM1 <= meanPM1)); //not necessary?

   Int_t call = 2;   //moved from below
   if (conditionC) {
      call = 2;  //'P'  //not necessary?
   } else if (conditionA && conditionB) {
      call = 0;  //'A'
   } else if ((conditionA && !conditionB) || (!conditionA && conditionB)) {
      call = 1;  //'M'
   }//if

// Return values
   value1 = call;
   value2 = 0;  //pvalue???
   return errNoErr;
}//Calculate

//______________________________________________________________________________
Double_t XMeanDifferenceCall::GetMean(Int_t npar, Double_t *pars, 
                              Int_t length, Double_t *array)
{
   // Calculate mean for present call 
if(kCSa) cout << "------XMeanDifferenceCall::GetMean------" << endl;

   Double_t mean = 0;
   Double_t var  = 0;
   Int_t    num  = 0;

   fMean->InitParameters(npar, pars); 
   fMean->InitArrays(length,fInten1,fStdev1,fNPix1,fInten2,fStdev2,fNPix2);
   fMean->SetArray(length, array);
   fMean->Calculate(mean, var, num);
   fMean->DeleteArray();
   fMean->DeleteParameters();

   return mean;
}//FGetMean


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDetectionCall                                                       //
//                                                                      //
// Detection call algorithm                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDetectionCall::XDetectionCall()
               :XCallDetector()
{
   // Default DetectionCall constructor
   if(kCS) cout << "---XDetectionCall::XDetectionCall(default)------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XDetectionCall::XDetectionCall(const char *name, const char *type)
               :XCallDetector(name, type)
{
   // Normal DetectionCall constructor
   if(kCS) cout << "---XDetectionCall::XDetectionCall------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XDetectionCall::~XDetectionCall()
{
   // DetectionCall destructor
   if(kCS) cout << "---XDetectionCall::~XDetectionCall------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XDetectionCall::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate detection call
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XDetectionCall::Calculate------" << endl;

// Test number of parameters
   if (TestNumParameters(4) != errNoErr) return errInitParameters;

// Get parameters
   Double_t tau      = fPars[0];
   Double_t alpha1   = fPars[1];
   Double_t alpha2   = fPars[2];
   Double_t saturate = (fPars[3] > 0.0) ?  46000.0 : 0.0;
//??   Double_t saturate = (fPars[3] > 0.0) ?  fTreeInfo->GetValue("fMaxInten") : 0.0;

//??? for Wilcox in TStatUtils
//   Bool_t exact   = (fNPar > 4) ? fPars[4] : 0;
//   Bool_t correct = (fNPar > 5) ? fPars[5] : 1;

// Correct for saturated MMs (fInten2)
   Int_t total = 0;
   Int_t last  = 0;
   if (saturate > 0.0) {
      Int_t *ignore = new Int_t[fLength];

      total = 0;
      for (Int_t i=0; i<fLength; i++) {
         if (fInten2[i] > saturate) {
	         ignore[i] = 1;
	         total++;
         } else {
            ignore[i] = 0;
         }//if
      }//for_i

      // ignore probes with saturated MMs unless all are saturated
      last = 0;
      if ((total > 0) & (total < fLength)) {
         for (Int_t i=0; i<fLength; i++) {
	         if (!ignore[i]) {
	            fInten1[last] = fInten1[i];  //PM
	            fInten2[last] = fInten2[i];  //MM
	            last++;
	         }//if
         }//for_i

         fLength = last;
      }//if

      delete [] ignore;
   }//if

   Double_t *score = new Double_t[fLength];
   for (Int_t i=0; i<fLength; i++){   
      score[i] = (fInten1[i] - fInten2[i]) / (fInten1[i] + fInten2[i]);
   }//for_i

// Detection p-value (value2)
   value2 = this->WilcoxTest(fLength, score, tau);
//   value2 = TStat::WilcoxTest(fLength, score, tau, exact, correct, ???);

// Present call (value1)
   if      (value2 < alpha1) value1 = 2.0;  //"P"
   else if (value2 < alpha2) value1 = 1.0;  //"M"
   else                      value1 = 0.0;  //"A"

   delete [] score;

   return errNoErr;
}//Calculate

//______________________________________________________________________________
Double_t XDetectionCall::WilcoxTest(Int_t n, Double_t *x, Double_t mu)
{
   // Translation of relevant bits of the wilcox.test method from R
   if(kCSa) cout << "------XDetectionCall::WilcoxTest------" << endl;

   Int_t err = errNoErr;

   Int_t    prev, ntie;
   Double_t stat, nties, z, sigma, pval;    

// Get ridd of zero values
   Int_t j = 0;
   for(Int_t i=0; i<n; i++) {
      x[j] = x[i] - mu;
      if (x[j] != 0) j++;
   }//for_i
   n = j;

   if (n == 0) return 1; //pval=1.0

// Init local arrays
   Double_t *rank  = 0;
   Double_t *absx  = 0;
   Int_t    *index = 0;
   if (!(rank  = new (nothrow) Double_t[n])) {err = errInitMemory; goto cleanup;}
   if (!(absx  = new (nothrow) Double_t[n])) {err = errInitMemory; goto cleanup;}
   if (!(index = new (nothrow) Int_t[n]))    {err = errInitMemory; goto cleanup;}

   for (Int_t i=0; i<n; i++) {
      absx[i] = TMath::Abs(x[i]);
   }//for_i

// Get index and rank for absx
   TStat::TrueRank1(n, absx, index, rank);
   for (Int_t i=0; i<n; i++) {
      rank[i] = (x[index[i]] > 0) ? rank[i] : -rank[i];
   }//for_i

   stat = 0.0;    
   for (Int_t i=0; i<n; i++) {
      if (rank[i] > 0) {
         stat += rank[i];
      }//if
   }//for_i

   nties = 0.0;
   prev  = 0;
   ntie  = 0;
   for (Int_t i=1; i<n; i++) {
      if (rank[prev] == rank[i]) {
         ntie++;
      } else {
         if (ntie > 1) {
            nties += ntie*ntie*ntie - ntie; 
         }//if
         ntie = 0;
         prev = i;
      }//if
   }//for_i

   // need to convert to double otherwise large numbers create nan
   z     = stat - (Double_t)n*(n + 1)/4.0;
   sigma = TMath::Sqrt((Double_t)n*(n + 1)*(2*n + 1)/24.0 - nties/48.0);
   pval  = TStat::PNormApprox(z/sigma);

// Cleanup
cleanup:
   if (index) delete [] index;
   if (absx)  delete [] absx;
   if (rank)  delete [] rank;

   return ((err == errNoErr) ? (1.0 - pval) : (Double_t)err);
}//WilcoxTest


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMAS4Call                                                            //
//                                                                      //
// Prsent call algorithm of MAS4                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMAS4Call::XMAS4Call()
          :XCallDetector()
{
   // Default MAS4Call constructor
   if(kCS) cout << "---XMAS4Call::XMAS4Call(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XMAS4Call::XMAS4Call(const char *name, const char *type)
          :XCallDetector(name, type)
{
   // Normal MAS4Call constructor
   if(kCS) cout << "---XMAS4Call::XMAS4Call------" << endl;

}//Constructor

//______________________________________________________________________________
XMAS4Call::~XMAS4Call()
{
   // MAS4Call destructor
   if(kCS) cout << "---XMAS4Call::~XMAS4Call------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XMAS4Call::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate present call of MAS4
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XMAS4Call::Calculate------" << endl;

   cout << "Error: MAS4 Present Call is not yet implemented." << endl;
   return errAbort;

//   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDABGCall                                                            //
//                                                                      //
// Detection Above BackGround call algorithm                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDABGCall::XDABGCall()
          :XCallDetector()
{
   // Default DABGCall constructor
   if(kCS) cout << "---XDABGCall::XDABGCall(default)------" << endl;

   fNDefPar = 3;
   fNMaxGC  = 0;
   fNMaxBin = 0;
   fNBins   = 0;
   fGCTable = 0;
}//Constructor

//______________________________________________________________________________
XDABGCall::XDABGCall(const char *name, const char *type)
          :XCallDetector(name, type)
{
   // Normal DABGCall constructor
   if(kCS) cout << "---XDABGCall::XDABGCall------" << endl;

   fNDefPar = 3;
   fNMaxGC  = kProbeLength+1;
   fNMaxBin = 0;
   fGCTable = 0;

   fNBins = new Int_t[fNMaxGC];
   for (Int_t k=0; k<fNMaxGC; k++) fNBins[k] = 0;
}//Constructor

//______________________________________________________________________________
XDABGCall::~XDABGCall()
{
   // DABGCall destructor
   if(kCS) cout << "---XDABGCall::~XDABGCall------" << endl;

   DeleteTable(fGCTable, fNMaxGC);
   if (fNBins) {delete [] fNBins; fNBins = 0;}
}//Destructor

//______________________________________________________________________________
Int_t XDABGCall::Calculate(Int_t n, Double_t *x, Int_t *msk)
{
   // Get background intensities and save in table fGCTable sorted for GC-content
   // Note: array msk contains the GC-content for the probes with:
   //       eINITMASK - 1 >= GC-content >= eINITMASK - probe-length - 1 for MM probes 
   if(kCS) cout << "------XDABGCall::Calculate(table)------" << endl;

// Fill bins with number of bgrd probes for each GC-content 
   for (Int_t k=0; k<fNMaxGC; k++) fNBins[k] = 0;
   for (Int_t i=0; i<n; i++) {
      if (msk[i] < eINITMASK) {
//no      if (msk[i] <= eINITMASK) {
//no         Int_t gcbin = eINITMASK - msk[i];
         //need to set msk+1, since negative msk was set to -(numGC+1)
         Int_t gcbin = eINITMASK - (msk[i] + 1);
         fNBins[gcbin]++;
      }//if
   }//for_i

// Get bin with minimum/maximum GC-content
   Int_t mingc = 0;               //minimum GC-content of probes
   Int_t maxgc = kProbeLength;    //maximum GC-content of probes
   for (Int_t k=0; k<fNMaxGC; k++) {
      if (fNBins[k] != 0) {
         mingc = k;
         break;
      }//if
   }//for_k
   for (Int_t k=mingc; k<fNMaxGC; k++) {
      if (fNBins[k] == 0) {
         maxgc = k - 1;
//?         break;
         continue;
      }//if
      maxgc = k;
   }//for_k

   if (XManager::fgVerbose) {
      cout << "      range background of GC-content: " << endl;
      cout << "         minimum GC-content is " << mingc << endl;
      cout << "         maximum GC-content is " << maxgc << endl;
   }//if

// Init memory for table containing intensities for each GC-content
   // get maximum number of probes to reserve memory for table
   fNMaxBin = TMath::MaxElement(fNMaxGC, fNBins);
   DeleteTable(fGCTable, fNMaxGC);
   if (!(fGCTable = CreateTable(fNMaxGC, fNMaxBin))) return errInitMemory;

// Fill table with bgrd intensity for each bin consecutively
   for (Int_t k=0; k<fNMaxGC; k++) fNBins[k] = 0;  //reset nbins
   for (Int_t i=0; i<n; i++) {
      if (msk[i] < eINITMASK) {
         Int_t gcbin = eINITMASK - (msk[i] + 1);
         fGCTable[gcbin][fNBins[gcbin]++] = x[i];
      }//if
   }//for_i

   return errNoErr;
}//Calculate

//______________________________________________________________________________
Int_t XDABGCall::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate detection above background call, where value1 returns the
   // detection call and value2 returns the DABG p-value
   // Note: array fInten1 contains intensities of the current probeset and
   //       array fNPix1 is used to store the corresponding GC-content
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XDABGCall::Calculate------" << endl;

// Test number of parameters
   if (TestNumParameters(3) != errNoErr) return errInitParameters;

// Get parameters
   Double_t cut    = fPars[0];
   Double_t alpha1 = fPars[1];
   Double_t alpha2 = fPars[2];

// Determine if to use Fisher or Percentile to calculate p-values
   Bool_t doPercentile = (cut >= 0.0) && (cut <= 1.0);

// Compute p-value
   if (doPercentile) {
      value2 = this->PValuePercentile(fLength, fNPix1, fInten1, cut);
   } else {
      value2 = this->PValueFisher(fLength, fNPix1, fInten1);
   }//if

// Present call (value1)
   if      (value2 < alpha1) value1 = 2.0;  //"P"
   else if (value2 < alpha2) value1 = 1.0;  //"M"
   else                      value1 = 0.0;  //"A"

   return errNoErr;
}//Calculate

//______________________________________________________________________________
Double_t XDABGCall::PValueFisher(Int_t n, Int_t *arrgc, Double_t *inten)
{
   // Calculate p-value using Fisher's chi-squared method
   // Note: Adapted from Affymetrix APT: file DABG.cpp
   if(kCSa) cout << "------XDABGCall::PValueFisher------" << endl;

   Double_t pvalue = 1.0;

// Calculate p-value product
   Double_t pvalpro = 0.0;
   for (Int_t i=0; i<n; i++) {
      pvalue  = this->Intensity2PValue(arrgc[i], inten[i]);
      pvalpro = pvalpro + log(pvalue);
   }//for_i

// Set minimum value: why?
   if (pvalpro == 0.0) pvalpro = 0.000001;

// Get chi-squared probability
   Double_t stat = -2*pvalpro;  //why -2*
   pvalue = this->ChiSqrProb(2*n, (Float_t)stat);

   return pvalue;
}//PValueFisher

//______________________________________________________________________________
Double_t XDABGCall::PValuePercentile(Int_t n, Int_t *arrgc, Double_t *inten, Double_t cut)
{
   // Calculate p-value using an individual probe p-value at a particular percentile
   // Note: Adapted from Affymetrix APT: file DABG.cpp
   if(kCSa) cout << "------XDABGCall::PValuePercentile------" << endl;

   Double_t pvalue = 1.0;

//TO DO: replace code with non-template code!!
// Fill p-values
   std::vector<Double_t> vecpval;
   vecpval.reserve(n);
   for (Int_t i=0; i<n; i++) {
      pvalue = this->Intensity2PValue(arrgc[i], inten[i]);
      vecpval.push_back(pvalue);
   }//for_i

   if (n == 1) {
      pvalue = vecpval[0];
   } else {
      Double_t pos = cut*(n - 1);
      Int_t    idx = (Int_t)pos;

      // lookup the pvalue (and interpolate)
      sort(vecpval.begin(), vecpval.end());
      pvalue = vecpval[idx] + ((pos - idx)*(vecpval[idx+1] - vecpval[idx]));
   }//if

   return pvalue;
}//PValuePercentile

//______________________________________________________________________________
Double_t XDABGCall::Intensity2PValue(Int_t gcbin, Double_t inten)
{
   // Calculate p-value for intensity inten with GC-content gcbin
   // Note: Adapted from Affymetrix APT: file DABG.cpp
   if(kCSa) cout << "------XDABGCall::Intensity2PValue------" << endl;

   Double_t pvalue = 1.0;
   Int_t    size   = fNBins[gcbin]; //number of probes with GC-content gcbin

   if (size == 0) {
      pvalue = 1.0;
   } else if (size == 1) {
      Double_t bgrd = fGCTable[gcbin][0];

      if (inten < bgrd) {
         pvalue = 1.0;
      } else if (inten > bgrd) {
         pvalue = 0.0;
      } else {
         pvalue = 0.5;
      }//if
   } else {
//TO DO: replace code with non-template code!!
      // for function lower_bound() need to convert array to template vector
      std::vector<Double_t> vecinten;
      vecinten.reserve(size);
      for (Int_t i=0; i<size; i++) {
         vecinten.push_back(fGCTable[gcbin][i]);
      }//for_i

      // since fGCTable is not sorted
      sort(vecinten.begin(),vecinten.end());

      // get position of inten in corresponding sorted bgrd intensities
      Int_t idx = lower_bound(vecinten.begin(),vecinten.end(),inten)
                - vecinten.begin();  //need to subtract start position

      if (idx == size) idx = size - 1;
      pvalue = 1.0 - ((Double_t)idx)/((Double_t)size);
   }//if

   return pvalue;
}//Intensity2PValue

//______________________________________________________________________________
template <typename T1> T1 XDABGCall::ChiSqrProb(Int_t n, T1 x)
{
   // Fisher's chi-squared probability
   // Note: Adapted from Affymetrix APT: file stats-distributions.h
   if(kCSa) cout << "------XDABGCall::ChiSqrProb------" << endl;

   Double_t p = 1.0;

   if (x <= 0) {
      p = 1.0;
   } else if (n > 100) {
      p = UProb((pow((Double_t)(x/n),(Double_t)(1.0/3.0))
        - (1.0 - 2.0/9.0/n))/sqrt(((2.0/9.0)/n)));
   } else if (x > 400) {
      p = 0.0;
   } else {  
      Double_t a, i, k;

      if ((n%2) != 0) {
         p = 2.0*UProb(sqrt(x));
         a = sqrt((Float_t)(2.0/TMath::Pi()))*exp(-x/2)/sqrt(x);
         k = 1;
      } else {
         p = a = exp(-x/2.0);
         k = 2;
      }//if

      for (i=k; i<=(n-2); i+=2) {
         a *= x/i;
         p += a;
      }//for_i
   }//if

   return (T1)p;
}//ChiSqrProb

//______________________________________________________________________________
template <typename T1> T1 XDABGCall::UProb(T1 x)
{
   // Upper probability of the u distribution (u=-0.85)
   // Note: Adapted from Affymetrix APT: file stats-distributions.h
   if(kCSa) cout << "------XDABGCall::UProb------" << endl;

   Double_t p    = 0.0;
   Double_t absx = fabs(x);

   if (absx < 1.9) {
      p = pow((1.0
        + absx*(0.049867347
              + absx*(0.0211410061
                    + absx*(0.0032776263
                          + absx*(0.0000380036
                                + absx*(0.0000488906
                                      + absx*0.000005383)))))), -16)/2; // ?
   } else if (absx <= 100.0) {
      for (Int_t i=18; i>=1; i--) {
         p = i/(absx + p);
      }//for_i
      p = exp(-0.5*absx*absx) / sqrt(2.0*TMath::Pi())/(absx + p);
   }//if

   if (x < 0) p = 1.0 - p;

   return (T1)p;
}//UProb


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XINICall                                                             //
//                                                                      //
// Informative call algorithm of FARMS                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XINICall::XINICall()
         :XCallDetector()
{
   // Default INICall constructor
   if(kCS) cout << "---XINICall::XINICall(default)------" << endl;

   fNDefPar   = 8;
   fMultichip = kTRUE;
}//Constructor

//______________________________________________________________________________
XINICall::XINICall(const char *name, const char *type)
         :XCallDetector(name, type)
{
   // Normal INICall constructor
   if(kCS) cout << "---XINICall::XINICall------" << endl;

   fNDefPar   = 8;
   fMultichip = kTRUE;
}//Constructor

//______________________________________________________________________________
XINICall::~XINICall()
{
   // INICall destructor
   if(kCS) cout << "---XINICall::~XINICall------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XINICall::SetArray(Int_t length, Double_t *array)
{
   // Set size of array for calculation to length doubles and set contents
   // and convert array to logarithm of base "log2"
   if(kCSa) cout << "------XINICall::SetArray------" << endl;

   Int_t err = XHybridizer::SetArray(length, array);
   if (err != errNoErr) return err;

   fArray = Array2Log(fLength, fArray, 1.0, "log2");

   return err;
}//SetArray

//______________________________________________________________________________
Int_t XINICall::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Calculate Informative/NonInformative call of FARMS
   // Note: Informative call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XINICall::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t nrow = (Int_t)(fLength / n);
   Int_t ncol = n; 

// Get parameters
   Int_t    version  = (Int_t)(fPars[0]);  // version of package farms 
   Double_t weight   = fPars[1];           // hyperparameter
   Double_t mu       = fPars[2];           // hyperparameter
   Double_t scale    = fPars[3];           // scaling parameter
   Double_t tol      = fPars[4];           // termination tolerance of EM
   Int_t    cyc      = (Int_t)(fPars[5]);  // maximum number of cycles of EM
   Double_t alpha1   = fPars[6];
   Double_t alpha2   = fPars[7];

   if (cyc <= 0) cyc = 2*ncol;

   if (version == 131) {
      err = this->DoFARMS131(nrow, ncol, fArray, x, y, weight, mu, scale, tol, cyc);
   } else if (version == 130) {
      err = this->DoFARMS130(nrow, ncol, fArray, x, y, weight, mu, scale, tol, cyc);
   } else {
      cerr << "Error: Version <" << version << "> is not supported." << endl;
      err = errAbort;
   }//if

// Present call (x)
   for (Int_t k=0; k<ncol; k++) {
      if      (y[k] < alpha1) x[k] = 2.0;  //"P"
      else if (y[k] < alpha2) x[k] = 1.0;  //"M"
      else                    x[k] = 0.0;  //"A"
   }//for_k
   
   return err;
}//Calculate

//______________________________________________________________________________
Int_t XINICall::DoFARMS130(Int_t nrow, Int_t ncol, Double_t *inten, Double_t *x, 
                Double_t *y, Double_t weight, Double_t mu, Double_t scale, 
                Double_t tol, Double_t cyc)
{
   // Calculate FARMS (Factor Analysis for Robust Microarray Summarization)
   // Reference: Hochreiter S, et al., Bioinformatics 2006 Apr 15;22(8):943-9.
   // Note: Adapted from R package farms_1.3: file farms.R (license GNU GPL >= 2)
   if(kCSa) cout << "------XINICall::DoFARMS130------" << endl;

   Int_t    err    = errNoErr;
   Double_t sum    = 0.0;
   Double_t sumL   = 0.0;
   Double_t alpha  = 0.0;
   Double_t bbeta  = 0.0;
   Double_t mean   = 0.0;
   Double_t a      = 0.0;
   Double_t ezz    = 0.0;

// Initialize local arrays
   Double_t *xmean  = 0;  // column means
   Double_t *X      = 0;  // matrix X = x - xmean
   Double_t *mm     = 0;  // temporary matrix
   Double_t *XX     = 0;  // matrix (positive definit)
   Double_t *diagXX = 0;  // diagonale(XX)
   Double_t *L      = 0;  // factor loadings
   Double_t *Lold   = 0;  // temporary L
   Double_t *Ph     = 0;  // diagonal uniqueness matrix
   Double_t *PsiL   = 0;  // 
   Double_t *beta   = 0;  // 
   Double_t *XXbeta = 0;  // 
   Double_t *TXbeta = 0;  // 
   Double_t *EZZ    = 0;  //

// Initialize memory for local arrays
   if (!(X  = new (nothrow) Double_t[nrow*ncol])) {err = errInitMemory; goto cleanup;}
   if (!(mm = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}
   if (!(XX = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}

   if (!(xmean  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(diagXX = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(L      = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Lold   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Ph     = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(PsiL   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(beta   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(XXbeta = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(TXbeta = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(EZZ    = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}

   // row means
   for (Int_t i=0; i<nrow; i++) {
      xmean[i] = 0.0;
      for (Int_t j=0; j<ncol; j++) {
         xmean[i] += inten[i*ncol + j];
      }//for_j
      xmean[i] = xmean[i]/ncol;
   }//for_i

   // center data: X = x - xmean
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = inten[i*ncol + j] - xmean[i];
      }//for_i
   }//for_j

   // mm = t(X)*X/ncol
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            sum += X[i*ncol + j]*X[k*ncol + j];
         }//for_j
         mm[i*nrow + k] = sum/ncol;
      }//for_k
   }//for_i

   // XX = (mm + t(mm))/2 but positive definit
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = (mm[i*nrow + k] + mm[k*nrow + i])/2.0;
         XX[i*nrow + k] = (sum > 0.0) ? sum : 0.0;
      }//for_k
   }//for_i

   // diagonale of XX
   mean = 0.0;
   for (Int_t i=0; i<nrow; i++) {
      diagXX[i] = XX[i*nrow + i];
      mean     += diagXX[i];
   }//for_i
   mean  = mean/nrow;
   alpha = weight/mean;
   bbeta = mu*alpha;

   // factor loadings: L
   for (Int_t i=0; i<nrow; i++) {
      L[i]    = 0.8*diagXX[i];
      Lold[i] = 0.0;
   }//for_i

   // diagonal uniqueness matrix: Ph
   for (Int_t i=0; i<nrow; i++) {
      Ph[i] = 0.05*diagXX[i];
   }//for_i

//  Expectation-Maximization (EM) algorithm
   while (cyc >0) {
   // E step
      a = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         PsiL[i] = L[i]/Ph[i];
         a      += L[i]*PsiL[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         beta[i] = PsiL[i]/a;
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         sum = 0.0;
         for (Int_t k=0; k<nrow; k++) {
            sum += XX[i*nrow + k]*beta[k];
         }//for_k
         XXbeta[i] = sum;
      }//for_i

      sum = 0.0;
      ezz = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         sum += beta[i]*L[i];
         ezz += beta[i]*XXbeta[i];
      }//for_i
      ezz = ezz - sum;

      for (Int_t i=0; i<nrow; i++) {
         TXbeta[i] = XXbeta[i] + bbeta*Ph[i];
         EZZ[i]    = ezz       + alpha*Ph[i];
      }//for_i

   // M step
      for (Int_t i=0; i<nrow; i++) {
         L[i] = TXbeta[i]/EZZ[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         Ph[i] = diagXX[i] - XXbeta[i]*L[i] +  alpha*Ph[i]*L[i]*(bbeta - L[i]);
      }//for_i

      sumL = sum = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         sumL += L[i]*L[i];
         sum  += Lold[i]*Lold[i];
      }//for_i
      sumL = TMath::Abs(TMath::Sqrt(sumL) - TMath::Sqrt(sum));

      if (sumL < tol) break;

      for (Int_t i=0; i<nrow; i++) Lold[i] = L[i];
      cyc--;
   }//while

// Final result x and stdev y
   for (Int_t j=0; j<ncol; j++) {

      y[j] = 1.0/a;
   }//for_j


// Cleanup
cleanup:
   if (EZZ)    {delete [] EZZ;    EZZ    = 0;}
   if (TXbeta) {delete [] TXbeta; TXbeta = 0;}
   if (XXbeta) {delete [] XXbeta; XXbeta = 0;}
   if (beta)   {delete [] beta;   beta   = 0;}
   if (PsiL)   {delete [] PsiL;   PsiL   = 0;}
   if (Ph)     {delete [] Ph;     Ph     = 0;}
   if (Lold)   {delete [] Lold;   Lold   = 0;}
   if (L)      {delete [] L;      L      = 0;}
   if (diagXX) {delete [] diagXX; diagXX = 0;}
   if (xmean)  {delete [] xmean;  xmean  = 0;}
   if (XX)     {delete [] XX;     XX     = 0;}
   if (mm)     {delete [] mm;     mm     = 0;}
   if (X)      {delete [] X;      X      = 0;}
   
   return err;
}//DoFARMS130

//______________________________________________________________________________
Int_t XINICall::DoFARMS131(Int_t nrow, Int_t ncol, Double_t *inten, 
                Double_t *x, Double_t *y, Double_t weight, Double_t mu, 
                Double_t scale, Double_t tol, Double_t cyc)
{
   // Calculate FARMS (Factor Analysis for Robust Microarray Summarization)
   // Reference: Hochreiter S, et al., Bioinformatics 2006 Apr 15;22(8):943-9.
   // Note: Adapted from R package farms_1.3.1: file farms.R (license GNU GPL >= 2)
   if(kCSa) cout << "------XINICall::DoFARMS131------" << endl;

   Int_t    err    = errNoErr;
   Double_t sum    = 0.0;
   Double_t var    = 0.0;
   Double_t alpha  = 0.0;
   Double_t bbeta  = 0.0;
   Double_t mean   = 0.0;
   Double_t a      = 0.0;
   Double_t ezz    = 0.0;

// Initialize local arrays
   Double_t *xmean  = 0;  // column means
   Double_t *xstdv  = 0;  // standard deviation
   Double_t *X      = 0;  // matrix X = x - xmean
   Double_t *mm     = 0;  // temporary matrix
   Double_t *XX     = 0;  // matrix (positive definit)
   Double_t *diagXX = 0;  // diagonale(XX)
   Double_t *L      = 0;  // factor loadings
   Double_t *Lold   = 0;  // temporary L
   Double_t *Ph     = 0;  // diagonal uniqueness matrix
   Double_t *PsiL   = 0;  // 
   Double_t *beta   = 0;  // 
   Double_t *arr1   = 0;  // 
   Double_t *arr2   = 0;  // 
   Double_t *EZZ    = 0;  //
   Double_t *c      = 0;  //

// Initialize memory for local arrays
   if (!(X  = new (nothrow) Double_t[nrow*ncol])) {err = errInitMemory; goto cleanup;}
   if (!(mm = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}
   if (!(XX = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}

   if (!(xmean  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(xstdv  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(diagXX = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(L      = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Lold   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Ph     = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(PsiL   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(beta   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(arr1   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(arr2   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(EZZ    = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(c      = new (nothrow) Double_t[ncol])) {err = errInitMemory; goto cleanup;}

// Normalize
   // calculate row means
   for (Int_t i=0; i<nrow; i++) {
      xmean[i] = 0.0;
      for (Int_t j=0; j<ncol; j++) {
         xmean[i] += inten[i*ncol + j];
      }//for_j
      xmean[i] = xmean[i]/ncol;
   }//for_i

   // center data: XX = inten - xmean
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = inten[i*ncol + j] - xmean[i];
      }//for_i
   }//for_j

   // mm = X*t(X)/ncol
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            sum += X[i*ncol + j]*X[k*ncol + j];
         }//for_j
         mm[i*nrow + k] = sum/ncol;
      }//for_k
   }//for_i

   // calculate standard deviation
   for (Int_t i=0; i<nrow; i++) {
      xstdv[i] = TMath::Sqrt(mm[i*nrow + i]);
   }//for_i

   // standardize X to variance 1: X = inten/stdev
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = inten[i*ncol + j]/xstdv[i];
      }//for_i
   }//for_j

// XX
   // row means
   for (Int_t i=0; i<nrow; i++) {
      xmean[i] = 0.0;
      for (Int_t j=0; j<ncol; j++) {
         xmean[i] += X[i*ncol + j];
      }//for_j
      xmean[i] = xmean[i]/ncol;
   }//for_i

   // center data: X = X - xmean
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = X[i*ncol + j] - xmean[i];
      }//for_i
   }//for_j

   // mm = t(X)*X/ncol
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            sum += X[i*ncol + j]*X[k*ncol + j];
         }//for_j
         mm[i*nrow + k] = sum/ncol;
      }//for_k
   }//for_i

   // XX = (mm + t(mm))/2 but positive definit
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = (mm[i*nrow + k] + mm[k*nrow + i])/2.0;
         XX[i*nrow + k] = (sum > 0.0) ? sum : 0.0;
      }//for_k
   }//for_i

   // diagonale of XX
   for (Int_t i=0; i<nrow; i++) {
      diagXX[i] = XX[i*nrow + i];
   }//for_i

   alpha = weight*nrow;
   bbeta = mu*alpha;

   // factor loadings: L
   for (Int_t i=0; i<nrow; i++) {
      L[i]    = TMath::Sqrt(0.8*diagXX[i]);
      Lold[i] = 0.0;
   }//for_i

   // diagonal uniqueness matrix: Ph
   for (Int_t i=0; i<nrow; i++) {
      Ph[i] = diagXX[i] - L[i]*L[i];
   }//for_i

//  Expectation-Maximization (EM) algorithm
   while (cyc >0) {
   // E step
      a = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         PsiL[i] = L[i]/Ph[i];
         a      += L[i]*PsiL[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         beta[i] = PsiL[i]/a;
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         sum = 0.0;
         for (Int_t k=0; k<nrow; k++) {
            sum += XX[i*nrow + k]*beta[k];
         }//for_k
         arr1[i] = sum;
      }//for_i

      sum = 0.0;
      ezz = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         sum += beta[i]*L[i];
         ezz += beta[i]*arr1[i];
      }//for_i
      ezz = ezz - sum;

      for (Int_t i=0; i<nrow; i++) {
         arr2[i] = arr1[i] + bbeta*Ph[i];
         EZZ[i]  = ezz     + alpha*Ph[i];
      }//for_i

   // M step
      for (Int_t i=0; i<nrow; i++) {
         L[i] = arr2[i]/EZZ[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         Ph[i] = diagXX[i] - arr1[i]*L[i] +  alpha*Ph[i]*L[i]*(bbeta - L[i]);
      }//for_i

      sum = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         Double_t LL = Lold[i] - L[i];
         sum  += TMath::Sqrt(LL*LL);
      }//for_i
      if (sum < tol) break;

      for (Int_t i=0; i<nrow; i++) Lold[i] = L[i];
      cyc--;
   }//while

   sum = 0.0;
   for (Int_t j=0; j<ncol; j++) {
      c[j] = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         c[j] += X[i*ncol + j]*beta[i];
      }//for_i
      sum  += c[j];
   }//for_j
   mean = sum/ncol;

   var = TStat::StDev(ncol, c, mean);

   for (Int_t i=0; i<nrow; i++) {
      L[i] = L[i]*var;
   }//for_i

   a = 1.0;
   for (Int_t i=0; i<nrow; i++) {
      PsiL[i] = L[i]/Ph[i];
      a      += L[i]*PsiL[i];
   }//for_i

// Final expression and stdev
   for (Int_t j=0; j<ncol; j++) {
      y[j] = 1.0/a;
   }//for_j

// Cleanup
cleanup:
   if (c)      {delete [] c;      c      = 0;}
   if (EZZ)    {delete [] EZZ;    EZZ    = 0;}
   if (arr2)   {delete [] arr2;   arr2   = 0;}
   if (arr1)   {delete [] arr1;   arr1   = 0;}
   if (beta)   {delete [] beta;   beta   = 0;}
   if (PsiL)   {delete [] PsiL;   PsiL   = 0;}
   if (Ph)     {delete [] Ph;     Ph     = 0;}
   if (Lold)   {delete [] Lold;   Lold   = 0;}
   if (L)      {delete [] L;      L      = 0;}
   if (diagXX) {delete [] diagXX; diagXX = 0;}
   if (xstdv)  {delete [] xstdv;  xstdv  = 0;}
   if (xmean)  {delete [] xmean;  xmean  = 0;}
   if (XX)     {delete [] XX;     XX     = 0;}
   if (mm)     {delete [] mm;     mm     = 0;}
   if (X)      {delete [] X;      X      = 0;}
   
   return err;
}//DoFARMS131


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XArithmeticMean                                                      //
//                                                                      //
// Arithmetic mean expression algorithm                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XArithmeticMean::XArithmeticMean()
                :XExpressor()
{
   // Default ArithmeticMean constructor
   if(kCS) cout << "---XArithmeticMean::XArithmeticMean(default)------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XArithmeticMean::XArithmeticMean(const char *name, const char *type)
                :XExpressor(name, type)
{
   // Normal ArithmeticMean constructor
   if(kCS) cout << "---XArithmeticMean::XArithmeticMean------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XArithmeticMean::~XArithmeticMean()
{
   // ArithmeticMean destructor
   if(kCS) cout << "---XArithmeticMean::~XArithmeticMean------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XArithmeticMean::CreateArray(Int_t length)
{
   // Initialize array: copy intensity
   if(kCSa) cout << "------XArithmeticMean::CreateArray------" << endl;

   if (fArray) DeleteArray();
   if (!(fArray = new (nothrow) Double_t[length])) {return errInitMemory;}
   fLength = length;

   memcpy(fArray, fInten1, length*sizeof(Double_t));

   return errNoErr;
}//CreateArray

//______________________________________________________________________________
Int_t XArithmeticMean::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate arithmetic trimmed mean
   // value1 = mean; value2 = variance; num = length of trimmed array
   if(kCSa) cout << "------XArithmeticMean::Calculate------" << endl;

//Better!!
//!!! TStat::Mean()
//!!! TMath::Mean()

// Get trim parameter (default = 0)
   Double_t trim = (fNPar > 0) ? fPars[0] : 0.0;  //trim vlaue

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[fLength])) {return errInitMemory;}
   TMath::Sort(fLength, fArray, index);

   // start-index and end-index
   Int_t start, end;
   if (trim < 0.5) {
      start = (Int_t)TMath::Floor((Double_t)fLength * trim);
      end   = fLength - start;
   } else {
      if ((fLength % 2) == 0){
         start = (Int_t)TMath::Floor((Double_t)fLength / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)TMath::Ceil((Double_t)fLength / 2.0);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 0;
   for (Int_t i=start; i<end; i++) {
      mean += fArray[index[i]];
   }//for_i
   mean /= trimlen;

// Calculate variance
   Double_t var = 0;
   for (Int_t i=start; i<end; i++) {
      var += (fArray[index[i]] - mean)*(fArray[index[i]] - mean);
   }//for_i
   var /= trimlen;
//   var /= (trimlen - 1); // degrees of freedom? test for null division
//   var /= (fLength - 1); // degrees of freedom? test for null division

   if (index) {delete [] index; index = 0;}

// Return values
   value1 = mean;
   value2 = var;
   num    = trimlen;
   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeometricMean                                                       //
//                                                                      //
// Geometric mean expression algorithm                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGeometricMean::XGeometricMean()
               :XArithmeticMean()
{
   // Default GeometricMean constructor
   if(kCS) cout << "---XGeometricMean::XGeometricMean(default)------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XGeometricMean::XGeometricMean(const char *name, const char *type)
               :XArithmeticMean(name, type)
{
   // Normal GeometricMean constructor
   if(kCS) cout << "---XGeometricMean::XGeometricMean------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XGeometricMean::~XGeometricMean()
{
   // GeometricMean destructor
   if(kCS) cout << "---XGeometricMean::~XGeometricMean------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGeometricMean::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate geometric trimmed mean
   // value1 = mean; value2 = variance; num = length of trimmed array
   if(kCSa) cout << "------XGeometricMean::Calculate------" << endl;

//Better!!
//!!! TStat::GeoMean()
//!!! TMath::GeomMean()

// Get trim parameter (default = 0)
   Double_t trim = (fNPar > 0) ? fPars[0] : 0.0;  //trim vlaue

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[fLength])) {return errInitMemory;}
   TMath::Sort(fLength, fArray, index);

   // start-index and end-index
   Int_t start, end;
   if (trim < 0.5) {
      start = (Int_t)TMath::Floor((Double_t)fLength * trim);
      end   = fLength - start;
   } else {
      if ((fLength % 2) == 0){
         start = (Int_t)TMath::Floor((Double_t)fLength / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)TMath::Ceil((Double_t)fLength / 2.0);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 1;
   for (Int_t i=start; i<end; i++) {
      mean *= fArray[index[i]];
   }//for_i
   mean = pow(mean, 1.0/(Double_t)trimlen);

// Calculate variance
   Double_t var = 0;
   for (Int_t i=start; i<end; i++) {
      var += (fArray[index[i]] - mean)*(fArray[index[i]] - mean);
   }//for_i
   var /= trimlen;
//   var /= (trimlen - 1); // degrees of freedom? test for null division
//   var /= (fLength - 1); // degrees of freedom? test for null division

   if (index) {delete [] index; index = 0;}

// Return values
   value1 = mean;
   value2 = var;
   num    = trimlen;
   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XWeightedMean                                                        //
//                                                                      //
// Weighted mean expression algorithm                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XWeightedMean::XWeightedMean()
              :XArithmeticMean()
{
   // Default WeightedMean constructor
   if(kCS) cout << "---XWeightedMean::XWeightedMean(default)------" << endl;

   fNDefPar = 1;
}//Constructor

//______________________________________________________________________________
XWeightedMean::XWeightedMean(const char *name, const char *type)
              :XArithmeticMean(name, type)
{
   // Normal WeightedMean constructor
   if(kCS) cout << "---XWeightedMean::XWeightedMean------" << endl;

   fNDefPar = 1;
}//Constructor

//______________________________________________________________________________
XWeightedMean::~XWeightedMean()
{
   // WeightedMean destructor
   if(kCS) cout << "---XWeightedMean::~XWeightedMean------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XWeightedMean::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate weighted trimmed mean
   // value1 = mean; value2 = variance; num = length of trimmed array
   if(kCSa) cout << "------XWeightedMean::Calculate------" << endl;

//Better!!
//!!! TStat::Mean(..,w)
//!!! TMath::Mean(..,w)

// Test number of parameters
   if (TestNumParameters(1) != errNoErr) return errInitParameters;

// Get parameters
   Double_t wmaxinten = fPars[0];
   Double_t maxinten  = fTreeInfo->GetValue("fMaxInten");
   Double_t maxpix    = fTreeInfo->GetValue("fMaxNPixels");

   Double_t w1, w2, w3, value;
   Double_t *weight = 0;
   if (!(weight = new (nothrow) Double_t[fLength])) return errInitMemory;
   
// Calculate weighted mean
   Double_t mean      = 0;
   Double_t sumweight = 0;
   Int_t    len       = 0;
   for (Int_t i=0; i<fLength; i++) {
   // weight for ratio MM/PM
      if (fInten2[i] > fInten1[i]) {
         weight[i] = 0;
         continue;
      } else if ((fInten2[i] >= maxinten) && (fInten1[i] >= maxinten)) {
         w1 = wmaxinten;
      } else {
         w1 = fInten1[i] ? (1 - fInten2[i]/fInten1[i]) : 0; //omit zero division
      }//if
   // weight for StdDev[%]
      w2 = fInten1[i] ? (1 - fStdev1[i]/fInten1[i]) : 0; //omit zero division     
   // weight for ratio NumPixels/MaxNumPixels
      w3 = fNPix1[i]/maxpix;
   // overall weight
      weight[i] = w1 * w2 * w3; 
   
      sumweight += weight[i];
      mean += fArray[i] * weight[i];
      len++;
   }//for_i

// Calculate variance
   Double_t var = 0;
   if (sumweight > 0) {
      mean /= sumweight;
   // variance: Netz; Formeln der Mathematik; p.562
      for (Int_t i=0; i<fLength; i++) {
         value = fArray[i] - mean;
         var += value * value * weight[i];
      }//for_i 

      if (len > 1) {
         var /= (len - 1) * sumweight;
      } else {
         var = 0; //var = SP[?]*SP[?];
      }//if
   } else {
      mean = -1;
      var  =  0; //std = sqrt(vVar)
   }//if

   if (weight) {delete [] weight; weight = 0;}

// Return values
   value1 = mean;
   value2 = var;
   num    = len;
   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCCorrectedMean                                                     //
//                                                                      //
// Weighted mean expression algorithm with correction for GC-content    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGCCorrectedMean::XGCCorrectedMean()
                 :XArithmeticMean()
{
   // Default GCCorrectedMean constructor
   if(kCS) cout << "---XGCCorrectedMean::XGCCorrectedMean(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XGCCorrectedMean::XGCCorrectedMean(const char *name, const char *type)
                 :XArithmeticMean(name, type)
{
   // Normal GCCorrectedMean constructor
   if(kCS) cout << "---XGCCorrectedMean::XGCCorrectedMean------" << endl;

}//Constructor

//______________________________________________________________________________
XGCCorrectedMean::~XGCCorrectedMean()
{
   // GCCorrectedMean destructor
   if(kCS) cout << "---XGCCorrectedMean::~XGCCorrectedMean------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XGCCorrectedMean::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate weighted trimmed mean containing correction for GC-content
   // value1 = mean; value2 = variance; num = length of trimmed array
   if(kCSa) cout << "------XGCCorrectedMean::Calculate------" << endl;
   
// Calculate corrected mean
   Double_t mean = 0;
   Int_t    len  = 0;

// Calculate variance
   Double_t var = 0;

cout << "Note: GCCorrectedMean is not implemented yet!" << endl;

// Return values
   value1 = mean;
   value2 = var;
   num    = len;
   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XWeightedDiff                                                        //
//                                                                      //
// Weighted difference expression algorithm                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XWeightedDiff::XWeightedDiff()
              :XExpressor()
{
   // Default WeightedDiff constructor
   if(kCS) cout << "---XWeightedDiff::XWeightedDiff(default)------" << endl;

   fNDefPar = 1;
}//Constructor

//______________________________________________________________________________
XWeightedDiff::XWeightedDiff(const char *name, const char *type)
              :XExpressor(name, type)
{
   // Normal WeightedDiff constructor
   if(kCS) cout << "---XWeightedDiff::XWeightedDiff------" << endl;

   fNDefPar = 1;
}//Constructor

//______________________________________________________________________________
XWeightedDiff::~XWeightedDiff()
{
   // WeightedDiff destructor
   if(kCS) cout << "---XWeightedDiff::~XWeightedDiff------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XWeightedDiff::CreateArray(Int_t length)
{
   // Initialize array: subtract background corrected MM from PM intensity
   if(kCSa) cout << "------XWeightedDiff::CreateArray------" << endl;

   if (fArray) DeleteArray();
   if (!(fArray = new (nothrow) Double_t[length])) return errInitMemory;
   fLength = length;

   for (Int_t i=0; i<fLength; i++) {
      fArray[i] = fInten1[i] - fInten2[i]; 
   }//for_i

   return errNoErr;
}//CreateArray

//______________________________________________________________________________
Int_t XWeightedDiff::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate weighted difference
   // value1 = mean; value2 = variance; num = length of trimmed array
   if(kCSa) cout << "------XWeightedDiff::Calculate------" << endl;

// Test number of parameters
   if (TestNumParameters(1) != errNoErr) return errInitParameters;

// Get parameters
   Double_t wmaxinten = fPars[0];
   Double_t maxinten  = fTreeInfo->GetValue("fMaxInten");
   Double_t maxpix    = fTreeInfo->GetValue("fMaxNPixels");

   Double_t w1, w2, w3, value;
   Double_t *weight = 0;
   if (!(weight = new (nothrow) Double_t[fLength])) return errInitMemory;
   for (Int_t i=0; i<fLength; i++) weight[i] = 1.0; 
   
// Calculate weighted mean
   Double_t mean      = 0;
   Double_t sumweight = 0;
   Int_t    len       = 0;
   for (Int_t i=0; i<fLength; i++) {
   // weight for ratio MM/PM
      if (fInten2[i] > fInten1[i]) {
         weight[i] = 0;
         continue;
      } else if ((fInten2[i] >= maxinten) && (fInten1[i] >= maxinten)) {
         weight[i] = wmaxinten;
         sumweight += weight[i];
         mean += fInten1[i] * weight[i];
//??         mean += (fInten1[i] - fBg1[i]) * weight[i];
         len++;
         continue;
      } else {
         w1 = 1 - fInten2[i]/fInten1[i];
      }//if
   // weight for StdDev[%]
      w2 = 1 - (fStdev1[i]/fInten1[i]) * (fStdev2[i]/fInten2[i]);      
   // weight for ratio NumPixels/MaxNumPixels
      w3 = (fNPix1[i]/maxpix) * (fNPix2[i]/maxpix);
   // overall weight
      weight[i] = w1 * w2 * w3; 
   
      sumweight += weight[i];
      mean += fArray[i] * weight[i];
      len++;
   }//for_i

// Calculate variance
   Double_t var = 0;
   if (sumweight > 0) {
      mean /= sumweight;
   // variance: Netz; Formeln der Mathematik; p.562
//      var = 0;
      for (Int_t i=0; i<fLength; i++) {
         value = fArray[i] - mean;
         var += value * value * weight[i];
      }//for_i 

      if (len > 1) {
         var /= (len - 1) * sumweight;
      } else {
         var = 0; //var = SP[?]*SP[?];
      }//if
   } else {
      mean = -1;
      var  =  0; //std = sqrt(vVar)
   }//if

   if (weight) {delete [] weight; weight = 0;}

// Return values
   value1 = mean;
   value2 = var;
   num    = len;
   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAvgDif                                                              //
//                                                                      //
// Average Difference expression algorithm (MAS4)                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAvgDif::XAvgDif()
        :XWeightedDiff()
{
   // Default AvgDif constructor
   if(kCS) cout << "---XAvgDif::XAvgDif(default)------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XAvgDif::XAvgDif(const char *name, const char *type)
        :XWeightedDiff(name, type)
{
   // Normal AvgDif constructor
   if(kCS) cout << "---XAvgDif::XAvgDif------" << endl;

   fNDefPar = 0;
}//Constructor

//______________________________________________________________________________
XAvgDif::~XAvgDif()
{
   // AvgDif destructor
   if(kCS) cout << "---XAvgDif::~XAvgDif------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XAvgDif::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate average difference as described in MAS4
   // value1 = mean; value2 = variance; num = length of array within STP
   if(kCSa) cout << "------XAvgDif::Calculate------" << endl;

   Double_t STP = (fNPar > 0) ? fPars[0] : 3.0;  //default=3
   Int_t    len = fLength;

// Calculate mean difference
   Double_t mean = 0;
   for (Int_t i=0; i<fLength; i++) {
      mean += fArray[i];
   }//for_i
   mean /= fLength;

// Calculate variance
   Double_t value;
   Double_t var = 0;
   for (Int_t i=0; i<fLength; i++) {
      value = fArray[i] - mean;
      var += value * value;
   }//for_i
//   var /= fLength;
//TEST WITH: !!! compare to MAS4!
   var /= (fLength - 1); // degrees of freedom? test for null division

// Super scoring if more than eight probe pairs
   if (fLength > 8) {
      Double_t range   = STP * sqrt(var);
      Double_t rangeHi = mean + range;
      Double_t rangeLo = mean - range;

   // mean difference
      mean = 0;
      len  = 0;
      for (Int_t i=0; i<fLength; i++) {
         if ((fArray[i] <= rangeHi) || (fArray[i] >= rangeLo)) {
//??         if ((fArray[i] < rangeHi) || (fArray[i] > rangeLo)) {
            mean += fArray[i];
            len++;
         } else {
            fArray[i] = -999999;
         }//if
      }//for_i
      if (len > 1) {
         mean /= len;
      }//if

   // variance
      var = 0;
      for (Int_t i=0; i<fLength; i++) {
         if (fArray[i] > -999999) {
            value = fArray[i] - mean;
            var += value * value;
         }//if
      }//for_i
      if (len > 1) {
         var /= len;
//??         var /= len - 1; // degrees of freedom?
      }//if
   }//if

// Return values
   value1 = mean;
   value2 = var;
   num    = len;
   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTukeyBiweight                                                       //
//                                                                      //
// Signal value expression based on tukey biweight algorithm (MAS5)     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTukeyBiweight::XTukeyBiweight()
               :XExpressor()
{
   // Default TukeyBiweight constructor
   if(kCS) cout << "---XTukeyBiweight::XTukeyBiweight(default)------" << endl;

   fNDefPar = 6;
}//Constructor

//______________________________________________________________________________
XTukeyBiweight::XTukeyBiweight(const char *name, const char *type)
               :XExpressor(name, type)
{
   // Normal TukeyBiweight constructor
   if(kCS) cout << "---XTukeyBiweight::XTukeyBiweight------" << endl;

   fNDefPar = 6;
}//Constructor

//______________________________________________________________________________
XTukeyBiweight::~XTukeyBiweight()
{
   // TukeyBiweight destructor
   if(kCS) cout << "---XTukeyBiweight::~XTukeyBiweight------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XTukeyBiweight::CreateArray(Int_t length)
{
   // Initialize array: subtract background corrected MM from PM intensity
   if(kCSa) cout << "------XTukeyBiweight::CreateArray------" << endl;

// Test number of parameters
   if (TestNumParameters(6) != errNoErr) return errInitParameters;

// Get parameters
   Double_t neglog = fPars[5];  //value used if value is negative

   if (fArray) DeleteArray();
   if (!(fArray = new (nothrow) Double_t[length])) return errInitMemory;
   fLength = length;

   Double_t *pm = 0;
   Double_t *mm = 0;
   if (!(pm = new (nothrow) Double_t[fLength])) {return errInitMemory;}
   if (!(mm = new (nothrow) Double_t[fLength])) {delete [] pm; return errInitMemory;}

   // PM and MM
   for (Int_t i=0; i<length; i++) {
      pm[i] = fInten1[i]; 
      mm[i] = fInten2[i];
   }//for_i

   // convert PM and MM to log
   pm = Array2Log(length, pm, neglog, fLogBase);
   mm = Array2Log(length, mm, neglog, fLogBase);

   // subtract PM - MM
   for (Int_t i=0; i<length; i++) {
      fArray[i] = pm[i] - mm[i]; 
   }//for_i

   if (mm) {delete [] mm; mm = 0;}
   if (pm) {delete [] pm; pm = 0;}

   return errNoErr;
}//CreateArray

//______________________________________________________________________________
Int_t XTukeyBiweight::Calculate(Double_t &value1, Double_t &value2, Int_t &num)
{
   // Calculate signal value as described in MAS5
   // value1 = mean; value2 = variance; num = length of array within STP
   if(kCSa) cout << "------XTukeyBiweight::Calculate------" << endl;

// Test number of parameters
   if (TestNumParameters(6) != errNoErr) return errInitParameters;

// Parameters
   Double_t ctau   = fPars[0];  //contrast tau
   Double_t stau   = fPars[1];  //scale tau
   Double_t delta  = fPars[2];  //delta
   Double_t c      = fPars[3];  //c for tukey biweight
   Double_t eps    = fPars[4];  //epsilon for tukey biweight
   Double_t neglog = fPars[5];  //value used if value is negative
   Double_t var    = 0.0;

// Tukey biweight
   Double_t sb = TStat::TukeyBiweight(fLength, fArray, var, c, eps);

   Double_t *ct = 0;
   // ct must be at least have length two, if fLength=1
   if (!(ct = new (nothrow) Double_t[fLength + 1])) return errInitMemory;

// Compute contrast value and subtract from PM
   for (Int_t i=0; i<fLength; i++) {
      ct[i] = fInten2[i];

      if (fInten2[i] >= fInten1[i]) {
         if (sb > ctau) {
            ct[i] = fInten1[i] / TMath::Power(2, sb);
         } else if (sb <= ctau) {
            ct[i] = fInten1[i] / TMath::Power(2, ctau/(1 + (ctau - sb)/stau));
         }//if
      }//if

      fArray[i] = TMath::Max((fInten1[i] - ct[i]), delta); 
   }//for_i

// Compute expression value and variance
   fArray = Array2Log(fLength, fArray, neglog, fLogBase);
   ct[0]  = TStat::TukeyBiweight(fLength, fArray, ct[1], c, eps);
//Problem: convert var or sqrt(var) from log2 back to pow(2,?) ??
   ct     = Array2Pow(2, ct, fLogBase);

// Return values
   value1 = ct[0];
   value2 = ct[1];
   num    = fLength;

   if (ct) {delete [] ct; ct = 0;}

   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMedianPolish                                                        //
//                                                                      //
// Median polish expression algorithm used for RMA                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMedianPolish::XMedianPolish()
              :XMultichipExpressor()
{
   // Default MedianPolish constructor
   if(kCS) cout << "---XMedianPolish::XMedianPolish(default)------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XMedianPolish::XMedianPolish(const char *name, const char *type)
              :XMultichipExpressor(name, type)
{
   // Normal MedianPolish constructor
   if(kCS) cout << "---XMedianPolish::XMedianPolish------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XMedianPolish::~XMedianPolish()
{
   // MedianPolish destructor
   if(kCS) cout << "---XMedianPolish::~XMedianPolish------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XMedianPolish::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Calculate median polish
   if(kCSa) cout << "------XMedianPolish::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t nrow = (Int_t)(fLength / n);
   Int_t ncol = n;

   if (fResiduals) {delete [] fResiduals; fResiduals = 0;}
   if (!(fResiduals = new (nothrow) Double_t[fLength])) return errInitMemory;

   Double_t *rowmed = 0;
   if (!(rowmed = new (nothrow) Double_t[nrow])) return errInitMemory;

   err = DoMedianPolish(nrow, ncol, fArray, x, rowmed, y, fResiduals, 0);

   if (rowmed) {delete [] rowmed; rowmed = 0;}
   
   return err;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFARMS                                                               //
//                                                                      //
// Factor Analysis for Robust Microarray Summarization (Hochreiter S)   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XFARMS::XFARMS()
       :XMultichipExpressor()
{
   // Default FARMS constructor
   if(kCS) cout << "---XFARMS::XFARMS(default)------" << endl;

   fNDefPar = 5;
}//Constructor

//______________________________________________________________________________
XFARMS::XFARMS(const char *name, const char *type)
       :XMultichipExpressor(name, type)
{
   // Normal FARMS constructor
   if(kCS) cout << "---XFARMS::XFARMS------" << endl;

   fNDefPar = 5;
}//Constructor

//______________________________________________________________________________
XFARMS::~XFARMS()
{
   // FARMS destructor
   if(kCS) cout << "---XFARMS::~XFARMS------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XFARMS::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Calculate FARMS (Factor Analysis for Robust Microarray Summarization)
   // Reference: Hochreiter S, et al., Bioinformatics 2006 Apr 15;22(8):943-9.
   if(kCSa) cout << "------XFARMS::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t nrow = (Int_t)(fLength / n);
   Int_t ncol = n; 

// Get parameters
   Int_t    version  = (Int_t)(fPars[0]);  // version of package farms 
   Double_t weight   = fPars[1];           // hyperparameter
   Double_t mu       = fPars[2];           // hyperparameter
   Double_t scale    = fPars[3];           // scaling parameter
   Double_t tol      = fPars[4];           // termination tolerance of EM
   Int_t    cyc      = (Int_t)(fPars[5]);  // maximum number of cycles of EM
   Int_t    weighted = (version == 131) ? (Int_t)(fPars[6]) : 0;  // weighted mean

   if (cyc <= 0) cyc = 2*ncol;

   if (version == 131) {
      err = this->DoFARMS131(nrow, ncol, fArray, x, y, weight, mu, scale,
                                  tol, cyc, weighted);
   } else if (version == 130) {
      err = this->DoFARMS130(nrow, ncol, fArray, x, y, weight, mu, scale, tol, cyc);
   } else {
      cerr << "Error: Version <" << version << "> is not supported." << endl;
      err = errAbort;
   }//if
   
   return err;
}//Calculate

//______________________________________________________________________________
Int_t XFARMS::DoFARMS130(Int_t nrow, Int_t ncol, Double_t *inten, Double_t *x, 
              Double_t *y, Double_t weight, Double_t mu, Double_t scale, 
              Double_t tol, Double_t cyc)
{
   // Calculate FARMS (Factor Analysis for Robust Microarray Summarization)
   // Reference: Hochreiter S, et al., Bioinformatics 2006 Apr 15;22(8):943-9.
   // Note: Adapted from R package farms_1.3: file farms.R (license GNU GPL >= 2)
   if(kCSa) cout << "------XFARMS::DoFARMS130------" << endl;

   Int_t    err    = errNoErr;
   Double_t sum    = 0.0;
   Double_t sumL   = 0.0;
   Double_t alpha  = 0.0;
   Double_t bbeta  = 0.0;
   Double_t mean   = 0.0;
   Double_t a      = 0.0;
   Double_t ezz    = 0.0;
   Double_t lambda = 0.0;

// Initialize local arrays
   Double_t *xmean  = 0;  // column means
   Double_t *X      = 0;  // matrix X = x - xmean
   Double_t *mm     = 0;  // temporary matrix
   Double_t *XX     = 0;  // matrix (positive definit)
   Double_t *diagXX = 0;  // diagonale(XX)
   Double_t *L      = 0;  // factor loadings
   Double_t *Lold   = 0;  // temporary L
   Double_t *Ph     = 0;  // diagonal uniqueness matrix
   Double_t *PsiL   = 0;  // 
   Double_t *beta   = 0;  // 
   Double_t *XXbeta = 0;  // 
   Double_t *TXbeta = 0;  // 
   Double_t *EZZ    = 0;  //
   Double_t *c      = 0;  //

// Initialize memory for local arrays
   if (!(X  = new (nothrow) Double_t[nrow*ncol])) {err = errInitMemory; goto cleanup;}
   if (!(mm = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}
   if (!(XX = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}

   if (!(xmean  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(diagXX = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(L      = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Lold   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Ph     = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(PsiL   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(beta   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(XXbeta = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(TXbeta = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(EZZ    = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(c      = new (nothrow) Double_t[ncol])) {err = errInitMemory; goto cleanup;}

   // row means
   for (Int_t i=0; i<nrow; i++) {
      xmean[i] = 0.0;
      for (Int_t j=0; j<ncol; j++) {
         xmean[i] += inten[i*ncol + j];
      }//for_j
      xmean[i] = xmean[i]/ncol;
   }//for_i

   // center data: X = x - xmean
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = inten[i*ncol + j] - xmean[i];
      }//for_i
   }//for_j

   // mm = t(X)*X/ncol
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            sum += X[i*ncol + j]*X[k*ncol + j];
         }//for_j
         mm[i*nrow + k] = sum/ncol;
      }//for_k
   }//for_i

   // XX = (mm + t(mm))/2 but positive definit
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = (mm[i*nrow + k] + mm[k*nrow + i])/2.0;
         XX[i*nrow + k] = (sum > 0.0) ? sum : 0.0;
      }//for_k
   }//for_i

   // diagonale of XX
   mean = 0.0;
   for (Int_t i=0; i<nrow; i++) {
      diagXX[i] = XX[i*nrow + i];
      mean     += diagXX[i];
   }//for_i
   mean  = mean/nrow;
   alpha = weight/mean;
   bbeta = mu*alpha;

   // factor loadings: L
   for (Int_t i=0; i<nrow; i++) {
      L[i]    = 0.8*diagXX[i];
      Lold[i] = 0.0;
   }//for_i

   // diagonal uniqueness matrix: Ph
   for (Int_t i=0; i<nrow; i++) {
      Ph[i] = 0.05*diagXX[i];
   }//for_i

//  Expectation-Maximization (EM) algorithm
   while (cyc >0) {
   // E step
      a = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         PsiL[i] = L[i]/Ph[i];
         a      += L[i]*PsiL[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         beta[i] = PsiL[i]/a;
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         sum = 0.0;
         for (Int_t k=0; k<nrow; k++) {
            sum += XX[i*nrow + k]*beta[k];
         }//for_k
         XXbeta[i] = sum;
      }//for_i

      sum = 0.0;
      ezz = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         sum += beta[i]*L[i];
         ezz += beta[i]*XXbeta[i];
      }//for_i
      ezz = ezz - sum;

      for (Int_t i=0; i<nrow; i++) {
         TXbeta[i] = XXbeta[i] + bbeta*Ph[i];
         EZZ[i]    = ezz       + alpha*Ph[i];
      }//for_i

   // M step
      for (Int_t i=0; i<nrow; i++) {
         L[i] = TXbeta[i]/EZZ[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         Ph[i] = diagXX[i] - XXbeta[i]*L[i] +  alpha*Ph[i]*L[i]*(bbeta - L[i]);
      }//for_i

      sumL = sum = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         sumL += L[i]*L[i];
         sum  += Lold[i]*Lold[i];
      }//for_i
      sumL = TMath::Abs(TMath::Sqrt(sumL) - TMath::Sqrt(sum));

      if (sumL < tol) break;

      for (Int_t i=0; i<nrow; i++) Lold[i] = L[i];
      cyc--;
   }//while

   for (Int_t j=0; j<ncol; j++) {
      c[j] = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         c[j] += X[i*ncol + j]*beta[i];
      }//for_i
   }//for_j

   lambda = mean = 0.0;
   for (Int_t i=0; i<nrow; i++) {
      lambda += L[i];
      mean   += xmean[i];
   }//for_i
   lambda = lambda/nrow;
   mean   = mean/nrow;

// Final result x and stdev y
   for (Int_t j=0; j<ncol; j++) {
      x[j] = scale*lambda*c[j] + mean;
      y[j] = 1.0/a;
   }//for_j

// Convert results
   x = Array2Pow(ncol, x, fLogBase);
//?   y = Array2Pow(ncol, y, fLogBase);

// Cleanup
cleanup:
   if (c)      {delete [] c;      c      = 0;}
   if (EZZ)    {delete [] EZZ;    EZZ    = 0;}
   if (TXbeta) {delete [] TXbeta; TXbeta = 0;}
   if (XXbeta) {delete [] XXbeta; XXbeta = 0;}
   if (beta)   {delete [] beta;   beta   = 0;}
   if (PsiL)   {delete [] PsiL;   PsiL   = 0;}
   if (Ph)     {delete [] Ph;     Ph     = 0;}
   if (Lold)   {delete [] Lold;   Lold   = 0;}
   if (L)      {delete [] L;      L      = 0;}
   if (diagXX) {delete [] diagXX; diagXX = 0;}
   if (xmean)  {delete [] xmean;  xmean  = 0;}
   if (XX)     {delete [] XX;     XX     = 0;}
   if (mm)     {delete [] mm;     mm     = 0;}
   if (X)      {delete [] X;      X      = 0;}
   
   return err;
}//DoFARMS130

//______________________________________________________________________________
Int_t XFARMS::DoFARMS131(Int_t nrow, Int_t ncol, Double_t *inten, 
              Double_t *x, Double_t *y, Double_t weight, Double_t mu, 
              Double_t scale, Double_t tol, Double_t cyc, Bool_t weighted)
{
   // Calculate FARMS (Factor Analysis for Robust Microarray Summarization)
   // Reference: Hochreiter S, et al., Bioinformatics 2006 Apr 15;22(8):943-9.
   // Note: Adapted from R package farms_1.3.1: file farms.R (license GNU GPL >= 2)
   if(kCSa) cout << "------XFARMS::DoFARMS131------" << endl;

   Int_t    err    = errNoErr;
   Double_t sum    = 0.0;
   Double_t var    = 0.0;
   Double_t alpha  = 0.0;
   Double_t bbeta  = 0.0;
   Double_t mean   = 0.0;
   Double_t a      = 0.0;
   Double_t ezz    = 0.0;

// Initialize local arrays
   Double_t *xmean  = 0;  // column means
   Double_t *xstdv  = 0;  // standard deviation
   Double_t *X      = 0;  // matrix X = x - xmean
   Double_t *mm     = 0;  // temporary matrix
   Double_t *XX     = 0;  // matrix (positive definit)
   Double_t *diagXX = 0;  // diagonale(XX)
   Double_t *L      = 0;  // factor loadings
   Double_t *Lold   = 0;  // temporary L
   Double_t *Ph     = 0;  // diagonal uniqueness matrix
   Double_t *PsiL   = 0;  // 
   Double_t *beta   = 0;  // 
   Double_t *arr1   = 0;  // 
   Double_t *arr2   = 0;  // 
   Double_t *EZZ    = 0;  //
   Double_t *c      = 0;  //

// Initialize memory for local arrays
   if (!(X  = new (nothrow) Double_t[nrow*ncol])) {err = errInitMemory; goto cleanup;}
   if (!(mm = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}
   if (!(XX = new (nothrow) Double_t[nrow*nrow])) {err = errInitMemory; goto cleanup;}

   if (!(xmean  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(xstdv  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(diagXX = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(L      = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Lold   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(Ph     = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(PsiL   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(beta   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(arr1   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(arr2   = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(EZZ    = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(c      = new (nothrow) Double_t[ncol])) {err = errInitMemory; goto cleanup;}

// Normalize
   // calculate row means
   for (Int_t i=0; i<nrow; i++) {
      xmean[i] = 0.0;
      for (Int_t j=0; j<ncol; j++) {
         xmean[i] += inten[i*ncol + j];
      }//for_j
      xmean[i] = xmean[i]/ncol;
   }//for_i

   // center data: XX = inten - xmean
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = inten[i*ncol + j] - xmean[i];
      }//for_i
   }//for_j

   // mm = X*t(X)/ncol
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            sum += X[i*ncol + j]*X[k*ncol + j];
         }//for_j
         mm[i*nrow + k] = sum/ncol;
      }//for_k
   }//for_i

   // calculate standard deviation
   for (Int_t i=0; i<nrow; i++) {
      xstdv[i] = TMath::Sqrt(mm[i*nrow + i]);
      if (xstdv[i] == 0) xstdv[i] = 1.0;
   }//for_i

   // standardize X to variance 1: X = inten/stdev
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = inten[i*ncol + j]/xstdv[i];
      }//for_i
   }//for_j

// XX
   // row means
   for (Int_t i=0; i<nrow; i++) {
      xmean[i] = 0.0;
      for (Int_t j=0; j<ncol; j++) {
         xmean[i] += X[i*ncol + j];
      }//for_j
      xmean[i] = xmean[i]/ncol;
   }//for_i

   // center data: X = X - xmean
   for (Int_t j=0; j<ncol; j++) {
      for (Int_t i=0; i<nrow; i++) {
         X[i*ncol + j] = X[i*ncol + j] - xmean[i];
      }//for_i
   }//for_j

   // mm = t(X)*X/ncol
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            sum += X[i*ncol + j]*X[k*ncol + j];
         }//for_j
         mm[i*nrow + k] = sum/ncol;
      }//for_k
   }//for_i

   // XX = (mm + t(mm))/2 but positive definit
   for (Int_t i=0; i<nrow; i++) {
      for (Int_t k=0; k<nrow; k++) {
         sum = (mm[i*nrow + k] + mm[k*nrow + i])/2.0;
         XX[i*nrow + k] = (sum > 0.0) ? sum : 0.0;
      }//for_k
   }//for_i

   // diagonale of XX
   for (Int_t i=0; i<nrow; i++) {
      diagXX[i] = XX[i*nrow + i];
   }//for_i

   alpha = weight*nrow;
   bbeta = mu*alpha;

   // factor loadings: L
   for (Int_t i=0; i<nrow; i++) {
      L[i]    = TMath::Sqrt(0.8*diagXX[i]);
      Lold[i] = 0.0;
   }//for_i

   // diagonal uniqueness matrix: Ph
   for (Int_t i=0; i<nrow; i++) {
      Ph[i] = diagXX[i] - L[i]*L[i];
   }//for_i

//  Expectation-Maximization (EM) algorithm
   while (cyc >0) {
   // E step
      a = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         PsiL[i] = L[i]/Ph[i];
         a      += L[i]*PsiL[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         beta[i] = PsiL[i]/a;
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         sum = 0.0;
         for (Int_t k=0; k<nrow; k++) {
            sum += XX[i*nrow + k]*beta[k];
         }//for_k
         arr1[i] = sum;
      }//for_i

      sum = 0.0;
      ezz = 1.0;
      for (Int_t i=0; i<nrow; i++) {
         sum += beta[i]*L[i];
         ezz += beta[i]*arr1[i];
      }//for_i
      ezz = ezz - sum;

      for (Int_t i=0; i<nrow; i++) {
         arr2[i] = arr1[i] + bbeta*Ph[i];
         EZZ[i]  = ezz     + alpha*Ph[i];
      }//for_i

   // M step
      for (Int_t i=0; i<nrow; i++) {
         L[i] = arr2[i]/EZZ[i];
      }//for_i

      for (Int_t i=0; i<nrow; i++) {
         Ph[i] = diagXX[i] - arr1[i]*L[i] +  alpha*Ph[i]*L[i]*(bbeta - L[i]);
         if (Ph[i] == 0) Ph[i] = 1.0;
      }//for_i

      sum = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         Double_t LL = Lold[i] - L[i];
         sum  += TMath::Sqrt(LL*LL);
      }//for_i
      if (sum < tol) break;

      for (Int_t i=0; i<nrow; i++) Lold[i] = L[i];
      cyc--;
   }//while

   sum = 0.0;
   for (Int_t j=0; j<ncol; j++) {
      c[j] = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         c[j] += X[i*ncol + j]*beta[i];
      }//for_i
      sum  += c[j];
   }//for_j
   mean = sum/ncol;

   var = TStat::StDev(ncol, c, mean);
   if (var == 0) var = 1.0;

   for (Int_t j=0; j<ncol; j++) {
      c[j] = c[j]/var;
   }//for_j

   for (Int_t i=0; i<nrow; i++) {
      L[i] = L[i]*var;
   }//for_i

   a = 1.0;
   for (Int_t i=0; i<nrow; i++) {
      PsiL[i] = L[i]/Ph[i];
      a      += L[i]*PsiL[i];
   }//for_i

   for (Int_t i=0; i<nrow; i++) {
      arr1[i] = L[i]*xstdv[i];
   }//for_i

   for (Int_t i=0; i<nrow; i++) {
      arr2[i] = xmean[i]*xstdv[i];
   }//for_i

   if (weighted) {
      sum = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         PsiL[i] = L[i]*L[i]/Ph[i];
         sum    += PsiL[i];
      }//for_i
      for (Int_t i=0; i<nrow; i++) {
         PsiL[i] = PsiL[i]/sum;
      }//for_i

      var = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         var +=arr1[i]* PsiL[i];
      }//for_i
   } else {
      var  = TStat::Median(nrow, arr1);
   }//if

   mean = TStat::Mean(nrow, arr2);

// Final expression and stdev
   for (Int_t j=0; j<ncol; j++) {
      x[j] = scale*var*c[j] + mean;
      y[j] = 1.0/a;
   }//for_j

// Convert results
   x = Array2Pow(ncol, x, fLogBase);
//?   y = Array2Pow(ncol, y, fLogBase);

// Cleanup
cleanup:
   if (c)      {delete [] c;      c      = 0;}
   if (EZZ)    {delete [] EZZ;    EZZ    = 0;}
   if (arr2)   {delete [] arr2;   arr2   = 0;}
   if (arr1)   {delete [] arr1;   arr1   = 0;}
   if (beta)   {delete [] beta;   beta   = 0;}
   if (PsiL)   {delete [] PsiL;   PsiL   = 0;}
   if (Ph)     {delete [] Ph;     Ph     = 0;}
   if (Lold)   {delete [] Lold;   Lold   = 0;}
   if (L)      {delete [] L;      L      = 0;}
   if (diagXX) {delete [] diagXX; diagXX = 0;}
   if (xstdv)  {delete [] xstdv;  xstdv  = 0;}
   if (xmean)  {delete [] xmean;  xmean  = 0;}
   if (XX)     {delete [] XX;     XX     = 0;}
   if (mm)     {delete [] mm;     mm     = 0;}
   if (X)      {delete [] X;      X      = 0;}
   
   return err;
}//DoFARMS131


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDFW                                                                 //
//                                                                      //
// Distribution Free Weighted Fold Change (Chen Z)                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDFW::XDFW()
     :XMultichipExpressor()
{
   // Default DFW constructor
   if(kCS) cout << "---XDFW::XDFW(default)------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XDFW::XDFW(const char *name, const char *type)
     :XMultichipExpressor(name, type)
{
   // Normal DFW constructor
   if(kCS) cout << "---XDFW::XDFW------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XDFW::~XDFW()
{
   // DFW destructor
   if(kCS) cout << "---XDFW::~XDFW------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XDFW::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Calculate DFW (Distribution Free Weighted Fold Change)
   // Reference: Chen Z, et al., Bioinformatics 2007 Feb 1;23(3):321-7.
   if(kCSa) cout << "------XDFW::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t nrow = (Int_t)(fLength / n);
   Int_t ncol = n; 

   Double_t max, min, WR, WSD;

// Get parameters or use default values
   Int_t    mm = (fNPar > 0) ? (Int_t)fPars[0] : 3;
   Int_t    nn = (fNPar > 1) ? (Int_t)fPars[1] : 1;
   Double_t cc = (fNPar > 2) ? fPars[2] : 0.01;

// Initialize local arrays
   Double_t *range  = 0;
   Double_t *stdev  = 0;
   Double_t *weight = 0;
   Double_t *arr    = 0;

// Initialize memory for local arrays
   if (!(range  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(stdev  = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(weight = new (nothrow) Double_t[nrow])) {err = errInitMemory; goto cleanup;}
   if (!(arr    = new (nothrow) Double_t[ncol])) {err = errInitMemory; goto cleanup;}

// Compute range and stdvn for each row
   for (Int_t i=0; i<nrow; i++) {
      if (ncol == 1) {
         range[i] = fArray[i];
         stdev[i] = 1.0;
      } else {
         Double_t mean = 0.0;
         for (Int_t j=0; j<ncol; j++) {
            arr[j] = fArray[i*ncol + j];
            mean  += arr[j];
         }//for_j
         mean = mean/ncol;

         range[i] = TStat::Max(ncol, arr) - TStat::Min(ncol, arr);
         stdev[i] = TStat::StDev(ncol, arr, mean);
      }//if
   }//for_i

// Weighted range (WR)
   weight = this->Weight(nrow, range, weight);
   for (Int_t j=0; j<ncol; j++) {
      arr[j] = 0.0;
      for (Int_t i=0; i<nrow; i++) {
         arr[j] += fArray[i*ncol + j] * weight[i];
      }//for_i
   }//for_j

   max = TStat::Max(ncol, arr);
   min = TStat::Min(ncol, arr);
   WR  = max - min;

// Final stdev y
   if (WR == 0.0) {
      for (Int_t j=0; j<ncol; j++) {
         y[j] = 0.0;
      }//for_j
   } else {
      for (Int_t j=0; j<ncol; j++) {
         y[j] = (arr[j] - min)/WR;
      }//for_j
   }//if 

// Weighted stdev (WSD)
   weight = this->Weight(nrow, stdev, weight);
   WSD = 0.0;
   for (Int_t i=0; i<nrow; i++) {
      WSD  += stdev[i] * weight[i];
   }//for_i

// Final result x
   for (Int_t j=0; j<ncol; j++) {
      x[j] = min + cc * y[j] * TMath::Power(WR, mm) * TMath::Power(WSD, nn);
   }//for_j

// Convert results
   x = Array2Pow(ncol, x, fLogBase);
   y = Array2Pow(ncol, y, fLogBase);

// Cleanup
cleanup:
   if (arr)    {delete [] arr;    arr    = 0;}
   if (weight) {delete [] weight; weight = 0;}
   if (stdev)  {delete [] stdev;  stdev  = 0;}
   if (range)  {delete [] range;  range  = 0;}
   
   return err;
}//Calculate

//______________________________________________________________________________
Double_t *XDFW::Weight(Int_t n, const Double_t *arr, Double_t *w)
{
   // Weighting function using Tukey biweight.
   if(kCSa) cout << "------XDFW::Weight------" << endl;

   Double_t *x   = new Double_t[n];
   Double_t  med = TStat::Median(n, arr);
   Double_t  max = 0.0;
   Double_t  sum = 0.0;

   for (Int_t i=0; i<n; i++) {
      x[i] = arr[i] - med;
      if  (TMath::Abs(x[i]) > max) max =  TMath::Abs(x[i]); 
   }//for_i

   if (max == 0.0) {
      for (Int_t i=0; i<n; i++) {
         w[i] = arr[i];
         sum += w[i];
      }//for_i
   } else {
      // Tukey biweight
      Double_t uu  = 0.0;
      Double_t uu2 = 0.0;
      for (Int_t i=0; i<n; i++) {
         uu   = x[i]/max;
         uu2  = (1.0 - uu*uu);
         w[i] = uu2*uu2;
         sum += w[i];
      }//for_i
   }//if

   // final weight
   sum = (sum == 0.0) ? n : sum;  //prevent zero division
   for (Int_t i=0; i<n; i++) {
      w[i] = w[i]/sum;
   }//for_i

   delete [] x;

   return w;
}//Weight


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFIRMA                                                               //
//                                                                      //
// Finding Isoforms using Robust Multichip Analysis (Purdom E)          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XFIRMA::XFIRMA()
       :XSpliceExpressor()
{
   // Default FIRMA constructor
   if(kCS) cout << "---XFIRMA::XFIRMA(default)------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XFIRMA::XFIRMA(const char *name, const char *type)
       :XSpliceExpressor(name, type)
{
   // Normal FIRMA constructor
   if(kCS) cout << "---XFIRMA::XFIRMA------" << endl;

   fNDefPar = 3;
}//Constructor

//______________________________________________________________________________
XFIRMA::~XFIRMA()
{
   // XFIRMA destructor
   if(kCS) cout << "---XFIRMA::~XFIRMA------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XFIRMA::Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z, 
              Double_t *dx, Double_t *dy, Int_t *idx, Int_t nsub)
{
   // Calculate summarization and alternative splicing
   // Reference: Purdom E, et al., Bioinformatics 2008 Aug 1;24(15):1707-14.
   // Source code is based on R function firma() extracted from aroma.affymetrix, see:
   // http://www.aroma-project.org/node/81
   if(kCSa) cout << "------XFIRMA::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t nrow = (Int_t)(fLength / n);
   Int_t ncol = n;
   Int_t cnt  = 0;

   Double_t se = 0.0;

// Initialize residuals array
   if (fResiduals) {delete [] fResiduals; fResiduals = 0;}
   if (!(fResiduals = new (nothrow) Double_t[fLength])) return errInitMemory;

// Initialize memory for arrays
   Double_t *inten   = 0;
   Double_t *score   = 0;
   Double_t *rowmed  = 0;
   Double_t *trLevel = 0;
   Double_t *trStdev = 0;
   Double_t *psLevel = 0;
   Double_t *psStdev = 0;
   Double_t *psScore = 0;
   if (!(inten   = new (nothrow) Double_t[fLength])) {err = errInitMemory; goto cleanup;}
   if (!(score   = new (nothrow) Double_t[fLength])) {err = errInitMemory; goto cleanup;}
   if (!(rowmed  = new (nothrow) Double_t[nrow]))    {err = errInitMemory; goto cleanup;}
   if (!(trLevel = new (nothrow) Double_t[ncol]))    {err = errInitMemory; goto cleanup;}
   if (!(trStdev = new (nothrow) Double_t[ncol]))    {err = errInitMemory; goto cleanup;}
   if (!(psLevel = new (nothrow) Double_t[ncol]))    {err = errInitMemory; goto cleanup;}
   if (!(psStdev = new (nothrow) Double_t[ncol]))    {err = errInitMemory; goto cleanup;}
   if (!(psScore = new (nothrow) Double_t[ncol]))    {err = errInitMemory; goto cleanup;}

// Calculate summarization for transcript intensities
   err = DoMedianPolish(nrow, ncol, fArray, trLevel, rowmed, trStdev, score, 0);
   if (err != errNoErr) goto cleanup;

   // affinities sum to zero (on log scale)
   rowmed[nrow-1] = 0.0;
   for (Int_t i=0; i<nrow-1; i++) rowmed[nrow-1] -= rowmed[i];
   rowmed = Array2Pow(nrow, rowmed, fLogBase);

// Compute residuals
   cnt = 0;
   if (strcmp(fLogBase, "0") == 0) {
      for (Int_t i=0; i<nrow; i++) {
         for (Int_t k=0; k<ncol; k++) {
            fResiduals[cnt] = fArray[cnt]/(rowmed[i]*trLevel[k]);
            cnt++;
         }//for_k
      }//for_i
   } else if (strcmp(fLogBase, "log2") == 0) {
      for (Int_t i=0; i<nrow; i++) {
         for (Int_t k=0; k<ncol; k++) {
            fResiduals[cnt] = TMath::Power(2, fArray[cnt])/(rowmed[i]*trLevel[k]);
            fResiduals[cnt] = TMath::Log2(fResiduals[cnt]);
            cnt++;
         }//for_k
      }//for_i
   } else if (strcmp(fLogBase, "log10") == 0) {
      for (Int_t i=0; i<nrow; i++) {
         for (Int_t k=0; k<ncol; k++) {
            fResiduals[cnt] = TMath::Power(10, fArray[cnt])/(rowmed[i]*trLevel[k]);
            fResiduals[cnt] = TMath::Log10(fResiduals[cnt]);
            cnt++;
         }//for_k
      }//for_i
   } else if (strcmp(fLogBase, "log") == 0) {
      for (Int_t i=0; i<nrow; i++) {
         for (Int_t k=0; k<ncol; k++) {
            fResiduals[cnt] = TMath::Power(TMath::E(), fArray[cnt])/(rowmed[i]*trLevel[k]);
            fResiduals[cnt] = TMath::Log(fResiduals[cnt]);
            cnt++;
         }//for_k
      }//for_i
   } else {
      cerr << "Error: LogBase <" << fLogBase << "> is not known." << endl;
      err = errInitParameters; goto cleanup;
   }//if

// Estimate standard error of residuals (center=0)
   se = TStat::MAD(fLength, fResiduals, 0.0, 1.4826, kFALSE, kFALSE);

// Calculate summarization and score for probeset/exon intensities
   cnt = 0;
   for (Int_t i=0; i<nsub; i++) {
      for (Int_t k=0; k<ncol; k++) {
         for (Int_t j=0; j<idx[i]; j++) {
            inten[j*ncol+k] = fArray[(j+cnt)*ncol + k];
            score[j]        = fResiduals[(j+cnt)*ncol + k]/se;
         }//for_j

         psScore[k] = TMath::Median(idx[i], score);
      }//for_k
      cnt += idx[i];

      err = DoMedianPolish(idx[i], ncol, inten, psLevel, rowmed, psStdev, score, 0);
      if (err != errNoErr) goto cleanup;

      for (Int_t k=0; k<ncol; k++) {
         x[i*ncol+k]  = trLevel[k];
         y[i*ncol+k]  = psLevel[k];
         z[i*ncol+k]  = psScore[k];
         dx[i*ncol+k] = trStdev[k];
         dy[i*ncol+k] = psStdev[k];
      }//for_k
   }//for_i

cleanup:
   if (psScore) {delete [] psScore; psScore = 0;}
   if (psStdev) {delete [] psStdev; psStdev = 0;}
   if (psLevel) {delete [] psLevel; psLevel = 0;}
   if (trStdev) {delete [] trStdev; trStdev = 0;}
   if (trLevel) {delete [] trLevel; trLevel = 0;}
   if (rowmed)  {delete [] rowmed;  rowmed  = 0;}
   if (score)   {delete [] score;   score   = 0;}
   if (inten)   {delete [] inten;   inten   = 0;}

   return err;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRMAQualifier                                                        //
//                                                                      //
// Class for GeneChip RMA quality control algorithms                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRMAQualifier::XRMAQualifier()
              :XQualifier()
{
   // Default RMAQualifier constructor
   if(kCS) cout << "---XRMAQualifier::XRMAQualifier(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XRMAQualifier::XRMAQualifier(const char *name, const char *type)
              :XQualifier(name, type)
{
   // Normal RMAQualifier constructor
   if(kCS) cout << "---XRMAQualifier::XRMAQualifier------" << endl;

}//Constructor

//______________________________________________________________________________
XRMAQualifier::~XRMAQualifier()
{
   // RMAQualifier destructor
   if(kCS) cout << "---XRMAQualifier::~XRMAQualifier------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPLMQualifier                                                        //
//                                                                      //
// Class for GeneChip PLM quality control algorithms                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPLMQualifier::XPLMQualifier()
              :XQualifier()
{
   // Default PLMQualifier constructor
   if(kCS) cout << "---XPLMQualifier::XPLMQualifier(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XPLMQualifier::XPLMQualifier(const char *name, const char *type)
              :XQualifier(name, type)
{
   // Normal PLMQualifier constructor
   if(kCS) cout << "---XPLMQualifier::XPLMQualifier------" << endl;

}//Constructor

//______________________________________________________________________________
XPLMQualifier::~XPLMQualifier()
{
   // PLMQualifier destructor
   if(kCS) cout << "---XPLMQualifier::~XPLMQualifier------" << endl;

}//Destructor

