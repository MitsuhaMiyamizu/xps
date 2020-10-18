// Author: Christian Stratowa 11/25/2002             last modified: 01/14/2011

/*
 *******************************************************************************
 ***********************  Statistics Package for ROOT  *************************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2011 Dr. Christian Stratowa
 *
 *  Written by: Christian Stratowa, Vienna, Austria <cstrato@aon.at>
 *
 *******************************************************************************
 * Algorithms for statistical functions partially adapted from:                *
 * R: A Computer Language for Statistical Data Analysis                        *
 * Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka                    *
 * Copyright (C) 1999-2001 The R Development Core Team                         *
 * R is free software, for licensing see the GNU General Public License        *
 * http:/cran.r-project.org/                                                   *
 *******************************************************************************
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

#ifndef __TStat__
#define __TStat__

#include <Riostream.h>
#include <cmath>

// ROOT includes
#include "TMath.h"
#include "TString.h"


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStat                                                                //
//                                                                      //
// Class containing statistical methods                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TStat {

   private:
      
   public :
      TStat() {}
      virtual ~TStat() {}

      static Bool_t    Ident(Int_t n, const Double_t *arr);
      static Double_t  Max(Int_t n, const Double_t *arr);
      static Double_t  Min(Int_t n, const Double_t *arr);
      static Double_t *PMax(Int_t n, const Double_t *arr1, const Double_t *arr2,
                          Double_t *max);
      static Double_t *PMax(Int_t n, Int_t m, const Double_t **table,
                          Double_t *max);
      static Double_t *PMin(Int_t n, const Double_t *arr1, const Double_t *arr2,
                          Double_t *min);
      static Double_t *PMin(Int_t n, Int_t m, const Double_t **table, 
                          Double_t *min);
      static Double_t *CumMax(Int_t n, const Double_t *arr, Double_t *max);
      static Double_t *CumMin(Int_t n, const Double_t *arr, Double_t *min);
      static Double_t  CumProd(Int_t n, const Double_t *arr);
      static Double_t  CumSum(Int_t n, const Double_t *arr);
      static Double_t  LnFact(Int_t n);
      static Double_t  BinomCoeff(Int_t n, Int_t k);
      static void      Rank(Int_t n, Double_t *arr, Int_t *index, Int_t *rank,
                          Bool_t down = kTRUE);
      static void      Rank(Int_t n, Double_t *arr, Double_t *rank);
      static Int_t    *TrueRank(Int_t n, Double_t *arr, Int_t *rank,
                          Bool_t down = kTRUE);
      static Double_t *TrueRank1(Int_t n, Double_t *arr, Int_t *index,
                          Double_t *rank);

      static Double_t GeoMean(Int_t n, const Double_t *arr);
      static Double_t GeoMean(Int_t n, const Double_t *arr, Int_t &len,
                         const Double_t na);
      static Double_t GeoMean(Int_t n, const Double_t *arr, const Double_t *w);
      static Double_t GeoMean(Int_t n, const Double_t *arr, const Double_t trim,
                         Double_t &var, Int_t &len);
      static Double_t Mean(Int_t n, const Double_t *arr);
      static Double_t Mean(Int_t begin, Int_t end, const Double_t *arr);
      static Double_t Mean(Int_t n, const Double_t *arr, Int_t &len,
                         const Double_t na);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t *w);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t trim);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t trim,
                         Double_t &var, Int_t &len);
      static Double_t Mean(Int_t n, const Double_t *arr, const Double_t *w,
                         const Double_t trim, Double_t &var, Int_t &len);
      static Double_t Median(Int_t n, const Double_t *arr);
      static Double_t Median(Int_t n, const Double_t *arr, const Int_t *index);
      static Double_t Median(Int_t n, const Double_t *arr, Bool_t low, Bool_t high);
      static Double_t Median(Int_t n, const Double_t *arr, UShort_t logbase);
      static Double_t MAD(Int_t n, const Double_t *arr, Float_t constant = 1.4826);
      static Double_t MAD(Int_t n, const Double_t *arr, const Double_t trim,
                         Float_t constant = 1.4826);
      static Double_t MAD(Int_t n, const Double_t *arr, const Double_t center,
                         Double_t constant, Bool_t low, Bool_t high);
      static Double_t StDev(Int_t n, const Double_t *arr, const Double_t mean);
      static Double_t Var(Int_t n, const Double_t *arr, const Double_t mean);
      static Double_t Var(Int_t n, const Double_t *arr, const Double_t mean,
                         Int_t &len, const Double_t na);

      static Double_t MedianPolish(Int_t nrow, Int_t ncol, Double_t *x,
                         Double_t *rowmed, Double_t *colmed, Double_t *residu,
                         Int_t maxiter = 10, Double_t eps = 0.01,
                         Bool_t verbose = kFALSE);
      static Double_t MedianPolishTranspose(Int_t nrow, Int_t ncol, Double_t *x,
                         Double_t *rowmed, Double_t *colmed, Double_t *residu,
                         Int_t maxiter = 10, Double_t eps = 0.01,
                         Bool_t verbose = kFALSE);
      static Double_t TukeyBiweight(Int_t n, const Double_t *x, Double_t &var,
                         Double_t c = 5.0, Double_t eps = 0.0001);

      static Double_t Quantile(Int_t n, const Double_t *arr, const Double_t q);
      static Double_t Quantile(Int_t n, const Double_t *arr, const Int_t *index,
                         const Double_t q);
      static Double_t IQR(Int_t n, const Double_t *arr, const Double_t qlo = 0.25,
                         const Double_t qhi = 0.75);
      static Double_t IQR(Int_t n, const Double_t *arr, const Int_t *index,
                         const Double_t qlo = 0.25, const Double_t qhi = 0.75);

      static Double_t *Quantiles(Int_t n, Double_t *arr, Int_t *index,
                          Int_t nquant, Double_t *q, Double_t *quant);

      static void     NextPerm(Int_t n, Int_t k, Int_t *grp);
      static void     NextPerm(Int_t n, Int_t k, Int_t *grp1, Int_t nk, Int_t *grp2);
      static void     Sample(Int_t n, Int_t k, Int_t *grp);
      static void     Sample(Int_t n, Int_t k, Int_t *grp1, Int_t nk, Int_t *grp2);

      static Double_t CorPearson(Int_t n, Double_t *x, Double_t *y);
      static Double_t CorSpearman(Int_t n, Double_t *x, Double_t *y);

      static void     MassDist(Int_t nx, Double_t *x, Double_t *w, Double_t xmin,
                         Double_t xmax, Int_t ny, Double_t *y);
      static void     TwiddleFactor4FFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im);
      static void     TwiddleFactor4IFFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im);
      static void     FFT(Int_t n, Double_t *f_re, Double_t *f_im);
      static void     IFFT(Int_t n, Double_t *f_re, Double_t *f_im);
      static Int_t    FFTDensityConvolve(Int_t n, Double_t *x_re, Double_t *y_re);
      static void     Kernelize(Int_t n, Double_t *x, Double_t bw, const char *kernel);
      static Double_t Bandwidth(Int_t n, Double_t *x, Double_t iqr);
      static void     LinearInterpolate(Double_t *xin, Double_t *yin,
                         Int_t nout, Double_t *xout, Double_t *yout);
      static Int_t    Density(Int_t n, Double_t *x, Double_t *w, Int_t nout,
                         Double_t *xout, Double_t *yout, const char *kernel = "gaussian");
      static Double_t MaxDensity(Int_t n, Double_t *x, Double_t *w, Int_t npts = 512,
                         const char *kernel = "gaussian");

      static Double_t PNormApprox(Double_t x);

      ClassDef(TStat,0) //Stat
};

//______________________________________________________________________________
inline Bool_t TStat::Ident(Int_t n, const Double_t *arr)
{
   // Return true if all data in array are equal

   Bool_t ident = kTRUE;
   for (Int_t i=0; i<n-1; i++) {
      if (arr[i] != arr[i+1]) {ident = kFALSE; break;}
   }//for

   return ident;
}//Max

//______________________________________________________________________________
inline Double_t TStat::Max(Int_t n, const Double_t *arr)
{
   // Maximum of array arr.

   return arr[TMath::LocMax(n,arr)];
}//Max

//______________________________________________________________________________
inline Double_t TStat::Min(Int_t n, const Double_t *arr)
{
   // Minimum of array arr.

   return arr[TMath::LocMin(n,arr)];
}//Min

//______________________________________________________________________________
inline Double_t *TStat::PMax(Int_t n, const Double_t *arr1, const Double_t *arr2,
                        Double_t *max)
{
  // Returns maximum of arr1 and arr2 at each index in array max.

   for (Int_t i=0; i<n; i++) max[i] = (arr1[i] > arr2[i]) ? arr1[i] : arr2[i];
   return max;
}//PMax

//______________________________________________________________________________
inline Double_t *TStat::PMin(Int_t n, const Double_t *arr1, const Double_t *arr2,
                        Double_t *min)
{
  // Returns minimum of arr1 and arr2 at each index in array min.

   for (Int_t i=0; i<n; i++) min[i] = (arr1[i] < arr2[i]) ? arr1[i] : arr2[i];
   return min;
}//PMin

//______________________________________________________________________________
inline Double_t *TStat::CumMax(Int_t n, const Double_t *arr, Double_t *max)
{
  // Returns cumulative maxima of arr in array max.

   max[0] = arr[0];
   for (Int_t i=1; i<n; i++) max[i] = (arr[i] > max[i-1]) ? arr[i] : max[i-1];
   return max;
}//CumMax

//______________________________________________________________________________
inline Double_t *TStat::CumMin(Int_t n, const Double_t *arr, Double_t *min)
{
  // Returns cumulative minima of arr in array min.

   min[0] = arr[0];
   for (Int_t i=1; i<n; i++) min[i] = (arr[i] < min[i-1]) ? arr[i] : min[i-1];
   return min;
}//CumMin

//______________________________________________________________________________
inline Double_t TStat::CumProd(Int_t n, const Double_t *arr)
{
   // Cumulative sum.

   Double_t prod = 0;
   for (Int_t i=0; i<n; i++) prod *= arr[i];
   return prod;
}//CumProd

//______________________________________________________________________________
inline Double_t TStat::CumSum(Int_t n, const Double_t *arr)
{
   // Cumulative sum.

   Double_t sum = 0;
   for (Int_t i=0; i<n; i++) sum += arr[i];
   return sum;
}//CumSum

//______________________________________________________________________________
inline Double_t TStat::LnFact(Int_t n)
{
   // ln(n!)

   if (n > 1) return TMath::LnGamma(n+1.0);
   else       return 0.0;
}//LnFact

//______________________________________________________________________________
inline Double_t TStat::BinomCoeff(Int_t n, Int_t k)
{
  // Binomial coefficient "n chooose k" as double.

  return floor(0.5 + exp(LnFact(n)-LnFact(k)-LnFact(n-k)));
}//CumProd

//______________________________________________________________________________
inline Double_t TStat::StDev(Int_t n, const Double_t *arr, const Double_t mean)
{
   // Calculate standard deviation for arithmetic mean

   return TMath::Sqrt(Var(n, arr, mean));
}//StDev


/*******************************************************************************
Class TMEstimator and its inherited classes implement "M-estimators", see:
http://en.wikipedia.org/wiki/Robust_statistics
http://en.wikipedia.org/wiki/M-estimator

Class methods Rho(), Psi(), Weight() are implemented as described in Table 1 of:
http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html

However, implementation of methods Psi(), Weight(), Derivative() is adapted from:
Bioconductor package "affyPLM" file "src/psi_fns.c"
created by: B. M. Bolstad <bolstad@stat.berkeley.edu>
*******************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMEstimator                                                          //
//                                                                      //
// Base class for M-estimators                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TMEstimator {

   protected:
      Double_t fConst;  // tuning constant
      TString  fName;   // name of estimator
      
   public:
   enum {
      kRho    = 0,   // rho function
      kPsi    = 1,   // psi function
      kDeriv  = 2,   // derivative of psi function
      kWeight = 3    // weight
   };    
      
   public:
      TMEstimator(): fConst(), fName()        {}
      TMEstimator(Double_t c, const char *name): fConst(c), fName(name) {}
      virtual ~TMEstimator()                  {}

      virtual Double_t Rho(Double_t x)        {return 0;}
      virtual Double_t Psi(Double_t x)        {return 0;}
      virtual Double_t Derivative(Double_t x) {return 0;}
      virtual Double_t Weight(Double_t x)     {return 0;}

      Double_t Calculate(Double_t x, Int_t type);
      const char *GetName()             const {return fName;}

      static TMEstimator *Estimator(const char *name, Double_t c = 0.0);
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THuberEstimator                                                      //
//                                                                      //
// Class for M-estimator Huber                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class THuberEstimator: public TMEstimator {

   public:
      THuberEstimator();
      THuberEstimator(Double_t c);
      virtual ~THuberEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t THuberEstimator::Rho(Double_t x)
{
   return (TMath::Abs(x) <= fConst) ? x*x/2.0 : fConst*(TMath::Abs(x) - fConst/2.0);
}//Rho

//______________________________________________________________________________
inline Double_t THuberEstimator::Psi(Double_t x)
{
   return (TMath::Abs(x) <= fConst) ? x : ((x < 0) ? -fConst : fConst);
}//Rho

//______________________________________________________________________________
inline Double_t THuberEstimator::Derivative(Double_t x)
{
   return (TMath::Abs(x) <= fConst) ? 1.0 : 0.0;
}//Rho

//______________________________________________________________________________
inline Double_t THuberEstimator::Weight(Double_t x)
{
   return (1 < fConst/TMath::Abs(x)) ?  1.0 : fConst/TMath::Abs(x);
}//Rho

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFairEstimator                                                       //
//                                                                      //
// Class for M-estimator Fair                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TFairEstimator: public TMEstimator {

   public:
      TFairEstimator();
      TFairEstimator(Double_t c);
      virtual ~TFairEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t TFairEstimator::Rho(Double_t x)
{
   return fConst*fConst*(TMath::Abs(x)/fConst - TMath::Log(1.0 + TMath::Abs(x)/fConst));
}//Rho

//______________________________________________________________________________
inline Double_t TFairEstimator::Psi(Double_t x)
{
   return x/(1.0 + TMath::Abs(x)/fConst);
}//Rho

//______________________________________________________________________________
inline Double_t TFairEstimator::Derivative(Double_t x)
{
   Double_t u = 1.0 + TMath::Abs(x)/fConst;
   return (x >= 0) ? 1.0/u - x/(fConst*u*u) : 1.0/u + x/(fConst*u*u);
}//Rho

//______________________________________________________________________________
inline Double_t TFairEstimator::Weight(Double_t x)
{
   return 1.0/(1.0 + TMath::Abs(x)/fConst);
}//Rho

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCauchyEstimator                                                     //
//                                                                      //
// Class for M-estimator Cauchy                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TCauchyEstimator: public TMEstimator {

   public:
      TCauchyEstimator();
      TCauchyEstimator(Double_t c);
      virtual ~TCauchyEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t TCauchyEstimator::Rho(Double_t x)
{
   return fConst*fConst*TMath::Log(1.0 + (x/fConst)*(x/fConst))/2.0;
}//Rho

//______________________________________________________________________________
inline Double_t TCauchyEstimator::Psi(Double_t x)
{
   return x/(1.0 + (x/fConst)*(x/fConst));
}//Rho

//______________________________________________________________________________
inline Double_t TCauchyEstimator::Derivative(Double_t x)
{
   Double_t u = fConst*fConst + x*x;
   return fConst*fConst*(fConst*fConst - x*x)/(u*u);
}//Rho

//______________________________________________________________________________
inline Double_t TCauchyEstimator::Weight(Double_t x)
{
   return 1.0/(1.0 + (x/fConst)*(x/fConst));
}//Rho

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGemanMcClureEstimator                                               //
//                                                                      //
// Class for M-estimator GemanMcClure                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TGemanMcClureEstimator: public TMEstimator {

   public:
      TGemanMcClureEstimator();
      TGemanMcClureEstimator(Double_t c);
      virtual ~TGemanMcClureEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t TGemanMcClureEstimator::Rho(Double_t x)
{
   return (x*x/2.0)/(1.0 + x*x);
}//Rho

//______________________________________________________________________________
inline Double_t TGemanMcClureEstimator::Psi(Double_t x)
{
   return x/((1.0 + x*x)*(1.0 + x*x));
}//Rho

//______________________________________________________________________________
inline Double_t TGemanMcClureEstimator::Derivative(Double_t x)
{
   Double_t u = 1.0 + x*x;
   return (1.0 - 3.0*x*x)/(u*u*u);
}//Rho

//______________________________________________________________________________
inline Double_t TGemanMcClureEstimator::Weight(Double_t x)
{
   return 1.0/((1.0 + x*x)*(1.0 + x*x));
}//Rho

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TWelschEstimator                                                     //
//                                                                      //
// Class for M-estimator Welsch                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TWelschEstimator: public TMEstimator {

   public:
      TWelschEstimator();
      TWelschEstimator(Double_t c);
      virtual ~TWelschEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t TWelschEstimator::Rho(Double_t x)
{
   return fConst*fConst*(1.0 - TMath::Exp(-(x/fConst)*(x/fConst)))/2.0;
}//Rho

//______________________________________________________________________________
inline Double_t TWelschEstimator::Psi(Double_t x)
{
   return x*TMath::Exp(-(x/fConst)*(x/fConst));
}//Rho

//______________________________________________________________________________
inline Double_t TWelschEstimator::Derivative(Double_t x)
{
   return TMath::Exp(-(x/fConst)*(x/fConst))*(1.0 - 2.0*(x*x)/(fConst*fConst));
}//Rho

//______________________________________________________________________________
inline Double_t TWelschEstimator::Weight(Double_t x)
{
   return TMath::Exp(-(x/fConst)*(x/fConst));
}//Rho

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTukeyEstimator                                                      //
//                                                                      //
// Class for M-estimator Tukey                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TTukeyEstimator: public TMEstimator {

   public:
      TTukeyEstimator();
      TTukeyEstimator(Double_t c);
      virtual ~TTukeyEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t TTukeyEstimator::Rho(Double_t x)
{
   Double_t u = 1.0 - (x/fConst)*(x/fConst);
   return (TMath::Abs(x) <= fConst) ? fConst*fConst*(1.0 - u*u*u)/6.0 : fConst*fConst/6.0;
}//Rho

//______________________________________________________________________________
inline Double_t TTukeyEstimator::Psi(Double_t x)
{
   Double_t u = 1.0 - (x/fConst)*(x/fConst);
   return (TMath::Abs(x) <= fConst) ? x*u*u : 0.0;
}//Rho

//______________________________________________________________________________
inline Double_t TTukeyEstimator::Derivative(Double_t x)
{
   Double_t u = (x/fConst)*(x/fConst);
   return (TMath::Abs(x) <= fConst) ? (1.0 - u)*(1.0 - 5.0*u) : 0.0;
}//Rho

//______________________________________________________________________________
inline Double_t TTukeyEstimator::Weight(Double_t x)
{
   Double_t u = 1.0 - (x/fConst)*(x/fConst);
   return (TMath::Abs(x) <= fConst) ?  u*u : 0.0;
}//Rho

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAndrewEstimator                                                     //
//                                                                      //
// Class for M-estimator Andrew                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TAndrewEstimator: public TMEstimator {

   public:
      TAndrewEstimator();
      TAndrewEstimator(Double_t c);
      virtual ~TAndrewEstimator();

      virtual Double_t Rho(Double_t x);
      virtual Double_t Psi(Double_t x);
      virtual Double_t Derivative(Double_t x);
      virtual Double_t Weight(Double_t x);
};

//______________________________________________________________________________
inline Double_t TAndrewEstimator::Rho(Double_t x)
{
//TO DO:
   return (TMath::Abs(x) <= fConst*TMath::Pi()) ? -99999999.9 : -99999999.9;
}//Rho

//______________________________________________________________________________
inline Double_t TAndrewEstimator::Psi(Double_t x)
{
   return (TMath::Abs(x) <= fConst*TMath::Pi()) ? fConst*TMath::Sin(x/fConst) : 0.0;
}//Rho

//______________________________________________________________________________
inline Double_t TAndrewEstimator::Derivative(Double_t x)
{
   return (TMath::Abs(x) <= fConst*TMath::Pi()) ? TMath::Cos(x/fConst) : 0.0;
}//Rho

//______________________________________________________________________________
inline Double_t TAndrewEstimator::Weight(Double_t x)
{
   return (TMath::Abs(x) <= fConst*TMath::Pi()) ?  TMath::Sin(x/fConst)/(x/fConst) : 0.0;
}//Rho

#endif
