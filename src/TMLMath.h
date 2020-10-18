// Author: Christian Stratowa 11/25/2002             last modified: 12/18/2005

/*
 *******************************************************************************
 ***********************  Statistics Package for ROOT  *************************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2008 Dr. Christian Stratowa
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

#ifndef __TMLMath__
#define __TMLMath__

//include to define IEEE_754, WORDS_BIGENDIAN, WIN32
#include "RconfigR.h"

#include <Riostream.h>
#include <cmath>

// R includes
#include "RmathR.h"

// R
extern double NA_REAL;
extern double R_PosInf;
extern double R_NegInf;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMLMath                                                              //
//                                                                      //
// Class containing special functions from R Mathlib                    //
// Mathlib : A C Library of Special Functions.                          //
// Copyright (C) 1998 Ross Ihaka                                        //
// Copyright (C) 2000 The R Development Core Team                       //
// http:/cran.r-project.org/                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TMLMath {
   private:
      static Double_t D1Mach(Int_t i);
      static Int_t    I1Mach(Int_t i);
      static Double_t PBetaRaw(Double_t x, Double_t pin, Double_t qin,
                         Int_t lower_tail);
      static void     PNormBoth(Double_t x, Double_t &cum, Double_t &ccum, Int_t i_tail, Int_t log_p);
      
   public :
      TMLMath() {}
      virtual ~TMLMath() {}

	// Misc
      static Double_t Expm1(Double_t x);
      static Double_t FMax2(Double_t x, Double_t y);
      static Double_t FMin2(Double_t x, Double_t y);
      static Double_t FMod(Double_t x1, Double_t x2);
      static Double_t FTrunc(Double_t x);
      static Double_t Log(Double_t x);
      static Double_t Log1p(Double_t x);
      static Double_t Pow(Double_t x, Double_t y);
      static Double_t PowDi(Double_t x, Int_t n);

	// Not a Number
      static Int_t    IsNaN(Double_t x);
      static Int_t    Finite(Double_t x);

	// Beta Distribution
      static Double_t PBeta(Double_t x, Double_t pin, Double_t qin,
                         Int_t lower_tail = kTRUE, Int_t log_p = kFALSE);

	// Normal Distribution
      static Double_t DNorm(Double_t x, Double_t mu = 0, Double_t sigma = 1,
                         Int_t give_log = kFALSE);
      static Double_t PNorm(Double_t x, Double_t mu = 0, Double_t sigma = 1,
                         Int_t lower_tail = kTRUE, Int_t log_p = kFALSE);
      static Double_t QNorm(Double_t p, Double_t mu = 0, Double_t sigma = 1,
                         Int_t lower_tail = kTRUE, Int_t log_p = kFALSE);

	// Student t Distribution
      static Double_t PT(Double_t x, Double_t n,
                         Int_t lower_tail = kTRUE, Int_t log_p = kFALSE);
      static Double_t QT(Double_t p, Double_t ndf,
                         Int_t lower_tail = kTRUE, Int_t log_p = kFALSE);

	// Gamma and Related Functions
      static Double_t GammaFn(Double_t x);
      static void     GammaLims(Double_t &xmin, Double_t &xmax);
      static Double_t LGammaCor(Double_t x);
      static Double_t LGammaFn(Double_t x);
      static Double_t LBeta(Double_t a, Double_t b);

   // Chebyshev Series
      static Int_t    ChebyshevInit(Double_t *dos, Int_t nos, Double_t eta);
      static Double_t ChebyshevEval(Double_t x, const Double_t *a, const Int_t n);

      ClassDef(TMLMath,0) //NMath
};

//______________________________________________________________________________
inline Double_t TMLMath::FMax2(Double_t x, Double_t y)
{
   // Return aximum of x and y.
   // Converted from: R-1.6.1/src/nmath/fmax2.c 

#ifdef IEEE_754
	if (IsNaN(x) || IsNaN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}//FMax2

//______________________________________________________________________________
inline Double_t TMLMath::FMin2(Double_t x, Double_t y)
{
   // Return aximum of x and y.
   // Converted from: R-1.6.1/src/nmath/fmin2.c 

#ifdef IEEE_754
	if (IsNaN(x) || IsNaN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}//FMin2

//______________________________________________________________________________
inline Double_t TMLMath::FMod(Double_t x1, Double_t x2)
{
   // Modulo.
   // Converted from: R-1.6.1/src/nmath/mlutils.c 

   Double_t q = x1 / x2;
   return x1 - floor(q) * x2;
}//FMod

//______________________________________________________________________________
inline Double_t TMLMath::FTrunc(Double_t x)
{
   // Truncation toward zero.

   if(x >= 0) return floor(x);
   else       return ceil(x);
}//FTrunc


#endif
