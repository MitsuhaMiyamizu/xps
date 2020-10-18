// Author: Christian Stratowa 11/25/2002             last modified: 06/25/2005

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

using namespace std;

// R includes:
#include "RnmathR.h"
#include "RdpqR.h"  

// ROOT
#include "TMath.h"
#include "TRandom.h"
#include "TList.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TFriendElement.h"

#include "TMLMath.h"

//debug: print class method names
const Bool_t  kCS  = 0;
const Bool_t  kCSa = 0;

// R
double NA_REAL  = ML_NAN;
double R_PosInf = ML_POSINF;
double R_NegInf = ML_NEGINF;

ClassImp(TMLMath);

Int_t gSignGam;   //defined in LGammaFn(x)


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMLMath                                                              //
//                                                                      //
// Class containing special functions from R Mathlib                    //
// This source is based on Mathlib : A C Library of Special Functions.  //
// Copyright (C) 1998 Ross Ihaka                                        //
// Copyright (C) 2000 The R Development Core Team                       //
// http:/cran.r-project.org/                                            //
// Code converted to C++ by C. Stratowa.                                //
//                                                                      //
// R and Mathlib is free software; you can redistribute it and/or modify//
// it under the terms of the GNU General Public License as published by //
// the Free Software Foundation; either version 2 of the License, or    //
// (at your option) any later version.                                  //
//                                                                      //
// R and Mathlib is distributed in the hope that it will be useful,     //
// but WITHOUT ANY WARRANTY; without even the implied warranty of       //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        //
// GNU General Public License for more details.                         //
//                                                                      //
// You should have received a copy of the GNU General Public License    //
// along with R and Mathlib; if not, write to:                          //
// The Free Software  Foundation, Inc.,                                 //
// 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
/*
//______________________________________________________________________________
Double_t TMLMath::XXX()
{
   // XXX
   if(kCSa) cout << "------TMLMath::XXX------" << endl;

}//XXX
*/

//______________________________________________________________________________
Double_t TMLMath::Expm1(Double_t x)
{
   // DESCRIPTION
   //
   //	Compute the Exponential minus 1:  exp(x) - 1
   // accurately also when x is close to zero, i.e. |x| << 1
   //
   // NOTES
   //	As log1p(), this is a standard function in some C libraries,
   //	particularly GNU and BSD (but is neither ISO/ANSI C nor POSIX).
   // Converted from: R-1.6.1/src/nmath/expm1.c 
   if(kCSa) cout << "------TMLMath::Expm1------" << endl;

   Double_t y = 0;
   Double_t a = fabs(x);

   if (a < DBL_EPSILON) return x;
   if (a > 0.697)       return exp(x) - 1;  //negligible cancellation

   if (a > 1e-8) y = exp(x) - 1;
   else          y = (x / 2 + 1) * x; //Taylor expansion, more accurate in this s range

   // Newton step for solving   log(1 + y) = x   for y :
   // WARNING: does not work for y ~ -1: bug in 1.5.0
   y -= (1 + y) * (Log1p (y) - x);

   return y;
}//Expm1

//______________________________________________________________________________
Double_t TMLMath::Log(Double_t x)
{
   // Compute logarithm
   // Converted from: R-1.6.1/src/nmath/mlutils.c 

   return (x > 0 ? log(x) : x < 0 ? ML_NAN : ML_NEGINF);
}//Log

//______________________________________________________________________________
Double_t TMLMath::Log1p(Double_t x)
{
   //Compute the relative error logarithm: log(1 + x)
   //
   // NOTES
   //    This code is a translation of the Fortran subroutine `dlnrel'
   //    written by W. Fullerton of Los Alamos Scientific Laboratory.
   // Converted from: R-1.6.1/src/nmath/log1p.c 
   if(kCSa) cout << "------TMLMath::Log1p------" << endl;

   // Series for log1p on the interval -.375 to .375
   // with weighted error   6.35e-32
   // log weighted error  31.20
   // significant figures required  30.93
   // decimal places required  32.01
   const Double_t alnrcs[43] = {
      +.10378693562743769800686267719098e+1,
      -.13364301504908918098766041553133e+0,
      +.19408249135520563357926199374750e-1,
      -.30107551127535777690376537776592e-2,
      +.48694614797154850090456366509137e-3,
      -.81054881893175356066809943008622e-4,
      +.13778847799559524782938251496059e-4,
      -.23802210894358970251369992914935e-5,
      +.41640416213865183476391859901989e-6,
      -.73595828378075994984266837031998e-7,
      +.13117611876241674949152294345011e-7,
      -.23546709317742425136696092330175e-8,
      +.42522773276034997775638052962567e-9,
      -.77190894134840796826108107493300e-10,
      +.14075746481359069909215356472191e-10,
      -.25769072058024680627537078627584e-11,
      +.47342406666294421849154395005938e-12,
      -.87249012674742641745301263292675e-13,
      +.16124614902740551465739833119115e-13,
      -.29875652015665773006710792416815e-14,
      +.55480701209082887983041321697279e-15,
      -.10324619158271569595141333961932e-15,
      +.19250239203049851177878503244868e-16,
      -.35955073465265150011189707844266e-17,
      +.67264542537876857892194574226773e-18,
      -.12602624168735219252082425637546e-18,
      +.23644884408606210044916158955519e-19,
      -.44419377050807936898878389179733e-20,
      +.83546594464034259016241293994666e-21,
      -.15731559416479562574899253521066e-21,
      +.29653128740247422686154369706666e-22,
      -.55949583481815947292156013226666e-23,
      +.10566354268835681048187284138666e-23,
      -.19972483680670204548314999466666e-24,
      +.37782977818839361421049855999999e-25,
      -.71531586889081740345038165333333e-26,
      +.13552488463674213646502024533333e-26,
      -.25694673048487567430079829333333e-27,
      +.48747756066216949076459519999999e-28,
      -.92542112530849715321132373333333e-29,
      +.17578597841760239233269760000000e-29,
      -.33410026677731010351377066666666e-30,
      +.63533936180236187354180266666666e-31,
   };
   const Double_t xmin = -1 + sqrt(1/DBL_EPSILON);

#ifdef NOMORE_FOR_THREADS
   static Int_t nlnrel = 0;

   if (nlnrel == 0) {
   // initialize chebychev coefficients
      nlnrel = ChebyshevInit(alnrcs, 43, DBL_EPSILON/20);
   }//if
#else
# define nlnrel 22
/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
#endif

   if (x == 0.) return 0.;/* speed */
   if (x == -1) return(ML_NEGINF);
   if (x  < -1) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

   if (fabs(x) <= .375) {
   // Improve on speed (only);
	// again give result accurate to IEEE double precision:
      if (fabs(x) < .5 * DBL_EPSILON) return x;

      if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
         return x * (1 - .5 * x);

      return x * (1 - x * ChebyshevEval(x / .375, alnrcs, nlnrel));
   }//if

   if (x < xmin) {
      // answer less than half precision because x too near -1
      ML_ERROR(ME_PRECISION);
   }//if

   return Log(1 + x);
}//Log1p

//______________________________________________________________________________
Double_t TMLMath::Pow(Double_t x, Double_t y)
{
   // Power x ^ y
   // Converted from: R-1.6.1/src/nmath/mlutils.c 
   if(kCSa) cout << "------TMLMath::Pow------" << endl;

   if(x == 1. || y == 0.)
      return(1.);

   if(x == 0.) {
      if(y > 0.) return(0.);
      /* y < 0 */return(ML_POSINF);
   }//if

   if (Finite(x) && Finite(y))
      return(pow(x,y));  //pow(), not Pow()!!!

   if (IsNaN(x) || IsNaN(y)) {
#ifdef IEEE_754
      return(x + y);
#else
      return(NA_REAL);
#endif
   }//if

   if(!Finite(x)) {
      if(x > 0){		/* Inf ^ y */
         return((y < 0.)? 0. : ML_POSINF);
      } else {			/* (-Inf) ^ y */
         if(Finite(y) && y == floor(y)) /* (-Inf) ^ n */
            return((y < 0.) ? 0. : ((Int_t)FMod(y,2.) ? x  : -x));
      }//if
   }//if

   if(!Finite(y)) {
      if(x >= 0) {
         if(y > 0)		/* y == +Inf */
            return((x >= 1)? ML_POSINF : 0.);
         else		/* y == -Inf */
            return((x < 1) ? ML_POSINF : 0.);
      }//if
   }//if

   return(ML_NAN);		/* all other cases: (-Inf)^{+-Inf,
				   non-int}; (neg)^{+-Inf} */
}//Pow

//______________________________________________________________________________
Double_t TMLMath::PowDi(Double_t x, Int_t n)
{
   // Power x ^ y
   // Converted from: R-1.6.1/src/nmath/mlutils.c 
   if(kCSa) cout << "------TMLMath::PowDi------" << endl;

   Double_t pow = 1.0;

   if (IsNaN(x)) return x;

   if (n != 0) {
      if (!Finite(x)) return Pow(x, (Double_t)n);
      if (n < 0) { n = -n; x = 1/x; }
      for (;;) {
         if(n & 01) pow *= x;
         if(n >>= 1) x *= x; else break;
      }//for
   }//if

   return pow;
}//PowDi

//______________________________________________________________________________
Int_t TMLMath::IsNaN(Double_t x)
{
   // Not a Number.
   // Converted from: R-1.6.1/src/nmath/mlutils.c 

#ifdef IEEE_754
   return (isnan(x) != 0);
#else /* not IEEE_754 */
# ifndef HAVE_ISNAN
   return (x == ML_NAN);
# else
   return (isnan(x) != 0 || x == ML_NAN);
# endif
#endif /* not IEEE_754 */
}//IsNA

//______________________________________________________________________________
Int_t TMLMath::Finite(Double_t x)
{
   // Number not finite
   // Converted from: R-1.6.1/src/nmath/mlutils.c 

#ifdef IEEE_754

/* Include the header file defining finite() */
#ifdef HAVE_IEEE754_H
# include <ieee754.h>		/* newer Linuxen */
#else
# ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>		/* others [Solaris 2.5.x], .. */
# endif
#endif
#if defined(Win32) && defined(_MSC_VER)
# include <float.h>
#endif

# ifdef Macintosh
   return isfinite(x);
# endif
# ifdef HAVE_WORKING_FINITE
   return finite(x);
# else
#  ifdef _AIX
#   include <fp.h>
    return FINITE(x);
#  else
    return (!isnan(x) & (x != ML_POSINF) & (x != ML_NEGINF));
#  endif
# endif

#else /* not IEEE_754 */

# ifndef HAVE_FINITE
   return (x !=  ML_NAN && x < ML_POSINF && x > ML_NEGINF);
# else
   Int_t finite(Double_t);
   return finite(x);
# endif

#endif /* not IEEE_754 */
}//Finite

//______________________________________________________________________________
Double_t TMLMath::PBeta(Double_t x, Double_t pin, Double_t qin,
         Int_t lower_tail, Int_t log_p)
{
   // DESCRIPTION
   // Returns distribution function of the beta distribution.
   // ( = The incomplete beta ratio I_x(p,q) ).
   // Converted from: R-1.6.1/src/nmath/pbeta.c 
   if(kCSa) cout << "------TMLMath::PBeta------" << endl;

#ifdef IEEE_754
   if (IsNaN(x) || IsNaN(pin) || IsNaN(qin)) return x + pin + qin;
#endif

   if (pin <= 0 || qin <= 0) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

   if (x <= 0) return R_DT_0;
   if (x >= 1) return R_DT_1;

   return R_D_val(PBetaRaw(x, pin, qin, lower_tail));
}//PBeta

//______________________________________________________________________________
Double_t TMLMath::DNorm(Double_t x, Double_t mu, Double_t sigma, Int_t give_log)
{
   // DESCRIPTION
   // Compute the density of the normal distribution.
   // Converted from: R-1.6.1/src/nmath/dnorm.c 
   if(kCSa) cout << "------TMLMath::DNorm------" << endl;

#ifdef IEEE_754
   if (IsNaN(x) || IsNaN(mu) || IsNaN(sigma)) return x + mu + sigma;
#endif
   if (sigma <= 0) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

   x = (x - mu) / sigma;

   // M_1_SQRT_2PI = 1 / sqrt(2 * pi)
   return (give_log ?
            -(M_LN_SQRT_2PI  +	0.5 * x * x + Log(sigma)) :
              M_1_SQRT_2PI * exp(-0.5 * x * x)  /	  sigma);
}//DNorm

//______________________________________________________________________________
Double_t TMLMath::PNorm(Double_t x, Double_t mu, Double_t sigma,
         Int_t lower_tail, Int_t log_p)
{
   //  DESCRIPTION
   //	The main computation evaluates near-minimax approximations derived
   //	from those in "Rational Chebyshev approximations for the error
   //	function" by W. J. Cody, Math. Comp., 1969, 631-637.  This
   //	transportable program uses rational functions that theoretically
   //	approximate the normal distribution function to at least 18
   //	significant decimal digits.  The accuracy achieved depends on the
   //	arithmetic system, the compiler, the intrinsic functions, and
   //	proper selection of the machine-dependent constants.
   //
   //  REFERENCE
   //	Cody, W. D. (1993).
   //	ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
   //	Special Function Routines and Test Drivers".
   //	ACM Transactions on Mathematical Software. 19, 22-32.
   // Converted from: R-1.6.1/src/nmath/pnorm.c 
   if(kCSa) cout << "------TMLMath::PNorm------" << endl;

// Note: The structure of these checks has been carefully thought through.
// For example, if x == mu and sigma == 0, we still get the correct answer.
#ifdef IEEE_754
   if (IsNaN(x) || IsNaN(mu) || IsNaN(sigma))
	   return x + mu + sigma;
#endif
   if (sigma < 0) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

   x = (x - mu) / sigma;
   if (!Finite(x)) {
	   if (IsNaN(x)) {return(ML_NAN);}
	   if(x < 0)     {return R_DT_0;}
	   else          {return R_DT_1;}
   }//if

   Double_t p  = 0.;
   Double_t cp = 0.;
   PNormBoth(x, p, cp, (lower_tail ? 0 : 1), log_p);

   return (lower_tail ? p : cp);
}//PNorm

//______________________________________________________________________________
Double_t TMLMath::QNorm(Double_t p, Double_t mu, Double_t sigma,
         Int_t lower_tail, Int_t log_p)
{
   // DESCRIPTION
   // Compute the quantile function for the normal distribution.
   // For small to moderate probabilities, algorithm referenced
   // below is used to obtain an initial approximation which is
   // polished with a final Newton step.
   // 
   // For very large arguments, an algorithm of Wichura is used.
   // 
   // REFERENCE
   // Beasley, J. D. and S. G. Springer (1977).
   // Algorithm AS 111: The percentage points of the normal distribution,
   // Applied Statistics, 26, 118-121.
   // 
   //      Wichura, M.J. (1988).
   //      Algorithm AS 241: The Percentage Points of the Normal Distribution.
   //      Applied Statistics, 37, 477-484.
   // Converted from: R-1.6.1/src/nmath/qnorm.c 
   if(kCSa) cout << "------TMLMath::QNorm------" << endl;

   Double_t p_, q, r, val;

#ifdef IEEE_754
   if (IsNaN(p) || IsNaN(mu) || IsNaN(sigma)) return p + mu + sigma;
#endif
   if (p == R_DT_0)	return ML_NEGINF;
   if (p == R_DT_1)	return ML_POSINF;
   R_Q_P01_check(p);

   if(sigma  < 0)	{ML_ERROR(ME_DOMAIN); return ML_NAN;}
   if(sigma == 0)	return mu;

   p_ = R_DT_qIv(p);  // real lower_tail prob. p
   q = p_ - 0.5;

#ifdef DEBUG_qnorm
   REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
            p,mu,sigma, lower_tail, log_p, q);
#endif

#ifdef OLD_qnorm
   /* --- use  AS 111 --- */
   if (fabs(q) <= 0.42) {
   // 0.08 <= p <= 0.92
      r = q * q;
      val = q * (((-25.44106049637 * r + 41.39119773534) * r
            - 18.61500062529) * r + 2.50662823884)
            / ((((3.13082909833 * r - 21.06224101826) * r
            + 23.08336743743) * r + -8.47351093090) * r + 1.0);
   } else {
	// p < 0.08 or p > 0.92, set r = min(p, 1 - p)
      if (q > 0) {r = R_DT_CIv(p);}  //1-p
      else       {r = p_;}           // = R_DT_Iv(p) ^=  p
#ifdef DEBUG_qnorm
	REprintf("\t 'middle p': r = %7g\n", r);
#endif
      if(r > DBL_EPSILON) {
         r = sqrt(- ((log_p && 
            ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : Log(r)));
#ifdef DEBUG_qnorm
         REprintf("\t new r = %7g ( =? sqrt(- log(r)) )\n", r);
#endif
         val = (((2.32121276858 * r + 4.85014127135) * r
               - 2.29796479134) * r - 2.78718931138)
               / ((1.63706781897 * r + 3.54388924762) * r + 1.0);
         if (q < 0) val = -val;
      } else if (r >= DBL_MIN) { /* r = p <= eps : Use Wichura */
         val = -2 * (log_p ? R_D_Lval(p) : Log(R_D_Lval(p)));
         r = Log(2 * M_PI * val);
#ifdef DEBUG_qnorm
         REprintf("\t DBL_MIN <= r <= DBL_EPS: val = %g, new r = %g\n", val, r);
#endif
         p = val * val;
         r = r/val + (2 - r)/p + (-14 + 6 * r - r * r)/(2 * p * val);
         val = sqrt(val * (1 - r));
         if (q < 0.0) val = -val;
         return mu + sigma * val;
      } else {
#ifdef DEBUG_qnorm
         REprintf("\t r < DBL_MIN : giving up (-> +- Inf \n");
#endif
         ML_ERROR(ME_RANGE);
         if(q < 0.0) {return ML_NEGINF;}
         else        {return ML_POSINF;}
      }//if
   }//if
/* FIXME: This could be improved when log_p or !lower_tail ?
 *	  (using p, not p_ , and a different derivative )
 */
#ifdef DEBUG_qnorm
   REprintf("\t before final step: val = %7g\n", val);
#endif
   /* Final Newton step: */
   val = val -
         (PNorm(val, 0., 1., /*lower*/TRUE, /*log*/FALSE) - p_) /
          DNorm(val, 0., 1., /*log*/FALSE);

#else
/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*    ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

      Produces the normal deviate Z corresponding to a given lower
      tail area of P; Z is accurate to about 1 part in 10**16.

      (original fortran code used PARAMETER(..) for the coefficients
      and provided hash codes for checking them...)
*/
   if (fabs(q) <= .425) {  // 0.075 <= p <= 0.925
      r = .180625 - q * q;
      val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
   } else {  // closer than 0.075 from {0,1} boundary
   // r = min(p, 1-p) < 0.075
      if (q > 0) {r = R_DT_CIv(p);} //1-p 
      else       {r = p_;}          // = R_DT_Iv(p) ^=  p

      r = sqrt(- ((log_p &&
         ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : Log(r)));
      // r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 )
#ifdef DEBUG_qnorm
   REprintf("\t close to 0 or 1: r = %7g\n", r);
#endif
      if (r <= 5.) { // <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11
         r += -1.6;
         val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
      } else {  // very close to  0 or 1
         r += -5.;
         val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
      }//if
   if(q < 0.0) val = -val;  // return (q >= 0.)? r : -r ;
   }//if

#endif
/*-- Switch of AS 111 <-> AS 241 --- */
   return mu + sigma * val;
}//QNorm

//______________________________________________________________________________
Double_t TMLMath::PT(Double_t x, Double_t n, Int_t lower_tail, Int_t log_p)
{
   // The Student t distribution
   // Return P[ T <= x ] where T ~ t_{n} (t distrib. with n degrees of freedom)
   // Converted from: R-1.6.1/src/nmath/pt.c 
   if(kCSa) cout << "------TMLMath::PT------" << endl;

   Double_t val = 0;
#ifdef IEEE_754
   if (IsNaN(x) || IsNaN(n)) return x + n;
#endif

   if (n <= 0.0) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

   if (!Finite(x)) return (x < 0) ? R_DT_0 : R_DT_1;

   if (!Finite(n)) return PNorm(x, 0.0, 1.0, lower_tail, log_p);

   if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
   // Approx. from	 Abramowitz & Stegun 26.7.8 (p.949)
      val = 1./(4.*n);
      return PNorm(x*(1.- val)/sqrt(1.+ x*x*2.*val),0.0,1.0,lower_tail, log_p);
   }//if

   val = PBeta(n / (n + x * x), n / 2.0, 0.5, /*lower_tail*/1, log_p);

   // Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):
   if (x <= 0.) lower_tail = !lower_tail;

   if (log_p) {
      if(lower_tail) {return Log1p(-0.5*exp(val));}
      else           {return val - M_LN2;}    // = log(.5* pbeta(....))
   } else {
      val /= 2.;
      return R_D_Cval(val);
   }//if
}//PT

//______________________________________________________________________________
Double_t TMLMath::QT(Double_t p, Double_t ndf, Int_t lower_tail, Int_t log_p)
{
   // DESCRIPTION
   // The "Student" t distribution quantile function.
   // 
   // NOTES
   // This is a C translation of the Fortran routine given in:
   // Algorithm 396: Student's t-quantiles by
   // G.W. Hill CACM 13(10), 619-620, October 1970
   // Converted from: R-1.6.1/src/nmath/qt.c 
   if(kCSa) cout << "------TMLMath::QT------" << endl;

   const Double_t eps = 1.e-12;

   Double_t a, b, c, d, p_, P, q, x, y;
   Bool_t   neg;

#ifdef IEEE_754
   if (IsNaN(p) || IsNaN(ndf)) return p + ndf;
#endif
   if (p == R_DT_0) return ML_NEGINF;
   if (p == R_DT_1) return ML_POSINF;
   R_Q_P01_check(p);

//??   if (ndf < 1) ML_ERR_return_NAN;  /* FIXME:  not yet treated here */
   if (ndf < 1) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

    /* FIXME: This test should depend on  ndf  AND p  !!
     * -----  and in fact should be replaced by
     * something like Abramowitz & Stegun 26.7.5 (p.949)
     */
   if (ndf > 1e20) return QNorm(p, 0., 1., lower_tail, log_p);

   p_ = R_D_qIv(p); /* note: exp(p) may underflow to 0; fix later */

   // 0 <= P <= 1  in all cases
   if((lower_tail && p_ > 0.5) || (!lower_tail && p_ < 0.5)) {
      neg = kFALSE;
      P = 2 * R_D_Cval(p_);
   } else {
      neg = kTRUE;
      P = 2 * R_D_Lval(p_);
   }//if 

   if (fabs(ndf - 2) < eps) {	//df ~= 2
      if (P > 0) {
         q = sqrt(2 / (P * (2 - P)) - 2);
      } else { //P = 0, but maybe = exp(p) !
         if (log_p) {q = M_SQRT2 * exp(- .5 * R_D_Lval(p));}
         else       {q = ML_POSINF;}
      }//if
   } else if (ndf < 1 + eps) { //df ~= 1  (df < 1 excluded above !)
      if (P > 0) {
         q = - tan((P+1) * M_PI_2);
      } else { //P = 0, but maybe p_ = exp(p) !
         if (log_p) q = M_1_PI * exp(-R_D_Lval(p)); //cot(e) ~ 1/e
         else q = ML_POSINF;
      }
   } else {	//-- usual case; including, e.g., df = 1.1
      a = 1 / (ndf - 0.5);
      b = 48 / (a * a);
      c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
      d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * M_PI_2) * ndf;
      if (P > 0 || !log_p) {
         y = Pow(d * P, 2 / ndf);
      } else {   //P = 0 && log_p;  P = 2*exp(p*)
         y = exp(2 / ndf * (Log(d) + M_LN2 + R_D_Lval(p)));
      }//if

      if (y > 0.05 + a) {
      // Asymptotic inverse expansion about normal
         if (P > 0 || !log_p) {
            x = QNorm(0.5 * P, 0., 1., kTRUE, kFALSE);
         } else {
         // P = 0 && log_p;  P = 2*exp(p') */
            x = QNorm( p, 0., 1., lower_tail, kTRUE);
         }//if

         y = x * x;
         if (ndf < 5) c += 0.3 * (ndf - 4.5) * (x + 0.6);
         c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
         y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1) * x;
         y = Expm1(a * y * y);
      } else {
         y = ((1 / (((ndf + 6) / (ndf * y) - 0.089 * d - 0.822)
             * (ndf + 2) * 3) + 0.5 / (ndf + 4))
             * y - 1) * (ndf + 1) / (ndf + 2) + 1 / y;
      }
      q = sqrt(ndf * y);
   }//if

   if (neg) q = -q;
   return q;
}//QT

//______________________________________________________________________________
Double_t TMLMath::GammaFn(Double_t x)
{
   // DESCRIPTION
   // This function computes the value of the gamma function.
   // 
   // NOTES
   // This function is a translation into C of a Fortran subroutine
   // by W. Fullerton of Los Alamos Scientific Laboratory.
   // 
   // The accuracy of this routine compares (very) favourably with
   // those of the Sun Microsystems portable mathematical library.
   // Converted from: R-1.6.1/src/nmath/gamma.c 
   if(kCSa) cout << "------TMLMath::GammaFn------" << endl;

   const Double_t gamcs[42] = {
      +.8571195590989331421920062399942e-2,
      +.4415381324841006757191315771652e-2,
      +.5685043681599363378632664588789e-1,
      -.4219835396418560501012500186624e-2,
      +.1326808181212460220584006796352e-2,
      -.1893024529798880432523947023886e-3,
      +.3606925327441245256578082217225e-4,
      -.6056761904460864218485548290365e-5,
      +.1055829546302283344731823509093e-5,
      -.1811967365542384048291855891166e-6,
      +.3117724964715322277790254593169e-7,
      -.5354219639019687140874081024347e-8,
      +.9193275519859588946887786825940e-9,
      -.1577941280288339761767423273953e-9,
      +.2707980622934954543266540433089e-10,
      -.4646818653825730144081661058933e-11,
      +.7973350192007419656460767175359e-12,
      -.1368078209830916025799499172309e-12,
      +.2347319486563800657233471771688e-13,
      -.4027432614949066932766570534699e-14,
      +.6910051747372100912138336975257e-15,
      -.1185584500221992907052387126192e-15,
      +.2034148542496373955201026051932e-16,
      -.3490054341717405849274012949108e-17,
      +.5987993856485305567135051066026e-18,
      -.1027378057872228074490069778431e-18,
      +.1762702816060529824942759660748e-19,
      -.3024320653735306260958772112042e-20,
      +.5188914660218397839717833550506e-21,
      -.8902770842456576692449251601066e-22,
      +.1527474068493342602274596891306e-22,
      -.2620731256187362900257328332799e-23,
      +.4496464047830538670331046570666e-24,
      -.7714712731336877911703901525333e-25,
      +.1323635453126044036486572714666e-25,
      -.2270999412942928816702313813333e-26,
      +.3896418998003991449320816639999e-27,
      -.6685198115125953327792127999999e-28,
      +.1146998663140024384347613866666e-28,
      -.1967938586345134677295103999999e-29,
      +.3376448816585338090334890666666e-30,
      -.5793070335782135784625493333333e-31
   };

   Int_t i, n;
   Double_t y;
   Double_t sinpiy, value;

#ifdef NOMORE_FOR_THREADS
   static Int_t    ngam = 0;
   static Double_t xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

// Initialize machine dependent constants, the first time gamma() is called.
/* FIXME for threads ! */
   if (ngam == 0) {
      ngam = ChebyshevInit(gamcs, 42, DBL_EPSILON/20);
      GammaLims(xmin, xmax);
      xsml = exp(FMax2(Log(DBL_MIN), -Log(DBL_MAX)) + 0.01);
      //   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE
      dxrel = sqrt(1/DBL_EPSILON);
   }
#else
// For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
// (xmin, xmax) are non-trivial, see ./gammalims.c 
// xsml = exp(.01)*DBL_MIN 
// dxrel = sqrt(1/DBL_EPSILON) = 2 ^ 26
//# define ngam 22
//# define xmin -170.5674972726612
//# define xmax  171.61447887182298
//# define xsml 2.2474362225598545e-308
//# define dxrel 67108864.
 const Int_t    ngam  =  22;
 const Double_t xmin  = -170.5674972726612;
 const Double_t xmax  =  171.61447887182298;
 const Double_t xsml  =  2.2474362225598545e-308;
 const Double_t dxrel =  67108864;
#endif

   if(IsNaN(x)) return x;

   y = fabs(x);

   if (y <= 10) {
	// Compute gamma(x) for -10 <= x <= 10
	// Reduce the interval and find gamma(1 + y) for 0 <= y < 1
	// first of all.
      n = (Int_t)x;
//      n = x;
      if (x < 0) --n;
      y = x - n;   // n = floor(x)  ==>	y in [ 0, 1 ) 
      --n;
      value = ChebyshevEval(y * 2 - 1, gamcs, ngam) + .9375;
      if (n == 0) return value;  // x = 1.dddd = 1+y 

		if (n < 0) {
		// compute gamma(x) for -10 <= x < 1
			// If argument is exactly zero or a negative integer then return NaN
			if (x == 0 || (x < 0 && x == n + 2)) {
				ML_ERROR(ME_RANGE);
				return ML_NAN;
			}//if

			// Answer is less than half precision because x too near a negative integer
			if (x < -0.5 && fabs(x - (Int_t)(x - 0.5) / x) < dxrel) {
				ML_ERROR(ME_PRECISION);
			}//if

			// The argument is so close to 0 that the result would overflow.
			if (y < xsml) {
				ML_ERROR(ME_RANGE);
				if (x > 0) {return ML_POSINF;}
				else       {return ML_NEGINF;}
			}//if

			n = -n;
			for (i=0; i<n; i++) value /= (x + i);
			return value;
		} else {
		// gamma(x) for 2 <= x <= 10 
			for (i=1; i<=n; i++) value *= (y + i);
			return value;
		}
	} else {
   // gamma(x) for  y = |x| > 10.
      if (x > xmax) {			// Overflow
         ML_ERROR(ME_RANGE);
         return ML_POSINF;
      }//if

      if (x < xmin) {			// Underflow 
         ML_ERROR(ME_UNDERFLOW);
         return ML_UNDERFLOW;
      }//if

      value = exp((y - 0.5) * Log(y) - y + M_LN_SQRT_2PI + LGammaCor(y));

      if (x > 0) return value;

      if (fabs((x - (Int_t)(x - 0.5))/x) < dxrel) {
      // The answer is less than half precision because 
      // the argument is too near a negative integer.
         ML_ERROR(ME_PRECISION);
      }//if

      sinpiy = sin(M_PI * y);
      if (sinpiy == 0) {    // Negative integer arg - overflow 
         ML_ERROR(ME_RANGE);
         return ML_POSINF;
      }//if

      return ((-M_PI) / (y * sinpiy * value));
   }//if
}//GammaFn

//______________________________________________________________________________
void TMLMath::GammaLims(Double_t &xmin, Double_t &xmax)
{
   // DESCRIPTION
   // This function calculates the minimum and maximum legal bounds
   // for x in GammaFn(x).  These are not the only bounds, but they
   // are the only non-trivial ones to calculate.
   // 
   // NOTES
   // This routine is a translation into C of a Fortran subroutine
   // by W. Fullerton of Los Alamos Scientific Laboratory.
   // Converted from: R-1.6.1/src/nmath/gammalims.c 
   if(kCSa) cout << "------TMLMath::GammaLims------" << endl;

// FIXME: Even better: If IEEE, #define these in RnmathR.h
//	       and don't call GammaLims() at all
#ifdef IEEE_754
   xmin = (Double_t)(-170.5674972726612);
   xmax = (Double_t)(171.61447887182298);  //(3 Intel/Sparc architectures)
//   xmin = -170.5674972726612;
//   xmax =  171.61447887182298;  //(3 Intel/Sparc architectures)
#else
   Double_t alnbig, alnsml, xln, xold;

   alnsml = Log(D1Mach(1));
   xmin = -alnsml;
   for (Int_t i=1; i<=10; ++i) {
      xold = xmin;
      xln  = Log(xmin);
      xmin -= xmin * ((xmin + .5) * xln - xmin - .2258 + alnsml) /
              (xmin * xln + .5);
      if (fabs(xmin - xold) < .005) {
         xmin = -(xmin) + .01;
         goto find_xmax;
      }//if
   }//for_i

   // unable to find xmin
   ML_ERROR(ME_NOCONV);
   xmin = xmax = ML_NAN;

find_xmax:
   alnbig = Log(D1Mach(2));
   xmax = alnbig;
   for (Int_t i=1; i<=10; ++i) {
      xold = xmax;
      xln = Log(xmax);
      xmax -= xmax * ((xmax - .5) * xln - xmax + .9189 - alnbig) /
              (xmax * xln - .5);
      if (fabs(xmax - xold) < .005) {
         xmax += -.01;
         goto done;
      }//if
   }//for_i

   // unable to find xmax
   ML_ERROR(ME_NOCONV);
   xmin = xmax = ML_NAN;

done:
   xmin = FMax2(xmin, -(xmax) + 1);
#endif
}//GammaLims

//______________________________________________________________________________
Double_t TMLMath::LGammaCor(Double_t x)
{
   // DESCRIPTION
   // Compute the log gamma correction factor for x >= 10 so that
   // log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
   // [ lgammacor(x) is called	Del(x)	in other contexts (e.g. dcdflib)]
   // 
   // NOTES
   // This routine is a translation into C of a Fortran subroutine
   // written by W. Fullerton of Los Alamos Scientific Laboratory.
   // Converted from: R-1.6.1/src/nmath/lgammacor.c 
   if(kCSa) cout << "------TMLMath::LGammaCor------" << endl;

   const Double_t algmcs[15] = {
      +.1666389480451863247205729650822e+0,
      -.1384948176067563840732986059135e-4,
      +.9810825646924729426157171547487e-8,
      -.1809129475572494194263306266719e-10,
      +.6221098041892605227126015543416e-13,
      -.3399615005417721944303330599666e-15,
      +.2683181998482698748957538846666e-17,
      -.2868042435334643284144622399999e-19,
      +.3962837061046434803679306666666e-21,
      -.6831888753985766870111999999999e-23,
      +.1429227355942498147573333333333e-24,
      -.3547598158101070547199999999999e-26,
      +.1025680058010470912000000000000e-27,
      -.3401102254316748799999999999999e-29,
      +.1276642195630062933333333333333e-30
   };

#ifdef NOMORE_FOR_THREADS
   static Int_t    nalgm = 0;
   static Double_t xbig = 0, xmax = 0;

// Initialize machine dependent constants, the first time gamma() is called.
//	FIXME for threads !
   if (nalgm == 0) {
   // For IEEE double precision : nalgm = 5
      nalgm = ChebyshevInit(algmcs, 15, DBL_EPSILON/2);
      xbig  = 1 / sqrt(DBL_EPSILON/2);   // ~ 94906265.6 for IEEE double
      xmax  = exp(FMin2(Log(DBL_MAX / 12), -Log(12 * DBL_MIN)));
      //    = DBL_MAX / 48 ~= 3.745e306 for IEEE double
   }//if
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
# define nalgm 5
# define xbig  94906265.62425156
# define xmax  3.745194030963158e306
#endif

   Double_t tmp;
   if (x < 10) {
//??      ML_ERR_return_NAN;
      ML_ERROR(ME_DOMAIN); return ML_NAN;
   } else if (x >= xmax) {
      ML_ERROR(ME_UNDERFLOW);
      return ML_UNDERFLOW;
   } else if (x < xbig) {
      tmp = 10 / x;
      return ChebyshevEval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
   } else {
      return 1 / (x * 12);
   }//if
}//LGammaCor

//______________________________________________________________________________
Double_t TMLMath::LGammaFn(Double_t x)
{
   // DESCRIPTION
   // This function computes log|gamma(x)|.  At the same time
   // the variable "signgam" is set to the sign of the gamma
   // function.
   // 
   // NOTES
   // This routine is a translation into C of a Fortran subroutine
   // by W. Fullerton of Los Alamos Scientific Laboratory.
   // 
   // The accuracy of this routine compares (very) favourably
   // with those of the Sun Microsystems portable mathematical
   // library.
   // Converted from: R-1.6.1/src/nmath/lgamma.c 
   if(kCSa) cout << "------TMLMath::LGammaFn------" << endl;

   Double_t ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
   static Double_t xmax  = 0.;
   static Double_t dxrel = 0.;
   if (xmax == 0) {// initialize machine dependent constants _ONCE_ 
      xmax  = D1Mach(2)/Log(D1Mach(2));   // = 2.533 e305  for IEEE double
      dxrel = sqrt(D1Mach(4));    // sqrt(Eps) ~ 1.49 e-8  for IEEE double 
   }//if
#else
// For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
// xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
// dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765696e-8
#endif

   gSignGam = 1;

#ifdef IEEE_754
   if (IsNaN(x)) return x;
#endif

   if (x <= 0 && x == (Int_t)x) { // Negative integer argument
      ML_ERROR(ME_RANGE);
      return ML_POSINF;   //+Inf, since lgamma(x) = log|gamma(x)|
   }//if

   y = fabs(x);

   if (y <= 10) return Log(fabs(GammaFn(x)));
   // ELSE  y = |x| > 10 ---------------------- 
   if (y > xmax) {
      ML_ERROR(ME_RANGE);
      return ML_POSINF;
   }//if

   if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
      if (x > 1e17)          {return (x*(Log(x) - 1.));}
      else if (x > 4934720.) {return (M_LN_SQRT_2PI + (x - 0.5) * Log(x) - x);}
      else
#endif
      return M_LN_SQRT_2PI + (x - 0.5) * Log(x) - x + LGammaCor(x);
   }//if
   // else: x < -10; y = -x
   sinpiy = fabs(sin(M_PI * y));

   if (sinpiy == 0) { // negative integer argument 
      MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
//??      ML_ERR_return_NAN;
      ML_ERROR(ME_DOMAIN); return ML_NAN;
   }//if

   ans = M_LN_SQRT_PId2 + (x - 0.5) * Log(y) - x - Log(sinpiy) - LGammaCor(y);

   if (fabs((x - (Int_t)(x - 0.5)) * ans / x) < dxrel) {
   // The answer is less than half precision because
   // the argument is too near a negative integer.
      ML_ERROR(ME_PRECISION);
   }//if

   if (x <= 0 && ((Int_t)(-x))%2 == 0) gSignGam = -1;

   return ans;
}//LGammaFn

//______________________________________________________________________________
Double_t TMLMath::LBeta(Double_t a, Double_t b)
{
   // DESCRIPTION
   // This function returns the value of the log beta function.
   // 
   // NOTES
   // This routine is a translation into C of a Fortran subroutine
   // by W. Fullerton of Los Alamos Scientific Laboratory.
   // Converted from: R-1.6.1/src/nmath/lbeta.c 
   if(kCSa) cout << "------TMLMath::LBeta------" << endl;

   Double_t p, q;
   p = q = a;
   if (b < p) p = b; // := min(a,b)
   if (b > q) q = b; // := max(a,b)

#ifdef IEEE_754
   if (IsNaN(a) || IsNaN(b)) return a + b;
#endif

// both arguments must be >= 0 
   if (p < 0)           {ML_ERROR(ME_DOMAIN); return ML_NAN;}
   else if (p == 0)     {return ML_POSINF;}
   else if (!Finite(q)) {return ML_NEGINF;}

   Double_t corr = 0;
   if (p >= 10) {
   // p and q are big.
      corr = LGammaCor(p) + LGammaCor(q) - LGammaCor(p + q);
      return Log(q) * -0.5 + M_LN_SQRT_2PI + corr
             + (p - 0.5) * Log(p / (p + q)) + q * Log1p(-p / (p + q));
   } else if (q >= 10) {
   // p is small, but q is big.
      corr = LGammaCor(q) - LGammaCor(p + q);
      return LGammaFn(p) + corr + p - p * Log(p + q)
             + (q - 0.5) * Log1p(-p / (p + q));
   } else {
   // p and q are small: p <= q > 10.
      return Log(GammaFn(p) * (GammaFn(q) / GammaFn(p + q)));
   }//if
}//LBeta

//______________________________________________________________________________
Int_t TMLMath::ChebyshevInit(Double_t *dos, Int_t nos, Double_t eta)
{
   // "ChebyshevInit" determines the number of terms for the
   // double precision orthogonal series "dos" needed to insure
   // the error is no larger than "eta".  Ordinarily eta will be
   // chosen to be one-tenth machine precision.
   //
   // NOTES
   //    This routine is a translation into C of Fortran routines
   //    by W. Fullerton of Los Alamos Scientific Laboratory.
   //    Based on the Fortran routine dcsevl by W. Fullerton.
   //    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
   // Converted from: R-1.6.1/src/nmath/chebyshev.c 
   if(kCSa) cout << "------TMLMath::ChebyshevInit------" << endl;

   if (nos < 1) return 0;

   Int_t    i   = 0;
   Double_t err = 0.0;
   for (Int_t ii=1; ii<=nos; ii++) {
      i = nos - ii;
      err += fabs(dos[i]);
      if (err > eta) return i;
   }//for

   return i;
}//ChebyshevInit

//______________________________________________________________________________
Double_t TMLMath::ChebyshevEval(Double_t x, const Double_t *a, const Int_t n)
{
   // "ChebyshevEval" evaluates the n-term Chebyshev series "a" at "x".
   //
   // NOTES
   //    This routine is a translation into C of Fortran routines
   //    by W. Fullerton of Los Alamos Scientific Laboratory.
   //    Based on the Fortran routine dcsevl by W. Fullerton.
   //    Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
   // Converted from: R-1.6.1/src/nmath/chebyshev.c 
   if(kCSa) cout << "------TMLMath::ChebyshevEval------" << endl;

   if (n < 1 || n > 1000)   {ML_ERROR(ME_DOMAIN); return ML_NAN;}
   if (x < -1.1 || x > 1.1) {ML_ERROR(ME_DOMAIN); return ML_NAN;}

   Double_t b2 = 0;
   Double_t b1 = 0;
   Double_t b0 = 0;
   Double_t twox = x * 2;
   for (Int_t i=1;i<=n;i++) {
      b2 = b1;
      b1 = b0;
      b0 = twox * b1 - b2 + a[n - i];
   }//for_i

   return (b0 - b2) * 0.5;
}//ChebyshevEval


//______________________________________________________________________________
Double_t TMLMath::D1Mach(Int_t i)
{
   // NaNs propagated correctly 
   // ??? FIXME:  Eliminate calls to these
   //  =====   o   from C code when
   //      o   it is only used to initialize "static" variables (threading)
   //  and use the DBL_... constants instead
   // Converted from: R-1.6.1/src/nmath/d1mach.c 
   if(kCSa) cout << "------TMLMath::D1Mach------" << endl;

   switch (i) {
      case 1:  return DBL_MIN;
      case 2:  return DBL_MAX;
      // = FLT_RADIX  ^ - DBL_MANT_DIG
      // for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON 
      case 3:  return Pow((Double_t)I1Mach(10), -(Double_t)I1Mach(14));
      // = FLT_RADIX  ^ (1- DBL_MANT_DIG) =
      // for IEEE:  = 2^52 = 4503599627370496 = 1/DBL_EPSILON
      case 4:  return Pow((Double_t)I1Mach(10), 1-(Double_t)I1Mach(14));
      case 5:  return log10(2.0);  // = M_LOG10_2 in Rmath.h
      default: return 0.0;
   }//switch
}//D1Mach

//______________________________________________________________________________
Int_t TMLMath::I1Mach(Int_t i)
{
   // Converted from: R-1.6.1/src/nmath/i1mach.c 
   if(kCSa) cout << "------TMLMath::I1Mach------" << endl;

   switch (i) {

      case  1: return 5;
      case  2: return 6;
      case  3: return 0;
      case  4: return 0;

      case  5: return CHAR_BIT * sizeof(Int_t);
      case  6: return sizeof(Int_t)/sizeof(char);

      case  7: return 2;
      case  8: return CHAR_BIT * sizeof(Int_t) - 1;
      case  9: return INT_MAX;

      case 10: return FLT_RADIX;

      case 11: return FLT_MANT_DIG;
      case 12: return FLT_MIN_EXP;
      case 13: return FLT_MAX_EXP;

      case 14: return DBL_MANT_DIG;
      case 15: return DBL_MIN_EXP;
      case 16: return DBL_MAX_EXP;

      default: return 0;
   }//switch
}//I1Mach

//______________________________________________________________________________
Double_t TMLMath::PBetaRaw(Double_t x, Double_t pin, Double_t qin, Int_t lower_tail)
{
   // DESCRIPTION
   // Returns distribution function of the beta distribution.
   // ( = The incomplete beta ratio I_x(p,q) ).
   // 
   // NOTES
   // This routine is a translation into C of a Fortran subroutine
   // by W. Fullerton of Los Alamos Scientific Laboratory.
   // 
   // REFERENCE
   // Bosten and Battiste (1974).
   // Remark on Algorithm 179, CACM 17, p153, (1974).
   // Converted from: R-1.6.1/src/nmath/pbeta.c 
   if(kCSa) cout << "------TMLMath::PBetaRaw------" << endl;

   Double_t ans, c, finsum, p, ps, p1, q, term, xb, xi, y;
   Int_t    n, i, ib, swap_tail;

   const Double_t eps   = .5 * DBL_EPSILON;
   const Double_t sml   = DBL_MIN;
   const Double_t lneps = Log(eps);
   const Double_t lnsml = Log(sml);

// swap tails if x is greater than the mean
   if (pin / (pin + qin) < x) {
      swap_tail = 1;
      y = 1 - x;
      p = qin;
      q = pin;
   } else {
      swap_tail = 0;
      y = x;
      p = pin;
      q = qin;
   }//if

   if ((p + q) * y / (p + 1) < eps) {
   // tail approximation
      xb = p * Log(FMax2(y, sml)) - Log(p) - LBeta(p, q);
      if (xb > lnsml && y != 0) {
         ans = (swap_tail == lower_tail) ? -Expm1(xb) : exp(xb);
      } else {
         ans = (swap_tail == lower_tail) ? 1. : 0;
      }//if
   } else {
   // evaluate the infinite sum first.  term will equal
   // y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i)

   /*___ FIXME ___:  This takes forever (or ends wrongly)
      when (one or) both p & q  are huge */

   ps = q - floor(q);
   if (ps == 0) ps = 1;
   xb = p * Log(y) - LBeta(ps, p) - Log(p);
   ans = 0;
   if (xb >= lnsml) {
      ans = exp(xb);
      term = ans * p;
      if (ps != 1) {
         n = (Int_t)FMax2(lneps/Log(y), 4.0);
         for(i=1 ; i<=n ; i++) {
            xi = i;
            term *= (xi - ps) * y / xi;
            ans += term / (p + xi);
         }//for_i
      }//if
   }//if

// now evaluate the finite sum, maybe.
   if (q > 1) {
      xb = p * Log(y) + q * Log1p(-y) - LBeta(p, q) - Log(q);
      ib = FMax2(xb / lnsml, 0.0);  //Int_t???
      term = exp(xb - ib * lnsml);
      c = 1 / (1 - y);
      p1 = q * c / (p + q - 1);

      finsum = 0;
      n = q;  //Int_t???
      if (q == n) n--;
      for(i=1 ; i<=n ; i++) {
         if (p1 <= 1 && term / eps <= finsum) break;
         xi = i;
         term = (q - xi + 1) * c * term / (p + q - xi);
         if (term > 1) {
            ib--;
            term *= sml;
         }//if
         if (ib == 0)
            finsum += term;
      }//for_i
         ans += finsum;
   }//if

   if (swap_tail == lower_tail)
      ans = 1 - ans;
      ans = FMax2(FMin2(ans, 1.), 0.);
   }//if

   return ans;
}//PBetaRaw

//______________________________________________________________________________
void TMLMath::PNormBoth(Double_t x, Double_t &cum, Double_t &ccum, Int_t i_tail, Int_t log_p)
{
   // i_tail in {0,1,2} means: "lower", "upper", or "both" :
   // if(lower) return  *cum := P[X <= x]
   // if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
   // Converted from: R-1.6.1/src/nmath/pnorm.c 
   if(kCSa) cout << "------TMLMath::PNormBoth------" << endl;

   const Double_t a[5] = {
		2.2352520354606839287,
		161.02823106855587881,
		1067.6894854603709582,
		18154.981253343561249,
		0.065682337918207449113
   };
   const Double_t b[4] = {
		47.20258190468824187,
		976.09855173777669322,
		10260.932208618978205,
		45507.789335026729956
   };
   const Double_t c[9] = {
		0.39894151208813466764,
		8.8831497943883759412,
		93.506656132177855979,
		597.27027639480026226,
		2494.5375852903726711,
		6848.1904505362823326,
		11602.651437647350124,
		9842.7148383839780218,
		1.0765576773720192317e-8
   };
   const Double_t d[8] = {
	   22.266688044328115691,
	   235.38790178262499861,
		1519.377599407554805,
		6485.558298266760755,
		18615.571640885098091,
		34900.952721145977266,
		38912.003286093271411,
		19685.429676859990727
   };
   const Double_t p[6] = {
	   0.21589853405795699,
	   0.1274011611602473639,
	   0.022235277870649807,
	   0.001421619193227893466,
	   2.9112874951168792e-5,
	   0.02307344176494017303
   };
   const Double_t q[5] = {
		1.28426009614491121,
		0.468238212480865118,
		0.0659881378689285515,
		0.00378239633202758244,
		7.29751555083966205e-5
   };

   Double_t xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
   Double_t min = DBL_MIN;
#endif
   Int_t i, lower, upper;

#ifdef IEEE_754
   if(IsNaN(x)) { cum = ccum = x; return; }
#endif

   // Consider changing these :
   eps = DBL_EPSILON * 0.5;

   // i_tail in {0,1,2} =^= {lower, upper, both}
   lower = i_tail != 1;
   upper = i_tail != 0;

   y = fabs(x);
   if (y <= 0.67448975) {
   // qnorm(3/4) = .6744.... -- earlier had 0.66291
	   if (y > eps) {
	      xsq = x * x;
	      xnum = a[4] * xsq;
	      xden = xsq;
	      for (i = 0; i < 3; ++i) {
		      xnum = (xnum + a[i]) * xsq;
		      xden = (xden + b[i]) * xsq;
	      }//for
	   } else {
         xnum = xden = 0.0;
      }//if

      temp = x * (xnum + a[3]) / (xden + b[3]);
      if (lower) cum  = 0.5 + temp;
      if (upper) ccum = 0.5 - temp;
      if (log_p) {
         if (lower) cum  = Log(cum);
         if (upper) ccum = Log(ccum);
      }//if
   } else if (y <= M_SQRT_32) {
   // Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657
      xnum = c[8] * y;
      xden = y;
      for (i = 0; i < 7; ++i) {
         xnum = (xnum + c[i]) * y;
         xden = (xden + d[i]) * y;
      }//for
      temp = (xnum + c[7]) / (xden + d[7]);

#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */
#define do_del(X)							\
	xsq = FTrunc(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if (log_p) {							\
	    cum = (-xsq * xsq * 0.5) + (-del * 0.5) + Log(temp);	\
	    if ((lower && x > 0.) || (upper && x <= 0.))			\
		  ccum = Log1p(-exp(-xsq * xsq * 0.5) * 		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    ccum = 1.0 - cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = cum; if(lower) cum = ccum; ccum = temp;	\
	}

      do_del(y);
      swap_tail;
   } else if ((-37.5193 < x) || (x < 8.2924)) { /* originally had y < 50 */
   // Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 8.29)
   xsq = 1.0 / (x * x);
   xnum = p[5] * xsq;
   xden = xsq;
   for (i = 0; i < 4; ++i) {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
   }//for
   temp = xsq * (xnum + p[4]) / (xden + q[4]);
   temp = (M_1_SQRT_2PI - temp) / y;

   do_del(x);
   swap_tail;
   } else { /* x < -37.5193  OR	8.2924 < x */
      if (log_p) {
      // be better than to just return log(0) or log(1)
         xsq = x*x;
         if( xsq * DBL_EPSILON < 1.)
            del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
         else
            del = 0.;
         cum = -.5*xsq - M_LN_SQRT_2PI - log(y) + Log1p(del);
         ccum = -0.;/*log(1)*/
         swap_tail;
      } else {
         if (x > 0) {cum = 1.; ccum = 0.;}
         else       {cum = 0.; ccum = 1.;}
      }//if
   }//if

#ifdef NO_DENORMS
   /* do not return "denormalized" -- needed ?? */
   if (log_p) {
      if (cum > -min)  cum  = -0.;
      if (ccum > -min) ccum = -0.;
   } else {
      if (cum < min)  cum  = 0.;
      if (ccum < min) ccum = 0.;
   }//if
#endif
   return;
}//PNormBoth


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Mathlib utility functions                                            //
//                                                                      //
// Code taken from file R-1.6.1/src/nmath/mlutils.c                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#ifndef IEEE_754

void ml_error(int n)
{
   switch (n) {

   case ME_NONE:
      errno = 0;
      break;

   case ME_DOMAIN:
   case ME_NOCONV:
      errno = EDOM;
      break;

   case ME_RANGE:
      errno = ERANGE;
   break;

   default:
      break;
   }//switch
}//ml_error

#endif



