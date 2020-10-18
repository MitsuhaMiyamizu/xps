// Author: Christian Stratowa 25/11/2002             last modified: 26/12/2002

/*
*******************************************************************************
***********************  Statistics Package for ROOT  *************************
*******************************************************************************
 *
 *  Copyright (C) 2000-2007 Dr. Christian Stratowa
 *
 *  Adapted from R by: Christian Stratowa, Vienna, Austria <cstrato@aon.at>
 *  
 *  Algorithms for statistical functions adapted for ROOT from:
 *  header file R-1.6.1/src/nmath/dpq.h
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2001  The R Development Core Team
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
 */

	/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p
							/* "DEFAULT" */
							/* --------- */
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */

#define R_D_Lval(p)	(lower_tail ? (p) : (1 - (p)))	/*  p  */
#define R_D_Cval(p)	(lower_tail ? (1 - (p)) : (p))	/*  1 - p */

//????????? log, exp, expm1, sqrt, floor, fabs
//#define R_D_val(x)	(log_p	? Log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_D_qIv(p)	(log_p	? exp(p) : (p))		/*  p  in qF(p,..) */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))	/* exp(x) */
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */

#define R_DT_val(x)	R_D_val(R_D_Lval(x))		/*  x  in pF */
#define R_DT_Cval(x)	R_D_val(R_D_Cval(x))		/*  1 - x */
/*#define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))		 *  p  in qF ! */
#define R_DT_qIv(p)	(log_p ? (lower_tail ? exp(p) : - Expm1(p)) \
			       : R_D_Lval(p))

/*#define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))		 *  1 - p in qF */
#define R_DT_CIv(p)	(log_p ? (lower_tail ? -Expm1(p) : exp(p)) \
			       : R_D_Cval(p))


#define R_DT_exp(x)	R_D_exp(R_D_Lval(x))		/* exp(x) */
#define R_DT_Cexp(x)	R_D_exp(R_D_Cval(x))		/* exp(1 - x) */

#define R_DT_log(p)	(lower_tail ? R_D_log(p) :		\
			 log1p(- (log_p ? exp(p) : p)))/* log(p)	in qF */

#define R_DT_Clog(p)	(lower_tail ?				\
			 Log1p(- (log_p ? exp(p) : p)) :	\
			 R_D_log(p))			/* log(1 - p)	in qF */

#define R_Q_P01_check(p)			\
    if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
   {ML_ERROR(ME_DOMAIN); return ML_NAN;}
//	ML_ERR_return_NAN


/* additions for density functions (C.Loader) */
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_D_nonint(x) 	  (fabs((x) - floor((x)+0.5)) > 1e-7)
#define R_D_notnnegint(x) (x < 0. || R_D_nonint(x))

#define R_D_nonint_check(x) 				\
   if(R_D_nonint(x)) {					\
	MATHLIB_WARNING("non-integer x = %f", x);	\
	return R_D__0;					\
   }
