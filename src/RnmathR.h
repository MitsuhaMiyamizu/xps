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
 *  header file R-1.6.1/src/nmath/nmath.h
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

/* Private header file for use during compilation of Mathlib */
#ifndef RNMATHR_H
#define RNMATHR_H

//include to define IEEE_754, WORDS_BIGENDIAN
#include "RconfigR.h"

#include <stdio.h>
# define MATHLIB_ERROR(fmt,x)	{ printf(fmt,x); exit(1); }
# define MATHLIB_WARNING(fmt,x)		printf(fmt,x)
# define MATHLIB_WARNING2(fmt,x,x2)	printf(fmt,x,x2)
# define MATHLIB_WARNING3(fmt,x,x2,x3)	printf(fmt,x,x2,x3)
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) printf(fmt,x,x2,x3,x4)


#ifdef IEEE_754
#define ML_POSINF	(1.0 / 0.0)
#define ML_NEGINF	((-1.0) / 0.0)
#define ML_NAN		(0.0 / 0.0)
#else
#define ML_POSINF	DBL_MAX
#define ML_NEGINF	(-DBL_MAX)
#define ML_NAN		(-DBL_MAX*(1-1e-15))
#endif


#ifdef  __cplusplus
extern "C" {
#endif

#ifdef IEEE_754
#define ML_ERROR(x)	/* nothing */
#define ML_UNDERFLOW	(DBL_MIN * DBL_MIN)
#define ML_VALID(x)	(!IsNaN(x))
//#define ML_VALID(x)	(!ISNAN(x))
#else /*--- NO IEEE: No +/-Inf, NAN,... ---*/
void ml_error(int n);
//extern  void ml_error(int n);
#define ML_ERROR(x)	ml_error(x)
#define ML_UNDERFLOW	0
#define ML_VALID(x)	(errno == 0)
#endif

#ifdef  __cplusplus
}
#endif

#define ME_NONE		0
/*	no error */
#define ME_DOMAIN	1
/*	argument out of domain */
#define ME_RANGE	2
/*	value out of range */
#define ME_NOCONV	4
/*	process did not converge */
#define ME_PRECISION	8
/*	does not have "full" precision */
#define ME_UNDERFLOW	16
/*	and underflow occured (important for IEEE)*/

//#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN); return ML_NAN; }

/* Wilcoxon Rank Sum Distribution */

#define WILCOX_MAX 50

/* Wilcoxon Signed Rank Distribution */

#define SIGNRANK_MAX 50


#endif /* RNMATHR_H */
