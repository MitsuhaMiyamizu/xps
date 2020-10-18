// Author: Christian Stratowa 11/25/2002             last modified: 04/11/2009

/*
 *******************************************************************************
 ***********************  Statistics Package for ROOT  *************************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2009 Dr. Christian Stratowa
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

#ifndef __StatUtils__
#define __StatUtils__

// ROOT includes
#include "TNamed.h"
#include "TTree.h"

extern Int_t NumSep(const char *name, const char *sep);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TUnivariateTest                                                      //
//                                                                      //
// Base class for univariate tests                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TUnivariateTest: public TNamed {

   protected:
      Double_t   fConfLevel;    //confidence level (<0 - do not calculate)
      Double_t   fConfHi;       //upper value of confidence interval
      Double_t   fConfLo;       //lower value of confidence interval
      Double_t   fMu;           //true value of the mean
      Double_t   fMean1;        //mean value for group1
      Double_t   fMean2;        //mean value for group2
      Double_t   fStat;         //value of the statistic used, e.g. Z, T
      Double_t   fPValue;       //p-value
      Double_t   fPChance;      //p-value of random sampling
      Int_t      fNPerm;        //number of permutations for p-chance
      Double_t   fDF;           //degrees of freedom
      Double_t   fSE;           //standard error s.e.
      TString    fAlternative;  //alternative hypothesis
      TString    fAdjustment;   //method used to adjust p-value, e.g. Bonferroni
      Double_t   fNA;           //value to be used as NA
      Bool_t     fTwoSample;    //TRUE if two-sample test
      Bool_t     fPaired;       //TRUE if paired test
      Bool_t     fHasNA;        //TRUE if data have missing values
      Bool_t     fAdjPVal;      //adjust to p-value(TRUE) or p-chance(FALSE)
      
   protected:
      
   public:
      TUnivariateTest();
      TUnivariateTest(const char *name, const char *type = "uvt");
      TUnivariateTest(const char *name, const char *type, Double_t na);
      virtual ~TUnivariateTest();

      void     Init(Int_t nperm = 0, Double_t mu = 0, Bool_t paired = kFALSE,
                  Double_t conflevel = 0.95, const char *alt = "twosided");

      virtual void     PrintInfo();

      virtual Double_t PValue(Double_t stat, Double_t df);
      virtual Double_t PValue(Double_t stat, Double_t df, Double_t se,
                          Double_t &lo, Double_t &hi);

      virtual Double_t Statistic(Int_t n, Double_t *grp, Double_t &mean,
                          Double_t &se, Double_t &df, Double_t mu = 0);
      virtual Double_t Statistic(Int_t n, Double_t *grp, Double_t &mean,
                          Double_t &se, Double_t &df, Double_t mu, Double_t na);
      virtual Double_t Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                          Double_t *grp2, Double_t &mean1, Double_t &mean2,
                          Double_t &se, Double_t &df, Double_t mu = 0);
      virtual Double_t Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                          Double_t *grp2, Double_t &mean1, Double_t &mean2,
                          Double_t &se, Double_t &df, Double_t mu, Double_t na);

      void     Test(Int_t n, Double_t *grp, Double_t mu = 0);
      void     Test(Int_t n1, Double_t *grp1, Int_t n2, Double_t *grp2,
                  Int_t nperm = -1, Double_t mu = 0);
      Int_t    Test(const char *infile, const char *outfile, const char *varlist,
                  Int_t nrows, Int_t nperm = -1, Double_t mu = 0,
                  const char *sepi = "\t", const char *sepo = "\t", char delim = '\n');
      Int_t    Test(Int_t n, Int_t *gid, TTree *intree, const char *leafname,
                  TTree *outtree, const char *varlist = "*", Int_t nperm = -1,
                  Double_t mu = 0);
      Int_t    Test(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                  TTree *outtree, const char *varlist = "*", Int_t nperm = -1,
                  Double_t mu = 0);

      Double_t  PChance(Int_t n1, Double_t *grp1, Int_t n2, Double_t *grp2,
                   Int_t nperm, Double_t tstat);
      Double_t *PChance(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                   Double_t *pcha, Int_t *nperm, Double_t *tstat);
      Double_t  Permute(Int_t n, Double_t *arr, Int_t n1, Int_t nperm,
                   Double_t tstat);
      Double_t *Permute(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                   Double_t *pcha, Int_t *nperm, Double_t *tstat);
      Double_t  Sample(Int_t n, Double_t *arr, Int_t n1, Int_t nperm,
                   Double_t tstat);
      Double_t *Sample(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                   Double_t *pcha, Int_t *nperm, Double_t *tstat);

      Double_t *PAdjustWY(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                   Double_t *pcha, Int_t *nperm, Double_t *tstat, Double_t *padj);
      Double_t *PermuteWY(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                   Double_t *pcha, Int_t *nperm, Double_t *tstat, Double_t *padj);
      Double_t *SampleWY(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                   Double_t *pcha, Int_t *nperm, Double_t *tstat, Double_t *padj);

      Double_t *PAdjust(Int_t n, Double_t *pval, Double_t *padj);
      Double_t *Bonferroni(Int_t n, Double_t *pval, Double_t *padj);
      Double_t *BY(Int_t n, Double_t *pval, Double_t *padj);
      Double_t *FDR(Int_t n, Double_t *pval, Double_t *padj);
      Double_t *Hochberg(Int_t n, Double_t *pval, Double_t *padj);
      Double_t *Holm(Int_t n, Double_t *pval, Double_t *padj);

      void     SetAdjustment(const char *adj, Bool_t adjpval = kTRUE)
                                {fAdjustment = adj; fAdjPVal = adjpval;}
      void     SetAlternative(const char *alt) {fAlternative = alt;}
      void     SetConfidenceLevel(Double_t cl) {fConfLevel   = cl;}
      void     SetIsPaired(Bool_t paired)      {fPaired      = paired;}
      void     SetNumPerm(Int_t nperm)         {fNPerm       = nperm;}
      void     SetTwoSample(Bool_t twos)       {fTwoSample   = twos;}
      void     SetNA(Double_t na)              {fNA = na; fHasNA = kTRUE;}

      TString  GetAdjustment()       const {return fAdjustment;}
      TString  GetAlternative()      const {return fAlternative;}
      Double_t GetConfidenceLevel()  const {return fConfLevel;}
      Bool_t   GetIsPaired()         const {return fPaired;}
      Double_t GetMu()               const {return fMu;}
      Bool_t   GetTwoSample()        const {return fTwoSample;}

      Double_t GetStatistic()        const {return fStat;}
      Double_t GetDegFreedom()       const {return fDF;}
      Double_t GetStdError()         const {return fSE;}
      Double_t GetPValue()           const {return fPValue;}
      Double_t GetPChance()          const {return fPChance;}
      Int_t    GetNumPerm()          const {return fNPerm;}

      ClassDef(TUnivariateTest,1) //UnivariateTest
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStudentTest                                                         //
//                                                                      //
// Class for Student's t-test                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TStudentTest: public TUnivariateTest {

   protected:
      Bool_t     fEqualVar;     //TRUE if variances treated as equal
      
   protected:
      
   public:
      TStudentTest();
      TStudentTest(const char *name, const char *type = "stt");
      TStudentTest(const char *name, const char *type, Double_t na);
      virtual ~TStudentTest();

      void     Init(Int_t nperm = 0, Double_t mu = 0, Bool_t paired = kFALSE,
                  Double_t conflevel = 0.95, Bool_t varequ = kFALSE,
                  const char *alt = "twosided");

      virtual void     PrintInfo();

      virtual Double_t PValue(Double_t stat, Double_t df);
      virtual Double_t PValue(Double_t stat, Double_t df, Double_t se,
                          Double_t &lo, Double_t &hi);

      virtual Double_t Statistic(Int_t n, Double_t *grp, Double_t &mean,
                          Double_t &se, Double_t &df, Double_t mu = 0);
      virtual Double_t Statistic(Int_t n, Double_t *grp, Double_t &mean,
                          Double_t &se, Double_t &df, Double_t mu, Double_t na);
      virtual Double_t Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                          Double_t *grp2, Double_t &mean1, Double_t &mean2,
                          Double_t &se, Double_t &df, Double_t mu = 0);
      virtual Double_t Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                          Double_t *grp2, Double_t &mean1, Double_t &mean2,
                          Double_t &se, Double_t &df, Double_t mu, Double_t na);

      void     SetEqualVariance(Bool_t var) {fEqualVar = var;}
      Bool_t   GetEqualVariance()     const {return fEqualVar;}

      ClassDef(TStudentTest,1) //StudentTest
};


#endif
