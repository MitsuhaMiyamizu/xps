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

using namespace std;

#include <new>  //needed for new (nothrow)

// ROOT
#include "TMath.h"
#include "TRandom.h"
#include "TList.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TFriendElement.h"

#include "StatUtils.h"
#include "TMLMath.h"
#include "TStat.h"

//debug: print class method names
const Bool_t  kCS  = 0;
const Bool_t  kCSa = 0;

ClassImp(TUnivariateTest);
ClassImp(TStudentTest);

const Int_t  kBufSize = 1024;    //bufsize for chars
const Int_t  kLineBuf = 16635;   //to read lines from infile
const char*  kGroup   = "Group"; //name for infile 2.row/1.col 


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TUnivariateTest                                                      //
//                                                                      //
// Base class for univariate tests                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TUnivariateTest::TUnivariateTest()
                :TNamed()
{
   // Default UnivariateTest constructor
   if(kCS) cout << "---TUnivariateTest::TUnivariateTest(default)------" << endl;

   fConfHi      = R_PosInf;
   fConfLo      = R_NegInf;
   fMean1       = NA_REAL;
   fMean2       = NA_REAL;
   fStat        = NA_REAL;
   fPValue      = -1;
   fPChance     = -1;
   fDF          = 0;
   fSE          = 0;
   fNA          = -1;
   fTwoSample   = kTRUE;
   fHasNA       = kFALSE;
   fAdjPVal     = kFALSE;
   fAlternative = "twosided";
   fAdjustment  = "none";

   this->Init();
}//Constructor

//______________________________________________________________________________
TUnivariateTest::TUnivariateTest(const char *name, const char *type)
                :TNamed(name, type)
{
   // Normal UnivariateTest constructor
   if(kCS) cout << "---TUnivariateTest::TUnivariateTest------" << endl;

   fConfHi      = R_PosInf;
   fConfLo      = R_NegInf;
   fMean1       = NA_REAL;
   fMean2       = NA_REAL;
   fStat        = NA_REAL;
   fPValue      = -1;
   fPChance     = -1;
   fDF          = 0;
   fSE          = 0;
   fNA          = -1;
   fTwoSample   = kTRUE;
   fHasNA       = kFALSE;
   fAdjPVal      = kFALSE;
   fAlternative = "twosided";
   fAdjustment  = "none";

   this->Init();
}//Constructor

//______________________________________________________________________________
TUnivariateTest::TUnivariateTest(const char *name, const char *type, Double_t na)
                :TNamed(name, type)
{
   // Normal UnivariateTest constructor
   // na is the value used in the data set to indicate missing values
   if(kCS) cout << "---TUnivariateTest::TUnivariateTest(na)------" << endl;

   fConfHi      = R_PosInf;
   fConfLo      = R_NegInf;
   fMean1       = NA_REAL;
   fMean2       = NA_REAL;
   fStat        = NA_REAL;
   fPValue      = na;
   fPChance     = na;
   fDF          = 0;
   fSE          = 0;
   fNA          = na;
   fTwoSample   = kTRUE;
   fHasNA       = kTRUE;
   fAlternative = "twosided";
   fAdjustment  = "none";

   this->Init();
}//Constructor

//______________________________________________________________________________
TUnivariateTest::~TUnivariateTest()
{
   // UnivariateTest destructor
   if(kCS) cout << "---TUnivariateTest::~TUnivariateTest------" << endl;

}//Destructor

//______________________________________________________________________________
void TUnivariateTest::Init(Int_t nperm, Double_t mu, Bool_t paired,
                      Double_t conflevel, const char *alt)
{
   // Initialize univariate test
   // not necessesary if default values are used
   // if conflevel < 0 then do not calculate confidence interval
   if(kCS) cout << "---TUnivariateTest::Init------" << endl;

   fNPerm  = nperm;
   fMu     = mu;
   fPaired = paired;

   if (conflevel <= 1.0) {
      fConfLevel = conflevel;
   } else {
      cout << "Warning: Confidence level >1, thus set to default 0.95" << endl;
      fConfLevel = 0.95;
   }//if

   fAlternative = alt;
}//Init

//______________________________________________________________________________
void TUnivariateTest::PrintInfo()
{
   // Print results of test
   if(kCS) cout << "------TUnivariateTest::PrintInfo------" << endl;

   cout << "==============================================================================" << endl;
   cout << endl;
   if (!fTwoSample) {
      cout << "         One Sample normal-test" << endl;
   } else if (fPaired) {
      cout << "         Paired normal-test" << endl;
   } else {
      cout << "         Two Sample normal-test" << endl;
   }//if
   cout << endl;

   cout << "z  = " << fStat << endl;
   cout << "p-value  = " << fPValue << endl;

   if (fNPerm > 0) {
      cout << "numperm  = " << fNPerm << endl;
      cout << "p-chance = " << fPChance << endl;
   }//if

   cout << "alternative hypothesis: true ";
   if (fPaired || fTwoSample) cout << "difference in means ";
   else                       cout << "mean ";
   if (strcmp(fAlternative, "greater") == 0)   cout << "is greater than ";
   else if (strcmp(fAlternative, "less") == 0) cout << "is less than ";
   else                                        cout << "is not equal to ";
   cout << fMu << endl;

   if (fConfLevel >= 0) {
      cout << 100*fConfLevel << " percent confidence interval:" << endl;
      cout << " [ " << fConfLo << " ,  " << fConfHi << " ]" << endl;
   }//if

   cout << "sample estimates: " << endl;
   cout << "mean(grp1)";
   if (fPaired || !fTwoSample) cout << endl;
   else                        cout << "      mean(grp2)" << endl;
   cout << "   " << fMean1;
   if (fPaired || !fTwoSample) cout << endl;
   else                        cout << "           " << fMean2 << endl;
   cout << endl;

   cout << "==============================================================================" << endl;
}//PrintInfo

//______________________________________________________________________________
Double_t TUnivariateTest::PValue(Double_t stat, Double_t df)
{
   // Calculate p-value only
   if(kCSa) cout << "------TUnivariateTest::PValue------" << endl;

   if (TMLMath::IsNaN(stat)) return NA_REAL;

   Double_t pval = 1;
   if (strcmp(fAlternative, "twosided") == 0) {
      pval = 2*TMLMath::PNorm(-TMath::Abs(stat), 0, 1);
   } else if (strcmp(fAlternative, "greater") == 0) {
      pval = TMLMath::PNorm(stat, 0, 1, kFALSE);
   } else if (strcmp(fAlternative, "less") == 0) {
      pval = TMLMath::PNorm(stat, 0 ,1, kTRUE);
   } else {
      cerr << "Error: Alternative not known" << endl;
      return NA_REAL;
   }//if

   return pval;
}//PValue

//______________________________________________________________________________
Double_t TUnivariateTest::PValue(Double_t stat, Double_t df, Double_t se,
                          Double_t &lo, Double_t &hi)
{
   // Calculate p-value and confidence interval [lo, hi]
   if(kCSa) cout << "------TUnivariateTest::PValue(cf)------" << endl;

   if (TMLMath::IsNaN(stat)) return NA_REAL;

   Double_t pval = 1;
   Double_t cint = 0;
   if (strcmp(fAlternative, "twosided") == 0) {
      pval = 2*TMLMath::PNorm(-TMath::Abs(stat), 0, 1);
      if (fConfLevel >= 0) {
         cint = TMLMath::QNorm(0.5 + fConfLevel/2.);
         hi   = fMu + (stat + cint)*se;
         lo   = fMu + (stat - cint)*se;
      }//if
   } else if (strcmp(fAlternative, "greater") == 0) {
      pval = TMLMath::PNorm(stat, 0, 1, kFALSE);
      if (fConfLevel >= 0) {
         cint = TMLMath::QNorm(fConfLevel);
         hi   = R_PosInf;
         lo   = fMu + (stat - cint)*se;
      }//if
   } else if (strcmp(fAlternative, "less") == 0) {
      pval = TMLMath::PNorm(stat, 0 ,1, kTRUE);
      if (fConfLevel >= 0) {
         cint = TMLMath::QNorm(fConfLevel);
         hi   = fMu + (stat + cint)*se;
         lo   = R_NegInf;
      }//if
   } else {
      cerr << "Error: Alternative not known" << endl;
      return NA_REAL;
   }//if

   return pval;
}//PValue

//______________________________________________________________________________
Double_t TUnivariateTest::Statistic(Int_t n, Double_t *grp, Double_t &mean,
                          Double_t &se, Double_t &df, Double_t mu)
{
   // Return z-statistic for one-sample test
   if(kCSa) cout << "------TUnivariateTest::Statistic(1)------" << endl;

   if (n < 2) {
      cerr << "Error: Less than two values in group" << endl;
      return NA_REAL;
   }//if

   Double_t mn, var, dfr, err, zs;
   mn  = TStat::Mean(n, grp);
   var = TStat::Var(n, grp, mn);
   dfr = n - 1;
   err = TMath::Sqrt(var/n);
   zs  = (mn - mu)/err;

   mean = mn;
   se   = err;
   df   = dfr;
   return zs;
}//Statistic

//______________________________________________________________________________
Double_t TUnivariateTest::Statistic(Int_t n, Double_t *grp, Double_t &mean,
                          Double_t &se, Double_t &df, Double_t mu, Double_t na)
{
   // Return z-statistic for one-sample test 
   // missing values are allowed but have to be set to a chosen value na
   if(kCSa) cout << "------TUnivariateTest::Statistic(1na)------" << endl;

   if (n < 2) {
      cerr << "Error: Less than two values in group" << endl;
      return NA_REAL;
   }//if

   Int_t len = n;  //length of array w/o na
   Double_t mn, var, dfr, err, zs;
   mn  = TStat::Mean(n, grp, len, na);
   var = TStat::Var(n, grp, mn, len, na);
   dfr = len - 1;
   err = TMath::Sqrt(var/len);
   zs  = (mn - mu)/err;

   mean = mn;
   se   = err;
   df   = dfr;
   return zs;
}//Statistic

//______________________________________________________________________________
Double_t TUnivariateTest::Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                          Double_t *grp2, Double_t &mean1, Double_t &mean2,
                          Double_t &se, Double_t &df, Double_t mu)
{
   // Return z-statistic for two-sample test
   if(kCSa) cout << "------TUnivariateTest::Statistic(2)------" << endl;

   Double_t mn1, mn2, err, dfr, zs;
   mn1 = err = dfr = zs = 0;
   mn2 = NA_REAL;
   if (fPaired) {
      if (n1 != n2) {
         cerr << "Error: Group1 and group2 must have paired values" << endl;
         return NA_REAL;
      }//if

      // paired grp = grp1 - grp2
      Double_t *grp = 0;
      if (!(grp = new (nothrow) Double_t[n1])) {
         cerr << "Error: Could not initialize memory!" << endl;
         return NA_REAL;
      }//if
      for (Int_t i=0; i<n1; i++) grp[i] = grp1[i] - grp2[i];

      // one-sample test
      zs = this->Statistic(n1, grp, mn1, err, dfr, mu);

      if (grp) delete [] grp;
   } else {
      if ((n1 < 2) || (n2 < 2)) {
         cerr << "Error: Less than two values in one of the groups" << endl;
         return NA_REAL;
      }//if

      Double_t var1, var2;
      mn1  = TStat::Mean(n1, grp1);
      mn2  = TStat::Mean(n2, grp2);
      var1 = TStat::Var(n1, grp1, mn1);
      var2 = TStat::Var(n2, grp2, mn2);

      dfr = n1 + n2 - 2;
      err = TMath::Sqrt(var1/n1 + var2/n2);
      zs  = (mn1 - mn2 - mu)/err;
   }//if

   mean1 = mn1;
   mean2 = mn2;
   se    = err;
   df    = dfr;
   return zs;
}//Statistic

//______________________________________________________________________________
Double_t TUnivariateTest::Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                          Double_t *grp2, Double_t &mean1, Double_t &mean2,
                          Double_t &se, Double_t &df, Double_t mu, Double_t na)
{
   // Return z-statistic for two-sample test 
   // missing values are allowed but have to be set to a chosen value na
   if(kCSa) cout << "------TUnivariateTest::Statistic(2na)------" << endl;

   Double_t mn1, mn2, err, dfr, zs;
   mn1 = err = dfr = zs = 0;
   mn2 = NA_REAL;
   if (fPaired) {
      if (n1 != n2) {
         cerr << "Error: Group1 and group2 must have paired values" << endl;
         return NA_REAL;
      }//if

      // paired grp = grp1 - grp2
      Double_t *grp = 0;
      if (!(grp = new (nothrow) Double_t[n1])) {
         cerr << "Error: Could not initialize memory!" << endl;
         return NA_REAL;
      }//if
      // use only complete cases
      Int_t m = n1;
      for (Int_t i=0; i<n1; i++) {
         if ((grp1[i] != na) && (grp2[i] != na)) grp[i] = grp1[i] - grp2[i];
         else  m--;
      }//for

      // one-sample test
      zs = this->Statistic(m, grp, mn1, err, dfr, mu);

      if (grp) delete [] grp;
   } else {
      if ((n1 < 2) || (n2 < 2)) {
         cerr << "Error: Less than two values in one of the groups" << endl;
         return NA_REAL;
      }//if

      Int_t len1 = n1;  //length of grp1 w/o na
      Int_t len2 = n2;  //length of grp2 w/o na
      Double_t var1, var2;
      mn1  = TStat::Mean(n1, grp1, len1, na);
      mn2  = TStat::Mean(n2, grp2, len2, na);
      var1 = TStat::Var(n1, grp1, mn1, len1, na);
      var2 = TStat::Var(n2, grp2, mn2, len2, na);

      if ((len1 < 2) || (len2 < 2)) {
         // prevent lots of error messages during random sampling
         if (fNPerm <= 0) {
            cerr << "Error: Less than 2 non-missing values in one of the groups"
                 << endl;
         }//if
         return NA_REAL;
      }//if

      dfr = len1 + len2 - 2;
      err = TMath::Sqrt(var1/len1 + var2/len2);
      zs  = (mn1 - mn2 - mu)/err;
   }//if

   mean1 = mn1;
   mean2 = mn2;
   se    = err;
   df    = dfr;
   return zs;
}//Statistic

//______________________________________________________________________________
void TUnivariateTest::Test(Int_t n, Double_t *grp, Double_t mu)
{
   // Do one-sample test
   if(kCS) cout << "------TUnivariateTest::Test(1)------" << endl;

   if (fHasNA) fStat = this->Statistic(n, grp, fMean1, fSE, fDF, mu, fNA);
   else        fStat = this->Statistic(n, grp, fMean1, fSE, fDF, mu);

   if (fConfLevel < 0) fPValue = this->PValue(fStat, fDF);
   else                fPValue = this->PValue(fStat, fDF, fSE, fConfLo, fConfHi);

   fMu = mu;
   fTwoSample = kFALSE;
}//Test

//______________________________________________________________________________
void TUnivariateTest::Test(Int_t n1, Double_t *grp1, Int_t n2, Double_t *grp2,
                      Int_t nperm, Double_t mu)
{
   // Do two-sample test
   if(kCS) cout << "------TUnivariateTest::Test(2)------" << endl;

   if (nperm >= 0) fNPerm = nperm;
   if (mu != 0)    fMu    = mu;
//   if (mu != fMu)  fMu    = mu;

// Check size of groups
   if ((n1 < 2) || (n2 < 2)) {
      cerr << "Error: Less than two values in one of the groups" << endl;
      return;
   }//if

   if (fHasNA) {
      fStat = this->Statistic(n1, grp1, n2, grp2, fMean1, fMean2, fSE, fDF, mu, fNA);
   } else {
      fStat = this->Statistic(n1, grp1, n2, grp2, fMean1, fMean2, fSE, fDF, mu);
   }//if

   if (fConfLevel < 0) fPValue = this->PValue(fStat, fDF);
   else                fPValue = this->PValue(fStat, fDF, fSE, fConfLo, fConfHi);

   if (fNPerm > 0) {
      fPChance = this->PChance(n1, grp1, n2, grp2, fNPerm, fStat);
   }//if

   fTwoSample = kTRUE;
}//Test

//______________________________________________________________________________
Int_t TUnivariateTest::Test(const char *infile, const char *outfile,
                       const char *varlist, Int_t nrows, Int_t nperm,
                       Double_t mu, const char *sepi, const char *sepo,
                       char delim)
{
   // Do multiple univariate tests for infile and export result as outfile.
   // infile must contain:  1st colum: row name or identifier
   //    1st row:   header row
   //    2nd row:   "Group" in 1st column and group id 1 or 2
   //               two-sample test: arbitrary distribution of groups allowed
   //               paired test of n pairs: pair is in col[i] and col[n+i]
   //                                   or: pair is in col[i] and col[i+1]
   //               use group id <=0 if you want to exclude columns
   //               (for paired test corresponding column pairs must be excluded)
   //               if 2nd row is missing, a one-sample test is done
   // varlist is the list of variables to calculate and to export:
   //    "stat" - univariate statistic
   //    "mn1"  - mean value of group1
   //    "mn2"  - mean value of group2  (not used for one-sample test)
   //    "se"   - standard error
   //    "df"   - degrees of freedom
   //    "pval" - p-value
   //    "nper" - number of permutations used to determine p-chance
   //    "pcha" - p-chance
   //    "padj" - p-value or p-chance adjusted for multiple comparisons
   //    e.g. for usual test, varlist = "stat:mn1:mn2:se:df:pval"
   //    "*" is equal to "stat:mn1:mn2:se:df:pval:nper:pcha:padj"
   // nrows is number of data rows (minus header and group)
   // nperm is number of permutations to be done
   if(kCS) cout << "------TUnivariateTest::Test------" << endl;

   if (nperm >= 0) fNPerm = nperm;
   else            nperm  = fNPerm;
   if (mu != 0)    fMu    = mu;
//   if (mu != fMu)  fMu    = mu;

   Int_t err = 0;
   Int_t idx = 0;
   Int_t nid = 0;
   Int_t n1  = 0;
   Int_t n2  = 0;

// Check for two-sample test
   Bool_t ok  = (fTwoSample && nperm > 0);  //permute two-sample/paired sample
   Bool_t ok2 = (fTwoSample && !fPaired);   //two-sample only

// Check for Westfall-Young adjustment
   Bool_t wy = (strcmp(fAdjustment, "wy") == 0);

// Open input file and output file
   ifstream input(infile, ios::in);
   if (!input.good()) {
      cerr << "Error: Could not create input <" << infile << ">" << endl;
      return 1;
   }//if
   ofstream output(outfile, ios::out);
   if (!output) {
      cerr << "Error: Could not create output <" << outfile << ">" << endl;
      return 1;
   }//if

// Init local arrays
   Int_t    *gid  = 0;
   Double_t *grp1 = 0;
   Double_t *grp2 = 0;

   TString  *arrName = 0;  //row names
   Double_t *arrStat = 0;  //statistic
   Double_t *arrMn1  = 0;  //mean of grp1
   Double_t *arrMn2  = 0;  //mean of grp2
   Double_t *arrSE   = 0;  //standard error
   Double_t *arrDF   = 0;  //degree of freedom
   Double_t *arrPVal = 0;  //p-value
   Int_t    *arrNPer = 0;  //number of permutations
   Double_t *arrPCha = 0;  //p-chance
   Double_t *arrPAdj = 0;  //p-adjusted
   Double_t **table  = 0;  //optional table to store infile data

// Decompose varlist
   Bool_t hasStat = kFALSE;
   Bool_t hasMn1  = kFALSE;
   Bool_t hasMn2  = kFALSE;
   Bool_t hasSE   = kFALSE;
   Bool_t hasDF   = kFALSE;
   Bool_t hasPVal = kFALSE;
   Bool_t hasNPer = kFALSE;
   Bool_t hasPCha = kFALSE;
   Bool_t hasPAdj = kFALSE;

   if (strcmp(varlist,"*") == 0) {
      hasStat = kTRUE;
      hasMn1  = kTRUE;
      hasMn2  = kTRUE;
      hasSE   = kTRUE;
      hasDF   = kTRUE;
      hasPVal = kTRUE;
      hasNPer = kTRUE;
      hasPCha = kTRUE;
      hasPAdj = kTRUE;
   } else {
      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;
      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname,"stat") == 0) {hasStat = kTRUE;}
         if (strcmp(vname,"mn1")  == 0) {hasMn1  = kTRUE;}
         if (strcmp(vname,"mn2")  == 0) {hasMn2  = kTRUE;}
         if (strcmp(vname,"se")   == 0) {hasSE   = kTRUE;}
         if (strcmp(vname,"df")   == 0) {hasDF   = kTRUE;}
         if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE;}
         if (strcmp(vname,"nper") == 0) {hasNPer = kTRUE;}
         if (strcmp(vname,"pcha") == 0) {hasPCha = kTRUE;}
         if (strcmp(vname,"padj") == 0) {hasPAdj = kTRUE;}
         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while
      delete [] dname;
   }//if

// Read header to get number of data columns 
   char nextline[kLineBuf];
   input.getline(nextline, kLineBuf, delim);
   if (!input.good()) return 1;
   Int_t numdata = NumSep(nextline, sepi); //+ 1; for total #cols

// Read group indices from second row
   input.getline(nextline, kLineBuf, delim);
   if (!input.good()) return 1;

//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
// Init group id
   if (!(gid = new (nothrow) Int_t[numdata])) {err = 1; goto cleanup;}
   for (Int_t i=0; i<numdata; i++) gid[i] = 0;

   // check for "Group" identifier kGroup in first column
   Bool_t hasRow2;
   char   name[kBufSize];
   strtok(strcpy(name,nextline),sepi);
   if (strcmp(name,kGroup)  == 0) {
      // read group indices
      for (Int_t i=0; i<numdata; i++) {
         gid[i] = atoi(strtok(NULL, sepi));
      }//for_i

//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
      // create arrays grp1 and grp2
      for (Int_t i=0; i<numdata; i++) {
         if (gid[i] == 1) {n1++; nid++;}
         if (gid[i] == 2) {n2++; nid++;}
      }//if
      if (!(grp1 = new (nothrow) Double_t[n1])) {err = 1; goto cleanup;} 
      if (!(grp2 = new (nothrow) Double_t[n2])) {err = 1; goto cleanup;} 

      // check if one/two-sample test
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
      fTwoSample = (n1 < nid) ? kTRUE : kFALSE;
      hasRow2    = kTRUE;
   } else {
      cout << "Warning: Group identifiers missing; assuming one-sample test."
           << endl;

      // create array grp1
      nid = numdata;
      n1  = numdata;
      if (!(grp1 = new (nothrow) Double_t[n1])) {err = 1; goto cleanup;} 
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
      for (Int_t i=0; i<numdata; i++) gid[i] = 1;

      fTwoSample = kFALSE;
      hasRow2    = kFALSE;
   }//if

// Check size of groups
   if (nid < 2) {
      cerr << "Error: Less than two data columns selected" << endl;
      return 1;
   }//if

   ok  = (fTwoSample && nperm > 0);  //permute two-sample/paired sample
   ok2 = (fTwoSample && !fPaired);   //two-sample only

// Create local arrays
   if (!(arrName = new (nothrow) TString[nrows]))  {err = 1; goto cleanup;}
   if (!(arrStat = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrMn1  = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrMn2  = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrSE   = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrDF   = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrPVal = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrNPer = new (nothrow) Int_t[nrows]))    {err = 1; goto cleanup;}
   if (!(arrPCha = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrPAdj = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}

// Init local arrays
   for (Int_t i=0; i<nrows; i++) {
      arrStat[i] = 1;
      arrMn1[i]  = fMu;
      arrMn2[i]  = fMu;
      arrSE[i]   = 0;
      arrDF[i]   = 0;
      arrPVal[i] = 1;
      arrNPer[i] = fNPerm;
      arrPCha[i] = 0;
      arrPAdj[i] = 0;
   }//for_i

// Create table necessary to permute data
   if (ok) {
      if (!(table = new (nothrow) Double_t*[nrows])) {ok = kFALSE; goto cleantable;} 
      for (Int_t i=0; i<nrows; i++) {
         table[i] = 0;
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         if (!(table[i] = new (nothrow) Double_t[nid])) {ok = kFALSE; goto cleantable;} 
      }//for_i
   cleantable:
      if (!ok) {
         cerr << "Error: Not enough memory to create table to store data:";
         cerr <<        "No permutation can be done to calculate p-chance."
              << endl;
         for (Int_t i=0; i<nrows; i++) {
            if (table[i]) {delete [] table[i]; table[i] = 0;}
         }//for_i
         if (table) delete [] table;
      }//if
   }//if

// Read data
   while (input.good()) {
      if (hasRow2) input.getline(nextline, kLineBuf, delim);
      if (input.fail() || (idx == nrows)) break;
      hasRow2 = kTRUE;  //necessary in case of missing group row

      // read row name
      strcpy(name,strtok(&nextline[0], sepi));
      arrName[idx] = name;

      Double_t value = 0;
      Double_t mn1, mn2, se, df, stat, pval, pcha, padj;
      mn1 = se = df = stat = 0;
      mn2 = NA_REAL;
      pval = pcha = padj = NA_REAL;
      if (fTwoSample) {
      //Two-sample test
         // read row data
         Int_t id1 = 0;
         Int_t id2 = 0;
         for (Int_t i=0; i<numdata; i++) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
            value = atof(strtok(NULL, sepi));
            if (gid[i] == 1) grp1[id1++] = value;
            if (gid[i] == 2) grp2[id2++] = value;
         }//for_i

         // do test
         if (fHasNA) {
            stat = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, mu, fNA);
         } else {
            stat = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, mu);
         }//if

         // fill table
         if (ok) {
            for (Int_t i=0; i<n1; i++) table[idx][i]    = grp1[i];
            for (Int_t i=0; i<n2; i++) table[idx][n1+i] = grp2[i];
         }//if
      } else {
      //One-sample test
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         Int_t id1 = 0;
         for (Int_t i=0; i<numdata; i++) {
            value = atof(strtok(NULL, sepi));
            if (gid[i] == 1) grp1[id1++] = value;
         }//for_i

         // do test
         if (fHasNA) {
            stat = this->Statistic(n1, grp1, mn1, se, df, mu, fNA);
         } else {
            stat = this->Statistic(n1, grp1, mn1, se, df, mu);
         }//if
      }//if
      arrStat[idx] = stat;
      arrMn1[idx]  = mn1;
      arrMn2[idx]  = mn2;
      arrSE[idx]   = se;
      arrDF[idx]   = df;
      arrPVal[idx] = this->PValue(stat, df);
      idx++;
   }//while

   if (idx != nrows) {
      cout << "Warning: Number of lines read <"    << idx  
           << "> is not equal to number of rows <" << nrows << ">"
           << endl;
   }//if
   nrows = idx;

// Permute data for p-chance
   if (ok && !wy) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
      arrPCha = this->PChance(nrows,nid,table,n1,arrPCha,arrNPer,arrStat);
//      arrPCha = this->PChance(nrows,numdata,table,n1,arrPCha,arrNPer,arrStat);
      if (arrPCha == 0) ok = kFALSE;
   }//if

// Adjust p-values to p-adjust
   if (hasPAdj) {
      if (!fAdjPVal && ok && wy) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         arrPAdj = this->PAdjustWY(nrows,nid,table,n1,arrPCha,arrNPer,arrStat,arrPAdj);
         if (arrPCha == 0) {hasPAdj = kFALSE; ok = kFALSE;}
      } else if (!fAdjPVal && ok) {
         arrPAdj = this->PAdjust(nrows, arrPCha, arrPAdj);
      } else {
         arrPAdj = this->PAdjust(nrows, arrPVal, arrPAdj);
      }//if
   }//if

// Write header
   output << "Name";
   if (hasStat)       output << sepo << "Statistic";
   if (hasMn1)        output << sepo << "Mean1";
   if (hasMn2 && ok2) output << sepo << "Mean2";
   if (hasSE)         output << sepo << "SE";
   if (hasDF)         output << sepo << "DF";
   if (hasPVal)       output << sepo << "P-Value";
   if (hasNPer && ok) output << sepo << "NumPerm";
   if (hasPCha && ok) output << sepo << "P-Chance";
   if (hasPAdj)       output << sepo << "P-Adjusted";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nrows; i++) {
      output << arrName[i];
      if (hasStat)       output << sepo << arrStat[i];
      if (hasMn1)        output << sepo << arrMn1[i];
      if (hasMn2 && ok2) output << sepo << arrMn2[i];
      if (hasSE)         output << sepo << arrSE[i];
      if (hasDF)         output << sepo << arrDF[i];
      if (hasPVal)       output << sepo << arrPVal[i];
      if (hasNPer && ok) output << sepo << arrNPer[i];
      if (hasPCha && ok) output << sepo << arrPCha[i];
      if (hasPAdj)       output << sepo << arrPAdj[i];
      output << endl;
   }//for_i

// Cleanup
cleanup:
   if (arrName) {delete [] arrName; arrName = 0;}
   if (arrStat) {delete [] arrStat; arrStat = 0;}
   if (arrMn1)  {delete [] arrMn1;  arrMn1  = 0;}
   if (arrMn2)  {delete [] arrMn2;  arrMn2  = 0;}
   if (arrSE)   {delete [] arrSE;   arrSE   = 0;}
   if (arrDF)   {delete [] arrDF;   arrDF   = 0;}
   if (arrPVal) {delete [] arrPVal; arrPVal = 0;}
   if (arrNPer) {delete [] arrNPer; arrNPer = 0;}
   if (arrPCha) {delete [] arrPCha; arrPCha = 0;}
   if (arrPAdj) {delete [] arrPAdj; arrPAdj = 0;}
   if (gid)  delete [] gid;
   if (grp1) delete [] grp1;
   if (grp2) delete [] grp2;

// Delete table
   if (ok) {
      for (Int_t i=0; i<nrows; i++) {
         if (table[i]) {delete [] table[i]; table[i] = 0;}
      }//for_i
      if (table) delete [] table;
   }//if

// Close files
   input.close();
   output.close();

   return err;
}//Test

//______________________________________________________________________________
Int_t TUnivariateTest::Test(Int_t n, Int_t *gid, TTree *intree, const char *leafname,
                       TTree *outtree, const char *varlist, Int_t nperm, Double_t mu)
{
   // Do multiple univariate tests for intree and store result in outtree.
   // intree must be a collection of n tree friends containing leaf leafname
   // gid is an array containing the group IDs for the n treefriends:
   //     gid=1:  for group1 
   //     gid=2:  for group2
   //     gid<=0: trees at corresponding gid position(s) will be ignored
   // varlist is the list of variables to calculate and to store in outtree:
   //    "stat" - univariate statistic
   //    "mn1"  - mean value of group1
   //    "mn2"  - mean value of group2  (not used for one-sample test)
   //    "se"   - standard error
   //    "df"   - degrees of freedom
   //    "pval" - p-value
   //    "nper" - number of permutations used to determine p-chance
   //    "pcha" - p-chance
   //    "padj" - p-value or p-chance adjusted for multiple comparisons
   //    e.g. for usual test, varlist = "stat:mn1:mn2:se:df:pval"
   //    "*" is equal to "stat:mn1:mn2:se:df:pval:nper:pcha:padj"
   // nperm is number of permutations to be done
   if(kCS) cout << "------TUnivariateTest::Test------" << endl;

   if (nperm >= 0) fNPerm = nperm;
   else            nperm  = fNPerm;
   if (mu != 0)    fMu    = mu;

   Int_t err = 0;
   Int_t n1  = 0;
   Int_t n2  = 0;
   Bool_t ok=kFALSE, ok2, wy;

   if (intree == 0 || outtree == 0) {
      cerr << "Error: Intree and/or outtree is missing." << endl;
      return 1;
   }//if

// Get list of tree friends
   TFriendElement *fe = 0;
   TList *friends = intree->GetListOfFriends();
//BETTER: if (!friends) {cerr << "..." << endl;}
   Int_t nfriends = friends->GetSize();
   Int_t nrows    = (Int_t)(intree->GetEntries());
   if (nfriends+1 != n) {
      cerr << "Error: Number of trees <" << (nfriends+1) 
           << "> is not equal to size of group array <" << n << ">."
           << endl;
      return 1;
   }//if

// Get parent tree
   fe = (TFriendElement*)friends->At(0);
   TTree   *tree0 = fe->GetParentTree();
   TLeaf   *leaf0 = tree0->FindLeaf(leafname);
   TBranch *brch0 = leaf0->GetBranch();
   TTree   *treej = 0;
   TBranch *brchj = 0;
   TLeaf   *leafj = 0;
//PROBLEM: need to check if leaf0 has correct leafname
//PROBLEM: need to check if leafj (for correct gid!) have correct leafname!!!!

// Init local arrays
   Double_t *grp1 = 0;
   Double_t *grp2 = 0;

   Double_t *arrStat = 0;  //statistic
   Double_t *arrMn1  = 0;  //mean of grp1
   Double_t *arrMn2  = 0;  //mean of grp2
   Double_t *arrSE   = 0;  //standard error
   Double_t *arrDF   = 0;  //degree of freedom
   Double_t *arrPVal = 0;  //p-value
   Int_t    *arrNPer = 0;  //number of permutations
   Double_t *arrPCha = 0;  //p-chance
   Double_t *arrPAdj = 0;  //p-adjusted
   Double_t **table  = 0;  //optional table to store infile data

// Decompose varlist
   Bool_t hasStat = kFALSE;
   Bool_t hasMn1  = kFALSE;
   Bool_t hasMn2  = kFALSE;
   Bool_t hasSE   = kFALSE;
   Bool_t hasDF   = kFALSE;
   Bool_t hasPVal = kFALSE;
   Bool_t hasNPer = kFALSE;
   Bool_t hasPCha = kFALSE;
   Bool_t hasPAdj = kFALSE;

   if (strcmp(varlist,"*") == 0) {
      hasStat = kTRUE;
      hasMn1  = kTRUE;
      hasMn2  = kTRUE;
      hasSE   = kTRUE;
      hasDF   = kTRUE;
      hasPVal = kTRUE;
      hasNPer = kTRUE;
      hasPCha = kTRUE;
      hasPAdj = kTRUE;
   } else {
      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;
      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname,"stat") == 0) {hasStat = kTRUE;}
         if (strcmp(vname,"mn1")  == 0) {hasMn1  = kTRUE;}
         if (strcmp(vname,"mn2")  == 0) {hasMn2  = kTRUE;}
         if (strcmp(vname,"se")   == 0) {hasSE   = kTRUE;}
         if (strcmp(vname,"df")   == 0) {hasDF   = kTRUE;}
         if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE;}
         if (strcmp(vname,"nper") == 0) {hasNPer = kTRUE;}
         if (strcmp(vname,"pcha") == 0) {hasPCha = kTRUE;}
         if (strcmp(vname,"padj") == 0) {hasPAdj = kTRUE;}
         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while
      delete [] dname;
   }//if

//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
// Create arrays grp1 and grp2
   Int_t nid = 0;
   for (Int_t i=0; i<n; i++) {
      if (gid[i] == 1) {n1++; nid++;}
      if (gid[i] == 2) {n2++; nid++;}
   }//if

// Check size of groups
   if (nid < 2) {
      cerr << "Error: Less than two trees selected" << endl;
      return 1;
   }//if

//   for (Int_t i=0; i<n; i++) (gid[i] == 1) ? n1++ : n2++;
   if (!(grp1 = new (nothrow) Double_t[n1])) {err = 1; goto cleanup;} 
   if (!(grp2 = new (nothrow) Double_t[n2])) {err = 1; goto cleanup;} 

// Check if one/two-sample test
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
   fTwoSample = (n1 < nid) ? kTRUE : kFALSE;
   ok  = (fTwoSample && nperm > 0);  //permute two-sample/paired sample
   ok2 = (fTwoSample && !fPaired);   //two-sample only

// Check for Westfall-Young adjustment
   wy = (strcmp(fAdjustment, "wy") == 0);

// Create local arrays
   if (!(arrStat = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrMn1  = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrMn2  = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrSE   = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrDF   = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrPVal = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrNPer = new (nothrow) Int_t[nrows]))    {err = 1; goto cleanup;}
   if (!(arrPCha = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrPAdj = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}

// Init local arrays
   for (Int_t i=0; i<nrows; i++) {
      arrStat[i] = 1;
      arrMn1[i]  = fMu;
      arrMn2[i]  = fMu;
      arrSE[i]   = 0;
      arrDF[i]   = 0;
      arrPVal[i] = 1;
      arrNPer[i] = fNPerm;
      arrPCha[i] = 0;
      arrPAdj[i] = 0;
   }//for_i

// Create table necessary to permute data
   if (ok) {
      if (!(table = new (nothrow) Double_t*[nrows])) {ok = kFALSE; goto cleantable;} 
      for (Int_t i=0; i<nrows; i++) {
         table[i] = 0;
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         if (!(table[i] = new (nothrow) Double_t[nid])) {ok = kFALSE; goto cleantable;} 
      }//for_i
   cleantable:
      if (!ok) {
         cerr << "Error: Not enough memory to create table to store data:";
         cerr <<        "No permutation can be done to calculate p-chance."
              << endl;
         for (Int_t i=0; i<nrows; i++) {
            if (table[i]) {delete [] table[i]; table[i] = 0;}
         }//for_i
         if (table) delete [] table;
      }//if
   }//if

// Read intree entries
   Double_t mn1, mn2, se, df, stat, pval, pcha, padj;
   Int_t perm;
   for (Int_t i=0; i<nrows; i++) {
      mn1 = se = df = stat = 0.0;
      mn2 = NA_REAL;
      pval = pcha = padj = NA_REAL;
      if (fTwoSample) {
      //Two-sample test
         // read entry i from tree plus friends
         Int_t id1 = 0;
         Int_t id2 = 0;

         brch0->GetEntry(i);
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         if (gid[0] == 1) grp1[id1++] = leaf0->GetValue();
         if (gid[0] == 2) grp2[id2++] = leaf0->GetValue();

         for (Int_t j=0; j<nfriends; j++) {
            if (gid[j+1] <= 0) continue;

            fe = (TFriendElement*)friends->At(j);
            treej = fe->GetTree();
            leafj = treej->FindLeaf(leafname);
            brchj = leafj->GetBranch();
            brchj->GetEntry(i);
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
            if (gid[j+1] == 1) grp1[id1++] = leafj->GetValue();
            if (gid[j+1] == 2) grp2[id2++] = leafj->GetValue();
         }//for_j

         // do test
         if (fHasNA) {
            stat = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, mu, fNA);
         } else {
            stat = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, mu);
         }//if

         // fill table
         if (ok) {
            for (Int_t k=0; k<n1; k++) table[i][k]    = grp1[k];
            for (Int_t k=0; k<n2; k++) table[i][n1+k] = grp2[k];
         }//if
      } else {
      //One-sample test
         // read entry i from tree plus friends
         Int_t id1 = 0;

         brch0->GetEntry(i);
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         if (gid[0] == 1) grp1[id1++] = leaf0->GetValue();

         for (Int_t j=0; j<nfriends; j++) {
            if (gid[j+1] <= 0) continue;

            fe = (TFriendElement*)friends->At(j);
            treej = fe->GetTree();
            leafj = treej->FindLeaf(leafname);
            brchj = leafj->GetBranch();
            brchj->GetEntry(i);
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
            if (gid[j+1] == 1) grp1[id1++] = leafj->GetValue();
         }//for_j

         // do test
         if (fHasNA) {
            stat = this->Statistic(n1, grp1, mn1, se, df, mu, fNA);
         } else {
            stat = this->Statistic(n1, grp1, mn1, se, df, mu);
         }//if
      }//if
      arrStat[i] = stat;
      arrMn1[i]  = mn1;
      arrMn2[i]  = mn2;
      arrSE[i]   = se;
      arrDF[i]   = df;
      arrPVal[i] = this->PValue(stat, df);
   }//for_i

// Permute data for p-chance
   if (ok && !wy) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
      arrPCha = this->PChance(nrows,nid,table,n1,arrPCha,arrNPer,arrStat);
      if (arrPCha == 0) ok = kFALSE;
   }//if

// Adjust p-values to p-adjust
   if (hasPAdj) {
      if (!fAdjPVal && ok && wy) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         arrPAdj = this->PAdjustWY(nrows,nid,table,n1,arrPCha,arrNPer,arrStat,arrPAdj);
         if (arrPCha == 0) {hasPAdj = kFALSE; ok = kFALSE;}
      } else if (!fAdjPVal && ok) {
         arrPAdj = this->PAdjust(nrows, arrPCha, arrPAdj);
      } else {
         arrPAdj = this->PAdjust(nrows, arrPVal, arrPAdj);
      }//if
   }//if

// Create branches for outtree
   mn1 = mn2 = se = df = stat = 0.0;
   pval = pcha = padj = 0.0;
   perm = 0;
   if (hasStat)       outtree->Branch("statBr", &stat, "stat/D");
   if (hasMn1)        outtree->Branch("mn1Br",  &mn1,  "mn1/D");
   if (hasMn2 && ok2) outtree->Branch("mn2Br",  &mn2,  "mn2/D");
   if (hasSE)         outtree->Branch("seBr",   &se,   "se/D");
   if (hasDF)         outtree->Branch("dfBr",   &df,   "df/D");
   if (hasPVal)       outtree->Branch("pvalBr", &pval, "pval/D");
   if (hasNPer && ok) outtree->Branch("permBr", &perm, "perm/D");
   if (hasPCha && ok) outtree->Branch("pchaBr", &pcha, "pcha/D");
   if (hasPAdj)       outtree->Branch("padjBr", &padj, "padj/D");

// Fill outtree
   for (Int_t i=0; i<nrows; i++) {
      if (hasStat)       stat = arrStat[i];
      if (hasMn1)        mn1  = arrMn1[i];
      if (hasMn2 && ok2) mn2  = arrMn2[i];
      if (hasSE)         se   = arrSE[i];
      if (hasDF)         df   = arrDF[i];
      if (hasPVal)       pval = arrPVal[i];
      if (hasNPer && ok) perm = arrNPer[i];
      if (hasPCha && ok) pcha = arrPCha[i];
      if (hasPAdj)       padj = arrPAdj[i];
      outtree->Fill();
   }//for_i

// Cleanup
cleanup:
   if (arrStat) {delete [] arrStat; arrStat = 0;}
   if (arrMn1)  {delete [] arrMn1;  arrMn1  = 0;}
   if (arrMn2)  {delete [] arrMn2;  arrMn2  = 0;}
   if (arrSE)   {delete [] arrSE;   arrSE   = 0;}
   if (arrDF)   {delete [] arrDF;   arrDF   = 0;}
   if (arrPVal) {delete [] arrPVal; arrPVal = 0;}
   if (arrNPer) {delete [] arrNPer; arrNPer = 0;}
   if (arrPCha) {delete [] arrPCha; arrPCha = 0;}
   if (arrPAdj) {delete [] arrPAdj; arrPAdj = 0;}
   if (grp1) delete [] grp1;
   if (grp2) delete [] grp2;

// Delete table
   if (ok) {
      for (Int_t i=0; i<nrows; i++) {
         if (table[i]) {delete [] table[i]; table[i] = 0;}
      }//for_i
      if (table) delete [] table;
   }//if

   return err;
}//Test

//______________________________________________________________________________
Int_t TUnivariateTest::Test(Int_t n, Int_t *gid, TTree **intree, const char *leafname,
                       TTree *outtree, const char *varlist, Int_t nperm, Double_t mu)
{
   // Do multiple univariate tests for intree and store result in outtree.
   // intree is an array of n trees containing leaf leafname
   // gid is an array containing the group IDs for the n trees:
   //     gid=1:  for group1 
   //     gid=2:  for group2
   //     gid<=0: trees at corresponding gid position(s) will be ignored
   // varlist is the list of variables to calculate and to store in outtree:
   //    "stat" - univariate statistic
   //    "mn1"  - mean value of group1
   //    "mn2"  - mean value of group2  (not used for one-sample test)
   //    "se"   - standard error
   //    "df"   - degrees of freedom
   //    "pval" - p-value
   //    "nper" - number of permutations used to determine p-chance
   //    "pcha" - p-chance
   //    "padj" - p-value or p-chance adjusted for multiple comparisons
   //    e.g. for usual test, varlist = "stat:mn1:mn2:se:df:pval"
   //    "*" is equal to "stat:mn1:mn2:se:df:pval:nper:pcha:padj"
   // nperm is number of permutations to be done
   if(kCS) cout << "------TUnivariateTest::Test------" << endl;

   if (nperm >= 0) fNPerm = nperm;
   else            nperm  = fNPerm;
   if (mu != 0)    fMu    = mu;

   Int_t err = 0;
   Int_t n1  = 0;
   Int_t n2  = 0;
   Bool_t ok=kFALSE, ok2, wy;

   if (intree == 0 || outtree == 0) {
      cerr << "Error: Intree and/or outtree is missing." << endl;
      return 1;
   }//if

   Int_t nrows = (Int_t)(intree[0]->GetEntries());

   TBranch **brch = new TBranch*[n];
   TLeaf   **leaf = new TLeaf*[n];
//PROBLEM: need to check if leaf has correct leafname
//PROBLEM: need to check if leaf (for correct gid!) have correct leafname!!!!

   for (Int_t j=0; j<n; j++) {
      leaf[j] = intree[j]->FindLeaf(leafname);      
      brch[j] = leaf[j]->GetBranch();
   }//for_j

// Init local arrays
   Double_t *grp1 = 0;
   Double_t *grp2 = 0;

   Double_t *arrStat = 0;  //statistic
   Double_t *arrMn1  = 0;  //mean of grp1
   Double_t *arrMn2  = 0;  //mean of grp2
   Double_t *arrSE   = 0;  //standard error
   Double_t *arrDF   = 0;  //degree of freedom
   Double_t *arrPVal = 0;  //p-value
   Int_t    *arrNPer = 0;  //number of permutations
   Double_t *arrPCha = 0;  //p-chance
   Double_t *arrPAdj = 0;  //p-adjusted
   Double_t **table  = 0;  //optional table to store infile data

// Decompose varlist
   Bool_t hasStat = kFALSE;
   Bool_t hasMn1  = kFALSE;
   Bool_t hasMn2  = kFALSE;
   Bool_t hasSE   = kFALSE;
   Bool_t hasDF   = kFALSE;
   Bool_t hasPVal = kFALSE;
   Bool_t hasNPer = kFALSE;
   Bool_t hasPCha = kFALSE;
   Bool_t hasPAdj = kFALSE;

   if (strcmp(varlist,"*") == 0) {
      hasStat = kTRUE;
      hasMn1  = kTRUE;
      hasMn2  = kTRUE;
      hasSE   = kTRUE;
      hasDF   = kTRUE;
      hasPVal = kTRUE;
      hasNPer = kTRUE;
      hasPCha = kTRUE;
      hasPAdj = kTRUE;
   } else {
      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;
      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname,"stat") == 0) {hasStat = kTRUE;}
         if (strcmp(vname,"mn1")  == 0) {hasMn1  = kTRUE;}
         if (strcmp(vname,"mn2")  == 0) {hasMn2  = kTRUE;}
         if (strcmp(vname,"se")   == 0) {hasSE   = kTRUE;}
         if (strcmp(vname,"df")   == 0) {hasDF   = kTRUE;}
         if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE;}
         if (strcmp(vname,"nper") == 0) {hasNPer = kTRUE;}
         if (strcmp(vname,"pcha") == 0) {hasPCha = kTRUE;}
         if (strcmp(vname,"padj") == 0) {hasPAdj = kTRUE;}
         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while
      delete [] dname;
   }//if

//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
// Create arrays grp1 and grp2
   Int_t nid = 0;
   for (Int_t i=0; i<n; i++) {
      if (gid[i] == 1) {n1++; nid++;}
      if (gid[i] == 2) {n2++; nid++;}
   }//if

// Check size of groups
   if (nid < 2) {
      cerr << "Error: Less than two trees selected" << endl;
      return 1;
   }//if

//   for (Int_t i=0; i<n; i++) (gid[i] == 1) ? n1++ : n2++;
   if (!(grp1 = new (nothrow) Double_t[n1])) {err = 1; goto cleanup;} 
   if (!(grp2 = new (nothrow) Double_t[n2])) {err = 1; goto cleanup;} 

// Check if one/two-sample test
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
   fTwoSample = (n1 < nid) ? kTRUE : kFALSE;
   ok  = (fTwoSample && nperm > 0);  //permute two-sample/paired sample
   ok2 = (fTwoSample && !fPaired);   //two-sample only

// Check for Westfall-Young adjustment
   wy = (strcmp(fAdjustment, "wy") == 0);

// Create local arrays
   if (!(arrStat = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrMn1  = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrMn2  = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrSE   = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrDF   = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrPVal = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrNPer = new (nothrow) Int_t[nrows]))    {err = 1; goto cleanup;}
   if (!(arrPCha = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}
   if (!(arrPAdj = new (nothrow) Double_t[nrows])) {err = 1; goto cleanup;}

// Init local arrays
   for (Int_t i=0; i<nrows; i++) {
      arrStat[i] = 1;
      arrMn1[i]  = fMu;
      arrMn2[i]  = fMu;
      arrSE[i]   = 0;
      arrDF[i]   = 0;
      arrPVal[i] = 1;
      arrNPer[i] = fNPerm;
      arrPCha[i] = 0;
      arrPAdj[i] = 0;
   }//for_i

// Create table necessary to permute data
   if (ok) {
      if (!(table = new (nothrow) Double_t*[nrows])) {ok = kFALSE; goto cleantable;} 
      for (Int_t i=0; i<nrows; i++) {
         table[i] = 0;
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         if (!(table[i] = new (nothrow) Double_t[nid])) {ok = kFALSE; goto cleantable;} 
      }//for_i
   cleantable:
      if (!ok) {
         cerr << "Error: Not enough memory to create table to store data:";
         cerr <<        "No permutation can be done to calculate p-chance."
              << endl;
         for (Int_t i=0; i<nrows; i++) {
            if (table[i]) {delete [] table[i]; table[i] = 0;}
         }//for_i
         if (table) delete [] table;
      }//if
   }//if

// Read intree entries
   Double_t mn1, mn2, se, df, stat, pval, pcha, padj;
   Int_t perm;
   for (Int_t i=0; i<nrows; i++) {
      mn1 = se = df = stat = 0.0;
      mn2 = NA_REAL;
      pval = pcha = padj = NA_REAL;
      if (fTwoSample) {
      //Two-sample test
         // read entry i from trees
         Int_t id1 = 0;
         Int_t id2 = 0;

         for (Int_t j=0; j<n; j++) {
            if (gid[j] <= 0) continue;

            brch[j]->GetEntry(i);
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
            if (gid[j] == 1) grp1[id1++] = leaf[j]->GetValue();
            if (gid[j] == 2) grp2[id2++] = leaf[j]->GetValue();
         }//for_j

         // do test
         if (fHasNA) {
            stat = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, mu, fNA);
         } else {
            stat = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, mu);
         }//if

         // fill table
         if (ok) {
            for (Int_t k=0; k<n1; k++) table[i][k]    = grp1[k];
            for (Int_t k=0; k<n2; k++) table[i][n1+k] = grp2[k];
         }//if
      } else {
      //One-sample test
         // read entry i from tree plus friends
         Int_t id1 = 0;

         for (Int_t j=0; j<n; j++) {
            if (gid[j] <= 0) continue;

            brch[j]->GetEntry(i);
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
            if (gid[j] == 1) grp1[id1++] = leaf[j]->GetValue();
         }//for_j

         // do test
         if (fHasNA) {
            stat = this->Statistic(n1, grp1, mn1, se, df, mu, fNA);
         } else {
            stat = this->Statistic(n1, grp1, mn1, se, df, mu);
         }//if
      }//if
      arrStat[i] = stat;
      arrMn1[i]  = mn1;
      arrMn2[i]  = mn2;
      arrSE[i]   = se;
      arrDF[i]   = df;
      arrPVal[i] = this->PValue(stat, df);
   }//for_i

// Permute data for p-chance
   if (ok && !wy) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
      arrPCha = this->PChance(nrows,nid,table,n1,arrPCha,arrNPer,arrStat);
      if (arrPCha == 0) ok = kFALSE;
   }//if

// Adjust p-values to p-adjust
   if (hasPAdj) {
      if (!fAdjPVal && ok && wy) {
//TO DO: check if gid=1, gid=2 or gid=-1 (ignore)
         arrPAdj = this->PAdjustWY(nrows,nid,table,n1,arrPCha,arrNPer,arrStat,arrPAdj);
         if (arrPCha == 0) {hasPAdj = kFALSE; ok = kFALSE;}
      } else if (!fAdjPVal && ok) {
         arrPAdj = this->PAdjust(nrows, arrPCha, arrPAdj);
      } else {
         arrPAdj = this->PAdjust(nrows, arrPVal, arrPAdj);
      }//if
   }//if

// Create branches for outtree
   mn1 = mn2 = se = df = stat = 0.0;
   pval = pcha = padj = 0.0;
   perm = 0;
   if (hasStat)       outtree->Branch("statBr", &stat, "stat/D");
   if (hasMn1)        outtree->Branch("mn1Br",  &mn1,  "mn1/D");
   if (hasMn2 && ok2) outtree->Branch("mn2Br",  &mn2,  "mn2/D");
   if (hasSE)         outtree->Branch("seBr",   &se,   "se/D");
   if (hasDF)         outtree->Branch("dfBr",   &df,   "df/D");
   if (hasPVal)       outtree->Branch("pvalBr", &pval, "pval/D");
   if (hasNPer && ok) outtree->Branch("permBr", &perm, "perm/D");
   if (hasPCha && ok) outtree->Branch("pchaBr", &pcha, "pcha/D");
   if (hasPAdj)       outtree->Branch("padjBr", &padj, "padj/D");

// Fill outtree
   for (Int_t i=0; i<nrows; i++) {
      if (hasStat)       stat = arrStat[i];
      if (hasMn1)        mn1  = arrMn1[i];
      if (hasMn2 && ok2) mn2  = arrMn2[i];
      if (hasSE)         se   = arrSE[i];
      if (hasDF)         df   = arrDF[i];
      if (hasPVal)       pval = arrPVal[i];
      if (hasNPer && ok) perm = arrNPer[i];
      if (hasPCha && ok) pcha = arrPCha[i];
      if (hasPAdj)       padj = arrPAdj[i];
      outtree->Fill();
   }//for_i

// Cleanup
cleanup:
   if (arrStat) {delete [] arrStat; arrStat = 0;}
   if (arrMn1)  {delete [] arrMn1;  arrMn1  = 0;}
   if (arrMn2)  {delete [] arrMn2;  arrMn2  = 0;}
   if (arrSE)   {delete [] arrSE;   arrSE   = 0;}
   if (arrDF)   {delete [] arrDF;   arrDF   = 0;}
   if (arrPVal) {delete [] arrPVal; arrPVal = 0;}
   if (arrNPer) {delete [] arrNPer; arrNPer = 0;}
   if (arrPCha) {delete [] arrPCha; arrPCha = 0;}
   if (arrPAdj) {delete [] arrPAdj; arrPAdj = 0;}
   if (grp1) delete [] grp1;
   if (grp2) delete [] grp2;

// Delete table
   if (ok) {
      for (Int_t i=0; i<nrows; i++) {
         if (table[i]) {delete [] table[i]; table[i] = 0;}
      }//for_i
      if (table) delete [] table;
   }//if

   delete [] leaf;
   delete [] brch;

   return err;
}//Test

//______________________________________________________________________________
Double_t TUnivariateTest::PChance(Int_t n1, Double_t *grp1, Int_t n2, Double_t *grp2,
                 Int_t nperm, Double_t stat)
{
   // Return p-chance as a measure how often permutation of the (n1+n2) data
   // into two random groups results in a statistic that is smaller than
   // the calculated stat value (depending on the given alternative).
   // If the number of permutation "nperm" is larger than the number of 
   //  all possible permutations, then all possible permutations are used.
   if(kCSa) cout << "------TUnivariateTest::PChance------" << endl;

// Check size of groups
   if ((n1 < 2) || (n2 < 2)) {
      cerr << "Error: Less than two values in one of the groups" << endl;
      return NA_REAL;
   }//if

// Check if stat has been calculated
   if (TMLMath::IsNaN(stat)) {
      cerr << "Error: Need to calculate statistic first!" << endl;
      return NA_REAL;
   }//if

// Determine number of permutations
   Int_t n  = n1 + n2;
   Int_t bc = (Int_t)TStat::BinomCoeff(n, n1);

// Init array for sampling
   Double_t *arr = 0;
   if (!(arr = new (nothrow) Double_t[n])) {
      cerr << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   for (Int_t i=0; i<n1; i++) arr[i]    = grp1[i];
   for (Int_t i=0; i<n2; i++) arr[n1+i] = grp2[i];

// Calculate p-chance
   Double_t pchance = 0;
   if (nperm < bc) {
      pchance = this->Sample(n, arr, n1, nperm, stat);
   } else {
      pchance = this->Permute(n, arr, n1, bc, stat);
   }//if

   delete [] arr;

   return pchance;
}//PChance

//______________________________________________________________________________
Double_t *TUnivariateTest::PChance(Int_t nrows, Int_t n, Double_t **table,
                           Int_t n1, Double_t *pcha, Int_t *nperm, Double_t *stat)
{
   // Return  array p-chance as a measure how often permutation of each row of
   // table[nrows,n] into two random groups of size n1 and n-n1 results in a
   // statistic that is smaller than the calculated stat value.
   // If the number of permutations fNPerm is larger than the number of 
   //  all possible permutations, then all possible permutations are used.
   // Array nperm returns the actual number of permutations for each row.
   if(kCS) cout << "------TUnivariateTest::PChance(t)------" << endl;

// Check size of groups
   if ((n < 2) || ((n - n1) < 2)) {
      cerr << "Error: Less than two values in one of the groups" << endl;
      return 0;
   }//if

// Determine number of permutations
   Int_t bc = (Int_t)TStat::BinomCoeff(n, n1);

// Calculate p-chance
   if (fNPerm < bc) {
      pcha = this->Sample(nrows, n, table, n1, pcha, nperm, stat);
   } else {
      fNPerm = bc;
      for (Int_t i=0; i<nrows; i++) nperm[i] = bc;
      pcha = this->Permute(nrows, n, table, n1, pcha, nperm, stat);
   }//if

   return pcha;
}//PChance

//______________________________________________________________________________
Double_t TUnivariateTest::Permute(Int_t n, Double_t *arr, Int_t n1, Int_t nperm,
                          Double_t stat)
{
   // Return p-chance as a measure how often permutation of the n data
   // into two random groups results in a statistic that is smaller than
   // the calculated stat value (depending on the given alternative).
   if(kCSa) cout << "------TUnivariateTest::Permute------" << endl;

// Create indices and initialize
   Int_t n2 = n - n1;
   Int_t *idx1 = new Int_t[n1];
   Int_t *idx2 = new Int_t[n2];
   for (Int_t i=0; i<n1; i++) idx1[i] = i;
   for (Int_t i=0; i<n2; i++) idx2[i] = n1 + i;

// Create arrays for permuted data
   Double_t *grp1 = new Double_t[n1];
   Double_t *grp2 = new Double_t[n2];

// Calculate stat(perm) and compare with stat
   Int_t    numperm = nperm;
   Int_t    nchance = 1;
   Double_t s       = 0;
   Double_t pchance = 0;
   if (strcmp(fAlternative, "twosided") == 0) {
      for (Int_t i=1; i<nperm; i++) {
         // fill grp1 and grp2 with permuted data
         TStat::NextPerm(n, n1, idx1, n2, idx2);
         for (Int_t i=0; i<n1; i++) grp1[i] = arr[idx1[i]];
         for (Int_t i=0; i<n2; i++) grp2[i] = arr[idx2[i]];

         // skip if data in one of the arrays are equal
         if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
            numperm--;
            continue;
         }//if

         // test if stat(perm) >= stat
         Double_t mn1, mn2, se, df;
         if (fHasNA) {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
         } else {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
         }//if
         if (TMLMath::IsNaN(s))          numperm--;
         else if (fabs(s) >= fabs(stat)) nchance++;
      }//for_i
   } else if (strcmp(fAlternative, "greater") == 0) {
      for (Int_t i=1; i<nperm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);
         for (Int_t k=0; k<n1; k++) grp1[k] = arr[idx1[k]];
         for (Int_t k=0; k<n2; k++) grp2[k] = arr[idx2[k]];

         if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
            numperm--;
            continue;
         }//if

         Double_t mn1, mn2, se, df;
         if (fHasNA) {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
         } else {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
         }//if
         if (TMLMath::IsNaN(s)) numperm--;
         else if (s < stat)     nchance++;
      }//for_i
   } else if (strcmp(fAlternative, "less") == 0) {
      for (Int_t i=1; i<nperm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);
         for (Int_t i=0; i<n1; i++) grp1[i] = arr[idx1[i]];
         for (Int_t i=0; i<n2; i++) grp2[i] = arr[idx2[i]];

         if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
            numperm--;
            continue;
         }//if

         Double_t mn1, mn2, se, df;
         if (fHasNA) {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
         } else {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
         }//if
         if (TMLMath::IsNaN(s)) numperm--;
         else if (s >= stat)    nchance++;
      }//for_i
   } else {
      cerr << "Error: Alternative not known" << endl;
      pchance = NA_REAL;
      goto cleanup;
   }//if

   fNPerm  = numperm;   //set to actual number of permutations
   pchance = (numperm ? ((Double_t)nchance / numperm) : 0);

// Cleanup
cleanup:
   delete [] idx1;
   delete [] idx2;
   delete [] grp1;
   delete [] grp2;

   return pchance;
}//Permute

//______________________________________________________________________________
Double_t *TUnivariateTest::Permute(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                           Double_t *pcha, Int_t *nperm, Double_t *stat)
{
   // Return  array p-chance as a measure how often permutation of each row of
   // table[nrows,n] into two random groups of size n1 and n-n1 results in a
   // statistic that is smaller than the calculated stat value (depending 
   // on the given alternative).
   // Array nperm returns the actual number of permutations for each row.
   if(kCS) cout << "------TUnivariateTest::Permute(t)------" << endl;

// Create indices for sampling and initialize
   Int_t n2 = n - n1;
   Int_t *idx1 = new Int_t[n1];
   Int_t *idx2 = new Int_t[n2];
   for (Int_t i=0; i<n1; i++) idx1[i] = i;
   for (Int_t i=0; i<n2; i++) idx2[i] = n1 + i;

// Create arrays for permuted data
   Double_t *grp1 = new Double_t[n1];
   Double_t *grp2 = new Double_t[n2];

// Calculate stat(perm) and compare with stat
   for (Int_t j=0; j<nrows; j++) pcha[j] = 1.0;
   Double_t s = 0;
   if (strcmp(fAlternative, "twosided") == 0) {
      for (Int_t i=1; i<fNPerm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);
         for (Int_t j=0; j<nrows; j++) {
            // fill grp1 and grp2 with permuted data
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            // skip if data in one of the arrays are equal
            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;  //????????
            }//if

            // test if s(perm) >= stat
            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s))             nperm[j]--;
            else if (fabs(s) >= fabs(stat[j])) pcha[j] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "greater") == 0) {
      for (Int_t i=1; i<fNPerm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);
         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s)) nperm[j]--;
            else if (s < stat[j])  pcha[j] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "less") == 0) {
      for (Int_t i=1; i<fNPerm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);
         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s)) nperm[j]--;
            else if (s >= stat[j]) pcha[j] += 1;
         }//for_j
      }//for_i
   } else {
      cerr << "Error: Alternative not known" << endl;
      pcha = 0;
      goto cleanup;
   }//if

   for (Int_t j=0; j<nrows; j++) {
      pcha[j] = (nperm[j] ? (pcha[j] / nperm[j]) : 0);
   }//for_j

// Cleanup
cleanup:
   delete [] idx1;
   delete [] idx2;
   delete [] grp1;
   delete [] grp2;

   return pcha;
}//Permute

//______________________________________________________________________________
Double_t TUnivariateTest::Sample(Int_t n, Double_t *arr, Int_t n1, Int_t nperm,
                          Double_t stat)
{
   // Return p-chance as a measure how often random sampling of the n data
   // into two random groups results in a statistic that is smaller than
   // the calculated stat value (depending on the given alternative).
   if(kCSa) cout << "------TUnivariateTest::Sample------" << endl;

// Create indices and initialize
   Int_t n2 = n - n1;
   Int_t *idx1 = new Int_t[n1];
   Int_t *idx2 = new Int_t[n2];
   for (Int_t i=0; i<n1; i++) idx1[i] = i;
   for (Int_t i=0; i<n2; i++) idx2[i] = n1 + i;

// Create arrays for permuted data
   Double_t *grp1 = new Double_t[n1];
   Double_t *grp2 = new Double_t[n2];

// Calculate tstat(perm) and compare with tstat
   Int_t    numperm = nperm;
   Int_t    nchance = 0;
   Double_t s       = 0;
   Double_t pchance = 0;
   if (strcmp(fAlternative, "twosided") == 0) {
      for (Int_t i=0; i<nperm; i++) {
         // fill grp1 and grp2 with randomized data
         TStat::Sample(n, n1, idx1, n2, idx2);
         for (Int_t i=0; i<n1; i++) grp1[i] = arr[idx1[i]];
         for (Int_t i=0; i<n2; i++) grp2[i] = arr[idx2[i]];

         // skip if data in one of the arrays are equal
         if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
            numperm--;
            continue;
         }//if

         // test if s(perm) >= tstat
         Double_t mn1, mn2, se, df;
         if (fHasNA) {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
         } else {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
         }//if
         if (TMLMath::IsNaN(s))          numperm--;
         else if (fabs(s) >= fabs(stat)) nchance++;
      }//for_i
   } else if (strcmp(fAlternative, "greater") == 0) {
      for (Int_t i=0; i<nperm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);
         for (Int_t i=0; i<n1; i++) grp1[i] = arr[idx1[i]];
         for (Int_t i=0; i<n2; i++) grp2[i] = arr[idx2[i]];

         if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
            numperm--;
            continue;
         }//if

         Double_t mn1, mn2, se, df;
         if (fHasNA) {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
         } else {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
         }//if
         if (TMLMath::IsNaN(s)) numperm--;
         else if (s < stat)     nchance++;
      }//for_i
   } else if (strcmp(fAlternative, "less") == 0) {
      for (Int_t i=0; i<nperm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);
         for (Int_t i=0; i<n1; i++) grp1[i] = arr[idx1[i]];
         for (Int_t i=0; i<n2; i++) grp2[i] = arr[idx2[i]];

         if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
            numperm--;
            continue;
         }//if

         Double_t mn1, mn2, se, df;
         if (fHasNA) {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
         } else {
            s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
         }//if
         if (TMLMath::IsNaN(s)) numperm--;
         else if (s >= stat)    nchance++;
      }//for_i
   } else {
      cerr << "Error: Alternative not known" << endl;
      pchance = NA_REAL;
      goto cleanup;
   }//if

   fNPerm  = numperm;   //set to actual number of permutations
   pchance = (numperm ? ((Double_t)nchance / numperm) : 0);

// Cleanup
cleanup:
   delete [] idx1;
   delete [] idx2;
   delete [] grp1;
   delete [] grp2;

   return pchance;
}//Sample

//______________________________________________________________________________
Double_t *TUnivariateTest::Sample(Int_t nrows, Int_t n, Double_t **table, Int_t n1,
                           Double_t *pcha, Int_t *nperm, Double_t *stat)
{
   // Return  array p-chance as a measure how often random sampling of each row
   // of table[nrows,n] into two random groups of size n1 and n-n1 results in a
   // statistic that is smaller than the calculated stat value (depending 
   // on the given alternative).
   if(kCS) cout << "------TUnivariateTest::Sample(t)------" << endl;

// Create indices and initialize
   Int_t n2 = n - n1;
   Int_t *idx1 = new Int_t[n1];
   Int_t *idx2 = new Int_t[n2];
   for (Int_t i=0; i<n1; i++) idx1[i] = i;
   for (Int_t i=0; i<n2; i++) idx2[i] = n1 + i;

// Create arrays for permuted data
   Double_t *grp1 = new Double_t[n1];
   Double_t *grp2 = new Double_t[n2];

// Calculate stat(perm) and compare with tstat
   for (Int_t j=0; j<nrows; j++) pcha[j] = 0;
   Double_t s = 0;
   if (strcmp(fAlternative, "twosided") == 0) {
      for (Int_t i=0; i<fNPerm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);
         for (Int_t j=0; j<nrows; j++) {
            // fill grp1 and grp2 with permuted data
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            // skip if data in one of the arrays are equal
            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;  //????????
            }//if

            // test if s(perm) >= tstat
            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s))             nperm[j]--;
            else if (fabs(s) >= fabs(stat[j])) pcha[j] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "greater") == 0) {
      for (Int_t i=0; i<fNPerm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);
         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s)) nperm[j]--;
            else if (s < stat[j])  pcha[j] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "less") == 0) {
      for (Int_t i=0; i<fNPerm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);
         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s)) nperm[j]--;
            else if (s >= stat[j]) pcha[j] += 1;
         }//for_j
      }//for_i
   } else {
      cerr << "Error: Alternative not known" << endl;
      pcha = 0;
      goto cleanup;
   }//if

   for (Int_t j=0; j<nrows; j++) {
      pcha[j] = (nperm[j] ? (pcha[j] / nperm[j]) : 0);
   }//for_j

// Cleanup
cleanup:
   delete [] idx1;
   delete [] idx2;
   delete [] grp1;
   delete [] grp2;

   return pcha;
}//Sample

//______________________________________________________________________________
Double_t *TUnivariateTest::PAdjustWY(Int_t nrows, Int_t n, Double_t **table,
                           Int_t n1, Double_t *pcha, Int_t *nperm,
                           Double_t *stat, Double_t *padj)
{
   // Return the Westfall and Young step-down adjusted p-chance values as
   // array p-adjust. Hereby p-chance is a measure how often permutation of 
   // each row of table[nrows,n] into two random groups of size n1 and n-n1
   // results in a statistic that is smaller than the calculated stat value.
   // If the number of permutations fNPerm is larger than the number of 
   // all possible permutations, then all possible permutations are used.
   // Array pcha returns the p-value obtained by permutation of the data.
   // Array nperm returns the actual number of permutations for each row.
   // References:
   // (1) Westfall P.H. and Young S.S. (1993)
   //     Resampling-based multiple testing:examples and methods for
   //     p-value adjustment.
   //     Wiley series in probability and mathematical statistics; Wiley.
   // (2) Dudoit S., Yang Y.H., Callow M.J., Speed T.P.  (2000) 
   //     Statistical methods for identifying differentially expressed genes 
   //     in replicated cDNA microarray experiments. 
   //     UC Berkeley, Technical report #578
   // (3) Manduchi E. (2000)
   //     tpWY, see: http://www.cbil.upenn.edu/tpWY/

   if(kCS) cout << "------TUnivariateTest::PAdjustWY------" << endl;

// Check size of groups
   if ((n < 2) || ((n - n1) < 2)) {
      cerr << "Error: Less than two values in one of the groups" << endl;
      return 0;
   }//if

// Determine number of permutations
   Int_t bc = (Int_t)TStat::BinomCoeff(n, n1);

// Calculate p-chance
   if (fNPerm < bc) {
      padj = this->SampleWY(nrows, n, table, n1, pcha, nperm, stat, padj);
   } else {
      fNPerm = bc;
      for (Int_t i=0; i<nrows; i++) nperm[i] = bc;
      padj = this->PermuteWY(nrows, n, table, n1, pcha, nperm, stat, padj);
   }//if

   return padj;
}//PAdjustWY

//______________________________________________________________________________
Double_t *TUnivariateTest::PermuteWY(Int_t nrows, Int_t n, Double_t **table,
                           Int_t n1, Double_t *pcha, Int_t *nperm,
                           Double_t *stat, Double_t *padj)
{
   // Return the Westfall and Young step-down adjusted p-chance values as
   // array p-adjust. Hereby p-chance is a measure how often permutation of 
   // each row of table[nrows,n] into two groups of size n1 and n-n1 results
   // in a statistic that is smaller than the calculated stat value (depending 
   // on the given alternative).
   // Array pcha returns the p-value obtained by permutation of the data.
   // Array nperm returns the actual number of permutations for each row.
   if(kCS) cout << "------TUnivariateTest::PermuteWY------" << endl;

   Int_t    *idx1 = 0;
   Int_t    *idx2 = 0;
   Int_t    *indx = 0;
   Double_t *grp1 = 0;
   Double_t *grp2 = 0;
   Double_t *s    = 0;
   Double_t *u    = 0;

   Int_t nr1 = nrows - 1;
   Int_t k   = nr1;

// Create indices for sampling and initialize
   Int_t n2 = n - n1;
   if (!(idx1 = new (nothrow) Int_t[n1])) {padj = 0; goto cleanup;}
   if (!(idx2 = new (nothrow) Int_t[n2])) {padj = 0; goto cleanup;}
   for (Int_t i=0; i<n1; i++) idx1[i] = i;
   for (Int_t i=0; i<n2; i++) idx2[i] = n1 + i;

// Create arrays for permuted data
   if (!(grp1 = new (nothrow) Double_t[n1])) {padj = 0; goto cleanup;}
   if (!(grp2 = new (nothrow) Double_t[n2])) {padj = 0; goto cleanup;}

// Create index to sort stat
   if (!(indx = new (nothrow) Int_t[nrows])) {padj = 0; goto cleanup;}

// Init stat(perm), pcha, padj
   if (!(s = new (nothrow) Double_t[nrows])) {padj = 0; goto cleanup;}
   if (!(u = new (nothrow) Double_t[nrows])) {padj = 0; goto cleanup;}
   for (Int_t j=0; j<nrows; j++) s[j]    = 1;
   for (Int_t j=0; j<nrows; j++) u[j]    = 1;
   for (Int_t j=0; j<nrows; j++) pcha[j] = 1;
   for (Int_t j=0; j<nrows; j++) padj[j] = 1;

// Calculate stat(perm) and compare with stat
   if (strcmp(fAlternative, "twosided") == 0) {
      // sort fabs(stat)
      for (Int_t j=0; j<nrows; j++) u[j] = fabs(stat[j]);
      TMath::Sort(nrows, u, indx, kTRUE);

      for (Int_t i=1; i<fNPerm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);

         for (Int_t j=0; j<nrows; j++) {
            // fill grp1 and grp2 with permuted data
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            // skip if data in one of the arrays are equal
            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               s[j] = -1;
               nperm[j]--;
               continue;  //????????
            }//if

            // test if s(perm) >= stat
            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s[j]))    {s[j] = -1; nperm[j]--;}
            else if (fabs(s[j]) >= fabs(stat[j])) {pcha[j] += 1;}
         }//for_j

         k = indx[nr1];
         u[nr1] = s[k];
         if (fabs(u[nr1]) >= fabs(stat[k])) {
            padj[k] += 1;
         }//if
         for (Int_t j=nr1-1; j>=0; j--) {
            k = indx[j];
            if (s[k] == -1) continue; //???????
            u[j] = (fabs(u[j+1]) >= fabs(s[k])) ? u[j+1] : s[k];
            if (fabs(u[j]) >= fabs(stat[k])) padj[k] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "greater") == 0) {
      TMath::Sort(nrows, stat, indx, kFALSE);

      for (Int_t i=1; i<fNPerm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);

         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s[j])) {s[j] = -1; nperm[j]--;}
            else if (s[j] < stat[j])  {pcha[j] += 1;}
         }//for_j

         k = indx[nr1];
         u[nr1] = s[k];
         if (u[nr1] < stat[k]) {
            padj[k] += 1;
         }//if
         for (Int_t j=nr1-1; j>=0; j--) {
            k = indx[j];
            if (s[k] == -1) continue;
            u[j] = (u[j+1] < s[k]) ? u[j+1] : s[k];
            if (u[j] < stat[k]) padj[k] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "less") == 0) {
      TMath::Sort(nrows, stat, indx, kTRUE);

      for (Int_t i=1; i<fNPerm; i++) {
         TStat::NextPerm(n, n1, idx1, n2, idx2);

         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s[j])) {s[j] = -1; nperm[j]--;}
            else if (s[j] >= stat[j]) {pcha[j] += 1;}
         }//for_j

         k = indx[nr1];
         u[nr1] = s[k];
         if (u[nr1] >= stat[k]) {
            padj[k] += 1;
         }//if
         for (Int_t j=nr1-1; j>=0; j--) {
            k = indx[j];
            if (s[k] == -1) continue;
            u[j] = (u[j+1] >= s[k]) ? u[j+1] : s[k];
            if (u[j] >= stat[k]) padj[k] += 1;
         }//for_j
      }//for_i
   } else {
      cerr << "Error: Alternative not known" << endl;
      pcha = 0;
      goto cleanup;
   }//if

   for (Int_t j=0; j<nrows; j++) {
      k = indx[j];
      pcha[j] = (nperm[j] ? (pcha[j] / nperm[j]) : 0);
      padj[k] = (nperm[k] ? (padj[k] / nperm[k]) : 0);
      if (j>=1) padj[k] = (padj[k]>=padj[indx[j-1]]) ? padj[k] : padj[indx[j-1]];
   }//for_j

// Cleanup
cleanup:
   if (s)    delete [] s;
   if (u)    delete [] u;
   if (indx) delete [] indx;
   if (idx1) delete [] idx1;
   if (idx2) delete [] idx2;
   if (grp1) delete [] grp1;
   if (grp2) delete [] grp2;

   return padj;
}//PermuteWY

//______________________________________________________________________________
Double_t *TUnivariateTest::SampleWY(Int_t nrows, Int_t n, Double_t **table,
                           Int_t n1, Double_t *pcha, Int_t *nperm,
                           Double_t *stat, Double_t *padj)
{
   // Return the Westfall and Young step-down adjusted p-chance values as
   // array p-adjust. Hereby p-chance is a measure how often random sampling of 
   // each row of table[nrows,n] into two groups of size n1 and n-n1 results
   // in a statistic that is smaller than the calculated stat value (depending 
   // on the given alternative).
   // Array pcha returns the p-value obtained by permutation of the data.
   // Array nperm returns the actual number of permutations for each row.
   if(kCS) cout << "------TUnivariateTest::SampleWY------" << endl;

   Int_t    *idx1 = 0;
   Int_t    *idx2 = 0;
   Int_t    *indx = 0;
   Double_t *grp1 = 0;
   Double_t *grp2 = 0;
   Double_t *s    = 0;
   Double_t *u    = 0;

   Int_t nr1 = nrows - 1;
   Int_t k   = nr1;

// Create indices for sampling and initialize
   Int_t n2 = n - n1;
   if (!(idx1 = new (nothrow) Int_t[n1])) {padj = 0; goto cleanup;}
   if (!(idx2 = new (nothrow) Int_t[n2])) {padj = 0; goto cleanup;}
   for (Int_t i=0; i<n1; i++) idx1[i] = i;
   for (Int_t i=0; i<n2; i++) idx2[i] = n1 + i;

// Create arrays for permuted data
   if (!(grp1 = new (nothrow) Double_t[n1])) {padj = 0; goto cleanup;}
   if (!(grp2 = new (nothrow) Double_t[n2])) {padj = 0; goto cleanup;}

// Create index to sort stat
   if (!(indx = new (nothrow) Int_t[nrows])) {padj = 0; goto cleanup;}

// Init stat(perm), pcha, padj
   if (!(s = new (nothrow) Double_t[nrows])) {padj = 0; goto cleanup;}
   if (!(u = new (nothrow) Double_t[nrows])) {padj = 0; goto cleanup;}
   for (Int_t j=0; j<nrows; j++) s[j]    = 1;
   for (Int_t j=0; j<nrows; j++) u[j]    = 1;
   for (Int_t j=0; j<nrows; j++) pcha[j] = 1;
   for (Int_t j=0; j<nrows; j++) padj[j] = 1;

// Calculate stat(perm) and compare with tstat
   if (strcmp(fAlternative, "twosided") == 0) {
      // sort fabs(stat)
      for (Int_t j=0; j<nrows; j++) u[j] = fabs(stat[j]);
      TMath::Sort(nrows, u, indx, kTRUE);

      for (Int_t i=0; i<fNPerm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);

         for (Int_t j=0; j<nrows; j++) {
            // fill grp1 and grp2 with permuted data
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            // skip if data in one of the arrays are equal
            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;  //????????
            }//if

            // test if s(perm) >= stat
            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s[j]))    {s[j] = -1; nperm[j]--;}
            else if (fabs(s[j]) >= fabs(stat[j])) {pcha[j] += 1;}
         }//for_j

         k = indx[nr1];
         u[nr1] = s[k];
         if (fabs(u[nr1]) >= fabs(stat[k])) {
            padj[k] += 1;
         }//if
         for (Int_t j=nr1-1; j>=0; j--) {
            k = indx[j];
            if (s[k] == -1) continue; //???????
            u[j] = (fabs(u[j+1]) >= fabs(s[k])) ? u[j+1] : s[k];
            if (fabs(u[j]) >= fabs(stat[k])) padj[k] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "greater") == 0) {
      TMath::Sort(nrows, stat, indx, kFALSE);

      for (Int_t i=0; i<fNPerm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);

         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s[j])) {s[j] = -1; nperm[j]--;}
            else if (s[j] < stat[j])  {pcha[j] += 1;}
         }//for_j

         k = indx[nr1];
         u[nr1] = s[k];
         if (u[nr1] < stat[k]) {
            padj[k] += 1;
         }//if
         for (Int_t j=nr1-1; j>=0; j--) {
            k = indx[j];
            if (s[k] == -1) continue;
            u[j] = (u[j+1] < s[k]) ? u[j+1] : s[k];
            if (u[j] < stat[k]) padj[k] += 1;
         }//for_j
      }//for_i
   } else if (strcmp(fAlternative, "less") == 0) {
      TMath::Sort(nrows, stat, indx, kTRUE);

      for (Int_t i=0; i<fNPerm; i++) {
         TStat::Sample(n, n1, idx1, n2, idx2);

         for (Int_t j=0; j<nrows; j++) {
            for (Int_t k=0; k<n1; k++) grp1[k] = table[j][idx1[k]];
            for (Int_t k=0; k<n2; k++) grp2[k] = table[j][idx2[k]];

            if (TStat::Ident(n1,grp1) || TStat::Ident(n2,grp2)) {
               nperm[j]--;
               continue;
            }//if

            Double_t mn1, mn2, se, df;
            if (fHasNA) {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu, fNA);
            } else {
               s[j] = this->Statistic(n1, grp1, n2, grp2, mn1, mn2, se, df, fMu);
            }//if
            if (TMLMath::IsNaN(s[j])) {s[j] = -1; nperm[j]--;}
            else if (s[j] >= stat[j]) {pcha[j] += 1;}
         }//for_j

         k = indx[nr1];
         u[nr1] = s[k];
         if (u[nr1] >= stat[k]) {
            padj[k] += 1;
         }//if
         for (Int_t j=nr1-1; j>=0; j--) {
            k = indx[j];
            if (s[k] == -1) continue;
            u[j] = (u[j+1] >= s[k]) ? u[j+1] : s[k];
            if (u[j] >= stat[k]) padj[k] += 1;
         }//for_j
      }//for_i
   } else {
      cerr << "Error: Alternative not known" << endl;
      pcha = 0;
      goto cleanup;
   }//if

   for (Int_t j=0; j<nrows; j++) {
      k = indx[j];
      pcha[j] = (nperm[j] ? (pcha[j] / nperm[j]) : 0);
      padj[k] = (nperm[k] ? (padj[k] / nperm[k]) : 0);
      if (j>=1) padj[k] = (padj[k]>=padj[indx[j-1]]) ? padj[k] : padj[indx[j-1]];
   }//for_j

// Cleanup
cleanup:
   if (s)    delete [] s;
   if (u)    delete [] u;
   if (indx) delete [] indx;
   if (idx1) delete [] idx1;
   if (idx2) delete [] idx2;
   if (grp1) delete [] grp1;
   if (grp2) delete [] grp2;

   return padj;
}//SampleWY

//______________________________________________________________________________
Double_t *TUnivariateTest::PAdjust(Int_t n, Double_t *pval, Double_t *padj)
{
   // Calculate adjusted p-value padj
   // Use method "SetAdjustment("adj.method") to select adjustment method.
   // Supported adjustment methods are:
   //    adj.method = "none"         no adjustment
   //    adj.method = "bonferroni"   Bonferroni adjustment
   //    adj.method = "by"           FDR adjustment (Benjamini & Yekutieli)
   //    adj.method = "bh" or "fdr"  FDR adjustment (Benjamini & Hochberg)
   //    adj.method = "hochberg"     Hochberg adjustment
   //    adj.method = "holm"         Holm adjustment
   // Converted from: R-1.6.1/src/library/base/R/p.adjust.R 
   // Additional adjustment methods are:
   //    adj.method = "wy"           Westfall-Young step-down adjustment
   //                                see method "PAdjustWY()"
   if(kCS) cout << "------TUnivariateTest::PAdjust------" << endl;

   if (strcmp(fAdjustment, "none") == 0) {
      for (Int_t i=0; i<n; i++) padj[i] = pval[i];
   } else if (strcmp(fAdjustment, "bonferroni") == 0) {
      padj = this->Bonferroni(n, pval, padj);
   } else if (strcmp(fAdjustment, "by") == 0) {
      padj = this->BY(n, pval, padj);
   } else if (strcmp(fAdjustment, "fdr") == 0 || strcmp(fAdjustment, "bh") == 0) {
      padj = this->FDR(n, pval, padj);
   } else if (strcmp(fAdjustment, "hochberg") == 0) {
      padj = this->Hochberg(n, pval, padj);
   } else if (strcmp(fAdjustment, "holm") == 0) {
      padj = this->Holm(n, pval, padj);
   } else {
      cerr << "Error: Adjustment method not known, using method <none>" << endl;
      for (Int_t i=0; i<n; i++) padj[i] = pval[i];
   }//if

   return padj;
}//PAdjust

//______________________________________________________________________________
Double_t *TUnivariateTest::Bonferroni(Int_t n, Double_t *pval, Double_t *padj)
{
   // Calculate Bonferroni adjustment
   if(kCS) cout << "------TUnivariateTest::Bonferroni------" << endl;

   for (Int_t i=0; i<n; i++) padj[i] = (n*pval[i] < 1) ? n*pval[i] : 1;

   return padj;
}//Bonferroni

//______________________________________________________________________________
Double_t *TUnivariateTest::BY(Int_t n, Double_t *pval, Double_t *padj)
{
   // Calculate FDR adjustment (Benjamini & Yekutieli)
   if(kCS) cout << "------TUnivariateTest::BY------" << endl;

   Double_t q = 0;

   Int_t    *index = 0;
   Int_t    *rank  = 0;
   Double_t *vec1  = 0;
   Double_t *tmp1  = 0;
   Double_t *tmp2  = 0;

   if (!(index = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(rank  = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(vec1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp2  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}

   TStat::Rank(n, pval, index, rank, kTRUE);
   for (Int_t i=0; i<n; i++) vec1[i] = 1;
   for (Int_t i=0; i<n; i++) q += 1.0/(Double_t)(i+1);
   for (Int_t i=0; i<n; i++) tmp1[i] = ((Double_t)n/(Double_t)(n-i))*q*pval[index[i]];
   tmp2 = TStat::CumMin(n, tmp1, tmp2);
   tmp1 = TStat::PMin(n, vec1, tmp2, tmp1);
   for (Int_t i=0; i<n; i++) padj[i] = tmp1[rank[i]];

cleanup:
   delete [] index;
   delete [] rank;
   delete [] vec1;
   delete [] tmp1;
   delete [] tmp2;

   return padj;
}//BY

//______________________________________________________________________________
Double_t *TUnivariateTest::FDR(Int_t n, Double_t *pval, Double_t *padj)
{
   // Calculate FDR adjustment (Benjamini & Hochberg)
   if(kCS) cout << "------TUnivariateTest::FDR------" << endl;

   Int_t    *index = 0;
   Int_t    *rank  = 0;
   Double_t *vec1  = 0;
   Double_t *tmp1  = 0;
   Double_t *tmp2  = 0;

   if (!(index = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(rank  = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(vec1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp2  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}

   TStat::Rank(n, pval, index, rank, kTRUE);
   for (Int_t i=0; i<n; i++) vec1[i] = 1;
   for (Int_t i=0; i<n; i++) tmp1[i] = ((Double_t)n/(Double_t)(n-i))*pval[index[i]];
   tmp2 = TStat::CumMin(n, tmp1, tmp2);
   tmp1 = TStat::PMin(n, vec1, tmp2, tmp1);
   for (Int_t i=0; i<n; i++) padj[i] = tmp1[rank[i]];

cleanup:
   delete [] index;
   delete [] rank;
   delete [] vec1;
   delete [] tmp1;
   delete [] tmp2;

   return padj;
}//FDR

//______________________________________________________________________________
Double_t *TUnivariateTest::Hochberg(Int_t n, Double_t *pval, Double_t *padj)
{
   // Calculate Hochberg adjustment
   if(kCS) cout << "------TUnivariateTest::Hochberg------" << endl;

   Int_t    *index = 0;
   Int_t    *rank  = 0;
   Double_t *vec1  = 0;
   Double_t *tmp1  = 0;
   Double_t *tmp2  = 0;

   if (!(index = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(rank  = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(vec1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp2  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}

   TStat::Rank(n, pval, index, rank, kTRUE);
   for (Int_t i=0; i<n; i++) vec1[i] = 1;
   for (Int_t i=0; i<n; i++) tmp1[i] = (i+1)*pval[index[i]];
   tmp2 = TStat::CumMin(n, tmp1, tmp2);
   tmp1 = TStat::PMin(n, vec1, tmp2, tmp1);
   for (Int_t i=0; i<n; i++) padj[i] = tmp1[rank[i]];

cleanup:
   delete [] index;
   delete [] rank;
   delete [] vec1;
   delete [] tmp1;
   delete [] tmp2;

   return padj;
}//Hochberg

//______________________________________________________________________________
Double_t *TUnivariateTest::Holm(Int_t n, Double_t *pval, Double_t *padj)
{
   // Calculate Holm adjustment
   if(kCS) cout << "------TTTest::Holm------" << endl;

   Int_t    *index = 0;
   Int_t    *rank  = 0;
   Double_t *vec1  = 0;
   Double_t *tmp1  = 0;
   Double_t *tmp2  = 0;

   if (!(index = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(rank  = new (nothrow) Int_t[n]))    {padj = pval; goto cleanup;}
   if (!(vec1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp1  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}
   if (!(tmp2  = new (nothrow) Double_t[n])) {padj = pval; goto cleanup;}

   TStat::Rank(n, pval, index, rank, kFALSE);
   for (Int_t i=0; i<n; i++) vec1[i] = 1;
   for (Int_t i=0; i<n; i++) tmp1[i] = (n-i)*pval[index[i]];
   tmp2 = TStat::CumMax(n, tmp1, tmp2);
   tmp1 = TStat::PMin(n, vec1, tmp2, tmp1);
   for (Int_t i=0; i<n; i++) padj[i] = tmp1[rank[i]];

cleanup:
   delete [] index;
   delete [] rank;
   delete [] vec1;
   delete [] tmp1;
   delete [] tmp2;

   return padj;
}//Holm


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStudentTest                                                         //
//                                                                      //
// Class for Student's t-test                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TStudentTest::TStudentTest()
             :TUnivariateTest()
{
   // Default StudentTest constructor
   if(kCS) cout << "---TStudentTest::TStudentTest(default)------" << endl;

   this->Init();
}//Constructor

//______________________________________________________________________________
TStudentTest::TStudentTest(const char *name, const char *type)
             :TUnivariateTest(name, type)
{
   // Normal StudentTest constructor
   if(kCS) cout << "---TStudentTest::TStudentTest------" << endl;

   this->Init();
}//Constructor

//______________________________________________________________________________
TStudentTest::TStudentTest(const char *name, const char *type, Double_t na)
             :TUnivariateTest(name, type, na)
{
   // Normal StudentTest constructor
   // na is the value used in the data set to indicate missing values
   if(kCS) cout << "---TStudentTest::TStudentTest(na)------" << endl;

   this->Init();
}//Constructor

//______________________________________________________________________________
TStudentTest::~TStudentTest()
{
   // StudentTest destructor
   if(kCS) cout << "---TStudentTest::~TStudentTest------" << endl;

}//Destructor

//______________________________________________________________________________
void TStudentTest::Init(Int_t nperm, Double_t mu, Bool_t paired, Double_t conflevel,
                   Bool_t varequ, const char *alt)
{
   // Initialize t-test
   // not necessesary if default values are used
   // if conflevel < 0 then do not calculate confidence interval
   if(kCS) cout << "---TStudentTest::Init------" << endl;

   fEqualVar = varequ;
   TUnivariateTest::Init(nperm, mu, paired, conflevel, alt);
}//Init

//______________________________________________________________________________
void TStudentTest::PrintInfo()
{
   // Print results of t-test
   if(kCS) cout << "------TStudentTest::PrintInfo------" << endl;

   cout << "==============================================================================" << endl;
   cout << endl;
   if (!fTwoSample) {
      cout << "         One Sample t-test" << endl;
   } else if (fPaired) {
      cout << "         Paired t-test" << endl;
   } else if (fEqualVar) {
      cout << "         Two Sample t-test" << endl;
   } else {
      cout << "         Welch Two Sample t-test" << endl;
   }//if
   cout << endl;

   cout << "t  = " << fStat << endl;
   cout << "df = " << fDF << endl;
   cout << "p-value  = " << fPValue << endl;

   if (fNPerm > 0) {
      cout << "numperm  = " << fNPerm << endl;
      cout << "p-chance = " << fPChance << endl;
   }//if

   cout << "alternative hypothesis: true ";
   if (fPaired || fTwoSample) cout << "difference in means ";
   else                       cout << "mean ";
   if (strcmp(fAlternative, "greater") == 0)   cout << "is greater than ";
   else if (strcmp(fAlternative, "less") == 0) cout << "is less than ";
   else                                        cout << "is not equal to ";
   cout << fMu << endl;

   if (fConfLevel >= 0) {
      cout << 100*fConfLevel << " percent confidence interval:" << endl;
      cout << " [ " << fConfLo << " ,  " << fConfHi << " ]" << endl;
   }//if

   cout << "sample estimates: " << endl;
   cout << "mean(grp1)";
   if (fPaired || !fTwoSample) cout << endl;
   else                        cout << "      mean(grp2)" << endl;
   cout << "   " << fMean1;
   if (fPaired || !fTwoSample) cout << endl;
   else                        cout << "           " << fMean2 << endl;
   cout << endl;

   cout << "==============================================================================" << endl;
}//PrintInfo

//______________________________________________________________________________
Double_t TStudentTest::PValue(Double_t stat, Double_t df)
{
   // Calculate p-value only
   if(kCSa) cout << "------TStudentTest::PValue------" << endl;

   if (TMLMath::IsNaN(stat)) return NA_REAL;

   Double_t pval = 1;
   if (strcmp(fAlternative, "twosided") == 0) {
      pval = 2*TMLMath::PT(-TMath::Abs(stat), df);
   } else if (strcmp(fAlternative, "greater") == 0) {
      pval = TMLMath::PT(stat, df, kFALSE);
   } else if (strcmp(fAlternative, "less") == 0) {
      pval = TMLMath::PT(stat, df, kTRUE);
   } else {
      cerr << "Error: Alternative not known" << endl;
      return NA_REAL;
   }//if

   return pval;
}//PValue

//______________________________________________________________________________
Double_t TStudentTest::PValue(Double_t stat, Double_t df, Double_t se, Double_t &lo,
                       Double_t &hi)
{
   // Calculate p-value and confidence interval [lo, hi]
   if(kCSa) cout << "------TStudentTest::PValue(cf)------" << endl;

   if (TMLMath::IsNaN(stat)) return NA_REAL;

   Double_t pval = 1;
   Double_t cint = 0;
   if (strcmp(fAlternative, "twosided") == 0) {
      pval = 2*TMLMath::PT(-TMath::Abs(stat), df);
      if (fConfLevel >= 0) {
         cint = TMLMath::QT(0.5 + fConfLevel/2.0, df);
         hi   = fMu + (stat + cint)*se;
         lo   = fMu + (stat - cint)*se;
      }//if
   } else if (strcmp(fAlternative, "greater") == 0) {
      pval = TMLMath::PT(stat, df, kFALSE);
      if (fConfLevel >= 0) {
         cint = TMLMath::QT(fConfLevel, df);
         hi   = R_PosInf;
         lo   = fMu + (stat - cint)*se;
      }//if
   } else if (strcmp(fAlternative, "less") == 0) {
      pval = TMLMath::PT(stat, df, kTRUE);
      if (fConfLevel >= 0) {
         cint = TMLMath::QT(fConfLevel, df);
         hi   = fMu + (stat + cint)*se;
         lo   = R_NegInf;
      }//if
   } else {
      cerr << "Error: Alternative not known" << endl;
      return NA_REAL;
   }//if

   return pval;
}//PValue

//______________________________________________________________________________
Double_t TStudentTest::Statistic(Int_t n, Double_t *grp, Double_t &mean,
                       Double_t &se, Double_t &df, Double_t mu)
{
   // Return t-statistic for one-sample t-test
   if(kCSa) cout << "------TStudentTest::Statistic(1)------" << endl;

   if (n < 2) {
      cerr << "Error: Less than two values in group" << endl;
      return NA_REAL;
   }//if

   Double_t mn, var, dfr, err, ts;
   mn  = TStat::Mean(n, grp);
   var = TStat::Var(n, grp, mn);
   dfr = n - 1;
   err = TMath::Sqrt(var/n);
   ts  = (mn - mu)/err;

   mean = mn;
   se   = err;
   df   = dfr;
   return ts;
}//Statistic

//______________________________________________________________________________
Double_t TStudentTest::Statistic(Int_t n, Double_t *grp, Double_t &mean,
                       Double_t &se, Double_t &df, Double_t mu, Double_t na)
{
   // Return t-statistic for one-sample t-test 
   // missing values are allowed but have to be set to a chosen value na
   if(kCSa) cout << "------TStudentTest::Statistic(1na)------" << endl;

   if (n < 2) {
      cerr << "Error: Less than two values in group" << endl;
      return NA_REAL;
   }//if

   Int_t len = n;  //length of array w/o na
   Double_t mn, var, dfr, err, ts;
   mn  = TStat::Mean(n, grp, len, na);
   var = TStat::Var(n, grp, mn, len, na);
   dfr = len - 1;
   err = TMath::Sqrt(var/len);
   ts  = (mn - mu)/err;

   mean = mn;
   se   = err;
   df   = dfr;
   return ts;
}//Statistic

//______________________________________________________________________________
Double_t TStudentTest::Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                       Double_t *grp2, Double_t &mean1, Double_t &mean2,
                       Double_t &se, Double_t &df, Double_t mu)
{
   // Return t-statistic for two-sample t-test
   if(kCSa) cout << "------TStudentTest::Statistic(2)------" << endl;

   Double_t mn1, mn2, err, dfr, ts;
   mn1 = err = dfr = ts = 0.0;
   mn2 = NA_REAL;
   if (fPaired) {
      if (n1 != n2) {
         cerr << "Error: Group1 and group2 must have paired values" << endl;
         return NA_REAL;
      }//if

      // paired grp = grp1 - grp2
      Double_t *grp = 0;
      if (!(grp = new Double_t[n1])) {
         cerr << "Error: Could not initialize memory!" << endl;
         return NA_REAL;
      }//if
      for (Int_t i=0; i<n1; i++) grp[i] = grp1[i] - grp2[i];

      // one-sample t-test
      ts = this->Statistic(n1, grp, mn1, err, dfr, mu);

      if (grp) delete [] grp;
   } else {
      if ((n1 < 2) || (n2 < 2)) {
         cerr << "Error: Less than two values in one of the groups" << endl;
         return NA_REAL;
      }//if

      Double_t var1, var2;
      mn1  = TStat::Mean(n1, grp1);
      mn2  = TStat::Mean(n2, grp2);
      var1 = TStat::Var(n1, grp1, mn1);
      var2 = TStat::Var(n2, grp2, mn2);

      if (fEqualVar) {
         dfr = n1 + n2 - 2;
         Double_t var = ((n1-1)*var1 + (n2-1)*var2) / dfr;
         err = TMath::Sqrt(var*(1./n1 + 1./n2));
      } else {
         Double_t v1 = var1 / n1;
         Double_t v2 = var2 / n2;
         Double_t v  = v1 + v2;
         err = TMath::Sqrt(v);
         dfr = v*v / (v1*v1/(n1-1) + v2*v2/(n2-1));
      }//if

      ts = (mn1 - mn2 - mu)/err;
   }//if

   mean1 = mn1;
   mean2 = mn2;
   se    = err;
   df    = dfr;
   return ts;
}//Statistic

//______________________________________________________________________________
Double_t TStudentTest::Statistic(Int_t n1, Double_t *grp1, Int_t n2,
                       Double_t *grp2, Double_t &mean1, Double_t &mean2,
                       Double_t &se, Double_t &df, Double_t mu, Double_t na)
{
   // Return t-statistic for two-sample t-test 
   // missing values are allowed but have to be set to a chosen value na
   if(kCSa) cout << "------TStudentTest::Statistic(2na)------" << endl;

   Double_t mn1, mn2, err, dfr, ts;
   mn1 = err = dfr = ts = 0.0;
   mn2 = NA_REAL;
   if (fPaired) {
      if (n1 != n2) {
         cerr << "Error: Group1 and group2 must have paired values" << endl;
         return NA_REAL;
      }//if

      // paired grp = grp1 - grp2
      Double_t *grp = 0;
      if (!(grp = new Double_t[n1])) {
         cerr << "Error: Could not initialize memory!" << endl;
         return NA_REAL;
      }//if
      // use only complete cases
      Int_t m = n1;
      for (Int_t i=0; i<n1; i++) {
         if ((grp1[i] != na) && (grp2[i] != na)) grp[i] = grp1[i] - grp2[i];
         else  m--;
      }//for

      // one-sample t-test
      ts = this->Statistic(m, grp, mn1, err, dfr, mu);

      if (grp) delete [] grp;
   } else {
      if ((n1 < 2) || (n2 < 2)) {
         cerr << "Error: Less than two values in one of the groups" << endl;
         return NA_REAL;
      }//if

      Int_t len1 = n1;  //length of grp1 w/o na
      Int_t len2 = n2;  //length of grp2 w/o na
      Double_t var1, var2;
      mn1  = TStat::Mean(n1, grp1, len1, na);
      mn2  = TStat::Mean(n2, grp2, len2, na);
      var1 = TStat::Var(n1, grp1, mn1, len1, na);
      var2 = TStat::Var(n2, grp2, mn2, len2, na);

      if ((len1 < 2) || (len2 < 2)) {
         // prevent lots of error messages during random sampling
         if (fNPerm <= 0) {
            cerr << "Error: Less than 2 non-missing values in one of the groups"
                 << endl;
         }//if
         return NA_REAL;
      }//if

      if (fEqualVar) {
         dfr = len1 + len2 - 2;
         Double_t var = ((len1-1)*var1 + (len2-1)*var2) / dfr;
         err = TMath::Sqrt(var*(1.0/len1 + 1.0/len2));
      } else {
         Double_t v1 = var1 / len1;
         Double_t v2 = var2 / len2;
         Double_t v  = v1 + v2;
         err = TMath::Sqrt(v);
         dfr = v*v / (v1*v1/(len1-1) + v2*v2/(len2-1));
      }//if

      ts = (mn1 - mn2 - mu)/err;
   }//if

   mean1 = mn1;
   mean2 = mn2;
   se    = err;
   df    = dfr;
   return ts;
}//Statistic


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Utility functions                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
Int_t NumSep(const char *name, const char *sep)
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
}//NumSep


