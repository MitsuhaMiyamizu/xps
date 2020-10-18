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

using namespace std;

#include <new>  //needed for new (nothrow)

// ROOT
#include "TMath.h"
#include "TRandom.h"

#include "TMLMath.h"
#include "TStat.h"

//debug: print class method names
const Bool_t  kCS  = 0;
const Bool_t  kCSa = 0;

ClassImp(TStat);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TStat                                                                //
//                                                                      //
// Class containing statistical methods                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
Double_t *TStat::PMax(Int_t n, Int_t m, const Double_t **table, Double_t *max)
{
  // Returns maximum in each row of table[1:n,1:m] in array max[1:n].

   for (Int_t i=0; i<n; i++) {
      Double_t rmax = table[i][0];
      for (Int_t j=1; j<m; j++) {
         rmax = (table[i][j] > rmax) ? table[i][j] : rmax;
      }//for_j
      max[i] = rmax;
   }//for_i
   return max;
}//PMax

//______________________________________________________________________________
Double_t *TStat::PMin(Int_t n, Int_t m, const Double_t **table, Double_t *min)
{
  // Returns minimum in each row of table[1:n,1:m] in array min[1:n].

   for (Int_t i=0; i<n; i++) {
      Double_t rmin = table[i][0];
      for (Int_t j=1; j<m; j++) {
         rmin = (table[i][j] < rmin) ? table[i][j] : rmin;
      }//for_j
      min[i] = rmin;
   }//for_i
   return min;
}//PMin

//______________________________________________________________________________
void TStat::Rank(Int_t n, Double_t *arr, Int_t *index, Int_t *rank, Bool_t down)
{
   // Return sort index and rank of array arr
   if(kCSa) cout << "------TStat::Rank(int)------" << endl;

   if (n <= 0) return;
   if (n == 1) {
      index[0] = 0;
      rank[0] = 0;
      return;
   }//if

   TMath::Sort(n,arr,index,down);

   for (Int_t i=0; i<n; i++) {
      rank[index[i]] = i;
   }//for_i
}//Rank

//______________________________________________________________________________
void TStat::Rank(Int_t n, Double_t *arr, Double_t *rank)
{
   // Return true rank of array arr analogously to the R function rank(), i.e.
   // equal values in array get identical rank, however, rank starts with rank=1
   // and can be a real number! Note that array arr must already be sorted!
   // Note: Array rank can be used as index, e.g. arr[(Int_t)floor(rank[i])-1]
   if(kCSa) cout << "------TStat::Rank(double)------" << endl;

   Int_t i = 0;
   Int_t j = 0;
   
   while (i < n) {
      j = i;

      while ((j < n-1) && (arr[j] == arr[j+1])) j++;

      if (i != j) {
         for (Int_t k=i; k<=j; k++) rank[k] = (i + j + 2) / 2.0;
      } else {
         rank[i] = i + 1;
      }//if

      i = j + 1;
   }//while
}//Rank

//______________________________________________________________________________
Int_t *TStat::TrueRank(Int_t n, Double_t *arr, Int_t *rank, Bool_t down)
{
   // Return true rank of array arr, i.e. equal values in array get
   // identical rank! Note that array arr must already be sorted!
   // Warning: TrueRank() can not be used as index, use Rank() or GetRanks()
   if(kCSa) cout << "------TStat::TrueRank------" << endl;

   if (n <= 0) return 0;
   if (n == 1) {
      rank[0] = 0;
      return rank;
   }//if

   Int_t *index = new Int_t[n];
   TMath::Sort(n,arr,index,down);

   Int_t k = 0;
   for (Int_t i=0; i<n; i++) {
      if ((i > 0) && (arr[index[i]] == arr[index[i-1]])) {
         rank[index[i]] = i - 1;
         k++;
      }//if
      rank[index[i]] = i - k;
   }//for_i

   delete [] index;
   return rank;
}//TrueRank

//______________________________________________________________________________
Double_t *TStat::TrueRank1(Int_t n, Double_t *arr, Int_t *index, Double_t *rank)
{
   // Return true rank of array arr, i.e. equal values in array get
   // identical rank, however, rank starts with rank=1 and can be a real number!
   // Warning:
   // TrueRank1() can not be used as index, for this purpose use Rank()
////////////
//PROBLEM: currently onyly true for already sorted arr!!!!!
////////////
   if(kCSa) cout << "------TStat::TrueRank1------" << endl;

   if (n <= 0) return 0;
   if (n == 1) {
      index[0] = 0;
      rank[0]  = 1.0;
      return rank;
   }//if

   Double_t *sorted = new Double_t[n];

   TMath::Sort(n, arr, index, kFALSE);
   for (Int_t i=0; i<n; i++) sorted[i] = arr[index[i]];

   Int_t rnk  = 1;
   Int_t sum  = 1;
   Int_t ntie = 1;
   Int_t prev = 0;

   rank[0] = 1.0;
   for(Int_t i=1; i<=n; i++) {
      if (i < n && sorted[i] == sorted[prev]) { 
         ntie++;
         rnk++;
         sum += rnk;
      } else {
         if (ntie > 1) {
            while (prev < i) {
               rank[prev] = (Double_t)sum / (Double_t)ntie;
               prev++;
            }//while
         }//if

         rnk++;
         sum = rnk;
         if (i < n) rank[i] = rnk;
         prev = i;
         ntie = 1;
      }//if
   }//for_i

   delete [] sorted;

   return rank;
}//TrueRank1

//______________________________________________________________________________
Double_t TStat::GeoMean(Int_t n, const Double_t *arr)
{
   // Calculate geometric mean
   if(kCSa) cout << "------TStat::GeoMean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

   Double_t mean = 0;
   for (Int_t i=0; i<n; i++) mean += TMath::Log10(arr[i]);
   mean = TMath::Power(10, mean/(Double_t)n);

   return mean;
}//GeoMean

//______________________________________________________________________________
Double_t TStat::GeoMean(Int_t n, const Double_t *arr, Int_t &len, const Double_t na)
{
   // Calculate geometric mean for array containing missing values set to na
   if(kCSa) cout << "------TStat::GeoMean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return ((arr[0] != na) ? arr[0] : NA_REAL);

   Double_t mean  = 0;
   Int_t    count = n;
   for (Int_t i=0; i<n; i++) {
      (!(arr[i] == na || TMLMath::IsNaN(arr[i]))) ? (mean += TMath::Log10(arr[i])) : count--;
   }//for
   mean = (count > 0) ? TMath::Power(10, mean/(Double_t)n) : NA_REAL;

   len = count;
   return mean;
}//GeoMean

//______________________________________________________________________________
Double_t TStat::GeoMean(Int_t n, const Double_t *arr, const Double_t *w)
{
   // Calculate geometric mean with weights w
   if(kCSa) cout << "------TStat::GeoMean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

   Double_t mean = 0.0;
   Double_t sumw = 0.0;
   for (Int_t i=0; i<n; i++) {;
      sumw += w[i];
      mean += TMath::Log10(arr[i] * w[i]);
   }//for_i

   if (sumw > 0.0) {
      mean = TMath::Power(10, mean/sumw);
   } else {
      cout << "Error: Sum of weights is null!" << endl;
      mean = 0.0;   //??
//      mean = DBL_MAX;   //better: Inf?
   }//if

   return mean;
}//GeoMean

//______________________________________________________________________________
Double_t TStat::GeoMean(Int_t n, const Double_t *arr, const Double_t trim,
                Double_t &var, Int_t &len)
{
   // Calculate geometric trimmed mean (trim >= 0.5 is equal to median)
   // var is variance and len is length of trimmed array
   if(kCSa) cout << "------TStat::GeoMean(trim)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) {var = 0.0; len = 1; return arr[0];}

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index);

   // start-index and end-index
   Int_t start, end;
   if (trim < 0.5) {
      start = (Int_t)floor(n * trim);
      end   = n - start;
   } else {
      if ((n % 2) == 0){
         start = (Int_t)floor(n / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)floor(n / 2.0);
//??         start = (Int_t)ceil(n / 2);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 0.0;
   for (Int_t i=start; i<end; i++) {
      mean += TMath::Log10(arr[index[i]]);
   }//for_i
   mean = TMath::Power(10, mean/(Double_t)trimlen);

// Calculate variance
   Double_t var1 = 0;
   for (Int_t i=start; i<end; i++) {
      var1 += (arr[index[i]] - mean)*(arr[index[i]] - mean);
   }//for_i
   if (trimlen > 1) var1 /= (trimlen - 1);
   else var1 = 0;

   delete [] index;

   var = var1;
   len = trimlen;
   return mean;
}//GeoMean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t begin, Int_t end, const Double_t *arr)
{
   // Calculate arithmetic mean for array[begin, end-1]
   if(kCSa) cout << "------TStat::Mean------" << endl;

   if (begin < 0)    return NA_REAL;
   if (begin == end) return arr[begin];

   Double_t mean = 0.0;
   for (Int_t i=begin; i<end; i++) mean += arr[i];

   return mean/(end - begin);
}//Mean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t n, const Double_t *arr)
{
   // Calculate arithmetic mean
   if(kCSa) cout << "------TStat::Mean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

   Double_t mean = 0.0;
   for (Int_t i=0; i<n; i++) mean += arr[i];

   return mean/n;
}//Mean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t n, const Double_t *arr, Int_t &len, const Double_t na)
{
   // Calculate arithmetic mean for array containing missing values set to na
   if(kCSa) cout << "------TStat::Mean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return ((arr[0] != na) ? arr[0] : NA_REAL);

   Double_t mean  = 0.0;
   Int_t    count = n;
   for (Int_t i=0; i<n; i++) {
      (!(arr[i] == na || TMLMath::IsNaN(arr[i]))) ? (mean += arr[i]) : count--;
   }//for

   len = count;
   return (count > 0) ? mean/count : NA_REAL;
}//Mean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t n, const Double_t *arr, const Double_t *w)
{
   // Calculate arithmetic mean with weights w
   if(kCSa) cout << "------TStat::Mean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

   Double_t mean = 0.0;
   Double_t sumw = 0.0;
   for (Int_t i=0; i<n; i++) {;
      sumw += w[i];
      mean += arr[i] * w[i];
   }//for_i

   if (sumw > 0.0) {
      mean /= sumw;
   } else {
      cout << "Error: Sum of weights is null!" << endl;
      mean = NA_REAL;
   }//if

   return mean;
}//Mean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t n, const Double_t *arr, const Double_t trim)
{
   // Calculate arithmetic trimmed mean (trim >= 0.5 is equal to median)
   if(kCSa) cout << "------TStat::Mean(trim)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index);

   // start-index and end-index
   Int_t start, end;
   if (trim < 0.5) {
      start = (Int_t)floor(n * trim);
      end   = n - start;
   } else {
      if ((n % 2) == 0){
         start = (Int_t)floor(n / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)floor(n / 2.0);
//??         start = (Int_t)ceil(n / 2);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 0.0;
   for (Int_t i=start; i<end; i++) {
      mean += arr[index[i]];
   }//for_i
   mean /= trimlen;

   delete [] index;

   return mean;
}//Mean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t n, const Double_t *arr, const Double_t trim,
                Double_t &var, Int_t &len)
{
   // Calculate arithmetic trimmed mean (trim >= 0.5 is equal to median)
   // var is variance and len is length of trimmed array
   if(kCSa) cout << "------TStat::Mean(trim)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) {var = 0.0; len = 1; return arr[0];}

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index);

   // start-index and end-index
   Int_t start, end;
   if (trim < 0.5) {
      start = (Int_t)floor(n * trim);
      end   = n - start;
   } else {
      if ((n % 2) == 0){
         start = (Int_t)floor(n / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)floor(n / 2.0);
//??         start = (Int_t)ceil(n / 2);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 0.0;
   for (Int_t i=start; i<end; i++) {
      mean += arr[index[i]];
   }//for_i
   mean /= trimlen;

// Calculate variance
   Double_t var1 = 0;
   Double_t tmp  = 0;
   for (Int_t i=start; i<end; i++) {
      tmp = arr[index[i]] - mean;
      var1 += tmp * tmp;
   }//for_i
   if (trimlen > 1) var1 /= (trimlen - 1);
   else var1 = 0;

   delete [] index;

   var = var1;
   len = trimlen;
   return mean;
}//Mean

//______________________________________________________________________________
Double_t TStat::Mean(Int_t n, const Double_t *arr, const Double_t *w,
                const Double_t trim, Double_t &var, Int_t &len)
{
   // Calculate arithmetic trimmed mean with weights w
   // var is variance and len is length of trimmed array
   if(kCSa) cout << "------TStat::Mean------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) {var = 0; len = 1; return arr[0];}

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index);

   // start-index and end-index
   Int_t start, end;
   if (trim < 0.5) {
      start = (Int_t)floor(n * trim);
      end   = n - start;
   } else {
      if ((n % 2) == 0){
         start = (Int_t)floor(n / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)floor(n / 2.0);
//??         start = (Int_t)ceil(n / 2);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 0.0;
   Double_t sumw    = 0.0;
   for (Int_t i=start; i<end; i++) {
      sumw += w[index[i]];
      mean += arr[index[i]] * w[index[i]];
   }//for_i

// Calculate variance
   Double_t var1 = 0.0;
   Double_t tmp  = 0.0;
   if (sumw > 0.0) {
      mean /= sumw;
      for (Int_t i=start; i<end; i++) {
         tmp = arr[index[i]] - mean;
         var1 += tmp * tmp * w[index[i]];
      }//for_i
//??      if (trimlen > 1) var1 /= (trimlen - 1) * sumw;
      if (trimlen > 1) var1 /= (trimlen - 1);
      else var1 = 0;
   } else {
      cout << "Error: Sum of weights is null!" << endl;
      mean = NA_REAL; 
      var1 = NA_REAL;
   }//if

   delete [] index;

   var = var1;
   len = trimlen;
   return mean;
}//Mean

//______________________________________________________________________________
Double_t TStat::Median(Int_t n, const Double_t *arr)
{
   // Calculate median
   if(kCSa) cout << "------TStat::Median(n)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index);

// Find median
   Int_t k;
   Double_t median = 0.0;
   if ((n % 2) == 0){
      k = (Int_t)floor(n / 2.0) - 1;
      median = (arr[index[k]] + arr[index[k+1]])/2.0;
   } else {
      k = (Int_t)floor(n / 2.0);
      median = arr[index[k]];
   }//if

   delete [] index;

   return median;
}//Median

//______________________________________________________________________________
Double_t TStat::Median(Int_t n, const Double_t *arr, const Int_t *index)
{
   // Calculate median of array "arr" with sort index "index"
   if(kCSa) cout << "------TStat::Median(n)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

   Int_t k;
   Double_t median = 0.0;
   if ((n % 2) == 0){
      k = (Int_t)floor(n / 2.0) - 1;
      median = (arr[index[k]] + arr[index[k+1]])/2.0;
   } else {
      k = (Int_t)floor(n / 2.0);
      median = arr[index[k]];
   }//if

   return median;
}//Median

//______________________________________________________________________________
Double_t TStat::Median(Int_t n, const Double_t *arr, Bool_t low, Bool_t high)
{
   // Calculate median
   if(kCSa) cout << "------TStat::Median(lo,hi)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return arr[0];

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index, kFALSE);

// Find median
   Int_t k;
   Double_t median = 0.0;
   if ((n % 2) == 0){
      k = (Int_t)floor(n / 2.0) - 1;

      if      (low && !high) median = arr[index[k]];
      else if (high && !low) median = arr[index[k+1]];
      else                   median = (arr[index[k]] + arr[index[k+1]])/2.0;
   } else {
      k = (Int_t)floor(n / 2.0);
      median = arr[index[k]];
   }//if

   delete [] index;

   return median;
}//Median

//______________________________________________________________________________
Double_t TStat::Median(Int_t n, const Double_t *arr, UShort_t logbase)
{
   // Calculate median
   // Convert first to logbase: 0 - linear, 1 - log, 2 - log2, 10 - log10
   if(kCSa) cout << "------TStat::Median(logbase)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) {
      if (logbase == 0)  {return arr[0];              } else
      if (logbase == 1)  {return TMath::Log(arr[0]);  } else
      if (logbase == 2)  {return TMath::Log2(arr[0]); } else
      if (logbase == 10) {return TMath::Log10(arr[0]);} 
   }//if

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index);

// Find median
   Int_t k;
   Double_t median = 0.0;
   if ((n % 2) == 0){
      k = (Int_t)floor(n / 2.0) - 1;
      if (logbase == 0)  {median =  (arr[index[k]]               + arr[index[k+1]])/2.0;              } else
      if (logbase == 1)  {median =  (TMath::Log(arr[index[k]])   + TMath::Log(arr[index[k+1]]))/2.0;  } else
      if (logbase == 2)  {median =  (TMath::Log2(arr[index[k]])  + TMath::Log2(arr[index[k+1]]))/2.0; } else
      if (logbase == 10) {median =  (TMath::Log10(arr[index[k]]) + TMath::Log10(arr[index[k+1]]))/2.0;} 
   } else {
      k = (Int_t)floor(n / 2.0);
      if (logbase == 0)  {median =  arr[index[k]];              } else
      if (logbase == 1)  {median =  TMath::Log(arr[index[k]]);  } else
      if (logbase == 2)  {median =  TMath::Log2(arr[index[k]]); } else
      if (logbase == 10) {median =  TMath::Log10(arr[index[k]]);} 
   }//if

   delete [] index;

   return median;
}//Median

//______________________________________________________________________________
Double_t TStat::MAD(Int_t n, const Double_t *arr, Float_t constant)
{
   // Calculate median absolute deviation
   if(kCSa) cout << "------TStat::MAD(n)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return 0;

   Double_t *x = 0;
   if (!(x = new (nothrow) Double_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if

   Double_t center = TStat::Median(n, arr);

   for (Int_t i=0; i<n; i++) {
      x[i] = TMath::Abs(arr[i] - center);
   }//for_i

   Double_t mad = TStat::Median(n, x);

   delete [] x;

   return constant * mad;
}//MAD

//______________________________________________________________________________
Double_t TStat::MAD(Int_t n, const Double_t *arr, const Double_t trim, Float_t constant)
{
   // Calculate mean/median absolute deviation
   if(kCSa) cout << "------TStat::MAD(trim)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return 0;

   Double_t *x = 0;
   if (!(x = new (nothrow) Double_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if

   Double_t center = TStat::Mean(n, arr, trim);

   for (Int_t i=0; i<n; i++) {
      x[i] = TMath::Abs(arr[i] - center);
   }//for_i

   Double_t mad = TStat::Median(n, x);

   delete [] x;

   return constant * mad;
}//MAD

//______________________________________________________________________________
Double_t TStat::MAD(Int_t n, const Double_t *arr, Double_t center, Double_t constant,
                    Bool_t low, Bool_t high)
{
   // Calculate mean/median absolute deviation
   if(kCSa) cout << "------TStat::MAD(center)------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return 0;

   Double_t *x = 0;
   if (!(x = new (nothrow) Double_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if

   for (Int_t i=0; i<n; i++) {
      x[i] = TMath::Abs(arr[i] - center);
   }//for_i

   Double_t mad = TStat::Median(n, x, low, high);

   delete [] x;

   return constant * mad;
}//MAD

//______________________________________________________________________________
Double_t TStat::Var(Int_t n, const Double_t *arr, const Double_t mean)
{
   // Calculate variance for arithmetic mean
   if(kCSa) cout << "------TStat::Var------" << endl;

   if (n <= 0) return NA_REAL;
   if (n == 1) return 0;

   Double_t var = 0.0;
   Double_t tmp = 0.0;
   for (Int_t i=0; i<n; i++) {
      tmp = arr[i] - mean;
      var += tmp * tmp;
   }//for_i

   return var / (n - 1);
}//Var

//______________________________________________________________________________
Double_t TStat::Var(Int_t n, const Double_t *arr, const Double_t mean,
                Int_t &len, const Double_t na)
{
   // Calculate variance for arithmetic mean
   // array arr can contain missing values set to na
   if(kCSa) cout << "------TStat::Var------" << endl;

   if (TMLMath::IsNaN(mean) || !TMLMath::Finite(mean) || (mean == na)) {
      return NA_REAL;
   }//if
   if (n == 1) return 0;

   Double_t var   = 0.0;
   Double_t tmp   = 0.0;
   Int_t    count = n;
   for (Int_t i=0; i<n; i++) {
      if (!(arr[i] == na || TMLMath::IsNaN(arr[i]))) {
         tmp = arr[i] - mean;
         var += tmp * tmp;
      } else {
         count--;
      }//if
   }//for_i

   len = count;
   return (count > 1) ? var/(count-1) : NA_REAL;
}//Var

//______________________________________________________________________________
Double_t TStat::MedianPolish(Int_t nrow, Int_t ncol, Double_t *x, Double_t *rowmed,
                Double_t *colmed, Double_t *residu, Int_t maxiter, Double_t eps,
                Bool_t verbose)
{
   // Compute median polish algorithm for a table of size nrow*ncol, whereby
   // the table is imported as array x[ij] with ij = i*ncol + j
   // Return value is the total median; 
   // row medians and column medians are returned in arrays rowmed and colmed
   // residuals are returned in array residu of size nrow*ncol
   // Median polish alternately removes the row and column medians until the 
   // proportional reduction in the sum of absolute residuals is less than eps
   // or until there have been maxiter iterations. 
   if(kCSa) cout << "------TStat::MedianPolish------" << endl;

   Int_t i, j, k, iter;
   Int_t size = nrow*ncol;

// Init return values
   for (i=0; i<nrow; i++) rowmed[i] = 0;
   for (j=0; j<ncol; j++) colmed[j] = 0;
   for (k=0; k<size; k++) residu[k] = x[k];

   Double_t totmed = 0.0;
   Double_t delta  = 0.0;
   Double_t oldsum = 0.0;
   Double_t newsum = 0.0;

   Double_t *array  = new Double_t[nrow >= ncol ? nrow : ncol]; 
   Double_t *rdelta = new Double_t[nrow]; 
   Double_t *cdelta = new Double_t[ncol]; 

// Loop until converged or maxiter is reached  
   for (iter=1; iter<=maxiter; iter++){
      // get medians of rows
      for (i=0; i<nrow; i++) { 
         for (j=0; j<ncol; j++) array[j] = residu[i*ncol + j];
//         rdelta[i] = TStat::Median(ncol, array);
         rdelta[i] = TMath::Median(ncol, array);
      }//for_i

      // subtract row medians from each row
      for (i=0; i<nrow; i++) { 
         for (j=0; j<ncol; j++) residu[i*ncol + j] -= rdelta[i];
      }//for_i

      // add rdelta to rowmed
      for (i=0; i<nrow; i++) rowmed[i] += rdelta[i];

      // subtract median of colmed from colmed
//      delta = TStat::Median(ncol, colmed);
      delta = TMath::Median(ncol, colmed);
      for (j=0; j<ncol; j++) colmed[j] -= delta;
      totmed = totmed + delta;

      // get medians of columns
      for (j=0; j<ncol; j++){
         for (i=0; i<nrow; i++) array[i] = residu[i*ncol + j];
//         cdelta[j] = TStat::Median(nrow, array);
         cdelta[j] = TMath::Median(nrow, array);
      }//for_j

      // subtract column medians from each column
      for (j=0; j<ncol; j++){
         for (i=0; i<nrow; i++) residu[i*ncol + j] -= cdelta[j];
      }//for_j

      // add cdelta to colmed
      for (j=0; j<ncol; j++) colmed[j] += cdelta[j];

      // subtract median of rowmed from rowmed
//      delta = TStat::Median(nrow, rowmed);
      delta = TMath::Median(nrow, rowmed);
      for (i=0; i<nrow; i++) rowmed[i] -= delta;
      totmed = totmed + delta;

      // sum of the absolute values of all elements
      newsum = 0.0;
//??      for (k=0; k<size; k++) newsum += residu[k];
      for (k=0; k<size; k++) newsum += TMath::Abs(residu[k]);

      // test for convergence
      if (newsum == 0.0 || fabs(newsum - oldsum) < eps*newsum) break;
      oldsum = newsum;
   }//for_iter

   if (verbose && (iter >= maxiter)) {
      cout << "Warning: MedianPolish did not converge in <" << maxiter
           << "> iterations." << endl;
   }//if

// Cleanup
   delete [] array;
   delete [] cdelta;
   delete [] rdelta;

   return totmed;
}//MedianPolish

//______________________________________________________________________________
Double_t TStat::MedianPolishTranspose(Int_t nrow, Int_t ncol, Double_t *x, 
                Double_t *rowmed, Double_t *colmed, Double_t *residu,
                Int_t maxiter, Double_t eps, Bool_t verbose)
{
   // Compute median polish algorithm for a table of size nrow*ncol, whereby
   // the table is imported as array x[ij] with ij = i*ncol + j
   // Return value is the total median; 
   // row medians and column medians are returned in arrays rowmed and colmed
   // residuals are returned in array residu of size nrow*ncol
   // Note for transposed algorithm:
   // Median polish alternately removes the column and row medians, starting
   // with removal of columns, until the proportional reduction in the sum of 
   // absolute residuals is less than eps or until there have been maxiter iterations. 
   if(kCSa) cout << "------TStat::MedianPolishTranspose------" << endl;

   Int_t i, j, k, iter;
   Int_t size = nrow*ncol;

// Init return values
   for (i=0; i<nrow; i++) rowmed[i] = 0;
   for (j=0; j<ncol; j++) colmed[j] = 0;
   for (k=0; k<size; k++) residu[k] = x[k];

   Double_t totmed = 0.0;
   Double_t delta  = 0.0;
   Double_t oldsum = 0.0;
   Double_t newsum = 0.0;

   Double_t *array  = new Double_t[nrow >= ncol ? nrow : ncol]; 
   Double_t *rdelta = new Double_t[nrow]; 
   Double_t *cdelta = new Double_t[ncol]; 

// Loop until converged or maxiter is reached  
   for (iter=1; iter<=maxiter; iter++){
      // get medians of columns
      for (j=0; j<ncol; j++){
         for (i=0; i<nrow; i++) array[i] = residu[i*ncol + j];
//         cdelta[j] = TStat::Median(nrow, array);
         cdelta[j] = TMath::Median(nrow, array);
      }//for_j

      // subtract column medians from each column
      for (j=0; j<ncol; j++){
         for (i=0; i<nrow; i++) residu[i*ncol + j] -= cdelta[j];
      }//for_j

      // add cdelta to colmed
      for (j=0; j<ncol; j++) colmed[j] += cdelta[j];

      // subtract median of rowmed from rowmed
//      delta = TStat::Median(nrow, rowmed);
      delta = TMath::Median(nrow, rowmed);
      for (i=0; i<nrow; i++) rowmed[i] -= delta;
      totmed = totmed + delta;

      // get medians of rows
      for (i=0; i<nrow; i++) { 
         for (j=0; j<ncol; j++) array[j] = residu[i*ncol + j];
//         rdelta[i] = TStat::Median(ncol, array);
         rdelta[i] = TMath::Median(ncol, array);
      }//for_i

      // subtract row medians from each row
      for (i=0; i<nrow; i++) { 
         for (j=0; j<ncol; j++) residu[i*ncol + j] -= rdelta[i];
      }//for_i

      // add rdelta to rowmed
      for (i=0; i<nrow; i++) rowmed[i] += rdelta[i];

      // subtract median of colmed from colmed
//      delta = TStat::Median(ncol, colmed);
      delta = TMath::Median(ncol, colmed);
      for (j=0; j<ncol; j++) colmed[j] -= delta;
      totmed = totmed + delta;

      // sum of the absolute values of all elements
      newsum = 0.0;
//??      for (k=0; k<size; k++) newsum += residu[k];
      for (k=0; k<size; k++) newsum += TMath::Abs(residu[k]);

      // test for convergence
      if (newsum == 0.0 || fabs(newsum - oldsum) < eps*newsum) break;
      oldsum = newsum;
   }//for_iter

   if (verbose && (iter >= maxiter)) {
      cout << "Warning: MedianPolish did not converge in <" << maxiter
           << "> iterations." << endl;
   }//if

// Cleanup
   delete [] array;
   delete [] cdelta;
   delete [] rdelta;

   return totmed;
}//MedianPolishTranspose

//______________________________________________________________________________
Double_t TStat::TukeyBiweight(Int_t n, const Double_t *x, Double_t &var,
         Double_t c, Double_t eps)
{
   // Calculate one-step Tukey Biweight for array x of length n
   // c is a tuning constant and eps is a small value to prevent zero division
   // Variable var will return the square of S(bi), i.e. of the measure of
   // uncertainty of the biweight. The t-distribution can be used to determine
   // the confidence interval: T(bi) +/- t*sqrt(var/n)
   if(kCSa) cout << "------TStat::TukeyBiweight------" << endl;

   if (n == 1) return x[0];

   Double_t m, s, cseps, uu2;
   Double_t sumwx = 0.0;
   Double_t sumw  = 0.0;
   Double_t tb    = 0.0;

// Init local arrays
   Double_t *u = 0;
   Double_t *w = 0;
   if (!(u = new (nothrow) Double_t[n])) goto cleanup;
   if (!(w = new (nothrow) Double_t[n])) goto cleanup;
   for (Int_t i=0; i<n; i++) u[i] = w[i] = 0.0;

   m = TStat::Median(n, x);
   for (Int_t i=0; i<n; i++) {
      u[i] = TMath::Abs(x[i] - m);
   }//for_i
   s = TStat::Median(n, u);

// Calculate one-step tukey biweight
   sumwx = 0.0;
   sumw  = 0.0;
   cseps = c*s + eps;
   for (Int_t i=0; i<n; i++) {
      u[i] = (x[i] - m) / cseps;
      uu2  = 1 - u[i]*u[i];
      w[i] = (TMath::Abs(u[i]) <= 1.0) ? uu2*uu2 : 0.0;
      sumwx += w[i]*x[i];
      sumw  += w[i];
   }//for_i
   tb = sumwx / sumw;

// Calculate variance for tukey biweight
   sumwx = 0.0;
   sumw  = 0.0;
   for (Int_t i=0; i<n; i++) {
      uu2  = (x[i] - tb)*w[i];
      sumwx += uu2*uu2;
      sumw  += TMath::Sqrt(w[i])*(1 - 5*u[i]*u[i]);
   }//for_i
   sumw = TMath::Abs(sumw);
   var = n*sumwx / (sumw*sumw);

// Cleanup
cleanup:
   if (w) {delete [] w; w = 0;}
   if (u) {delete [] u; u = 0;}

   return tb;
}//TukeyBiweight

//______________________________________________________________________________
Double_t TStat::Quantile(Int_t n, const Double_t *arr, const Double_t q)
{
   // Return value of array arr at quantile q
   if(kCSa) cout << "------TStat::Quantile------" << endl;

   if (n == 1) return arr[0];

   if (q < 0.0 || q > 1.0) {
      cout << "Error: Quantile q is not within [0,1]!" << endl;
      return NA_REAL;
   }//if

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index, kFALSE);

// Find quantile
   Double_t qu = (n - 1) * q;
   Int_t    lo = (Int_t)floor(qu);
   Int_t    hi = (Int_t)ceil(qu);

   Double_t ql = arr[index[lo]];
   Double_t qh = arr[index[hi]];
   Double_t qq = (ql == qh) ? 0.0 : qh - ql;

   delete [] index;

   return ql + qq * (qu - lo);
}//Quantile

//______________________________________________________________________________
Double_t TStat::Quantile(Int_t n, const Double_t *arr, const Int_t *index, const Double_t q)
{
   // Return value of array "arr" with sort index "index" at quantile q
   if(kCSa) cout << "------TStat::Quantile------" << endl;

   if (n == 1) return arr[0];

   if (q < 0.0 || q > 1.0) {
      cout << "Error: Quantile q is not within [0,1]!" << endl;
      return NA_REAL;
   }//if

   Double_t qu = (n - 1) * q;
   Int_t    lo = (Int_t)floor(qu);
   Int_t    hi = (Int_t)ceil(qu);
   Double_t ql = arr[index[lo]];
   Double_t qh = arr[index[hi]];
   Double_t qq = (ql == qh) ? 0.0 : qh - ql;

   return ql + qq * (qu - lo);
}//Quantile

//______________________________________________________________________________
Double_t TStat::IQR(Int_t n, const Double_t *arr, const Double_t qlo, const Double_t qhi)
{
   // Return interquantile range of array arr between quantiles qlo and qhi
   // Default setting returns interquartile range!
   if(kCSa) cout << "------TStat::IQR------" << endl;

   if (n == 1) return 0.0;

   if (qlo < 0.0 || qlo > 1.0 || qhi < 0.0 || qhi > 1.0) {
      cout << "Error: Quantile qlo or qhi is not within [0,1]!" << endl;
      return NA_REAL;
   }//if

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return NA_REAL;
   }//if
   TMath::Sort(n, arr, index, kFALSE);

// Find quantile
   Double_t qu, ql, qh, qq, qs;
   Int_t    lo, hi;

   qu = (n - 1) * qhi;
   lo = (Int_t)floor(qu);
   hi = (Int_t)ceil(qu);
   ql = arr[index[lo]];
   qh = arr[index[hi]];
   qq = (ql == qh) ? 0.0 : qh - ql;
   qs = ql + qq * (qu - lo);

   qu = (n - 1) * qlo;
   lo = (Int_t)floor(qu);
   hi = (Int_t)ceil(qu);
   ql = arr[index[lo]];
   qh = arr[index[hi]];
   qq = (ql == qh) ? 0.0 : qh - ql;
   qs = qs - ql - qq * (qu - lo);

   delete [] index;

   return qs;
}//IQR

//______________________________________________________________________________
Double_t TStat::IQR(Int_t n, const Double_t *arr, const Int_t *index,
                    const Double_t qlo, const Double_t qhi)
{
   // Return interquantile range of array "arr" with sort index "index" between 
   // quantiles qlo and qhi. Default setting returns interquartile range!
   if(kCSa) cout << "------TStat::IQR------" << endl;

   if (n == 1) return 0.0;

   if (qlo < 0.0 || qlo > 1.0 || qhi < 0.0 || qhi > 1.0) {
      cout << "Error: Quantile qlo or qhi is not within [0,1]!" << endl;
      return NA_REAL;
   }//if

   Double_t qu, ql, qh, qq, qs;
   Int_t    lo, hi;

   qu = (n - 1) * qhi;
   lo = (Int_t)floor(qu);
   hi = (Int_t)ceil(qu);
   ql = arr[index[lo]];
   qh = arr[index[hi]];
   qq = (ql == qh) ? 0.0 : qh - ql;
   qs = ql + qq * (qu - lo);

   qu = (n - 1) * qlo;
   lo = (Int_t)floor(qu);
   hi = (Int_t)ceil(qu);
   ql = arr[index[lo]];
   qh = arr[index[hi]];
   qq = (ql == qh) ? 0.0 : qh - ql;
   qs = qs - ql - qq * (qu - lo);

   return qs;
}//IQR

//______________________________________________________________________________
Double_t *TStat::Quantiles(Int_t n, Double_t *arr, Int_t *index,
                           Int_t nquant, Double_t *q, Double_t *quant)
{
   // Return array quant of nquant quantiles q for array arr of size n
   if(kCSa) cout << "------TStat::Quantiles------" << endl;

   TMath::Sort(n, arr, index, kFALSE);

   for (Int_t i=0; i<nquant; i++) {
      quant[i] = TStat::Quantile(n, arr, index, q[i]);
   }//for_i

   return quant;
}//Quantiles

//______________________________________________________________________________
void TStat::NextPerm(Int_t n, Int_t k, Int_t *grp) 
{
   // Return next permutation of array grp
   // Input:  ordered subset "a" (=grp) {a0, a1, ..., a(k-1)} of size k 
   //         from set {0,1,...,n-1} with a0 < a1 < ... < a(k-1)
   // Output: ordered subset "b" (=grp) {b0, b1, ..., b(k-1)} of size k,  
   //         which immediately follows subset a lexicographic order
   // Can be used in for-loop to determine all subsets of size k
   if(kCSa) cout << "------TStat::NextPerm------" << endl;

   Int_t idx = 0;
  
   Int_t i = k - 1;
   Int_t m = n - 1;
   while (grp[i] == m) {
      i--;
      m--;
      idx++;
   }//while

   i = k - 1;
   m = n - 1;
   while (grp[i] == m) {
      grp[i] = grp[i-idx] + 1 + idx;
      i--;
      m--;
      idx--;
   }//while
   
   if (i >= 0) grp[i] += 1;
}//NextPerm

//______________________________________________________________________________
void TStat::NextPerm(Int_t n, Int_t k, Int_t *grp1, Int_t nk, Int_t *grp2) 
{
   // Return next permutation of array grp1 and its complement grp2
   // Can be used in for-loop to determine all subsets of size k
   if(kCSa) cout << "------TStat::NextPerm------" << endl;

// Next permutation for grp1
   NextPerm(n, k, grp1);

// Place complement in grp2
   Int_t idx = 0; 
   for (Int_t j=0; j<n; j++) {
      Int_t tmp = 1;
      for (Int_t i=0; i<k; i++) {
         if (grp1[i] == j) {
            tmp = 0;
            break;
         }//if
      }//for_i

      if (tmp == 1) {
         grp2[idx] = j;
         idx++;
      }//if
   }//for_j

// Check if array grp2 is of correct size
   if (idx != nk) {
      cerr << "Error: Array grp2 is not of size <" << idx << ">!" << endl;
   }//if
}//NextPerm

//______________________________________________________________________________
void TStat::Sample(Int_t n, Int_t k, Int_t *grp) 
{
   // Takes a sample of the specified size k from set {0, 1,..., n-1}
   // and places it into grp
   if(kCSa) cout << "------TStat::Sample------" << endl;

   if (n <= 1) {grp[0] = 1; return;}

// Create set
   Int_t *set = 0;
   if (!(set = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return;
   }//if
   for (Int_t i=0; i<n; i++) set[i] = i;

// Create random subset and place in grp
   Int_t i, tmp; 
   for (Int_t j=0; j<k; j++) {
      i = n - j;
      while (i == n-j) {
         i = (Int_t)floor((n-j)*(gRandom->Rndm()));
      }//while
      grp[j] = set[i];
      tmp = set[n-j-1];
      set[n-j-1] = set[i];
      set[i] = tmp;
   }//for_j

   delete [] set;
}//Sample

//______________________________________________________________________________
void TStat::Sample(Int_t n, Int_t k, Int_t *grp1, Int_t nk, Int_t *grp2) 
{
   // Take a sample of the specified size k from set {0, 1,..., n-1}
   // and places it into grp1 and its complement in grp2
   if(kCSa) cout << "------TStat::Sample------" << endl;

   if (n <= 2) {grp1[0] = 0; grp2[0] = 1; return;}

// Create set
   Int_t *set = 0;
   if (!(set = new (nothrow) Int_t[n])) {
      cout << "Error: Could not initialize memory!" << endl;
      return;
   }//if
   for (Int_t i=0; i<n; i++) set[i] = i;

// Create random subset and place in grp1
   Int_t i, tmp; 
   for (Int_t j=0; j<k; j++) {
      i = n - j;
      while (i == n-j) {
         i = (Int_t)floor((n-j)*(gRandom->Rndm()));
      }//while
      grp1[j] = set[i];
      tmp = set[n-j-1];
      set[n-j-1] = set[i];
      set[i] = tmp;
   }//for_j

// Place complement in grp2
   Int_t idx = 0; 
   for (Int_t j=0; j<n; j++) {
      tmp = 1;
      for (Int_t i=0; i<k; i++) {
         if (grp1[i] == j) {
            tmp = 0;
            break;
         }//if
      }//for_i

      if (tmp == 1) {
         grp2[idx] = j;
         idx++;
      }//if
   }//for_j

// Check if array grp2 is of correct size
   if (idx != nk) {
      cerr << "Error: Array grp2 is not of size <" << idx << ">!" << endl;
   }//if

   delete [] set;
}//Sample

//______________________________________________________________________________
Double_t TStat::CorPearson(Int_t n, Double_t *x, Double_t *y)
{
   // Pearson correlation coefficient between vectors x and y

   Double_t sx, sy, sxx, syy, sxy;
   sx = sy = sxx = syy = sxy = 0.0;
   for (Int_t i=0; i<n; i++) {
      sx  += x[i];
      sy  += y[i];
      sxx += x[i] * x[i];
      syy += y[i] * y[i];
      sxy += x[i] * y[i];
   }//for

   return (sxy - sx*sy/n)/TMath::Sqrt((sxx - sx*sx/n)*(syy - sy*sy/n));
}//CorPearson

//______________________________________________________________________________
Double_t TStat::CorSpearman(Int_t n, Double_t *x, Double_t *y)
{
   // Spearman correlation coefficient between vectors x and y

   Int_t *index = new Int_t[n];
   Int_t *rankx = new Int_t[n];
   Int_t *ranky = new Int_t[n];

   TStat::Rank(n, x, index, rankx, kTRUE);
   TStat::Rank(n, y, index, ranky, kTRUE);

   Double_t d, sdd;
   d = sdd = 0.0;
   for (Int_t i=0; i<n; i++) {
      d = rankx[i] - ranky[i];
      sdd += d * d;
   }//for

   delete [] ranky;
   delete [] rankx;
   delete [] index;

   return (1 - 6*sdd/(n*(n*n - 1)));
}//CorSpearman

//______________________________________________________________________________
void TStat::MassDist(Int_t nx, Double_t *x, Double_t *w, Double_t xmin, Double_t xmax,
                     Int_t ny, Double_t *y)
{
   // Discretize data x with xmin and xmax minimum and maximum value of x
   // Array w of length nx is weight for x
   // Output array y contains discretation scheme of data
   // see Applied Statistics R50 and Applied Statistics 176
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::MassDist------" << endl;
  
   Int_t    ixmin = 0;
   Int_t    ixmax = ny - 2;
   Double_t mass  = 0.0;
   Double_t delta = (xmax - xmin) / (ny - 1);
  
   for (Int_t i=0; i<ny; i++) y[i] = 0.0;

   for (Int_t i=0; i<nx; i++) mass += w[i];
   mass = 1.0/mass;

   for(Int_t i=0; i<nx; i++) {
      if (TMath::Finite(x[i])) {
         Double_t pos = (x[i] - xmin) / delta;
         Int_t    ix  = (Int_t)TMath::Floor(pos);
         Double_t fx  = pos - ix;

         if ((ixmin <= ix) && (ix <= ixmax)) {
            y[ix]     += w[i]*(1 - fx);
            y[ix + 1] += w[i]*fx;
         } else if (ix == -1) {
            y[0] += w[i]*fx;
         } else if (ix == ixmax + 1) {
            y[ix] += w[i]*(1 - fx);
         }//if
      }//if
   }//for_i
  
   for (Int_t i=0; i<ny; i++) y[i] *= mass;
}//MassDist

//______________________________________________________________________________
void TStat::TwiddleFactor4FFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im)
{
   // Twiddle factor for FFT
   // n     - length of data series
   // tf_re - on output contains real part of twiddle factor
   // tf_im - on output contains imaginary part of twiddle factor
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::TwiddleFactor4FFT------" << endl;

   if (i ==0 ) {
      tf_re = 1.0;
      tf_im = 0.0;
   } else {
      tf_re =  TMath::Cos(2*TMath::Pi()*(Double_t)i / (Double_t)n);  
      tf_im = -TMath::Sin(2*TMath::Pi()*(Double_t)i / (Double_t)n); 
   }//if
}//TwiddleFactor4FFT

//______________________________________________________________________________
void TStat::TwiddleFactor4IFFT(Int_t n, Int_t i, Double_t &tf_re, Double_t &tf_im)
{
   // Twiddle factor for Inverse FFT
   // n     - length of data series
   // tf_re - on output contains real part of twiddle factor
   // tf_im - on output contains imaginary part of twiddle factor
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::TwiddleFactor4IFFT------" << endl;

   if (i == 0) {
      tf_re = 1.0;
      tf_im = 0.0;
   } else {
      tf_re = TMath::Cos(2*TMath::Pi()*(Double_t)i / (Double_t)n);  
      tf_im = TMath::Sin(2*TMath::Pi()*(Double_t)i / (Double_t)n); 
   }//if
}//TwiddleFactor4IFFT

//______________________________________________________________________________
void TStat::FFT(Int_t n, Double_t *f_re, Double_t *f_im)
{
   // Fast Fourier Transform in space, result is in reverse bit order.
   // compute FFT using Decimation In Frequency of a data sequence of length 2^n
   // f_re - real component of data series
   // f_im - imaginary component of data series
   // n    -  where 2^n is length of data series
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::FFT------" << endl;
  
   Int_t    baseE, baseO, blocks, points, points2;
   Double_t even_re, even_im, odd_re, odd_im, tf_re, tf_im;

   blocks = 1;
   points = 1 << n;

   for (Int_t i=0; i<n; i++) {
      points2 = points >> 1;
      baseE = 0;

      for (Int_t j=0; j<blocks; j++) {
         baseO = baseE + points2;

         for (Int_t k=0; k<points2; k++) {
            even_re = f_re[baseE + k] + f_re[baseO + k]; 
            even_im = f_im[baseE + k] + f_im[baseO + k];

            TStat::TwiddleFactor4FFT(points, k, tf_re, tf_im);

            odd_re = (f_re[baseE + k] - f_re[baseO + k])*tf_re
                   - (f_im[baseE + k] - f_im[baseO + k])*tf_im;
            odd_im = (f_re[baseE + k] - f_re[baseO + k])*tf_im
                   + (f_im[baseE + k] - f_im[baseO + k])*tf_re; 

            f_re[baseE + k] = even_re;
            f_im[baseE + k] = even_im;
            f_re[baseO + k] = odd_re;
            f_im[baseO + k] = odd_im;
         }//for_k

         baseE = baseE + points;
      }//for_j 
                   
      blocks = blocks << 1; 
      points = points >> 1;
   }//for_i
}//FFT

//______________________________________________________________________________
void TStat::IFFT(Int_t n, Double_t *f_re, Double_t *f_im)
{
   // Inverse Fast Fourier Transform in space, where input is in reverse bit
   // order and output is in normal order.
   // compute IFFT using Decimation In Time of a data sequence of length 2^n
   // f_re - real component of data series
   // f_im - imaginary component of data series
   // n    -  where 2^n is length of data series
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::IFFT------" << endl;

   Int_t    blocks, points, points2, baseB, baseT;
   Double_t top_re, top_im, bot_re, bot_im, tf_re, tf_im;

   blocks = 1 << (n-1);
   points = 2;  
   for (Int_t i=0; i<n; i++) {
      points2 = points >> 1;
      baseT = 0;

      for (Int_t j=0; j<blocks; j++) {
         baseB = baseT + points2;

         for (Int_t k=0; k<points2; k++) {
            top_re = f_re[baseT + k];
            top_im = f_im[baseT + k];
	
            TStat::TwiddleFactor4IFFT(points, k, tf_re, tf_im);

            bot_re = f_re[baseB + k]*tf_re - f_im[baseB + k]*tf_im;
            bot_im = f_re[baseB + k]*tf_im + f_im[baseB + k]*tf_re;

            f_re[baseT + k] = top_re + bot_re;
            f_im[baseT + k] = top_im + bot_im;
            f_re[baseB + k] = top_re - bot_re;   
            f_im[baseB + k] = top_im - bot_im; 
         }//for_k  
  
         baseT = baseT + points; 
      }//for_j

      blocks = blocks >> 1;
      points = points << 1;
   }//for_i
}//IFFT

//______________________________________________________________________________
Int_t TStat::FFTDensityConvolve(Int_t n, Double_t *x_re, Double_t *y_re)
{
   // Fast Fourier Transform density convolve
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::FFTDensityConvolve------" << endl;

   Int_t err = 0;

   // to stop rounding problems 
   Int_t nlog2 = (Int_t)(TMath::Log((Double_t)n) / TMath::Log(2.0) + 0.5);

// Init local arrays
   Double_t *x_im = 0;
   Double_t *y_im = 0;
   Double_t *c_re = 0;
   Double_t *c_im = 0;

// Create local arrays
   if (!(x_im = new (nothrow) Double_t[n])) {err = 1; goto cleanup;}
   if (!(y_im = new (nothrow) Double_t[n])) {err = 1; goto cleanup;}
   if (!(c_re = new (nothrow) Double_t[n])) {err = 1; goto cleanup;}
   if (!(c_im = new (nothrow) Double_t[n])) {err = 1; goto cleanup;}
   // important to initialize to zero
   for (Int_t i=0; i<n; i++) x_im[i] = y_im[i] = c_re[i] = c_im[i] = 0.0;

   TStat::FFT(nlog2, y_re, y_im);
   TStat::FFT(nlog2, x_re, x_im);
  
   for (Int_t i=0; i<n; i++){
      c_re[i] = y_re[i]*x_re[i] + y_im[i]*x_im[i];
      c_im[i] = y_im[i]*x_re[i] - y_re[i]*x_im[i];
   }//for_i
  
   TStat::IFFT(nlog2, c_re, c_im);

   for (Int_t i=0; i<n; i++){
      x_re[i] = c_re[i];
   }//for_i

cleanup:
   delete [] c_im;
   delete [] c_re;
   delete [] y_im;
   delete [] x_im;

   return err;
}//FFTDensityConvolve

//______________________________________________________________________________
void TStat::Kernelize(Int_t n, Double_t *x, Double_t bw, const char *kernel)
{
   // Kernelize arry x of length n with bandwidth bw
   // kernel is the type of kernel to use
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ and expanded for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::Kernelize------" << endl;
  
   Double_t a, ax;

   if (strcmp(kernel, "gaussian") == 0) {
   // Gaussian Kernel
      for (Int_t i=0; i<n; i++) {
         x[i] = TMath::Gaus(x[i], 0, bw, kTRUE);
      }//for_i
   } else if (strcmp(kernel, "epanechnikov") == 0) {
   // Epanechnikov Kernel
      a = bw * TMath::Sqrt(5.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = 3.0/(4.0*a)*(1.0 - (ax/a)*(ax/a));
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "rectangular") == 0) {
   // Rectangular Kernel
      a = bw * TMath::Sqrt(3.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = 0.5 / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "triangular") == 0) {
   // Triangular Kernel
      a = bw * TMath::Sqrt(6.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = (1.0 - ax/a) / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "biweight") == 0) {
   // Biweight Kernel
      a = bw * TMath::Sqrt(7.0);
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = 15.0/16.0 * (1.0 - (ax/a)*(ax/a)) * (1.0 - (ax/a)*(ax/a)) / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "cosine") == 0) {
   // Cosine Kernel
      a = bw / TMath::Sqrt(1.0/3.0 - 2.0/(TMath::Pi()*TMath::Pi()));
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = (1.0 + TMath::Cos(TMath::Pi()*x[i]/a)) / (2*a);
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   } else if (strcmp(kernel, "optcosine") == 0) {
   // Optcosine Kernel
      a = bw / TMath::Sqrt(1.0 - 8.0/(TMath::Pi()*TMath::Pi()));
      for (Int_t i=0; i<n; i++) {
         ax = TMath::Abs(x[i]);
         if (ax < a) {
            x[i] = (TMath::Pi()/4.0)*TMath::Cos(TMath::Pi()*x[i]/(2*a)) / a;
         } else {
            x[i] = 0.0;
         }//if
      }//for_i
   }//if
}//Kernelize

//______________________________________________________________________________
Double_t TStat::Bandwidth(Int_t n, Double_t *x, Double_t iqr)
{
   // Compute kernel bandwidth for array x of size n
   // iqr is the interquartile range of x
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::Bandwidth------" << endl;

   Double_t hi, lo;

   // standard deviation of array x  
//   hi = TMath::RMS(n, x); // (x-m)/n not (x-m)/(n-1) =>is not stdev!!
   hi = TMath::Sqrt(TStat::Var(n, x, TStat::Mean(n, x)));
  
   if (hi > iqr) lo = iqr/1.34;
   else          lo = hi;

   if (lo == 0) {
      if      (hi != 0)               lo = hi;
      else if (TMath::Abs(x[1]) != 0) lo = TMath::Abs(x[1]);
      else                            lo = 1.0;
   }//if

   return (0.9*lo*TMath::Power((Double_t)n, -0.2));
}//Bandwidth 

//______________________________________________________________________________
void TStat::LinearInterpolate(Double_t *xin, Double_t *yin,
                              Int_t nout, Double_t *xout, Double_t *yout)
{
   // Given xin and yin, interpolate linearly at xout and put results in yout
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::LinearInterpolate------" << endl;

   Int_t i, j, ij;

   for (Int_t k=0 ; k<nout; k++) {
      i = 0;
      j = nout - 1;
    
      if(xout[k] < xin[i]) {yout[k] = yin[0];      continue;}
      if(xout[k] > xin[j]) {yout[k] = yin[nout-1]; continue;}
 
      // find the correct interval by bisection
      while (i < j - 1) {  // xin[i] <= xout[k] <= xin[j]
         ij = (i + j)/2;   // i+1 <= ij <= j-1
         if (xout[k] < xin[ij]) j = ij;
         else                   i = ij; // still i < j
      }//while
  
      if (xout[k] == xin[j]) {yout[k] = yin[j]; continue;}
      if (xout[k] == xin[i]) {yout[k] = yin[i]; continue;}
  
      yout[k] = yin[i] + (yin[j] - yin[i])*((xout[k] - xin[i])/(xin[j] - xin[i]));
   }//for_k
}//LinearInterpolate

//______________________________________________________________________________
Int_t TStat::Density(Int_t n, Double_t *x, Double_t *w, Int_t nout,  Double_t *xout,
                     Double_t *yout, const char *kernel)
{
   // Kernel density of array x with weight w of size n
   // yout is the array of density values of size nout
   // xout is the array of the corresponding x coordinates
   // Note: size nout of output should be a power of two, preferably 512 or above
   // Important note:
   // Adapted from file "weightedkerneldensity.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::Density------" << endl;

   Int_t err   = 0;
   Int_t nout2 = 2*nout;
   Double_t lo, hi, iqr, bw, from, to;

// Init local arrays
   Int_t    *indx = 0;
   Double_t *sort = 0;
   Double_t *xin  = 0;
   Double_t *yin  = 0;
   Double_t *yord = 0;

// Create local arrays
   if (!(indx = new (nothrow) Int_t[n]))        {err = 1; goto cleanup;}
   if (!(sort = new (nothrow) Double_t[n]))     {err = 1; goto cleanup;}
   if (!(xin  = new (nothrow) Double_t[nout]))  {err = 1; goto cleanup;}
   if (!(yin  = new (nothrow) Double_t[nout2])) {err = 1; goto cleanup;}
   if (!(yord = new (nothrow) Double_t[nout2])) {err = 1; goto cleanup;}
   for (Int_t i=0; i<nout;  i++) xin[i] = 0.0;
   for (Int_t i=0; i<nout2; i++) yin[i] = yord[i] = 0.0;

// Sort x to get lo, hi and interquartile range
   TMath::Sort(n, x, indx, kFALSE);
   for (Int_t i=0; i<n; i++) sort[i] = x[indx[i]];
  
   lo  = sort[0];
   hi  = sort[n-1];
   iqr = sort[(Int_t)(0.75*n + 0.5)] - sort[(Int_t)(0.25*n + 0.5)];

// Bandwidth  
   bw = TStat::Bandwidth(n, x, iqr);
  
   lo = lo - 7.0*bw;
   hi = hi + 7.0*bw;

   for (Int_t i=0; i<=nout; i++) {
      yin[i] = (Double_t)i/(Double_t)(nout2 - 1)*2*(hi - lo);
   }//for_i

   for (Int_t i=nout+1; i<nout2; i++) {
      yin[i] = -yin[nout2 - i];
   }//for_i

   TStat::Kernelize(nout2, yin, bw, kernel);

   TStat::MassDist(n, x, w, lo, hi, nout, yord);

   err = TStat::FFTDensityConvolve(nout2, yin, yord);
   if (err != 0) goto cleanup;

   // corrections to get correct output range 
   to   = hi - 4.0*bw;
   from = lo + 4.0*bw;

   for (Int_t i=0; i<nout; i++) {
      xin[i]  = (Double_t)i / (Double_t)(nout - 1)*(hi - lo)   + lo;
      xout[i] = (Double_t)i / (Double_t)(nout - 1)*(to - from) + from;
   }//for_i

   for (Int_t i=0; i<nout; i++) {
      yin[i] = yin[i]/nout2;
   }//for_i

   TStat::LinearInterpolate(xin, yin, nout, xout, yout);

cleanup:
   delete [] xin;
   delete [] yord;
   delete [] yin;
   delete [] sort;
   delete [] indx;

   return err;
}//Density

//______________________________________________________________________________
Double_t TStat::MaxDensity(Int_t n, Double_t *x, Double_t *w, Int_t npts,
                           const char *kernel)
{
   // Compute density of array x with weight w and return maximal density
   // Note: number of points npts should be a power of two, preferably 512 or above
   // Important note:
   // Adapted from file "rma_background3.c" of program "RMAExpress", see:
   // http://stat-www.berkeley.edu/users/bolstad/RMAExpress/RMAExpress.html
   // Created by Ben Bolstad  <bolstad@stat.berkeley.edu>
   // License: GPL V2 or later
   // Converted to C++ for XPS by Christian Stratowa 
   if(kCSa) cout << "------TStat::MaxDensity------" << endl;
 
   Int_t   i     = 0;
   Double_t xmax = 0;
   Double_t ymax = 0;

// Init local arrays
   Double_t *arr  = 0;
   Double_t *xden = 0;
   Double_t *yden = 0;
   if (!(arr  = new (nothrow) Double_t[n]))    goto cleanup;
   if (!(xden = new (nothrow) Double_t[npts])) goto cleanup;
   if (!(yden = new (nothrow) Double_t[npts])) goto cleanup;

   for (i=0; i<n;    i++) arr[i]  = x[i];
   for (i=0; i<npts; i++) xden[i] = yden[i] = 0;
  
   TStat::Density(n, arr, w, npts, xden, yden, kernel);

//   for (i=0; i<npts; i++) ymax = (yden[i] > ymax) ? yden[i] : ymax;
   ymax = TStat::Max(npts, yden);

// Get index for max density of xden   
   i = 0;
   while(1) {
      if (yden[i] == ymax) break;
      i++;
   }//while
   xmax = xden[i];

// Cleanup
cleanup:
   if (yden) {delete [] yden; yden = 0;}
   if (xden) {delete [] xden; xden = 0;}
   if (arr)  {delete [] arr;  arr  = 0;}

  return xmax;
}//MaxDensity

//______________________________________________________________________________
Double_t TStat::PNormApprox(Double_t x)
{
   // Numerical approximation to the normal distribution as described in   
   // Abramowitz and Stegun: Handbook of Mathematical functions
   //                        page 931: 26.2.1, 932:26.2.17
   // Note: taken from simpleaffy2.c in the Bioconductor simpleaffy package
   // Copyright (C) 2004 Crispin Miller
   if(kCSa) cout << "------TStat::PNormApprox------" << endl;

   if (x >  6.0) return 1.0;
   if (x < -6.0) return 0.0;

   Double_t b1 =  0.31938153; 
   Double_t b2 = -0.356563782; 
   Double_t b3 =  1.781477937;
   Double_t b4 = -1.821255978;
   Double_t b5 =  1.330274429; 
   Double_t p  =  0.2316419; 
   Double_t c2 =  0.3989423;

   Double_t a = TMath::Abs(x); 
   Double_t t = 1.0 / (1.0 + a*p); 
   Double_t b = c2 * TMath::Exp((-x)*(x/2.0)); 
   Double_t n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t; 

   n = 1.0 - b*n;
   if (x < 0.0) n = 1.0 - n; 

   return n; 
}//PNormApprox


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMEstimator                                                          //
//                                                                      //
// Base class for M-estimators                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
Double_t TMEstimator::Calculate(Double_t x, Int_t type)
{
   if (type == TMEstimator::kRho) {
      return Rho(x);
   } else if (type == TMEstimator::kPsi) {
      return Psi(x);
   } else if (type == TMEstimator::kDeriv) {
      return Derivative(x);
   } else if (type == TMEstimator::kWeight) {
      return Weight(x);
   }//if

   return 0;
}//Calculate

//______________________________________________________________________________
TMEstimator *TMEstimator::Estimator(const char *name, Double_t c)
{
   TMEstimator *estimator = 0;
   if (strcmp(name, "huber") == 0) {
      if (c == 0.0) estimator = new THuberEstimator();
      else          estimator = new THuberEstimator(c);
   } else if (strcmp(name, "fair") == 0) {
      if (c == 0.0) estimator = new TFairEstimator();
      else          estimator = new TFairEstimator(c);
   } else if (strcmp(name, "cauchy") == 0) {
      if (c == 0.0) estimator = new TCauchyEstimator();
      else          estimator = new TCauchyEstimator(c);
   } else if (strcmp(name, "gemanmcclure") == 0) {
      if (c == 0.0) estimator = new TGemanMcClureEstimator();
      else          estimator = new TGemanMcClureEstimator(c);
   } else if (strcmp(name, "welsch") == 0) {
      if (c == 0.0) estimator = new TWelschEstimator();
      else          estimator = new TWelschEstimator(c);
   } else if (strcmp(name, "tukey") == 0) {
      if (c == 0.0) estimator = new TTukeyEstimator();
      else          estimator = new TTukeyEstimator(c);
   } else if (strcmp(name, "andrew") == 0) {
      if (c == 0.0) estimator = new TAndrewEstimator();
      else          estimator = new TAndrewEstimator(c);
   }//if

   return estimator;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THuberEstimator                                                      //
//                                                                      //
// Class for M-estimator Huber                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
THuberEstimator::THuberEstimator()
                :TMEstimator()
{
   // Default HuberEstimator constructor
   if(kCS) cout << "---THuberEstimator::THuberEstimator(default)------" << endl;

   fConst = 1.345;
   fName  = "Huber";
}//Constructor

//______________________________________________________________________________
THuberEstimator::THuberEstimator(Double_t c)
                :TMEstimator(c, "Huber")
{
   // Normal HuberEstimator constructor
   if(kCS) cout << "---THuberEstimator::THuberEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
THuberEstimator::~THuberEstimator()
{
   // HuberEstimator destructor
   if(kCS) cout << "---THuberEstimator::~THuberEstimator------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFairEstimator                                                       //
//                                                                      //
// Class for M-estimator Fair                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TFairEstimator::TFairEstimator()
               :TMEstimator()
{
   // Default FairEstimator constructor
   if(kCS) cout << "---TFairEstimator::TFairEstimator(default)------" << endl;

   fConst = 1.3998;
   fName  = "Fair";
}//Constructor

//______________________________________________________________________________
TFairEstimator::TFairEstimator(Double_t c)
               :TMEstimator(c, "Fair")
{
   // Normal FairEstimator constructor
   if(kCS) cout << "---TFairEstimator::TFairEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
TFairEstimator::~TFairEstimator()
{
   // FairEstimator destructor
   if(kCS) cout << "---TFairEstimator::~TFairEstimator------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCauchyEstimator                                                     //
//                                                                      //
// Class for M-estimator Cauchy                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TCauchyEstimator::TCauchyEstimator()
                 :TMEstimator()
{
   // Default CauchyEstimator constructor
   if(kCS) cout << "---TCauchyEstimator::TCauchyEstimator(default)------" << endl;

   fConst = 2.3849;
   fName  = "Cauchy";
}//Constructor

//______________________________________________________________________________
TCauchyEstimator::TCauchyEstimator(Double_t c)
                 :TMEstimator(c, "Cauchy")
{
   // Normal CauchyEstimator constructor
   if(kCS) cout << "---TCauchyEstimator::TCauchyEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
TCauchyEstimator::~TCauchyEstimator()
{
   // CauchyEstimator destructor
   if(kCS) cout << "---TCauchyEstimator::~TCauchyEstimator------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGemanMcClureEstimator                                               //
//                                                                      //
// Class for M-estimator GemanMcClure                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TGemanMcClureEstimator::TGemanMcClureEstimator()
                       :TMEstimator()
{
   // Default GemanMcClureEstimator constructor
   if(kCS) cout << "---TGemanMcClureEstimator::TGemanMcClureEstimator(default)------" << endl;

   fConst = 1.0;
   fName  = "GemanMcClure";
}//Constructor

//______________________________________________________________________________
TGemanMcClureEstimator::TGemanMcClureEstimator(Double_t c)
                       :TMEstimator(c, "GemanMcClure")
{
   // Normal GemanMcClureEstimator constructor
   if(kCS) cout << "---TGemanMcClureEstimator::TGemanMcClureEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
TGemanMcClureEstimator::~TGemanMcClureEstimator()
{
   // GemanMcClureEstimator destructor
   if(kCS) cout << "---TGemanMcClureEstimator::~TGemanMcClureEstimator------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TWelschEstimator                                                     //
//                                                                      //
// Class for M-estimator Welsch                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TWelschEstimator::TWelschEstimator()
                 :TMEstimator()
{
   // Default WelschEstimator constructor
   if(kCS) cout << "---TWelschEstimator::TWelschEstimator(default)------" << endl;

   fConst = 2.9846;
   fName  = "Welsch";
}//Constructor

//______________________________________________________________________________
TWelschEstimator::TWelschEstimator(Double_t c)
                 :TMEstimator(c, "Welsch")
{
   // Normal WelschEstimator constructor
   if(kCS) cout << "---TWelschEstimator::TWelschEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
TWelschEstimator::~TWelschEstimator()
{
   // WelschEstimator destructor
   if(kCS) cout << "---TWelschEstimator::~TWelschEstimator------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTukeyEstimator                                                      //
//                                                                      //
// Class for M-estimator Tukey                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TTukeyEstimator::TTukeyEstimator()
                :TMEstimator()
{
   // Default TukeyEstimator constructor
   if(kCS) cout << "---TTukeyEstimator::TTukeyEstimator(default)------" << endl;

   fConst = 4.6851;
   fName  = "Tukey";
}//Constructor

//______________________________________________________________________________
TTukeyEstimator::TTukeyEstimator(Double_t c)
                :TMEstimator(c, "Tukey")
{
   // Normal TukeyEstimator constructor
   if(kCS) cout << "---TTukeyEstimator::TTukeyEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
TTukeyEstimator::~TTukeyEstimator()
{
   // TukeyEstimator destructor
   if(kCS) cout << "---TTukeyEstimator::~TTukeyEstimator------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TAndrewEstimator                                                     //
//                                                                      //
// Class for M-estimator Andrew                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
TAndrewEstimator::TAndrewEstimator()
                 :TMEstimator()
{
   // Default AndrewEstimator constructor
   if(kCS) cout << "---TAndrewEstimator::TAndrewEstimator(default)------" << endl;

   fConst = 1.339;
   fName  = "Andrew";
}//Constructor

//______________________________________________________________________________
TAndrewEstimator::TAndrewEstimator(Double_t c)
                 :TMEstimator(c, "Andrew")
{
   // Normal AndrewEstimator constructor
   if(kCS) cout << "---TAndrewEstimator::TAndrewEstimator------" << endl;

}//Constructor

//______________________________________________________________________________
TAndrewEstimator::~TAndrewEstimator()
{
   // AndrewEstimator destructor
   if(kCS) cout << "---TAndrewEstimator::~TAndrewEstimator------" << endl;

}//Destructor

