// File created: 12/16/2002                          last modified: 07/03/2010
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2010 Dr. Christian Stratowa
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

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <cstdlib>

#include "XPSFilter.h"

#include "XPSSchemes.h"
#include "XPSData.h"
#include "XPSUtils.h"

#include "TBranch.h"
#include "TLeaf.h"
#include "TTree.h"
//#include "TList.h"

const Int_t kLineBuf = 16635;   //to read lines

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XFilter);
ClassImp(XPreFilter);
ClassImp(XUniFilter);
ClassImp(XMultiFilter);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFilter                                                              //
//                                                                      //
// Base class for filter                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XFilter::XFilter()
        :XAlgorithm()
{
   // Default Filter constructor
   if(kCS) cout << "---XFilter::XFilter(default)------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XFilter::XFilter(const char *name, const char *type)
        :XAlgorithm(name, type)
{
   // Normal Filter constructor
   if(kCS) cout << "---XFilter::XFilter------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XFilter::XFilter(const char *name, const char *type, Double_t na)
        :XAlgorithm(name, type)
{
   // Normal Filter constructor
   if(kCS) cout << "---XFilter::XFilter------" << endl;

   this->Init();
 
   fNA    = na;
   fHasNA = kTRUE;
}//Constructor

//______________________________________________________________________________
XFilter::~XFilter()
{
   // Filter destructor
   if(kCS) cout << "---XFilter::~XFilter------" << endl;

   if (fMask) {delete [] fMask; fMask = 0;}
}//Destructor

//______________________________________________________________________________
void XFilter::Init()
{
   // Initialize parameters
   if(kCS) cout << "------XFilter::Init------" << endl;

   fNData      = 0;
   fSorted     = 0; 
   fMean       = 0;
   fTrim       = 0.0;
   fVar        = 0;
   fMin        = 0;
   fMax        = 0;
   fEpsilon    = 0.001;
   fMinFilters = 1;
   fNMask      = 0;
   fMask       = 0;
}//Init

//______________________________________________________________________________
Int_t XFilter::Initialize(Int_t min, Bool_t /*reset*/)
{
   // Re-initialize filters:
   // - min: minimum number of selected filter methods to satisfy
   if(kCS) cout << "------XFilter::Initialize------" << endl;

   fMinFilters = min;

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XFilter::MeanVarMinMax(Int_t n, Double_t *arr)
{
   // Calculate mean, variance, minimum, maximum for array arr of data
   // Note: this is a helper-function to sort array arr and calculate values  
   //       used in different filter methods, only once per row
   if(kCSa) cout << "------XFilter::MeanVarMinMax------" << endl;

//TO DO!!! if(fHasNA)

// Init array sorted only once
   if (!fSorted) if (!(fSorted = new Double_t[n])) return errInitMemory;

   if (n == 1) {
      fVar = 0;
      fMin = fMax = fMean = fSorted[0] = arr[0];
      return errNoErr;
   }//if

// Create index and sort array
   Int_t *index = 0;
   if (!(index = new Int_t[n])) return errInitMemory;
   TMath::Sort(n, arr, index, kFALSE);

   // start-index and end-index
   Int_t start, end;
   if (fTrim < 0.5) {
      start = (Int_t)floor(n * fTrim);
      end   = n - start;
   } else {
      if ((n % 2) == 0){
         start = (Int_t)floor(n / 2.0) - 1;
         end   = start + 2;
      } else {
         start = (Int_t)floor(n / 2.0);
         end   = start + 1;
     }//if
   }//if

// Calculate trimmed mean
   Int_t    trimlen = end - start;
   Double_t mean    = 0;
   for (Int_t i=0; i<n; i++) {
      fSorted[i] = arr[index[i]];
      mean += (i>=start && i<end) ? arr[index[i]] : 0;
   }//for_i
   mean /= trimlen;

// Calculate variance
   Double_t var = 0;
   Double_t tmp = 0;
   for (Int_t i=start; i<end; i++) {
      tmp = arr[index[i]] - mean;
      var += tmp * tmp;
   }//for_i
   if (trimlen > 1) var /= (trimlen - 1);
   else var = 0;

   delete [] index;

/////////////
//PROBLEM: need to test fEspilon=0!!!
/////////////
// add fEpsilon to prevent zero division or to set fMean to 1 (for fEpsilon=0)
   fMean = (fEpsilon > 0) ? ((mean != 0) ? mean : fEpsilon) : 1;
   fVar  = var;
   fMin  = (fSorted[0] != 0) ? fSorted[0] : ((fEpsilon > 0) ? fEpsilon : 1);
   fMax  = fSorted[n-1];

   return errNoErr;
}//MeanVarMinMax

//______________________________________________________________________________
Int_t XFilter::FillMaskTree(TTree *unittree, TTree *masktree, Int_t n, Int_t *arr)
{
   // Fill masktree with unit branch from unittree and with mask array arr
   if(kCS) cout << "------XFilter::FillMaskTree------" << endl;

   if ((unittree == 0) || (masktree == 0)) return errGetTree;

// Check if trees have equal entries
   Int_t nentries = (Int_t)(unittree->GetEntries());
   if (nentries != n) {
      cerr << "Error: Tree <" << unittree->GetName() << " has not <"
           << n << "> entries." << endl;
      return errAbort;
   }//if

// Get leaf and branch from unittree
   TLeaf   *leaf = unittree->FindLeaf("fUnitID");
   TBranch *brch = leaf->GetBranch();

// Create leaf and branch for masktree
   Int_t split   = 99;
   Int_t buffer  = 64000; //other size????
   XUnitID *unit = new XUnitID();
   XMask   *mask = new XMask();

   masktree->Branch("UnitBranch", "XUnitID", &unit, buffer, split);
   masktree->Branch("MaskBranch", "XMask", &mask, buffer, split);

// Fill masktree
   for (Int_t i=0; i<n; i++) {
      // read entry i from unittree
      brch->GetEntry(i);
      Int_t id = (Int_t)leaf->GetValue();

      unit->SetUnitID(id);
      mask->SetFlag(arr[i]);

      masktree->Fill();
   }//for_i

   return errNoErr;
}//FillMaskTree


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreFilter                                                           //
//                                                                      //
// Class for nonspecific filter                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPreFilter::XPreFilter()
           :XFilter()
{
   // Default PreFilter constructor
   if(kCS) cout << "---XPreFilter::XPreFilter(default)------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XPreFilter::XPreFilter(const char *name, const char *type)
           :XFilter(name, type)
{
   // Normal PreFilter constructor
   if(kCS) cout << "---XPreFilter::XPreFilter------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XPreFilter::XPreFilter(const char *name, const char *type, Double_t na)
           :XFilter(name, type, na)
{
   // Normal Filter constructor
   if(kCS) cout << "---XPreFilter::XPreFilter------" << endl;

   this->Init();
}//Constructor

//______________________________________________________________________________
XPreFilter::~XPreFilter()
{
   // Filter destructor
   if(kCS) cout << "---XPreFilter::~XPreFilter------" << endl;

}//Destructor

//______________________________________________________________________________
void XPreFilter::Init()
{
   // Initialize parameters
   if(kCS) cout << "------XPreFilter::Init------" << endl;

   fMAD           = 0.0;
   fCov2mn        = 0.0;
   fVar2mn        = 0.0;
   fDif2mn        = 0.0;
   fMax2min       = 0.0;
   fGap2mn        = 0.0;
   fWindow        = 0.0;
   fLoThreshold   = 0.0;
   fLoCondition   = "percent";
   fLowerID       = 1;
   fLoSamples     = 100.0;
   fUpThreshold   = 0.0;
   fUpCondition   = "percent";
   fUpperID       = 1;
   fUpSamples     = 100.0;
   fQRatio        = 1.0;
   fLoQ           = 0.05;
   fHiQ           = 0.95;
   fEntropy       = 1.0;
   fNQuantiles    = 10;
   fCallCondition = "percent";
//??   fCallID        = 1;
   fCallSamples   = 100.0;
   fCallPValue    = 1.0;
   fNCall         = 0;

   fHasMAD = fHasCov = fHasVar = fHasDif = fHasM2m = fHasGap = 0;
   fHasLoT = fHasUpT = fHasQua = fHasEnt = fHasCal = kFALSE;
}//Init

//______________________________________________________________________________
Int_t XPreFilter::Initialize(Int_t min, Bool_t reset)
{
   // Re-initialize filters:
   // - min: minimum number of selected filter methods to satisfy
   //        min = 1:  equivalent to "OR", i.e. satisfy at least one filter
   //        min = 10: equivalent to "AND", i.e. satisfy all 10 selected filters
   //        e.g. 5 filters selected and min=3: 3 of 5 filters must be satisfied
   // - reset:  clear filters, i.e. reset filter methods to kFALSE
   //   Note:  after resetting, filter methods must be initialized explicitely
   if(kCS) cout << "------XPreFilter::Initialize------" << endl;

   fMinFilters = min;

   // reset all filters
   if (reset) {
      fHasMAD = fHasCov = fHasVar = fHasDif = fHasM2m = fHasGap = 0;
      fHasLoT = fHasUpT = fHasQua = fHasEnt = fHasCal = kFALSE;
   }//if

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XPreFilter::InitType(const char *type, Option_t *options,
                  Int_t npars, Double_t *pars)
{
   // Initialize filter type
   if(kCS) cout << "------XPreFilter::InitType------" << endl;

   if (strcmp(type, "variation") == 0) {
      return this->InitVariation(options, npars, pars);
   } else if (strcmp(type, "lowerthreshold") == 0) {
      return this->InitLowerThreshold(options, npars, pars);
   } else if (strcmp(type, "upperthreshold") == 0) {
      return this->InitUpperThreshold(options, npars, pars);
   } else if (strcmp(type, "quantile") == 0) {
      return this->InitQuantile(npars, pars);
   } else if (strcmp(type, "entropy") == 0) {
      return this->InitEntropy(npars, pars);
   } else if (strcmp(type, "call") == 0) {
      return this->InitCall(options, npars, pars);
   } else {
      cerr << "Error: PreFilter algorithm <" << type << "> not known" << endl;
   }//if

   return errInitSetting;
}//InitType

//______________________________________________________________________________
Int_t XPreFilter::InitVariation(const char *varlist, Int_t npars, Double_t *pars)
{
   // Initialize variation filters defined in varlist:
   //    "mad"     - median absolute deviation
   //    "cov2mn"  - coefficient of variation sqrt(var)/mean
   //    "var2mn"  - variance variance/mean
   //    "dif2mn"  - difference (max - min)/mean
   //    "max2min" - ratio max/min
   //    "gap2mn"  - gap between levels gap/mean
   // e.g. varlist = "var2mn:gap2mn"
   // npars parameters must be given in order of varlist, i.e. first the
   // cutoff levels for each filter, and at last followed by;
   // - window size for gap if gap2mn
   // - trim value for mean if cov2mn or var2mn or dif2mn or gap2mn
   // - epsilon as value to replace mean and min:
   //      epsilon > 0: replace mean=0 and min=0 (to prevent zero division)
   //      epsilon = 0: set always mean=1 (and for min=0 set min=1)
   // e.g. for "var2mn:gap2mn": 5, 0.35, 0.4, 0.0(trim), 0.05(window), 0.01(epsilon) 
   if(kCS) cout << "------XPreFilter::InitVariation------" << endl;

   fHasMAD = fHasCov = fHasVar = fHasDif = fHasM2m = fHasGap = 0;

   Int_t idx = 0;
   char *vname = new char[strlen(varlist) + 1];
   char *dname = vname;
   vname = strtok(strcpy(vname,varlist),":");
   while(vname) {
      if (strcmp(vname,"mad")     == 0) {fHasMAD = ++idx;}
      if (strcmp(vname,"cov2mn")  == 0) {fHasCov = ++idx;}
      if (strcmp(vname,"var2mn")  == 0) {fHasVar = ++idx;}
      if (strcmp(vname,"dif2mn")  == 0) {fHasDif = ++idx;}
      if (strcmp(vname,"max2min") == 0) {fHasM2m = ++idx;}
      if (strcmp(vname,"gap2mn")  == 0) {fHasGap = ++idx;}
      vname = strtok(NULL, ":");
      if (vname == 0) break;
   }//while
   delete [] dname;

   Short_t n = idx;
   if (fHasGap) idx++;  //for window size
   if (fHasCov || fHasVar || fHasDif || fHasGap) idx++;  //for trim value
   idx++; //for fEpsilon

   if (idx != npars) {
      cerr << "Error: Number of parameters for varlist <"      << idx  
           << "> is not equal to given number of parameters <" << npars << ">."
           << endl;
      return errAbort;
   }//if

   idx = 0;
   for (Short_t j=1; j<=n; j++) {
      if (fHasMAD == j) fMAD     = pars[idx++];
      if (fHasCov == j) fCov2mn  = pars[idx++];
      if (fHasVar == j) fVar2mn  = pars[idx++];
      if (fHasDif == j) fDif2mn  = pars[idx++];
      if (fHasM2m == j) fMax2min = pars[idx++];
      if (fHasGap == j) fGap2mn  = pars[idx++];
   }//for_j
   if (fHasGap) fWindow = pars[idx++];
   if (fHasCov || fHasVar || fHasDif || fHasGap) fTrim = pars[idx++];

   fEpsilon = (pars[idx] > 0.0) ? pars[idx] : 0.0;
//??   fEpsilon = (pars[idx] > 0.0) ? pars[idx] : fEpsilon;

   return errNoErr;
}//InitVariation

//______________________________________________________________________________
Int_t XPreFilter::InitLowerThreshold(const char *condition, Int_t npars,
                  Double_t *pars)
{
   // Initialize lower threshold filter
   // pars[0] is lower threshold
   // condition = "percent":    pars[1] is percent of samples
   // condition = "samples":    pars[1] is number of samples
   // condition = "mean":       pars[1] is mean of samples
   // condition = "percentile": pars[1] is percentile of samples
   if(kCS) cout << "------XPreFilter::InitLowerThreshold------" << endl;

   if (npars != 2) return errInitSetting;

   fLoCondition = condition;
   fLoThreshold = pars[0];
   fLoSamples   = pars[1];

   fHasLoT = kTRUE;
   return errNoErr;
}//InitLowerThreshold

//______________________________________________________________________________
Int_t XPreFilter::InitUpperThreshold(const char *condition, Int_t npars,
                  Double_t *pars)
{
   // Initialize upper threshold filter
   // pars[0] is upper threshold
   // condition = "percent":    pars[1] is percent of samples
   // condition = "samples":    pars[1] is number of samples
   // condition = "mean":       pars[1] is mean of samples
   // condition = "percentile": pars[1] is percentile of samples
   if(kCS) cout << "------XPreFilter::InitUpperThreshold------" << endl;

   if (npars != 2) return errInitSetting;

   fUpCondition = condition;
   fUpThreshold = pars[0];
   fUpSamples   = pars[1];

   fHasUpT = kTRUE;
   return errNoErr;
}//InitUpperThreshold

//______________________________________________________________________________
Int_t XPreFilter::InitQuantile(Int_t npars, Double_t *pars)
{
   // Initialize quantile filter
   // pars[0] is ratio of high quantile to low quantile
   // pars[1] is low quantile
   // pars[2] is high quantile
   if(kCS) cout << "------XPreFilter::InitQuantile------" << endl;

   if (npars != 3) return errInitSetting;

   fQRatio = pars[0];
   fLoQ    = pars[1];
   fHiQ    = pars[2];

   fHasQua = kTRUE;
   return errNoErr;
}//InitQuantile

//______________________________________________________________________________
Int_t XPreFilter::InitEntropy(Int_t npars, Double_t *pars)
{
   // Initialize entropy filter
   // pars[0] is cutoff entropy
   // pars[1] is number of quantiles
   if(kCS) cout << "------XPreFilter::InitEntropy------" << endl;

   if (npars != 2) return errInitSetting;

   fEntropy    = pars[0];
   fNQuantiles = (Int_t)pars[1];

   fHasEnt = kTRUE;
   return errNoErr;
}//InitEntropy

//______________________________________________________________________________
Int_t XPreFilter::InitCall(const char *condition , Int_t npars, Double_t *pars)
{
   // Initialize present call filter
   // pars[0]: range for detection p-value is [0,1]
   //          set pars[0]=1 to use "call" otherwise use "detection p-value"
   // condition = "percent": pars[1] is percent of samples less than pvalue 
   // condition = "samples": pars[1] is minimum number of samples less than pvalue 
   if(kCS) cout << "------XPreFilter::InitCall------" << endl;

   if (npars != 2) return errInitSetting;

   fCallCondition = condition;
   fCallPValue    = pars[0];
   fCallSamples   = pars[1];

   fHasCal = kTRUE;
   return errNoErr;
}//InitCall

//______________________________________________________________________________
void XPreFilter::InitThresholdConditions()
{
   // Initialize threshold conditions
   if(kCS) cout << "------XPreFilter::InitThresholdConditions------" << endl;

// Lower threshold
   Double_t numsamples = fLoSamples;
   if (strcmp(fLoCondition.Data(),"percent") == 0) {
      fLowerID = 1;
      numsamples = ceil(fNData*fLoSamples/100.0);
   } else if (strcmp(fLoCondition.Data(),"samples") == 0) {
      fLowerID = 2;
   } else if (strcmp(fLoCondition.Data(),"mean") == 0) {
      fLowerID = 3;
   } else if (strcmp(fLoCondition.Data(),"percentile") == 0) {
      fLowerID = 4;
   }//if
   fLoSamples = (numsamples < fNData) ? numsamples : fNData;

// Upper threshold
   numsamples = fUpSamples;
   if (strcmp(fUpCondition.Data(),"percent") == 0) {
      fUpperID = 1;
      numsamples = ceil(fNData*fUpSamples/100.0);
   } else if (strcmp(fUpCondition.Data(),"samples") == 0) {
      fUpperID = 2;
   } else if (strcmp(fUpCondition.Data(),"mean") == 0) {
      fUpperID = 3;
   } else if (strcmp(fUpCondition.Data(),"percentile") == 0) {
      fUpperID = 4;
   }//if
   fUpSamples = (numsamples < fNData) ? numsamples : fNData;

   return;
}//InitThresholdConditions

//______________________________________________________________________________
void XPreFilter::InitCallConditions()
{
   // Initialize call conditions
   if(kCS) cout << "------XPreFilter::InitCallConditions------" << endl;

// Present call
   Double_t numsamples = fCallSamples;
   if (strcmp(fCallCondition.Data(),"percent") == 0) {
      numsamples = ceil(fNCall*fCallSamples/100.0);
   }//if
   fCallSamples = (numsamples < fNCall) ? numsamples : fNCall;

   return;
}//InitCallConditions

//______________________________________________________________________________
Int_t XPreFilter::Calculate(const char *infile, const char *outfile, 
                  const char *varlist, Int_t nrows, const char *sepi,
                  const char *sepo, char delim)
{
   // Calculate filter
   // varlist is varlist of infile, e.g. "expr:call:pval"
   // outfile exports only those rows of infile which pass the filter conditions
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   // Note: since outfile has same columns as infile, sepo is equal to sepi
   if(kCS) cout << "------XPreFilter::Calculate------" << endl;

//////////////
// PROBLEM: also for log,log2,log10!!???
// maybe infile must have log data?
// BETTER: create new method Transform() in class X??? to create log-table!
//////////////

   Int_t err  = errNoErr;
   Int_t idx  = 0;
   Int_t okay = 0;
   Int_t nvar = 0;
   Int_t nhas = 0;

// Decompose varlist
   Bool_t hasExpr = kFALSE;
   Bool_t hasCall = kFALSE;
   Bool_t hasPVal = kFALSE;

   char *vname = new char[strlen(varlist) + 1];
   char *dname = vname;
   vname = strtok(strcpy(vname,varlist),":");
   while(vname) {
      if (strcmp(vname,"expr") == 0) {hasExpr = kTRUE; nhas++;}
      if (strcmp(vname,"call") == 0) {hasCall = kTRUE; nhas++;}
      if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE; nhas++;}
      nvar++;
      vname = strtok(NULL, ":");
      if (vname == 0) break;
      nvar++;
   }//while
   delete [] dname;

   if (nhas != nvar) {
      cerr << "Error: Only <expr>, <call>, <pval> are allowed in varlist" << endl;
      return errAbort;
   }//if

   if (!hasExpr) {
      cerr << "Error: varlist must have at least <expr> as variable" << endl;
      cerr << "       and <expr> must be the first variable in varlist." << endl;
      return errAbort;
   }//if

   if (!(hasCall || hasPVal)) fHasCal = kFALSE;

// Ensure that fMinFilters is not larger than number of selected filters
   fMinFilters = this->SetMinFilters(fMinFilters);

   char *ch;
   char  name[kBufSize];
   char  nextline[kLineBuf];

// Open input file and output file
   ifstream input(infile, ios::in);
   if (!input.good()) {
      cerr << "Error: Could not create input <" << infile << ">" << endl;
      return errReadingInput;
   }//if

   ofstream output(outfile, ios::out);
   if (!output) {
      cerr << "Error: Could not create output <" << outfile << ">" << endl;
      return errOpenOutput;
   }//if

// Init local arrays
   TString  *arrName = 0;  //row names
   Double_t *arrExpr = 0;  //expression levels
   Int_t    *arrCall = 0;  //present call
   Double_t *arrPVal = 0;  //detection p-value

// Read header to get number of data columns
   input.getline(nextline, kLineBuf, delim);
   if (!input.good()) {err = errReadingInput; goto cleanup;}

   fNData = (Int_t)(NumSeparators(nextline, sepi)/nvar); //numdata = numcolumns/numvariables
   fNCall = fHasCal ? fNData : 0;

// Create array fMask
   fNMask = nrows;
   if (!(fMask = new Int_t[nrows])) {err = errInitMemory; goto cleanup;}
   // init values to NA (= -1)
   for (Int_t i=0; i<nrows; i++) fMask[i] = -1;

// Create local arrays
   if (!(arrName = new TString[nrows]))   {err = errInitMemory; goto cleanup;}
   if (!(arrExpr = new Double_t[fNData])) {err = errInitMemory; goto cleanup;}
   if (!(arrCall = new Int_t[fNData]))    {err = errInitMemory; goto cleanup;}
   if (!(arrPVal = new Double_t[fNData])) {err = errInitMemory; goto cleanup;}

   // init values to NA (= -1)
   for (Int_t i=0; i<fNData; i++) {
      arrExpr[i] = -1;
      arrCall[i] = -1;
      arrPVal[i] = -1;
   }//for_i

// Initialize filter conditions
   this->InitThresholdConditions();
   this->InitCallConditions();

// Read data and compute mask
   while (input.good()) {
      input.getline(nextline, kLineBuf, delim);
      if (input.fail() || (idx == nrows)) break;

      // read row name
      strcpy(name,strtok(&nextline[0], sepi));
      arrName[idx] = name;

      // read row data
      for (Int_t i=0; i<fNData; i++) {
         arrExpr[i] = atof(strtok(NULL, sepi));

         if (hasCall) {
            ch = strtok(NULL, sepi);
            if      (ch[0] == 'P') arrCall[i] = 2;
            else if (ch[0] == 'A') arrCall[i] = 0;
            else if (ch[0] == 'M') arrCall[i] = 1;
            else                   arrCall[i] = -1;
         }//if

         if (hasPVal) {arrPVal[i] = atof(strtok(NULL, sepi));}
      }//for_i

      if ((err = MeanVarMinMax(fNData, arrExpr)) != errNoErr) goto cleanup;

      Short_t flag = 0;
      if (fHasMAD) flag += this->MAD();
      if (fHasCov) flag += this->CoefficientOfVariation();
      if (fHasVar) flag += this->Variance2Mean();
      if (fHasDif) flag += this->Difference2Mean();
      if (fHasM2m) flag += this->RatioMax2Min();
      if (fHasGap) flag += this->Gap2Mean();
      if (fHasLoT) flag += this->LowerThreshold();
      if (fHasUpT) flag += this->UpperThreshold();
      if (fHasQua) flag += this->QuantileHi2Lo();
      if (fHasCal) flag += this->PresentCall(arrCall, arrPVal);

      fMask[idx] = flag;
      idx++;
   }//while

   if (idx != nrows) {
      cout << "Warning: Number of lines read <"    << idx  
           << "> is not equal to number of rows <" << nrows << ">"
           << endl;
   }//if
   nrows = idx;

// Entropy filter
   if (fHasEnt) {
//      err = this->Entropy(nrows, fNData, table, arrMask);
   }//if

// Re-read infile and write only lines passing filter(s) to outfile
   input.clear();
   input.close();
   input.open(infile, ios::in);

   // read header
   input.getline(nextline, kLineBuf, delim);
   output << nextline << delim;

   // read data and write filtered data
   idx = 0;
   while (input.good()) {
      input.getline(nextline, kLineBuf, delim);
      if (input.fail() || (idx == nrows)) break;

      if (fMask[idx] >= fMinFilters) {output << nextline << delim; okay++;}
      idx++;
   }//while

   if (XManager::fgVerbose) {
      cout << "Prefilter: <" << okay << "> genes of <" << idx 
           << "> genes fulfill filter criteria." << endl;
   }//if

// Cleanup
cleanup:
   if (arrPVal) {delete [] arrPVal; arrPVal = 0;}
   if (arrCall) {delete [] arrCall; arrCall = 0;}
   if (arrExpr) {delete [] arrExpr; arrExpr = 0;}
   if (arrName) {delete [] arrName; arrName = 0;}

// Close files
   input.close();
   output.close();

   return err;
}//Calculate

//______________________________________________________________________________
Int_t XPreFilter::Calculate(Int_t n, TTree **intree, const char *leafname,
                  TTree *outtree, const char *varlist)
{
   // Calculate prefilter for n intrees with leaf leafname and save in outtree
   // n        - number of input trees
   // intree   - array of n input trees
   // leafname - name of intree leaf containing data to be used
   // outtree  - name of output filter tree 
   // varlist  - filters to be used as defined in kPreFltr, e.g. "cv:ratio:call"
   //            default is "*": all initialized filters will be used 
   if(kCS) cout << "------XPreFilter::Calculate(tree)------" << endl;

   Int_t err  = errNoErr;
   Int_t flag = 0;
   Int_t okay = 0;

   if (intree == 0 || outtree == 0) {
      cerr << "Error: Intree and/or outtree is missing." << endl;
      return errGetTree;
   }//if

// Decompose varlist (only if not "*")
   if (strcmp(varlist, "*") != 0) {
      Short_t hasMAD, hasCov, hasVar, hasDif, hasM2m, hasGap;
      Bool_t  hasLoT, hasHiT, hasQua, hasEnt, hasCal;
      hasMAD = hasCov = hasVar = hasDif = hasM2m = hasGap = 0;
      hasLoT = hasHiT = hasQua = hasEnt = hasCal = kFALSE;

      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;
      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname, kPreFltr[0])  == 0) {hasMAD = 1;}
         if (strcmp(vname, kPreFltr[1])  == 0) {hasCov = 1;}
         if (strcmp(vname, kPreFltr[2])  == 0) {hasVar = 1;}
         if (strcmp(vname, kPreFltr[3])  == 0) {hasDif = 1;}
         if (strcmp(vname, kPreFltr[4])  == 0) {hasM2m = 1;}
         if (strcmp(vname, kPreFltr[5])  == 0) {hasGap = 1;}
         if (strcmp(vname, kPreFltr[6])  == 0) {hasLoT = kTRUE;}
         if (strcmp(vname, kPreFltr[7])  == 0) {hasHiT = kTRUE;}
         if (strcmp(vname, kPreFltr[8])  == 0) {hasQua = kTRUE;}
         if (strcmp(vname, kPreFltr[9])  == 0) {hasEnt = kTRUE;}
         if (strcmp(vname, kPreFltr[10]) == 0) {hasCal = kTRUE;}
         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while
      delete [] dname;

      // if varlist defines filters to use then set fHasXXX (else use InitXXX())
      if (hasMAD || hasCov || hasVar || hasDif || hasM2m || hasGap || 
          hasLoT || hasHiT || hasQua || hasEnt || hasCal) {
         fHasMAD = hasMAD;
         fHasCov = hasCov;
         fHasVar = hasVar;
         fHasDif = hasDif;
         fHasM2m = hasM2m;
         fHasGap = hasGap;
         fHasLoT = hasLoT;
         fHasUpT = hasHiT;
         fHasQua = hasQua;
         fHasEnt = hasEnt;
         fHasCal = hasCal;
      }//if
   }//if

// Check leaves of intree
   if (intree[0]->FindLeaf(leafname)  == 0) {
      cerr << "Error: Tree does not have leaf <" << leafname << ">." << endl;
      return errAbort;
   }//if

   if (intree[0]->FindLeaf("fUnitID") == 0) {
      cerr << "Error: Tree does not have leaf <fUnitID>." << endl;
      return errAbort;
   }//if

// Ensure that fMinFilters is not larger than number of selected filters
   fMinFilters = this->SetMinFilters(fMinFilters);

// Get number of samples and number of entries
   Int_t nentries = (Int_t)(intree[0]->GetEntries());
   fNData = n;

   // test for call trees if InitCall() is initialized
   if (fHasCal && (fNCall == 0)) {
      cout << "Error: No call trees have been selected." << endl;
      return errAbort;
   }//if

   // usually number of call trees should be equal to number of expr trees
   if (fHasCal && (fNCall != fNData)) {
      cout << "Warning: Number of call trees used is not equal to number of expression trees."
           << endl;
   }//if

//   TLeaf   *leafUnit = intree[0]->FindLeaf("fUnitID");
//   TBranch *brchUnit = leafUnit->GetBranch();

   TBranch **brchj = new TBranch*[fNData];
   TLeaf   **leafj = new TLeaf*[fNData];
   for (Int_t j=0; j<fNData; j++) {
      leafj[j] = intree[j]->FindLeaf(leafname);      
      brchj[j] = leafj[j]->GetBranch();
   }//for_j

// Create expr array
   Double_t *arrExpr = 0;  //expression levels
   if (!(arrExpr = new Double_t[fNData])) {err = errInitMemory; goto cleanup;}
   // init values to NA (= -1)
   for (Int_t i=0; i<fNData; i++) arrExpr[i] = -1;

// Create array fMask only if it does not exist already:
// if CallFlag() was called before then fMask contains the results
   if (!fMask) {
      fNMask = nentries;
      if (!(fMask = new Int_t[nentries])) {err = errInitMemory; goto cleanup;}
      for (Int_t i=0; i<nentries; i++) fMask[i] = 0;
   }//if

// Initialize filter conditions
   this->InitThresholdConditions();

// Read intree entries, calculate mask and fill outtree
   for (Int_t i=0; i<nentries; i++) {
//      brchUnit->GetEntry(i);
//      Int_t id = (Int_t)leafUnit->GetValue();

      // read entry i from tree plus friends
      for (Int_t j=0; j<fNData; j++) {
//no         brchj[j]->GetEntry(id);
         brchj[j]->GetEntry(i);
         arrExpr[j] = leafj[j]->GetValue();
      }//for_j

      if ((err = MeanVarMinMax(fNData, arrExpr)) != errNoErr) goto cleanup;

      flag = 0;
      if (fHasMAD) flag += this->MAD();
      if (fHasCov) flag += this->CoefficientOfVariation();
      if (fHasVar) flag += this->Variance2Mean();
      if (fHasDif) flag += this->Difference2Mean();
      if (fHasM2m) flag += this->RatioMax2Min();
      if (fHasGap) flag += this->Gap2Mean();
      if (fHasLoT) flag += this->LowerThreshold();
      if (fHasUpT) flag += this->UpperThreshold();
      if (fHasQua) flag += this->QuantileHi2Lo();
      if (fHasCal) flag += GetMask(i);  //add mask from CallFlag()

//      flag = (flag >= fMinFilters);
//BETTER?? since Entropy must be calculated later:
      fMask[i] = flag;

//ev.BETTER?? separate arrMsk for each type of filter:
//store number of selected(remaining) rows for each filter in fNxxxFlags
//in cout print result for each filter
//in XPSApp display result for each filter in GUI
   }//for_i

// Entropy filter
   if (fHasEnt) {
//      err = this->Entropy(nentries, fNData, table, fMask);
//better?      err = this->Entropy();
   }//if

   for (Int_t i=0; i<nentries; i++) {
      fMask[i] = (Int_t)(fMask[i] >= fMinFilters);
      if (fMask[i]) okay++;
   }//for_i

// Fill outtree with mask
   err = FillMaskTree(intree[0], outtree, nentries, fMask);

   if (XManager::fgVerbose) {
      cout << "Prefilter: <" << okay << "> genes of <" << nentries 
           << "> genes fulfill filter criteria." << endl;
   }//if

// Cleanup
cleanup:
   if (arrExpr) {delete [] arrExpr; arrExpr = 0;}

   delete [] leafj;
   delete [] brchj;

   return err;
}//Calculate

//______________________________________________________________________________
Int_t XPreFilter::CallFlag(Int_t n, TTree **intree, const char *varlist,
                  TTree *outtree)
{
   // Calculate call flag for n intrees with leaf leafname and save in outtree
   // n        - number of input trees
   // intree   - array of n input trees
   // varlist  - filters to be used as defined in kPreFltr, i.e. "call", "pval"
   //            default is "*": if intree has leaf fPValue then "pval" will be
   //                            used for filtering otherwise "call" will be used
   // outtree  - name of output filter tree or null to prevent filling outtree
   //            note: if CallFlag() is called before Calculate() set outtree=0
   if(kCS) cout << "------XPreFilter::CallFlag------" << endl;

   Int_t err  = errNoErr;
   Int_t okay = 0;

   if (intree == 0) {
      cerr << "Error: Intree is missing." << endl;
      return errGetTree;
   }//if

// Create branches for outtree
   XUnitID *unit = 0;
   XMask   *mask = 0;
   if (outtree) {
      Int_t split  = 99;
      Int_t buffer = 64000; //other size????
      unit = new XUnitID();
      mask = new XMask();

      outtree->Branch("UnitBranch", "XUnitID", &unit, buffer, split);
      outtree->Branch("MaskBranch", "XMask", &mask, buffer, split);
   }//if

// Decompose varlist
   Bool_t hasCall = kFALSE;
   Bool_t hasPVal = kFALSE;

   if (strcmp(varlist,"*") == 0) {
      hasCall = kTRUE;
      hasPVal = kTRUE;
   } else {
      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;
      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname,"call") == 0) {hasCall = kTRUE;}
         if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE;}
         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while
      delete [] dname;
   }//if

// Check leaves of intree
   if (intree[0]->FindLeaf("fUnitID") == 0) {
      cerr << "Error: Tree does not have leaf <fUnitID>." << endl;
      return errAbort;
   }//if
   if (intree[0]->FindLeaf("fCall")   == 0) hasCall = kFALSE;
   if (intree[0]->FindLeaf("fPValue") == 0) hasPVal = kFALSE;

// Check varlist setting
   if (!(hasCall || hasPVal)) {
      cout << "Warning: Cannot calculate call filter: no call tree(s) or wrong varlist."
           << endl;
      fHasCal = kFALSE;
      return errAbort;
   }//if

   // make sure that present call flags will be used for calculation
   if (hasCall && !hasPVal) fCallPValue = 1.0;

// Get number of samples and number of entries
   Int_t nentries = (Int_t)(intree[0]->GetEntries());
   fNCall = n;

   XPCall **callj = new XPCall*[fNCall];
   for (Int_t j=0; j<fNCall; j++) callj[j] = 0;

   for (Int_t j=0; j<fNCall; j++) {
      intree[j]->SetBranchAddress("CallBranch", &callj[j]);
   }//for_j

   TLeaf   *leafUnit = intree[0]->FindLeaf("fUnitID");
   TBranch *brchUnit = leafUnit->GetBranch();

// Init call arrays
   Int_t    *arrCall = 0;  //present call
   Double_t *arrPVal = 0;  //detection p-value

// Create call arrays
   if (!(arrCall = new Int_t[fNCall]))    {err = errInitMemory; goto cleanup;}
   if (!(arrPVal = new Double_t[fNCall])) {err = errInitMemory; goto cleanup;}

   // init values to NA (= -1)
   for (Int_t i=0; i<fNCall; i++) {
      arrCall[i] = -1;
      arrPVal[i] = -1;
   }//for_i

// Create array fMask
   fNMask = nentries;
   if (fMask) {delete [] fMask; fMask = 0;}
   if (!(fMask = new Int_t[nentries])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<nentries; i++) fMask[i] = 0;

// Initialize call conditions
   this->InitCallConditions();

// Read intree entries, calculate mask and fill outtree
   for (Int_t i=0; i<nentries; i++) {
      brchUnit->GetEntry(i);
      Int_t id = (Int_t)leafUnit->GetValue();

      // read entry i from tree plus friends
      for (Int_t j=0; j<fNCall; j++) {
//no         intree[j]->GetEntry(id);
         intree[j]->GetEntry(i);
         arrCall[j] = callj[j]->GetCall();
         arrPVal[j] = callj[j]->GetPValue();
      }//for_j

      fMask[i] = this->PresentCall(arrCall, arrPVal);
      okay    += fMask[i];

      if (outtree) {
         unit->SetUnitID(id);
         mask->SetFlag(fMask[i]);
         outtree->Fill();
      }//if
   }//for_i

   if (XManager::fgVerbose) {
      cout << "Call filter: <" << okay << "> genes of <" << nentries 
           << "> genes have present call for at least <" << fCallSamples
           << "> samples." << endl;
   }//if

// Cleanup
cleanup:
   if (arrPVal) {delete [] arrPVal; arrPVal = 0;}
   if (arrCall) {delete [] arrCall; arrCall = 0;}

   delete [] callj;

   return err;
}//CallFlag

//______________________________________________________________________________
Short_t XPreFilter::Gap2Mean()
{
   // gap filter gap/mean
   if(kCSa) cout << "------XPreFilter::Gap2Mean------" << endl;

   // start-index and end-index
   Int_t start = (Int_t)floor(fNData * fWindow);
   Int_t end   = fNData - start - 1; //since [i,i+1] checked in for-loop

   Int_t count = 0;
   for (Int_t i=start; i<end; i++) {
      count += ((fSorted[i+1] - fSorted[i])/fMean >= fGap2mn) ? 1 : 0;
   }//for_i

   return (Short_t)(count >= 1);
}//Gap2Mean

//______________________________________________________________________________
Short_t XPreFilter::LowerThreshold()
{
   // Lower threshold filter
   if(kCSa) cout << "------XPreFilter::LowerThreshold------" << endl;

   Short_t flag = 0;

   switch (fLowerID) {
      case 1:   //percent of samples
      case 2: { //number of samples
         Int_t count = 0;
         for (Int_t i=0; i<fNData; i++) {
            count += (fSorted[i] >= fLoThreshold) ? 1 : 0;
         }//for_i

         flag = (Short_t)(count >= fLoSamples);
         break;
      }//case

      case 3: { //trimmed mean
         // start-index and end-index for trimmed mean
         Int_t start, end;
         if (fLoSamples < 0.5) {
            start = (Int_t)floor(fNData * fLoSamples);
            end   = fNData - start;
         } else {
            if ((fNData % 2) == 0){
               start = (Int_t)floor(fNData / 2.0) - 1;
               end   = start + 2;
            } else {
               start = (Int_t)floor(fNData / 2.0);
               end   = start + 1;
           }//if
         }//if

         // calculate trimmed mean
         Int_t    trimlen = end - start;
         Double_t mean    = 0;
         for (Int_t i=start; i<end; i++) {
            mean += fSorted[i];
         }//for_i
         mean /= trimlen;

         // return flag
         flag = (Short_t)(mean >= fLoThreshold);
         break;
      }//case

      case 4: { //percentile
         Double_t qu = (fNData - 1) * fLoSamples;
         Int_t    lo = (Int_t)floor(qu);
         Int_t    hi = (Int_t)ceil(qu);

         Double_t ql = fSorted[lo];
         Double_t qh = fSorted[hi];
         Double_t qq = (ql == qh) ? 0 : qh - ql;
         Double_t qs = ql + qq * (qu - lo);

         // return flag
         flag = (Short_t)(qs >= fLoThreshold);
         break;
      }//case

      default:
         break;
   }//switch

   return flag;
}//LowerThreshold

//______________________________________________________________________________
Short_t XPreFilter::UpperThreshold()
{
   // Upper threshold filter
   if(kCSa) cout << "------XPreFilter::UpperThreshold------" << endl;

   Short_t flag = 0;

   switch (fUpperID) {
      case 1:   //percent of samples
      case 2: { //number of samples
         Int_t count = 0;
         for (Int_t i=0; i<fNData; i++) {
            count += (fSorted[i] <= fUpThreshold) ? 1 : 0;
         }//for_i

         flag = (Short_t)(count >= fUpSamples);
         break;
      }//case

      case 3: { //trimmed mean
         // start-index and end-index for trimmed mean
         Int_t start, end;
         if (fUpSamples < 0.5) {
            start = (Int_t)floor(fNData * fUpSamples);
            end   = fNData - start;
         } else {
            if ((fNData % 2) == 0){
               start = (Int_t)floor(fNData / 2.0) - 1;
               end   = start + 2;
            } else {
               start = (Int_t)floor(fNData / 2.0);
               end   = start + 1;
           }//if
         }//if

         // calculate trimmed mean
         Int_t    trimlen = end - start;
         Double_t mean    = 0;
         for (Int_t i=start; i<end; i++) {
            mean += fSorted[i];
         }//for_i
         mean /= trimlen;

         // return flag
         flag = (Short_t)(mean <= fUpThreshold);
         break;
      }//case

      case 4: { //percentile
         Double_t qu = (fNData - 1) * fUpSamples;
         Int_t    lo = (Int_t)floor(qu);
         Int_t    hi = (Int_t)ceil(qu);

         Double_t ql = fSorted[lo];
         Double_t qh = fSorted[hi];
         Double_t qq = (ql == qh) ? 0 : qh - ql;
         Double_t qs = ql + qq * (qu - lo);

         // return flag
         flag = (Short_t)(qs <= fUpThreshold);
         break;
      }//case

      default:
         break;
   }//switch

   return flag;
}//UpperThreshold

//______________________________________________________________________________
Short_t XPreFilter::QuantileHi2Lo()
{
   // Quantile filter QHi/QLo
   if(kCSa) cout << "------XPreFilter::QuantileHi2Lo------" << endl;

// Find quantile
   Double_t qu, ql, qh, qq, qs, qt;
   Int_t    lo, hi;

   qu = (fNData - 1) * fHiQ;
   lo = (Int_t)floor(qu);
   hi = (Int_t)ceil(qu);
   ql = fSorted[lo];
   qh = fSorted[hi];
   qq = (ql == qh) ? 0 : qh - ql;
   qs = ql + qq * (qu - lo);

   qu = (fNData - 1) * fLoQ;
   lo = (Int_t)floor(qu);
   hi = (Int_t)ceil(qu);
   ql = fSorted[lo];
   qh = fSorted[hi];
   qq = (ql == qh) ? 0 : qh - ql;
   qt = ql + qq * (qu - lo);
   qs = (qt != 0) ? qs/qt : 999999.9; //??allowed???

   return (Short_t)(qs >= fQRatio);
}//QuantileHi2Lo

//______________________________________________________________________________
Int_t XPreFilter::Entropy(Double_t **table)
{
   // Low entropy filter
   if(kCSa) cout << "------XPreFilter::Entropy------" << endl;

   cout << "Note: Entropy filter not yet implemented" << endl;
   return 0;
}//Entropy

//______________________________________________________________________________
Short_t XPreFilter::PresentCall(Int_t *call, Double_t *pval)
{
   // Present call filter
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XPreFilter::PresentCall------" << endl;

   Int_t count = 0;

   if (fCallPValue >= 1.0) { //call
      for (Int_t i=0; i<fNCall; i++) {
         count += (call[i] >= fCallPValue) ? 1 : 0;
      }//for_i
   } else { //pval
      for (Int_t i=0; i<fNCall; i++) {
         count += (pval[i] <= fCallPValue) ? 1 : 0;
      }//for_i
   }//if

   return (Short_t)(count >= fCallSamples);
}//PresentCall

//______________________________________________________________________________
Int_t XPreFilter::SetMinFilters(Int_t min)
{
   // Ensure that fMinFilters is not larger than number of selected filters
   if(kCS) cout << "------XPreFilter::SetMinFilters------" << endl;


   Int_t selftrs = fHasMAD + fHasCov + fHasVar + fHasDif + fHasM2m + fHasGap 
                 + (Int_t)fHasLoT + (Int_t)fHasUpT + (Int_t)fHasQua
                 + (Int_t)fHasEnt + (Int_t)fHasCal;

   return (min <= selftrs) ? min : selftrs;
}//SetMinFilters


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniFilter                                                           //
//                                                                      //
// Class for univariate filter                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUniFilter::XUniFilter()
           :XFilter()
{
   // Default UniFilter constructor
   if(kCS) cout << "---XUniFilter::XUniFilter(default)------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XUniFilter::XUniFilter(const char *name, const char *type)
           :XFilter(name, type)
{
   // Normal UniFilter constructor
   if(kCS) cout << "---XUniFilter::XUniFilter------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XUniFilter::XUniFilter(const char *name, const char *type, Double_t na)
           :XFilter(name, type, na)
{
   // Normal Filter constructor
   if(kCS) cout << "---XUniFilter::XUniFilter------" << endl;

   this->Init();
}//Constructor

//______________________________________________________________________________
XUniFilter::~XUniFilter()
{
   // Filter destructor
   if(kCS) cout << "---XUniFilter::~XUniFilter------" << endl;

}//Destructor

//______________________________________________________________________________
void XUniFilter::Init()
{
   // Initialize parameters
   if(kCS) cout << "------XUniFilter::Init------" << endl;

   fFCValue        = 0.0;
   fFCDirection    = 0;
   fStat           = 0.0;
   fPValue         = 0.0;
   fPChance        = 0.0;
   fPAdjust        = 0.0;
   fCallPValue     = 1.0;
   fCallCondition1 = "percent";
   fCallSamples1   = 100.0;
   fNCall1         = 0;
   fCallCondition2 = "percent";
   fCallSamples2   = 100.0;
   fNCall2         = 0;

   fHasStat = fHasPVal = fHasPCha = fHasPAdj = kFALSE;
   fHasFdCh = fHasUniT = fHasCall = kFALSE;
}//Init

//______________________________________________________________________________
Int_t XUniFilter::Initialize(Int_t min, Bool_t reset)
{
   // Re-initialize filters:
   // - min: minimum number of selected filter methods to satisfy
   //        min = 1: equivalent to "OR", i.e. satisfy at least one filter
   //        min = 3: equivalent to "AND", i.e. satisfy all selected filters
   //        e.g. 3 filters selected and min=2: 2 of 3 filters must be satisfied
   // - reset: clear filters, i.e. reset filter methods to kFALSE
   //   Note:  after resetting, filter methods must be initialized explicitely
   if(kCS) cout << "------XUniFilter::Initialize------" << endl;

   fMinFilters = min;

   // reset all filters
   if (reset) {
      fHasStat = fHasPVal = fHasPCha = fHasPAdj = kFALSE;
      fHasFdCh = fHasUniT = fHasCall = kFALSE;
   }//if

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XUniFilter::InitType(const char *type, Option_t *options, Int_t npars,
                  Double_t *pars)
{
   // Initialize filter type
   if(kCS) cout << "------XUniFilter::InitType------" << endl;

   if (strcmp(type, "foldchange") == 0) {
      return this->InitFoldChange(npars, pars);
   } else if (strcmp(type, "unitest") == 0) {
      return this->InitUniTest(options, npars, pars);
   } else if (strcmp(type, "call") == 0) {
      return this->InitCall(options, npars, pars);
   } else {
      cerr << "Error: UniFilter algorithm <" << type << "> not known" << endl;
   }//if

   return errInitSetting;
}//InitType

//______________________________________________________________________________
Int_t XUniFilter::InitFoldChange(Int_t npars, Double_t *pars)
{
   // Initialize fold change filter with fc = pars[0] and direction = pars[1]
   // fc: fold change value
   // direction = 0: select genes with fold change > fc AND < 1/fc
   // direction > 0: select genes with fold change > fc   - upregulated genes only
   // direction < 0: select genes with fold change < 1/fc - downregulated genes only
   if(kCS) cout << "------XUniFilter::InitFoldChange------" << endl;

   if (npars != 2) return errInitSetting;

   fFCValue     = pars[0];
   fFCDirection = (Int_t)pars[1];

   fHasFdCh = kTRUE;
   return errNoErr;
}//InitFoldChange

//______________________________________________________________________________
Int_t XUniFilter::InitUniTest(const char *varlist, Int_t npars, Double_t *pars)
{
   // Initialize UniTest filter
   // varlist is the list of variables used for filtering:
   //    "stat" - univariate statistic
   //    "pval" - p-value
   //    "pcha" - p-chance
   //    "padj" - p-value or p-chance adjusted for multiple comparisons
   //     e.g. varlist = "stat:pcha",
   //     usually, only one parameter is used, e.g. varlist = "pval"
   if(kCS) cout << "------XUniFilter::InitUniTest------" << endl;

   fHasStat = fHasPVal = fHasPCha = fHasPAdj = kFALSE;

   Int_t idx = 0;
   char *vname = new char[strlen(varlist) + 1];
   char *dname = vname;
   vname = strtok(strcpy(vname,varlist),":");
   while(vname) {
      if (strcmp(vname,"stat") == 0) {fHasStat = kTRUE; idx++;}
      if (strcmp(vname,"pval") == 0) {fHasPVal = kTRUE; idx++;}
      if (strcmp(vname,"pcha") == 0) {fHasPCha = kTRUE; idx++;}
      if (strcmp(vname,"padj") == 0) {fHasPAdj = kTRUE; idx++;}
      vname = strtok(NULL, ":");
      if (vname == 0) break;
   }//while
   delete [] dname;

   if (idx != npars) {
      cerr << "Error: Number of parameters for varlist <"      << idx  
           << "> is not equal to given number of parameters <" << npars << ">."
           << endl;
      return errAbort;
   }//if

   idx = 0;
   if (fHasStat) fStat    = pars[idx++];
   if (fHasPVal) fPValue  = pars[idx++];
   if (fHasPCha) fPChance = pars[idx++];
   if (fHasPAdj) fPAdjust = pars[idx++];

   fHasUniT = (fHasStat || fHasPVal || fHasPCha || fHasPAdj);

   return errNoErr;
}//InitUniTest

//______________________________________________________________________________
Int_t XUniFilter::InitCall(Option_t *options, Int_t npars, Double_t *pars)
{
   // Initialize present call filter
   // pars[0]: range for detection p-value is [0,1]
   //          pars[0]<1 use "detection p-value" to get all pval <= pars[0]
   //          pars[0]=1 use "detection call" to get 'P' only
   //          pars[0]=2 use "detection call" to get 'P' and 'M'
   // options must be options = "condition1:condition2"
   // condition1 = "percent": pars[1] is percent of control samples less than pvalue 
   // condition1 = "samples": pars[1] is minimum number of controls less than pvalue 
   // condition2 = "percent": pars[2] is percent of experiment samples less than pvalue 
   // condition2 = "samples": pars[2] is minimum number of experiment samples less than pvalue 
   if(kCS) cout << "------XUniFilter::InitCall------" << endl;

   if (npars != 3) return errInitSetting;

	TString optcpy  = options;
   char   *opt     = (char*)optcpy.Data();
   fCallCondition1 = strtok(opt, ":");
   fCallCondition2 = strtok(NULL, ":");

   fCallPValue   = pars[0];
   fCallSamples1 = pars[1];
   fCallSamples2 = pars[2];

   fHasCall = kTRUE;
   return errNoErr;
}//InitCall

//______________________________________________________________________________
void XUniFilter::InitCallConditions()
{
   // Initialize call conditions
   if(kCS) cout << "------XUniFilter::InitCallConditions------" << endl;

   Double_t numsamples = fCallSamples1;
   if (strcmp(fCallCondition1.Data(),"percent") == 0) {
      numsamples = ceil(fNCall1*fCallSamples1/100.0);
   }//if
   fCallSamples1 = (numsamples < fNCall1) ? numsamples : fNCall1;

   numsamples = fCallSamples2;
   if (strcmp(fCallCondition2.Data(),"percent") == 0) {
      numsamples = ceil(fNCall2*fCallSamples2/100.0);
   }//if
   fCallSamples2 = (numsamples < fNCall2) ? numsamples : fNCall2;

   return;
}//InitCallConditions

//______________________________________________________________________________
Int_t XUniFilter::Calculate(const char *infile, const char *outfile, 
                  const char *varlist, Int_t nrows, const char *sepi,
                  const char *sepo, char delim)
{
   // Calculate filter
//???   // varlist is varlist of infile, e.g. "expr:call:pval"
   // outfile exports only those rows of infile which pass the filter conditions
   // Note: since outfile has same columns as infile, sepo is equal to sepi
   if(kCS) cout << "------XUniFilter::Calculate------" << endl;

//////////////
// PROBLEM: also for log,log2,log10!!???
// maybe infile must have log data?
// BETTER: create new method Transform() in class X??? to create log-table!
//////////////

   Int_t err = errNoErr;

//TO DO!!

// Close files
//   input.close();
//   output.close();

   return err;
}//Calculate

//______________________________________________________________________________
Int_t XUniFilter::Calculate(TTree *intree, const char *leafname, TTree *outtree,
                  const char *varlist, Int_t base)
{
   // Calculate filter
   if(kCS) cout << "------XUniFilter::Calculate------" << endl;

   Int_t err  = errNoErr;
   Int_t flag = 0;
   Int_t okay = 0;

   if (intree == 0 || outtree == 0) {
      cerr << "Error: Intree and/or outtree is missing." << endl;
      return errGetTree;
   }//if

// Decompose varlist
   Bool_t hasStat, hasPVal, hasPCha, hasPAdj, hasFdCh, hasCall;
   hasStat = hasPVal = hasPCha = hasPAdj = hasFdCh = hasCall = kFALSE;

   char *vname = new char[strlen(varlist) + 1];
   char *dname = vname;
   vname = strtok(strcpy(vname,varlist),":");
   while(vname) {
      if (strcmp(vname, kUniFltr[0]) == 0) {hasStat = kTRUE;}
      if (strcmp(vname, kUniFltr[1]) == 0) {hasPVal = kTRUE;}
      if (strcmp(vname, kUniFltr[2]) == 0) {hasPCha = kTRUE;}
      if (strcmp(vname, kUniFltr[3]) == 0) {hasPAdj = kTRUE;}
      if (strcmp(vname, kUniFltr[4]) == 0) {hasFdCh = kTRUE;}
      if (strcmp(vname, kUniFltr[5]) == 0) {hasCall = kTRUE;}
      vname = strtok(NULL, ":");
      if (vname == 0) break;
   }//while
   delete [] dname;

   // if varlist defines filters to use then set fHasXXX (else use InitXXX())
   if (hasStat || hasPVal || hasPCha || hasPAdj || hasFdCh || hasCall) {
      fHasStat = hasStat;
      fHasPVal = hasPVal;
      fHasPCha = hasPCha;
      fHasPAdj = hasPAdj;
      fHasFdCh = hasFdCh;
      fHasCall = hasCall;
   }//if

// Check leaves of intree
   Bool_t hasLeafMn1  = kFALSE;
   Bool_t hasLeafMn2  = kFALSE;
   Bool_t hasLeafStat = kFALSE;
   Bool_t hasLeafPVal = kFALSE;
   Bool_t hasLeafPCha = kFALSE;
   Bool_t hasLeafPAdj = kFALSE;

   TLeaf *leaf = 0;
   if ((leaf = intree->FindLeaf("mn1"))  != 0) {hasLeafMn1  = kTRUE;}
   if ((leaf = intree->FindLeaf("mn2"))  != 0) {hasLeafMn2  = kTRUE;}
   if ((leaf = intree->FindLeaf("stat")) != 0) {hasLeafStat = kTRUE;}
   if ((leaf = intree->FindLeaf("pval")) != 0) {hasLeafPVal = kTRUE;}
   if ((leaf = intree->FindLeaf("pcha")) != 0) {hasLeafPCha = kTRUE;}
   if ((leaf = intree->FindLeaf("padj")) != 0) {hasLeafPAdj = kTRUE;}

   if (!(hasLeafMn1 || hasLeafMn2 || hasLeafStat || 
         hasLeafPVal || hasLeafPCha || hasLeafPAdj)) {
      cerr << "Error: Tree <" << intree->GetName() << "> does not have correct leaves."
           << endl;
      return errAbort;
   }//if

// Set fHasXXX to kFALSE if corresponding leaf is missing
   if (fHasFdCh && !(hasLeafMn1 && hasLeafMn2)) {
      cout << "Warning: Cannot apply fold-change filter since tree <"
           << intree->GetName() << "> has no leaf <mn1> and/or <mn2>."  << endl;
      fHasFdCh = kFALSE;
   }//if

   if (fHasStat && !hasLeafStat) {
      cout << "Warning: Cannot apply statistic filter since tree <"
           << intree->GetName() << "> has no leaf <stat>."  << endl;
      fHasStat = kFALSE;
   }//if

   if (fHasPVal && !hasLeafPVal) {
      cout << "Warning: Cannot apply pvalue filter since tree <"
           << intree->GetName() << "> has no leaf <pval>."  << endl;
      fHasPVal = kFALSE;
   }//if

   if (fHasPCha && !hasLeafPCha) {
      cout << "Warning: Cannot apply pchance filter since tree <"
           << intree->GetName() << "> has no leaf <pcha>."  << endl;
      fHasPCha = kFALSE;
   }//if

   if (fHasPAdj && !hasLeafPAdj) {
      cout << "Warning: Cannot apply padjust filter since tree <"
           << intree->GetName() << "> has no leaf <padj>."  << endl;
      fHasPAdj = kFALSE;
   }//if

   // test for call trees if InitCall() is initialized
   if (fHasCall && !(fNCall1 || fNCall2)) {
      cout << "Warning: No call trees have been selected, ignoring call filter."
           << endl;
      fHasCall = kFALSE;
   }//if

// Ensure that fMinFilters is not larger than number of selected filters
   fMinFilters = this->SetMinFilters(fMinFilters);

// Set branch status for intree
   intree->SetBranchStatus("*",  0);
   if (hasLeafMn1)  intree->SetBranchStatus("mn1Br",  1);
   if (hasLeafMn2)  intree->SetBranchStatus("mn2Br",  1);
   if (hasLeafStat) intree->SetBranchStatus("statBr", 1);
   if (hasLeafPVal) intree->SetBranchStatus("pvalBr", 1);
   if (hasLeafPCha) intree->SetBranchStatus("pchaBr", 1);
   if (hasLeafPAdj) intree->SetBranchStatus("padjBr", 1);

// Set branch address for intree
   Double_t mn1, mn2, stat, pval, pcha, padj;
   mn1 = mn2 = stat = pval = pcha = padj = 0.0;
   if (hasLeafMn1)  intree->SetBranchAddress("mn1Br",  &mn1);
   if (hasLeafMn2)  intree->SetBranchAddress("mn2Br",  &mn2);
   if (hasLeafStat) intree->SetBranchAddress("statBr", &stat);
   if (hasLeafPVal) intree->SetBranchAddress("pvalBr", &pval);
   if (hasLeafPCha) intree->SetBranchAddress("pchaBr", &pcha);
   if (hasLeafPAdj) intree->SetBranchAddress("padjBr", &padj);

   // needed for FillMaskTree(intree,...)
   intree->SetBranchStatus("fUnitID*",  1);
   XUnitID *uid = 0;
   intree->SetBranchAddress("UnitBranch", &uid);
   TBranch *brchUnit = intree->GetBranch("fUnitID"); //prevent crash SafeDelete(fDataFile)

// Get number of entries
   Int_t nentries = (Int_t)(intree->GetEntries());

// Create array fMask only if it does not exist already:
// if CallFlag() was called before then fMask contains the results
   if (!fMask) {
      fNMask = nentries;
      if (!(fMask = new Int_t[nentries])) {err = errInitMemory; goto cleanup;}
//??      for (Int_t i=0; i<nentries; i++) fMask[i] = 0;
      for (Int_t i=0; i<nentries; i++) fMask[i] = 1;
   }//if

// Read intree entries, calculate mask and fill outtree
   for (Int_t i=0; i<nentries; i++) {
//      brchUnit->GetEntry(i);
//      Int_t id = uid->GetUnitID();

      // read entry i from tree
      intree->GetEntry(i);
//no      intree->GetEntry(id);

      flag = 0;
      if (fHasStat) flag += this->Statistic(stat);
      if (fHasPVal) flag += this->PValue(pval);
      if (fHasPCha) flag += this->PChance(pcha);
      if (fHasPAdj) flag += this->PAdjust(padj);
      if (fHasFdCh) flag += this->FoldChange(mn1, mn2, base);
      if (fHasCall) flag += GetMask(i);  //add mask from CallFlag()

      fMask[i] = (Int_t)(flag >= fMinFilters);
      if (fMask[i]) okay++;
   }//for_i

// Fill outtree
   err = FillMaskTree(intree, outtree, nentries, fMask);

   if (XManager::fgVerbose) {
      cout << "Unifilter: <" << okay << "> genes of <" << nentries 
           << "> genes fulfill filter criteria." << endl;
   }//if

// Cleanup
cleanup:
//??

   return err;
}//Calculate

//______________________________________________________________________________
Int_t XUniFilter::CallFlag(Int_t n, Int_t *gid, TTree **intree,
                  const char *varlist, TTree *outtree)
{
   // Calculate call flag
   if(kCS) cout << "------XUniFilter::CallFlag------" << endl;

   Int_t err = errNoErr;

   if (intree == 0) {
      cerr << "Error: Intree is missing." << endl;
      return errGetTree;
   }//if

//?? not necessary??
// Decompose varlist
   Bool_t hasCall = kFALSE;
   Bool_t hasPVal = kFALSE;

   if (strcmp(varlist,"*") == 0) {
      hasCall = kTRUE;
      hasPVal = kTRUE;
   } else {
      char *vname = new char[strlen(varlist) + 1];
      char *dname = vname;
      vname = strtok(strcpy(vname,varlist),":");
      while(vname) {
         if (strcmp(vname,"call") == 0) {hasCall = kTRUE;}
         if (strcmp(vname,"pval") == 0) {hasPVal = kTRUE;}
         vname = strtok(NULL, ":");
         if (vname == 0) break;
      }//while
      delete [] dname;
   }//if

// Check leaves of intree
   if (intree[0]->FindLeaf("fUnitID") == 0) {
      cerr << "Error: Tree does not have leaf <fUnitID>." << endl;
      return errAbort;
   }//if
   if (intree[0]->FindLeaf("fCall")   == 0) hasCall = kFALSE;
   if (intree[0]->FindLeaf("fPValue") == 0) hasPVal = kFALSE;

// Check varlist setting
   if (!(hasCall || hasPVal)) {
      cout << "Warning: Cannot calculate call filter: no call tree(s) or wrong varlist."
           << endl;
      fHasCall = kFALSE;
      return errAbort;
   }//if

   // make sure that present call flags will be used for calculation
   if (hasCall && !hasPVal) fCallPValue = 1.0;

// Init trees
   Int_t nentries = (Int_t)(intree[0]->GetEntries());

//   TLeaf   *leafUnit  = 0;
//   TBranch *brchUnit  = 0;
//   leafUnit = intree[0]->FindLeaf("fUnitID");
//   brchUnit = leafUnit->GetBranch();

   XPCall **call = new XPCall*[n];
   for (Int_t j=0; j<n; j++) {
      call[j] = 0;
      intree[j]->SetBranchAddress("CallBranch", &call[j]);
   }//for_j

// Get number of trees in each group
   Int_t nid = 0;
   Int_t n1  = 0;
   Int_t n2  = 0;
   for (Int_t i=0; i<n; i++) {
      if (gid[i] == 1) {n1++; nid++;}
      if (gid[i] == 2) {n2++; nid++;}
   }//if
   fNCall1 = n1;
   fNCall2 = n2;

   if (n1 == 0 || n2 == 0) {
      cerr << "Error: Two groups are needed for call filter." << endl;
      delete [] call;
      return errAbort;
   }//if

   if (n1 + n2 < n) {
      cout << "Warning: Number of trees for group1 + group2 is less than number of call trees."
           << endl;
//??      return errAbort;
   }//if

   Double_t  cp   = 0.0;
   Double_t *grp1 = 0;
   Double_t *grp2 = 0;

// Create arrays grp1 and grp2
   if (!(grp1 = new Double_t[n1])) {err = errInitMemory; goto cleanup;} 
   if (!(grp2 = new Double_t[n2])) {err = errInitMemory; goto cleanup;} 

// Create array fMask
   fNMask = nentries;
   if (fMask) {delete [] fMask; fMask = 0;}
   if (!(fMask = new Int_t[nentries])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<nentries; i++) fMask[i] = 0;

// Initialize call conditions
   this->InitCallConditions();

// Read intree entries, calculate mask and fill outtree
   for (Int_t i=0; i<nentries; i++) {
//      brchUnit->GetEntry(i);
//      Int_t id = (Int_t)leafUnit->GetValue();

      // read entry i from tree plus friends
      Int_t id1 = 0;
      Int_t id2 = 0;
      for (Int_t j=0; j<n; j++) {
         if (gid[j] <= 0) continue;

         intree[j]->GetEntry(i);
//no         intree[j]->GetEntry(id);

         // get either call or detection p-value
         cp = (fCallPValue >= 1.0) ? call[j]->GetCall() : call[j]->GetPValue();

         if (gid[j] == 1) grp1[id1++] = cp;
         if (gid[j] == 2) grp2[id2++] = cp;
      }//for_j

      fMask[i] = this->PresentCall(id1, grp1, id2, grp2);
   }//for_i

// Cleanup
cleanup:
   if (grp2) {delete [] grp2; grp2 = 0;}
   if (grp1) {delete [] grp1; grp1 = 0;}

   delete [] call;

   return err;
}//CallFlag

//______________________________________________________________________________
Short_t XUniFilter::FoldChange(Double_t value1, Double_t value2, Int_t base)
{
   // Fold change filter
   if(kCSa) cout << "------XUniFilter::FoldChange------" << endl;

   // convert values from logarithm to base "base"
   value1 = (base > 0) ? ((base > 1) ? ((base > 2) ? TMath::Power(10, value1) : TMath::Power(2, value1)) : TMath::Power(TMath::E(), value1)) : value1;
   value2 = (base > 0) ? ((base > 1) ? ((base > 2) ? TMath::Power(10, value2) : TMath::Power(2, value2)) : TMath::Power(TMath::E(), value2)) : value2;

   // fold change
   Double_t fc = (value2 >= value1) ? value2/value1 : -value1/value2;

   // get filter flag for specified direction
   Short_t flag = 0;
   if      (fFCDirection == 0) flag = (Short_t)(TMath::Abs(fc) >= fFCValue);
   else if (fFCDirection >  0) flag = (Short_t)(fc >= fFCValue);
   else if (fFCDirection <  0) flag = (Short_t)(fc <= -fFCValue);

   return flag;
}//FoldChange

//______________________________________________________________________________
Short_t XUniFilter::PresentCall(Int_t n1, Double_t *grp1, Int_t n2, Double_t *grp2)
{
   // Conditions for present call filter
   // Control group grp1:    less than fCallSamples1 should be present
   // Experiment group grp2: more than fCallSamples2 should be present
   // Both coditions must be met
   // Note: Present call data will be stored as: 'P'=2, 'M'=1, 'A'=0
   if(kCSa) cout << "------XUniFilter::PresentCall------" << endl;

   Int_t count1 = 0;
   Int_t count2 = 0;

   if (fCallPValue >= 1.0) { //call
      for (Int_t i=0; i<n1; i++) {
         count1 += (grp1[i] >= fCallPValue) ? 1 : 0;
      }//for_i
      for (Int_t i=0; i<n2; i++) {
         count2 += (grp2[i] >= fCallPValue) ? 1 : 0;
      }//for_i
   } else { //pval
      for (Int_t i=0; i<n1; i++) {
         count1 += (grp1[i] <= fCallPValue) ? 1 : 0;
      }//for_i
      for (Int_t i=0; i<n2; i++) {
         count2 += (grp2[i] <= fCallPValue) ? 1 : 0;
      }//for_i
   }//if

   return (Short_t)((count1 <= fCallSamples1)  && (count2 >= fCallSamples2));
}//PresentCall

//______________________________________________________________________________
Int_t XUniFilter::SetMinFilters(Int_t min)
{
   // Ensure that fMinFilters is not larger than number of selected filters
   if(kCS) cout << "------XUniFilter::SetMinFilters------" << endl;

   Int_t selftrs = (Int_t)fHasStat + (Int_t)fHasPVal + (Int_t)fHasPCha
                 + (Int_t)fHasPAdj + (Int_t)fHasFdCh + (Int_t)fHasCall;

   return (min <= selftrs) ? min : selftrs;
}//SetMinFilters


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiFilter                                                         //
//                                                                      //
// Class for multivariate filter                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMultiFilter::XMultiFilter()
             :XFilter()
{
   // Default MultiFilter constructor
   if(kCS) cout << "---XMultiFilter::XMultiFilter(default)------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XMultiFilter::XMultiFilter(const char *name, const char *type)
             :XFilter(name, type)
{
   // Normal MultiFilter constructor
   if(kCS) cout << "---XMultiFilter::XMultiFilter------" << endl;

   this->Init();   
}//Constructor

//______________________________________________________________________________
XMultiFilter::XMultiFilter(const char *name, const char *type, Double_t na)
             :XFilter(name, type, na)
{
   // Normal Filter constructor
   if(kCS) cout << "---XMultiFilter::XMultiFilter------" << endl;

   this->Init();
}//Constructor

//______________________________________________________________________________
XMultiFilter::~XMultiFilter()
{
   // Filter destructor
   if(kCS) cout << "---XMultiFilter::~XMultiFilter------" << endl;

}//Destructor

//______________________________________________________________________________
void XMultiFilter::Init()
{
   // Initialize parameters
   if(kCS) cout << "------XMultiFilter::Init------" << endl;

}//Init

