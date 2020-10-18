// File created: 12/16/2002                          last modified: 11/11/2007
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2008 Dr. Christian Stratowa
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

#ifndef __XPSFilter__
#define __XPSFilter__

#include <Riostream.h>
#include <cmath>

#include "TMath.h"
#include "TNamed.h"
#include "TString.h"
#include "TTree.h"

#include "TStat.h"

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif

extern const char *kPreFltr[];
extern const char *kUniFltr[];
extern const char *kMultiFltr[];
/*//?? NOTE: EV BETTER?? since it is easier to add new filter classes!
class XVariationFilter: public XFilter 
class XLowerThresholdFilter: public XFilter 
class XUpperThresholdFilter: public XFilter 
class XQuantileFilter: public XFilter 
class XEntropyFilter: public XFilter 
class XCallFilter: public XFilter 
class XFoldChangeFilter: public XFilter 
class XUniTestFilter: public XFilter
etc
in class ??? use GetMask() of different filters to get final filter tree
// or: filter tree contains masks of different filters as branches?
// Replace XPreFilterSet with general XFilterSet where different filters
// can be added, e.g. XVariationFilter and XUniTestFilter, however, need
// to check if XUniTestFilter can be applied with current data!!
//BUT PROBLEM: e.g. XUniTestFilter can only be applied to tree.stt while
// e.g. XLowerThresholdFilter can only be applied to tree friends tree.sup!
*/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFilter                                                              //
//                                                                      //
// Base class for filter                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XFilter: public XAlgorithm {

   protected:
      Int_t      fNData;         //number of samples
      Double_t  *fSorted;        //[fNData] sorted array of expr for current row 
      Double_t   fMean;          //mean value of current row
      Double_t   fTrim;          //trim value for mean
      Double_t   fVar;           //variance of current row
      Double_t   fMin;           //minimum value of current row
      Double_t   fMax;           //maximum value of current row
      Double_t   fEpsilon;       //value to be added to prevent zero division
      Int_t      fMinFilters;    //minimum number of filters to satisfy
      Int_t      fNMask;         //size of final filter array fMask
      Int_t     *fMask;          //[fNMask] mask array of filter result
      
   protected:
      Int_t MeanVarMinMax(Int_t n, Double_t *arr);
      Int_t FillMaskTree(TTree *unittree, TTree *masktree, Int_t n, Int_t *arr);
      
   public:
      XFilter();
      XFilter(const char *name, const char *type = "ftr");
      XFilter(const char *name, const char *type, Double_t na);
      virtual ~XFilter();

      void    Init(); //called in constructor => called in TBrowser if stored in TFile
      virtual Int_t Initialize(Int_t min, Bool_t reset = kFALSE);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(const char * /*infile*/, const char * /*outfile*/, 
                       const char * /*varlist*/, Int_t /*nrows*/,
                       const char * /*sepi*/, const char * /*sepo*/,
                       char /*delim*/)                                {return 0;}
      virtual Int_t Calculate(TTree * /*intree*/, const char * /*leafname*/,
                       TTree * /*outtree*/, const char * /*varlist*/) {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, TTree ** /*intree*/, const char * /*leafname*/,
                       TTree * /*outtree*/, const char * /*varlist*/) {return 0;}
      virtual Int_t CallFlag(Int_t /*n*/, TTree ** /*intree*/, const char * /*varlist*/,
                       TTree * /*outtree*/)                           {return 0;}
      virtual Int_t CallFlag(Int_t /*n*/, Int_t * /*gid*/, TTree ** /*intree*/,
                       const char * /*varlist*/, TTree * /*outtree*/) {return 0;}

      Int_t  GetMask(Int_t i) const {return fMask ? fMask[i] : -1;}
      Int_t *GetMask()        const {return fMask;}

      ClassDef(XFilter,1) //Filter
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreFilter                                                           //
//                                                                      //
// Class for nonspecific filter                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPreFilter: public XFilter {

   protected:
      Double_t   fMAD;           //cutoff median absolute deviation
      Double_t   fCov2mn;        //cutoff stdev/mean
      Double_t   fVar2mn;        //cutoff variance/mean
      Double_t   fDif2mn;        //cutoff (max-min)/mean
      Double_t   fMax2min;       //cutoff max/min
      Double_t   fGap2mn;        //cutoff gap/mean
      Double_t   fWindow;        //window size of gap filter
      TString    fLoCondition;   //condition for lower threshold
      Int_t      fLowerID;       //id for lower condition
      Double_t   fLoThreshold;   //lower threshold
      Double_t   fLoSamples;     //number/percent of samples for lower threshold
      TString    fUpCondition;   //condition for upper threshold
      Int_t      fUpperID;       //id for upper condition
      Double_t   fUpThreshold;   //upper threshold
      Double_t   fUpSamples;     //number/percent of samples for upper threshold
      Double_t   fLoQ;           //lowest percentile
      Double_t   fHiQ;           //highest percentile
      Double_t   fQRatio;        //cutoff ratio HiQ/LoQ
      Double_t   fEntropy;       //cutoff entropy
      Int_t      fNQuantiles;    //number of quantiles
      TString    fCallCondition; //condition for present call
//??      Int_t      fCallID;        //id for present call condition
      Double_t   fCallPValue;    //detection p-value
      Double_t   fCallSamples;   //number/percent of samples for present call
      Int_t      fNCall;         //number of present call data
      Short_t    fHasMAD;        //greater 0 if apply MAD filter
      Short_t    fHasCov;        //greater 0 if apply coefficient-of-variation filter
      Short_t    fHasVar;        //greater 0 if apply variance filter
      Short_t    fHasDif;        //greater 0 if apply (max-min) filter
      Short_t    fHasM2m;        //greater 0 if apply max/min filter
      Short_t    fHasGap;        //greater 0 if apply gap filter
      Bool_t     fHasLoT;        //TRUE if apply lower threshold filter
      Bool_t     fHasUpT;        //TRUE if apply upper threshold filter
      Bool_t     fHasQua;        //TRUE if apply quantile filter
      Bool_t     fHasEnt;        //TRUE if apply entropy filter
      Bool_t     fHasCal;        //TRUE if apply present call filter
      
   protected:
      Int_t InitVariation(const char *varlist, Int_t npars, Double_t *pars);
      Int_t InitLowerThreshold(const char *condition, Int_t npars, Double_t *pars);
      Int_t InitUpperThreshold(const char *condition, Int_t npars, Double_t *pars);
      Int_t InitQuantile(Int_t npars, Double_t *pars);
      Int_t InitEntropy(Int_t npars, Double_t *pars);
      Int_t InitCall(const char *condition, Int_t npars, Double_t *pars);

      void    InitThresholdConditions();
      void    InitCallConditions();
      Short_t MAD();                    //inline
      Short_t CoefficientOfVariation(); //inline
      Short_t Variance2Mean();          //inline
      Short_t RatioMax2Min();           //inline
      Short_t Difference2Mean();        //inline
      Short_t Gap2Mean();
      Short_t LowerThreshold();
      Short_t UpperThreshold();
      Short_t QuantileHi2Lo();
      Int_t   Entropy(Double_t **table);
      Short_t PresentCall(Int_t *call, Double_t *pval);
      Int_t   SetMinFilters(Int_t min);
      
   public:
      XPreFilter();
      XPreFilter(const char *name, const char *type = "pfr");
      XPreFilter(const char *name, const char *type, Double_t na);
      virtual ~XPreFilter();

      void  Init(); //called in constructor => called in TBrowser if stored in TFile
      virtual Int_t Initialize(Int_t min, Bool_t reset = kFALSE);
      virtual Int_t InitType(const char *type, Option_t *options,
                       Int_t npars, Double_t *pars);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t Calculate(Int_t n, TTree **intree, const char *leafname,
                       TTree *outtree, const char *varlist);

      using XFilter::CallFlag;
      virtual Int_t CallFlag(Int_t n, TTree **intree, const char *varlist,
                       TTree *outtree);

      Double_t GetMedAbsDev()               const {return fMAD;}
      Double_t GetCov2Mean()                const {return fCov2mn;}
      Double_t GetVar2Mean()                const {return fVar2mn;}
      Double_t GetDif2Mean()                const {return fDif2mn;}
      Double_t GetMax2Min()                 const {return fMax2min;}
      Double_t GetGap2Mean()                const {return fGap2mn;}
      Double_t GetTrimValue()               const {return fTrim;}
      Double_t GetWindowSize()              const {return fWindow;}
      TString  GetLowerCondition()          const {return fLoCondition;}
      Double_t GetLowerThreshold()          const {return fLoThreshold;}
      Double_t GetLowerSamples()            const {return fLoSamples;}
      TString  GetUpperCondition()          const {return fUpCondition;}
      Double_t GetUpperThreshold()          const {return fUpThreshold;}
      Double_t GetUpperSamples()            const {return fUpSamples;}
      Double_t GetLowQuantile()             const {return fLoQ;}
      Double_t GetHighQuantile()            const {return fHiQ;}
      Double_t GetQuantileRatio()           const {return fQRatio;}
      Double_t GetEntropy()                 const {return fEntropy;}
      Int_t    GetNumberQuantiles()         const {return fNQuantiles;}
      TString  GetCallCondition()           const {return fCallCondition;}
      Double_t GetCallPValue()              const {return fCallPValue;}
      Double_t GetCallSamples()             const {return fCallSamples;}
      Bool_t   HasMADFilter()               const {return fHasMAD;}
      Bool_t   HasCoefOfVarFilter()         const {return fHasCov;}
      Bool_t   HasVarianceFilter()          const {return fHasVar;}
      Bool_t   HasDifferenceFilter()        const {return fHasDif;}
      Bool_t   HasRatioFilter()             const {return fHasM2m;}
      Bool_t   HasGapFilter()               const {return fHasGap;}
      Bool_t   HasLoThresholdFilter()       const {return fHasLoT;}
      Bool_t   HasUpThresholdFilter()       const {return fHasUpT;}
      Bool_t   HasQuantileFilter()          const {return fHasQua;}
      Bool_t   HasEntropyFilter()           const {return fHasEnt;}
      Bool_t   HasCallFilter()              const {return fHasCal;}

      ClassDef(XPreFilter,1) //PreFilter
};

//______________________________________________________________________________
inline Short_t XPreFilter::MAD()
{
   // Median absolute deviation filter

   return (Short_t)(TStat::MAD(fNData, fSorted) >= fMAD);
}//MAD

//______________________________________________________________________________
inline Short_t XPreFilter::CoefficientOfVariation()
{
   // Coefficient of variation filter

   return (Short_t)(TMath::Sqrt(fVar)/fMean >= fCov2mn);
}//CoefficientOfVariation

//______________________________________________________________________________
inline Short_t XPreFilter::Variance2Mean()
{
   // Variance filter

   return (Short_t)(fVar/fMean >= fVar2mn);
}//Variance2Mean

//______________________________________________________________________________
inline Short_t XPreFilter::RatioMax2Min()
{
   // Ratio filter max/min

   return (Short_t)(fMax/fMin >= fMax2min);
}//RatioMax2Min

//______________________________________________________________________________
inline Short_t XPreFilter::Difference2Mean()
{
   // Difference filter (max-min)/mean

   return (Short_t)((fMax - fMin)/fMean >= fDif2mn);
}//Difference2Mean


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniFilter                                                           //
//                                                                      //
// Class for univariate filter                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUniFilter: public XFilter {

   protected:
      Double_t   fFCValue;        //cutoff for fold change value
      Int_t      fFCDirection;    //direction of fold change
      Double_t   fStat;           //cutoff for univariate statistic
      Double_t   fPValue;         //cutoff for p-value
      Double_t   fPChance;        //cutoff for p-chance
      Double_t   fPAdjust;        //cutoff for adjusted p-value
      Double_t   fCallPValue;     //detection p-value
      TString    fCallCondition1; //condition for present call of group1
      Double_t   fCallSamples1;   //number/percent of samples1 for present call
      Int_t      fNCall1;         //number of present call data for group1
      TString    fCallCondition2; //condition for present call of group2
      Double_t   fCallSamples2;   //number/percent of samples2 for present call
      Int_t      fNCall2;         //number of present call data for group2
      Bool_t     fHasStat;        //TRUE if cutoff is univariate statistic
      Bool_t     fHasPVal;        //TRUE if cutoff is p-value
      Bool_t     fHasPCha;        //TRUE if cutoff is p-chance
      Bool_t     fHasPAdj;        //TRUE if cutoff is adjusted p-value
      Bool_t     fHasFdCh;        //TRUE if apply fold change filter
      Bool_t     fHasUniT;        //TRUE if apply unitest filter
      Bool_t     fHasCall;        //TRUE if apply present call filter
      
   protected:
      Int_t InitFoldChange(Int_t npars, Double_t *pars);
      Int_t InitUniTest(const char *varlist, Int_t npars, Double_t *pars);
      Int_t InitCall(Option_t *options, Int_t npars, Double_t *pars);

      void    InitCallConditions();
      Short_t Statistic(Double_t stat); //inline
      Short_t PValue(Double_t pval);    //inline
      Short_t PChance(Double_t pcha);   //inline
      Short_t PAdjust(Double_t padj);   //inline
      Short_t FoldChange(Double_t value1, Double_t value2, Int_t base);
      Short_t PresentCall(Int_t n1, Double_t *grp1, Int_t n2, Double_t *grp2);
      Int_t   SetMinFilters(Int_t min);
      
   public:
      XUniFilter();
      XUniFilter(const char *name, const char *type = "ufr");
      XUniFilter(const char *name, const char *type, Double_t na);
      virtual ~XUniFilter();

      void  Init(); //called in constructor => called in TBrowser if stored in TFile
      virtual Int_t Initialize(Int_t min, Bool_t reset = kFALSE);
      virtual Int_t InitType(const char *type, Option_t *options,
                       Int_t npars, Double_t *pars);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t Calculate(TTree *intree, const char *leafname,
                       TTree *outtree, const char *varlist, Int_t base = 0);

      using XFilter::CallFlag;
      virtual Int_t CallFlag(Int_t n, Int_t *gid, TTree **intree,
                       const char *varlist, TTree *outtree);

      Double_t GetFCValue()                 const {return fFCValue;}
      Int_t    GetFCDirection()             const {return fFCDirection;}
      Double_t GetStatistic()               const {return fStat;}
      Double_t GetPValue()                  const {return fPValue;}
      Double_t GetPChance()                 const {return fPChance;}
      Double_t GetPAdjust()                 const {return fPAdjust;}
      Double_t GetCallPValue()              const {return fCallPValue;}
      TString  GetCallCondition1()          const {return fCallCondition1;}
      Double_t GetCallSamples1()            const {return fCallSamples1;}
      TString  GetCallCondition2()          const {return fCallCondition2;}
      Double_t GetCallSamples2()            const {return fCallSamples2;}
      Bool_t   HasStatistic()               const {return fHasStat;}
      Bool_t   HasPValue()                  const {return fHasPVal;}
      Bool_t   HasPChance()                 const {return fHasPCha;}
      Bool_t   HasPAdjust()                 const {return fHasPAdj;}
      Bool_t   HasFoldChangeFilter()        const {return fHasFdCh;}
      Bool_t   HasUniTestFilter()           const {return fHasUniT;}
      Bool_t   HasCallFilter()              const {return fHasCall;}

      ClassDef(XUniFilter,1) //UniFilter
};

//______________________________________________________________________________
inline Short_t XUniFilter::PValue(Double_t pval)
{
   // PValue filter

   return (Short_t)(pval <= fPValue);
}//PValue

//______________________________________________________________________________
inline Short_t XUniFilter::PChance(Double_t pcha)
{
   // PChance filter

   return (Short_t)(pcha <= fPChance);
}//PChance

//______________________________________________________________________________
inline Short_t XUniFilter::PAdjust(Double_t padj)
{
   // PAdjust filter

   return (Short_t)(padj <= fPAdjust);
}//PAdjust

//______________________________________________________________________________
inline Short_t XUniFilter::Statistic(Double_t stat)
{
   // Statistic filter

   return (Short_t)(stat >= fStat);
}//Statistic


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiFilter                                                         //
//                                                                      //
// Class for multivariate filter                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMultiFilter: public XFilter {

   protected:
      
   protected:
      
   public:
      XMultiFilter();
      XMultiFilter(const char *name, const char *type = "mfr");
      XMultiFilter(const char *name, const char *type, Double_t na);
      virtual ~XMultiFilter();

      void  Init(); //called in constructor

      ClassDef(XMultiFilter,1) //MultiFilter
};

#endif

