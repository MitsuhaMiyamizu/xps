// File created: 12/16/2002                          last modified: 02/26/2010
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

#ifndef __XPSAnalysis__
#define __XPSAnalysis__

#include "XPSSchemes.h"

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif

class XFilter;
class XPreFilter;
class XUniFilter;
class XMultiFilter;
class XAnalyser;
class XUniTester;
class XMultiTester;
class XClusterizer;
class XRegressor;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalysisManager                                                     //
//                                                                      //
// Class for statistical analysis of microarray data                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAnalysisManager: public XProcessManager {

   private:

   protected:
      virtual Int_t     ImportDefaults(const char *infile);
      virtual Int_t     InitDefaults();
      virtual XSetting *NewSetting(const char *type, const char *infile);
      virtual XTreeSet *NewTreeSet(const char *type);
      virtual XPlot    *NewPlotter(const char *name, const char *title = "");

   public:
      XAnalysisManager();
      XAnalysisManager(const char *name, const char *type = "",
                       Int_t verbose = kTRUE);
      virtual ~XAnalysisManager();

      Int_t InitSetting(Int_t min, const char *logbase = "0", Double_t neglog = 1.0);
      Int_t Analyse(const char *setname, const char *leafname, const char *outtree,
               const char *varlist = "*");
      Int_t Analyse(const char *infile, const char *outfile, const char *varlist,
               Int_t nrows, const char *sepi = "\t", const char *sepo = "\t",
               char delim = '\n');

      ClassDef(XAnalysisManager,1) //AnalysisManager
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalySetting                                                        //
//                                                                      //
// Class for initialization of parameter settings for analysis          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAnalySetting: public XProcesSetting {

   protected:
      XFilter    *fFilter;     //filter
      XAnalyser  *fAnalyser;   //analyser
      Int_t       fMinFilters; //minimum number of filters to satisfy
      TString     fLogBase;    //base of logarithm to convert tree data
      Double_t    fNegLog;     //replacement for negative data in log

   protected:
      Int_t InitPreFilter(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitUniFilter(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitMultiFilter(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitUniTest(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitMultiTest(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitClusterizer(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitRegressor(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);

   public:
      XAnalySetting();
      XAnalySetting(const char *type, const char *infile);
      virtual ~XAnalySetting();

      virtual Int_t InitAlgorithm(const char *name, const char *type,
                       Option_t *options, const char *filename,
                       Int_t npars, Double_t *pars);
      virtual void  ResetAlgorithm(const char *name, const char *type);

      Int_t   InitFilters(Int_t min, const char *logbase, Double_t neglog);

      XFilter   *GetFilter()   const {return fFilter;}
      XAnalyser *GetAnalyser() const {return fAnalyser;}
      TString    GetLogBase()  const {return fLogBase;}
      Double_t   GetNegLog()   const {return fNegLog;}

      ClassDef(XAnalySetting,1) //AnalySetting
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreFilterHeader                                                     //
//                                                                      //
// Class for storing nonspecific filter information in tree header      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XPreFilterHeader: public XTreeHeader {

   protected:
//??      Int_t      fMinFilters;    //minimum number of filters to satisfy
      Double_t   fMAD;           //cutoff median absolute deviation
      Double_t   fCov2mn;        //cutoff stdev/mean
      Double_t   fVar2mn;        //cutoff variance/mean
      Double_t   fDif2mn;        //cutoff (max-min)/mean
      Double_t   fMax2min;       //cutoff max/min
      Double_t   fGap2mn;        //cutoff gap/mean
      Double_t   fTrim;          //trim value for mean
      Double_t   fWindow;        //window size of gap filter
      TString    fLoCondition;   //condition for lower threshold
      Double_t   fLoThreshold;   //lower threshold
      Double_t   fLoSamples;     //number/percent of samples for lower threshold
      TString    fUpCondition;   //condition for upper threshold
      Double_t   fUpThreshold;   //upper threshold
      Double_t   fUpSamples;     //number/percent of samples for upper threshold
      Double_t   fLoQ;           //lowest percentile
      Double_t   fHiQ;           //highest percentile
      Double_t   fQRatio;        //cutoff ratio HiQ/LoQ
      Double_t   fEntropy;       //cutoff entropy
      Int_t      fNQuantiles;    //number of quantiles
      TString    fCallCondition; //condition for present call
      Double_t   fCallPValue;    //detection p-value
      Double_t   fCallSamples;   //number/percent of samples for present call
      Short_t    fHasMAD;        //greater 0 if apply MAD filter
      Short_t    fHasCov;        //greater 0 if apply coefficient-of-variation filter
      Short_t    fHasVar;        //greater 0 if apply variance filter
      Short_t    fHasDif;        //greater 0 if apply (max-min) filter
      Short_t    fHasM2m;        //greater 0 if apply max/min filter
      Short_t    fHasGap;        //greater 0 if apply gap filter
      Bool_t     fHasLoT;        //TRUE if apply lower threshold filter
      Bool_t     fHasUpT;        //TRUE if apply upper threshold filter
      Bool_t     fHasQua;        //TRUE if apply qunatile filter
      Bool_t     fHasEnt;        //TRUE if apply entropy filter
      Bool_t     fHasCal;        //TRUE if apply present call filter
      
   public :
      XPreFilterHeader();
      XPreFilterHeader(const char *str, Int_t treeid = 0);
      virtual ~XPreFilterHeader();

      void     SetMedAbsDev(Double_t mad)         {fMAD           = mad;}
      void     SetCov2Mean(Double_t var)          {fCov2mn        = var;}
      void     SetVar2Mean(Double_t var)          {fVar2mn        = var;}
      void     SetDif2Mean(Double_t dif)          {fDif2mn        = dif;}
      void     SetMax2Min(Double_t m2m)           {fMax2min       = m2m;}
      void     SetGap2Mean(Double_t gap)          {fGap2mn        = gap;}
      void     SetTrimValue(Double_t trim)        {fTrim          = trim;}
      void     SetWindowSize(Double_t win)        {fWindow        = win;}
      void     SetLowerCondition(const char *con) {fLoCondition   = con;}
      void     SetLowerThreshold(Double_t val)    {fLoThreshold   = val;}
      void     SetLowerSamples(Double_t val)      {fLoSamples     = val;}
      void     SetUpperCondition(const char *con) {fUpCondition   = con;}
      void     SetUpperThreshold(Double_t val)    {fUpThreshold   = val;}
      void     SetUpperSamples(Double_t val)      {fUpSamples     = val;}
      void     SetLowQuantile(Double_t loQ)       {fLoQ           = loQ;}
      void     SetHighQuantile(Double_t hiQ)      {fHiQ           = hiQ;}
      void     SetQuantileRatio(Double_t h2l)     {fQRatio        = h2l;}
      void     SetEntropy(Double_t val)           {fEntropy       = val;}
      void     SetNumberQuantiles(Int_t val)      {fNQuantiles    = val;}
      void     SetCallCondition(const char *con)  {fCallCondition = con;}
      void     SetCallPValue(Double_t pval)       {fCallPValue    = pval;}
      void     SetCallSamples(Double_t val)       {fCallSamples   = val;}
      void     SetMADFilter(Short_t val)          {fHasMAD        = val;}
      void     SetCoefOfVarFilter(Short_t val)    {fHasCov        = val;}
      void     SetVarianceFilter(Short_t val)     {fHasVar        = val;}
      void     SetDifferenceFilter(Short_t val)   {fHasDif        = val;}
      void     SetRatioFilter(Short_t val)        {fHasM2m        = val;}
      void     SetGapFilter(Short_t val)          {fHasGap        = val;}
      void     SetLoThresholdFilter(Bool_t val)   {fHasLoT        = val;}
      void     SetUpThresholdFilter(Bool_t val)   {fHasUpT        = val;}
      void     SetQuantileFilter(Bool_t val)      {fHasQua        = val;}
      void     SetEntropyFilter(Bool_t val)       {fHasEnt        = val;}
      void     SetCallFilter(Bool_t val)          {fHasCal        = val;}

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

      ClassDef(XPreFilterHeader,1) //FilterHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniFilterHeader                                                     //
//                                                                      //
// Class for storing univariate filter information in tree header       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XUniFilterHeader: public XTreeHeader {

   protected:
//??      Int_t      fMinFilters;    //minimum number of filters to satisfy
      Double_t   fFCValue;        //cutoff for fold change value
      Int_t      fFCDirection;    //direction of fold change
      Double_t   fStat;           //cutoff for univariate statistic
      Double_t   fPValue;         //cutoff for p-value
      Double_t   fPChance;        //cutoff for p-chance
      Double_t   fPAdjust;        //cutoff for adjusted p-value
      Double_t   fCallPValue;     //detection p-value
      TString    fCallCondition1; //condition for present call of group1
      Double_t   fCallSamples1;   //number/percent of samples1 for present call
      TString    fCallCondition2; //condition for present call of group2
      Double_t   fCallSamples2;   //number/percent of samples2 for present call
      Bool_t     fHasStat;        //TRUE if cutoff is univariate statistic
      Bool_t     fHasPVal;        //TRUE if cutoff is p-value
      Bool_t     fHasPCha;        //TRUE if cutoff is p-chance
      Bool_t     fHasPAdj;        //TRUE if cutoff is adjusted p-value
      Bool_t     fHasFdCh;        //TRUE if apply fold change filter
      Bool_t     fHasUniT;        //TRUE if apply unitest filter
      Bool_t     fHasCall;        //TRUE if apply present call filter
      
   public :
      XUniFilterHeader();
      XUniFilterHeader(const char *str, Int_t treeid = 0);
      virtual ~XUniFilterHeader();

      void     SetFCValue(Double_t val)           {fFCValue        = val;}
      void     SetFCDirection(Int_t dir)          {fFCDirection    = dir;}
      void     SetStatistic(Double_t val)         {fStat           = val;}
      void     SetPValue(Double_t val)            {fPValue         = val;}
      void     SetPChance(Double_t val)           {fPChance        = val;}
      void     SetPAdjust(Double_t val)           {fPAdjust        = val;}
      void     SetCallPValue(Double_t val)        {fCallPValue     = val;}
      void     SetCallCondition1(const char *con) {fCallCondition1 = con;}
      void     SetCallSamples1(Double_t val)      {fCallSamples1   = val;}
      void     SetCallCondition2(const char *con) {fCallCondition2 = con;}
      void     SetCallSamples2(Double_t val)      {fCallSamples2   = val;}
      void     SetHasStatistic(Bool_t val)        {fHasStat        = val;}
      void     SetHasPValue(Bool_t val)           {fHasPVal        = val;}
      void     SetHasPChance(Bool_t val)          {fHasPCha        = val;}
      void     SetHasPAdjust(Bool_t val)          {fHasPAdj        = val;}
      void     SetFoldChangeFilter(Bool_t val)    {fHasFdCh        = val;}
      void     SetUniTestFilter(Bool_t val)       {fHasUniT        = val;}
      void     SetCallFilter(Bool_t val)          {fHasCall        = val;}

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

      ClassDef(XUniFilterHeader,1) //UniFilterHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiFilterHeader                                                   //
//                                                                      //
// Class for storing multivariate filter information in tree header     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XMultiFilterHeader: public XTreeHeader {

   protected:
      
   public :
      XMultiFilterHeader();
      XMultiFilterHeader(const char *str, Int_t treeid = 0);
      virtual ~XMultiFilterHeader();

      ClassDef(XMultiFilterHeader,1) //MultiFilterHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUniTestHeader                                                       //
//                                                                      //
// Class for storing univariate-test information in tree header         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XUniTestHeader: public XTreeHeader {

   protected:
      Double_t   fConfLevel;    //confidence level
      Double_t   fMu;           //true value of the mean
      Int_t      fNPerm;        //number of permutations for p-chance
      TString    fAlternative;  //alternative hypothesis
      TString    fAdjustment;   //method used to adjust p-value, e.g. Bonferroni
      Bool_t     fTwoSample;    //TRUE if two-sample test
      Bool_t     fPaired;       //TRUE if paired test
      Bool_t     fEqualVar;     //TRUE if variances treated as equal
      
   public :
      XUniTestHeader();
      XUniTestHeader(const char *str, Int_t treeid = 0);
      virtual ~XUniTestHeader();

      void     SetAdjustment(const char *adj)  {fAdjustment  = adj;}
      void     SetAlternative(const char *alt) {fAlternative = alt;}
      void     SetConfidenceLevel(Double_t cl) {fConfLevel   = cl;}
      void     SetEqualVariance(Bool_t var)    {fEqualVar    = var;}
      void     SetIsPaired(Bool_t paired)      {fPaired      = paired;}
      void     SetMu(Double_t mu)              {fMu          = mu;}
      void     SetNumPerm(Int_t nperm)         {fNPerm       = nperm;}
      void     SetTwoSample(Bool_t twos)       {fTwoSample   = twos;}

      TString  GetAdjustment()       const {return fAdjustment;}
      TString  GetAlternative()      const {return fAlternative;}
      Double_t GetConfidenceLevel()  const {return fConfLevel;}
      Bool_t   GetEqualVariance()    const {return fEqualVar;}
      Bool_t   GetIsPaired()         const {return fPaired;}
      Double_t GetMu()               const {return fMu;}
      Int_t    GetNumPerm()          const {return fNPerm;}
      Bool_t   GetTwoSample()        const {return fTwoSample;}

      ClassDef(XUniTestHeader,1) //UniTestInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultiTestHeader                                                     //
//                                                                      //
// Class for storing multivariate-test information in tree header       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XMultiTestHeader: public XTreeHeader {

   protected:
      
   public :
      XMultiTestHeader();
      XMultiTestHeader(const char *str, Int_t treeid = 0);
      virtual ~XMultiTestHeader();

      ClassDef(XMultiTestHeader,1) //MultiTestHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XClusterHeader                                                       //
//                                                                      //
// Class for storing cluster analysis information in tree header        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XClusterHeader: public XTreeHeader {

   protected:
      
   public :
      XClusterHeader();
      XClusterHeader(const char *str, Int_t treeid = 0);
      virtual ~XClusterHeader();

      ClassDef(XClusterHeader,1) //ClusterHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRegressionHeader                                                    //
//                                                                      //
// Class for storing regression-test information in tree header         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
 class XRegressionHeader: public XTreeHeader {

   protected:
      
   public :
      XRegressionHeader();
      XRegressionHeader(const char *str, Int_t treeid = 0);
      virtual ~XRegressionHeader();

      ClassDef(XRegressionHeader,1) //RegressionHeader
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnalySet                                                            //
//                                                                      //
// Class for analysis of microarray data                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAnalySet: public XProcesSet {

   protected:
//?? not necessary?? fCalls, fFilters ??
      TList     *fCalls;         //! list of selected call trees
      TList     *fFilters;       //! list of selected filter trees
      TTree     *fFilterTree;    //! current filter tree
      Int_t      fMinFilters;    //minimum number of filters to satisfy
      TString    fLogBase;       //base of logarithm to convert tree data
      Double_t   fNegLog;        //replacement for negative data in log
      TString    fSchemeType;    //chip type

   protected:
      using XProcesSet::HandleOption;
      virtual Int_t HandleOption(TTree *tree, Option_t *opt);
      virtual Int_t GetCallMask(Int_t ntree, TTree **tree, Int_t n, Int_t *msk);
      virtual Int_t GetFilterMask(Int_t ntree, TTree **tree, Int_t n, Int_t *msk);
      virtual Int_t ExportFilterTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);

      Bool_t        IsFilterTree(TTree *tree);
      Int_t         LogBase();
      Int_t         SelectCall(TTree *tree, Int_t id = 9999);
      Int_t         SelectFilter(TTree *tree, Int_t id = -1);

   public:
      XAnalySet();
      XAnalySet(const char *name, const char *type);
      virtual ~XAnalySet();

      virtual Int_t Analyse(const char * /*leafname*/, const char * /*outtree*/,
                       const char * /*varlist*/) {return 0;}
      virtual Int_t Analyse(const char * /*infile*/, const char * /*outfile*/, 
                       const char * /*varlist*/, Int_t /*nrows*/, const char * /*sepi*/,
                       const char * /*sepo*/, char /*delim*/)
                                                 {return 0;}
//?? in XProcesSet??
//      virtual TTree *CopyExprTrees(TTree *tree, Int_t n, Short_t *msk, Int_t base = 0);
//      virtual TTree *CopyCallTrees(TTree *tree, Int_t n, Short_t *msk);
      virtual Int_t CopyExprTrees(Int_t ntree, TTree **fromtree, TTree **totree,
                       Int_t n, Int_t *msk, Int_t base = 0, Bool_t save = kTRUE);
      virtual Int_t CopyCallTrees(Int_t ntree, TTree **fromtree, TTree **totree,
                       Int_t n, Int_t *msk, Bool_t save = kTRUE);

      void          SetFilterTree(TTree *tree)      {fFilterTree = tree;}
      void          SetMinFilters(Int_t min)        {fMinFilters = min;}
      void          SetLogBase(const char *base)    {fLogBase    = base;}
      void          SetSchemeType(const char *type) {fSchemeType = type;}
      TString       GetLogBase()              const {return fLogBase;}
      Double_t      GetNegLog()               const {return fNegLog;}
      TString       GetSchemeType()           const {return fSchemeType;}

      ClassDef(XAnalySet,1) //AnalySet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPreFilterSet                                                        //
//                                                                      //
// Class for nonspecific filtering of microarray data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPreFilterSet: public XAnalySet {

   protected:
//?? should fFilter be stored with XTreeSet?? (then XFilter::Init() must be public)
      XPreFilter   *fFilter;     //! nonspecific filter

   protected:
//      virtual Int_t ExportFilterTrees(Int_t n, TString *names, const char *varlist,
//                       ofstream &output, const char *sep);

   public:
      XPreFilterSet();
      XPreFilterSet(const char *name, const char *type);
      virtual ~XPreFilterSet();

      using XTreeSet::AddTreeHeader;
      virtual void  AddTreeHeader(const char *treename, Int_t treeid);
      virtual Int_t Analyse(const char *leafname, const char *outtree,
                       const char *varlist);
      virtual Int_t Analyse(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");

      ClassDef(XPreFilterSet,1) //PreFilterSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnivarSet                                                           //
//                                                                      //
// Class for univariate analysis of microarray data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUnivarSet: public XAnalySet {

   protected:
      XUniTester   *fAnalyser;     //! univariate analyser
      XUniFilter   *fFilter;       //! univariate filter

   private:
      void  AddHeader(const char *treename, Int_t treeid);
      Int_t UniTest(Int_t n, TTree **tree, const char *leafname,
               const char *outtree, const char *varlist, Option_t *option);
      Int_t Filter(Int_t n, TTree **tree, const char *leafname, Int_t nc,
               TTree **calltree, const char *outtree, const char *varlist,
               Option_t *option, Int_t base);

   protected:
      virtual Int_t ExportUnivarTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
//      virtual Int_t ExportFilterTrees(Int_t n, TString *names, const char *varlist,
//                       ofstream &output, const char *sep);

   public:
      XUnivarSet();
      XUnivarSet(const char *name, const char *type);
      virtual ~XUnivarSet();

      using XTreeSet::AddTreeHeader;
      virtual void  AddTreeHeader(const char *treename, Int_t treeid);

      virtual Int_t Analyse(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t Analyse(const char *leafname, const char *outtree,
                       const char *varlist);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");

      ClassDef(XUnivarSet,1) //UnivarSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultivarSet                                                         //
//                                                                      //
// Class for multivariate analysis of microarray data                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMultivarSet: public XAnalySet {

   protected:
      XMultiTester   *fAnalyser;     //! multivariate analyser
      XMultiFilter   *fFilter;       //! multivariate filter

   private:

   protected:
      virtual Int_t ExportMultivarTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);

   public:
      XMultivarSet();
      XMultivarSet(const char *name, const char *type);
      virtual ~XMultivarSet();

      virtual Int_t Analyse(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t Analyse(const char *leafname, const char *outtree,
                       const char *varlist);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");

      ClassDef(XMultivarSet,1) //MultivarSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XClusterSet                                                          //
//                                                                      //
// Class for cluster analysis of microarray data                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XClusterSet: public XAnalySet {

   protected:
      XClusterizer   *fAnalyser;     //! cluster analyser

   private:

   protected:
      virtual Int_t ExportClusterTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);

   public:
      XClusterSet();
      XClusterSet(const char *name, const char *type);
      virtual ~XClusterSet();

      virtual Int_t Analyse(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t Analyse(const char *leafname, const char *outtree,
                       const char *varlist);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");

      ClassDef(XClusterSet,1) //ClusterSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRegressionSet                                                       //
//                                                                      //
// Class for regression analysis of microarray data                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XRegressionSet: public XAnalySet {

   protected:
      XRegressor   *fAnalyser;     //! regression analyser

   private:

   protected:
      virtual Int_t ExportRegressionTrees(Int_t n, TString *names,
                       const char *varlist, ofstream &output, const char *sep);

   public:
      XRegressionSet();
      XRegressionSet(const char *name, const char *type);
      virtual ~XRegressionSet();

      virtual Int_t Analyse(const char *infile, const char *outfile, 
                       const char *varlist, Int_t nrows, const char *sepi,
                       const char *sepo, char delim);
      virtual Int_t Analyse(const char *leafname, const char *outtree,
                       const char *varlist);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");

      ClassDef(XRegressionSet,1) //RegressionSet
};

//?????? NOT NECESSARY??????????

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XScore                                                               //
//                                                                      //
// Class containing score of univariate test                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XScore: public TObject {

   protected:
      Double32_t    fStat;         //univariate statistic
      Double32_t    fDF;           //degrees of freedom
      Double_t      fPValue;       //p-value

   public:
      XScore() {}
      virtual ~XScore() {}

      void SetStatistic(Double_t stat) {fStat   = stat;}
      void SetDF(Double_t df)          {fDF     = df;}
      void SetPValue(Double_t pval)    {fPValue = pval;}

      Double_t GetStatistic()    const {return fStat;}
      Double_t GetDF()           const {return fDF;}
      Double_t GetPValue()       const {return fPValue;}

      ClassDef(XScore,1) //Score
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGrpMn                                                               //
//                                                                      //
// Class containing group means                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGrpMn: public TObject {

   protected:
      Double32_t    fMean1;        //mean value of group1
      Double32_t    fMean2;        //mean value of group2
      Double32_t    fSE;           //standard error

   public:
      XGrpMn() {}
      virtual ~XGrpMn() {}

      void SetMean1(Double_t mn1) {fMean1 = mn1;}
      void SetMean2(Double_t mn2) {fMean2 = mn2;}
      void SetStdErr(Double_t se) {fSE   = se;}

      Double_t GetMean1()   const {return fMean1;}
      Double_t GetMean2()   const {return fMean2;}
      Double_t GetStdErr()  const {return fSE;}

      ClassDef(XGrpMn,1) //GrpMn
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XChance                                                              //
//                                                                      //
// Class containing result of random permutation                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XChance: public TObject {

   protected:
      Double_t      fPChance;        //p-chance
      Int_t         fNPerm;          //number of permutations

   public:
      XChance() {}
      virtual ~XChance() {}

      void SetPChance(Double_t pchance)    {fPChance = pchance;}
      void SetNumPermutations(Int_t nperm) {fNPerm   = nperm;}

      Double_t GetPChance()          const {return fPChance;}
      Int_t    GetNumPermutations()  const {return fNPerm;}

      ClassDef(XChance,1) //Chance
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAdjust                                                              //
//                                                                      //
// Class containing adjusted p-values                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAdjust: public TObject {

   protected:
      Double_t      fPAdjust;      //p-value adjusted for multiple comparisons

   public:
      XAdjust() {}
      virtual ~XAdjust() {}

      void SetPAdjust(Double_t padj) {fPAdjust = padj;}

      Double_t GetPAdjust()    const {return fPAdjust;}

      ClassDef(XAdjust,1) //Adjust
};

#endif

