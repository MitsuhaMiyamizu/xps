// File created: 08/05/2002                          last modified: 12/29/2010
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2011 Dr. Christian Stratowa
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

#ifndef __XPSHybridizer__
#define __XPSHybridizer__

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif

class TMEstimator;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybridizer                                                          //
//                                                                      //
// Base class for microarray hybridization algorithms                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XHybridizer: public XAlgorithm {

   protected:
      Int_t       fLength;      //length of arrays
      Double_t   *fInten1;      //[fLength] Array of intensities
      Double_t   *fStdev1;      //[fLength] Array of StdDevs
      Int_t      *fNPix1;       //[fLength] Array of pixel number/percentage
      Double_t   *fInten2;      //[fLength] Array of intensities
      Double_t   *fStdev2;      //[fLength] Array of StdDevs
      Int_t      *fNPix2;       //[fLength] Array of pixel number/percentage
      Double_t   *fArray;       //[fLength] Array to be used for calculation
      Int_t       fNDefPar;     //number of default parameters
      Bool_t      fMultichip;   //TRUE if multichip algorithm

   public:
      XHybridizer();
      XHybridizer(const char *name, const char *type);
      virtual ~XHybridizer();

      virtual Int_t CreateArray(Int_t /*length*/) {return 0;}
      virtual void  DeleteArray();
      virtual Int_t SetArray(Int_t length, Double_t *array);

      void   InitArrays(Int_t length, Double_t *inten1, Double_t *stdev1,
                Int_t *npix1, Double_t *inten2, Double_t *stdev2, Int_t *npix2);
      void   InitTreeInfo(XTreeInfo *info) {fTreeInfo = info;}

      Int_t     GetLength()    const {return fLength;}
      Double_t *GetArray()     const {return fArray;}
      Bool_t    IsMultichip()  const {return fMultichip;}

      ClassDef(XHybridizer,1) //Hybridizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBackgrounder                                                        //
//                                                                      //
// Base class for microarray background analysis algorithm              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XBackgrounder: public XHybridizer {

   protected:
      Int_t        fNRows;       //number of rows of chip-matrix
      Int_t        fNCols;       //number of columns of chip-matrix

   public:
      XBackgrounder();
      XBackgrounder(const char *name, const char *type);
      virtual ~XBackgrounder();

      virtual Double_t *AdjustIntensity(Int_t n, Double_t *inten, Double_t *bgrd,
                           Double_t *stdv);
      virtual Double_t *AdjustError(Int_t n, Double_t *sd1, Double_t *sd2);
      virtual void      Smooth(const Double_t * /*arrIn*/, Double_t * /*arrOut*/, 
                           Int_t /*numrows*/, Int_t /*numcols*/)  {}

      void  SetNumRows(Int_t numrows) {fNRows = numrows;}
      void  SetNumCols(Int_t numcols) {fNCols = numcols;}
      Int_t GetNumRows()        const {return fNRows;}
      Int_t GetNumColumns()     const {return fNCols;}

      ClassDef(XBackgrounder,1) //Backgrounder
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCallDetector                                                        //
//                                                                      //
// Base class for microarray analysis present call algorithm            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCallDetector: public XHybridizer {

   protected:
      TString     fDataOpt;     //option to apply to data type
      TString     fBgrdOpt;     //option for background subtraction

   public:
      XCallDetector();
      XCallDetector(const char *name, const char *type);
      virtual ~XCallDetector();

      virtual void  SetOptions(Option_t *opt);

      void    SetDataOption(Option_t *opt) {fDataOpt = opt; fDataOpt.ToLower();}
      void    SetBgrdOption(Option_t *opt) {fBgrdOpt = opt; fBgrdOpt.ToLower();}
      TString GetDataOption()        const {return fDataOpt;}
      TString GetBgrdOption()        const {return fBgrdOpt;}

      ClassDef(XCallDetector,1) //CallDetector
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExpressor                                                           //
//                                                                      //
// Base class for microarray expression analysis algorithm              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExpressor: public XHybridizer {

   protected:
      TString      fBgrdOpt;    //option for background subtraction
      TString      fLogBase;    //logbase: 0, log, log2, log10
      TMEstimator *fEstimator;  //! optional M-estimator
      Bool_t       fSplicexpr;  //TRUE if spliceexpressor algorithm

   public:
      XExpressor();
      XExpressor(const char *name, const char *type);
      virtual ~XExpressor();

///////////
//TO DO PM correct, e.g. pmonly
//BETTER?: XExpressor::CreateArray() {if(pmonly) else if (idealmm) etc}
//=> ev XWeightedDiff and XAvgDif not necessary, since "subtractmm"
///////////
//      virtual Int_t CreateArray(Int_t length);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Double_t *dx, Double_t *dy, Int_t *idx, Int_t nsub)
                                                            {return 0;}

      virtual Int_t SetArray(Int_t length, Double_t *array);
      virtual void  SetOptions(Option_t *opt);

      void    SetBgrdOption(Option_t *opt) {fBgrdOpt = opt; fBgrdOpt.ToLower();}
      void    SetLogBase(const char *lgb)  {fLogBase = lgb; fLogBase.ToLower();}
      void    SetOption(Option_t *opt)     {fOption  = opt; fOption.ToLower();}
      TString GetBgrdOption()        const {return fBgrdOpt;}
      TString GetLogBase()           const {return fLogBase;}

      Bool_t  IsSpliceXpressor()     const {return fSplicexpr;}

      ClassDef(XExpressor,2) //Expressor
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMultichipExpressor                                                  //
//                                                                      //
// Base class for multichip expression analysis algorithm               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMultichipExpressor: public XExpressor {

   protected:
      Double_t   *fResiduals;   //[fLength] Array of residuals

   public:
      XMultichipExpressor();
      XMultichipExpressor(const char *name, const char *type);
      virtual ~XMultichipExpressor();

      Int_t     DoMedianPolish(Int_t nrow, Int_t ncol, Double_t *inten,
                   Double_t *level, Double_t *rowmed, Double_t *colmed,
                   Double_t *residu, Double_t *weight);
//      Int_t     DoPLM(Int_t nrow, Int_t ncol, Double_t *inten,
//                   Double_t *level, Double_t *rowmed, Double_t *colmed,
//                   Double_t *residu, Double_t *weight);
      Double_t *PseudoError(Int_t nrow, Int_t ncol, Double_t *residu, Double_t *se);
      Double_t *PseudoWeight(Int_t nrow, Int_t ncol, Double_t *residu, Double_t *weight);
//      Double_t *StandardError(Int_t nrow, Int_t ncol, Double_t *residu, Double_t *se);
//      Double_t *StandardWeight(Int_t nrow, Int_t ncol, Double_t *residu, Double_t *weight);

      Double_t *GetResiduals() {return fResiduals;}

      ClassDef(XMultichipExpressor,1) //MultichipExpressor
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSpliceExpressor                                                     //
//                                                                      //
// Algorithms used for summarization and splice detection               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSpliceExpressor: public XMultichipExpressor {

   protected:

   public:
      XSpliceExpressor();
      XSpliceExpressor(const char *name, const char *type);
      virtual ~XSpliceExpressor();

      ClassDef(XSpliceExpressor,1) //SpliceExpressor
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XQualifier                                                           //
//                                                                      //
// Base class for microarray quality control algorithms                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XQualifier: public XAlgorithm {

   protected:
      Double_t             fCoefHi;        //coefficient for high cutoff
      Double_t             fCoefLo;        //coefficient for low cutoff
//      TString              fDataOpt;       //option to apply to data type
      TString              fQualOpt;       //option  to apply to data type
//?      Int_t                fNQuantiles;   //length of quantiles array
//?      Double_t            *fQuantiles;    //[fNQuantiles] Array of quantiles
      XMultichipExpressor *fExpressor;    //! expressor used

   public:
      XQualifier();
      XQualifier(const char *name, const char *type);
      virtual ~XQualifier();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Double_t *w, Int_t *msk);
      virtual Int_t SetArray(Int_t length, Double_t *array);
      virtual void  SetOptions(Option_t *opt);

      Double_t  MeanBorder(Int_t begin, Int_t end, Double_t *arr);
      void      HiLoBorder(Int_t begin, Int_t end, Double_t *arr, Short_t *msk,
                       Double_t mean, Double_t &meanhi, Double_t &meanlo);

//      void      SetDataOption(Option_t *opt) {fDataOpt = opt; fDataOpt.ToLower();}
      void      SetQualOption(Option_t *opt) {fQualOpt = opt; fQualOpt.ToLower();}
//      TString   GetDataOption()        const {return fDataOpt;}
      TString   GetQualOption()        const {return fQualOpt;}
//      Int_t     GetLength()        const {return fExpressor->GetLength();}
      Option_t *GetOption()            const {return fExpressor->GetOption();}

      ClassDef(XQualifier,1) //Qualifier
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSectorBackground                                                    //
//                                                                      //
// Class for sector background algorithm                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSectorBackground: public XBackgrounder {

   protected:

   public:
      XSectorBackground();
      XSectorBackground(const char *name, const char *type);
      virtual ~XSectorBackground();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Int_t *msk);
      virtual void  Smooth(const Double_t *arrIn, Double_t *arrOut, 
                       Int_t numrows, Int_t numcols);

      ClassDef(XSectorBackground,1) //SectorBackground
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XWeightedBackground                                                  //
//                                                                      //
// Class for weighted sector background algorithm from MAS5             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XWeightedBackground: public XSectorBackground {

   protected:

   public:
      XWeightedBackground();
      XWeightedBackground(const char *name, const char *type);
      virtual ~XWeightedBackground();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Int_t *msk);

      ClassDef(XWeightedBackground,1) //WeightedBackground
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRMABackground                                                       //
//                                                                      //
// Class for global RMA background algorithm                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XRMABackground: public XBackgrounder {

   protected:
      TString     fKernel;      //kernel for density

   private:
      Int_t ComputeParameters(Int_t n, Double_t *x, Double_t *w, Double_t *pars);
      Int_t ComputeParameters(Int_t n1, Double_t *x1, Double_t *w1,
               Int_t n2, Double_t *x2, Double_t *w2, Double_t *pars);
      void  Adjust(Int_t n, Double_t *x, Double_t *pars);

   public:
      XRMABackground();
      XRMABackground(const char *name, const char *type);
      virtual ~XRMABackground();

      virtual void  SetOptions(Option_t *opt);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Int_t *msk);

      ClassDef(XRMABackground,1) //RMABackground
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCBackground                                                        //
//                                                                      //
// Class for background algorithm based on similar GC content           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCBackground: public XBackgrounder {

   protected:

   public:
      XGCBackground();
      XGCBackground(const char *name, const char *type);
      virtual ~XGCBackground();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Int_t *msk);

      ClassDef(XGCBackground,1) //GCBackground
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMeanDifferenceCall                                                  //
//                                                                      //
// Present call algorithm based on difference of means                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMeanDifferenceCall: public XCallDetector {

   protected:
      XHybridizer *fMean;     //! algorithm used for calculation

   public:
      XMeanDifferenceCall();
      XMeanDifferenceCall(const char *name, const char *type);
      virtual ~XMeanDifferenceCall();

      using XAlgorithm::Calculate;
      virtual Int_t    Calculate(Double_t &value1, Double_t &value2, Int_t &num);
      virtual Double_t GetMean(Int_t numPar, Double_t *pars, 
                          Int_t length, Double_t *array);

      ClassDef(XMeanDifferenceCall,1) //PresentCall
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDetectionCall                                                       //
//                                                                      //
// Detection call algorithm                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDetectionCall: public XCallDetector {

   protected:
// need to be put in TStatUtils as TWilcoxTest!!!
      Double_t WilcoxTest(Int_t n, Double_t *x, Double_t mu);

   public:
      XDetectionCall();
      XDetectionCall(const char *name, const char *type);
      virtual ~XDetectionCall();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XDetectionCall,1) //DetectionCall
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMAS4Call                                                            //
//                                                                      //
// Present call algorithm of MAS4                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMAS4Call: public XCallDetector {

   protected:

   public:
      XMAS4Call();
      XMAS4Call(const char *name, const char *type);
      virtual ~XMAS4Call();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XMAS4Call,1) //MAS4Call
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDABGCall                                                            //
//                                                                      //
// Detection Above BackGround call algorithm                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDABGCall: public XCallDetector {

   protected:
      Int_t        fNMaxGC;      //! maximum GC-content, i.e. kProbeLength+1
      Int_t        fNMaxBin;     //! maximum probesize for any GC-content
      Int_t       *fNBins;       //[fNMaxGC] number of probes with certain GC-content
      Double_t   **fGCTable;     //! Background values sorted for GC-content

   protected:
      Double_t PValueFisher(Int_t n, Int_t *arrgc, Double_t *inten);
      Double_t PValuePercentile(Int_t n, Int_t *arrgc, Double_t *inten, Double_t cut);

      Double_t Intensity2PValue(Int_t gcbin, Double_t inten);

      template <typename T1> T1 ChiSqrProb(int n, T1 x);
      template <typename T1> T1 UProb(T1 x);

   public:
      XDABGCall();
      XDABGCall(const char *name, const char *type);
      virtual ~XDABGCall();

      virtual Int_t Calculate(Int_t n, Double_t *x, Int_t *msk);
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);
      virtual Int_t Calculate(Int_t /*n*/, Int_t * /*x*/, Int_t * /*msk*/)
                                                            {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Int_t * /*msk*/)
                                                            {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                       Int_t * /*msk*/)                     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/,
                       Int_t * /*idx*/, Int_t * /*msk*/)     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
                                                            {return 0;}

      ClassDef(XDABGCall,1) //DABGCall
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XINICall                                                             //
//                                                                      //
// Informative call algorithm of FARMS                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XINICall: public XCallDetector {

   protected:

   private:
      Int_t DoFARMS130(Int_t nrow, Int_t ncol, Double_t *inten, Double_t *x, 
               Double_t *y, Double_t weight = 8.0, Double_t mu = 0.0, 
               Double_t scale = 2.0, Double_t tol = 0.00001, Double_t cyc = 0.0);
      Int_t DoFARMS131(Int_t nrow, Int_t ncol, Double_t *inten, 
               Double_t *x, Double_t *y, Double_t weight = 0.5, Double_t mu = 0.0, 
               Double_t scale = 1.0, Double_t tol = 0.00001, Double_t cyc = 0.0);

   public:
      XINICall();
      XINICall(const char *name, const char *type);
      virtual ~XINICall();

      virtual Int_t SetArray(Int_t length, Double_t *array);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);

      ClassDef(XINICall,1) //INICall
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XArithmeticMean                                                      //
//                                                                      //
// Arithmetic mean expression algorithm                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XArithmeticMean: public XExpressor {

   protected:

   public:
      XArithmeticMean();
      XArithmeticMean(const char *name, const char *type);
      virtual ~XArithmeticMean();

      virtual Int_t CreateArray(Int_t length);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XArithmeticMean,1) //ArithmeticMean
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeometricMean                                                       //
//                                                                      //
// Geometric mean expression algorithm                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGeometricMean: public XArithmeticMean {

   protected:

   public:
      XGeometricMean();
      XGeometricMean(const char *name, const char *type);
      virtual ~XGeometricMean();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XGeometricMean,1) //GeometricMean
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XWeightedMean                                                        //
//                                                                      //
// Weighted mean expression algorithm                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XWeightedMean: public XArithmeticMean {

   protected:

   public:
      XWeightedMean();
      XWeightedMean(const char *name, const char *type);
      virtual ~XWeightedMean();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XWeightedMean,1) //WeightedMean
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCCorrectedMean                                                     //
//                                                                      //
// Weighted mean expression algorithm with correction for GC-content    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCCorrectedMean: public XArithmeticMean {

   protected:

   public:
      XGCCorrectedMean();
      XGCCorrectedMean(const char *name, const char *type);
      virtual ~XGCCorrectedMean();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XGCCorrectedMean,1) //GCCorrectedMean
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XWeightedDiff                                                        //
//                                                                      //
// Weighted difference expression algorithm                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XWeightedDiff: public XExpressor {

   protected:

   public:
      XWeightedDiff();
      XWeightedDiff(const char *name, const char *type);
      virtual ~XWeightedDiff();

//BETTER?: XExpressor::CreateArray() {if(pmonly) else if (idealmm) etc}
      virtual Int_t CreateArray(Int_t length);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XWeightedDiff,1) //WeightedDiff
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAvgDif                                                              //
//                                                                      //
// Average Difference expression algorithm (MAS4)                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAvgDif: public XWeightedDiff {

   protected:

   public:
      XAvgDif();
      XAvgDif(const char *name, const char *type);
      virtual ~XAvgDif();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XAvgDif,1) //AvgDif
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTukeyBiweight                                                       //
//                                                                      //
// Signal value expression based on tukey biweight algorithm (MAS5)     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTukeyBiweight: public XExpressor {

   protected:

   public:
      XTukeyBiweight();
      XTukeyBiweight(const char *name, const char *type);
      virtual ~XTukeyBiweight();

//BETTER?: XExpressor::CreateArray() {if(pmonly) else if (idealmm) etc}
      virtual Int_t CreateArray(Int_t length);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Double_t &value1, Double_t &value2, Int_t &num);

      ClassDef(XTukeyBiweight,1) //TukeyBiweight
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMedianPolish                                                        //
//                                                                      //
// Median polish expression algorithm used for RMA                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMedianPolish: public XMultichipExpressor {

   protected:

   public:
      XMedianPolish();
      XMedianPolish(const char *name, const char *type);
      virtual ~XMedianPolish();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);

      ClassDef(XMedianPolish,1) //MedianPolish
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFARMS                                                               //
//                                                                      //
// Factor Analysis for Robust Microarray Summarization (Hochreiter S)   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XFARMS: public XMultichipExpressor {

   protected:

   private:
      Int_t DoFARMS130(Int_t nrow, Int_t ncol, Double_t *inten, Double_t *x, 
               Double_t *y, Double_t weight = 8.0, Double_t mu = 0.0, 
               Double_t scale = 2.0, Double_t tol = 0.00001, Double_t cyc = 0.0);
      Int_t DoFARMS131(Int_t nrow, Int_t ncol, Double_t *inten, 
               Double_t *x, Double_t *y, Double_t weight = 0.5, Double_t mu = 0.0, 
               Double_t scale = 1.0, Double_t tol = 0.00001, Double_t cyc = 0.0,
               Bool_t weighted = kTRUE);

   public:
      XFARMS();
      XFARMS(const char *name, const char *type);
      virtual ~XFARMS();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);

      ClassDef(XFARMS,1) //FARMS
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDFW                                                                 //
//                                                                      //
// Distribution Free Weighted Fold Change (Chen Z)                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDFW: public XMultichipExpressor {

   protected:

   private:
      Double_t *Weight(Int_t n, const Double_t *arr, Double_t *w);

   public:
      XDFW();
      XDFW(const char *name, const char *type);
      virtual ~XDFW();

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);

      ClassDef(XDFW,1) //DFW
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFIRMA                                                               //
//                                                                      //
// Finding Isoforms using Robust Multichip Analysis (Purdom E)          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XFIRMA: public XSpliceExpressor {

   protected:

   public:
      XFIRMA();
      XFIRMA(const char *name, const char *type);
      virtual ~XFIRMA();

//      using XAlgorithm::Calculate;
      using XExpressor::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Double_t *z,
                       Double_t *dx, Double_t *dy, Int_t *idx, Int_t nsub);

      ClassDef(XFIRMA,1) //FIRMA
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRMAQualifier                                                        //
//                                                                      //
// Class for GeneChip RMA quality control algorithms                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XRMAQualifier: public XQualifier {

   protected:

   public:
      XRMAQualifier();
      XRMAQualifier(const char *name, const char *type);
      virtual ~XRMAQualifier();

      ClassDef(XRMAQualifier,1) //RMAQualifier
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPLMQualifier                                                        //
//                                                                      //
// Class for GeneChip PLM quality control algorithms                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPLMQualifier: public XQualifier {

   protected:

   public:
      XPLMQualifier();
      XPLMQualifier(const char *name, const char *type);
      virtual ~XPLMQualifier();

      ClassDef(XPLMQualifier,1) //PLMQualifier
};

#endif

