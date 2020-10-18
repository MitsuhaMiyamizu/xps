// File created: 08/05/2002                          last modified: 07/26/2009
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2009 Dr. Christian Stratowa
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

#ifndef __XPSNormalizer__
#define __XPSNormalizer__

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormalizer                                                          //
//                                                                      //
// Base class for normalization                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormalizer: public XAlgorithm {

   protected:
      TString     fDataOpt;     //option to apply to data type
      TString     fBgrdOpt;     //option for background subtraction
      TString     fLogBase;     //logbase: 0, log, log2, log10
      TString     fMethod;      //approx method
      TString     fTies;        //approx ties
      Int_t       fNApar;       //number of approx parameters
      Double_t   *fApars;       //[fNApar] array of approx parameters
//      TString     fSUnits;       //units used for scaling: all, sel
//      TString     fSMethod;     //scaling method
//      Int_t       fSApar;       //number of scaling parameters
//      Double_t   *fSpars;       //[fNApar] array of scaling parameters
      Bool_t      fInitApprox;  //TRUE if approx is initialized

   protected:
      virtual Int_t DoNormalize(Int_t /*nin*/, const Double_t * /*xin*/,
                       const Double_t * /*yin*/, Int_t /*nout*/,
                       Double_t * /*xout*/, Double_t * /*yout*/)    {return 0;}

   public:
      XNormalizer();
      XNormalizer(const char *name, const char *type);
      virtual ~XNormalizer();

      virtual Int_t     AddArray(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/, 
                           const char * /*name*/ = "")              {return 0;}
      virtual Double_t *GetArray(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/, 
                           const char * /*name*/ = "")              {return 0;}

      virtual void  SetOptions(Option_t *opt);
      virtual Int_t InitApprox(const char *method, const char *ties,
                       Int_t npar, Double_t *pars);

      using XAlgorithm::Calculate;
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);

      void    SetDataOption(Option_t *opt) {fDataOpt = opt; fDataOpt.ToLower();}
      void    SetBgrdOption(Option_t *opt) {fBgrdOpt = opt; fBgrdOpt.ToLower();}
      void    SetLogBase(const char *lgb)  {fLogBase = lgb; fLogBase.ToLower();}
      TString GetDataOption()        const {return fDataOpt;}
      TString GetBgrdOption()        const {return fBgrdOpt;}
      TString GetLogBase()           const {return fLogBase;}

      ClassDef(XNormalizer,1) //Normalizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMeanNormalizer                                                      //
//                                                                      //
// Class using trimmed mean algorithm for normalization                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMeanNormalizer: public XNormalizer {

   protected:

   protected:
      virtual Int_t DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout);

   public:
      XMeanNormalizer();
      XMeanNormalizer(const char *name, const char *type);
      virtual ~XMeanNormalizer();

      ClassDef(XMeanNormalizer,1) //MeanNormalizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMedianNormalizer                                                    //
//                                                                      //
// Class using median as algorithm for normalization                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMedianNormalizer: public XNormalizer {

   protected:

   protected:
      virtual Int_t DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout);

   public:
      XMedianNormalizer();
      XMedianNormalizer(const char *name, const char *type);
      virtual ~XMedianNormalizer();

      ClassDef(XMedianNormalizer,1) //MedianNormalizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XKernelNormalizer                                                    //
//                                                                      //
// Class using kernel smoother as algorithm for normalization           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XKernelNormalizer: public XNormalizer {

   protected:

   protected:
      virtual Int_t DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout);

   public:
      XKernelNormalizer();
      XKernelNormalizer(const char *name, const char *type);
      virtual ~XKernelNormalizer();

      ClassDef(XKernelNormalizer,1) //KernelNormalizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLowessNormalizer                                                    //
//                                                                      //
// Class using lowess as algorithm for normalization                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XLowessNormalizer: public XNormalizer {

   protected:

   protected:
      virtual Int_t DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout);

   public:
      XLowessNormalizer();
      XLowessNormalizer(const char *name, const char *type);
      virtual ~XLowessNormalizer();

      ClassDef(XLowessNormalizer,1) //LowessNormalizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSuperNormalizer                                                     //
//                                                                      //
// Class using super smoother as algorithm for normalization            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSuperNormalizer: public XNormalizer {

   protected:

   protected:
      virtual Int_t DoNormalize(Int_t nin, const Double_t *xin, const Double_t *yin,
                       Int_t nout, Double_t *xout, Double_t *yout);

   public:
      XSuperNormalizer();
      XSuperNormalizer(const char *name, const char *type);
      virtual ~XSuperNormalizer();

      ClassDef(XSuperNormalizer,1) //SuperNormalizer
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XQuantileNormalizer                                                  //
//                                                                      //
// Class using quantile algorithm for normalization                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XQuantileNormalizer: public XNormalizer {

   protected:
      Int_t       fNData;       //number of data
      Int_t       fCount;       //counter
//not necessary?
      Int_t       fNMean1;      //number of mean values
      Int_t       fNMean2;      //number of mean values
      Double_t   *fMean1;       //[fNMean1] array of entry means
      Double_t   *fMean2;       //[fNMean2] array of entry means
      TFile      *fTmpFile;     //! temporary file 

   private:
      Int_t DoMean(Int_t n, Double_t *x);
      Int_t DoTrimmedMean(Int_t n, Double_t *x, Double_t trim);

   public:
      XQuantileNormalizer();
      XQuantileNormalizer(const char *name, const char *type);
      virtual ~XQuantileNormalizer();

      virtual Int_t     AddArray(Int_t n, Double_t *x, Int_t *msk, 
                           const char *name = "");
      virtual Double_t *GetArray(Int_t n, Double_t *x, Int_t *msk, 
                           const char *name = "");

      using XAlgorithm::Calculate;
      virtual Int_t     Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk);

      ClassDef(XQuantileNormalizer,1) //QuantileNormalizer
};

#endif

