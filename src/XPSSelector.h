// File created: 08/05/2002                          last modified: 09/14/2007
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

#ifndef __XPSSelector__
#define __XPSSelector__

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSelector                                                            //
//                                                                      //
// Base class for selection of units used for normalization             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSelector: public XAlgorithm {

   protected:
      TString        fOption;      //option

   public:
      XSelector();
      XSelector(const char *name, const char *type);
      XSelector(const XSelector &selector);
      XSelector& operator=(const XSelector& rhs);
      virtual ~XSelector();

      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk); 
      virtual Int_t Calculate(Double_t &/*value1*/, Double_t &/*value2*/, Int_t &/*num*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Int_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                       Int_t * /*msk*/)                      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/,
                       Int_t * /*idx*/, Int_t * /*msk*/)     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
                                                             {return 0;}

      virtual void   SetOptions(Option_t *opt);
      virtual Int_t *SetMask(Int_t n, Int_t *arr);

      Int_t     SetFlag(Int_t n, Int_t *arr);
      void      SetOption(Option_t *opt) {fOption = opt; fOption.ToLower();}
      Option_t *GetOption()        const {return fOption.Data();}

      ClassDef(XSelector,1) //Selector
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRankSelector                                                        //
//                                                                      //
// Class for selection of units based on ranking algorithm              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XRankSelector: public XSelector {

   protected:

   protected:
      Double_t Cutoff(Int_t n, const Double_t *arr, Double_t max,
                  Bool_t showGraph = kFALSE);

   public:
      XRankSelector();
      XRankSelector(const char *name, const char *type);
      virtual ~XRankSelector();

      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk); 
      virtual Int_t Calculate(Double_t &/*value1*/, Double_t &/*value2*/, Int_t &/*num*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Int_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                       Int_t * /*msk*/)                      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/,
                       Int_t * /*idx*/, Int_t * /*msk*/)     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
                                                             {return 0;}

      ClassDef(XRankSelector,1) //RankSelector
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProbeSelector                                                       //
//                                                                      //
// Class for selection of probes                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProbeSelector: public XSelector {

   protected:

   public:
      XProbeSelector();
      XProbeSelector(const char *name, const char *type);
      XProbeSelector(const XProbeSelector &selector);
      XProbeSelector& operator=(const XProbeSelector& rhs);
      virtual ~XProbeSelector();

      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk); 
      virtual Int_t Calculate(Double_t &/*value1*/, Double_t &/*value2*/, Int_t &/*num*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Int_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                       Int_t * /*msk*/)                      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/,
                       Int_t * /*idx*/, Int_t * /*msk*/)     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
                                                             {return 0;}

      virtual Int_t *SetGenomeMask(Int_t n, Int_t *arr);
      virtual Int_t *SetExonMask(Int_t n, Int_t *arr);

      ClassDef(XProbeSelector,1) //ProbeSelector
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnitSelector                                                        //
//                                                                      //
// Class for selection of units                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//old class XUnitSelector: public XSelector {
class XUnitSelector: public XProbeSelector {

   protected:

   public:
      XUnitSelector();
      XUnitSelector(const char *name, const char *type);
      XUnitSelector(const XUnitSelector &selector);
      XUnitSelector& operator=(const XUnitSelector& rhs);
      virtual ~XUnitSelector();

      virtual Int_t Calculate(Int_t n, Int_t *x, Int_t *msk); 

      virtual Int_t Calculate(Double_t &/*value1*/, Double_t &/*value2*/, Int_t &/*num*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                       Int_t * /*msk*/)                      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/,
                       Int_t * /*idx*/, Int_t * /*msk*/)     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
                                                             {return 0;}

      ClassDef(XUnitSelector,1) //UnitSelector
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUserSelector                                                        //
//                                                                      //
// Class for selection of units by importing a mask file                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUserSelector: public XSelector {

   protected:
      TString        fInput;      //input option

   protected:
      Int_t Import(const char *infile, Int_t n, Int_t *msk, char delim = '\n');

   public:
      XUserSelector();
      XUserSelector(const char *name, const char *type);
      virtual ~XUserSelector();

      virtual void  SetOptions(Option_t *opt);
      virtual Int_t Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk); 
      virtual Int_t Calculate(Double_t &/*value1*/, Double_t &/*value2*/, Int_t &/*num*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Int_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/)
                                                             {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                       Int_t * /*msk*/)                      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/,
                       Int_t * /*idx*/, Int_t * /*msk*/)     {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/,
                      Int_t * /*idx*/, Int_t * /*msk*/)      {return 0;}
      virtual Int_t Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
                                                             {return 0;}

      ClassDef(XUserSelector,1) //UserSelector
};

#endif

