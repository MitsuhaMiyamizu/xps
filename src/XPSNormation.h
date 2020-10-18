// File created: 08/05/2002                          last modified: 02/06/2011
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

#ifndef __XPSNormation__
#define __XPSNormation__

#ifndef __XPSProcessing__
#include "XPSProcessing.h"
#endif

class XSelector;
class XNormalizer;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormationManager                                                    //
//                                                                      //
// Class for normalization of microarray data                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormationManager: public XProcessManager {

   protected:

   protected:
      virtual Int_t     ImportDefaults(const char *infile);
      virtual Int_t     InitDefaults();
      virtual XSetting *NewSetting(const char *type, const char *infile);
      virtual XTreeSet *NewTreeSet(const char *type);

   public:
      XNormationManager();
      XNormationManager(const char *name, const char *arraytype = "",
                        Int_t verbose = kTRUE);
      virtual ~XNormationManager();

      Int_t Normalize(const char *setname, const char *method = "normalize");
      
      ClassDef(XNormationManager,1) //NormationManager
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormationSetting                                                    //
//                                                                      //
// Class for initialization of normalization parameter settings         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormationSetting: public XProcesSetting {

   protected:
      XSelector      *fSelector;    //unit selection algorithm
      XNormalizer    *fNormalizer;  //normalization algorithm

   protected:
      Int_t InitApprox(Option_t *options, Int_t npars, Double_t *pars);
      Int_t InitNormalizer(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);
      Int_t InitSelector(const char *type, Option_t *options,
               Int_t npars, Double_t *pars);

   public:
      XNormationSetting();
      XNormationSetting(const char *arraytype, const char *infile);
      virtual ~XNormationSetting();

      virtual Int_t InitAlgorithm(const char *name, const char *type,
                       Option_t *options, const char *filename,
                       Int_t npars, Double_t *pars);

      XSelector   *GetSelector()   const {return fSelector;}
      XNormalizer *GetNormalizer() const {return fNormalizer;}

      ClassDef(XNormationSetting,1) //NormationSetting
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedSet                                                           //
//                                                                      //
// Base class for normalization of microarray data                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormedSet: public XProcesSet {

   protected:
      XSelector    *fSelector;    //! unit selection algorithm
      XNormalizer  *fNormalizer;  //! normalization algorithm
//??      Double_t      fSF;          //factor to scale mean/median to targetinten

   protected:
      virtual void  AddMaskTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nunits, Int_t nflags);

      virtual Int_t ExportTreeInfo(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);

      virtual Int_t ExportExprTrees(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return 0;}
      virtual Int_t ExportMaskTrees(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return 0;}
      virtual Int_t ExportMaskTreeInfo(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return 0;}

      Int_t FillExprArray(TTree *tree, Int_t n, Int_t *idx, Double_t *arr);
      Int_t FillExprTree(const char *name, Int_t n, Int_t *idx, Double_t *arr);
      Int_t FillMaskArray(const char *name, Int_t n, Int_t *arr);
      Int_t FillMaskTree(const char *name, Int_t n, Int_t *idx, Int_t *arr);
      Int_t MeanReference(Int_t numexpr, TTree **exprtree, Int_t n, Int_t *idx, Double_t *arr);
      Int_t MedianReference(Int_t numexpr, TTree **exprtree, Int_t n, Int_t *idx, Double_t *arr);

   public:
      XNormedSet();
      XNormedSet(const char *name, const char *type);
      virtual ~XNormedSet();

      virtual Int_t Initialize(TFile *file, XSetting *setting,
                       const char *infile = "", const char *treename = "");
      virtual Int_t Normalize(const char *method) = 0;

      ClassDef(XNormedSet,1) //NormedSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedGCSet                                                         //
//                                                                      //
// Class for normalization of GeneChip oligonucleotide arrays           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormedGCSet: public XNormedSet {

   protected:

   protected:
      virtual Int_t ExportExprTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportMaskTrees(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);

   public:
      XNormedGCSet();
      XNormedGCSet(const char *name, const char *type);
      virtual ~XNormedGCSet();

      virtual Int_t Normalize(const char *method);

      ClassDef(XNormedGCSet,1) //NormedGCSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedGenomeSet                                                     //
//                                                                      //
// Class for normalization of GenomeChip oligonucleotide arrays         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormedGenomeSet: public XNormedGCSet {

   protected:

   public:
      XNormedGenomeSet();
      XNormedGenomeSet(const char *name, const char *type);
      virtual ~XNormedGenomeSet();

      ClassDef(XNormedGenomeSet,1) //NormedGenomeSet
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XNormedExonSet                                                       //
//                                                                      //
// Class for normalization of ExonChip oligonucleotide arrays           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XNormedExonSet: public XNormedGenomeSet {

   protected:

   public:
      XNormedExonSet();
      XNormedExonSet(const char *name, const char *type);
      virtual ~XNormedExonSet();

      ClassDef(XNormedExonSet,1) //NormedExonSet
};

#endif

