// File created: 11/16/2007                          last modified: 11/24/2007
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

#ifndef __XProjectHandler__
#define __XProjectHandler__

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

//classes for database
class XLoginInfo;
class XProjectInfo;
class XAuthorInfo;
class XDatasetInfo;
class XSourceInfo;
class XSampleInfo;
class XCellLineInfo;
class XPrimaryCellInfo;
class XTissueInfo;
class XBiopsyInfo;
class XArrayInfo;
class XHybInfo;
class XHybridizationList;
class XTreatmentInfo;
class XTreatmentList;

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHandler                                                             //
//                                                                      //
// Base class for handling datatypes for database                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class XHandler {

   protected:
      TNamed   *fNamed;    //! to allow multiple inheritance
      TList    *fList;     //! list containing datatypes

   protected:
      void  Add(TObject *obj);

   public:
      XHandler();
      XHandler(const char *name, const char *title = "");
      virtual ~XHandler();

      virtual Int_t BeginTransaction(const char *name = "") {return 0;}
      virtual Int_t CommitTransaction()                     {return 0;}

      void  SetHandlerName(const char *name)      {fNamed->SetName(name);}
      void  SetHandlerTitle(const char *title="") {fNamed->SetTitle(title);}     

      const char *GetHandlerName()  const {return fNamed->GetName();}
      const char *GetHandlerTitle() const {return fNamed->GetTitle();}

      ClassDef(XHandler,1) //Handler
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProjectHandler                                                      //
//                                                                      //
// Class for handling datatypes for database                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class XProjectHandler: public XHandler {

   protected:
      XHybridizationList *fHybridizations;   //List of hybridizations
      XTreatmentList     *fTreatments;       //List of treatments

   protected:
      void  AddHybridization(XHybInfo *info);
      void  AddTreatment(XTreatmentInfo *info);

   public:
      XProjectHandler();
      XProjectHandler(const char *name, const char *title = "");
      virtual ~XProjectHandler();

      virtual TString LoginInfo(XLoginInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    LoginInfo(const char *userID, const char *password,
                         Bool_t replace = kFALSE);
      virtual TString ProjectInfo(XProjectInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    ProjectInfo(const char *name, Long_t date, const char *type = "",
                         const char *description = "", const char *comment = "",
                         Bool_t replace = kFALSE);
      virtual TString AuthorInfo(XAuthorInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    AuthorInfo(const char *lastname, const char *firstname,
                         const char *type, const char *company, const char *department = "",
                         const char *mail ="", const char *phone ="",
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString DatasetInfo(XDatasetInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    DatasetInfo(const char *name, const char *type,
                         const char *sample, const char *submitter, Long_t date,
                         const char *description = "", const char *comment = "",
                         Bool_t replace = kFALSE);
      virtual TString SourceInfo(XSourceInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    SourceInfo(const char *name, const char *type,
                         const char *species, const char *subspecies,
                         const char *description = "", const char *comment = "",
                         Bool_t replace = kFALSE);
      virtual TString SampleInfo(XSampleInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    SampleInfo(const char *name, const char *type, const char *sex,
                         const char *pheno ="", const char *geno = "",
                         const char *extract = "", Bool_t isXeno = kFALSE,
                         const char *xenostrain = "", const char *xenosex = "",
                         Double_t xenoage = 0.0, const char *xageunits = "",
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString CellLineInfo(XCellLineInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    CellLineInfo(const char *name, const char *type,
                         const char *parent = "", const char *atcc = "",
                         const char *mod = "", const char *sex = "",
                         const char *pheno = "", const char *geno = "",
                         const char *extract = "", Bool_t isXeno = kFALSE,
                         const char *xenostrain = "", const char *xenosex = "",
                         Double_t xenoage = 0, const char *xageunits = "years",
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString PrimaryCellInfo(XPrimaryCellInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    PrimaryCellInfo(const char *name, const char *type, Long_t date, 
                         const char *description = "", const char *sex = "",
                         const char *pheno = "", const char *geno = "",
                         const char *extract = "", Bool_t isXeno = kFALSE,
                         const char *xenostrain = "", const char *xenosex = "",
                         Double_t xenoage = 0, const char *xageunits = "years",
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString TissueInfo(XTissueInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    TissueInfo(const char *name, const char *type = "",
                         const char *development = "", const char *morphology = "",
                         const char *disease = "", const char *stage = "",
                         Double_t donorage = 0, const char *ageunits = "",
                         const char *status = "", const char *sex = "",
                         const char *pheno = "", const char *geno = "",
                         const char *extract = "", Bool_t isXeno = kFALSE,
                         const char *xenostrain = "", const char *xenosex = "",
                         Double_t xenoage = 0, const char *xageunits = "years",
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString BiopsyInfo(XBiopsyInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    BiopsyInfo(const char *name, const char *type = "",
                         const char *morphology = "", const char *disease = "",
                         const char *stage = "", Double_t donorage = 0,
                         const char *ageunits = "", const char *status = "",
                         const char *sex = "", const char *pheno = "", const char *geno = "",
                         const char *extract = "", Bool_t isXeno = kFALSE,
                         const char *xenostrain = "", const char *xenosex = "",
                         Double_t xenoage = 0, const char *xageunits = "years",
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString ArrayInfo(XArrayInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    ArrayInfo(const char *name, const char *type,
                         const char *description = "", const char *comment = "",
                         Bool_t replace = kFALSE);
      virtual TString HybridizationInfo(XHybInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    HybridizationInfo(const char *hybname, const char *type = "",
                         const char *input = "", Long_t date = 0, const char *prep = "",
                         const char *protocol = "", const char *repname = "", Int_t replica = 1,
                         const char *comment = "", Bool_t replace = kFALSE);
      virtual TString TreatmentInfo(XTreatmentInfo *info,
                         Bool_t copy = kFALSE, Bool_t replace = kFALSE);
      virtual void    TreatmentInfo(const char *name, const char *type = "",
                         Double_t conc = 0.0, const char *concunit = "",
                         Double_t time = 0.0, const char *timeunit = "",
                         const char *admin = "", const char *comment = "",
                         Bool_t replace = kFALSE);
      
      ClassDef(XProjectHandler,1) //ProjectHandler
};

#endif
