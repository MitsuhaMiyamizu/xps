// File created: 11/16/2007                          last modified: 07/04/2009
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

/******************************************************************************
* Major Revision History:
* Oct 2001 - Initial versions of xps finished
* Nov 2007 - Add classes XHandler, XProjectHandler
*
******************************************************************************/

using namespace std;

#include <Riostream.h>

#include "XPSProjectHandler.h"
#include "XPSDataTypes.h"

const Bool_t  kCS  = 0; //debug: print function names
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XHandler);
ClassImp(XProjectHandler);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHandler                                                             //
//                                                                      //
// Base class for handling datatypes for database                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XHandler::XHandler()
{
   // Default Handler constructor
   if(kCS) cout << "---XHandler::XHandler(default)------" << endl;

   fNamed = new TNamed();
   fList  = 0;
}//Constructor

//______________________________________________________________________________
XHandler::XHandler(const char *name, const char *title)
{
   // Normal Handler constructor
   // Note: fNamed is used instead of inheritance XHandler: public TNamed
   //       in order to allow multiple inheritance avoiding "dreaded diamond"
   if(kCS) cout << "---XHandler::XHandler------" << endl;

   fNamed = new TNamed(name, title);
   fList  = 0;
}//Constructor

//______________________________________________________________________________
XHandler::~XHandler()
{
   // Handler destructor
   if(kCS) cout << "---XHandler::~XHandler------" << endl;

   SafeDelete(fNamed);
   SafeDelete(fList);
}//Destructor

//______________________________________________________________________________
void XHandler::Add(TObject *obj)
{
   // Add object to list
   if(kCS) cout << "------XHandler::Add------" << endl;

   if (fList == 0) fList = new TList();

   fList->Add(obj);
}//Add


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProjectHandler                                                      //
//                                                                      //
// Class for handling datatypes for database                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProjectHandler::XProjectHandler()
                :XHandler()
{
   // Default ProjectHandler constructor
   if(kCS) cout << "---XProjectHandler::XProjectHandler(default)------" << endl;

   fHybridizations = 0;
   fTreatments     = 0;
}//Constructor

//______________________________________________________________________________
XProjectHandler::XProjectHandler(const char *name, const char *title)
                :XHandler(name, title)
{
   // Normal ProjectHandler constructor
   if(kCS) cout << "---XProjectHandler::XProjectHandler------" << endl;

   fHybridizations = 0;
   fTreatments     = 0;
}//Constructor

//______________________________________________________________________________
XProjectHandler::~XProjectHandler()
{
   // ProjectHandler destructor
   if(kCS) cout << "---XProjectHandler::~XProjectHandler------" << endl;

   fHybridizations = 0;  //do not delete since contained in fList
   fTreatments     = 0;  //do not delete since contained in fList
}//Destructor

//______________________________________________________________________________
TString XProjectHandler::LoginInfo(XLoginInfo *info, Bool_t copy, Bool_t replace)
{
   // Store login information in content folder and return userID
   // If copy is kTRUE then create new XLoginInfo before adding to content
   // If replace is kTRUE then replace existing XLoginInfo in content
   if(kCS) cout << "------XProjectHandler::LoginInfo------" << endl;

   if (info == 0) return TString(0);

   XLoginInfo *loginfo = info;
   if (copy) {
      loginfo = new XLoginInfo(*info);
   }//if

   loginfo->SetReplace(replace);
   Add(loginfo);
   return loginfo->GetUserID();
}//LoginInfo

//______________________________________________________________________________
void XProjectHandler::LoginInfo(const char *userID, const char *password,
                      Bool_t replace)
{
   // Store login information in content folder
   // If replace is kTRUE then replace existing XLoginInfo in content
   if(kCS) cout << "------XProjectHandler::LoginInfo------" << endl;

   XLoginInfo *info = new XLoginInfo(userID, password);

   info->SetReplace(replace);
   Add(info);
}//LoginInfo

//______________________________________________________________________________
TString XProjectHandler::ProjectInfo(XProjectInfo *info, Bool_t copy, Bool_t replace)
{
   // Store project information in content folder and return project name
   // If copy is kTRUE then create new XProjectInfo before adding to content
   // If replace is kTRUE then replace existing XProjectInfo in content
   if(kCS) cout << "------XProjectHandler::ProjectInfo------" << endl;

   if (info == 0) return TString(0);

   XProjectInfo *projectinfo = info;
   if (copy) {
      projectinfo = new XProjectInfo(*info);
   }//if

   projectinfo->SetReplace(replace);
   Add(projectinfo);
   return projectinfo->GetProjectName();
}//ProjectInfo

//______________________________________________________________________________
void XProjectHandler::ProjectInfo(const char *name, Long_t date, const char *type,
                      const char *description, const char *comment, Bool_t replace)
{
   // Store project information in content folder
   // If replace is kTRUE then replace existing XProjectInfo in content
   if(kCS) cout << "------XProjectHandler::ProjectInfo------" << endl;

   XProjectInfo *info = new XProjectInfo(name, date, type, description, comment);

   info->SetReplace(replace);
   Add(info);
}//ProjectInfo

//______________________________________________________________________________
TString XProjectHandler::AuthorInfo(XAuthorInfo *info, Bool_t copy, Bool_t replace)
{
   // Store author information in content folder and return author name
   // If copy is kTRUE then create new XAuthorInfo before adding to content
   // If replace is kTRUE then replace existing XAuthorInfo in content
   if(kCS) cout << "------XProjectHandler::AuthorInfo------" << endl;

   if (info == 0) return TString(0);

   XAuthorInfo *authorinfo = info;
   if (copy) {
      authorinfo = new XAuthorInfo(*info);
   }//if

   authorinfo->SetReplace(replace);
   Add(authorinfo);
   return authorinfo->GetLastName();
}//AuthorInfo

//______________________________________________________________________________
void XProjectHandler::AuthorInfo(const char *lastname, const char *firstname,
                      const char *type, const char *company, const char *department,
                      const char *mail, const char *phone, const char *comment,
                      Bool_t replace)
{
   // Store author information in content folder
   // If replace is kTRUE then replace existing XAuthorInfo in content
   if(kCS) cout << "------XProjectHandler::AuthorInfo------" << endl;

   XAuthorInfo *info = new XAuthorInfo(lastname, firstname, type, company,
                                       department, mail, phone, comment);
   info->SetReplace(replace);
   Add(info);
}//AuthorInfo

//______________________________________________________________________________
TString XProjectHandler::DatasetInfo(XDatasetInfo *info, Bool_t copy, Bool_t replace)
{
   // Store dataset information in content folder and return datset name
   // If copy is kTRUE then create new XDatasetInfo before adding to content
   // If replace is kTRUE then replace existing XDatasetInfo in content
   if(kCS) cout << "------XProjectHandler::DatasetInfo------" << endl;

   if (info == 0) return TString(0);

   XDatasetInfo *setinfo = info;
   if (copy) {
      setinfo = new XDatasetInfo(*info);
   }//if

   setinfo->SetReplace(replace);
   Add(setinfo);
   return setinfo->GetDatasetName();
}//DatasetInfo

//______________________________________________________________________________
void XProjectHandler::DatasetInfo(const char *name, const char *type,
                      const char *sample, const char *submitter, Long_t date,
                      const char *description, const char *comment, Bool_t replace)
{
   // Store dataset information in content folder
   // If replace is kTRUE then replace existing XDatasetInfo in content
   if(kCS) cout << "------XProjectHandler::DatasetInfo------" << endl;

   XDatasetInfo *setinfo = new XDatasetInfo(name, type, sample, submitter, date,
                                            description, comment);

   // set data type for data contained in dataset
//?   setinfo->SetDataType(fDataType);

   setinfo->SetReplace(replace);
   Add(setinfo);
}//DatasetInfo

//______________________________________________________________________________
TString XProjectHandler::SourceInfo(XSourceInfo *info, Bool_t copy, Bool_t replace)
{
   // Store source information in content folder and return source name
   // If copy is kTRUE then create new XSourceInfo before adding to content
   // If replace is kTRUE then replace existing XSourceInfo in content
   if(kCS) cout << "------XProjectHandler::SourceInfo------" << endl;

   if (info == 0) return TString(0);

   XSourceInfo *sourceinfo = info;
   if (copy) {
      sourceinfo = new XSourceInfo(*info);
   }//if

   sourceinfo->SetReplace(replace);
   Add(sourceinfo);
   return sourceinfo->GetSourceName();
}//SourceInfo

//______________________________________________________________________________
void XProjectHandler::SourceInfo(const char *name, const char *type, const char *species,
                      const char *subspecies, const char *description,
                      const char *comment, Bool_t replace)
{
   // Store source information in content folder
   // If replace is kTRUE then replace existing XSourceInfo in content
   if(kCS) cout << "------XProjectHandler::SourceInfo------" << endl;

   XSourceInfo *info = new XSourceInfo(name, type, species, subspecies,
                                       description, comment);

   info->SetReplace(replace);
   Add(info);
}//SourceInfo

//______________________________________________________________________________
TString XProjectHandler::SampleInfo(XSampleInfo *info, Bool_t copy, Bool_t replace)
{
   // Store sample information in content folder and return sample name
   // If copy is kTRUE then create new XSampleInfo before adding to content
   // If replace is kTRUE then replace existing XSampleInfo in content
   if(kCS) cout << "------XProjectHandler::SampleInfo------" << endl;

   if (info == 0) return TString(0);

   XSampleInfo *sampleinfo = info;
   if (copy) {
      sampleinfo = new XSampleInfo(*info);
   }//if

   sampleinfo->SetReplace(replace);
   Add(sampleinfo);
   return sampleinfo->GetSampleName();
}//SampleInfo

//______________________________________________________________________________
void XProjectHandler::SampleInfo(const char *name, const char *type, const char *sex,
                      const char *pheno, const char *geno, const char *extract,
                      Bool_t isXeno, const char *xenostrain, const char *xenosex,
                      Double_t xenoage, const char *xageunits, const char *comment,
                      Bool_t replace)
{
   // Store sample information in content folder
   // If replace is kTRUE then replace existing XSampleInfo in content
   if(kCS) cout << "------XProjectHandler::SampleInfo------" << endl;

   XSampleInfo *info = new XSampleInfo(name, type, sex, pheno, geno, extract,
                           isXeno, xenostrain, xenosex, xenoage, xageunits, comment);

   info->SetReplace(replace);
   Add(info);
}//SampleInfo

//______________________________________________________________________________
TString XProjectHandler::CellLineInfo(XCellLineInfo *info, Bool_t copy, Bool_t replace)
{
   // Store celline information in content folder and return celline name
   // If copy is kTRUE then create new XCellLineInfo before adding to content
   // If replace is kTRUE then replace existing XCellLineInfo in content
   if(kCS) cout << "------XProjectHandler::CellLineInfo------" << endl;

   if (info == 0) return TString(0);

   XCellLineInfo *cellinfo = info;
   if (copy) {
      cellinfo = new XCellLineInfo(*info);
   }//if

   cellinfo->SetReplace(replace);
   Add(cellinfo);
   return cellinfo->GetSampleName();
}//CellLineInfo

//______________________________________________________________________________
void XProjectHandler::CellLineInfo(const char *name, const char *type, const char *parent, 
                      const char *atcc, const char *mod, const char *sex, const char *pheno,
                      const char *geno, const char *extract, Bool_t isXeno,
                      const char *xenostrain, const char *xenosex, Double_t xenoage,
                      const char *xageunits, const char *comment, Bool_t replace)
{
   // Store celline information in content folder
   // If replace is kTRUE then replace existing XCellLineInfo in content
   if(kCS) cout << "------XProjectHandler::CellLineInfo------" << endl;

   XCellLineInfo *info = new XCellLineInfo(name, type, parent, atcc, mod, sex,
                             pheno, geno, extract, isXeno, xenostrain, xenosex,
                             xenoage, xageunits, comment);

   info->SetReplace(replace);
   Add(info);
}//CellLineInfo

//______________________________________________________________________________
TString XProjectHandler::PrimaryCellInfo(XPrimaryCellInfo *info, Bool_t copy, Bool_t replace)
{
   // Store cell information in content folder and return cell name
   // If copy is kTRUE then create new XPrimaryCellInfo before adding to content
   // If replace is kTRUE then replace existing XPrimaryCellInfo in content
   if(kCS) cout << "------XProjectHandler::PrimaryCellInfo------" << endl;

   if (info == 0) return TString(0);

   XPrimaryCellInfo *cellinfo = info;
   if (copy) {
      cellinfo = new XPrimaryCellInfo(*info);
   }//if

   cellinfo->SetReplace(replace);
   Add(cellinfo);
   return cellinfo->GetSampleName();
}//PrimaryCellInfo

//______________________________________________________________________________
void XProjectHandler::PrimaryCellInfo(const char *name, const char *type,  
                      Long_t date, const char *description, const char *sex,
                      const char *pheno, const char *geno, const char *extract,
                      Bool_t isXeno, const char *xenostrain, const char *xenosex,
                      Double_t xenoage, const char *xageunits, const char *comment,
                      Bool_t replace)
{
   // Store sample information in content folder
   // If replace is kTRUE then replace existing XSampleInfo in content
   if(kCS) cout << "------XProjectHandler::PrimaryCellInfo------" << endl;

   XPrimaryCellInfo *info = new XPrimaryCellInfo(name, type, date, description,
                                sex, pheno, geno, extract, isXeno, xenostrain,
                                xenosex, xenoage, xageunits, comment);

   info->SetReplace(replace);
   Add(info);
}//PrimaryCellInfo

//______________________________________________________________________________
TString XProjectHandler::TissueInfo(XTissueInfo *info, Bool_t copy, Bool_t replace)
{
   // Store tissue information in content folder and return tissue name
   // If copy is kTRUE then create new XTissueInfo before adding to content
   // If replace is kTRUE then replace existing XTissueInfo in content
   if(kCS) cout << "------XProjectHandler::TissueInfo------" << endl;

   if (info == 0) return TString(0);

   XTissueInfo *tissueinfo = info;
   if (copy) {
      tissueinfo = new XTissueInfo(*info);
   }//if

   tissueinfo->SetReplace(replace);
   Add(tissueinfo);
   return tissueinfo->GetSampleName();
}//TissueInfo

//______________________________________________________________________________
void XProjectHandler::TissueInfo(const char *name, const char *type,
                      const char *development, const char *morphology, const char *disease,
                      const char *stage, Double_t donorage, const char *ageunits,
                      const char *status, const char *sex, const char *pheno,
                      const char *geno, const char *extract, Bool_t isXeno,
                      const char *xenostrain, const char *xenosex, Double_t xenoage,
                      const char *xageunits, const char *comment, Bool_t replace)
{
   // Store tissue information in content folder
   // If replace is kTRUE then replace existing XTissueInfo in content
   if(kCS) cout << "------XProjectHandler::TissueInfo------" << endl;

   XTissueInfo *info = new XTissueInfo(name, type, development, morphology, disease,
                           stage, donorage, ageunits, status, sex, pheno, geno,
                           extract, isXeno, xenostrain, xenosex, xenoage, xageunits,
                           comment);

   info->SetReplace(replace);
   Add(info);
}//TissueInfo

//______________________________________________________________________________
TString XProjectHandler::BiopsyInfo(XBiopsyInfo *info, Bool_t copy, Bool_t replace)
{
   // Store biopsy information in content folder and return biopsy name
   // If copy is kTRUE then create new XBiopsyInfo before adding to content
   // If replace is kTRUE then replace existing XBiopsyInfo in content
   if(kCS) cout << "------XProjectHandler::BiopsyInfo------" << endl;

   if (info == 0) return TString(0);

   XBiopsyInfo *biopsyinfo = info;
   if (copy) {
      biopsyinfo = new XBiopsyInfo(*info);
   }//if

   biopsyinfo->SetReplace(replace);
   Add(biopsyinfo);
   return biopsyinfo->GetSampleName();
}//BiopsyInfo

//______________________________________________________________________________
void XProjectHandler::BiopsyInfo(const char *name, const char *type,
                      const char *morphology, const char *disease, const char *stage,
                      Double_t donorage, const char *ageunits, const char *status,
                      const char *sex, const char *pheno, const char *geno,
                      const char *extract, Bool_t isXeno, const char *xenostrain,
                      const char *xenosex, Double_t xenoage, const char *xageunits,
                      const char *comment, Bool_t replace)
{
   // Store sample information in content folder
   // If replace is kTRUE then replace existing XBiopsyInfo in content
   if(kCS) cout << "------XProjectHandler::BiopsyInfo------" << endl;

   XBiopsyInfo *info = new XBiopsyInfo(name, type, morphology, disease, stage,
                           donorage, ageunits, status, sex, pheno, geno, extract,
                           isXeno, xenostrain, xenosex, xenoage, xageunits, comment);

   info->SetReplace(replace);
   Add(info);
}//BiopsyInfo

//______________________________________________________________________________
TString XProjectHandler::ArrayInfo(XArrayInfo *info, Bool_t copy, Bool_t replace)
{
   // Store array information in content folder and return array name
   // If copy is kTRUE then create new XArrayInfo before adding to content
   // If replace is kTRUE then replace existing XArrayInfo in content
   if(kCS) cout << "------XProjectHandler::ArrayInfo------" << endl;

   if (info == 0) return TString(0);

   XArrayInfo *arrayinfo = info;
   if (copy) {
      arrayinfo = new XArrayInfo(*info);
   }//if

   arrayinfo->SetReplace(replace);
   Add(arrayinfo);
   return arrayinfo->GetArrayName();
}//ArrayInfo

//______________________________________________________________________________
void XProjectHandler::ArrayInfo(const char *name, const char *type,
                      const char *description, const char *comment, Bool_t replace)
{
   // Store array information in content folder
   // If replace is kTRUE then replace existing XArrayInfo in content
   if(kCS) cout << "------XProjectHandler::ArrayInfo------" << endl;

   XArrayInfo *info = new XArrayInfo(name, type, description, comment);

   info->SetReplace(replace);
   Add(info);
}//ArrayInfo

//______________________________________________________________________________
TString XProjectHandler::HybridizationInfo(XHybInfo *info, Bool_t copy,
                         Bool_t replace)
{
   // Store hybridization information in fHybridizations and return hybridization name
   // If copy is kTRUE then create new XHybInfo before adding to fHybridizations
   // If replace is kTRUE then replace existing XHybInfo
   if(kCS) cout << "------XProjectHandler::HybridizationInfo------" << endl;

   if (info == 0) return TString(0);

   XHybInfo *hyb = info;
   if (copy) {
      hyb = new XHybInfo(*info);
   }//if

   hyb->SetReplace(replace);
   this->AddHybridization(hyb);
   return hyb->GetHybName();
}//HybridizationInfo

//______________________________________________________________________________
void XProjectHandler::HybridizationInfo(const char *hybname, const char *type,
                      const char *input, Long_t date, const char *prep ,
                      const char *protocol, const char *repname, Int_t replica,
                      const char *comment, Bool_t replace)
{
   // Store hybridization information in list fHybridizations
   // If replace is kTRUE then replace existing XHybInfo
   if(kCS) cout << "------XProjectHandler::HybridizationInfo------" << endl;

   XHybInfo *hyb = new XHybInfo(hybname, type, input, date, prep, protocol,
                                repname, replica, comment);

   hyb->SetReplace(replace);
   this->AddHybridization(hyb);
}//HybridizationInfo

//______________________________________________________________________________
TString XProjectHandler::TreatmentInfo(XTreatmentInfo *info, Bool_t copy,
                         Bool_t replace)
{
   // Store treatment information in fTreatments and return hybridization name
   // If copy is kTRUE then create new XTreatmentInfo before adding to fTreatments
   // If replace is kTRUE then replace existing XTreatmentInfo
   if(kCS) cout << "------XProjectHandler::TreatmentInfo------" << endl;

   if (info == 0) return TString(0);

   XTreatmentInfo *treatment = info;
   if (copy) {
      treatment = new XTreatmentInfo(*info);
   }//if

   treatment->SetReplace(replace);
   this->AddTreatment(treatment);
   return treatment->GetTreatmentName();
}//TreatmentInfo

//______________________________________________________________________________
void XProjectHandler::TreatmentInfo(const char *name, const char *type,
                      Double_t conc, const char *concunit, Double_t time,
                      const char *timeunit, const char *admin, const char *comment,
                      Bool_t replace)
{
   // Store treatment information in list fTreatments
   // If replace is kTRUE then replace existing XTreatmentInfo
   if(kCS) cout << "------XProjectHandler::TreatmentInfo------" << endl;

   XTreatmentInfo *treatment = new XTreatmentInfo(name, type, conc, concunit,
                                                  time, timeunit, admin, comment);

   treatment->SetReplace(replace);
   this->AddTreatment(treatment);
}//TreatmentInfo

//______________________________________________________________________________
void XProjectHandler::AddHybridization(XHybInfo *info)
{
   // Add hybridization information to fHybridizations
   if(kCS) cout << "------XProjectHandler::AddHybridization------" << endl;

   if (fHybridizations == 0) {
      fHybridizations = new XHybridizationList(info->GetHybName(),info->GetComment());
      Add(fHybridizations);  //add to fList
   }//if

   // set replace only if info should replace an old info
   if (info->Replace() == kTRUE) {
      fHybridizations->SetReplace(info->Replace());
   }//if

   fHybridizations->Add(info);
}//AddHybridization

//______________________________________________________________________________
void XProjectHandler::AddTreatment(XTreatmentInfo *info)
{
   // Add treatment information to fTreatments
   if(kCS) cout << "------XProjectHandler::AddTreatment------" << endl;

   if (fTreatments == 0) {
      fTreatments = new XTreatmentList(info->GetTreatmentName(),"");
      Add(fTreatments);  //add to fList
   }//if

   // set replace only if info should replace an old info
   if (info->Replace() == kTRUE) {
      fTreatments->SetReplace(info->Replace());
   }//if

   fTreatments->Add(info);
}//AddTreatment
