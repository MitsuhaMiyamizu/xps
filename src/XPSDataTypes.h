// File created: 10/27/2001                          last modified: 11/23/2007
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

#ifndef __XPSDataTypes__
#define __XPSDataTypes__

#include "TNamed.h"
#include "TList.h"
#include "TString.h"
#include "TTree.h"

class XLoginInfo;


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTypeInfo                                                        //
//                                                                      //
// Base class for data types                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDataTypeInfo: public TNamed {

   protected:
      TString   fDataType;    //data type
//??      Bool_t    fHasInfo;     //has datatype info
      Bool_t    fHasInfo;     //! has datatype info
      Bool_t    fReplace;     //! replace current datatype info
      TString   fComment;     //Comment

   public:
      XDataTypeInfo();
      XDataTypeInfo(const char *name, const char *title = "",
         const char *type = "", const char *comment = "");
      XDataTypeInfo(const XDataTypeInfo &info);
      XDataTypeInfo& operator=(const XDataTypeInfo& rhs);
      virtual ~XDataTypeInfo();

      void SetDataType(const char *type) {fDataType = type;}
      void SetHasInfo(Bool_t hasInfo)    {fHasInfo  = hasInfo;}
      void SetReplace(Bool_t replace)    {fReplace  = replace;}
      void SetComment(const char *name)  {fComment = name;}

      TString GetDataType() const {return fDataType;}
      Bool_t  HasInfo()     const {return fHasInfo;}
      Bool_t  Replace()     const {return fReplace;}
      TString GetComment()  const {return fComment;}

      ClassDef(XDataTypeInfo,2) //DataTypeInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTypeList                                                        //
//                                                                      //
// Class containing list of datatype infos                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDataTypeList: public XDataTypeInfo {

   protected:
      TList       *fList;            //List of datatype infos

   public:
      XDataTypeList();
      XDataTypeList(const char *name, const char *title = "",
         const char *type = "", const char *comment = "");
      XDataTypeList(const XDataTypeList &list);
      XDataTypeList& operator=(const XDataTypeList& rhs);
      virtual ~XDataTypeList();

      void           Add(XDataTypeInfo *info);
      void           AddAt(XDataTypeInfo *info, Int_t idx);
      Int_t          Remove(const char *name);
      XDataTypeInfo *At(Int_t idx);
      XDataTypeInfo *FindDataTypeInfo(const char *name);

      TList  *GetList()             const {return fList;}
      Int_t   GetSize()             const {return fList ? fList->GetSize() : 0;}

      ClassDef(XDataTypeList,1) //DataTypeList
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDatabaseInfo                                                        //
//                                                                      //
// Class describing XPS database                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDatabaseInfo: public XDataTypeInfo {

   protected:
      TList     *fLoginList;  //List to store userID and password
      TString    fAdminID;    //UserID of administrator
      TString    fDirectory;  //Directory

   public:
      XDatabaseInfo();
      XDatabaseInfo(const char *name, const char *directory, const char *type = "",
         const char *comment = "");
      XDatabaseInfo(const XDatabaseInfo &info);
      XDatabaseInfo& operator=(const XDatabaseInfo& rhs);
      virtual ~XDatabaseInfo();

      void    AddLoginInfo(XLoginInfo *login);
      void    RemoveLoginInfo(const char *userID);
      TString GetPassword(const char *userID);
      Bool_t  IsPresentID(const char *userID);

      void    SetDBName(const char *name)   {SetTitle(name);}
      void    SetAdminID(const char *id)    {fAdminID   = id;}
      void    SetDirectory(const char *dir) {fDirectory = dir;}
      TString GetDBName()             const {return GetTitle();}
      TString GetAdminID()	           const {return fAdminID;}
      TString GetDirectory()	        const {return fDirectory;}

      ClassDef(XDatabaseInfo,1) //DatabaseInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProjectInfo                                                         //
//                                                                      //
// Class describing XPS project                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProjectInfo: public XDataTypeInfo {

   protected:
      Long_t     fDate;          //Date of creation
      TString    fDescription;   //Description
//??      XAuthorInfo *fSubmitter;    //Submitter information
//??      XLoginInfo  *fLogin;        //Login information (user ID, password)

   public:
      XProjectInfo();
      XProjectInfo(const char *name, Long_t date = 0, const char *type = "",
         const char *description = "", const char *comment = "");
      XProjectInfo(const XProjectInfo &info);
      XProjectInfo& operator=(const XProjectInfo& rhs);
      virtual ~XProjectInfo();

      void    SetProjectName(const char *name) {SetTitle(name);}
      void    SetDate(Long_t date)             {fDate        = date;}
      void    SetDescription(const char *name) {fDescription = name;}
      TString GetProjectName()           const {return GetTitle();}
      Long_t  GetDate()                  const {return fDate;}
      TString GetDescription()           const {return fDescription;}

      ClassDef(XProjectInfo,1) //ProjectInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAuthorInfo                                                          //
//                                                                      //
// Class describing author, e.g. administrator, submitter               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAuthorInfo: public XDataTypeInfo {

   protected:
      TString      fFirstName;
      TString      fCompany;
      TString      fLabName;
      TString      fMail;
      TString      fPhone;

   public:
      XAuthorInfo();
      XAuthorInfo(const char *last, const char *first = "", const char *type = "",
         const char *company = "", const char *lab = "", const char *mail = "",
         const char *phone = "", const char *comment = "");
      XAuthorInfo(const XAuthorInfo &info);
      XAuthorInfo& operator=(const XAuthorInfo& rhs);
      virtual ~XAuthorInfo();

      void SetLastName(const char *name)    {SetTitle(name);   }
      void SetFirstName(const char *name)   {fFirstName = name;}
      void SetCompany(const char *name)     {fCompany   = name;}
      void SetLabName(const char *name)     {fLabName   = name;}
      void SetMailAddress(const char *name) {fMail      = name;}
      void SetPhoneNumber(const char *name) {fPhone     = name;}
      TString GetLastName()           const {return GetTitle();}
      TString GetFirstName()          const {return fFirstName;}
      TString GetCompany()            const {return fCompany;  }
      TString GetLabName()            const {return fLabName;  }
      TString GetMailAddress()        const {return fMail;     }
      TString GetPhoneNumber()        const {return fPhone;    }

      ClassDef(XAuthorInfo,1) //AuthorInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLoginInfo                                                           //
//                                                                      //
// Class to store userID and password (and encrpyt/decrypt password)    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XLoginInfo: public XDataTypeInfo {

   protected:
      TString      fPassword;  //Password, needs to be encrypted

   public:
      XLoginInfo();
      XLoginInfo(const char *userID, const char *password = "", const char *comment = "");
      XLoginInfo(const XLoginInfo &info);
      XLoginInfo& operator=(const XLoginInfo& rhs);
      virtual ~XLoginInfo();

      void    SetPassword(const char *name);
      TString GetPassword();

//?? encrypt userID?
      void    SetUserID(const char *name) {SetTitle(name);}
      TString GetUserID()           const {return GetTitle();}

      ClassDef(XLoginInfo,1) //LoginInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDatasetInfo                                                         //
//                                                                      //
// Class describing XPS dataset                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDatasetInfo: public XDataTypeInfo {

   protected:
      TString      fSample;        //Sample type
      TString      fSubmitter;     //Submitter name
      Long_t       fDate;          //Date of creation
      TString      fDescription;   //Description

   public:
      XDatasetInfo();
      XDatasetInfo(const char *name, const char *type = "", const char *sample = "",
         const char *submitter = "", Long_t date = 0, const char *description = "",
         const char *comment = "");
      XDatasetInfo(const XDatasetInfo &info);
      XDatasetInfo& operator=(const XDatasetInfo& rhs);
      virtual ~XDatasetInfo();

      void SetDatasetName(const char *name) {SetTitle(name);}
      void SetSampleType(const char *name)  {fSample      = name;}
      void SetSubmitter(const char *name)   {fSubmitter   = name;}
      void SetDate(Long_t date)             {fDate        = date;}
      void SetDescription(const char *name) {fDescription = name;}
      TString GetDatasetName()        const {return GetTitle();}
      TString GetSampleType()         const {return fSample;}
      TString GetSubmitter()          const {return fSubmitter;}
      Long_t  GetDate()               const {return fDate;}
      TString GetDescription()        const {return fDescription;}

      ClassDef(XDatasetInfo,1) //DatasetInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSourceInfo                                                          //
//                                                                      //
// Class describing sample source                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSourceInfo: public XDataTypeInfo {

   protected:
      TString      fSpecies;       //Species
      TString      fSubSpecies;    //Subspecies or strain type
      TString      fDescription;   //Description

   public:
      XSourceInfo();
      XSourceInfo(const char *source, const char *type = "",
         const char *species = "", const char *subspecies = "",
         const char *description = "", const char *comment = "");
      XSourceInfo(const XSourceInfo &info);
      XSourceInfo& operator=(const XSourceInfo& rhs);
      virtual ~XSourceInfo();

      void SetSourceName(const char *name)  {SetTitle(name);}
      void SetSpecies(const char *name)     {fSpecies     = name;}
      void SetSubSpecies(const char *name)  {fSubSpecies  = name;}
      void SetDescription(const char *name) {fDescription = name;}
      TString GetSourceName()         const {return GetTitle();}
      TString GetSpecies()            const {return fSpecies;}
      TString GetSubSpecies()         const {return fSubSpecies;}
      TString GetDescription()        const {return fDescription;}

      ClassDef(XSourceInfo,1) //SourceInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XArrayInfo                                                           //
//                                                                      //
// Class describing microarray used for hybridization                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XArrayInfo: public XDataTypeInfo {

   protected:
      TString      fDescription;   //Description

   public:
      XArrayInfo();
      XArrayInfo(const char *name, const char *type = "",
         const char *description = "", const char *comment = "");
      XArrayInfo(const XArrayInfo &info);
      XArrayInfo& operator=(const XArrayInfo& rhs);
      virtual ~XArrayInfo();

      void SetArrayName(const char *name)   {SetTitle(name);}
      void SetDescription(const char *name) {fDescription = name;}
      TString GetArrayName()          const {return GetTitle();}
      TString GetDescription()        const {return fDescription;}

      ClassDef(XArrayInfo,1) //ArrayInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybInfo                                                             //
//                                                                      //
// Class describing hybridization parameters                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XHybInfo: public XDataTypeInfo {

   protected:
      TString      fInput;         //Full name of input data
      TString      fPreparation;   //Hybridization preparation
      TString      fProtocol;      //Hybridization protocol
      TString      fReplicaName;   //Name for group of replica
      Long_t       fDate;          //Date of hybridization
      Int_t        fReplica;       //Replica number

   public:
      XHybInfo();
      XHybInfo(const char *hybname, const char *type = "", const char *input = "",
          Long_t date = 0, const char *prep = "", const char *protocol = "",
          const char *repname = "", Int_t replica = 1, const char *comment = "");
      XHybInfo(const XHybInfo &info);
      XHybInfo& operator=(const XHybInfo& rhs);
      virtual ~XHybInfo();

      void SetHybName(const char *name)     {SetTitle(name);}
      void SetInput(const char *name)       {fInput       = name;}
      void SetHybPrep(const char *name)     {fPreparation = name;}
      void SetProtocol(const char *name)    {fProtocol    = name;}
      void SetReplicaName(const char *name) {fReplicaName = name;}
      void SetDate(Long_t date)             {fDate        = date;}
      void SetReplicaNumber(Int_t num)      {fReplica     = num;}
      TString GetHybName()            const {return GetTitle();}
      TString GetInput()              const {return fInput;}
      TString GetHybPrep()            const {return fPreparation;}
      TString GetProtocol()           const {return fProtocol;}
      TString GetReplicaName()        const {return fReplicaName;}
      Long_t  GetDate()               const {return fDate;}
      Int_t   GetReplicaNumber()      const {return fReplica;}

      ClassDef(XHybInfo,1) //HybInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybridizationList                                                   //
//                                                                      //
// Class describing sample hybridizations                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XHybridizationList: public XDataTypeList {

   protected:

   public:
      XHybridizationList();
      XHybridizationList(const char *type, const char *comment = "");
      XHybridizationList(const XHybridizationList &list);
      XHybridizationList& operator=(const XHybridizationList& rhs);
      virtual ~XHybridizationList();

      ClassDef(XHybridizationList,1) //HybridizationList
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSampleInfo                                                          //
//                                                                      //
// Class describing sample used for hybridization                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSampleInfo: public XDataTypeInfo {

   protected:
      TString      fSex;           //Sex
      TString      fPhenotype;     //Phenotype
      TString      fGenotype;      //Genotype
      TString      fExtraction;    //Extraction
      TString      fXenoStrain;    //Strain for xenograft
      TString      fXenoSex;       //Sex of acceptor strain
      TString      fAgeUnits;      //Age units
      Double_t     fXenoAge;       //Age of acceptor strain
      Bool_t       fIsXenograft;   //True if sample is xenograft

   public:
      XSampleInfo();
      XSampleInfo(const char *name, const char *type = "", const char *sex = "",
         const char *pheno = "", const char *geno = "", const char *extract = "",
         Bool_t isXeno = kFALSE, const char *xenostrain = "", const char *xenosex = "",
         Double_t xenoage = 0.0, const char *xageunits = "", const char *comment = "");
      XSampleInfo(const XSampleInfo &info);
      XSampleInfo& operator=(const XSampleInfo& rhs);
      virtual ~XSampleInfo();

      void SetSampleName(const char *name) {SetTitle(name);}
      void SetSex(const char *name)        {fSex        = name;}
      void SetPhenotype(const char *name)  {fPhenotype  = name;}
      void SetGenotype(const char *name)   {fGenotype   = name;}
      void SetExtraction(const char *name) {fExtraction = name;}
      void SetXenoStrain(const char *name) {fXenoStrain = name;}
      void SetXenoSex(const char *name)    {fXenoSex    = name;}
      void SetAgeUnits(const char *name)   {fAgeUnits   = name;}
      void SetXenoAge(Double_t age)        {fXenoAge    = age;}
      TString  GetSampleName()       const {return GetTitle();}
      TString  GetSex()              const {return fSex;}
      TString  GetPhenotype()        const {return fPhenotype;}
      TString  GetGenotype()         const {return fGenotype;}
      TString  GetExtraction()       const {return fExtraction;}
      TString  GetXenoStrain()       const {return fXenoStrain;}
      TString  GetXenoSex()          const {return fXenoSex;}
      TString  GetAgeUnits()         const {return fAgeUnits;}
      Double_t GetXenoAge()          const {return fXenoAge;}

      ClassDef(XSampleInfo,1) //SampleInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCellLineInfo                                                        //
//                                                                      //
// Class describing sample of type cell line                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XCellLineInfo: public XSampleInfo {

   protected:
      TString      fParentCell;    //Parent cell line
      TString      fATCC;          //ATCC number
      TString      fModification;  //Modification

   public:
      XCellLineInfo();
      XCellLineInfo(const char *name, const char *type, const char *parent = "", 
         const char *atcc = "", const char *mod = "", const char *sex = "",
         const char *pheno = "", const char *geno = "", const char *extract = "",
         Bool_t isXeno = kFALSE, const char *xenostrain = "", const char *xenosex = "",
         Double_t xenoage = 0, const char *xageunits = "years", const char *comment = "");
      XCellLineInfo(const XCellLineInfo &info);
      XCellLineInfo& operator=(const XCellLineInfo& rhs);
      virtual ~XCellLineInfo();

      void SetParentCell(const char *name)   {fParentCell   = name;}
      void SetATCCNumber(const char *name)   {fATCC         = name;}
      void SetModification(const char *name) {fModification = name;}
      TString GetParentCell()          const {return fParentCell;}
      TString GetATCCNumber()          const {return fATCC;}
      TString GetModification()        const {return fModification;}

      ClassDef(XCellLineInfo,1) //CellLineInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPrimaryCellInfo                                                     //
//                                                                      //
// Class describing sample of type primary cell                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XPrimaryCellInfo: public XSampleInfo {

   protected:
      Long_t       fIsolationDate; //Date of isolation of primary cell
      TString      fDescription;   //Description

   public:
      XPrimaryCellInfo();
      XPrimaryCellInfo(const char *name, const char *type, Long_t date, 
         const char *description = "", const char *sex = "",
         const char *pheno = "", const char *geno = "", const char *extract = "",
         Bool_t isXeno = kFALSE, const char *xenostrain = "", const char *xenosex = "",
         Double_t xenoage = 0, const char *xageunits = "years", const char *comment = "");
      XPrimaryCellInfo(const XPrimaryCellInfo &info);
      XPrimaryCellInfo& operator=(const XPrimaryCellInfo& rhs);
      virtual ~XPrimaryCellInfo();

      void SetIsolationDate(Long_t date)    {fIsolationDate = date;}
      void SetDescription(const char *name) {fDescription   = name;}
      Long_t  GetIsolationDate()      const {return fIsolationDate;}
      TString GetDescription()        const {return fDescription;}

      ClassDef(XPrimaryCellInfo,1) //PrimaryCellInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTissueInfo                                                          //
//                                                                      //
// Class describing sample of type tissue                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTissueInfo: public XSampleInfo {

   protected:
      TString      fDevelopment;   //Development stage
      TString      fMorphology;    //Morphology
      TString      fDisease;       //Disease
      TString      fDiseaseStage;  //Disease stage
      Double_t     fDonorAge;      //Donor age
      TString      fAgeUnits;      //Age units
      TString      fDonorStatus;   //Donor status

   public:
      XTissueInfo();
      XTissueInfo(const char *name, const char *type, const char *development = "", 
         const char *morphology = "", const char *disease = "", const char *stage = "",
         Double_t donorage = 0, const char *ageunits = "", const char *status = "",
         const char *sex = "", const char *pheno = "", const char *geno = "",
         const char *extract = "", Bool_t isXeno = kFALSE, const char *xenostrain = "",
         const char *xenosex = "", Double_t xenoage = 0, const char *xageunits = "years",
         const char *comment = "");
      XTissueInfo(const XTissueInfo &info);
      XTissueInfo& operator=(const XTissueInfo& rhs);
      virtual ~XTissueInfo();

      void SetDevelopmentStage(const char *name) {fDevelopment  = name;}
      void SetMorphology(const char *name)       {fMorphology   = name;}
      void SetDisease(const char *name)          {fDisease      = name;}
      void SetDiseaseStage(const char *name)     {fDiseaseStage = name;}
      void SetDonorAge(Double_t age)             {fDonorAge     = age;}
      void SetAgeUnits(const char *name)         {fAgeUnits     = name;}
      void SetDonorStatus(const char *name)      {fDonorStatus  = name;}
      TString  GetDevelopmentStage()       const {return fDevelopment;}
      TString  GetMorphology()             const {return fMorphology;}
      TString  GetDisease()                const {return fDisease;}
      TString  GetDiseaseStage()           const {return fDiseaseStage;}
      Double_t GetDonorAge()               const {return fDonorAge;}
      TString  GetAgeUnits()               const {return fAgeUnits;}
      TString  GetDonorStatus()            const {return fDonorStatus;}

      ClassDef(XTissueInfo,1) //TissueInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBiopsyInfo                                                          //
//                                                                      //
// Class describing sample of type biopsy                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XBiopsyInfo: public XSampleInfo {

   protected:
      TString      fMorphology;    //Morphology
      TString      fDisease;       //Disease
      TString      fDiseaseStage;  //Disease stage
      Double_t     fDonorAge;      //Donor age
      TString      fAgeUnits;      //Age units
      TString      fDonorStatus;   //Donor status

   public:
      XBiopsyInfo();
      XBiopsyInfo(const char *name, const char *type, const char *morphology = "",
         const char *disease = "", const char *stage = "", Double_t donorage = 0,
         const char *ageunits = "", const char *status = "", const char *sex = "",
         const char *pheno = "", const char *geno = "", const char *extract = "",
         Bool_t isXeno = kFALSE, const char *xenostrain = "", const char *xenosex = "",
         Double_t xenoage = 0, const char *xageunits = "years", const char *comment = "");
      XBiopsyInfo(const XBiopsyInfo &info);
      XBiopsyInfo& operator=(const XBiopsyInfo& rhs);
      virtual ~XBiopsyInfo();

      void SetMorphology(const char *name)   {fMorphology   = name;}
      void SetDisease(const char *name)      {fDisease      = name;}
      void SetDiseaseStage(const char *name) {fDiseaseStage = name;}
      void SetDonorAge(Double_t age)         {fDonorAge     = age;}
      void SetAgeUnits(const char *name)     {fAgeUnits     = name;}
      void SetDonorStatus(const char *name)  {fDonorStatus  = name;}
      TString  GetMorphology()         const {return fMorphology;}
      TString  GetDisease()            const {return fDisease;}
      TString  GetDiseaseStage()       const {return fDiseaseStage;}
      Double_t GetDonorAge()           const {return fDonorAge;}
      TString  GetAgeUnits()           const {return fAgeUnits;}
      TString  GetDonorStatus()        const {return fDonorStatus;}

      ClassDef(XBiopsyInfo,1) //BiopsyInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreatmentInfo                                                       //
//                                                                      //
// Class describing a single treatment                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTreatmentInfo: public XDataTypeInfo {

   protected:
      Double_t     fConcentration;  //Concentration
      TString      fConcUnit;       //Concentration units
      Double_t     fTime;           //Time
      TString      fTimeUnit;       //Time units
      TString      fAdministration; //Administration

   public:
      XTreatmentInfo();
      XTreatmentInfo(const char *name, const char *type = "", Double_t conc = 0.0,
         const char *concunit = "", Double_t time = 0.0, const char *timeunit = "",
         const char *admin = "", const char *comment = "");
      XTreatmentInfo(const XTreatmentInfo &info);
      XTreatmentInfo& operator=(const XTreatmentInfo& rhs);
      virtual ~XTreatmentInfo();

      void SetTreatmentName(const char *name)  {SetTitle(name);}
      void SetConcentration(Double_t conc)     {fConcentration  = conc;}
      void SetConcUnit(const char *name)       {fConcUnit       = name;}
      void SetTime(Double_t time)              {fTime           = time;}
      void SetTimeUnit(const char *name)       {fTimeUnit       = name;}
      void SetAdministration(const char *name) {fAdministration = name;}
      TString  GetTreatmentName()        const {return GetTitle();}
      Double_t GetConcentration()        const {return fConcentration;}
      TString  GetConcUnit()             const {return fConcUnit;}
      Double_t GetTime()                 const {return fTime;}
      TString  GetTimeUnit()             const {return fTimeUnit;}
      TString  GetAdministration()       const {return fAdministration;}

      ClassDef(XTreatmentInfo,1) //TreatmentInfo
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreatmentList                                                       //
//                                                                      //
// Class describing sample treatments                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTreatmentList: public XDataTypeList {

   protected:

   public:
      XTreatmentList();
      XTreatmentList(const char *type, const char *comment = "");
      XTreatmentList(const XTreatmentList &list);
      XTreatmentList& operator=(const XTreatmentList& rhs);
      virtual ~XTreatmentList();

      ClassDef(XTreatmentList,1) //TreatmentList
};

#endif

