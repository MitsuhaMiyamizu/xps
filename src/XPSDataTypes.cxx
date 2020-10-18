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


/******************************************************************************
* Major Revision History:
* Oct 2001 - Initial versions finished
* Jun 2005 - Derive classes from base class XDataTypeInfo.
* Oct 2005 - Add member fDataType to base class XDataTypeInfo.
* Nov 2007 - Add member fReplace to base class XDataTypeInfo and set ClassDef=2.
*          - Add classes XDataTypeList, XHybridizationList, XTreatmentList.
*
******************************************************************************/

using namespace std;

#include <Riostream.h>
#include "TError.h"

#include "XPSDataTypes.h"

const Bool_t kCS  = 0; //debug: print function names
const Bool_t kCSa = 0; //debug: print function names in loops

ClassImp(XDataTypeInfo);
ClassImp(XDataTypeList);
ClassImp(XDatabaseInfo);
ClassImp(XProjectInfo);
ClassImp(XAuthorInfo);
ClassImp(XLoginInfo);
ClassImp(XDatasetInfo);
ClassImp(XSourceInfo);
ClassImp(XArrayInfo);
ClassImp(XHybInfo);
ClassImp(XHybridizationList);
ClassImp(XSampleInfo);
ClassImp(XCellLineInfo);
ClassImp(XPrimaryCellInfo);
ClassImp(XTissueInfo);
ClassImp(XBiopsyInfo);
ClassImp(XTreatmentInfo);
ClassImp(XTreatmentList);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTypeInfo                                                        //
//                                                                      //
// Base class for data types                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDataTypeInfo::XDataTypeInfo()
              :TNamed()
{
   // Default datatype information constructor
   if(kCS) cout << "---XDataTypeInfo::XDataTypeInfo(default)------" << endl;

   fDataType = "";
   fHasInfo  = kFALSE;
   fReplace  = kFALSE;
   fComment  = "";
}//Constructor

//______________________________________________________________________________
XDataTypeInfo::XDataTypeInfo(const char *name, const char *title, const char *type,
               const char *comment)
              :TNamed(name, title)
{
   // Datatype information constructor
   if(kCS) cout << "---XDataTypeInfo::XDataTypeInfo------" << endl;

   fDataType = type;
   fComment  = comment;
   fHasInfo  = kFALSE;
   fReplace  = kFALSE;
}//Constructor

//______________________________________________________________________________
XDataTypeInfo::XDataTypeInfo(const XDataTypeInfo &info) 
              :TNamed(info),fHasInfo(info.fHasInfo)
{
   // Datatype information copy constructor
   if(kCS) cout << "---XDataTypeInfo::XDataTypeInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XDataTypeInfo& XDataTypeInfo::operator=(const XDataTypeInfo& rhs)
{
   // Datatype information assignment operator.
   if(kCS) cout << "---XDataTypeInfo::operator=------" << endl;

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fDataType = rhs.fDataType;
      fComment  = rhs.fComment;
      fHasInfo  = rhs.fHasInfo;
      fReplace  = rhs.fReplace;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDataTypeInfo::~XDataTypeInfo()
{
   // Datatype information destructor
   if(kCS) cout << "---XDataTypeInfo::~XDataTypeInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDataTypeList                                                        //
//                                                                      //
// Class containing list of datatype infos                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDataTypeList::XDataTypeList()
              :XDataTypeInfo()
{
   // Default datatype information constructor
   if(kCS) cout << "---XDataTypeList::XDataTypeList(default)------" << endl;

   fList = new TList();
}//Constructor

//______________________________________________________________________________
XDataTypeList::XDataTypeList(const char *name, const char *title, const char *type,
               const char *comment)
              :XDataTypeInfo(name, title, type, comment)
{
   // Treatment datatype information constructor
   if(kCS) cout << "---XDataTypeList::XDataTypeList------" << endl;

   fList = new TList();
}//Constructor

//______________________________________________________________________________
XDataTypeList::XDataTypeList(const XDataTypeList &list) 
              :XDataTypeInfo(list)
{
   // Datatype information copy constructor
   if(kCS) cout << "---XDataTypeList::XDataTypeList(copy)------" << endl;

//to do: test if correct?
   fList = 0;
   if (list.fList != 0) {
      fList = new TList();
      for (Int_t i=0; i<(list.fList)->GetSize(); i++) {
         fList->AddAt((list.fList)->At(i), i);
      }//for_i
   }//if
}//CopyConstructor

//______________________________________________________________________________
XDataTypeList& XDataTypeList::operator=(const XDataTypeList& rhs)
{
   // Datatype information assignment operator.
   if(kCS) cout << "---XDataTypeList::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fList = 0;
      if (rhs.fList != 0) {
         fList = new TList();
         for (Int_t i=0; i<(rhs.fList)->GetSize(); i++) {
            fList->AddAt((rhs.fList)->At(i), i);
         }//for_i
      }//if
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDataTypeList::~XDataTypeList()
{
   // Datatype information destructor
   if(kCS) cout << "---XDataTypeList::~XDataTypeList------" << endl;

   if(fList) {fList->Delete(); delete fList; fList = 0;}
}//Destructor

//______________________________________________________________________________
void XDataTypeList::Add(XDataTypeInfo *info)
{
   // Add datatype info to datatype list
   if(kCSa) cout << "------XDataTypeList::Add------" << endl;

   fList->Add(info);
   fHasInfo = kTRUE;
}//Add

//______________________________________________________________________________
void XDataTypeList::AddAt(XDataTypeInfo *info, Int_t idx)
{
   // Add datatype info to datatype list
   if(kCSa) cout << "------XDataTypeList::AddAt------" << endl;

   fList->AddAt(info, idx);
   fHasInfo = kTRUE;
}//AddAt

//______________________________________________________________________________
Int_t XDataTypeList::Remove(const char *name)
{
   // Remove datatype info name from datatype list 
   if(kCS) cout << "------XDataTypeList::Remove------" << endl;

   Int_t size = fList->GetSize();
   if (size == 0) {fHasInfo = kFALSE; return 0;}

// Loop over treatment list
   TIter next(fList);
   XDataTypeInfo *info = 0;
   while ((info = (XDataTypeInfo*)next())) {
      TString infoname = info->GetName();
      if (strcmp(name,  infoname.Data())  == 0) {
         fList->Remove(info);
         size--;
      }//if
   }//while

   fHasInfo = (size > 0) ? kTRUE: kFALSE;
   return size;
}//Remove

//______________________________________________________________________________
XDataTypeInfo *XDataTypeList::At(Int_t idx)
{
   // Return datatype info at idx of list
   if(kCSa) cout << "------XDataTypeList::At------" << endl;

   return (XDataTypeInfo*)(fList->At(idx));
}//AddAt

//______________________________________________________________________________
XDataTypeInfo *XDataTypeList::FindDataTypeInfo(const char *name)
{
   // Find datatype info with name stored in list fList
   if(kCS) cout << "------XDataTypeList::FindDataTypeInfo------" << endl;

   return (XDataTypeInfo*)fList->FindObject(name);
}//FindDataTypeInfo


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDatabaseInfo                                                        //
//                                                                      //
// Class describing XPS database                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDatabaseInfo::XDatabaseInfo()
              :XDataTypeInfo()
{
   // Default database information constructor
   if(kCS) cout << "---XDatabaseInfo::XDatabaseInfo(default)------" << endl;

   fDirectory = "";
   fAdminID   = "";
   fLoginList = new TList();
}//Constructor

//______________________________________________________________________________
XDatabaseInfo::XDatabaseInfo(const char *name, const char *directory, const char *type,
               const char *comment)
              :XDataTypeInfo("Database", name, type, comment)
{
   // Database information constructor
   if(kCS) cout << "---XDatabaseInfo::XDatabaseInfo------" << endl;

   fDirectory = directory;
   fAdminID   = "";
   fLoginList = new TList();
   fHasInfo   = kTRUE;
}//Constructor

//______________________________________________________________________________
XDatabaseInfo::XDatabaseInfo(const XDatabaseInfo &info) 
              :XDataTypeInfo(info),fAdminID(info.fAdminID),fDirectory(info.fDirectory)
{
   // Database information copy constructor
   if(kCS) cout << "---XDatabaseInfo::XDatabaseInfo(copy)------" << endl;

//to do: test if correct?
   fLoginList = 0;
   if (info.fLoginList != 0) {
      fLoginList = new TList();
      for (Int_t i=0; i<(info.fLoginList)->GetSize(); i++) {
         fLoginList->AddAt((info.fLoginList)->At(i), i);
      }//for_i
   }//if
}//CopyConstructor

//______________________________________________________________________________
XDatabaseInfo& XDatabaseInfo::operator=(const XDatabaseInfo& rhs)
{
   // Database information assignment operator.
   if(kCS) cout << "---XDatabaseInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fAdminID   = rhs.fAdminID;
      fDirectory = rhs.fDirectory;
//to do: need to copy fLoginList!!!
   cout << "Error: Copy of fLoginList not yet implemented." << endl;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDatabaseInfo::~XDatabaseInfo()
{
   // Database information destructor
   if(kCS) cout << "---XDatabaseInfo::~XDatabaseInfo" << endl;

   if (fLoginList) {fLoginList->Delete(); delete fLoginList; fLoginList = 0;}
}//Destructor

//______________________________________________________________________________
void XDatabaseInfo::AddLoginInfo(XLoginInfo *login)
{
   // Add Login information to Login list
   if(kCS) cout << "------XDatabaseInfo::AddLoginInfo------" << endl;

   fLoginList->Add(login); 
}//AddLoginInfo

//______________________________________________________________________________
void XDatabaseInfo::RemoveLoginInfo(const char *userID)
{
   // Remove Login information to Login list
   if(kCS) cout << "------XDatabaseInfo::RemoveLoginInfo------" << endl;

   XLoginInfo *login = 0;
   login = (XLoginInfo*)fLoginList->FindObject(userID); 
   fLoginList->Remove(login); 
}//RemoveLoginInfo

//______________________________________________________________________________
TString XDatabaseInfo::GetPassword(const char *userID)
{
   // Return password for userID
   if(kCS) cout << "------XDatabaseInfo::GetPassword------" << endl;

   TString password = "";

   TIter next(fLoginList);
   XLoginInfo *login = 0;
   while ((login = (XLoginInfo*)next())) {
      if (strcmp(login->GetName(),userID) == 0) {
         password = login->GetPassword();
         break;
      }//if
   }//while

   return password;
}//GetPassword

//______________________________________________________________________________
Bool_t XDatabaseInfo::IsPresentID(const char *userID)
{
   // Return true if userID is in login list
   if(kCS) cout << "------XDatabaseInfo::IsPresentID------" << endl;

   Bool_t isPresent = kFALSE;

   TIter next(fLoginList);
   XLoginInfo *login = 0;
   while ((login = (XLoginInfo*)next())) {
      if (strcmp(login->GetName(),userID) == 0) {
         isPresent = kTRUE;
         break;
      }//if
   }//while

   return isPresent;
}//IsPresentID


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProjectInfo                                                         //
//                                                                      //
// Class describing XPS project                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProjectInfo::XProjectInfo()
             :XDataTypeInfo()
{
   // Default Project information constructor
   if(kCS) cout << "---XProjectInfo::XProjectInfo(default)------" << endl;

   fDate        = 0;
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XProjectInfo::XProjectInfo(const char *name, Long_t date, const char *type,
              const char *description, const char *comment)
             :XDataTypeInfo("Project", name, type, comment)
{
   // Project information constructor
   if(kCS) cout << "---XProjectInfo::XProjectInfo------" << endl;

   fDate        = date;
   fDescription = description;
//?? check if necessary infos name and date are set!
//also in other datatype infos!!!
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XProjectInfo::XProjectInfo(const XProjectInfo &info) 
             :XDataTypeInfo(info),fDescription(info.fDescription)
{
   // Project information copy constructor
   if(kCS) cout << "---XProjectInfo::XProjectInfo(copy)------" << endl;

   fDate = info.fDate;
}//CopyConstructor

//______________________________________________________________________________
XProjectInfo& XProjectInfo::operator=(const XProjectInfo& rhs)
{
   // Project information assignment operator.
   if(kCS) cout << "---XProjectInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fDate        = rhs.fDate;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XProjectInfo::~XProjectInfo()
{
   // Project information destructor
   if(kCS) cout << "---XProjectInfo::~XProjectInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAuthorInfo                                                          //
//                                                                      //
// Class describing author, e.g. administrator, submitter               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAuthorInfo::XAuthorInfo()
            :XDataTypeInfo()
{
   // Default author information constructor
   if(kCS) cout << "---XAuthorInfo::XAuthorInfo(default)------" << endl;

   fFirstName = "";
   fCompany   = "";
   fLabName   = "";
   fMail      = "";
   fPhone     = "";
}//Constructor

//______________________________________________________________________________
XAuthorInfo::XAuthorInfo(const char *last, const char *first, const char *type,
             const char *company, const char *lab, const char *mail,
             const char *phone, const char *comment)
            :XDataTypeInfo("Author", last, type, comment)
{
   // Author information constructor
   if(kCS) cout << "---XAuthorInfo::XAuthorInfo------" << endl;

   fFirstName = first;
   fCompany   = company;
   fLabName   = lab;
   fMail      = mail;
   fPhone     = phone;
   fHasInfo   = kTRUE;
}//Constructor

//______________________________________________________________________________
XAuthorInfo::XAuthorInfo(const XAuthorInfo &info) 
            :XDataTypeInfo(info),fFirstName(info.fFirstName),fCompany(info.fCompany),
             fLabName(info.fLabName),fMail(info.fMail),fPhone(info.fPhone)
{
   // Author information copy constructor
   if(kCS) cout << "---XAuthorInfo::XAuthorInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XAuthorInfo& XAuthorInfo::operator=(const XAuthorInfo& rhs)
{
   // Author information assignment operator.
   if(kCS) cout << "---XAuthorInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fFirstName = rhs.fFirstName;
      fCompany   = rhs.fCompany;
      fLabName   = rhs.fLabName;
      fMail      = rhs.fMail;
      fPhone     = rhs.fPhone;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XAuthorInfo::~XAuthorInfo()
{
   // Author information destructor
   if(kCS) cout << "---XAuthorInfo::~XAuthorInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLoginInfo                                                           //
//                                                                      //
// Class to store userID and password (and encrpyt/decrypt password)    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XLoginInfo::XLoginInfo()
           :XDataTypeInfo()
{
   //Default Login information constructor
   if(kCS) cout << "---XLoginInfo::XLoginInfo(default)------" << endl;

   fPassword = "";
}//Constructor

//______________________________________________________________________________
XLoginInfo::XLoginInfo(const char *userID, const char *password, const char *comment)
           :XDataTypeInfo("Login", userID, "", comment)
{
   // Login information constructor
   if(kCS) cout << "---XLoginInfo::XLoginInfo------" << endl;

   SetName(userID);
   this->SetPassword(password);

   fHasInfo = kTRUE;
}//Constructor

//______________________________________________________________________________
XLoginInfo::XLoginInfo(const XLoginInfo &info) 
           :XDataTypeInfo(info),fPassword(info.fPassword)
{
   // Login information copy constructor
   if(kCS) cout << "---XLoginInfo::XLoginInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XLoginInfo& XLoginInfo::operator=(const XLoginInfo& rhs)
{
   // Login information assignment operator.
   if(kCS) cout << "---XLoginInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fPassword = rhs.fPassword;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XLoginInfo::~XLoginInfo()
{
   // Login information destructor
   if(kCS) cout << "---XLoginInfo::~XLoginInfo------" << endl;

   fPassword = "";
}//Destructor

//______________________________________________________________________________
void XLoginInfo::SetPassword(const char *name)
{
   // Encrypt and store encrypted password
   if(kCS) cout << "------XLoginInfo::SetPassword------" << endl;

   TString password = name; 

// encrypt password
   if (password != "") {
      for (Int_t i = 0; i < password.Length(); i++) {
         char inv = ~password(i);
         password.Replace(i, 1, inv);
      }//for
   }//if

   fPassword = password; 
}//SetPassword

//______________________________________________________________________________
TString XLoginInfo::GetPassword()
{
   // Decrypt and return decrypted password 
   if(kCS) cout << "------XLoginInfo::GetPassword------" << endl;

   TString password = fPassword; 

// decrypt password
   if (password != "") {
      for (Int_t i = 0; i < password.Length(); i++) {
         char inv = ~password(i);
         password.Replace(i, 1, inv);
      }//for
   }//if

   return password;
}//GetPassword


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDatasetInfo                                                         //
//                                                                      //
// Class describing XPS dataset                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDatasetInfo::XDatasetInfo()
             :XDataTypeInfo()
{
   // Default dataset information constructor
   if(kCS) cout << "---XDatasetInfo::XDatasetInfo(default)------" << endl;

   fSample      = "";
   fSubmitter   = "";
   fDate        = 0;
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XDatasetInfo::XDatasetInfo(const char *name, const char *type, const char *sample,
              const char *submitter, Long_t date, const char *description,
              const char *comment)
             :XDataTypeInfo("Dataset", name, type, comment)
{
   // Dataset information constructor
   // type: dataset type: UD, TS, DR, MC
   if(kCS) cout << "---XDatasetInfo::XDatasetInfo------" << endl;

   fSample      = sample;
   fSubmitter   = submitter;
   fDate        = date;
   fDescription = description;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XDatasetInfo::XDatasetInfo(const XDatasetInfo &info) 
             :XDataTypeInfo(info), fSample(info.fSample), fSubmitter(info.fSubmitter),
              fDescription(info.fDescription)
{
   // Dataset information copy constructor
   if(kCS) cout << "---XDatasetInfo::XDatasetInfo(copy)------" << endl;

   fDate = info.fDate;
}//CopyConstructor

//______________________________________________________________________________
XDatasetInfo& XDatasetInfo::operator=(const XDatasetInfo& rhs)
{
   // Dataset information assignment operator.
   if(kCS) cout << "---XDatasetInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fSample      = rhs.fSample;
      fSubmitter   = rhs.fSubmitter;
      fDate        = rhs.fDate;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XDatasetInfo::~XDatasetInfo()
{
   // Dataset information destructor
   if(kCS) cout << "---XDatasetInfo::~XDatasetInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSourceInfo                                                          //
//                                                                      //
// Class describing sample source                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSourceInfo::XSourceInfo()
            :XDataTypeInfo()
{
   // Default source information constructor
   if(kCS) cout << "---XSourceInfo::XSourceInfo(default)------" << endl;

   fSpecies     = "";
   fSubSpecies  = "";
   fDescription = "";
}//Constructor

//______________________________________________________________________________
XSourceInfo::XSourceInfo(const char *source, const char *type, const char *species,
             const char *subspecies, const char *description, const char *comment)
            :XDataTypeInfo("Source", source, type, comment)
{
   // Source information constructor
   if(kCS) cout << "---XSourceInfo::XSourceInfo------" << endl;

   fSpecies     = species;
   fSubSpecies  = subspecies;
   fDescription = description;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XSourceInfo::XSourceInfo(const XSourceInfo &info) 
            :XDataTypeInfo(info),fSpecies(info.fSpecies),fSubSpecies(info.fSubSpecies),
             fDescription(info.fDescription)
{
   // Source information copy constructor
   if(kCS) cout << "---XSourceInfo::XSourceInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XSourceInfo& XSourceInfo::operator=(const XSourceInfo& rhs)
{
   // Source information assignment operator.
   if(kCS) cout << "---XSourceInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fSpecies     = rhs.fSpecies;
      fSubSpecies  = rhs.fSubSpecies;
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XSourceInfo::~XSourceInfo()
{
   // Source information destructor
   if(kCS) cout << "---XSourceInfo::~XSourceInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XArrayInfo                                                           //
//                                                                      //
// Class describing microarray used for hybridization                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XArrayInfo::XArrayInfo()
           :XDataTypeInfo()
{
   // Default array information constructor
   if(kCS) cout << "---XArrayInfo::XArrayInfo(default)------" << endl;

   fDescription = "";
}//Constructor

//______________________________________________________________________________
XArrayInfo::XArrayInfo(const char *name, const char *type, const char *description,
            const char *comment)
           :XDataTypeInfo("Array", name, type, comment)
{
   // Array information constructor
   if(kCS) cout << "---XArrayInfo::XArrayInfo------" << endl;

   fDescription = description;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XArrayInfo::XArrayInfo(const XArrayInfo &info) 
           :XDataTypeInfo(info), fDescription(info.fDescription)
{
   // Array information copy constructor
   if(kCS) cout << "---XArrayInfo::XArrayInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XArrayInfo& XArrayInfo::operator=(const XArrayInfo& rhs)
{
   // Array information assignment operator.
   if(kCS) cout << "---XArrayInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fDescription = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XArrayInfo::~XArrayInfo()
{
   // Array information destructor
   if(kCS) cout << "---XArrayInfo::~XArrayInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybInfo                                                             //
//                                                                      //
// Class describing hybridization parameters for one sample             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XHybInfo::XHybInfo()
         :XDataTypeInfo()
{
   // Default hybridization information constructor
   if(kCS) cout << "---XHybInfo::XHybInfo(default)------" << endl;

   fInput       = "";
   fPreparation = "";
   fProtocol    = "";
   fReplicaName = "";
   fDate        = 0;
   fReplica     = 0;
}//Constructor

//______________________________________________________________________________
XHybInfo::XHybInfo(const char *hybname, const char *type, const char *input,
          Long_t date, const char *prep, const char *protocol, const char *repname,
          Int_t replica, const char *comment)
         :XDataTypeInfo("Hybridization", hybname, type, comment)
{
   // Hybridization information constructor
   if(kCS) cout << "---XHybInfo::XHybInfo------" << endl;

//?? in title full input?? then fInput not needed!!
//   TString name = Path2Name(input,"/",".");
//   SetTitle(name);
   fInput       = input;
   fPreparation = prep;
   fProtocol    = protocol;
   fReplicaName = repname;
   fDate        = date;
   fReplica     = replica;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XHybInfo::XHybInfo(const XHybInfo &info) 
         :XDataTypeInfo(info),fInput(info.fInput),fPreparation(info.fPreparation),
          fProtocol(info.fProtocol),fReplicaName(info.fReplicaName)
{
   // Hybridization information copy constructor
   if(kCS) cout << "---XHybInfo::XHybInfo(copy)------" << endl;

   fDate    = info.fDate;
   fReplica = info.fReplica;
}//CopyConstructor

//______________________________________________________________________________
XHybInfo& XHybInfo::operator=(const XHybInfo& rhs)
{
   // Hybridization information assignment operator.
   if(kCS) cout << "---XHybInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fInput       = rhs.fInput;
      fPreparation = rhs.fPreparation;
      fProtocol    = rhs.fProtocol;
      fReplicaName = rhs.fReplicaName;
      fDate        = rhs.fDate;
      fReplica     = rhs.fReplica;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XHybInfo::~XHybInfo()
{
   // Hybridization information destructor
   if(kCS) cout << "---XHybInfo::~XHybInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XHybridizationList                                                   //
//                                                                      //
// Class describing sample hybridizations                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XHybridizationList::XHybridizationList()
                   :XDataTypeList()
{
   // Default hybridization list constructor
   if(kCS) cout << "---XHybridizationList::XHybridizationList(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XHybridizationList::XHybridizationList(const char *type, const char *comment)
                   :XDataTypeList("Hybridizations", "HybridizationList", type, comment)
{
   // Treatment hybridization list constructor
   if(kCS) cout << "---XHybridizationList::XHybridizationList------" << endl;

}//Constructor

//______________________________________________________________________________
XHybridizationList::XHybridizationList(const XHybridizationList &list) 
                   :XDataTypeList(list)
{
   // Hybridization list copy constructor
   if(kCS) cout << "---XHybridizationList::XHybridizationList(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XHybridizationList& XHybridizationList::operator=(const XHybridizationList& rhs)
{
   // Hybridization list assignment operator.
   if(kCS) cout << "---XHybridizationList::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeList::operator=(rhs);
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XHybridizationList::~XHybridizationList()
{
   // Hybridization information destructor
   if(kCS) cout << "---XHybridizationList::~XHybridizationList------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSampleInfo                                                          //
//                                                                      //
// Class describing sample used for hybridization                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSampleInfo::XSampleInfo()
            :XDataTypeInfo()
{
   // Default sample information constructor
   if(kCS) cout << "---XSampleInfo::XSampleInfo(default)------" << endl;

   fSex         = "";
   fPhenotype   = "";
   fGenotype    = "";
   fExtraction  = "";
   fXenoStrain  = "";
   fXenoSex     = "";
   fAgeUnits    = "";
   fXenoAge     = 0;
   fIsXenograft = kFALSE;
}//Constructor

//______________________________________________________________________________
XSampleInfo::XSampleInfo(const char *name, const char *type, const char *sex,
             const char *pheno, const char *geno, const char *extract, Bool_t isXeno,
             const char *xenostrain, const char *xenosex, Double_t xenoage,
             const char *ageunits, const char *comment)
            :XDataTypeInfo("Sample", name, type, comment)
{
   // Sample information constructor
   if(kCS) cout << "---XSampleInfo::XSampleInfo------" << endl;

   fSex         = sex;
   fPhenotype   = pheno;
   fGenotype    = geno;
   fExtraction  = extract;
   fXenoStrain  = xenostrain;
   fXenoSex     = xenosex;
   fAgeUnits    = ageunits;
   fXenoAge     = xenoage;
   fIsXenograft = isXeno;
   fHasInfo     = kTRUE;
}//Constructor

//______________________________________________________________________________
XSampleInfo::XSampleInfo(const XSampleInfo &info) 
            :XDataTypeInfo(info), fSex(info.fSex), fPhenotype(info.fPhenotype),
             fGenotype(info.fGenotype), fExtraction(info.fExtraction),
             fXenoStrain(info.fXenoStrain), fXenoSex(info.fXenoSex),
             fAgeUnits(info.fAgeUnits)
{
   // Sample information copy constructor
   if(kCS) cout << "---XSampleInfo::XSampleInfo(copy)------" << endl;

   fXenoAge     = info.fXenoAge;
   fIsXenograft = info.fIsXenograft;
}//CopyConstructor

//______________________________________________________________________________
XSampleInfo& XSampleInfo::operator=(const XSampleInfo& rhs)
{
   // Sample information assignment operator.
   if(kCS) cout << "---XSampleInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fSex         = rhs.fSex;
      fPhenotype   = rhs.fPhenotype;
      fGenotype    = rhs.fGenotype;
      fExtraction  = rhs.fExtraction;
      fXenoStrain  = rhs.fXenoStrain;
      fXenoSex     = rhs.fXenoSex;
      fAgeUnits    = rhs.fAgeUnits;
      fXenoAge     = rhs.fXenoAge;
      fIsXenograft = rhs.fIsXenograft;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XSampleInfo::~XSampleInfo()
{
   // Sample information destructor
   if(kCS) cout << "---XSampleInfo::~XSampleInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XCellLineInfo                                                        //
//                                                                      //
// Class describing sample of type cell line                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XCellLineInfo::XCellLineInfo()
              :XSampleInfo()
{
   // Default cell line information constructor
   if(kCS) cout << "---XCellLineInfo::XCellLineInfo(default)------" << endl;

   fParentCell   = "";
   fATCC         = "";
   fModification = "";
}//Constructor

//______________________________________________________________________________
XCellLineInfo::XCellLineInfo(const char *name, const char *type, const char *parent, 
               const char *atcc, const char *mod, const char *sex, const char *pheno,
               const char *geno, const char *extract, Bool_t isXeno,
               const char *xenostrain, const char *xenosex, Double_t xenoage,
               const char *ageunits, const char *comment)
              :XSampleInfo(name, type, sex, pheno, geno, extract, isXeno,
               xenostrain, xenosex, xenoage, ageunits, comment)
{
   // Cell line information constructor
   if(kCS) cout << "---XCellLineInfo::XCellLineInfo------" << endl;

   fParentCell   = parent;
   fATCC         = atcc;
   fModification = mod;
}//Constructor

//______________________________________________________________________________
XCellLineInfo::XCellLineInfo(const XCellLineInfo &info) 
              :XSampleInfo(info), fParentCell(info.fParentCell),
               fATCC(info.fATCC), fModification(info.fModification)
{
   // Cell line information copy constructor
   if(kCS) cout << "---XCellLineInfo::XCellLineInfo(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XCellLineInfo& XCellLineInfo::operator=(const XCellLineInfo& rhs)
{
   // Cell line information assignment operator.
   if(kCS) cout << "---XProjectInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fParentCell   = rhs.fParentCell;
      fATCC         = rhs.fATCC;
      fModification = rhs.fModification;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XCellLineInfo::~XCellLineInfo()
{
   // Cell line information destructor
   if(kCS) cout << "---XCellLineInfo::~XCellLineInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPrimaryCellInfo                                                     //
//                                                                      //
// Class describing sample of type primary cell                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XPrimaryCellInfo::XPrimaryCellInfo()
                 :XSampleInfo()
{
   // Default primary cell information constructor
   if(kCS) cout << "---XPrimaryCellInfo::XPrimaryCellInfo(default)------" << endl;

   fIsolationDate = 0;
   fDescription   = "";
}//Constructor

//______________________________________________________________________________
XPrimaryCellInfo::XPrimaryCellInfo(const char *name, const char *type,  
                  Long_t date, const char *description, const char *sex,
                  const char *pheno, const char *geno, const char *extract,
                  Bool_t isXeno, const char *xenostrain, const char *xenosex,
                  Double_t xenoage, const char *xageunits, const char *comment)
                 :XSampleInfo(name, type, sex, pheno, geno, extract, isXeno,
                  xenostrain, xenosex, xenoage, xageunits, comment)
{
   // Primary cell information constructor
   if(kCS) cout << "---XPrimaryCellInfo::XPrimaryCellInfo------" << endl;

   fIsolationDate = date;
   fDescription   = description;
}//Constructor

//______________________________________________________________________________
XPrimaryCellInfo::XPrimaryCellInfo(const XPrimaryCellInfo &info) 
                 :XSampleInfo(info), fDescription(info.fDescription)
{
   // Primary cell information copy constructor
   if(kCS) cout << "---XPrimaryCellInfo::XPrimaryCellInfo(copy)------" << endl;

   fIsolationDate = info.fIsolationDate;
}//CopyConstructor

//______________________________________________________________________________
XPrimaryCellInfo& XPrimaryCellInfo::operator=(const XPrimaryCellInfo& rhs)
{
   // Primary cell information assignment operator.
   if(kCS) cout << "---XPrimaryCellInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fIsolationDate = rhs.fIsolationDate;
      fDescription   = rhs.fDescription;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XPrimaryCellInfo::~XPrimaryCellInfo()
{
   // Primary cell information destructor
   if(kCS) cout << "---XPrimaryCellInfo::~XPrimaryCellInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTissueInfo                                                          //
//                                                                      //
// Class describing sample of type tissue                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTissueInfo::XTissueInfo()
            :XSampleInfo()
{
   // Default tissue information constructor
   if(kCS) cout << "---XTissueInfo::XTissueInfo(default)------" << endl;

   fDevelopment  = "";
   fMorphology   = "";
   fDisease      = "";
   fDiseaseStage = "";
   fDonorAge     = 0;
   fAgeUnits     = "";
   fDonorStatus  = "";
}//Constructor

//______________________________________________________________________________
XTissueInfo::XTissueInfo(const char *name, const char *type, const char *development, 
             const char *morphology, const char *disease, const char *stage,
             Double_t donorage, const char *ageunits, const char *status,
             const char *sex, const char *pheno, const char *geno, const char *extract,
             Bool_t isXeno, const char *xenostrain, const char *xenosex,
             Double_t xenoage, const char *xageunits, const char *comment)
            :XSampleInfo(name, type, sex, pheno, geno, extract, isXeno,
             xenostrain, xenosex, xenoage, xageunits, comment)
{
   // Tissue information constructor
   if(kCS) cout << "---XTissueInfo::XTissueInfo------" << endl;

   fDevelopment  = development;
   fMorphology   = morphology;
   fDisease      = disease;
   fDiseaseStage = stage;
   fDonorAge     = donorage;
   fAgeUnits     = ageunits;
   fDonorStatus  = status;
}//Constructor

//______________________________________________________________________________
XTissueInfo::XTissueInfo(const XTissueInfo &info) 
            :XSampleInfo(info), fDevelopment(info.fDevelopment),
             fMorphology(info.fMorphology), fDisease(info.fDisease),
             fDiseaseStage(info.fDiseaseStage), fAgeUnits(info.fAgeUnits),
             fDonorStatus(info.fDonorStatus)
{
   // Tissue information copy constructor
   if(kCS) cout << "---XTissueInfo::XTissueInfo(copy)------" << endl;

   fDonorAge = info.fDonorAge;
}//CopyConstructor

//______________________________________________________________________________
XTissueInfo& XTissueInfo::operator=(const XTissueInfo& rhs)
{
   // Tissue information assignment operator.
   if(kCS) cout << "---XTissueInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fDevelopment  = rhs.fDevelopment;
      fMorphology   = rhs.fMorphology;
      fDisease      = rhs.fDisease;
      fDiseaseStage = rhs.fDiseaseStage;
      fDonorAge     = rhs.fDonorAge;
      fAgeUnits     = rhs.fAgeUnits;
      fDonorStatus  = rhs.fDonorStatus;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XTissueInfo::~XTissueInfo()
{
   // Tissue information destructor
   if(kCS) cout << "---XTissueInfo::~XTissueInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XBiopsyInfo                                                          //
//                                                                      //
// Class describing sample of type biopsy                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XBiopsyInfo::XBiopsyInfo()
            :XSampleInfo()
{
   // Default biopsy information constructor
   if(kCS) cout << "---XBiopsyInfo::XBiopsyInfo(default)------" << endl;

   fMorphology   = "";
   fDisease      = "";
   fDiseaseStage = "";
   fDonorAge     = 0;
   fAgeUnits     = "";
   fDonorStatus  = "";
}//Constructor

//______________________________________________________________________________
XBiopsyInfo::XBiopsyInfo(const char *name, const char *type,
             const char *morphology, const char *disease, const char *stage,
             Double_t donorage, const char *ageunits, const char *status,
             const char *sex, const char *pheno, const char *geno, const char *extract,
             Bool_t isXeno, const char *xenostrain, const char *xenosex,
             Double_t xenoage, const char *xageunits, const char *comment)
            :XSampleInfo(name, type, sex, pheno, geno, extract, isXeno,
             xenostrain, xenosex, xenoage, xageunits, comment)
{
   // Biopsy information constructor
   if(kCS) cout << "---XBiopsyInfo::XBiopsyInfo------" << endl;

   fMorphology   = morphology;
   fDisease      = disease;
   fDiseaseStage = stage;
   fDonorAge     = donorage;
   fAgeUnits     = ageunits;
   fDonorStatus  = status;
}//Constructor

//______________________________________________________________________________
XBiopsyInfo::XBiopsyInfo(const XBiopsyInfo &info) 
            :XSampleInfo(info),fMorphology(info.fMorphology),fDisease(info.fDisease),
             fDiseaseStage(info.fDiseaseStage),fAgeUnits(info.fAgeUnits),
             fDonorStatus(info.fDonorStatus)
{
   // Biopsy information copy constructor
   if(kCS) cout << "---XBiopsyInfo::XBiopsyInfo(copy)------" << endl;

   fDonorAge = info.fDonorAge;
}//CopyConstructor

//______________________________________________________________________________
XBiopsyInfo& XBiopsyInfo::operator=(const XBiopsyInfo& rhs)
{
   // Biopsy information assignment operator.
   if(kCS) cout << "---XBiopsyInfo::operator=------" << endl;

   if (this != &rhs) {
      XSampleInfo::operator=(rhs);
      fMorphology   = rhs.fMorphology;
      fDisease      = rhs.fDisease;
      fDiseaseStage = rhs.fDiseaseStage;
      fDonorAge     = rhs.fDonorAge;
      fAgeUnits     = rhs.fAgeUnits;
      fDonorStatus  = rhs.fDonorStatus;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XBiopsyInfo::~XBiopsyInfo()
{
   // Biopsy information destructor
   if(kCS) cout << "---XBiopsyInfo::~XBiopsyInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreatmentInfo                                                       //
//                                                                      //
// Class describing a single treatment                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreatmentInfo::XTreatmentInfo()
               :XDataTypeInfo()
{
   // Default treatment info constructor
   if(kCS) cout << "---XTreatmentInfo::XTreatmentInfo(default)------" << endl;

   fConcentration  = 0;
   fConcUnit       = "";
   fTime           = 0;
   fTimeUnit       = "";
   fAdministration = "";
}//Constructor

//______________________________________________________________________________
XTreatmentInfo::XTreatmentInfo(const char *name, const char *type, Double_t conc,
                const char *concunit, Double_t time, const char *timeunit,
                const char *admin, const char *comment)
               :XDataTypeInfo("Treatment", name, type, comment)
{
   // Treatment info constructor
   if(kCS) cout << "---XTreatmentInfo::XTreatmentInfo------" << endl;

   fConcentration  = conc;
   fConcUnit       = concunit;
   fTime           = time;
   fTimeUnit       = timeunit;
   fAdministration = admin;
}//Constructor

//______________________________________________________________________________
XTreatmentInfo::XTreatmentInfo(const XTreatmentInfo &info) 
               :XDataTypeInfo(info), fConcUnit(info.fConcUnit),
                fTimeUnit(info.fTimeUnit), fAdministration(info.fAdministration)
{
   // Treatment info copy constructor
   if(kCS) cout << "---XTreatmentInfo::XTreatmentInfo(copy)------" << endl;

   fConcentration = info.fConcentration;
   fTime          = info.fTime;
}//CopyConstructor

//______________________________________________________________________________
XTreatmentInfo& XTreatmentInfo::operator=(const XTreatmentInfo& rhs)
{
   // Treatment info assignment operator.
   if(kCS) cout << "---XTreatmentInfo::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeInfo::operator=(rhs);
      fConcentration  = rhs.fConcentration;
      fConcUnit       = rhs.fConcUnit;
      fTime           = rhs.fTime;
      fTimeUnit       = rhs.fTimeUnit;
      fAdministration = rhs.fAdministration;
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XTreatmentInfo::~XTreatmentInfo()
{
   // Treatment info destructor
   if(kCS) cout << "---XTreatmentInfo::~XTreatmentInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreatmentList                                                       //
//                                                                      //
// Class describing sample treatments                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreatmentList::XTreatmentList()
               :XDataTypeList()
{
   // Default treatment list constructor
   if(kCS) cout << "---XTreatmentList::XTreatmentList(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XTreatmentList::XTreatmentList(const char *type, const char *comment)
               :XDataTypeList("Treatments", "TreatmentList", type, comment)
{
   // Treatment list constructor
   if(kCS) cout << "---XTreatmentList::XTreatmentList------" << endl;

}//Constructor

//______________________________________________________________________________
XTreatmentList::XTreatmentList(const XTreatmentList &list) 
               :XDataTypeList(list)
{
   // Treatment list copy constructor
   if(kCS) cout << "---XTreatmentList::XTreatmentList(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XTreatmentList& XTreatmentList::operator=(const XTreatmentList& rhs)
{
   // Treatment list assignment operator.
   if(kCS) cout << "---XTreatmentList::operator=------" << endl;

   if (this != &rhs) {
      XDataTypeList::operator=(rhs);
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XTreatmentList::~XTreatmentList()
{
   // Treatment list destructor
   if(kCS) cout << "---XTreatmentList::~XTreatmentList------" << endl;

}//Destructor
