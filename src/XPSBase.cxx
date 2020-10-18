// File created: 05/18/2002                          last modified: 01/24/2011
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

/******************************************************************************
* Major Revision History:
* May 2002 - Initial versions finished
* Jan 2003 - Support to read data in XML format.
*          - Add function AddTree() to XManager.
* Nov 2004 - Add classes XTreeInfo and XAlgorithm.
* Apr 2005 - Support to read binary data files.
* May 2005 - Replace class XContent with class XFolder (designed in Oct 2001).
*            (necessary to store database relevant information)
*          - Add treename to method Import()
*          - Fuse libraries libXPSUtils and libXPSBase => library libXPSBase.
* Jul 2005 - Replace TTree::AddFriend() with TList fTrees.
*          - Store trees for every XTreeSet in separate TDirectory.
* Oct 2005 - Add fDataType to XManager; set title of TDirectory to fDataType.
* Jan 2006 - To decrease size of saved trees replace  Int_t with Short_t for 
*            mask/flag in XAlgorithm and derived classes
*          - Add parameter filename to InitAlgorithm() and field fFile and
*            method NewFile() to XAlgorithm to allow saving certain data such
*            as e.g. temporary trees in a (temporary) file
* Apr 2006 - Add fManager to XTreeSet to allow use of XManager::HandleError().
* May 2007 - Add fOption to XManager to allow updating root file from "R".
* Sep 2007 - Change all msk arrays from "Short_t *msk" to "Int_t *msk".
*          - Add support for  verbose messages in XManager constructor
* Feb 2008 - Add methods GetTreeHeader() and GetTreeInfo() to XManager.
*          - Adapt source code to compile with MS VC++ on WinXP
* Jul 2009 - Add methods Set/GetBufSize() to change bufsize of tree branches
* Jan 2010 - Add further methods GetOption() to XAlgorithm.
*          - Move methods from XHybridizer to XAlgorithm.
* Jan 2011 - Add method ExportTreeUserInfo() to export tree user info.
*
******************************************************************************/

/******************************************************************************
* Note on memory management:
* Since ROOT does not use exception handling I have decided also not to use it.
* However, to be able to control memory management, especially when defining a
* series of large arrays, I decided to use the following memory initialization
* throughout my program: 
*    if (!(arr  = new (nothrow) Float_t[size])) {err = errInitMemory; goto cleanup;}
*  cleanup:
*    if (arr) {delete [] arr; arr = 0;}
*    etc
* 1. "new (nothrow)" is necessary otherwise a memory error will not return zero.
* 2. Yes, I decided to use "goto cleanup", since doing cleanup at the end of a
*    method turned out to be the most effective way to cleanup before exiting
*    a method, and doing cleanup before exiting the program.
* 3. "if (arr)" ensures that only already initialized arrays (objects) will be
*    deleted.
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <new>  //needed for new (nothrow)

#include "TDatime.h"
#include "TError.h"
#include "TFriendElement.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

#include "XPSBase.h"

ClassImp(XIdxString);
ClassImp(XLdxString);
ClassImp(XFolder);
ClassImp(XSetting);
ClassImp(XTreeInfo);
ClassImp(XTreeHeader);
ClassImp(XTreeSet);
ClassImp(XAlgorithm);
ClassImp(XManager);

Int_t  XManager::fgBufSize = 32000;
Int_t  XManager::fgVerbose = 1;
Bool_t XTreeSet::fgPrintHeader = kTRUE;

//XFolders
const int kMAXDEPTH = 64;
static const char *gFolderD[kMAXDEPTH];
static Int_t gFolderLevel = -1;
static char  gFolderPath[512];
const  char* kRootFolder = "Content"; //??
///////////
// NOTE: In root file containing database import all content XFolders
//       from all xxx_cel.root files as subfolders
///////////

const char* kContent = "Content";

//debug: print function names
const Bool_t kCS  = 0; 
const Bool_t kCSa = 0; 


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XFolder                                                              //
//                                                                      //
// Class containing info for database about objects stored in TFile     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XFolder::XFolder()
        :TFolder()
{
   // Default Folder constructor 
   if(kCS) cout << "---XFolder::XFolder(default)------" << endl;

   fType = "";
}//Constructor

//______________________________________________________________________________
XFolder::XFolder(const char *name, const char *title, const char *type,
         Bool_t isPublic, Bool_t newlist)
        :TFolder(name, title)
{
   // Folder constructor
   if(kCS) cout << "---XFolder::XFolder------" << endl;

   fType     = type;
   fIsPublic = isPublic;
   if (newlist) fFolders = new TList();
}//Constructor

//______________________________________________________________________________
XFolder::~XFolder()
{
   // Folder destructor
   if(kCS) cout << "---XFolder::~XFolder------" << endl;

}//Destructor

//______________________________________________________________________________
XFolder *XFolder::AddFolder(const char *name, const char *title, const char *type,
                  TCollection *collection)
{
   // Create new folder and add it to list of folders
   if(kCS) cout << "------XFolder::AddFolder------" << endl;

   // adapted from TFolder::AddFolder
   if (strchr(name, '/')) {
      ::Error("TFolder::TFolder", "folder name cannot contain a slash: %s", name); 
      return 0;
   }
   if (strlen(GetName()) == 0) {
      ::Error("XFolder::XFolder", "folder name cannot be \"\"");
      return 0;
   }
   XFolder *folder = new XFolder();
   folder->SetName(name);
   folder->SetTitle(title);
   folder->SetType(type);
   if (!fFolders) fFolders = new TList();
   fFolders->Add(folder);

   if (collection) folder->fFolders = collection;
   else            folder->fFolders = new TList();
   return folder;
}//AddFolder

//______________________________________________________________________________
const char *XFolder::FindFullPathName(const char *name) const
{
   // Return the full pathname for name
   if(kCS) cout << "------XFolder::FindFullPathName------" << endl;

   // adapted from TFolder::FindFullPathName
   TObject *obj = this->FindObject(name);
   if (obj || !fFolders) {
      gFolderLevel++;
      gFolderD[gFolderLevel] = GetName();
      gFolderPath[0] = '/';
      gFolderPath[1] = 0;
      for (Int_t l=0; l<=gFolderLevel; l++) {
         strcat(gFolderPath, "/");
         strcat(gFolderPath, gFolderD[l]);
      }
      strcat(gFolderPath, "/");
      strcat(gFolderPath, name);
      gFolderLevel = -1;
      return gFolderPath;
   }

   if (name[0] == '/') return 0;

   TIter next(fFolders);
   XFolder *folder;
   const char *found;
   gFolderLevel++;
   gFolderD[gFolderLevel] = GetName();
   while ((obj=next())) {
      if (!obj->InheritsFrom(XFolder::Class())) continue;
      if (obj->InheritsFrom(TClass::Class())) continue;
      folder = (XFolder*)obj;
      found = folder->FindFullPathName(name);
      if (found) return found;
   }
   gFolderLevel--;
   return 0;
}//FindFullPathName

//______________________________________________________________________________
TObject *XFolder::FindObject(const char *name) const
{
   // Search object in tree of folders inside current folder
   // Name may be of the forms:
   //   A, specify full pathname starting at RootFolder of file
   //      e.g. //RootFolder/xxx/yyy/name 
   //
   //   B, specify full pathname starting at top folder of current directory
   //      e.g. /xxx/yyy/name 
   //
   //   C, Specify a pathname relative to this folder
   //      e.g. xxx/yyy/name
   //     name
   if(kCS) cout << "------XFolder::FindObject------" << endl;

   // adapted from TFolder::FindObject
   XFolder *folder = 0;
   if (!fFolders) return 0;
   if (name == 0) return 0;
   if (name[0] == '/') {
      if (name[1] == '/') {
         TString str = "//" + TString(kRootFolder) + "/";
         if (!strstr(name,str.Data())) return 0;
         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name + str.Length());
      } else {
         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name+1);
      }//if
   }//if

   char cname[kBuf4096];
   strcpy(cname,name);
   TObject *obj;
   char *slash = strchr(cname,'/');
   if (slash) {
      *slash = 0;
      obj = fFolders->FindObject(cname);
      if (!obj) return 0;
      return obj->FindObject(slash+1);
   } else {
      return fFolders->FindObject(name);
   }//if
}//FindObject

//______________________________________________________________________________
TObject *XFolder::FindObject(const char *name, Bool_t ignoreCase) const
{
   // Search object in tree of folders inside current folder
   // Name may be of the forms:
   //   A, specify full pathname starting at RootFolder of file
   //      e.g. //RootFolder/xxx/yyy/name 
   //
   //   B, specify full pathname starting at top folder of current directory
   //      e.g. /xxx/yyy/name 
   //
   //   C, Specify a pathname relative to this folder
   //      e.g. xxx/yyy/name
   //     name
   // If ignoreCase == kTRUE then ignore case of name (but not of path)
   if(kCS) cout << "------XFolder::FindObject(ignore)------" << endl;

   // adapted from TFolder::FindObject
   XFolder *folder = 0;
   if (!fFolders) return 0;
   if (name == 0) return 0;
   if (name[0] == '/') {
      if (name[1] == '/') {
         TString str = "//" + TString(kRootFolder) + "/";
         if (!strstr(name,str.Data())) return 0;

         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name + str.Length());
      } else {
         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name+1);
      }//if
   }//if

   char cname[kBuf4096];
   strcpy(cname,name);
   TObject *obj;
   char *slash = strchr(cname,'/');
   if (slash) {
      *slash = 0;
      obj = fFolders->FindObject(cname);
      if (!obj) return 0;
      return obj->FindObject(slash+1);
   } else if (ignoreCase) {
      TString str1 = TString(name);
      str1.ToLower();

      TIter next(fFolders);
      TObject *obj;
      while ((obj = next())) {
         TString str2 = TString(obj->GetName());
         str2.ToLower();

         if (!strcmp(str1.Data(), str2.Data())) return obj;
      }//while
   } else {
      return fFolders->FindObject(name);
   }//if

   return 0;
}//FindObject

//______________________________________________________________________________
TObject *XFolder::FindObject(const char *name, const char *classname) const
{
   // Search object inherited from class classname in tree of folders inside
   // the current folder
   // Name may be of the forms:
   //   A, specify full pathname starting at RootFolder of file
   //      e.g. //RootFolder/xxx/yyy/name 
   //
   //   B, specify full pathname starting at top folder of current directory
   //      e.g. /xxx/yyy/name 
   //
   //   C, Specify a pathname relative to this folder
   //      e.g. xxx/yyy/name
   //     name
   if(kCS) cout << "------XFolder::FindObject(n,c)------" << endl;

   // adapted from TFolder::FindObject
   XFolder *folder = 0;
   if (!fFolders) return 0;
   if (name == 0) return 0;
   if (name[0] == '/') {
      if (name[1] == '/') {
         TString str = "//" + TString(kRootFolder) + "/";
         if (!strstr(name,str.Data())) return 0;

         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name + str.Length(), classname);
      } else {
         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name+1, classname);
      }//if
   }//if

   char cname[kBuf4096];
   strcpy(cname,name);
   TObject *obj;
   char *slash = strchr(cname,'/');
   if (slash) {
      *slash = 0;
      obj = fFolders->FindObject(cname);
      if (!obj) return 0;
      if (obj->InheritsFrom("XFolder")) {
         return ((XFolder*)obj)->FindObject(slash+1, classname);
      } else {
         return obj->FindObject(slash+1);
      }//if
   } else {
      TIter next(fFolders);
      TObject *obj;
      while ((obj = next())) {
         if (obj->InheritsFrom(classname)) {
            if (strcmp(name, obj->GetName()) == 0) return obj;
         }//if
      }//while
   }//if

   return 0;
}//FindObject

//______________________________________________________________________________
TObject *XFolder::FindObject(const char *name, const char *title,
                  const char *classname) const
{
   // Search object inherited from class classname in tree of folders inside
   // the current folder
   // Name may be of the forms:
   //   A, specify full pathname starting at RootFolder of file
   //      e.g. //RootFolder/xxx/yyy/name 
   //
   //   B, specify full pathname starting at top folder of current directory
   //      e.g. /xxx/yyy/name 
   //
   //   C, Specify a pathname relative to this folder
   //      e.g. xxx/yyy/name
   //     name
   if(kCS) cout << "------XFolder::FindObject(n,t,c)------" << endl;

   // adapted from TFolder::FindObject
   XFolder *folder = 0;
   if (!fFolders) return 0;
   if (name == 0) return 0;
   if (name[0] == '/') {
      if (name[1] == '/') {
         TString str = "//" + TString(kRootFolder) + "/";
         if (!strstr(name,str.Data())) return 0;

         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name + str.Length(), title, classname);
      } else {
         folder = (XFolder*)(gDirectory->Get(kRootFolder));
         if (!folder) return 0;
         return folder->FindObject(name+1, title, classname);
      }//if
   }//if

   char cname[kBuf4096];
   strcpy(cname,name);
   TObject *obj;
   char *slash = strchr(cname,'/');
   if (slash) {
      *slash = 0;
      obj = fFolders->FindObject(cname);
      if (!obj) return 0;
      if (obj->InheritsFrom("XFolder")) {
         return ((XFolder*)obj)->FindObject(slash+1, title, classname);
      } else {
         return obj->FindObject(slash+1);
      }//if
   } else {
      TIter next(fFolders);
      TObject *obj;
      while ((obj = next())) {
         if (obj->InheritsFrom(classname)) {
            if ((strcmp(name,  obj->GetName())  == 0) &&
                (strcmp(title, obj->GetTitle()) == 0)) return obj;
         }//if
      }//while
   }//if

   return 0;
}//FindObject

//______________________________________________________________________________
TObject *XFolder::FindObjectAny(const char *name) const
{
   // Find first object starting at this folder 
   if(kCS) cout << "------XFolder::FindObjectAny------" << endl;

   // adapted from TFolder::FindObjectAny
   TObject *obj = this->FindObject(name);
   if (obj || !fFolders) return obj;
   if (name[0] == '/') return 0;

   TIter next(fFolders);
   XFolder *folder;
   TObject *found;
   if (gFolderLevel >= 0) gFolderD[gFolderLevel] = GetName();
   while ((obj=next())) {
      if (!obj->InheritsFrom(XFolder::Class())) continue;
      if (obj->IsA() == TClass::Class()) continue;
      folder = (XFolder*)obj;
      found = folder->FindObjectAny(name);
      if (found) return found;
   }

   return 0;
}//FindObjectAny


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreeInfo                                                            //
//                                                                      //
// Base class containing info about tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreeInfo::XTreeInfo() 
          :TNamed()
{
   // Default TreeInfo constructor
   if(kCS) cout << "---XTreeInfo::XTreeInfo(default)------" << endl;

   fSetName  = "";
   fSetClass = "";
   fOption   = "";
}//Constructor

//______________________________________________________________________________
XTreeInfo::XTreeInfo(const char *name, const char *title) 
          :TNamed(name, title)
{
   // Normal TreeInfo constructor
   if(kCS) cout << "---XTreeInfo::XTreeInfo------" << endl;

   fSetName  = "";
   fSetClass = "";
   fOption   = "";
}//Constructor

//______________________________________________________________________________
XTreeInfo::~XTreeInfo()
{
   // TreeInfo destructor
   if(kCS) cout << "---XTreeInfo::~XTreeInfo------" << endl;

}//Destructor


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreeHeader                                                          //
//                                                                      //
// Class storing tree information in header                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreeHeader::XTreeHeader() 
            :TObjString("")
{
   // Default TreeHeader constructor
   if(kCS) cout << "---XTreeHeader::XTreeHeader(default)------" << endl;

   fInfile = "";
   fType   = "";
   fID     = 0;
   fNPar   = 0;
   fPars   = 0;
}//Constructor

//______________________________________________________________________________
XTreeHeader::XTreeHeader(const char *str, Int_t treeid) 
            :TObjString(str)
{
   // Normal TreeHeader constructor
   if(kCS) cout << "---XTreeHeader::XTreeHeader------" << endl;

   fInfile = "";
   fType   = "";
   fID     = treeid;
   fDatime.Set();
   fNPar   = 0;
   fPars   = 0;
}//Constructor

//______________________________________________________________________________
XTreeHeader::~XTreeHeader()
{
   // TreeHeader destructor
   if(kCS) cout << "---XTreeHeader::~XTreeHeader------" << endl;

   if (fPars) {delete [] fPars; fPars = 0;}
   fNPar = 0;
}//Destructor

//______________________________________________________________________________
void XTreeHeader::SetParameters(Int_t npar, Double_t *pars)
{
   // Set parameters used by algorithm
   if(kCS) cout << "------XTreeHeader::SetParameters------" << endl;

   if (npar == 0) return;
   if (pars == 0) return;

   fNPar = npar;
   fPars = new Double_t[fNPar];
   memcpy(fPars, pars, npar*sizeof(Double_t));
}//SetParameters


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSetting                                                             //
//                                                                      //
// Class for initialization of parameter settings                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSetting::XSetting()
         :TNamed()
{
   // Default Setting constructor
   if(kCS) cout << "---XSetting::XSetting(default)------" << endl;

   fNA    = -1;
   fHasNA = kFALSE;
}//Constructor

//______________________________________________________________________________
XSetting::XSetting(const char *type, const char *infile)
         :TNamed(type, infile)
//?? infile not here?
{
   // Normal Setting constructor
   if(kCS) cout << "---XSetting::XSetting------" << endl;

   fNA    = -1;
   fHasNA = kFALSE;
}//Constructor

//______________________________________________________________________________
XSetting::~XSetting()
{
   // Setting destructor
   if(kCS) cout << "---XSetting::~XSetting------" << endl;

}//Destructor

//______________________________________________________________________________
void XSetting::ResetAlgorithm(const char * /*name*/, const char * /*type*/)
{
   // Reset algorithm with name and optional with type
   if(kCS) cout << "------XSetting::ResetAlgorithm------" << endl;

   cout << "ResetAlgorithm() is not implemented." << endl;
}//ResetAlgorithm


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTreeSet                                                             //
//                                                                      //
// Base class containing info about tree sets stored in TFile           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XTreeSet::XTreeSet() 
         :TNamed()
{
   // Default TreeSet constructor
   if(kCS) cout << "---XTreeSet::XTreeSet(default)------" << endl;

   fHeaders    = new TList(); //since default constructor called for permanent obj 
   fTrees      = new TList();
   fSelections = 0;
   fTrash      = 0;
   fFile       = 0;
   fTree       = 0;
   fSetting    = 0;
   fManager    = 0;
   fInfile     = "";
   fTreeName   = "";
   fAsXML      = kFALSE;
}//Constructor

//______________________________________________________________________________
XTreeSet::XTreeSet(const char *name, const char *title) 
         :TNamed(name, title)
{
   // Normal TreeSet constructor
   if(kCS) cout << "---XTreeSet::XTreeSet------" << endl;

   fHeaders    = new TList();
   fTrees      = new TList();
   fSelections = 0;  //created on demand only
   fTrash      = new TList();
   fFile       = 0;
   fTree       = 0;
   fSetting    = 0;
   fManager    = 0;
   fInfile     = "";
   fTreeName   = "";
   fAsXML      = kFALSE;
}//Constructor

//______________________________________________________________________________
XTreeSet::XTreeSet(const XTreeSet &treeset) 
         :TNamed(treeset)
{
   // TreeSet copy constructor
   if(kCS) cout << "---XTreeSet::XTreeSet(copy)------" << endl;

//to do: need to copy fHeaders!!! and its contents?
}//CopyConstructor

//______________________________________________________________________________
XTreeSet& XTreeSet::operator=(const XTreeSet& rhs)
{
   // TreeSet assignment operator.
   if(kCS) cout << "---XTreeSet::operator=------" << endl;

   if (this != &rhs) {
      TNamed::operator=(rhs);
//to do: need to copy fHeaders!!!
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XTreeSet::~XTreeSet()
{
   // TreeSet destructor
   if(kCS) cout << "---XTreeSet::~XTreeSet------" << endl;

   if (fHeaders)    {fHeaders->Delete();    delete fHeaders;    fHeaders    = 0;}
   if (fSelections) {fSelections->Delete(); delete fSelections; fSelections = 0;}
   if (fTrash)      {fTrash->Delete();      delete fTrash;      fTrash      = 0;}

   fTrees->Clear("nodelete");  //do not delete trees
   SafeDelete(fTrees);

   fManager = 0;
   fSetting = 0;
   fTree    = 0;
   fFile    = 0;
}//Destructor

//______________________________________________________________________________
void XTreeSet::AddTree(TTree *tree)
{
   // Add tree to list fTrees
   if(kCS) cout << "------XTreeSet::AddTree------" << endl;

   fTrees->Add(tree);
}//AddTree

//______________________________________________________________________________
void XTreeSet::AddTreeInfo(TTree *tree, XTreeInfo *info)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XTreeSet::AddTreeInfo------" << endl;

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   tree->GetUserInfo()->Add(info);
}//AddTreeInfo

//______________________________________________________________________________
void XTreeSet::AddTreeInfo(TTree *tree, const char *name, Option_t *option)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XTreeSet::AddTreeInfo------" << endl;

   XTreeInfo *info = new XTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   tree->GetUserInfo()->Add(info);
}//AddTreeInfo

//______________________________________________________________________________
void XTreeSet::AddTreeHeader(const char *treename, Int_t treeid)
{
   // Add tree header with treename to list of trees
   // Note: Need to call AddTreeHeader() after calling WriteTree()
   // Note: For treeid > 0 treename will only be added to selections.
   if(kCS) cout << "------XTreeSet::AddTreeHeader------" << endl;

   if (treeid > 0) {
      this->Select(treename, treeid);
   } else {
      this->Select(treename, treeid);

      TString exten = Path2Name(treename,".","");

      XTreeHeader *header = new XTreeHeader(treename, treeid);
      header->SetInfile(fInfile);
      header->SetType(exten.Data());
      fHeaders->Add(header);
   }//if
}//AddTreeHeader

//______________________________________________________________________________
void XTreeSet::AddTreeHeader(const char *treename, const char *treetype,
                 Int_t treeid, Int_t npar, Double_t *pars)
{
   // Add tree header with treename to list of trees
   // Note: Need to call AddTreeHeader() after calling WriteTree()
   if(kCS) cout << "------XTreeSet::AddTreeHeader------" << endl;

   if (treeid > 0) {
      this->Select(treename, treeid);
   } else {
      TString type = treetype;
      if (strcmp(treetype, "") == 0) type = Path2Name(treename,".","");

      XTreeHeader *header = new XTreeHeader(treename, treeid);
      header->SetInfile(fInfile);
      header->SetType(type.Data());
      header->SetParameters(npar, pars);
      fHeaders->Add(header);
   }//if
}//AddTreeHeader

//______________________________________________________________________________
void XTreeSet::RemoveTreeHeader(XTreeHeader *header)
{
   // Remove tree header from list of trees
   if(kCS) cout << "------XTreeSet::RemoveTreeHeader------" << endl;

   fHeaders->Remove(header);

   SafeDelete(header);
}//RemoveTreeHeader

//______________________________________________________________________________
void XTreeSet::RemoveTreeHeaders()
{
   // Remove all tree headers from list of trees
   if(kCS) cout << "------XTreeSet::RemoveTreeHeaders------" << endl;

   fHeaders->Delete();
}//RemoveTreeHeaders

//______________________________________________________________________________
Int_t XTreeSet::Export(const char *exten, const char *varlist,
                       ofstream &output, const char *sep)
{
   // Export selected variables given in varlist which are stored in tree
   // in output, separated by sep
   if(kCS) cout << "------XTreeSet::Export(tree)------" << endl;

   Int_t err     = errNoErr;
   Int_t numsels = 0;

   if (fSelections) numsels = fSelections->GetSize();
   if (numsels == 0) return errGetTree;

// Fill array with tree names
   TString *arrName = 0;
   if (!(arrName = new (nothrow) TString[numsels])) return errInitMemory;

   Int_t numtrees = 0;
   for (Int_t i=0; i<numsels; i++) {
      TObjString *selstr = (TObjString*)(fSelections->At(i));
      TString treename = selstr->GetString();
      TString treexten = Path2Name(treename,".",";");
      if (strcmp(treexten.Data(), exten) == 0) {
         arrName[numtrees] = selstr->GetString();
         numtrees++;
      }//if
   }//for_i

// Export tree(s) with extension exten
   err = this->ExportTree(exten, numtrees, arrName, varlist, output, sep);

// Clean up
   delete [] arrName;

   return err;
}//Export

//______________________________________________________________________________
Int_t XTreeSet::ExportTrees(const char *exten, const char *varlist, 
                ofstream &output, const char *sep)
{
   // Export variables from varlist for all trees in current tree set
   // with extension exten
   if(kCS) cout << "------XTreeSet::ExportTrees------" << endl;

   Int_t err = errNoErr;
   Int_t n   = fHeaders->GetSize();

// Initialize array of treenames
   TString *arrName = 0;
   if (!(arrName = new (nothrow) TString[n])) return errInitMemory;

// Clear tree list first 
   fTrees->Clear();

// Add trees to list fTrees
   n = 0;
   TIter next(fHeaders);
   XTreeHeader *header;
   while ((header = (XTreeHeader*)next())) {
      TString treename = header->GetString();
      TString treexten = Path2Name(treename,".",";");
      if (strcmp(treexten.Data(), exten) == 0) {
         arrName[n] = header->GetString();
         TTree *tree = (TTree*)gDirectory->Get(arrName[n]);
         if (!tree) {delete [] arrName; return errGetTree;}
         fTrees->Add(tree);
         n++;
      }//if
   }//while

   if (n == 0) {
      cerr << "Error: Could not get tree(s) with extension <" << exten << ">."
           << endl;
      return errGetTree;
   }//if

// Export tree(s) with extension exten
   err = this->ExportTree(exten, n, arrName, varlist, output, sep);

// Clean up
   delete [] arrName;

   return err;
}//ExportTrees

//______________________________________________________________________________
Int_t XTreeSet::ExportTree(const char *exten, Int_t n, TString *names,  
                const char *varlist, ofstream &output, const char *sep)
{
   // Export selected variables given in varlist which are stored in trees
   // with treenames names in output
   // If varlist begins with "userinfo" then export tree user info variables
   // listed in varlist after "userinfo", e.g. varlist="userinfo:fMin:fMax"
   if(kCS) cout << "------XTreeSet::ExportTree------" << endl;

   Int_t err = errNoErr;

   if(strcmp(SubString(varlist,":",0).Data(), "userinfo") == 0) {
      err = this->ExportTreeInfo(exten, n, names, varlist, output, sep);
   } else {
      if (fAsXML) err = this->ExportTreeXML(exten, n, names, varlist, output, sep);
      else        err = this->ExportTreeType(exten, n, names, varlist, output, sep);
   }//if

   return err;
}//ExportTree

//______________________________________________________________________________
Int_t XTreeSet::Import(ifstream &input, Option_t *option, const char *sep,
                char delim, Int_t split)
{
   // Import data from input. If option = "UPDATE" replace existing data.
   // If IsBinaryFile() != 0 then read binary header and data
   // Note: Since IsBinaryFile() returns an integer, an inherited method Import()
   //       would allow to read multiple differnt binary files
   if(kCS) cout << "------XTreeSet::Import------" << endl;

   Int_t err = errNoErr;

   if (this->IsBinaryFile(input)) {
      if (err == errNoErr) err = this->ReadBinaryHeader(input, sep, delim);
      if (err == errNoErr) err = this->ReadBinaryData(input, option, sep, delim, split);
   } else {
      input.close();
      input.open(fInfile.Data(), ios::in);  // reopen as text file
      if (!input) return errOpenFile;

      if (err == errNoErr) err = this->ReadHeader(input, sep, delim);
      if (err == errNoErr) err = this->ReadData(input, option, sep, delim, split);
   }//if

   return err;
}//Import

//______________________________________________________________________________
Int_t XTreeSet::ImportXML(ifstream &input, Option_t *option, const char *sep,
                char delim, Int_t split)
{
   // Import data from input in XML-format.
   // If option = "UPDATE" replace existing data. 
   if(kCS) cout << "------XTreeSet::ImportXML------" << endl;

   Int_t err = errNoErr;

   if (!err) err = this->ReadXMLHeader(input, sep, delim);
   if (!err) err = this->ReadXMLData(input, option, sep, delim, split);

   return err;
}//Import

//______________________________________________________________________________
Int_t XTreeSet::Initialize(TFile *file, XSetting *setting, const char *infile,
                const char *treename)
{
   // Initialize setting
   if(kCS) cout << "------XTreeSet::Initialize------" << endl;

   if (setting == 0) return errAbort;

   fFile     = file;
   fSetting  = setting;
   fInfile   = infile;
   fTreeName = treename;

   return errNoErr;
}//Initialize

//______________________________________________________________________________
void XTreeSet::Select(const char *name, Int_t id)
{
   // Add selected tree to list fSelections
   if(kCS) cout << "------XTreeSet::Select------" << endl;

   if (!fSelections) fSelections = new TList();

   XIdxString *str = new XIdxString(id, name);
   fSelections->Add(str);
}//Select

//______________________________________________________________________________
Double_t **XTreeSet::CreateTable(Int_t nrow, Int_t ncol)
{
   // Create table
   if(kCS) cout << "------XTreeSet::CreateTable------" << endl;

   Double_t **table = 0;

   if (!(table = new (nothrow) Double_t*[nrow])) return 0;
   for (Int_t k=0; k<nrow; k++) {
      table[k] = 0;
      if (!(table[k] = new (nothrow) Double_t[ncol])) return 0; 
   }//for_i

   return table;
}//CreateTable

//______________________________________________________________________________
void XTreeSet::DeleteTable(Double_t **table, Int_t nrow)
{
   // Delete table
   if(kCS) cout << "------XTreeSet::DeleteTable------" << endl;

   if (table == 0) return;

   for (Int_t k=0; k<nrow; k++) {
      if (table[k]) {delete [] table[k]; table[k] = 0;}
   }//for_k
   delete [] table;
}//DeleteTable

//______________________________________________________________________________
Int_t XTreeSet::WriteTree(TTree *tree, Int_t option, Int_t bufsize)
{
   // Write tree to current directory
   // Note: Need to call AddTreeHeader() after calling WriteTree()
   if(kCS) cout << "------XTreeSet::WriteTree------" << endl;

   Int_t err = errNoErr;

// Check if old tree header exists and optionally remove
   if (option == TObject::kOverwrite) {
      // remove old header from tree list
      TIter next(fHeaders);
      XTreeHeader *header = 0;
      while ((header = (XTreeHeader*)next())) {
         TString oldtree = Path2Name(header->GetString(), dSEP, ";");

         if (strcmp(tree->GetName(), oldtree.Data()) == 0) {
            this->RemoveTreeHeader(header);
            if (XManager::fgVerbose) {
               cout << "Tree name <" << oldtree.Data() << "> is removed from header."
                    << endl;
            }//if
         }//if
      }//while
   }//if

// Write tree to file
   if (tree->Write("", option, bufsize) > 0) {
      err = errNoErr;
   } else {
      cerr << "Error: Could not write tree <" << tree->GetName()
           << "> to directory <" << gDirectory->GetName() << ">." << endl;
      err = errWriteObject;
   }//if

   return err;
}//WriteTree

//______________________________________________________________________________
Int_t XTreeSet::DeleteTree(const char *name, const char *exten, const char *cycle)
{
   // Delete tree name.exten;cycle 
   if(kCS) cout << "------XTreeSet::DeleteTree------" << endl;

// Set size to number of trees to delete
   Int_t size = 1;
   if ((strcmp(name, "*") == 0) || (strcmp(exten, "*") == 0)) {
      size = fHeaders->GetSize();
      if (size == 0) return errGetTree;
   }//if

// Loop over tree list
   TIter next(fHeaders);
   XTreeHeader *header;
   while ((header = (XTreeHeader*)next())) {
      TString treename = header->GetString();
      TString namepart = Path2Name(treename, dSEP, ".");
      TString xtenpart = Path2Name(treename, ".", ";");

      if (((strcmp(name,  namepart.Data())  == 0) || (strcmp(name, "*") == 0)) &&
          ((strcmp(exten, xtenpart.Data()) == 0) || (strcmp(exten, "*") == 0))) {
         if (strcmp(cycle, "") != 0) {
            this->RemoveTreeHeader(header);
            treename = treename + ";" + TString(cycle);
            cout << "Tree <" << treename.Data() << "> is deleted from file." << endl;
         } else {
            cout << "Tree <" << treename.Data() << "> is deleted from memory." << endl;
         }//if

         gDirectory->Delete(treename);
         size--;
      }//if
   }//while

   return size;
}//DeleteTree

//______________________________________________________________________________
TString XTreeSet::FindTree(const char *name)
{
   // Find tree name stored in list fHeaders and return string treename 
   // If name is not found return 0.
   if(kCS) cout << "------XTreeSet::FindTree------" << endl;

   XTreeHeader *header = 0;
   header = (XTreeHeader*)fHeaders->FindObject(name);
   if (header) return header->GetString();
   return TString(0);
}//FindTree

//______________________________________________________________________________
XTreeInfo *XTreeSet::GetTreeInfo(const char *name, TTree *tree)
{
   // Find tree info with name stored in list fUserInfo of tree
   if(kCS) cout << "------XTreeSet::GetTreeInfo------" << endl;

   XTreeInfo *info = 0;
   info = (XTreeInfo*)tree->GetUserInfo()->FindObject(name);
   if (info) return info;
   return 0;
}//GetTreeInfo

//______________________________________________________________________________
XTreeHeader *XTreeSet::GetTreeHeader(const char *name)
{
   // Find tree header stored in list fHeaders
   if(kCS) cout << "------XTreeSet::GetTreeHeader------" << endl;

   XTreeHeader *header = 0;
   header = (XTreeHeader*)fHeaders->FindObject(name);
   if (header) return header;
   return 0;
}//GetTreeHeader


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAlgorithm                                                           //
//                                                                      //
// Base class for algorithms                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XAlgorithm::XAlgorithm()
           :TNamed()
{
   // Default Algorithm constructor
   if(kCS) cout << "---XAlgorithm::XAlgorithm(default)------" << endl;

   fFile        = 0;
   fTreeInfo    = 0;
   fOption      = "";
   fNPar        = 0;
   fPars        = 0;
   fNA          = -1;
   fHasNA       = kFALSE;
   fIsFileOwner = kFALSE;
}//Constructor

//______________________________________________________________________________
XAlgorithm::XAlgorithm(const char *name, const char *type)
           :TNamed(name, type)
{
   // Normal Algorithm constructor
   if(kCS) cout << "---XAlgorithm::XAlgorithm------" << endl;

   fFile        = 0;
   fTreeInfo    = 0;
   fOption      = "";
   fNPar        = 0;
   fPars        = 0;
   fNA          = -1;
   fHasNA       = kFALSE;
   fIsFileOwner = kFALSE;
}//Constructor

//______________________________________________________________________________
XAlgorithm::XAlgorithm(const XAlgorithm &algorithm) 
           :TNamed(algorithm)
{
   // Algorithm copy constructor
   if(kCS) cout << "---XAlgorithm::XAlgorithm(copy)------" << endl;

   fFile        = algorithm.fFile;
   fTreeInfo    = algorithm.fTreeInfo;
   fOption      = algorithm.fOption;
   fNPar        = algorithm.fNPar;
   fNA          = algorithm.fNA;
   fHasNA       = algorithm.fHasNA;
   fIsFileOwner = algorithm.fIsFileOwner;

   fPars = new Double_t[fNPar];
   memcpy(fPars, algorithm.fPars, fNPar*sizeof(Double_t));
}//CopyConstructor

//______________________________________________________________________________
XAlgorithm& XAlgorithm::operator=(const XAlgorithm& rhs)
{
   // Algorithm assignment operator.
   if(kCS) cout << "---XAlgorithm::operator=------" << endl;

   if (this != &rhs) {
      TNamed::operator=(rhs);

      fFile        = rhs.fFile;
      fTreeInfo    = rhs.fTreeInfo;
      fOption      = rhs.fOption;
      fNPar        = rhs.fNPar;
      fNA          = rhs.fNA;
      fHasNA       = rhs.fHasNA;
      fIsFileOwner = rhs.fIsFileOwner;

      if (fPars) {delete [] fPars; fPars = 0;}
      fPars = new Double_t[fNPar];
      memcpy(fPars, rhs.fPars, fNPar*sizeof(Double_t));
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XAlgorithm::~XAlgorithm()
{
   // Algorithm destructor
   if(kCS) cout << "---XAlgorithm::~XAlgorithm------" << endl;

   this->DeleteParameters();
   fTreeInfo = 0;

   if (fIsFileOwner) SafeDelete(fFile);
   fFile = 0;
}//Destructor

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Double_t &/*value1*/, Double_t &/*value2*/, Int_t &/*num*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Int_t * /*x*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Double_t * /*x*/, Int_t * /*idx*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Int_t * /*idx*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*n*/, Double_t * /*x*/, Double_t * /*y*/, Double_t * /*z*/, Int_t * /*idx*/, Int_t * /*msk*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Int_t XAlgorithm::Calculate(Int_t /*nrow*/, Int_t /*ncol*/, Double_t ** /*table*/)
{
   return 0;
}//Calculate

//______________________________________________________________________________
Double_t **XAlgorithm::CreateTable(Int_t nrow, Int_t ncol)
{
   // Create table
   if(kCS) cout << "------XAlgorithm::CreateTable------" << endl;

   Double_t **table = 0;

   if (!(table = new (nothrow) Double_t*[nrow])) return 0;
   for (Int_t k=0; k<nrow; k++) {
      table[k] = 0;
      if (!(table[k] = new (nothrow) Double_t[ncol])) return 0; 
   }//for_i

   return table;
}//CreateTable

//______________________________________________________________________________
void XAlgorithm::DeleteTable(Double_t **table, Int_t nrow)
{
   // Delete table
   if(kCS) cout << "------XAlgorithm::DeleteTable------" << endl;

   if (table == 0) return;

   for (Int_t k=0; k<nrow; k++) {
      if (table[k]) {delete [] table[k]; table[k] = 0;}
   }//for_k
   delete [] table;
}//DeleteTable

//______________________________________________________________________________
TFile *XAlgorithm::NewFile(const char *name, const char *exten)
{
   // Create new root file "name" and return pointer to TFile
   // If name is "tmp" or "tmp_exten", i.e. "tmp.root" or "tmp_exten.root"
   // a temporary file will be created using RECREATE
   if(kCS) cout << "------XAlgorithm::NewFile------" << endl;

   fIsFileOwner = kFALSE;

// If no name is given, do not create file
   if (!name || (strcmp(name, "") == 0)) return 0;

// Add extension to file name
   TString filename = gSystem->BaseName(name);
//   TString dirname  = gSystem->DirName(name);
   TString dirname  = Name2Path(name, sSEP);

   filename = Path2Name(filename, "", ".");
   filename = dirname + dSEP + filename;
   filename = filename + "_" + TString(exten) + ".root";
   if (strcmp(dirname.Data(), "") == 0) {
      dirname = gSystem->WorkingDirectory();
   }//if

// For "tmp.root" create temporary file (can be overwritten)
   TString tmp = Path2Name(name, dSEP, ".root");
   tmp = Path2Name(tmp.Data(), dSEP, "_"); // if e.g. "tmp_abc.root"
   tmp.ToLower();
   if (strcmp(tmp.Data(), "tmp") == 0) {
      fFile = new TFile(filename, "RECREATE", dirname);

      if (!fFile || fFile->IsZombie()) {
         cerr << "Error: Could not create file <" << filename.Data() << ">" << endl;
         SafeDelete(fFile);
         return 0;
      } else if (fFile->IsOpen()) {
         if (XManager::fgVerbose) {
            cout << "Creating new temporary file <" << filename.Data() << "> for <"
                 << GetName() << ">..." << endl;
         }//if
         fIsFileOwner = kTRUE;
         return fFile;
      }//if
   }//if

// Create new root file
   const char *fname;
   if ((fname = gSystem->ExpandPathName(filename.Data()))) {
      fFile = gROOT->GetFile(fname);
      if (fFile) {
         cerr << "Error: File <" << filename.Data() << "> does already exist" << endl;
         delete [] (char*)fname;
         return 0;
      }//if

      if (gSystem->AccessPathName(filename.Data())) {
         fFile = new TFile(filename, "CREATE", dirname);
      } else {
         cerr << "Error: File <" << filename.Data() << "> does already exist" << endl;
         delete [] (char*)fname;
         return 0;
      }//if

      if (!fFile || fFile->IsZombie()) {
         /*fail*/;
      } else if (fFile->IsOpen()) {
         if (XManager::fgVerbose) {
            cout << "Creating new file <" << filename.Data() << "> for <"
                 << GetName() << ">..." << endl;
         }//if
         delete [] (char*)fname;
         fIsFileOwner = kTRUE;
         return fFile;
      }//if

      delete [] (char*)fname;
   }//if

   cerr << "Error: Could not create file <" << filename.Data() << ">" << endl;
   SafeDelete(fFile);

   return 0;
}//NewFile

//______________________________________________________________________________
Int_t XAlgorithm::InitParameters(Int_t npar, Double_t *pars)
{
   // Initialize parameters for algorithm
   if(kCSa) cout << "------XAlgorithm::InitParameters------" << endl;

//??   if (fPars && fNPar != npar) {delete [] fPars; fPars = 0;}

   if (npar == 0) return errNoErr;
   if (pars == 0) return errNoErr;

   fNPar = npar;
   if (!fPars) fPars = new (nothrow) Double_t[fNPar];
   if (!fPars) return errInitMemory;
   memcpy(fPars, pars, npar*sizeof(Double_t));

   return errNoErr;
}//InitParameters

//______________________________________________________________________________
void XAlgorithm::DeleteParameters()
{
   // Delete parameters
   if(kCSa) cout << "------XAlgorithm::DeleteParameters------" << endl;

   if (fPars) {delete [] fPars; fPars = 0;}
   fNPar = 0;
}//DeleteParameters

//______________________________________________________________________________
Option_t *XAlgorithm::GetOptions(const char *sep)
{
   // Extract substring before sep
   if(kCSa) cout << "------XAlgorithm::GetOptions------" << endl;

   TString option = SubString(fOption.Data(), sep, 0);

   return (Option_t*)option.Data();
}//GetOptions

//______________________________________________________________________________
Option_t *XAlgorithm::GetOptions(TString &suboption, const char *sep)
{
   // Return substring before sep as option
   // Pass substring after sep to parameter suboption
   // Note: suboption must be created/deleted in calling function as:
   //       TString *subopt = new TString[1]; delete [] subopt
   if(kCSa) cout << "------XAlgorithm::GetOptions------" << endl;

   TString option;

   option    = SubString(fOption.Data(), sep, 0);
   suboption = SubString(fOption.Data(), sep, 1);
   if (strcmp(suboption.Data(), option.Data()) == 0) {
      suboption = "";
   }//if

   return (Option_t*)option.Data();
}//GetOptions

//______________________________________________________________________________
Double_t *XAlgorithm::Array2Log(Int_t n, Double_t *x, Double_t neglog, const char *base)
{
   // Convert array x to logarithm of base and return converted x
   // Negative values will not be converted but set to neglog
   if(kCSa) cout << "------XAlgorithm::Array2Log------" << endl;

   if (n == 0 || x == 0) return 0;

   if (strcmp(base, "0") == 0) {
      return x;
   } else if (strcmp(base, "log2") == 0) {
      for (Int_t i=0; i<n; i++) { 
         x[i] = (x[i] > 0) ? TMath::Log2(x[i]) : neglog;
      }//for_i
   } else if (strcmp(base, "log10") == 0) {
      for (Int_t i=0; i<n; i++) { 
         x[i] = (x[i] > 0) ? TMath::Log10(x[i]) : neglog;
      }//for_i
   } else if (strcmp(base, "log") == 0) {
      for (Int_t i=0; i<n; i++) { 
         x[i] = (x[i] > 0) ? TMath::Log(x[i]) : neglog;
      }//for_i
   } else {
      cout << "Warning: LogBase <" << base
           << "> is not known, using LogBase = 0." << endl;
      base = "0";
   }//if

   return x;
}//Array2Log

//______________________________________________________________________________
Double_t *XAlgorithm::Array2Pow(Int_t n, Double_t *x, const char *base)
{
   // Convert array x from logarithm of base and return converted  x
   if(kCSa) cout << "------XAlgorithm::Array2Pow------" << endl;

   if (n == 0 || x == 0) return 0;

   if (strcmp(base, "0") == 0) {
      return x;
   } else if (strcmp(base, "log2") == 0) {
      for (Int_t i=0; i<n; i++) { 
         x[i] = TMath::Power(2, x[i]);
      }//for_i
   } else if (strcmp(base, "log10") == 0) {
      for (Int_t i=0; i<n; i++) { 
         x[i] = TMath::Power(10, x[i]);
      }//for_i
   } else if (strcmp(base, "log") == 0) {
      for (Int_t i=0; i<n; i++) { 
         x[i] = TMath::Power(TMath::E(), x[i]);
      }//for_i
   }//if

   return x;
}//Array2Pow

//______________________________________________________________________________
Int_t XAlgorithm::TestNumParameters(Int_t npar)
{
   // Internal method to test if fNPar is (at least) equal to npar
   if(kCSa) cout << "------XAlgorithm::TestNumParameters------" << endl;

   if (fNPar < npar) {
      cerr << "Error: At least <" << npar 
           << ">parameters are neeeded for algorithm of type <" << fTitle.Data()
           << ">." << endl;
      return errInitParameters;
   }//if

   return errNoErr;
}//TestNumParameters


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XManager                                                             //
//                                                                      //
// Base class for managing files                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XManager::XManager()
         :TNamed()
{
   // Default Manager constructor
   if(kCS) cout << "---XManager::XManager(default)------" << endl;

   fFile        = 0;
   fSetting     = 0;
   fContent     = 0;
   fTreeSet     = 0;
   fPlotter     = 0;
   fTrash       = 0;
   fIsFileOwner = kFALSE;
   fAbort       = kFALSE;
   fInterrupt   = kFALSE;
   fInitFlag    = kFALSE;
}//Constructor

//______________________________________________________________________________
XManager::XManager(const char *name, const char *title, Int_t verbose)
         :TNamed(name, title)
{
   // Normal Manager constructor
   if(kCS) cout << "---XManager::XManager------" << endl;

   fFile        = 0;
   fSetting     = 0;
   fContent     = 0;
   fTreeSet     = 0;
   fPlotter     = 0;
   fTrash       = new TList();
   fIsFileOwner = kFALSE;
   fAbort       = kFALSE;
   fInterrupt   = kFALSE;
   fInitFlag    = kFALSE;
   fgVerbose    = verbose;
}//Constructor

//______________________________________________________________________________
XManager::~XManager()
{
   // Manager destructor
   if(kCS) cout << "---XManager::~XManager------" << endl;

   this->Close();

   SafeDelete(fPlotter);
   SafeDelete(fSetting);
   if(fTrash) {fTrash->Delete(); delete fTrash; fTrash = 0;}
}//Destructor

//______________________________________________________________________________
Int_t XManager::Initialize(const char *type, const char *data, Option_t *option,
                const char *infile, Bool_t hasPlotter)
{
   // Initialize default types: 
   // type:   source type, e.g. GeneChip or GenePix
   // data:   data type, e.g. "Data" (for *.CEL) or "Schemes" (for *.CDF)
   // option: option, e.g. "R" for R memory management
   // If infile is given, import default arguments from infile
   // Note: Parameters defined in subclasses of XManager should be
   //       defined or reset in InitDefaults() or ImportDefaults()
   if(kCS) cout << "------XManager::Initialize------" << endl;

// Reset
   fAbort     = kFALSE;
   fInterrupt = kFALSE;
   fInitFlag  = kFALSE; //?? reset here???

// Check if type is present
   if (strcmp(type, "") != 0) {
      SetTitle(type);
   } else if (strcmp(GetTitle(), "") != 0) {
      type = GetTitle();
   } else {
      cerr << "Error: Type is not initialized." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

   fDataType = TString(data);
   fOption   = TString(option);

// Initialize setting
   SafeDelete(fSetting);
   fSetting = this->NewSetting(type, infile);
   if (!fSetting) {
      cout << "Error: Could not initialize setting." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

// Initialize defaults
   Int_t err = errNoErr;
   if (strcmp(infile,"") != 0) {
      err = this->ImportDefaults(infile);
      if (err != errNoErr) {
         cerr << "Error: Could not import defaults." << endl;
         fAbort = kTRUE;
         return errAbort;
      }//if
   } else {
      err = this->InitDefaults();
      if (err != errNoErr) {
         cerr << "Error: Could not initialize default settings." << endl;
         fAbort = kTRUE;
         return errAbort;
      }//if
   }//if

// Initialize plotter
   if (hasPlotter) {
      SafeDelete(fPlotter);
      fPlotter = this->NewPlotter("Plotter", type);
      if (!fPlotter) {
         cerr << "Error: Could not initialize plotter." << endl;
         fAbort = kTRUE;
         return errAbort;
      }//if
   }//if

   return errNoErr;
}//Initialize

//______________________________________________________________________________
Int_t XManager::InitNA(Double_t na)
{
   // Initialize value used to indicate missing values
   // Note: must be called immediately after method Initialize() or New()
   if(kCS) cout << "------XManager::InitNA------" << endl;

   if (!fSetting) return this->HandleError(errInitSetting);

   fSetting->InitNA(na);

   return errNoErr;
}//InitNA

//______________________________________________________________________________
Int_t XManager::InitAlgorithm(const char *name, const char *type,
                Option_t *options, const char *filename, Int_t npars,
                Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5,
                Double_t p6, Double_t p7, Double_t p8, Double_t p9, Double_t p10)
//                Int_t npars, ...)
{
   // Initialize algorithm with name and type and with options,
   // and initialize the corresponding npars parameters
   // For no options set options=""
   // Multiple options are separated by ":", e.g. options="opt1:opt2:opt2"
   // Note: if filename != "", then a temporary file for the algorithm with
   //       name filename will be created to store data, e.g. trees
   if(kCS) cout << "------XManager::InitAlgorithm------" << endl;

   if (fAbort) return errAbort;
   if (!fSetting) return this->HandleError(errInitSetting);

   Int_t err = errNoErr;

// Get parameters
/*
// cannot be handled by CINT on PPC MacOS X!
   Double_t *pars = new Double_t[npars];
   va_list argPtr;
   va_start(argPtr, npars);
   for (Int_t i=0; i<npars; i++) {
      pars[i] = va_arg(argPtr, Double_t);
   }//for_i
   va_end(argPtr);
*/
   if (npars > 10) {
      npars = 10;
      cout << "Warning: Maximum number of allowed parameters is ten." << endl;
   }//if
   Double_t *pars = new Double_t[npars];
   if (npars > 0) pars[0] = p1;
   if (npars > 1) pars[1] = p2;
   if (npars > 2) pars[2] = p3;
   if (npars > 3) pars[3] = p4;
   if (npars > 4) pars[4] = p5;
   if (npars > 5) pars[5] = p6;
   if (npars > 6) pars[6] = p7;
   if (npars > 7) pars[7] = p8;
   if (npars > 8) pars[8] = p9;
   if (npars > 9) pars[9] = p10;

// To lower
   TString sname = name;    sname.ToLower();
   TString stype = type;    stype.ToLower();
   TString sopts = options; sopts.ToLower();

// Initialize algorithm
   if (!fSetting) {
      cout << "Error: Setting is not initialized." << endl;
      return errAbort;
   }//if
   err = fSetting->InitAlgorithm(sname, stype, sopts, filename, npars, pars);
   //note: need to copy Double_t *pars in: (since delete pars below!)
   if (err) {
      cerr << "Error: Algorithm of type <" << type << "> is not known." << endl;
      fAbort = kTRUE;
   }//if

   if (pars) {delete [] pars; pars = 0;}
   return err;
}//InitAlgorithm

//______________________________________________________________________________
void XManager::ResetAlgorithm(const char *name, const char *type)
{
   // Reset algorithm with name and optional with type
   if(kCS) cout << "------XManager::ResetAlgorithm------" << endl;

// To lower
   TString sname = name; sname.ToLower();
   TString stype = type; stype.ToLower();

   fSetting->ResetAlgorithm(sname.Data(), stype.Data());
}//ResetAlgorithm

//______________________________________________________________________________
Int_t XManager::New(const char *name, const char *dir, const char *type,
                Option_t *data, const char *option)
{
   // Create new root file with name in directory dir for type and kind
   // type:   source type, e.g. GeneChip or GenePix
   // data:   data type, e.g. "Data" (for *.CEL) or "Schemes" (for *.CDF)
   // option: option, e.g. "R" for R memory management
   // Note: for name is "tmp" or "tmp_abc", i.e. "tmp.root" or "tmp_abc.root"
   //       a temporary file will be created using RECREATE
   if(kCS) cout << "------XManager::New------" << endl;

//?? check first that manger->Initialize() is called first!!!??
   if (fAbort) return errAbort;
   fAbort = kTRUE;

// Get current directory
   if (strcmp(dir, "") == 0) {
      dir = gSystem->WorkingDirectory();
      if (XManager::fgVerbose) {
         cout << "Note: No directory given to store root file:"  << endl;
         cout << "      Using working directory <" << dir << ">" << endl;
      }//if
   }//if

// New file
   TString fullname = FullName(dir, name, dSEP) + ".root";
   TString dirname  = Name2Path(fullname, sSEP);
   fFile = this->NewFile(fullname.Data(), dirname.Data());
   if (!fFile) return errCreateFile;
   fIsFileOwner = kTRUE;

// Change dir to file
   fFile->cd();

// Create file content list
   SetTitle(type);
   fDataType = (strcmp(fDataType.Data(),"") == 0) ? TString(data)   : fDataType;
   fOption   = (strcmp(fOption.Data(),  "") == 0) ? TString(option) : fOption;

   fContent = this->NewContent(kContent, data, type);
   if (!fContent) {
      cerr << "Error: Could not create content list for <" << name << ">" << endl;
      return errAbort;
   }//if

// Initialize default setting (if not already initialized)
   if (!fSetting) {
      fSetting = this->NewSetting(type, "");
      if (!fSetting) {
         cout << "Error: Could not initialize setting." << endl;
         return errAbort;
      }//if
   }//if

   fAbort    = kFALSE;
   fInitFlag = kTRUE;
   return errNoErr;
}//New

//______________________________________________________________________________
Int_t XManager::Open(const char *fullname, const char *data, Option_t *option,
                Option_t *fileopt)
{
   // Open root file with fullname "/dir/subdir/name.root" in fileopt mode
   // (default fileopt is "READ")
   // data: data type, e.g. "Data" (for *.CEL) or "Schemes" (for *.CDF)
   if(kCS) cout << "------XManager::Open------" << endl;

   return this->Update(fullname, data, option, "", "", fileopt);
}//Open

//______________________________________________________________________________
Int_t XManager::Update(const char *fullname, const char *data, Option_t *option,
                const char *userID, const char *password, Option_t *fileopt)
{
   // Open root file with fullname "/dir/subdir/name.root" in fileopt mode
   // (default fileopt is "UPDATE")
   // data:   data type, e.g. "Data" (for *.CEL) or "Schemes" (for *.CDF)
   // option: option, e.g. "R" for R memory management
   // Note: fullname must always end with ".root"
   if(kCS) cout << "------XManager::Update------" << endl;

// NOT YET USED - to prevent compiler warnings:
   userID = 0; password = 0;

   if (fAbort) return errAbort;

//   fDataType = TString(data);
   fDataType = (strcmp(fDataType.Data(),"") == 0) ? TString(data) : fDataType;
   fOption   = (strcmp(fOption.Data(),  "") == 0) ? TString(option) : fOption;

   // convert fileopt toupper
   TString opt = Path2Name(fileopt, dSEP, ".");
   opt.ToUpper();

// Check userID and password
   if (strcmp(opt.Data(), "UPDATE") == 0) {
//to do: see XPSManager::FUpdateDB
   }//if

// Open file
   Bool_t isOwner = kFALSE;
   fFile = this->OpenFile(fullname, opt.Data(), isOwner);
   if (!fFile) {fAbort = kTRUE; return errCreateFile;}
   // assure that manager remains owner if it calls OpenFile multiple times
   if (!fIsFileOwner) fIsFileOwner = isOwner;

// Change dir to file
   fFile->cd();

// Get Content
   fFile->GetObject(kContent, fContent);
   if (!fContent) {
      cerr << "Error: File index <" << kContent << "> is missing" << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

// Initialize default setting with type (if not already initialized)
   if (!fSetting) {
      fSetting = this->NewSetting(GetTitle(), "");
      if (!fSetting) {
         cout << "Error: Could not initialize setting." << endl;
         fAbort = kTRUE;
         return errAbort;
      }//if
   }//if

   fInitFlag = kTRUE;
   return errNoErr;
}//Update

//______________________________________________________________________________
void XManager::Close(Option_t *option)
{
   // Close root file
   if(kCS) cout << "------XManager::Close------" << endl;

   if (fFile) {
      if (fFile->IsWritable() && this->Save() == kFALSE) {
         cerr << "Could not save content to file <" << fFile->GetName() << ">."
              << endl;
      }//if

// NOT YET USED - to prevent compiler warnings:
   option = 0;
/*////////////////////
//?? TO DO: BETTER:
      TString opt = option;
      opt.ToLower();
      if (opt.Contains("save")) this->Save();
*/////////////////////

      SafeDelete(fContent);
      fTreeSet = 0;  //deleted in fContent

//      if (fIsFileOwner) {SafeDelete(fFile);}
//      else              {fFile = 0;}  //??
      if (fIsFileOwner) {
         fFile->Close("R");
         SafeDelete(fFile);
      } else {
         fFile->Close("R");
         fFile = 0;
      } //if
   }//if

   fInitFlag = kFALSE;
}//Close

//______________________________________________________________________________
Bool_t XManager::Save()
{
   // Save root file
   if(kCS) cout << "------XManager::Save------" << endl;

   Bool_t isSaved = kFALSE;

   if (fFile) {
      // check if file is writable
      if (!fFile->IsWritable()) return kFALSE;

      // write content to file only if new file or file is updated
      fFile->cd();
      if (fContent) {
         fContent->Write("", TObject::kWriteDelete);
      }//if

      // write file to disk
      fFile->Flush();

      isSaved = kTRUE;
   }//if

   return isSaved;
}//Save

//______________________________________________________________________________
Int_t XManager::AddTree(const char *setname, const char *intree,
                Int_t treeid, Option_t *option)
{
   // Add intree to treeset setname.
   // Argument intree can be: /path/filename.root/treeset/treename.exten
   // Note: If no root filename is given, the current gDirectory is used as file
   // For treename is "*.exten", all trees with exten are added from root file
   // Parameters treeid and option are optional parameters, which can be used
   // to group trees and/or set options.
   // Note: If tree headers should not be added to list fHeaders of class
   //       XTreeSet the derived class need to overwrite method AddTreeHeader()
   //       and e.g. do: if (treeid > 0) do something else AddTreeHeader()
   if(kCS) cout << "------XManager::AddTree------" << endl;

   if (fAbort || !fFile || !fContent) return errAbort;

   Int_t   err = errNoErr;
   TString opt = "";
   if (strcmp(option, "") != 0) {opt = option; opt.ToLower();}

/////////////
//?? PROBLEM: memory leak if fTreeSet with different setname exists??
// need to delete in fContent ?? only in case of export where fFile is READ only
////////////..

// Get treeset
   fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   if (!fTreeSet) {
      // Create new treeset of type
      fTreeSet = this->NewTreeSet(GetTitle());
      if (!fTreeSet) {
         cerr << "Error: Could not create tree set." << endl;
         fAbort = kTRUE;
         return errAbort;
      }//if

      // Init treeset 
      err = fTreeSet->Initialize(fFile, fSetting, "", "");
      if (err != errNoErr) {
         cerr << "Error: Could not initialize tree set." << endl;
         fAbort = kTRUE;
         return errAbort;
      }//if

      fTreeSet->SetName(setname);
      fTreeSet->SetManager(this);
      fContent->Add(fTreeSet);
   }//if

///////////////////////////////
// TO DO: recognize "./" as current working directory!!!
///////////////////////////////

// Extract tree name from intree
   TString inname = Path2Name(intree, dSEP, "");
   if (strstr(inname.Data(), ".root")) {
      inname = "";
   }//if
   if (strcmp(inname.Data(), "") == 0) {
      cerr << "Error: Treename for intree is missing." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

// Extract root filename from intree
   TFile * file = 0;
   TString filename = "";
   Bool_t  isOwner  = kFALSE;
   if (strstr(intree, ".root")) {
      filename = GetROOTName(intree) + ".root";
      file = this->OpenFile(filename.Data(), "READ", isOwner);
      if (!file) return perrOpenFile;
      file->cd();
   } else {
//?? not necessary?
      filename = gDirectory->GetName();
   }//if

   TDirectory *savedir = gDirectory;

// Get name of treeset and change directory
   TString sname  = "";
   if (strstr(intree,".root")) {
      TString substr = SubString(intree,'.',sSEP, kFALSE);
      if (substr) sname = Path2Name(substr.Data(), dSEP, "");
      if (sname.Contains("root")) sname = "";
   } else if (strstr(intree, dSEP)) {
      sname = Path2Name(intree,"", dSEP);
   }//if

   if (!gDirectory->cd(sname)) return this->HandleError(errGetDir, sname);

// Add trees to treeset
   TString name  = Path2Name(intree, dSEP, ".");
   TString exten = Path2Name(intree, ".", "");
   if (strcmp(name.Data(), "*") == 0) {
   // Loop over all trees with extension exten
      TKey *key = 0;
      TIter next(gDirectory->GetListOfKeys());
      while ((key = (TKey*)next())) {
         TString xten  = Path2Name(key->GetName(), ".", ";");
         TString kname = Path2Name(key->GetName(), "", ".");
         if (strcmp(xten.Data(), exten) == 0) {
            TTree* tree = (TTree*)gDirectory->Get(key->GetName());
            if (!tree) {
               cerr << "Error: Could not get tree <" << inname.Data() << ">." << endl;
               fAbort = kTRUE;
               return errGetTree;
            }//if

            fTreeSet->AddTree(tree);
            fTreeSet->AddTreeHeader(key->GetName(), treeid);

            // handle option
            err = fTreeSet->HandleOption(tree, opt.Data());
            if (err != errNoErr) {
               cerr << "Error: Could not handle option <" << option << ">." << endl;
               fAbort = kTRUE;
            }//if
         }//if
      }//while
   } else {
   // Add intree with name inname
      TTree* tree = (TTree*)gDirectory->Get(inname);
      if (!tree) {
         cerr << "Error: Could not get tree <" << inname.Data() << ">." << endl;
         fAbort = kTRUE;
         return errGetTree;
      }//if

      fTreeSet->AddTree(tree);
      fTreeSet->AddTreeHeader(inname.Data(), treeid);

      // handle option
      err = fTreeSet->HandleOption(tree, opt.Data());
      if (err != errNoErr) {
         cerr << "Error: Could not handle option <" << option << ">." << endl;
         fAbort = kTRUE;
      }//if
   }//if

   savedir->cd();

   return err;
}//AddTree

//______________________________________________________________________________
Int_t XManager::ExportSet(const char *setname, const char *exten, 
                const char *varlist, const char *outfile, const char *sep)
{
   // Export selected data as given in varlist, which are stored in treeset
   // setname as sep-delimited file outfile
   // Note: use only after selecting trees with exten with method AddTree()
//?? allow w/o AddTree()?? need to add all trees from known treeset?
   if(kCS) cout << "------XManager::ExportSet------" << endl;

   if (fAbort) return errAbort;

   Int_t err = errNoErr;

// Check for tree extension
   if (strcmp(exten,"") == 0) {
      cerr << "Error: Tree extension is missing." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

// Set directory to root file or to setname if it is a file directory
   if (!fFile) {fAbort = kTRUE; return errGetFile;}
   TDirectory *dir = fFile->GetDirectory(setname);
   if (dir) fFile->cd(setname);
   else     fFile->cd();

// Check for empty varlist
   if (strcmp(varlist, "") == 0) {
      varlist = "*";
      cout << "Warning: No varlist given, exporting all variables." << endl;
   }//if

// Test if outfile is of type "XML"
   TString xml = Path2Name(outfile,".","");
   xml.ToUpper();
   Bool_t asXML = (strcmp(xml.Data(), "XML") == 0);

// Set outname in case no outfile name is given
   TString outname = TString(outfile);
   if (strcmp(outname.Data(), "") == 0) {
      char tmp[kBuf4096]; 
      strcpy(tmp, setname);
      strcat(tmp, "_"); strcat(tmp, exten);
      outname = tmp;
   }//if

   // test outfile for exten: if no exten then add .txt or .csv.  
   TString basename = (strcmp(outfile, "") == 0) ? "" : gSystem->BaseName(outfile);
   if (strstr(basename.Data(),".") == 0) {
      if ((strcmp(sep, ",") == 0) ||
          (strcmp(sep, ";") == 0)) outname += ".csv";
      else outname += ".txt";
   }//if

// Find tree set to export tree
   fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   if (fTreeSet) {
      fTreeSet->AsXML(asXML);
      fTreeSet->SetManager(this);

      ofstream output(outname, ios::out);
      if (!output) {
         cerr << "Error: Could not create output <" << outname.Data() << ">" << endl;
         return errOpenOutput;
      }//if

      if (XManager::fgVerbose) {
         cout << "Exporting data from treeset <" << setname
              << "> to file <" << outname.Data() << ">..." << endl;
      }//if

      if (!err) err = fTreeSet->Initialize(fFile, fSetting, "", "");
      if (!err) err = fTreeSet->Export(exten, varlist, output, sep);

      output.close();
      return this->HandleError(err, setname);
   } else {
      cerr << "Error: Tree set <" << setname 
           << "> could not be found in file content" << endl;
      err = errGetTreeSet;
   }//if

   return err;
}//ExportSet

//______________________________________________________________________________
Int_t XManager::Export(const char *treename, const char *varlist, 
                const char *outfile, const char *sep)
{
   // Export selected data as given in varlist, which are stored in tree
   // treename as sep-delimited file outfile. Argument treename can be:
   // "setname.treename.exten" - treename from set setname with exten is exported
   // "setname.*.exten" - all trees from set setname with exten are exported
   if(kCS) cout << "------XManager::Export------" << endl;

   if (fAbort) return errAbort;

   Int_t err = errNoErr;

// Extract tree extension
   TString tname = Path2Name(treename, dSEP, ";");
   TString exten = Path2Name(tname.Data(), ".", "");
   if ((strcmp(exten.Data(), "") == 0) || (strcmp(exten.Data(), "root") == 0)) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return errAbort;
   }//if

// Extract set name and tree name
   TString setname = "";
   TString trename = "";
   Int_t   numsep  = NumSeparators(tname.Data(), ".");
   if (numsep == 2) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 1);
   } else if (numsep == 1) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 0);
   } else if (numsep == 0) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return errGetTree;
   }//if

// Extract root filename
   TString filename = "";
   if (strstr(treename, ".root")) {
      filename = GetROOTName(treename) + ".root";
      this->Open(filename.Data());
   }//if
   if (!fFile) {fAbort = kTRUE; return errGetFile;}

// Change directory
   if (!fFile->cd(setname)) return this->HandleError(errGetDir, setname);

// Check for empty varlist
   if (strcmp(varlist, "") == 0) {
      varlist = "*";
      cout << "Warning: No varlist given, exporting all variables." << endl;
   }//if

// Test if outfile is of type "XML"
   TString xml = Path2Name(outfile, ".", "");
   xml.ToUpper();
   Bool_t asXML = (strcmp(xml.Data(), "XML") == 0);

// Set outname in case no outfile name is given
   TString outname = TString(outfile);
   if (strcmp(outname.Data(), "") == 0) {
      char tmp[kBuf4096]; 
      if (strcmp(trename.Data(), "*") == 0) {
         if ((strcmp(setname.Data(), "*") == 0) ||
             (strcmp(setname.Data(), "")  == 0)) strcpy(tmp, "AllTrees");
         else strcpy(tmp, setname.Data());
      } else {
         strcpy(tmp, trename.Data());
      }//if
      strcat(tmp, "_"); strcat(tmp, exten.Data());
      outname = tmp;
   }//if

   // test outfile for exten: if no exten then add .txt or .csv.  
   TString basename = (strcmp(outfile, "") == 0) ? "" : gSystem->BaseName(outfile);
   if (strstr(basename.Data(),".") == 0) {
      if ((strcmp(sep, ",") == 0) ||
          (strcmp(sep, ";") == 0)) outname += ".csv";
      else outname += ".txt";
   }//if

// If treename is "*.exten", export all trees with exten
//??  not necessary? only strcmp(trename.Data(), "*") == 0) ??
   if ((strcmp(setname.Data(), trename.Data()) == 0) &&
       (strcmp(trename.Data(), "*") == 0)) {
      if ((strcmp(exten.Data(), "") == 0) || (strcmp(exten.Data(), "*") == 0)) {
         cerr << "Error: Missing extension, name should be <*.exten>" << endl;
         return errAbort;
      }//if

      ofstream output(outname, ios::out);
      if (!output) {
         cerr << "Error: Could not create output <" << outname.Data() << ">" << endl;
         return errOpenOutput;
      }//if

      if (XManager::fgVerbose) {
         cout << "Exporting data from all trees with extension <" << exten
              << "> to file <" << outname.Data() << ">..." << endl;
      }//if
      err = this->ExportTrees(exten.Data(), varlist, output, asXML, sep);

      output.close();
      return this->HandleError(err, treename);
   }//if

// Else find tree set to export tree
   // find treeset with setname
   fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   // else find treeset containing tree with trename
   if (!fTreeSet) {
      trename += "." + exten;  //treename.exten
      fTreeSet = (XTreeSet*)fContent->FindObject(trename, "XTreeSet");
   }//if
   if (!fTreeSet) {
      cerr << "Error: Tree set <" << setname.Data() 
           << "> could not be found in file content" << endl;
      return errGetTreeSet;
   }//if

   fTreeSet->AsXML(asXML);
   fTreeSet->SetManager(this);

   err = fTreeSet->Initialize(fFile, fSetting, "", "");
   if (err) return err;

// Change directory
   if (!fFile->cd(fTreeSet->GetName())) return this->HandleError(errGetDir, fTreeSet->GetName());

   ofstream output(outname, ios::out);
   if (!output) {
      cerr << "Error: Could not create output <" << outname.Data() << ">" << endl;
      return errOpenOutput;
   }//if

   if (XManager::fgVerbose) {
      cout << "Exporting data from tree <" << trename.Data()
           << "> to file <" << outname.Data() << ">..." << endl;
   }//if
   if (strcmp(trename.Data(), "*") == 0) {
      // if treename is "*.exten", export all trees with exten
      err = fTreeSet->ExportTrees(exten.Data(), varlist, output, sep);
   } else {
      trename += "." + exten;  //treename.exten
      err = fTreeSet->ExportTree(exten.Data(), 1, &trename, varlist, output, sep);
   }//if

   output.close();

   return this->HandleError(err, trename);
}//Export

//______________________________________________________________________________
Int_t XManager::Import(const char *setname, const char *infile, const char *treename,
                Option_t *option, const char *sep, char delim, Int_t split)
{
   // Import data belonging to dataset setname from file infile. Option should
   // be of type "option.exten" where exten describes the type of tree.
   // If option = "UPDATE" replace existing tree with new infile.
   // Note: If setname = "", infile name is used as setname
   // Note: If treename = "", infile name is used as treename
   // Note: If infile = "path/*.exten" then all infiles in "path" with 
   //       extension "exten" will be imported (treename will be infile name)
   // Note: infile of type XML needs to have extension "xml", i.e. "infile.xml", 
   if(kCS) cout << "------XManager::Import------" << endl;

   Int_t err = errNoErr;

   fInterrupt = kFALSE;
   if (fAbort) {fInterrupt = kTRUE; return errGeneral;}

// Test for fFile and cd
   if (!fFile) return this->HandleError(errGetFile, "*.root");
   fFile->cd();

// Disect infile
   TString name = Path2Name(infile, dSEP, ".");
   TString xten = Path2Name(infile, ".", "");

// Append directory to infile
   TString iname = TString(infile);
   if (!iname.Contains(dSEP)) {
      iname  = TString(gSystem->WorkingDirectory()) + TString(dSEP) + iname;
   }//if

   const char *fullname = gSystem->ExpandPathName(iname.Data());
   TString path = Name2Path(fullname, sSEP);

// Change system directory to path
   TString savedir = gSystem->WorkingDirectory();
   if (!gSystem->ChangeDirectory(path.Data())) {
      cerr << "Error: Path <" << path.Data() << "> is not known." << endl;
      return errGeneral;
   }//if

   void *dirp = 0;
   if ((dirp = gSystem->OpenDirectory(".")) == 0) {
      gSystem->ChangeDirectory(savedir.Data());
      return errGeneral;
   }//if

// Add infile(s) to list infiles
   TList *infiles = new TList();
   TObjString *objstr = 0;
   if (strcmp(name.Data(), "*") == 0) {
      const char *entry;
      while ((entry = gSystem->GetDirEntry(dirp)) != 0) {
         if (strcmp(entry, ".") && strcmp(entry, "..") &&
             (strcmp(Path2Name(entry, ".", "").Data(), xten.Data()) == 0)) {
            iname = path + TString(dSEP) + TString(entry);
            objstr = new TObjString(iname);
            infiles->Add(objstr);
         }//if
      }//while
   } else {
      objstr = new TObjString(fullname);
      infiles->Add(objstr);
   }//if

// Change to old system directory (to prevent <TFile::GetSize>: "cannot stat the file")
   gSystem->FreeDirectory(dirp);
   gSystem->ChangeDirectory(savedir.Data());

// Test if infile is of type "XML"
   xten.ToUpper();
   Bool_t isXML = (strcmp(xten.Data(), "XML") == 0);

// Set default setname to first infile
   if (strcmp(setname, "") == 0) {
      objstr  = (TObjString*)(infiles->At(0));
//?? memerr?      setname = objstr->GetName();
      setname = strcpy((char*)setname, objstr->GetName());
   }//if

// Append extension to option
   TString opt   = Path2Name(option, "", ".");
   TString exten = Path2Name(option, ".", "");
   opt.ToUpper();
   if (strcmp(exten.Data(), "") == 0) {
      cout << "Warning: No extension given for option, using default .def."
           << endl;
      exten = "def";
   }//if
   TString optext = opt + "." + exten;

// Test infile for exten: if sep="" use sep for .txt or .csv.  
   if (strcmp(sep,"") == 0) {
      if (strstr(infile,".csv") != 0) sep = ",";
      else                            sep = "\t";
   }//if

// Import infile(s)
   TIter next(infiles);
   while ((objstr = (TObjString*)next())) {
      TString objname = Path2Name(objstr->GetName(), dSEP, ".");

   // Set objname to treename if infile is not "*"
      if ((strcmp(treename, "") != 0) && (infiles->GetSize() == 1)) {
         objname = TString(treename);
      }//if

   // Open infile containing data
      // MS VC++ requires mode "ios::binary" for binary files
      ifstream input(objstr->GetName(), ios::in | ios::binary);
      if (!input) {
         cerr << "Error: File <" << objstr->GetName() << "> does not exist."
              << endl;
         fInterrupt = kTRUE;
         return errGeneral;
      }//if

   // Import data for existing or newly created tree set
      if (XManager::fgVerbose) {
         cout << "Importing <" << objstr->GetName() << "> as <"
              << (objname + "." + exten).Data() << ">..." << endl;
      }//if

      fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
      if (fTreeSet) {
         fTreeSet->SetManager(this);

         // Create directory with name of treeset and cd to directory
         TDirectory *dir = fFile->GetDirectory(setname);
         if (!dir) fFile->mkdir(setname, fDataType);
         fFile->cd(setname);

         // Update treeset
         if (strcmp(opt.Data(), "UPDATE") == 0) {
            // delete existing tree(s) with "name.exten;*"
            TString tname = objname + "." + exten + ";" + "*";
            err = fTreeSet->DeleteTree(objname, exten, "*");
            if (err != errNoErr) {
               cerr << "Error: Could not delete <" << err
                    << "> number of trees with name <" << tname.Data() << ">."
                    << endl;
               fInterrupt = kTRUE;
               input.close();
               return err;
            }//if
//??            fContent->Remove(fTreeSet);

            // update content
            err = fTreeSet->Initialize(fFile, fSetting, objstr->GetName(), objname);
            if (isXML) {
               err = fTreeSet->ImportXML(input, optext.Data(), sep, delim, split);
            } else {
               err = fTreeSet->Import(input, optext.Data(), sep, delim, split);
            }//if
            if (err == errNoErr) {
               fContent->Add(fTreeSet);
               if (XManager::fgVerbose) {
                  cout << "Existing dataset <" << setname << "> is updated..." << endl;
               }//if
            }//if
         } else if ((strcmp(opt.Data(), "CREATE") == 0) ||
                    (strcmp(opt.Data(), "NEW")    == 0) ) {
            TString tname = fTreeSet->FindTree(objname + "." + exten);
            if (strcmp(tname.Data(), "") != 0) {
               // Interrupt if tree for infile does already exist
               cerr << "Error: Data for <" << tname.Data() << "> exist already." << endl;
               fInterrupt = kTRUE;
               input.close();
               return errGeneral;
            }//if

         // Import infile 
            err = fTreeSet->Initialize(fFile, fSetting, objstr->GetName(), objname);
            if (isXML) {
               err = fTreeSet->ImportXML(input, optext.Data(), sep, delim, split);
            } else {
               err = fTreeSet->Import(input, optext.Data(), sep, delim, split);
            }//if
         }//if
      } else if ((strcmp(opt.Data(), "CREATE") == 0) ||
                 (strcmp(opt.Data(), "NEW")    == 0) ) {
//?? SafeDelete(fTreeSet);
//but: see XManager::Close(): fTreeSet deleted in fContent!

      // Create new treeset of type
         fTreeSet = this->NewTreeSet(GetTitle());
         if (!fTreeSet) {
            fAbort = kTRUE;
            input.close();
            return errAbort;
         }//if
         fTreeSet->SetName(setname);
         fTreeSet->SetManager(this);

      // Create directory with name of treeset and cd to directory
         fFile->mkdir(setname, fDataType);
         fFile->cd(setname);

      // Import new infile 
         err = fTreeSet->Initialize(fFile, fSetting, objstr->GetName(), objname);
         if (isXML) {
            err = fTreeSet->ImportXML(input, optext.Data(), sep, delim, split);
         } else {
            err = fTreeSet->Import(input, optext.Data(), sep, delim, split);
         }//if
         if (err == errNoErr) {
            fContent->Add(fTreeSet);
            if (XManager::fgVerbose) {
               cout << "New dataset <" << setname << "> is added to Content..."
                    << endl;
            }//if
         }//if
      }//if

      input.close();
   }//while

   err = this->HandleError(err, infile);

   delete [] (char*)fullname;

   infiles->Delete();
   delete infiles;

   return err;
}//Import

//______________________________________________________________________________
void XManager::Delete(const char *name)
{
   // Delete object with name
   // Delete treeset:
   //    name = "name" - delete treeset "name" and all trees in treeset
   //    name = "*"    - delete all treesets and all trees
   // Delete tree:
   //    name = "tree.exten;cycle" - delete tree with "tree.exten;cycle"
   //    name = "tree.exten"       - delete all cycles "tree.exten;*"
   // Delete tree in specified treeset:
   //    name = "treeset/tree.exten;cycle" - delete tree with "tree.exten;cycle"
   //    name = "treeset/tree.exten"       - delete all cycles "tree.exten;*"
   if(kCS) cout << "------XManager::Delete------" << endl;

   if (fAbort) {fInterrupt = kTRUE; return;}

   TString sname = Path2Name(name, "", dSEP);
   TString exten = Path2Name(name, ".", ";");
   TString cycle = Path2Name(name, ";", "");
   if (sname.Contains("."))    sname = "";
   if (strstr(name, ".") == 0) exten = "";
   if (strstr(name, ";") == 0) cycle = "";

   if (strcmp(exten.Data(), "") == 0) {
   // Delete treeset:
      // Loop over treesets
      fFile->cd();

      TIter next(fContent->GetListOfFolders());
      TObject *obj = 0;
      while ((obj = next())) {
         if (obj->InheritsFrom(XTreeSet::Class())) {
            TString setname = obj->GetName();
            if ((strcmp(sname.Data(), setname.Data()) == 0) ||
                (strcmp(sname.Data(), "*") == 0)) {
               this->DeleteTreeSet(setname);
               if (XManager::fgVerbose) {
                  cout << "Treeset <" << setname.Data() << "> has been deleted." << endl;
               }//if
            }//if
         }//if
      }//while
   } else {
   // Delete tree:
      this->DeleteTree(name);
   }//if
}//Delete

//______________________________________________________________________________
Int_t XManager::DeleteTree(const char *namecycle)
{
   // Delete tree with name "setname/treename.exten;cycle"
   if(kCS) cout << "------XManager::DeleteTree------" << endl;

   if (fAbort) {fInterrupt = kTRUE; return errAbort;}

   Int_t err = errNoErr;

   TString sname = Path2Name(namecycle, "", dSEP);
   TString tname = Path2Name(namecycle, dSEP, ".");
   TString exten = Path2Name(namecycle, ".", ";");
   TString cycle = Path2Name(namecycle, ";", "");
   if (strstr(namecycle, ".") == 0) exten = "";
   if (strstr(namecycle, ";") == 0) cycle = "";

   fTreeSet = (XTreeSet*)fContent->FindObject(sname, "XTreeSet");
   if (fTreeSet) {
      fTreeSet->SetManager(this);

      if (!fFile->cd(sname)) return this->HandleError(errGetDir, sname);

      err = fTreeSet->DeleteTree(tname, exten, cycle);
      if (err > 0) {
         cerr << "Warning: Did not delete <"  << err
              << "> trees of set <" << sname.Data() << ">."
              << endl;
         fInterrupt = kTRUE;
         return errGeneral;
      } else if (err < 0) {
         fAbort = kTRUE;
         return this->HandleError(err, sname);
      }//if

      // Remove set if all trees are deleted from file
      if (((strcmp(exten.Data(), "*")             == 0)  &&
           (strstr(namecycle,    ";")             != 0)) || 
           (fTreeSet->GetTreeHeaders()->GetSize() == 0)) {
         fContent->Remove(fTreeSet);

         this->DeleteTreeSetInfo(sname);
         this->DeleteDirectory(sname, "*");
      }//if
   } else {
      cerr << "Error: Tree set <" << sname.Data() 
           << "> could not be found in file content" << endl;
      err = errGetTreeSet;
   }//if

   return err;
}//DeleteTree

//______________________________________________________________________________
Int_t XManager::DeleteTreeSet(const char *setname)
{
   // Delete tree set with name setname
   if(kCS) cout << "------XManager::DeleteTreeSet------" << endl;

   if (fAbort) {fInterrupt = kTRUE; return errAbort;}

   Int_t     err = errNoErr;
   XTreeSet *set = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   if (set) {
      set->SetManager(this);

      if (!fFile->cd(setname)) return this->HandleError(errGetDir, setname);

      err = set->DeleteTree("*", "*", "*");
      if (err) {
         cerr << "Error: Could not delete <"  << err
              << "> trees of set <" << setname << ">." << endl;
         fInterrupt = kTRUE;
         err = errGeneral;
      }//if
      fContent->Remove(set);

      this->DeleteTreeSetInfo(setname);
      this->DeleteDirectory(setname, "*");

      SafeDelete(set);
   } else {
      cerr << "Error: Tree set <" << setname 
           << "> could not be found in file content" << endl;
      err = errGetTreeSet;
   }//if

   return err;
}//DeleteTreeSet

//______________________________________________________________________________
Int_t XManager::HandleError(Int_t err, const char *name1, const char *name2)
{
   // Handle error messages
   if(kCS) cout << "------XManager::HandleError------" << endl;

   switch (err) {
      case errNoErr:
         break;

      case errFatal:
         cerr << "Fatal Error has occured: Exiting program!" << endl;
         gSystem->Exit(1);
         break;

      case errAbort:
         cerr << "An error has occured: Need to abort current process." << endl;
         fAbort = kTRUE;
         break;

      case errGeneral:
         cerr << "Error: General error." << endl;
         fInterrupt = kTRUE;
         break;

      case errInitMemory:
         cerr << "Error: Could not initialize memory to read file <" 
              << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errCreateFile:
         cerr << "Error: Could not create file <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errCreateDir:
         cerr << "Error: Could not make directory <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errCreateTree:
         cerr << "Error: Could not create tree <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errCreateTreeSet:
         cerr << "Error: Could not create treeset <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errGetFile:
         cerr << "Error: Could not get file <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errGetDir:
         cerr << "Error: Could not get directory <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errGetTree:
         cerr << "Error: Could not get tree <" << name1 << ">." << endl;
         fAbort = kTRUE;
//??         fInterrupt = kTRUE; //instead of fAbort??
         break;

      case errWriteObject:
         cerr << "Error: Could not write object <" << name1
           << "> to directory <" << gDirectory->GetName() << ">." << endl;
         fAbort = kTRUE;
//??         fInterrupt = kTRUE; //instead of fAbort??
         break;

      case errGetTreeSet:
         cerr << "Error: Could not get treeset <" << name1 << ">." << endl;
         fAbort = kTRUE;
         break;

      case errGetTreeInfo:
         cerr << "Error: Could not get user info for tree <" << name1 << ">."
              << endl;
         fAbort = kTRUE;
         break;

      case errPrematureEOF:
         cerr << "Error: Premature end of file <" << name1 << "> reached." << endl;
         fInterrupt = kTRUE;
         break;

      case errHeaderLine:
         cerr << "Error: Header line of file <" << name1 << "> is not correct."
              << endl;
         fInterrupt = kTRUE;
         break;

      case errMissingColumn:
         cerr << "Error: Essential column of file <" << name1 << "> is missing."
              << endl;
         fInterrupt = kTRUE;
         break;

      case errMissingLine:
         cerr << "Error: A line of file <" << name1 << "> is missing."
              << endl;
         fInterrupt = kTRUE;
         break;

      case errReadingInput:
         cerr << "Error when reading file <" << name1 << ">."
              << endl;
         fInterrupt = kTRUE;
         break;

      case errOpenOutput:
         cerr << "Error: Could not create output <" << name1 << ">."
              << endl;
         fInterrupt = kTRUE;
         break;

      case errInitContent:
         cerr << "Error: Index <" << name1 << "> is not initialized." << endl;
         fAbort = kTRUE;
         break;

      case errMissingContent:
         cerr << "Error: <" << name1 << "> index <" << name2 << "> is missing."
              << endl;
         fAbort = kTRUE;
         break;

      case errInitSetting:
         cerr << "Error: Setting is not initialized." << endl;
         fAbort = kTRUE;
         break;

      case errInitTreeSet:
         cerr << "Error: Tree set is not initialized." << endl;
         fAbort = kTRUE;
         break;

      case errInitPlotter:
         cerr << "Error: Plotter is not initialized" << endl;
         fInterrupt = kTRUE;
         break;

      case errNumTreeEntries:
         cerr << "Error: Number of entries in tree <" << name1
              << "> is not equal to <" << name2 << ">" << endl;
         fAbort = kTRUE;
         break;

      case errEQTreeEntries:
         cerr << "Error: Trees <" << name1 << "> and <" << name2
              << "> have not equal number of entries." << endl;
         fAbort = kTRUE;
         break;

      case errClassTreeSet:
         cerr << "Error: Treeset <" << name1 << "> is inherited from class <"
              << name2 << ">." << endl;
         fAbort = kTRUE;
//??       fInterrupt = kTRUE;
         break;

      case errAlgorithm:
         cerr << "Error: <" << name1 << "> algorithm  <" << name2
              << "> is not known." << endl;
         fAbort = kTRUE;
         break;

      case errUnknownType:
         cerr << "Error: <" << name1 << "> type  <" << name2
              << "> is not known." << endl;
         fAbort = kTRUE;
         break;

      default:
         cerr << "An error has occured, aborting process" << endl;
         fAbort = kTRUE;
         break;
   }//switch

   return err;
}//HandleError

//______________________________________________________________________________
void XManager::NewCanvas(const char *name, const char *title, Int_t wtopx,
               Int_t wtopy, Int_t ww , Int_t wh, Int_t nx, Int_t ny)
{
   // Create new canvas
   if(kCS) cout << "------XManager::NewCanvas------" << endl;

   if (fAbort || !fPlotter) return;

   fPlotter->NewCanvas(name, title, wtopx, wtopy, ww , wh, nx, ny);

   if (fFile) fPlotter->SetFile(fFile);
}//NewCanvas

//______________________________________________________________________________
void XManager::CloseCanvas(Option_t *opt)
{
   // Close canvas
   if(kCS) cout << "------XManager::CloseCanvas------" << endl;

   if (fAbort || !fPlotter) return;

   fPlotter->CloseCanvas(opt);
}//CloseCanvas

//______________________________________________________________________________
Int_t XManager::Draw(const char *canvasname, const char *treename,
                const char *varlist, const char *logbases, const char *type,
                Option_t *opt, Double_t minX, Double_t maxX, Double_t minY,
                Double_t maxY, Double_t minZ, Double_t maxZ, const char *var2sort,
                Bool_t down)
{
   // Draw at most three variables from varlist for tree treename
   // Each variable is converted to logbase given in list logbases
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // The type can be: graph, multigraph, hist, profile
   // with the corresponding option opt (TGraph::PaintGraph, THistPainter::Paint)
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Variable var2sort determines the variable to sort for
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XManager::Draw(1)------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->Draw(canvasname, treename, varlist, logbases, type, opt,
                        minX, maxX, minY, maxY, minZ, maxZ, var2sort, down);
   return err;
}//Draw

//______________________________________________________________________________
Int_t XManager::Draw(const char *canvasname, const char *treename1,
                const char *treename2,  const char *varlist, const char *logbases,
                const char *type, Option_t *opt, Double_t minX, Double_t maxX,
                Double_t minY, Double_t maxY, Int_t sort, Bool_t down)
{
   // Draw variable from varlist for corresponding tree treename1 or treename2
   // Each variable is converted to logbase given in list logbases
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // type: graph, hist, mvaplot, profile - with corresponding option opt
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XManager::Draw(2)------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->Draw(canvasname, treename1, treename2, varlist, logbases, type,
                        opt, minX, maxX, minY, maxY, sort, down);
   return err;
}//Draw

//______________________________________________________________________________
Int_t XManager::Draw(const char *canvasname, const char *treename1,
                const char *treename2, const char *treename3, const char *varlist,
                const char *logbases, const char *type, Option_t *opt,
                Double_t minX, Double_t maxX, Double_t minY, Double_t maxY,
                Double_t minZ, Double_t maxZ, Int_t sort, Bool_t down)
{
   // Draw variable from varlist for corresponding tree treename1, 2 or 3
   // Each variable is converted to logbase given in list logbases
   // logbase can be: 0(linear), log(ln), log10, log2, e.g. "0:log10"
   // type can be: graph, multigraph, hist, profile with corresponding option opt
   // Axis range is given by min and  max:
   //   min = max = -1111: range is calculated for each axis separately (default) 
   //   min = max = 0: all axes have same range, min and max are calculated
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XManager::Draw(3)------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->Draw(canvasname, treename1, treename2, treename3, varlist, logbases,
                        type, opt, minX, maxX, minY, maxY, minZ, maxZ, sort, down);
   return err;
}//Draw

//______________________________________________________________________________
Int_t XManager::DrawImage(const char *canvasname, const char *treename,
                const char *varlist, const char *logbase, Option_t *opt,
                Double_t min, Double_t max, Option_t *orientation)
{
   // Draw image for tree containing data with (x,y)-coordinates
   // varlist must be of the form "fX:fY:fData" (or equivalent leaf names)
   // logbase: 0(linear), log(ln), log10, log2
   // orientation: U(up), D(down), L(rotate left), R(rotate right)
   //    +M: mirror image, e.g. RM (mirror of rotate right)
   if(kCS) cout << "------XManager::DrawImage------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->DrawImage(canvasname, treename, varlist, logbase, opt, 
                             min, max, orientation);
   return err;
}//DrawImage

//______________________________________________________________________________
Int_t XManager::DrawTree(const char *canvasname, const char *treename,
                const char *varexp, const char *selection, Option_t *opt)
{
   // Draw variable expressions for tree(s)
   if(kCS) cout << "------XManager::DrawTree------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->DrawTree(canvasname, treename, varexp, selection, opt);
   return err;
}//DrawTree

//______________________________________________________________________________
Int_t XManager::DrawEntries(const char *canvasname, const char *leafname,
                Int_t n, Int_t *entrylist, const char *logbase, const char *type,
                Option_t *opt, Double_t min, Double_t max, Int_t sort, Bool_t down)
{
   // Draw leaf with leafname for all trees listed in fHeaders and for n entries 
   // stored in array entrylist. 
   // logbase can be: 0(linear), log(ln), log10, log2
   // type can be: multigraph, hist(box), boxplot with corresponding option opt
   // Axis range is given by min and  max:
   //   min < max:  range is given by user 
   //   min = max:  all axes have same range, min and max are calculated
   // sort options: sort = 0   entries not sorted before drawing
   //               sort = -1  each entry is sorted individually before drawing
   //               sort = k+1 entries are sorted for entry k (first entry: k=0)
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XManager::DrawEntries------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->DrawEntries(canvasname, leafname, n, entrylist, logbase, type,
                               opt, min, max, sort, down);
   return err;
}//DrawEntries

//______________________________________________________________________________
Int_t XManager::DrawLeaves(const char *canvasname, const char *leafname,
                const char *logbase, const char *type, Option_t *opt,
                Double_t min, Double_t max, Int_t sort, Bool_t down)
{
   // Draw leaf with leafname for all trees listed in fHeaders. 
   // logbase can be: 0(linear), log(ln), log10, log2
   // type can be: multigraph, boxplot with corresponding option opt
   // Axis range is given by min and  max:
   //   min < max:  range is given by user 
   //   min = max:  all axes have same range, min and max are calculated
   // sort options: sort = 0   leafs not sorted before drawing
   //               sort = -1  each leafs is sorted individually before drawing
   //               sort = k   leafs are sorted for leaf of tree k
   // Note: If canvas is created by NewCanvas(), set canvasname = "".
   if(kCS) cout << "------XManager::DrawLeaves------" << endl;

   if (fAbort)    return errAbort;
   if (!fPlotter) return errInitPlotter;

   Int_t err = errNoErr;
   err = fPlotter->DrawLeaves(canvasname, leafname, logbase, type, opt,
                              min, max, sort, down);
   return err;
}//DrawLeaves

//______________________________________________________________________________
TTree *XManager::GetTree(const char *treename)
{
   // Get tree for tree with treename. Argument treename can be:
   // "setname.treename.exten" - treename from set setname with exten is exported
   // "filename.root/setname.treename.exten" - as above but including filename
   if(kCS) cout << "------XManager::GetTree------" << endl;

   if (fAbort) return 0;

// Extract tree extension
   TString tname = Path2Name(treename, dSEP, ";");
   TString exten = Path2Name(tname.Data(), ".", "");
   if ((strcmp(exten.Data(), "") == 0) || (strcmp(exten.Data(), "root") == 0)) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return 0;
   }//if

// Extract set name and tree name
   TString setname = "";
   TString trename = "";
   Int_t   numsep  = NumSeparators(tname.Data(), ".");
   if (numsep == 2) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 1);
   } else if (numsep == 1) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 0);
   } else if (numsep == 0) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return 0;
   }//if
   trename += "." + exten;

// Extract root filename
   TString filename = "";
   if (strstr(treename,".root")) {
      filename = GetROOTName(treename) + ".root";
      this->Open(filename.Data());
   }//if
   if (!fFile) {fAbort = kTRUE; return 0;}

   if (!fFile->cd(setname)) {
//      this->HandleError(errGetDir, setname);
      cerr << "Error: Tree set <" << setname.Data() 
           << "> could not be found in file content" << endl;
      return 0;
   }//if

   return (TTree*)gDirectory->Get(trename);
}//GetTree

//______________________________________________________________________________
XTreeHeader *XManager::GetTreeHeader(const char *treename)
{
   // Get tree header for tree with treename. Argument treename can be:
   // "setname.treename.exten" - treename from set setname with exten is exported
   // "filename.root/setname.treename.exten" - as above but including filename
   if(kCS) cout << "------XManager::GetTreeHeader------" << endl;

   if (fAbort) return 0;

// Extract tree extension
   TString tname = Path2Name(treename, dSEP, ";");
   TString exten = Path2Name(tname.Data(), ".", "");
   if ((strcmp(exten.Data(), "") == 0) || (strcmp(exten.Data(), "root") == 0)) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return 0;
   }//if

// Extract set name and tree name
   TString setname = "";
   TString trename = "";
   Int_t   numsep  = NumSeparators(tname.Data(), ".");
   if (numsep == 2) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 1);
   } else if (numsep == 1) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 0);
   } else if (numsep == 0) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return 0;
   }//if
   trename += "." + exten;

// Extract root filename
   TString filename = "";
   if (strstr(treename, ".root")) {
      filename = GetROOTName(treename) + ".root";
      this->Open(filename.Data());
   }//if
   if (!fFile) {fAbort = kTRUE; return 0;}
   fFile->cd();
//?   fFile->cd(setname);

   XTreeHeader *header = 0;
   fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   if (fTreeSet) {
      header = fTreeSet->GetTreeHeader(trename);
   } else {
      cerr << "Error: Tree set <" << setname.Data() 
           << "> could not be found in file content" << endl;
   }//if

   return header;
}//GetTreeHeader

//______________________________________________________________________________
XTreeInfo *XManager::GetTreeInfo(const char *treename)
{
   // Get tree info for tree with treename. Argument treename can be:
   // "setname.treename.exten" - treename from set setname with exten is exported
   // "filename.root/setname.treename.exten" - as above but including filename
   if(kCS) cout << "------XManager::GetTreeInfo------" << endl;

   if (fAbort) return 0;

// Extract tree extension
   TString tname = Path2Name(treename, dSEP, ";");
   TString exten = Path2Name(tname.Data(), ".", "");
   if ((strcmp(exten.Data(), "") == 0) || (strcmp(exten.Data(), "root") == 0)) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return 0;
   }//if

// Extract set name and tree name
   TString setname = "";
   TString trename = "";
   Int_t   numsep  = NumSeparators(tname.Data(), ".");
   if (numsep == 2) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 1);
   } else if (numsep == 1) {
      setname = SubString(tname.Data(), ".", 0);
      trename = SubString(tname.Data(), ".", 0);
   } else if (numsep == 0) {
      cerr << "Error: Tree name is missing." << endl;
      fAbort = kTRUE;
      return 0;
   }//if
   trename += "." + exten;

// Extract root filename
   TString filename = "";
   if (strstr(treename,".root")) {
      filename = GetROOTName(treename) + ".root";
      this->Open(filename.Data());
   }//if
   if (!fFile) {fAbort = kTRUE; return 0;}
   fFile->cd(setname);

   XTreeInfo *info = 0;
   TTree     *tree = 0;
   fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   if (fTreeSet) {
      tree = (TTree*)gDirectory->Get(trename);
      info = (tree) ? fTreeSet->GetTreeInfo(trename, tree) : 0;
   } else {
      cerr << "Error: Tree set <" << setname.Data() 
           << "> could not be found in file content" << endl;
   }//if

   return info;
}//GetTreeInfo

//______________________________________________________________________________
void XManager::PrintContents(const char *setname)
{
   // Print contents of tree set with name setname
   if(kCS) cout << "------XManager::PrintContents------" << endl;

   if (fAbort) return;
   fFile->cd();
//??   fFile->cd(setname);

   fTreeSet = (XTreeSet*)fContent->FindObject(Path2Name(setname, dSEP, "."), "XTreeSet");
   if (fTreeSet) {
      XTreeSet::SetPrintHeader(kTRUE);
      fTreeSet->PrintInfo();
   } else {
      cerr << "Error: Tree set <" << setname 
           << "> could not be found in file content" << endl;
   }//if
}//PrintContents

//______________________________________________________________________________
void XManager::PrintContents()
{
   // Print contents of file, i.e. information about the trees
   if(kCS) cout << "------XManager::PrintContents------" << endl;

   if (fAbort) return;
   fFile->cd();
//??   fFile->cd(setname); //how to get setname?

   XTreeSet::SetPrintHeader(kTRUE);

   TIter next(fContent->GetListOfFolders());
   TObject *obj = 0;
   while ((obj = next())) {
      if (obj->InheritsFrom(XTreeSet::Class())) {
         ((XTreeSet*)obj)->PrintInfo();
      }//if
   }//while
}//PrintContents

//______________________________________________________________________________
void XManager::SetMaxFileSize(Long64_t maxsize)
{
   // Set maximum tree size, where maxsize is defined as KB, i.e. maxsize=100 means 100KB
   // Default tree file size is 1.9 GB, default SetMaxFileSize() is 1.9 TB
   if(kCS) cout << "------XManager::SetMaxFileSize------" << endl;

   TTree::SetMaxTreeSize(1000*maxsize);
}//SetMaxFileSize

//______________________________________________________________________________
void XManager::SetBufSize(Int_t bufsize)
{
   // Set buffer size for certain tree baskets
   if(kCS) cout << "------XManager::SetBufSize------" << endl;

   fgBufSize = (bufsize < 100) ? 100 : bufsize;
}//SetBufSize

//______________________________________________________________________________
void XManager::SetVerbose(Int_t verbose)
{
   // Set verbose output
   if(kCS) cout << "------XManager::SetVerbose------" << endl;

   fgVerbose = verbose;
}//SetVerbose

//______________________________________________________________________________
Int_t XManager::GetBufSize(Int_t n, Int_t m)
{
   // Get buffer size for certain tree baskets
   // For default bufsize=32000 decrease bufsize for n trees in steps of m
   if(kCS) cout << "------XManager::GetBufSize------" << endl;

   if (fgBufSize == 32000) return (fgBufSize/(1 + n/m));

   return fgBufSize;
}//GetBufSize

//______________________________________________________________________________
Int_t XManager::GetVerbose()
{
   // Get verbose output
   if(kCS) cout << "------XManager::GetVerbose------" << endl;

   return fgVerbose;
}//GetVerbose

//______________________________________________________________________________
Bool_t XManager::IsOpen(TFile *file, const char *filename)
{
   // Check if root file is already open
   // Note: If new file is identical to old file, a warning will be given that the
   //       file is already open but return value is kFALSE to prevent closing file
   if(kCS) cout << "------XManager::IsOpen------" << endl;

   if (file) {
      TString oldname = file->GetName();
      TString newname = Path2Name(filename, dSEP, ".") + ".root";

      TString xpaname;
//ccc char memory problem??
      const char *fname;
      if ((fname = gSystem->ExpandPathName(filename))) {
         xpaname = TString(fname);
         delete [] (char*)fname;
      }//if

      if ((strcmp(oldname.Data(), newname.Data()) == 0) ||
          (strcmp(oldname.Data(), xpaname.Data()) == 0)) {
         cout << "Warning: File <" << oldname.Data() << "> is already open." << endl;
         return kFALSE;
      }//if

      return kTRUE;
   }//if

   return kFALSE;
}//IsOpen

//______________________________________________________________________________
XSetting *XManager::NewSetting(const char *type, const char *infile)
{
   // Create new setting
   if(kCS) cout << "------XManager::NewSetting------" << endl;

   XSetting *setting = new XSetting(type, infile);
   return setting;
}//NewSetting

//______________________________________________________________________________
XPlot *XManager::NewPlotter(const char *name, const char *title)
{
   // Create new plotter
   if(kCS) cout << "------XManager::NewPlotter------" << endl;

   XPlot *plotter = new XPlot(name, title);
   return plotter;
}//NewPlotter

//______________________________________________________________________________
Int_t XManager::ExportTrees(const char *exten, const char *varlist, 
                ofstream &output, Bool_t asXML, const char *sep)
{
   // Export variables from varlist for all trees with extension exten
   if(kCS) cout << "------XManager::ExportTrees------" << endl;

   Int_t err = errNoErr;

// Count number of trees with extension exten to init array arrName
   Int_t numtrees = 0;
   TKey *key = 0;
   TIter next(gDirectory->GetListOfKeys());
   while ((key = (TKey*)next())) {
      TString xten = Path2Name(key->GetName(),".",";");
      if (strcmp(xten.Data(), exten) == 0) {
         numtrees++;
      }//if
   }//while

   TString *arrName = 0;
   if (!(arrName = new (nothrow) TString[numtrees])) return errInitMemory;

// Create temporary treeset of type
   XTreeSet *set = this->NewTreeSet(GetTitle());
   if (!set) {delete [] arrName; return errCreateTreeSet;}
   set->SetManager(this);

// Add trees to treeset
   Int_t i = 0;
   key = 0;
   next.Reset();
   while ((key = (TKey*)next())) {
      TString xten = Path2Name(key->GetName(),".",";");
      if (strcmp(xten.Data(), exten) == 0) {
//?? PROBLEM? key name;cycle??
         arrName[i] = key->GetName();
         TTree *tree = (TTree*)gDirectory->Get(arrName[i]);
         if (!tree) {delete [] arrName; delete set; return errGetTree;}
         set->AddTree(tree);
         i++;
      }//if
   }//while

// Export all trees with exten
   if ((err = set->Initialize(fFile, fSetting, "", "")) == 0) {
      set->AsXML(asXML);
      err = set->ExportTree(exten, numtrees, arrName, varlist, output, sep);
   }//if

// Clean up
   delete [] arrName;
   SafeDelete(set);

   return err;
}//ExportTrees

//______________________________________________________________________________
Int_t XManager::HandleOption(const char *setname, const char *treename, 
                Option_t *option)
{
   // Find tree set setname to handle option for treename
   // In XManager HandleOption() is called from AddTree() only.
   if(kCS) cout << "------XManager::HandleOption------" << endl;

   fTreeSet = (XTreeSet*)fContent->FindObject(setname, "XTreeSet");
   if (fTreeSet) {
      return fTreeSet->HandleOption(treename, option);
   } else {
      cerr << "Error: Treeset <" << setname << "> not found." << endl;
   }//if

   return errAbort;
}//HandleOption

//______________________________________________________________________________
TFile *XManager::NewFile(const char *name, const char *title)
{
   // Create new root file and return pointer to TFile
   // If name is "tmp" or "tmp_abc", i.e. "tmp.root" or "tmp_abc.root"
   // a temporary file will be created using RECREATE
   // For "tmpdt" or "tmpd" or "tmpt" a temporary file containing
   // date (d) and/or time (t) will be created
   // If name starts with "td_" or "d_" or "t_" a permanent file containing
   // date (d) and/or time (t) will be created, i.e. "name_date_time.root"
   if(kCS) cout << "------XManager::NewFile------" << endl;

//?? better?? fInterrupt = kTRUE;
   fAbort = kTRUE;

   TFile *file = 0;
//to do: check name for TNetFile!!

// Get date and time
   TDatime *datime = new TDatime();
   Int_t    date   = datime->GetDate();
   Int_t    time   = datime->GetTime();
   delete datime;

// Dissect name
   TString pathname = SubString(name, '\"', sSEP, kTRUE);
   TString tmp      = Path2Name(name, dSEP, ".");
   TString filename = SubString(tmp.Data(), '_', '.', kTRUE);

   tmp = Path2Name(tmp.Data(), dSEP, "_"); // if e.g. "tmp_abc.root"
   tmp.ToLower();

// For "tmp.root" create temporary file (can be overwritten)
   if ((strcmp(tmp.Data(), "tmp")   == 0) ||
       (strcmp(tmp.Data(), "tmpd")  == 0) ||
       (strcmp(tmp.Data(), "tmpt")  == 0) ||
       (strcmp(tmp.Data(), "tmpdt") == 0)) {

      filename = pathname + dSEP + filename;
      if (strcmp(tmp.Data(), "tmpd")  == 0) {
         filename += "_"; filename += date;
      } else if (strcmp(tmp.Data(), "tmpt") == 0) {
         filename += "_"; filename += time;
      } else if (strcmp(tmp.Data(), "tmpdt") == 0) {
         filename += "_"; filename += date;
         filename += "_"; filename += time;
      }//if
      filename += ".root";

      file = new TFile(filename, "RECREATE", title);

      if (!file || file->IsZombie()) {
         cerr << "Error: Could not create file <" << name << ">" << endl;
         SafeDelete(file);
         fAbort = kTRUE;
         return 0;
      } else if (file->IsOpen()) {
         if (XManager::fgVerbose) {
            cout << "Creating new temporary file <" << filename.Data() << ">..." << endl;
         }//if
         fAbort = kFALSE;
         return file;
      }//if
   }//if

// Create new root file
   const char *fname;
   if ((fname = gSystem->ExpandPathName(name))) {
      file = gROOT->GetFile(fname);
      if (file) {
         cerr << "Error: File <" << name << "> does already exist" << endl;
         delete [] (char*)fname;
         return 0;
      }//if

      if (gSystem->AccessPathName(fname)) {
         filename = pathname + dSEP + filename;
         if (strcmp(tmp.Data(), "d")  == 0) {
            filename += "_"; filename += date; filename += ".root";
         } else if (strcmp(tmp.Data(), "t") == 0) {
            filename += "_"; filename += time; filename += ".root";
         } else if  (strcmp(tmp.Data(), "dt") == 0) {
            filename += "_"; filename += date;
            filename += "_"; filename += time;
            filename += ".root";
         } else {
            filename = name;
         }//if
         file = new TFile(filename, "CREATE", title);
      } else {
         cerr << "Error: File <" << name << "> does already exist" << endl;
         delete [] (char*)fname;
         return 0;
      }//if

      if (!file || file->IsZombie()) {
         fAbort = kTRUE;
      } else if (file->IsOpen()) {
         if (XManager::fgVerbose) {
            cout << "Creating new file <" << filename.Data() << ">..." << endl;
         }//if
         fAbort = kFALSE;
         delete [] (char*)fname;
         return file;
      }//if

      delete [] (char*)fname;
   }//if

   cerr << "Error: Could not create file <" << name << ">" << endl;
   SafeDelete(file);

   return 0;
}//NewFile

//______________________________________________________________________________
TFile *XManager::OpenFile(const char *name, Option_t *option, Bool_t &isOwner)
{
   // Open root file and return pointer to TFile
   // Note: If file with name is already open, isOwner is set to kFALSE.
   //       This is important when used in GUI where only one manager can
   //       be the owner of file.
   // Note: fOption="R" to allow updating root files from within R due to
   //       the strange R memory management
   if(kCS) cout << "------XManager::OpenFile------" << endl;

   isOwner = kFALSE;

   // convert option toupper
   TString opt = TString(option);
//   TString opt = Path2Name(option,"",":");
//   TString cmd = Path2Name(option,":","");
   opt.ToUpper();

// Abort if option is RECREATE file
   if (strcmp(opt.Data(), "RECREATE") == 0) {
      cerr << "Error: Trying to recreate existing file <" << name << ">" << endl;
      return 0;
   }//if

   TFile *file = 0;
   const char *fname;
   if ((fname = gSystem->ExpandPathName(name))) {
      if (strcmp(fOption.Data(), "R") == 0) {
//?         file = TFile::Open(name, opt.Data());
         file = TFile::Open(fname, opt.Data());
         isOwner = kTRUE;
      } else {
         file = gROOT->GetFile(fname);

         // need to set correct option (why?)
         if (file) {
            if (strcmp(file->GetOption(), "UPDATE") == 0) {
               file->ReOpen(file->GetOption());
//should this be done?? seems to have no effect??
               isOwner = this->GetFileOwner();
            } else {
               file->ReOpen(opt.Data());
            }//if
         }//if
//wrong?      file->SetOption(opt.Data());
//wrong?      file->SetWritable(kTRUE);

         if (!file) { //do not use "else"
            file = TFile::Open(name, opt.Data());
            isOwner = kTRUE;
         }//if
      }//if

      delete [] (char*)fname;
   }//if

   if (!file || file->IsZombie()) {
      fAbort = kTRUE;
   } else if (file->IsOpen()) {
      if (isOwner) {  //if file is opened for the first time
         if (XManager::fgVerbose) {
            cout << "Opening file <" << name << "> in <" << option << "> mode..."
                 << endl;
         }//if
      }//if
      return file;
   }//if

   cerr << "Error: Could not open file <" << name << ">" << endl;
   SafeDelete(file);
   fAbort = kTRUE;
   return 0;
}//OpenFile

//______________________________________________________________________________
void XManager::CloseFile(TFile *file, Option_t *option)
{
   // Close root file
   if(kCS) cout << "------XManager::CloseFile------" << endl;

   if (file) {
      file->Close(option);
   }//if
}//CloseFile

//______________________________________________________________________________
void XManager::DeleteFile(TFile *file)
{
   // Delete root file
   if(kCS) cout << "------XManager::DeleteFile------" << endl;

//?? Delete root file from hard disk?
   SafeDelete(file);
}//DeleteFile

//______________________________________________________________________________
void XManager::DeleteDirectory(const char *name, const char *cycle)
{
   // Delete directory name in root file
   if(kCS) cout << "------XManager::DeleteDirectory------" << endl;

   if (!fFile) return;

   TString dir = TString(name);
   if (strcmp(cycle, "") != 0) dir = dir + ";" + TString(cycle);

   fFile->cd();
   fFile->Delete(dir);
}//DeleteDirectory

//______________________________________________________________________________
XFolder *XManager::NewContent(const char *name, const char *title, const char *type)
{
   // Create new file content
   if(kCS) cout << "------XManager::NewContent------" << endl;

   XFolder *content = new XFolder(name, title, type);
   // set owner to allow deletion of treesets in content
   content->SetOwner(kTRUE);
   //note: cannot add XFolder! gROOT->GetListOfBrowsables()->Add(content, name);

   return content;
}//NewContent

