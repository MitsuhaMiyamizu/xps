// File created: 08/05/2002                          last modified: 12/30/2010
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

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <cstdlib>
#include <new>  //needed for new (nothrow)

#include "TF1.h"
#include "TGraph.h"
#include "TGraphSmooth.h"

//for Benchmark test only:
#include <TBenchmark.h>

#include "XPSSelector.h"

//debug: print function names
const Bool_t  kCS  = 0; 
const Bool_t  kCSa = 0; //debug: print function names in loops

ClassImp(XSelector);
ClassImp(XRankSelector);
ClassImp(XProbeSelector);
ClassImp(XUnitSelector);
ClassImp(XUserSelector);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSelector                                                            //
//                                                                      //
// Class for selection of units used for normalization (default)        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSelector::XSelector()
          :XAlgorithm()
{
   // Default Selector constructor
   if(kCS) cout << "---XSelector::XSelector(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSelector::XSelector(const char *name, const char *type)
          :XAlgorithm(name, type)
{
   // Normal Selector constructor
   if(kCS) cout << "---XSelector::XSelector------" << endl;

}//Constructor

//______________________________________________________________________________
XSelector::XSelector(const XSelector &selector) 
          :XAlgorithm(selector)
{
   // Selector copy constructor
   if(kCS) cout << "---XSelector::XSelector(copy)------" << endl;

   fOption = selector.fOption;
}//CopyConstructor

//______________________________________________________________________________
XSelector& XSelector::operator=(const XSelector& rhs)
{
   // Selector assignment operator.
   if(kCS) cout << "---XSelector::operator=------" << endl;

   if (this != &rhs) {
      XAlgorithm::operator=(rhs);

      fOption = rhs.fOption;
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XSelector::~XSelector()
{
   // Selector destructor
   if(kCS) cout << "---XSelector::~XSelector------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XSelector::Calculate(Int_t n, Double_t * /*x*/, Double_t * /*y*/, Int_t *msk)
{
   // Leave array msk unchanged (none) or set mask to one (all)
   if(kCS) cout << "------XSelector::Calculate------" << endl;

   if (strcmp(fOption.Data(), "none") == 0) {
      // do not change msk
   } else if (strcmp(fOption.Data(), "all") == 0) {
      // set all entries of mask array to 1
      for (Int_t i=0; i<n; i++) msk[i] = 1;
   } else {
      cerr << "Error: Default selector does not have option <" << fOption.Data()
           << ">! Aborting execution." << endl;
      return errAbort;
   }//if

   return errNoErr;
}//Select

//______________________________________________________________________________
void XSelector::SetOptions(Option_t *opt)
{
   // Set selector option
   if(kCS) cout << "------XSelector::SetOptions------" << endl;

   fOption = opt; fOption.ToLower();
}//SetOptions

//______________________________________________________________________________
Int_t *XSelector::SetMask(Int_t /*n*/, Int_t *arr)
{
   // Return array arr
   if(kCSa) cout << "------XSelector::SetMask------" << endl;

   return arr;
}//SetMask

//______________________________________________________________________________
Int_t XSelector::SetFlag(Int_t n, Int_t *arr)
{
   // Return one if sum of flags for selected entry[i] is at least numflags
   if(kCSa) cout << "------XSelector::SetFlag------" << endl;

   Int_t numflags = (Int_t)fPars[3];
   numflags = (numflags > n) ? n : numflags;

   Int_t flags = 0;
   for (Int_t i=0; i<n; i++) {
      flags += arr[i];
//??      flags += (arr[i] >= 0) ? arr[i] : 0;  //NA = -1
   }//for_i

   return (flags >= numflags) ? 1 : 0;
}//SetFlag


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XRankSelector                                                        //
//                                                                      //
// Class for selection of units based on ranking algorithm              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XRankSelector::XRankSelector()
              :XSelector()
{
   // Default RankSelector constructor
   if(kCS) cout << "---XRankSelector::XRankSelector(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XRankSelector::XRankSelector(const char *name, const char *type)
              :XSelector(name, type)
{
   // Normal RankSelector constructor
   if(kCS) cout << "---XRankSelector::XRankSelector------" << endl;

}//Constructor

//______________________________________________________________________________
XRankSelector::~XRankSelector()
{
   // RankSelector destructor
   if(kCS) cout << "---XRankSelector::~XRankSelector------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XRankSelector::Calculate(Int_t n, Double_t *x, Double_t *y, Int_t *msk)
{
   // Select units used for normalization:
   // Fill msk with mask to select non-varying genes
   if(kCS) cout << "------XRankSelector::Calculate------" << endl;

   Int_t err    = errNoErr;
   Int_t length = 0;        
   Int_t niter  = 0;        

// Check percentage parameter
   if ((fPars[2] > 1.0) || (fPars[2] < 0.0)) {
      cout << "Warning: Percentage parameter to select <" << fPars[2]
           << "> must be within [0,1)! Setting to default value 0.25." << endl;
      fPars[2] = 0.25;
   } else if (fPars[2] == 1.0) {
      // set all entries of mask array to 1
      for (Int_t i=0; i<n; i++) msk[i] = 1;

      return errNoErr;
   }//if

   Double_t max      = 0;
   Double_t cutoff   = fPars[0];
   Double_t minunits = n*fPars[2];

// Initialize local arrays
   Int_t    *idxX  = 0;
   Int_t    *idxY  = 0;
   Int_t    *rankX = 0;
   Int_t    *rankY = 0;
   Double_t *diff  = 0;
   if (!(idxX  = new (nothrow) Int_t[n]))    {err = errInitMemory; goto cleanup;} 
   if (!(idxY  = new (nothrow) Int_t[n]))    {err = errInitMemory; goto cleanup;} 
   if (!(rankX = new (nothrow) Int_t[n]))    {err = errInitMemory; goto cleanup;} 
   if (!(rankY = new (nothrow) Int_t[n]))    {err = errInitMemory; goto cleanup;} 
   if (!(diff  = new (nothrow) Double_t[n])) {err = errInitMemory; goto cleanup;} 

// Sort and rank arrays x and y 
   TGraphSmooth::Rank(n, x, idxX, rankX, kFALSE);
   TGraphSmooth::Rank(n, y, idxY, rankY, kFALSE);

// Calculate absolute difference of ranks
   for (Int_t i=0; i<n; i++) {
      diff[i] = TMath::Abs(rankX[i] - rankY[i]);
      max = (max > diff[i]) ? max : diff[i];
   }//for_i

// Set cutoff = 1 if variance vector is null (max=0)
   if (max == 0) {
      cout << "Warning: Variance vector is null! Setting cutoff = 1." << endl;
      cutoff = 1;
   }//if

// Calculate cutoff value if initialized to 0
   if (cutoff == 0) {
      cutoff = this->Cutoff(n, diff, max, kFALSE);
      if (cutoff < 0) {err = (Int_t)cutoff; goto cleanup;}  //use cutoff as err
   }//if

// Increase cutoff until number of selected units > minunits
   do {
      if (XManager::fgVerbose) {
         cout << "Cutoff value is <" << cutoff << ">" << endl;
      }//if
      // calculate mask to select non-varying genes
      length = 0;        
      for (Int_t i=0; i<n; i++) {
         msk[i] = (diff[i] < cutoff) ? 1 : 0;
         length += msk[i];
      }//for_i
      cutoff += 0.2 * cutoff;  //increase in percent
      if (XManager::fgVerbose) {
         cout << "Number of selected genes is <" << length << ">";
         if (niter > 0) cout << " after <" << niter << "> iterations.";
         cout << endl;
      }//if

      // need to prevent infinite loop!
//?? as parameter??    if (niter++ ==  fPars[4]) {
      if (niter++ == 10) {
         cout << "Warning: Could not reach minimum #units after <"
              << niter << "> iterations!" << endl;
         break;
      }//if
   } while (length < minunits);

// Delete local variables
cleanup:
   delete [] diff;
   delete [] rankY;
   delete [] rankX;
   delete [] idxY;
   delete [] idxX;

   return errNoErr;
}//Calculate

//______________________________________________________________________________
Double_t XRankSelector::Cutoff(Int_t n, const Double_t *arr, Double_t max, 
                        Bool_t showGraph)
{
   // Calculate cutoff value
   if(kCS) cout << "------XRankSelector::Cutoff------" << endl;

   Int_t err = errNoErr;
   Int_t len = 0;
/*
// Special case: variance vector is null!
   if (max == 0) {
      cerr << "Error: Variance vector is null! Aborting normalization." << endl;
      cout << "       Please set <cutoff> value manually." << endl;
      return errAbort;
   }//if
*/

   Double_t d = 0;
   Double_t trim, minX, maxX;
   Double_t par0, par1;
   TF1     *f1 = 0;
   TGraph  *gr = 0;

// Initialize local arrays
   Int_t    *index = 0;
   Double_t *dX    = 0;
   Double_t *dY    = 0;
   if (!(index = new (nothrow) Int_t[n]))    {err = errInitMemory; goto cleanup;} 
   if (!(dX    = new (nothrow) Double_t[n])) {err = errInitMemory; goto cleanup;}
   if (!(dY    = new (nothrow) Double_t[n])) {err = errInitMemory; goto cleanup;} 

// Sort difference array arr
   TMath::Sort(n, arr, index, kFALSE);
   for (Int_t i=0; i<n; i++) {
      d = arr[index[i]];
      if (d > 0) {
         dX[len] = len;
         dY[len] = TMath::Log10(d);
         len++;
      }//if
   }//for_i
   len--;

// Special case: len = 1  //is this allowed????
   if (len == 1) {
      cout << "Warning: Array <difference of ranks> has only ONE entry!" << endl;
      par0 = dY[0];
      goto cleanup;
   }//if

// Set function for fitting
   trim = fPars[1];
   minX = 0;
   maxX = len;  //=1?
   if (trim > 0 && trim < 0.45) {
      minX = len * trim;
      maxX = len * (1 - trim); 
   }//if
   f1 = new TF1("f1", "pol1", minX, maxX);

// Fit graph with function
   gr = new TGraph(len, dX, dY);
   gr->Fit("f1", "RQ"); //"R"-use function range
   par0 = f1->GetParameter(0);
   par1 = f1->GetParameter(1);

// Delete locally created variables
cleanup:
//   delete gr;
//   delete f1;
   SafeDelete(gr);
   SafeDelete(f1);
   if (dY)    delete [] dY;
   if (dX)    delete [] dX;
   if (index) delete [] index;

// Return cutoff value
//   return TMath::Power(10, par0);
   return ((err == errNoErr) ? TMath::Power(10, par0) : (Double_t)err);
}//Cutoff


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProbeSelector                                                       //
//                                                                      //
// Class for selection of probes                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProbeSelector::XProbeSelector()
               :XSelector()
{
   // Default ProbeSelector constructor
   if(kCS) cout << "---XProbeSelector::XProbeSelector(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XProbeSelector::XProbeSelector(const char *name, const char *type)
               :XSelector(name, type)
{
   // Normal ProbeSelector constructor
   if(kCS) cout << "---XProbeSelector::XProbeSelector------" << endl;

}//Constructor

//______________________________________________________________________________
XProbeSelector::XProbeSelector(const XProbeSelector &selector) 
               :XSelector(selector)
{
   // ProbeSelector copy constructor
   if(kCS) cout << "---XProbeSelector::XProbeSelector(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XProbeSelector& XProbeSelector::operator=(const XProbeSelector& rhs)
{
   // ProbeSelector assignment operator.
   if(kCS) cout << "---XProbeSelector::operator=------" << endl;

   if (this != &rhs) {
      XSelector::operator=(rhs);
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XProbeSelector::~XProbeSelector()
{
   // ProbeSelector destructor
   if(kCS) cout << "---XProbeSelector::~XProbeSelector------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XProbeSelector::Calculate(Int_t n, Double_t * /*x*/, Double_t * /*y*/, Int_t *msk)
{
   // Fill msk with mask to select probes
   if(kCS) cout << "------XProbeSelector::Calculate------" << endl;

   if (strcmp(fOption.Data(), "none") == 0) {
      // do not change msk
   } else if (strcmp(fOption.Data(), "all") == 0) {
      // set all entries of mask array to 1
      for (Int_t i=0; i<n; i++) msk[i] = 1;
   } else if (strcmp(fOption.Data(), "both") == 0) {
      // set msk to one for msk=1 or msk=0 (PM and MM probes for genes only)
      for (Int_t i=0; i<n; i++) msk[i] = (msk[i] == 1) ? 1 : ((msk[i] == 0) ? 1 : 0);
   } else if (strcmp(fOption.Data(), "pmonly") == 0) {
      // set msk to one for msk=1 (PM probes only)
      for (Int_t i=0; i<n; i++) msk[i] = (msk[i] == 1) ? 1 : 0;
   } else if (strcmp(fOption.Data(), "mmonly") == 0) {
      // set msk to one for msk=0 (MM probes only)
      for (Int_t i=0; i<n; i++) msk[i] = (msk[i] == 0) ? 1 : 0;
   } else if (strcmp(fOption.Data(), "genome") == 0) {
      msk = this->SetGenomeMask(n, msk);
      if (msk == 0) return errInitParameters;
   } else if (strcmp(fOption.Data(), "exon") == 0) {
      msk = this->SetExonMask(n, msk);
      if (msk == 0) return errInitParameters;
   } else {
      cerr << "Error: Probe selector does not have option <" << fOption.Data()
           << ">! Aborting execution." << endl;
      return errAbort;
   }//if

   return errNoErr;
}//Calculate

//______________________________________________________________________________
Int_t *XProbeSelector::SetGenomeMask(Int_t n, Int_t *arr)
{
   // Return modified array arr, depending on parameter settings
   if(kCS) cout << "------XProbeSelector::SetGenomeMask------" << endl;

// Test number of parameters (at least one)
   if (TestNumParameters(1) != errNoErr) return 0;

   Int_t typepm = (Int_t)fPars[0];
   Int_t typemm = 0;
   // check if parameter for typemm exists
   if (fNPar > 1) typemm = (Int_t)fPars[1];

   // convert negative typemm to int
   typemm = (typemm >= 0) ? typemm : (abs(typemm) << 15);

// Set bit mask to level
   XBitSet bitmsk;
   bitmsk.ResetBit(XBitSet::kBitMask);
   bitmsk.SetBit(typepm);

// Set mask=1 if typepm is combination of levels
   for (Int_t i=0; i<n; i++) {
      Int_t x = arr[i];

      if ((x == eMMAT || x == eMMST) && (bitmsk.TestBit(x) == kTRUE)) {
            arr[i] = 0;
      } else if (x > 0 && bitmsk.TestBit(x) == kTRUE) {
         arr[i] = 1;
      } else {
         arr[i] = ((fNPar > 1) && (x == typemm)) ? 0 : eINITMASK;
      }//if
   }//for_i

   if (XManager::fgVerbose) {
      cout << "      setting selector mask for typepm <" << typepm << ">" << endl;
   }//if

   return arr;
}//SetGenomeMask

//______________________________________________________________________________
Int_t *XProbeSelector::SetExonMask(Int_t n, Int_t *arr)
{
   // Return modified array arr, depending on parameter settings
   if(kCS) cout << "------XProbeSelector::SetExonMask------" << endl;

   return this->SetGenomeMask(n, arr);
}//SetExonMask


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnitSelector                                                        //
//                                                                      //
// Class for selection of units                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUnitSelector::XUnitSelector()
              :XProbeSelector()
{
   // Default UnitSelector constructor
   if(kCS) cout << "---XUnitSelector::XUnitSelector(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XUnitSelector::XUnitSelector(const char *name, const char *type)
              :XProbeSelector(name, type)
{
   // Normal UnitSelector constructor
   if(kCS) cout << "---XUnitSelector::XUnitSelector------" << endl;

}//Constructor

//______________________________________________________________________________
XUnitSelector::XUnitSelector(const XUnitSelector &selector) 
//old              :XSelector(selector)
              :XProbeSelector(selector)
{
   // UnitSelector copy constructor
   if(kCS) cout << "---XUnitSelector::XUnitSelector(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XUnitSelector& XUnitSelector::operator=(const XUnitSelector& rhs)
{
   // UnitSelector assignment operator.
   if(kCS) cout << "---XUnitSelector::operator=------" << endl;

   if (this != &rhs) {
      XSelector::operator=(rhs);
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XUnitSelector::~XUnitSelector()
{
   // UnitSelector destructor
   if(kCS) cout << "---XUnitSelector::~XUnitSelector------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XUnitSelector::Calculate(Int_t n, Int_t * /*x*/, Int_t *msk)
{
   // Fill msk with mask to select units
   if(kCS) cout << "------XUnitSelector::Calculate------" << endl;

   if (strcmp(fOption.Data(), "none") == 0) {
      // do not change msk
   } else if (strcmp(fOption.Data(), "all") == 0) {
      // set all entries of mask array to 1
      for (Int_t i=0; i<n; i++) msk[i] = 1;
   } else if (strcmp(fOption.Data(), "gene") == 0) {
      // do not change msk
   } else if (strcmp(fOption.Data(), "genome") == 0) {
      msk = SetGenomeMask(n, msk);
      if (msk == 0) return errInitParameters;
   } else if (strcmp(fOption.Data(), "exon") == 0) {
      msk = SetExonMask(n, msk);
      if (msk == 0) return errInitParameters;
   } else {
      cerr << "Error: Unit selector does not have option <" << fOption.Data()
           << ">! Aborting execution." << endl;
      return errAbort;
   }//if

   return errNoErr;
}//Calculate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUserSelector                                                        //
//                                                                      //
// Class for selection of units by importing a mask file                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUserSelector::XUserSelector()
              :XSelector()
{
   // Default UserSelector constructor
   if(kCS) cout << "---XUserSelector::XUserSelector(default)------" << endl;

   fInput = "";
}//Constructor

//______________________________________________________________________________
XUserSelector::XUserSelector(const char *name, const char *type)
              :XSelector(name, type)
{
   // Normal UserSelector constructor
   if(kCS) cout << "---XUserSelector::XUserSelector------" << endl;

   fInput = "";
}//Constructor

//______________________________________________________________________________
XUserSelector::~XUserSelector()
{
   // UserSelector destructor
   if(kCS) cout << "---XUserSelector::~XUserSelector------" << endl;

}//Destructor

//______________________________________________________________________________
void XUserSelector::SetOptions(Option_t *opt)
{
   // Fill msk with user-defined mask to select non-varying genes
   if(kCS) cout << "------XUserSelector::SetOptions------" << endl;

//	TString optcpy = opt;
//   char *options = (char*)optcpy.Data();

   char *options = new char[strlen(opt) + 1];
   char *doption = options;
   options = strcpy(options, opt);

   if (NumSeparators(opt, ":") == 0) {
      fInput  = strtok(options, ":");
      fOption = ""; 
   } else {
      fInput  = strtok(options, ":");
      fOption = strtok(NULL, ":");
   }//if

   delete [] doption;
}//SetOptions

//______________________________________________________________________________
Int_t XUserSelector::Calculate(Int_t n, Double_t * /*x*/, Double_t * /*y*/,
      Int_t *msk)
{
   // Fill msk with user-defined mask to select non-varying genes
   if(kCS) cout << "------XUserSelector::Calculate------" << endl;

   if (strcmp(fInput.Data(), "all") == 0) {
      // set all entries of mask array to 1
      for (Int_t i=0; i<n; i++) {
         msk[i] = 1;
      }//for_i

      return errNoErr;
   } else {
      // initialize mask array
      for (Int_t i=0; i<n; i++) {
         msk[i] = 0;
      }//for_i

      // import mask from fInput
      return this->Import(fInput.Data(), n, msk);
   }//if

   return errAbort;
}//Calculate

//______________________________________________________________________________
Int_t XUserSelector::Import(const char *infile, Int_t n, Int_t *msk, char delim)
{
   // Import user-defined mask from infile
   if(kCS) cout << "------XUserSelector::Import------" << endl;

// Open infile containing mask
   ifstream input(infile, ios::in);
   if (!input) {
      cerr << "Error: File <" << infile << "> does not exist." << endl;
      return errReadingInput;
   }//if

// Check for header line UNIT_ID
//ccc char memory problem??
   char  nextline[kBufSize];
   input.getline(nextline, kBufSize, delim);
   if (strncmp("UNIT_ID", nextline, 7) != 0) {
      cerr << "Error: First column of file <" << infile << "> must be UNIT_ID."
           << endl;
      return errHeaderLine;
   }//if

// Import mask file
   Int_t unitID, flag;
   Int_t numunits = 0;
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.fail() || (numunits > n)) break;
      sscanf(nextline,"%i %i", &unitID, &flag);
//      sscanf(nextline,"%i%*c%i", &unitID, &flag);  //to allow ,; as separator?
      if (unitID < n) {
         msk[unitID] = flag;
      } else {
         cerr << "Error: UNIT_ID <" << unitID << "> is larger than <." << n
              << ">." << endl;
         input.close();
         return errReadingInput;
      }//if
      numunits++;
   }//while
   input.close();

// Check percentage parameter
   if ((fPars[0] >= 1.0) || (fPars[0] < 0.0)) {
      cout << "Warning: Percentage parameter to select <" << fPars[0]
           << "> must be witin [0,1)! Setting to default value 0.25." << endl;
      fPars[0] = 0.25;
   }//if

// Check if number of selected units > minunits
   Int_t minunits = (Int_t)(n*fPars[0]);
   if (numunits < minunits) {
      cerr << "Error: Number of imported units <" << numunits << "> " 
           << "is less than minimum requested units <" << minunits << ">."
           << endl;
      return errGeneral;
   }//if

   return errNoErr;
}//Import

