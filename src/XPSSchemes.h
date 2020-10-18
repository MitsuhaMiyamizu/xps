// File created: 05/18/2002                          last modified: 10/17/2010
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

#ifndef __XPSSchemes__
#define __XPSSchemes__

#include "XPSBase.h"

// Error messages for XPSSchemes
enum EErrorSchemes {
   errExtension  = -101,
   errCDFVersion = -102,
   errAnnVersion = -103,
};

// Exon level of annotation
enum ELevel {
//set to 2^n to allow combinations (exception no level)
   eMETACORE     = 8192,
   eMETAEXTENDED = 4096,
   eMETAFULL     = 2048,
   ePARACORE     = 1024,
   ePARAEXTENDED = 512,
   ePARAFULL     = 256,
   eAMBIGUOUS    = 128,
   eFREE         = 64,

   eNOLEVEL      = -64,
};

// Probe type
enum EProbeType {
   ePMAT = 16, ePMST = 32,
   eMMAT =  4, eMMST =  8,
};

// Exon probeset type: fUnitType for ExonArray:
// - eMAIN is subdivided into subtypes ELevel 
// - eMAINEXON is added to ELevel for "exon:main"
enum EProbesetType {
   eMAINEXON    =  16384, //only "exon:main"

   eMAIN        =  2,
   eCONTROLAFFX =  1,
   eCONTROLCHIP =  0,

   eGENOMIC     = -1,  //32768:  bitshift << 15
   eANTIGENOMIC = -2,  //65536
   eINTRON      = -4,  //131072
   eEXON        = -8,  //262144 //only "exon", not "exon:main"
   eUNMAPPED    = -16, //524288
//   reserverd    = -32, //1048576
   eUNKNOWNTYPE = -64, //2097152

   eINITMASK    = -16384, //init mask != 0 since mask will be set to 1 or 0
};

// Exon controlchip type
enum EControlchipType {
//   eGENOMICST  = -1,  eANTIGENOMICST = -2,
//   eINTRONST   = -4,  eEXONST        = -8,
//   eUNMAPPEDST = -16,
   eUNKNOWN   = -64,
   eJUMBOAT   = -100, eJUMBOST   = -101,
   eTHERMOAT  = -102, eTHERMOST  = -103,
   eTRIGRIDAT = -104, eTRIGRIDST = -105,
   eGENERICAT = -106, eGENERICST = -107,
   eBLANK     = -108,
   eCTRLPMST  = -109,
};

// Exon bounded info
enum EBounded {
   eBND_T_NOBND_T = 3, //bounded==TRUE && NoBoundedEvidence==TRUE
   eBND_F_NOBND_T = 2, //bounded==FALSE && NoBoundedEvidence==TRUE
   eBND_T_NOBND_F = 1, //bounded==TRUE && NoBoundedEvidence==FALSE
   eBND_F_NOBND_F = 0, //bounded==FALSE && NoBoundedEvidence==FALSE
};

//length of oligoprobe, i.e. #nucleotides
const Int_t kProbeLength = 25;

// tree scheme extensions and types
extern const char *kExtenScheme[];
extern const char *kTypeScheme[];
extern const char *kExtenAnnot[];
extern const char *kTypeAnnot[];


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSchemeManager                                                       //
//                                                                      //
// Class for managing import of microarray scheme definitions           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSchemeManager: public XManager {

   private:

   protected:
      virtual XTreeSet *NewTreeSet(const char *settype); 

   public:
      XSchemeManager();
      XSchemeManager(const char *name, const char *title = "",
                     Int_t verbose = kTRUE);
      virtual ~XSchemeManager();

      using XManager::ExportTrees;
      virtual Int_t ExportTrees(const char *exten, const char *varlist, 
                       ofstream &output, const char *sep);

      Int_t NewScheme(const char *name, const char *infile,
               const char *type = "scheme",
               const char *sep = "", char delim = '\n', Int_t split = 99);
      Int_t UpdateScheme(const char *name, const char *infile,
               const char *type = "scheme",
               const char *sep = "", char delim = '\n', Int_t split = 99);
      void  DeleteScheme(const char *name);

      Int_t NewAnnotation(const char *name, const char *infile,
               const char *type = "transcript",
               const char *sep = "", char delim = '\n', Int_t split = 99);
      Int_t UpdateAnnotation(const char *name, const char *infile,
               const char *type = "transcript",
               const char *sep = "", char delim = '\n', Int_t split = 99);

      virtual Int_t New(const char *name, const char *dir, const char *type,
                       const char *data = "Schemes", Option_t *option = "");
      virtual Int_t HandleError(Int_t err, const char *name1 = "",
                       const char *name2 = "");

      ClassDef(XSchemeManager,1) //SchemeManager
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLayoutTreeInfo                                                      //
//                                                                      //
// Class containing info about layout tree stored in fUserInfo of tree  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XLayoutTreeInfo: public XTreeInfo {

   protected:
      Int_t    fNRows;      //number of rows of chip-matrix (=TotalX)
      Int_t    fNCols;      //number of columns of chip-matrix (=TotalY)
      Short_t  fSequential; //indicates that layout data are sequential

   public :
      XLayoutTreeInfo();
      XLayoutTreeInfo(const char *name, const char *title);
      virtual ~XLayoutTreeInfo();

      virtual void     AddUserInfo(XTreeSet *set);
      virtual void     AddUserInfo(Int_t nrows, Int_t ncols, Short_t seq);
      virtual Double_t GetValue(const char *name);

      ClassDef(XLayoutTreeInfo,1) //LayoutTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProbeTreeInfo                                                       //
//                                                                      //
// Class containing info about probe tree stored in fUserInfo of tree   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProbeTreeInfo: public XTreeInfo {

   protected:
      Int_t      fMinGC;      //minimum GC content of oligos on array
      Int_t      fMaxGC;      //maximum GC content of oligos on array
      Double_t   fMinTm;      //minimum melting temperature of oligos on array
      Double_t   fMaxTm;      //maximum melting temperature of oligos on array

   public :
      XProbeTreeInfo();
      XProbeTreeInfo(const char *name, const char *title);
      virtual ~XProbeTreeInfo();

      virtual void     AddUserInfo(XTreeSet *set);
      virtual void     AddUserInfo(Int_t minGC, Int_t maxGC, Double_t minTm,
                          Double_t maxTm);
      virtual Double_t GetValue(const char *name);

      ClassDef(XProbeTreeInfo,1) //ProbeTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSchemeTreeInfo                                                      //
//                                                                      //
// Class containing info about scheme tree stored in fUserInfo of tree  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSchemeTreeInfo: public XTreeInfo {

   protected:
      Int_t    fNRows;      //number of rows of chip-matrix (=TotalX)
      Int_t    fNCols;      //number of columns of chip-matrix (=TotalY)
      Int_t    fNProbes;    //number of annotated probes on array
      Int_t    fNControls;  //number of controls on chip
      Int_t    fNGenes;     //number of genes on chip
      Int_t    fNUnits;     //number of units on chip

   public :
      XSchemeTreeInfo();
      XSchemeTreeInfo(const char *name, const char *title);
      virtual ~XSchemeTreeInfo();

      virtual void     AddUserInfo(XTreeSet *set);
      virtual void     AddUserInfo(Int_t nrows, Int_t ncols, Int_t nprobes,
                          Int_t nctrls, Int_t ngenes, Int_t nunits);
      virtual Double_t GetValue(const char *name);

      ClassDef(XSchemeTreeInfo,1) //SchemeTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnitTreeInfo                                                        //
//                                                                      //
// Class containing info about unit tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUnitTreeInfo: public XTreeInfo {

   protected:
      Int_t    fNControls;    //number of controls on chip
      Int_t    fNAffx;        //number of AFFX controls on array
      Int_t    fNGenes;       //number of genes on chip
      Int_t    fMinNPairs;    //minimum number of PM/MM pairs
      Int_t    fMaxNPairs;    //maximum number of PM/MM pairs
      Int_t    fMinNCells;    //number of cells with minimum number of PM/MM pairs
      Int_t    fMaxNCells;    //number of cells with maximum number of PM/MM pairs

   public :
      XUnitTreeInfo();
      XUnitTreeInfo(const char *name, const char *title);
      virtual ~XUnitTreeInfo();

      virtual void     AddUserInfo(XTreeSet *set);
      virtual void     AddUserInfo(Int_t nctrls, Int_t naffx, Int_t ngenes,
                          Int_t nmin, Int_t min, Int_t nmax, Int_t max);
      virtual Double_t GetValue(const char *name);

      ClassDef(XUnitTreeInfo,1) //UnitTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeTreeInfo                                                      //
//                                                                      //
// Class containing info about exon tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeTreeInfo: public XTreeInfo {

   protected:
      Int_t    fNControls;  //number of controls on chip
      Int_t    fNAffx;      //number of AFFX controls on array
      Int_t    fNSubunits;  //number of subunits (exons, probesets) on chip
      Int_t    fMinNCells;  //minimum number of subunit cells
      Int_t    fMaxNCells;  //maximum number of subunit cells

   public :
      XGenomeTreeInfo();
      XGenomeTreeInfo(const char *name, const char *title);
      virtual ~XGenomeTreeInfo();

      virtual void     AddUserInfo(XTreeSet *set);
      virtual void     AddUserInfo(Int_t nctrls, Int_t naffx, Int_t nsubs,
                          Int_t mincells, Int_t maxcells);
      virtual Double_t GetValue(const char *name);

      ClassDef(XGenomeTreeInfo,1) //GenomeTreeInfo
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDNAChip                                                             //
//                                                                      //
// Base class for DNA-chips                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XDNAChip: public XTreeSet {

   protected:
      TString   fCxyTreeName;  //name of layout tree
      TString   fScmTreeName;  //name of scheme tree
      TString   fIdxTreeName;  //name of unit tree
      TString   fPrbTreeName;  //name of probe info tree
      TString   fAnnTreeName;  //name of annotation tree
      TString   fVersionAnnot; //annotation version or build
      Int_t     fNRows;        //number of rows of chip-matrix (=TotalX)
      Int_t     fNCols;        //number of columns of chip-matrix (=TotalY)
      Int_t     fNProbes;      //number of annotated probes on array
      Int_t     fNControls;    //number of controls on chip
      Int_t     fNGenes;       //number of genes on chip
      Int_t     fNUnits;       //number of units on chip
      Short_t   fSequential;   //indicates that layout data are sequential

   protected:
      virtual Int_t ExportTrees(const char *exten, const char *varlist, 
                       ofstream &output, const char *sep);
      virtual Int_t ExportLayoutTree(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportProbeTree(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitTree(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportTransAnnotTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportExonAnnotTree(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportProbesetAnnotTree(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportControlAnnotTree(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportLayoutXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportProbeXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportSchemeXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportUnitXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportTransAnnotXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportExonAnnotXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportProbesetAnnotXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ExportControlAnnotXML(Int_t /*n*/, TString * /*names*/, const char * /*varlist*/,
                       ofstream &/*output*/, const char * /*sep*/) {return errNoErr;}
      virtual Int_t ImportLayout(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split)        = 0;
      virtual Int_t ImportLayoutXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split)        = 0;
      virtual Int_t ImportProbeInfo(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split)        = 0;
      virtual Int_t ImportProbeInfoXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split)        = 0;
      virtual Int_t ImportScheme(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split)        = 0;
      virtual Int_t ImportSchemeXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split)        = 0;
      virtual Int_t ImportTransAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportTransAnnotationXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportExonAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportExonAnnotationXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbesetAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbesetAnnotationXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportControlAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportControlAnnotationXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);

      Int_t XY2Index(Int_t x, Int_t y) {return (x + y*fNCols);  }
      Int_t Index2X(Int_t index)       {return (index % fNCols);}
      Int_t Index2Y(Int_t index)       {return (index / fNCols);}

   public :
      XDNAChip();
      XDNAChip(const char *name, const char *title);
      XDNAChip(const XDNAChip &chip);
      XDNAChip& operator=(const XDNAChip& rhs);
      virtual ~XDNAChip();
		
      virtual Int_t ExportTreeType(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t ExportTreeXML(const char *exten, Int_t n, TString *names,  
                       const char *varlist, ofstream &output, const char *sep);
      virtual Int_t Import(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
		
      TString GetLayoutTree()  const {return fCxyTreeName;}
      TString GetSchemeTree()  const {return fScmTreeName;}
      TString GetUnitTree()    const {return fIdxTreeName;}
      TString GetProbeTree()   const {return fPrbTreeName;}
      TString GetAnnotTree()   const {return fAnnTreeName;}
      Int_t   GetNumRows()     const {return fNRows;}
      Int_t   GetNumColumns()  const {return fNCols;}
      Int_t   GetNumProbes()   const {return fNProbes;}
      Int_t   GetNumControls() const {return fNControls;}
      Int_t   GetNumGenes()    const {return fNGenes;}
      Int_t   GetNumUnits()    const {return fNUnits;}
      Short_t GetSequential()  const {return fSequential;}

      ClassDef(XDNAChip,1) //DNAChip
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMicroArray                                                          //
//                                                                      //
// Class for cDNA microarrays                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XMicroArray: public XDNAChip {

   private:

   protected:
      virtual Int_t ExportProbeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportProbeXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ImportLayout(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportLayoutXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfo(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfoXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportScheme(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportSchemeXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      
   public :
      XMicroArray();
      XMicroArray(const char *name, const char *title);
      XMicroArray(const XMicroArray &chip);
      XMicroArray& operator=(const XMicroArray& rhs);
      virtual ~XMicroArray();
		
      virtual void  PrintInfo();

      ClassDef(XMicroArray,1) //Microarray
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XOligoArray                                                          //
//                                                                      //
// Class for oligonucleotide arrays                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XOligoArray: public XDNAChip {

   private:

   protected:
      virtual Int_t ExportLayoutTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportProbeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportLayoutXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportProbeXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ImportLayout(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportLayoutXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfo(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfoXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportScheme(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportSchemeXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      
   public :
      XOligoArray();
      XOligoArray(const char *name, const char *title);
      XOligoArray(const XOligoArray &chip);
      XOligoArray& operator=(const XOligoArray& rhs);
      virtual ~XOligoArray();
		
      ClassDef(XOligoArray,1) //OligoArray
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChip                                                            //
//                                                                      //
// Class for Affymetrix GeneChip arrays                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGeneChip: public XOligoArray {

   protected:
      TString    fAnpTreeName; //name of probeset annotation tree
      Int_t      fNProbesets;  //number of probe_sets on chip
      Int_t      fNAffx;       //number of AFFX controls on array
      Int_t      fMinGC;       //minimum GC content of oligos on array
      Int_t      fMaxGC;       //maximum GC content of oligos on array
      Double_t   fMinTm;       //minimum melting temperature of oligos on array
      Double_t   fMaxTm;       //maximum melting temperature of oligos on array

   protected:
      virtual void  AddLayoutTreeInfo(TTree *tree, const char *name, Option_t *option = "");
      virtual void  AddProbeTreeInfo(TTree *tree, const char *name, Option_t *option = "");
      virtual void  AddSchemeTreeInfo(TTree *tree, const char *name, Option_t *option = "");
      virtual void  AddUnitTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nctrls, Int_t naffx, Int_t ngenes,
                       Int_t nmin, Int_t min, Int_t nmax, Int_t max);
      virtual Int_t ContentGC(const char *sequence, Double_t &Tm, const char *method);
      virtual Double_t MeltingTemperature(Int_t numGC, Int_t length, const char *method);
      virtual Int_t ExportProbeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportTransAnnotTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportProbeXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitXML(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t IsBinaryFile(ifstream &input);
      virtual Int_t ImportProbeInfo(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfoXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      virtual Int_t ReadBinaryHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadBinaryData(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportScheme(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportTransAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportSchemeXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);

      Short_t ProbeType(const char *type);
      
   public:
      XGeneChip();
      XGeneChip(const char *name, const char *title);
      XGeneChip(const XGeneChip &chip);
      XGeneChip& operator=(const XGeneChip& rhs);
      virtual ~XGeneChip();

      virtual void  PrintInfo();

      TString  GetProbesetAnnotTree() const {return fAnpTreeName;}
      Int_t    GetNumProbesets()      const {return fNProbesets;}
      Int_t    GetNumAffx()           const {return fNAffx;}
      Int_t    GetMinContentGC()      const {return fMinGC;}
      Int_t    GetMaxContentGC()      const {return fMaxGC;}
      Double_t GetMinMeltingT()       const {return fMinTm;}
      Double_t GetMaxMeltingT()       const {return fMaxTm;}
 
      ClassDef(XGeneChip,1) //GeneChip
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChip                                                             //
//                                                                      //
// Class for Affymetrix Mapping arrays                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XSNPChip: public XGeneChip {

   private:

   protected:
//??      virtual Int_t ExportProbeTree(Int_t n, TString *names, const char *varlist,
//??      virtual Int_t ImportAnnotation(ifstream &input, Option_t *option,
//??                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfo(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfoXML(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
//??      virtual Int_t ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      
   public:
      XSNPChip();
      XSNPChip(const char *name, const char *title);
      XSNPChip(const XSNPChip &chip);
      XSNPChip& operator=(const XSNPChip& rhs);
      virtual ~XSNPChip();
 
      ClassDef(XSNPChip,1) //SNPChip
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChip                                                          //
//                                                                      //
// Class for Affymetrix HuGene arrays                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeChip: public XGeneChip {

   protected:
      TList     *fAffxNames;   //! list of control->affx names

   protected:
      virtual void  AddGenomeTreeInfo(TTree *tree, const char *name, Option_t *option,
                       Int_t nctrls, Int_t naffx, Int_t nsubs, Int_t mincells, Int_t maxcells);
      virtual Int_t ExportLayoutTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportTransAnnotTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ImportLayout(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbeInfo(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ReadHeader(ifstream &input, const char *sep, char delim);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      virtual Int_t ImportTransAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);

      Short_t ProbesetType(const char *type);
      TString ProbesetTypeID2Name(Short_t id);
      Int_t   SchemeMask(Int_t type, Int_t subtype);
      Int_t   SchemeMask(Short_t xtype);
      
   public:
      XGenomeChip();
      XGenomeChip(const char *name, const char *title);
      XGenomeChip(const XGenomeChip &chip);
      XGenomeChip& operator=(const XGenomeChip& rhs);
      virtual ~XGenomeChip();
 
      ClassDef(XGenomeChip,1) //GenomeChip
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChip                                                            //
//                                                                      //
// Class for Affymetrix Exon arrays                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonChip: public XGenomeChip {

   private:
      TString    fExnTreeName; //name of exon unit tree
      TString    fPbsTreeName; //name of probeset unit tree
      TString    fAnxTreeName; //name of exon annotation tree
      Int_t      fNExonUnits;  //number of exon units on chip
      Int_t      fNExons;      //number of exons on chip

   protected:
//?      virtual Int_t ExportExonTree(Int_t n, TString *names, const char *varlist,
//?                       ofstream &output, const char *sep);
      virtual Int_t ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportUnitTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportExonAnnotTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ExportProbesetAnnotTree(Int_t n, TString *names, const char *varlist,
                       ofstream &output, const char *sep);
      virtual Int_t ReadData(ifstream &input, Option_t *option, const char *sep,
                       char delim, Int_t split);
      virtual Int_t ImportTransAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportExonAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportProbesetAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      virtual Int_t ImportControlAnnotation(ifstream &input, Option_t *option,
                       const char *sep, char delim, Int_t split);
      
      Int_t   LayoutToX(Int_t index);
      Int_t   LayoutToY(Int_t index);
      Int_t   ControlchipType(const char *type);
      TString LevelID2Level(Short_t id);
      Int_t   ProbesetLevel(const char *level, Short_t type);
      Int_t   SchemeMask(Short_t xtype, Short_t level, Short_t bound);

   public:
      XExonChip();
      XExonChip(const char *name, const char *title);
      XExonChip(const XExonChip &chip);
      XExonChip& operator=(const XExonChip& rhs);
      virtual ~XExonChip();

      TString GetExonUnitTree()     const {return fExnTreeName;}
      TString GetProbesetUnitTree() const {return fPbsTreeName;}
      TString GetExonAnnotTree()    const {return fAnxTreeName;}
      Int_t   GetNumExonUnits()     const {return fNExonUnits;}
      Int_t   GetNumExons()         const {return fNExons;}
 
      ClassDef(XExonChip,1) //ExonChip
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XPosition                                                            //
//                                                                      //
// Class for the (x,y)-position of a single microarray cell             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
//Note: tree wastes memory since TObject has fUniqueID, fBits!!!
//      => do not derive from TObject
class XPosition {

   protected:
      Int_t      fX;            //x-coordinate of feature cell
      Int_t      fY;            //y-coordinate of feature cell

   public:
      XPosition() {}
      virtual ~XPosition() {}

      void  SetX(Int_t x) {fX = x;}
      void  SetY(Int_t y) {fY = y;}
      Int_t GetX()  const {return fX;}
      Int_t GetY()  const {return fY;}

      ClassDef(XPosition,1) //Position
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLayout                                                              //
//                                                                      //
// Class for the layout of a single microarray cell                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XLayout: public XPosition {
   private:
      Int_t   fProbeID;      //probe_id

   public:
      XLayout() {}
      virtual ~XLayout() {}

      void  SetProbeID(Int_t id) {fProbeID = id;}
      Int_t GetProbeID()   const {return fProbeID;}

      ClassDef(XLayout,1) //Layout
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XScheme                                                              //
//                                                                      //
// General class describing the array layout                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XScheme: public XPosition {

   protected:
      Int_t       fUnitID;         //internal unit index
      Int_t       fProbeLen;       //length in bases of probe
      Int_t       fMask;           //e.g. for oligo: PM=1, MM=0, masked=-1
      
   public :
      XScheme()          {}
      virtual ~XScheme() {}

      void SetUnitID(Int_t unitID)   {fUnitID   = unitID;}
      void SetProbeLength(Int_t len) {fProbeLen = len;}
      void SetMask(Int_t mask)       {fMask     = mask;}

      Int_t GetUnitID()      const {return fUnitID;}
      Int_t GetProbeLength() const {return fProbeLen;}
      Int_t GetMask()        const {return fMask;}

      ClassDef(XScheme,1) //Scheme
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCScheme                                                            //
//                                                                      //
// Class describing the GeneChip array layout                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCScheme: public XScheme {

   protected:
      Int_t         fAtomNo;         //relative position of probe within unit
      char          fPBase;          //probe base at substitution position
      char          fTBase;          //target base at substitution position
      
   public :
      XGCScheme()          {}
      virtual ~XGCScheme() {}

      void SetAtomNumber(Int_t atom) {fAtomNo = atom;}
      void SetProbeBase(char pbase)  {fPBase = pbase;}
      void SetTargetBase(char tbase) {fTBase = tbase;}

      Int_t GetAtomNumber()    const {return fAtomNo;}
      char  GetProbeBase()     const {return fPBase;}
      char  GetTargetBase()    const {return fTBase;}

      ClassDef(XGCScheme,1) //GeneChipScheme
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonScheme                                                          //
//                                                                      //
// Class describing the Exon array layout                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonScheme: public XScheme {

   private:
      Int_t         fExonID;         //internal exon index
      Int_t         fProbesetID;     //internal probeset index
      
   public :
      XExonScheme()          {}
      virtual ~XExonScheme() {}

      void  SetExonID(Int_t id)     {fExonID     = id;}
      void  SetProbesetID(Int_t id) {fProbesetID = id;}
      Int_t GetExonID()       const {return fExonID;}
      Int_t GetProbesetID()   const {return fProbesetID;}

      ClassDef(XExonScheme,1) //ExonScheme
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProbe                                                               //
//                                                                      //
// Class describing properties of the probe                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProbe: public XPosition {

   protected:
      TString       fProbeName;   //name of probe, e.g. IMAGE clone name 
      TString       fSequence;    //sequence of probe 

   public :
      XProbe()          {}
      virtual ~XProbe() {}

      void SetProbeName(const char *name) {fProbeName = name;}
      void SetSequence(const char *seq)   {fSequence  = seq;}

      const char* GetProbeName() const {return fProbeName;}
      const char* GetSequence()  const {return fSequence;}

      ClassDef(XProbe,1) //Probe
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCProbe                                                             //
//                                                                      //
// Class describing properties of GeneChip probes                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCProbe: public XProbe {

   private:
      Int_t         fPosition;    //position of nucleotide in target sequence
      Int_t         fNumberGC;    //number of G/C in probe
      Double32_t    fTm;          //melting temperature of probe
      Short_t       fProbeType;   //probe type

   public :
      XGCProbe()          {}
      virtual ~XGCProbe() {}

      void SetPosition(Int_t pos)     {fPosition  = pos;}
      void SetNumberGC(Int_t gc)      {fNumberGC  = gc;}
      void SetTMelting(Double_t tm)   {fTm        = tm;}
      void SetProbeType(Short_t id)   {fProbeType = id;}

      Int_t       GetPosition()    const {return fPosition;}
      Int_t       GetNumberGC()    const {return fNumberGC;}
      Double_t    GetTMelting()    const {return fTm;}
      Short_t     GetProbeType()   const {return fProbeType;}

      ClassDef(XGCProbe,1) //GCProbe
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnitID                                                              //
//                                                                      //
// Class for the ID of a single microarray unit                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUnitID {

   private:
      Int_t    fUnitID;      //unit-ID

   public:
      XUnitID() {}
      virtual ~XUnitID() {}

      void  SetUnitID(Int_t unitID) {fUnitID = unitID;}
      Int_t GetUnitID()       const {return fUnitID;}

      ClassDef(XUnitID,1) //UnitID
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnit                                                                //
//                                                                      //
// Class describing the units on the array                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XUnit {

   protected:
      TString       fUnitName;    //name of unit: gene identifier or QC 
      Int_t         fUnitID;      //unit-ID
      Int_t         fUnitType;    //unit type: control<=0, expression=3
      Int_t         fNCells;      //number of cells for this unit

   public :
      XUnit()          {}
      virtual ~XUnit() {}

      void SetUnitName(const char *name) {fUnitName = name;}
      void SetUnitID(Int_t unitID)       {fUnitID   = unitID;}
      void SetUnitType(Int_t type)       {fUnitType = type;}
      void SetNumCells(Int_t numcells)   {fNCells   = numcells;}

      const char* GetUnitName() const {return fUnitName;}
      Int_t       GetUnitID()   const {return fUnitID;}
      Int_t       GetUnitType() const {return fUnitType;}
      Int_t       GetNumCells() const {return fNCells;}

      ClassDef(XUnit,1) //Unit
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGCUnit                                                              //
//                                                                      //
// Class describing the units on the GeneChip array                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGCUnit: public XUnit {

   protected:
      Int_t         fNAtoms;      //number of atoms per cell

   public :
      XGCUnit()          {}
      virtual ~XGCUnit() {}

      void  SetNumAtoms(Int_t numatoms) {fNAtoms = numatoms;}
      Int_t GetNumAtoms()         const {return fNAtoms;}

      ClassDef(XGCUnit,1) //GCUnit
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonUnit                                                            //
//                                                                      //
// Class describing the units on the GeneChip exon array                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonUnit: public XGCUnit {

   private:
      Int_t         fSubUnitID;    //probeset_id, exon_id, transcript_cluster_id
      Int_t         fNSubunits;    //number of subunits per cell

   public :
      XExonUnit()          {}
      virtual ~XExonUnit() {}

      void  SetSubUnitID(Int_t subid)    {fSubUnitID = subid;}
      void  SetNumSubunits(Int_t numsub) {fNSubunits = numsub;}
      Int_t GetSubUnitID()         const {return fSubUnitID;}
      Int_t GetNumSubunits()       const {return fNSubunits;}

      ClassDef(XExonUnit,1) //ExonUnit
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XAnnotation                                                          //
//                                                                      //
// Class describing the annotation of the units on the array            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XAnnotation {

   protected:
      TString       fName;        //name of unit: gene name
      TString       fSymbol;      //symbol of unit: gene symbol
      TString       fChromosome;  //chromosome: NA if not known or no gene 
      TString       fCytoBand;    //cytoband: NA if not known or no gene
      TString       fSequence;    //sequence: NA if not known or no gene

   public :
      XAnnotation()          {}
      virtual ~XAnnotation() {}

      void SetName(const char *name)         {fName       = name;}
      void SetSymbol(const char *symbol)     {fSymbol     = symbol;}
      void SetChromosome(const char *chromo) {fChromosome = chromo;}
      void SetCytoBand(const char *cytoband) {fCytoBand   = cytoband;}
      void SetSequence(const char *seq)      {fSequence   = seq;}

      const char* GetName()       const {return fName;}
      const char* GetSymbol()     const {return fSymbol;}
      const char* GetChromosome() const {return fChromosome;}
      const char* GetCytoBand()   const {return fCytoBand;}
      const char* GetSequence()   const {return fSequence;}

      ClassDef(XAnnotation,1) //Annotation
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XTransAnnotation                                                     //
//                                                                      //
// Class describing the transcript annotation of the units on GeneChips //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XTransAnnotation {

   private:
      Int_t         fUnitID;       //internal unit ID
      TString       fTranscriptID; //probe set ID, transcript_cluster_id
      TString       fName;         //name of unit: gene name
      TString       fSymbol;       //symbol of unit: gene symbol
      TString       fAccession;    //accession number of unit
      Int_t         fEntrezID;     //Entrez ID
      TString       fChromosome;   //chromosome: NA if not known or no gene 
      TString       fCytoBand;     //cytoband: NA if not known or no gene
      Int_t         fStart;        //start position on chromosome
      Int_t         fStop;         //stop position on chromosome
      char          fStrand;       //chromosome strand +/-

   public :
      XTransAnnotation()           {}
      virtual ~XTransAnnotation()  {}

      void SetUnitID(Int_t id)     {fUnitID = id;}
      void SetTranscriptID(const char *name) {fTranscriptID = name;}
      void SetName(const char *name)         {fName       = name;}
      void SetSymbol(const char *symbol)     {fSymbol     = symbol;}
      void SetAccession(const char *access)  {fAccession  = access;}
      void SetEntrezID(Int_t id)             {fEntrezID   = id;}
      void SetChromosome(const char *chromo) {fChromosome = chromo;}
      void SetCytoBand(const char *cytoband) {fCytoBand   = cytoband;}
      void SetStart(Int_t start)             {fStart      = start;}
      void SetStop(Int_t stop)               {fStop       = stop;}
      void SetStrand(char strand)            {fStrand     = strand;}

      Int_t       GetUnitID()       const {return fUnitID;}
      const char* GetTranscriptID() const {return fTranscriptID;}
      const char* GetName()         const {return fName;}
      const char* GetSymbol()       const {return fSymbol;}
      const char* GetAccession()    const {return fAccession;}
      Int_t       GetEntrezID()     const {return fEntrezID;}
      const char* GetChromosome()   const {return fChromosome;}
      const char* GetCytoBand()     const {return fCytoBand;}
      Int_t       GetStart()        const {return fStart;}
      Int_t       GetStop()         const {return fStop;}
      char        GetStrand()       const {return fStrand;}

      ClassDef(XTransAnnotation,1) //TransAnnotation
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeAnnotation                                                    //
//                                                                      //
// Class describing the annotation of the transcripts on the array      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XGenomeAnnotation: public XTransAnnotation {

   private:
      Short_t       fCrossHybType; //cross-hybridization type
      Short_t       fProbesetType; //probeset_type

   public :
      XGenomeAnnotation()          {}
      virtual ~XGenomeAnnotation() {}

      void SetCrossHybType(Short_t id) {fCrossHybType = id;}
      void SetProbesetType(Short_t id) {fProbesetType = id;}

      Short_t GetCrossHybType()  const {return fCrossHybType;}
      Short_t GetProbesetType()  const {return fProbesetType;}

      ClassDef(XGenomeAnnotation,1) //GenomeAnnotation
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonAnnotation                                                      //
//                                                                      //
// Class describing the annotation of the exons on the array            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XExonAnnotation: public XGenomeAnnotation {

   private:
      Int_t         fTranscriptID; //transcript_cluster_id (integer)
      Int_t         fExonID;       //exon_id

   public :
      XExonAnnotation()          {}
      virtual ~XExonAnnotation() {}

      void SetTranscriptID(Int_t id)  {fTranscriptID = id;}
      void SetExonID(Int_t id)        {fExonID = id;}

      Int_t   GetTranscriptID() const {return fTranscriptID;}
      Int_t   GetTranscriptIX() const {return fTranscriptID;}
      Int_t   GetExonID()       const {return fExonID;}

      ClassDef(XExonAnnotation,1) //ExonAnnotation
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProbesetAnnotation                                                  //
//                                                                      //
// Class describing the annotation of the probesets on the array        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class XProbesetAnnotation: public XExonAnnotation {

   private:
      Int_t         fProbesetID;   //probeset_id
      Int_t         fNProbes;      //probe_count
      Short_t       fLevelID;      //level, converted to id
      Short_t       fBounded;      //bounded and noBoundedEvidence
      Short_t       fFull;         //fl, i.e. flag for putative full-length mRNA

   public :
      XProbesetAnnotation()          {}
      virtual ~XProbesetAnnotation() {}

      void SetProbesetID(Int_t id)     {fProbesetID = id;}
      void SetNumProbes(Int_t id)      {fNProbes = id;}
      void SetLevelID(Short_t id)      {fLevelID = id;}
      void SetBounded(Short_t id)      {fBounded = id;}
      void SetFullFlag(Short_t fl)     {fFull = fl;}

      Int_t   GetProbesetID()    const {return fProbesetID;}
      Int_t   GetNumProbes()     const {return fNProbes;}
      Short_t GetLevelID()       const {return fLevelID;}
      Short_t GetBounded()       const {return fBounded;}
      Short_t GetFullFlag()      const {return fFull;}

      ClassDef(XProbesetAnnotation,1) //ProbesetAnnotation
};

#endif
