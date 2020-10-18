// File created: 05/18/2002                          last modified: 03/30/2013
// Author: Christian Stratowa 06/18/2000

/*
 *******************************************************************************
 *********************  XPS - eXpression Profiling System  *********************
 *******************************************************************************
 *
 *  Copyright (C) 2000-2013 Dr. Christian Stratowa
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
* Jun 2002 - Support to UpdateProbeInfo and UpdateAnnotation.
* Dec 2002 - Reorganize XPSSchemes.h to be dependent on XPSBase.h.
*          - Support to read  *.CDF files in XML format.
* Apr 2005 - Support to read binary *.CDF files.
* May 2005 - Add class XSNPChip to support SNPChips (Mapping arrays)
* Jan 2006 - To decrease size of saved trees replace Double_t with Double32_t
*            for fTm and Int_t with Short_t for fMask
* Mar 2006 - Add support for alternative schemes, i.e. CDF files.
* Apr 2006 - Add class XUnitTreeInfo.
* Aug 2006 - Add support for exon arrays, classes XExonChip, XLayout,
*            XLayoutTreeInfo, XExonScheme, XExonAnnotation.
* Oct 2006 - Expand support for exon arrays, add class XProbesetAnnotation
* Jan 2007 - Add support for new exon probeset annotation file (version 1.8)
* Aug 2007 - Add support for whole genome arrays, classes XGenomeChip,
*            XGenomeAnnotation
* Sep 2007 - Change fMask from class XScheme from "Short_t to "Int_t"
* Feb 2009 - Add support for whole genome array probeset annotation file (na27)
* Apr 2009 - Add support for both versions of probe info file "xxx_probe.tab"
* Dec 2012 - Protect against tabs in Affymetrix annotation files
*
******************************************************************************/

using namespace std;

//#ifndef ROOT_Varargs
#include "Varargs.h"
//#endif

#include <cstdlib>
#include <new>  //needed for new (nothrow)

#include "THashTable.h"

#include "XPSSchemes.h"


// Tree extensions and types for schemes:
// scheme trees
const char *kExtenScheme[] = { "scm", "cxy", "idx", "prb", "exn", "pbs", ""};
const char *kTypeScheme[]  = { "scheme",
                               "layout",
                               "unit",
                               "probe",
                               "exon",
                               "probeset",
                               ""};

// annotation trees
const char *kExtenAnnot[] = { "ann", "anx", "anp", "anc", ""};
const char *kTypeAnnot[]  = { "transcript",
                              "exon",
                              "probeset",
                              "control",
                              ""};


// Header line for probe info file "xxx_probe.tab": two versions
const Int_t kNPrbCols = 14;
const char *kPrbHeader[14] = { "Probe Set Name",               "Probe Set ID", 
                               "Serial Order",                 "serial order", 
                               "Probe X",                      "probe x", 
                               "Probe Y",                      "probe y", 
                               "Probe Interrogation Position", "probe interrogation position", 
                               "Probe Sequence",               "probe sequence", 
                               "Target Strandedness",          "target strandedness"};

// Header line for annotation file "xxx_annot4xps.txt"
const Int_t kNAnnCols = 6;
const char *kAnnHeader[6] = { "UnitName",
                              "GeneSymbol", 
                              "GeneName", 
                              "Chromosome", 
                              "CytoBand", 
                              "TargetSequence" };

// Header line for annotation file of GeneChip array
const Int_t kNAnnotCols = 43;
const char *kAnnotHeader[43] = {"Probe Set ID",
                                "GeneChip Array",
                                "Species Scientific Name",
                                "Annotation Date",
                                "Sequence Type",
                                "Sequence Source",
                                "Transcript ID(Array Design)",
                                "Target Description",
                                "Representative Public ID",
                                "Archival UniGene Cluster",
                                "UniGene ID",
                                "Genome Version",
                                "Alignments",
                                "Gene Title",
                                "Gene Symbol",
                                "Chromosomal Location",
                                "Unigene Cluster Type",
                                "Ensembl",
                                "Entrez Gene",
                                "SwissProt",
                                "EC",
                                "OMIM",
                                "RefSeq Protein ID",
                                "RefSeq Transcript ID",
                                "FlyBase",
                                "AGI",
                                "WormBase",
                                "MGI Name",
                                "RGD Name",
                                "SGD accession number",
                                "Gene Ontology Biological Process",
                                "Gene Ontology Cellular Component",
                                "Gene Ontology Molecular Function",
                                "Pathway",
                                "Protein Families",
                                "Protein Domains",
                                "InterPro",
                                "Trans Membrane",
                                "QTL",
                                "Annotation Description",
                                "Annotation Transcript Cluster",
                                "Transcript Assignments",
                                "Annotation Notes" };


// Header line for annotation file of exon array
const Int_t kNProbesetCols = 39;
const char *kProbesetHeader[39] = {"probeset_id",
                                   "seqname",
                                   "strand",
                                   "start",
                                   "stop",
                                   "probe_count",
                                   "transcript_cluster_id",
                                   "exon_id",
                                   "psr_id",
                                   "gene_assignment",
                                   "mrna_assignment",
//version 1.6                                   "probeset_type",
                                   "crosshyb_type",
                                   "number_independent_probes",
                                   "number_cross_hyb_probes",
                                   "number_nonoverlapping_probes",
                                   "level",
                                   "bounded",
                                   "noBoundedEvidence",
                                   "has_cds",
                                   "fl",
                                   "mrna",
                                   "est",
                                   "vegaGene",
                                   "vegaPseudoGene",
                                   "ensGene",
                                   "sgpGene",
                                   "exoniphy",
                                   "twinscan",
                                   "geneid",
                                   "genscan",
                                   "genscanSubopt",
                                   "mouse_fl",
                                   "mouse_mrna",
                                   "rat_fl",
                                   "rat_mrna",
                                   "microRNAregistry",
                                   "rnaGene",
                                   "mitomap",
                                   "probeset_type" };

// Header line for transcript annotation file of genome/exon arrays
const Int_t kNTranscriptCols = 18;
const char *kTranscriptHeader[18] = {"transcript_cluster_id",
                                     "probeset_id",
                                     "seqname",
                                     "strand",
                                     "start",
                                     "stop",
                                     "total_probes",
                                     "gene_assignment",
                                     "mrna_assignment",
                                     "swissprot",
                                     "unigene",
                                     "GO_biological_process",
                                     "GO_cellular_component",
                                     "GO_molecular_function",
                                     "pathway",
                                     "protein_domains",
                                     "crosshyb_type",
                                     "category" };


// Header line for control annotation file of exon array
const Int_t kNControlCols = 3;
const char *kControlHeader[3] = {"probeset_id",
                                 "group_name",
                                 "probeset_name" };


ClassImp(XSchemeManager);
ClassImp(XLayoutTreeInfo);
ClassImp(XProbeTreeInfo);
ClassImp(XSchemeTreeInfo);
ClassImp(XUnitTreeInfo);
ClassImp(XGenomeTreeInfo);
ClassImp(XDNAChip);
ClassImp(XMicroArray);
ClassImp(XOligoArray);
ClassImp(XGeneChip);
ClassImp(XSNPChip);
ClassImp(XGenomeChip);
ClassImp(XExonChip);
ClassImp(XPosition);
ClassImp(XLayout);
ClassImp(XScheme);
ClassImp(XGCScheme);
ClassImp(XExonScheme);
ClassImp(XProbe);
ClassImp(XGCProbe);
ClassImp(XUnitID);
ClassImp(XUnit);
ClassImp(XGCUnit);
ClassImp(XExonUnit);
ClassImp(XAnnotation);
ClassImp(XTransAnnotation);
ClassImp(XGenomeAnnotation);
ClassImp(XExonAnnotation);
ClassImp(XProbesetAnnotation);

const Int_t kCDFMagicNumber   = 67;
const Int_t kCDFVersionNumber = 1;
const Int_t kUnitNameLength   = 64; //maximal unit name length
const Int_t kUnitNameLength1  = 65; //maximal unit name length + 1

//version number of annotation file(s)
const Int_t kVersionAnnot22 = 22;

//buffer to read annotation file
const Int_t kAnnBuf = 256000;

//separators for multipart entries
const char *kSepSl2 = " // ";
const char *kSepSl3 = " /// ";
const Int_t kNumSl2 = 4;
const Int_t kNumSl3 = 5;
const Int_t kNumAlg = 3;    //maximum number of alignment substrings separated by kSepSl2
const Int_t kNumSub = 1000; //maximum number of alignment substrings separated by kSepSl3

//debug: print function names
const Bool_t kCS  = 0; 
const Bool_t kCSa = 0; //debug: print function names in loops


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSchemeManager                                                       //
//                                                                      //
// Class for managing import of microarray scheme definitions           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSchemeManager::XSchemeManager()
               :XManager()
{
   // Default SchemeManager constructor
   if(kCS) cout << "---XSchemeManager::XSchemeManager(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSchemeManager::XSchemeManager(const char *name, const char *title, Int_t verbose)
               :XManager(name, title, verbose)
{
   // Normal XSchemeManager constructor
   if(kCS) cout << "---XSchemeManager::XSchemeManager------" << endl;

}//Constructor

//______________________________________________________________________________
XSchemeManager::~XSchemeManager()
{
   // SchemeManager destructor
   if(kCS) cout << "---XSchemeManager::~XSchemeManager------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XSchemeManager::ExportTrees(const char *exten, const char * /*varlist*/, 
                      ofstream &/*output*/, const char * /*sep*/)
{
   // Export variables from varlist for all trees in current file
   // with extension exten
   // Note: Not applicable for scheme trees
   if(kCS) cout << "------XSchemeManager::ExportTrees------" << endl;

   cout << "Note: The option <*." << exten << "> is not available for schemes" << endl;
   cout << "      since scheme trees have different number of entries." << endl;

   return errNoErr;
}//ExportTrees

//______________________________________________________________________________
Int_t XSchemeManager::NewScheme(const char *name, const char *infile,
                      const char *type, const char *sep, char delim, Int_t split)
{
   // Read new scheme from infile. Scheme can exist of following types of infile:
   // type = layout: infile containing (x,y) coordinates for probes (optional)
   // type = scheme: infile containing scheme information (optional including layout)
   // type = probe:  infile containing information about probes on array (optional)
   if(kCS) cout << "------XSchemeManager::NewScheme------" << endl;

   if (fAbort) return errAbort;

   TString exten  = Type2Extension(type, kTypeScheme, kExtenScheme);
   TString option = TString("CREATE.") + exten;

// Abort if scheme tree does already exist
   TTree  *stree = 0; 
   TString sname = Path2Name(name, dSEP, ".") + "." + exten;
   stree = (TTree*)fFile->Get(sname); 
   if (stree) {
      cerr << "Error: Scheme tree <" << sname.Data() << "> does already exist."
           << endl;
      return errGetTree;
   }//if

   return Import(name, infile, name, option, sep, delim, split);
}//NewScheme

//______________________________________________________________________________
Int_t XSchemeManager::UpdateScheme(const char *name, const char *infile, 
                      const char *type, const char *sep, char delim, Int_t split)
{
   // Update scheme file with name with data from infile of type(s):
   // type = layout: infile containing (x,y) coordinates for probes (optional)
   // type = scheme: infile containing scheme information (optional including layout)
   // type = probe:  infile containing information about probes on array (optional)
   if(kCS) cout << "------XSchemeManager::UpdateScheme------" << endl;

   if (fAbort) return errAbort;

   TString exten  = Type2Extension(type, kTypeScheme, kExtenScheme);
   TString option = TString("UPDATE.") + exten;

   return Import(name, infile, name, option, sep, delim, split);
}//UpdateScheme

//______________________________________________________________________________
void XSchemeManager::DeleteScheme(const char *name)
{
   // Delete scheme with name from file
   if(kCS) cout << "------XSchemeManager::DeleteScheme------" << endl;

   DeleteTreeSet(name);
}//DeleteScheme

//______________________________________________________________________________
Int_t XSchemeManager::NewAnnotation(const char *name, const char *infile,
                      const char *type, const char *sep, char delim, Int_t split)
{
   // Read the probe annotation from infile and store the information in new tree
   // The infile can have annotation of the following type: 
   // type = transcript: infile contains transcript annotations (default)
   // type = exon:       infile contains exon annotations (for exon arrays only)
   // type = probeset:   infile contains probeset annotations (exon arrays only)
   // type = control:    infile contains AFFX annotations (for exon arrays only)
//TO DO: EV???:
   // type = xps:        infile contains annotation file in xps-format
//                         (e.g necessary for annotation of alternative CDFs!!!)

   if(kCS) cout << "------XSchemeManager::NewAnnotation------" << endl;

   if (fAbort) return errAbort;
//??   if (fInterrupt) return errGeneral;

   TString exten  = Type2Extension(type, kTypeAnnot, kExtenAnnot);
   TString option = TString("CREATE.") + exten;

// Abort if annotation tree does already exist
   TTree  *atree = 0; 
   TString aname = Path2Name(name, dSEP, ".") + "." + exten;
   atree = (TTree*)fFile->Get(aname); 
   if (atree) {
      cerr << "Error: Scheme tree <" << aname.Data() << "> does already exist."
           << endl;
      return errGetTree;
   }//if

   return Import(name, infile, name, option, sep, delim, split);
}//NewAnnotation

//______________________________________________________________________________
Int_t XSchemeManager::UpdateAnnotation(const char *name, const char *infile,
                      const char *type, const char *sep, char delim, Int_t split)
{
   // Read the probe annotation from infile and store the information in tree
   // The infile can have annotation of the following type: 
   // type = transcript: infile contains transcript annotations (default)
   // type = exon:       infile contains exon annotations (for exon arrays only)
   // type = probeset:   infile contains probeset annotations (exon arrays only)
   // type = control:    infile contains AFFX annotations (for exon arrays only)
//TO DO: EV???:
   // type = xps:        infile contains annotation file in xps-format
//                         (e.g necessary for annotation of alternative CDFs!!!)
   if(kCS) cout << "------XSchemeManager::UpdateAnnotation------" << endl;

   if (fAbort) return errAbort;
//   if (fInterrupt) return errGeneral;

   TString exten  = Type2Extension(type, kTypeAnnot, kExtenAnnot);
   TString option = TString("UPDATE.") + exten;

   return Import(name, infile, name, option, sep, delim, split);
}//UpdateAnnotation

//______________________________________________________________________________
Int_t XSchemeManager::New(const char *name, const char *dir, const char *type,
                      const char *data, Option_t *option)
{
   // Create new root file with name in directory dir for type and data type
   if(kCS) cout << "------XSchemeManager::New------" << endl;

// Close existing file
   if (fFile) {
      if (XManager::fgVerbose) {
         cout << "Closing existing file <" << fFile->GetName() << ">..." << endl;
      }//if
      Close();  //??
   }//if

   return XManager::New(name, dir, type, data, option);
}//New

//______________________________________________________________________________
Int_t XSchemeManager::HandleError(Int_t err, const char *name1, const char *name2)
{
   // Handle error messages
   if(kCS) cout << "------XSchemeManager::HandleError------" << endl;

   switch (err) {
      case errExtension:
         cerr << "Error: Tree(s) with extension <" << name1 << "> not known" << endl;
         fInterrupt = kTRUE;
         break;

      case errCDFVersion:
         cerr << "Error: CDF with version/magic number <" << name1 
              << "> is not supported." << endl;
         fAbort = kTRUE;
         break;

      case errAnnVersion:
         cerr << "Error: Wrong version <" << name1 << "> of annotation file." << endl;
         fAbort = kTRUE;
         break;

      default:
         err = XManager::HandleError(err, name1, name2);
         break;
   }//switch

   return err;
}//HandleError

//______________________________________________________________________________
XTreeSet *XSchemeManager::NewTreeSet(const char *type)
{
   // Create new tree set of type:
   // "MicroArray", "OligoArray", "GeneChip", "SNPChip", "ExonChip"
   if(kCS) cout << "------XSchemeManager::NewTreeSet------" << endl;

   XTreeSet *set = 0;
   if (strcmp(type,"GeneChip") == 0) {
      set = new XGeneChip("untitled", type);
   } else if (strcmp(type,"SNPChip") == 0) {
      set = new XSNPChip("untitled", type);
   } else if (strcmp(type,"GenomeChip") == 0) {
      set = new XGenomeChip("untitled", type);
   } else if (strcmp(type,"ExonChip") == 0) {
      set = new XExonChip("untitled", type);
   } else if (strcmp(type,"OligoArray") == 0) {
      set = new XOligoArray("untitled", type);
   } else if (strcmp(type,"MicroArray") == 0) {
      set = new XMicroArray("untitled", type);
   } else {
      cerr << "Error: chip type <" << type << "> not known" << endl;
   }//if

   return set;
}//NewTreeSet


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XLayoutTreeInfo                                                      //
//                                                                      //
// Class containing info about layout tree stored in fUserInfo of tree  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XLayoutTreeInfo::XLayoutTreeInfo() 
                :XTreeInfo()
{
   // Default LayoutTreeInfo constructor
   if(kCS) cout << "---XLayoutTreeInfo::XLayoutTreeInfo(default)------" << endl;

   fNRows      = 0; 
   fNCols      = 0;
   fSequential = 1;
}//Constructor

//______________________________________________________________________________
XLayoutTreeInfo::XLayoutTreeInfo(const char *name, const char *title) 
                :XTreeInfo(name, title)
{
   // Normal LayoutTreeInfo constructor
   if(kCS) cout << "---XLayoutTreeInfo::XLayoutTreeInfo------" << endl;

   fNRows      = 0; 
   fNCols      = 0;
   fSequential = 1;
}//Constructor

//______________________________________________________________________________
XLayoutTreeInfo::~XLayoutTreeInfo()
{
   // LayoutTreeInfo destructor
   if(kCS) cout << "---XLayoutTreeInfo::~XLayoutTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XLayoutTreeInfo::AddUserInfo(XTreeSet *set)
{
   // Add user info from tree set
   if(kCS) cout << "------XLayoutTreeInfo::AddUserInfo------" << endl;

   XDNAChip *dcset = (XDNAChip*)set;
   fNRows      = dcset->GetNumRows();
   fNCols      = dcset->GetNumColumns();
   fSequential = dcset->GetSequential();
}//AddUserInfo

//______________________________________________________________________________
void XLayoutTreeInfo::AddUserInfo(Int_t nrows, Int_t ncols, Short_t seq)
{
   // Add user info from tree set
   if(kCS) cout << "------XLayoutTreeInfo::AddUserInfo------" << endl;

   fNRows      = nrows;
   fNCols      = ncols;
   fSequential = seq;
}//AddUserInfo

//______________________________________________________________________________
Double_t XLayoutTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCSa) cout << "------XLayoutTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNRows") == 0) {
      return fNRows;
   } else if (strcmp(name, "fNCols") == 0) {
      return fNCols;
   } else if (strcmp(name, "fSequential") == 0) {
      return fSequential;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XProbeTreeInfo                                                       //
//                                                                      //
// Class containing info about probe tree stored in fUserInfo of tree   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XProbeTreeInfo::XProbeTreeInfo() 
               :XTreeInfo()
{
   // Default ProbeTreeInfo constructor
   if(kCS) cout << "---XProbeTreeInfo::XProbeTreeInfo(default)------" << endl;

   fMinGC = -1;
   fMaxGC = -1;
   fMinTm = -1.0;
   fMaxTm = -1.0;
}//Constructor

//______________________________________________________________________________
XProbeTreeInfo::XProbeTreeInfo(const char *name, const char *title) 
               :XTreeInfo(name, title)
{
   // Normal ProbeTreeInfo constructor
   if(kCS) cout << "---XProbeTreeInfo::XProbeTreeInfo------" << endl;

   fMinGC = -1;
   fMaxGC = -1;
   fMinTm = -1.0;
   fMaxTm = -1.0;
}//Constructor

//______________________________________________________________________________
XProbeTreeInfo::~XProbeTreeInfo()
{
   // ProbeTreeInfo destructor
   if(kCS) cout << "---XProbeTreeInfo::~XProbeTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XProbeTreeInfo::AddUserInfo(XTreeSet *set)
{
   // Add user info from tree set
   if(kCS) cout << "------XProbeTreeInfo::AddUserInfo------" << endl;

   XGeneChip *gcset = (XGeneChip*)set;
   fMinGC = gcset->GetMinContentGC();
   fMaxGC = gcset->GetMaxContentGC();
   fMinTm = gcset->GetMinMeltingT();
   fMaxTm = gcset->GetMaxMeltingT();
}//AddUserInfo

//______________________________________________________________________________
void XProbeTreeInfo::AddUserInfo(Int_t minGC, Int_t maxGC, Double_t minTm,
                     Double_t maxTm)
{
   // Add user info from tree set
   if(kCS) cout << "------XProbeTreeInfo::AddUserInfo------" << endl;

   fMinGC = minGC;
   fMaxGC = maxGC;
   fMinTm = minTm;
   fMaxTm = maxTm;
}//AddUserInfo

//______________________________________________________________________________
Double_t XProbeTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCSa) cout << "------XProbeTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fMinGC") == 0) {
      return fMinGC;
   } else if (strcmp(name, "fMaxGC") == 0) {
      return fMaxGC;
   } else if (strcmp(name, "fMinTm") == 0) {
      return fMinTm;
   } else if (strcmp(name, "fMaxTm") == 0) {
      return fMaxTm;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSchemeTreeInfo                                                      //
//                                                                      //
// Class containing info about scheme tree stored in fUserInfo of tree  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSchemeTreeInfo::XSchemeTreeInfo() 
                :XTreeInfo()
{
   // Default SchemeTreeInfo constructor
   if(kCS) cout << "---XSchemeTreeInfo::XSchemeTreeInfo(default)------" << endl;

   fNRows     = 0; 
   fNCols     = 0;
   fNProbes   = 0;
   fNControls = 0;
   fNGenes    = 0;
   fNUnits    = 0;
}//Constructor

//______________________________________________________________________________
XSchemeTreeInfo::XSchemeTreeInfo(const char *name, const char *title) 
                :XTreeInfo(name, title)
{
   // Normal SchemeTreeInfo constructor
   if(kCS) cout << "---XSchemeTreeInfo::XSchemeTreeInfo------" << endl;

   fNRows     = 0; 
   fNCols     = 0;
   fNProbes   = 0;
   fNControls = 0;
   fNGenes    = 0;
   fNUnits    = 0;
}//Constructor

//______________________________________________________________________________
XSchemeTreeInfo::~XSchemeTreeInfo()
{
   // SchemeTreeInfo destructor
   if(kCS) cout << "---XSchemeTreeInfo::~XSchemeTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XSchemeTreeInfo::AddUserInfo(XTreeSet *set)
{
   // Add user info from tree set
   if(kCS) cout << "------XSchemeTreeInfo::AddUserInfo------" << endl;

   XDNAChip *dcset = (XDNAChip*)set;
   fNRows     = dcset->GetNumRows();
   fNCols     = dcset->GetNumColumns();
   fNProbes   = dcset->GetNumProbes();
   fNControls = dcset->GetNumControls();
   fNGenes    = dcset->GetNumGenes();
   fNUnits    = dcset->GetNumUnits();
}//AddUserInfo

//______________________________________________________________________________
void XSchemeTreeInfo::AddUserInfo(Int_t nrows, Int_t ncols, Int_t nprobes,
                      Int_t nctrls, Int_t ngenes, Int_t nunits)
{
   // Add user info from tree set
   if(kCS) cout << "------XSchemeTreeInfo::AddUserInfo------" << endl;

   fNRows     = nrows;
   fNCols     = ncols;
   fNProbes   = nprobes;
   fNControls = nctrls;
   fNGenes    = ngenes;
   fNUnits    = nunits;
}//AddUserInfo

//______________________________________________________________________________
Double_t XSchemeTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCSa) cout << "------XSchemeTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNRows") == 0) {
      return fNRows;
   } else if (strcmp(name, "fNCols") == 0) {
      return fNCols;
   } else if (strcmp(name, "fNProbes") == 0) {
      return fNProbes;
   } else if (strcmp(name, "fNControls") == 0) {
      return fNControls;
   } else if (strcmp(name, "fNGenes") == 0) {
      return fNGenes;
   } else if (strcmp(name, "fNUnits") == 0) {
      return fNUnits;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XUnitTreeInfo                                                        //
//                                                                      //
// Class containing info about unit tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XUnitTreeInfo::XUnitTreeInfo() 
              :XTreeInfo()
{
   // Default UnitTreeInfo constructor
   if(kCS) cout << "---XUnitTreeInfo::XUnitTreeInfo(default)------" << endl;

   fNControls = 0;
   fNAffx     = 0;
   fNGenes    = 0;
   fMinNPairs = 0; 
   fMaxNPairs = 0;
   fMinNCells = 0; 
   fMaxNCells = 0;
}//Constructor

//______________________________________________________________________________
XUnitTreeInfo::XUnitTreeInfo(const char *name, const char *title) 
              :XTreeInfo(name, title)
{
   // Normal UnitTreeInfo constructor
   if(kCS) cout << "---XUnitTreeInfo::XUnitTreeInfo------" << endl;

   fNControls = 0;
   fNAffx     = 0;
   fNGenes    = 0;
   fMinNPairs = 0; 
   fMaxNPairs = 0;
   fMinNCells = 0; 
   fMaxNCells = 0;
}//Constructor

//______________________________________________________________________________
XUnitTreeInfo::~XUnitTreeInfo()
{
   // UnitTreeInfo destructor
   if(kCS) cout << "---XUnitTreeInfo::~XUnitTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XUnitTreeInfo::AddUserInfo(XTreeSet *set)
{
   // Add user info from tree set
   if(kCS) cout << "------XUnitTreeInfo::AddUserInfo------" << endl;

   XGeneChip *gcset = (XGeneChip*)set;
   fNControls = gcset->GetNumControls();
   fNAffx     = gcset->GetNumAffx();
   fNGenes    = gcset->GetNumGenes();

// not available
   fMinNPairs = -1;
   fMaxNPairs = -1;
   fMinNCells = -1;
   fMaxNCells = -1;
}//AddUserInfo

//______________________________________________________________________________
void XUnitTreeInfo::AddUserInfo(Int_t nctrls, Int_t naffx, Int_t ngenes,
                    Int_t nmin, Int_t min, Int_t nmax, Int_t max)
{
   // Add user info from tree set
   if(kCS) cout << "------XUnitTreeInfo::AddUserInfo------" << endl;

   fNControls = nctrls;
   fNAffx     = naffx;
   fNGenes    = ngenes;
   fMinNPairs = min;
   fMaxNPairs = max;
   fMinNCells = nmin;
   fMaxNCells = nmax;
}//AddUserInfo

//______________________________________________________________________________
Double_t XUnitTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCSa) cout << "------XUnitTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNControls")        == 0) {
      return fNControls;
   } else if (strcmp(name, "fNAffx")     == 0) {
      return fNAffx;
   } else if (strcmp(name, "fNGenes")    == 0) {
      return fNGenes;
   } else if (strcmp(name, "fMinNPairs") == 0) {
      return fMinNPairs;
   } else if (strcmp(name, "fMaxNPairs") == 0) {
      return fMaxNPairs;
   } else if (strcmp(name, "fMinNCells") == 0) {
      return fMinNCells;
   } else if (strcmp(name, "fMaxNCells") == 0) {
      return fMaxNCells;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeTreeInfo                                                      //
//                                                                      //
// Class containing info about exon tree stored in fUserInfo of tree    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenomeTreeInfo::XGenomeTreeInfo() 
                :XTreeInfo()
{
   // Default GenomeTreeInfo constructor
   if(kCS) cout << "---XGenomeTreeInfo::XGenomeTreeInfo(default)------" << endl;

   fNControls = 0;
   fNAffx     = 0;
   fNSubunits = 0;
   fMinNCells = 0;
   fMaxNCells = 0;
}//Constructor

//______________________________________________________________________________
XGenomeTreeInfo::XGenomeTreeInfo(const char *name, const char *title) 
                :XTreeInfo(name, title)
{
   // Normal GenomeTreeInfo constructor
   if(kCS) cout << "---XGenomeTreeInfo::XGenomeTreeInfo------" << endl;

   fNControls = 0;
   fNAffx     = 0;
   fNSubunits = 0;
   fMinNCells = 0;
   fMaxNCells = 0;
}//Constructor

//______________________________________________________________________________
XGenomeTreeInfo::~XGenomeTreeInfo()
{
   // GenomeTreeInfo destructor
   if(kCS) cout << "---XGenomeTreeInfo::~XGenomeTreeInfo------" << endl;

}//Destructor

//______________________________________________________________________________
void XGenomeTreeInfo::AddUserInfo(XTreeSet *set)
{
   // Add user info from tree set
   if(kCS) cout << "------XGenomeTreeInfo::AddUserInfo------" << endl;

   XGenomeChip *ecset = (XGenomeChip*)set;
   fNControls = ecset->GetNumControls();
   fNAffx     = ecset->GetNumAffx();
//???
   fNSubunits = ecset->GetNumUnits();

// not available
   fMinNCells = -1;
   fMaxNCells = -1;
}//AddUserInfo

//______________________________________________________________________________
void XGenomeTreeInfo::AddUserInfo(Int_t nctrls, Int_t naffx, Int_t nsubs,
                                  Int_t mincells, Int_t maxcells)
{
   // Add user info from tree set
   if(kCS) cout << "------XGenomeTreeInfo::AddUserInfo------" << endl;

   fNControls = nctrls;
   fNAffx     = naffx;
   fNSubunits = nsubs;
   fMinNCells = mincells;
   fMaxNCells = maxcells;
}//AddUserInfo

//______________________________________________________________________________
Double_t XGenomeTreeInfo::GetValue(const char *name)
{
   // Return value for class member field name
   if(kCSa) cout << "------XGenomeTreeInfo::GetValue------" << endl;

   if (strcmp(name, "fNControls") == 0) {
      return fNControls;
   } else if (strcmp(name, "fNAffx") == 0) {
      return fNAffx;
   } else if (strcmp(name, "fNSubunits") == 0) {
      return fNSubunits;
   } else if (strcmp(name, "fMinNCells") == 0) {
      return fMinNCells;
   } else if (strcmp(name, "fMaxNCells") == 0) {
      return fMaxNCells;
   }//if
   return 0;
}//GetValue


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XDNAChip                                                             //
//                                                                      //
// Base class for DNA-chips                                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XDNAChip::XDNAChip()
         :XTreeSet()
{
   // Default DNAChip constructor
   if(kCS) cout << "---XDNAChip::XDNAChip(default)------" << endl;

   fCxyTreeName = "";
   fScmTreeName = "";
   fIdxTreeName = "";
   fPrbTreeName = "";
   fAnnTreeName = "";
}//Constructor

//______________________________________________________________________________
XDNAChip::XDNAChip(const char *name, const char *title)
         :XTreeSet(name, title)
{
   // Normal DNAChip constructor
   if(kCS) cout << "---XDNAChip::XDNAChip------" << endl;

   fCxyTreeName  = "NA";
   fScmTreeName  = "NA";
   fIdxTreeName  = "NA";
   fPrbTreeName  = "NA";
   fAnnTreeName  = "NA";
   fVersionAnnot = "NA";
   fNRows        = 0;
   fNCols        = 0;
   fNControls    = 0;
   fNGenes       = 0;
   fNUnits       = 0;
   fSequential   = 1;
}//Constructor

//______________________________________________________________________________
XDNAChip::XDNAChip(const XDNAChip &chip) 
         :XTreeSet(chip)
{
   // DNAChip copy constructor
   if(kCS) cout << "---XDNAChip::XDNAChip(copy)------" << endl;

   fCxyTreeName  = chip.fCxyTreeName.Copy();
   fScmTreeName  = chip.fScmTreeName.Copy();
   fIdxTreeName  = chip.fIdxTreeName.Copy();
   fPrbTreeName  = chip.fPrbTreeName.Copy();
   fAnnTreeName  = chip.fAnnTreeName.Copy();
   fVersionAnnot = chip.fVersionAnnot.Copy();
   fNRows        = chip.fNRows;
   fNCols        = chip.fNCols;
   fNControls    = chip.fNControls;
   fNGenes       = chip.fNGenes;
   fNUnits       = chip.fNUnits;
   fSequential   = chip.fSequential;
}//CopyConstructor

//______________________________________________________________________________
XDNAChip& XDNAChip::operator=(const XDNAChip& rhs)
{
   // DNAChip assignment operator.

   if (this != &rhs) {
      XTreeSet::operator=(rhs);

      fCxyTreeName  = rhs.fCxyTreeName.Copy();
      fScmTreeName  = rhs.fScmTreeName.Copy();
      fIdxTreeName  = rhs.fIdxTreeName.Copy();
      fPrbTreeName  = rhs.fPrbTreeName.Copy();
      fAnnTreeName  = rhs.fAnnTreeName.Copy();
      fVersionAnnot = rhs.fVersionAnnot.Copy();
      fNRows        = rhs.fNRows;
      fNCols        = rhs.fNCols;
      fNControls    = rhs.fNControls;
      fNGenes       = rhs.fNGenes;
      fNUnits       = rhs.fNUnits;
      fSequential   = rhs.fSequential;
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XDNAChip::~XDNAChip()
{
   // DNAChip destructor
   if(kCS) cout << "---XDNAChip::~XDNAChip------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XDNAChip::ExportTreeType(const char *exten, Int_t n, TString *names, 
                const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output
   if(kCS) cout << "------XDNAChip::ExportTreeType------" << endl;

   if ((strcmp(exten, kExtenScheme[2]) == 0) ||
       (strcmp(exten, kExtenScheme[4]) == 0) ||
       (strcmp(exten, kExtenScheme[5]) == 0)) {
      return this->ExportUnitTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenScheme[1]) == 0) {
      return this->ExportLayoutTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenScheme[0]) == 0) {
      return this->ExportSchemeTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenScheme[3]) == 0) {
      return this->ExportProbeTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[0]) == 0) {
      return this->ExportTransAnnotTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[1]) == 0) {
      return this->ExportExonAnnotTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[2]) == 0) {
      return this->ExportProbesetAnnotTree(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[3]) == 0) {
      return this->ExportControlAnnotTree(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//ExportTreeType

//______________________________________________________________________________
Int_t XDNAChip::ExportTreeXML(const char *exten, Int_t n, TString *names, 
                const char *varlist, ofstream &output, const char *sep)
{
   // Export data stored in tree treename to file output as XML-file
   if(kCS) cout << "------XDNAChip::ExportTreeXML------" << endl;

   if ((strcmp(exten, kExtenScheme[2]) == 0) ||
       (strcmp(exten, kExtenScheme[4]) == 0) ||
       (strcmp(exten, kExtenScheme[5]) == 0)) {
      return this->ExportUnitXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenScheme[1]) == 0) {
      return this->ExportLayoutXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenScheme[0]) == 0) {
      return this->ExportSchemeXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenScheme[3]) == 0) {
      return this->ExportProbeXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[0]) == 0) {
      return this->ExportTransAnnotXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[1]) == 0) {
      return this->ExportExonAnnotXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[2]) == 0) {
      return this->ExportProbesetAnnotXML(n, names, varlist, output, sep);
   } else if (strcmp(exten, kExtenAnnot[3]) == 0) {
      return this->ExportControlAnnotXML(n, names, varlist, output, sep);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//ExportTreeXML

//______________________________________________________________________________
Int_t XDNAChip::Import(ifstream &input, Option_t *option, const char *sep,
                char delim, Int_t split)
{
   // Import data from input. Option consists of "option.type" with type
   // determining type of tree, e.g. "CREATE.scm" imports new scheme
   // If option contains "UPDATE" replace existing data
   if(kCS) cout << "------XDNAChip::Import------" << endl;

   TString exten = Path2Name(option,".","");
   if (strcmp(exten.Data(), kExtenScheme[0]) == 0) {           //scheme
     return this->ImportScheme(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenScheme[1]) == 0) {    //layout
      return this->ImportLayout(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenScheme[3]) == 0) {    //probe
      return this->ImportProbeInfo(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[0]) == 0) {     //transcript
      return this->ImportTransAnnotation(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[1]) == 0) {     //exon
      return this->ImportExonAnnotation(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[2]) == 0) {     //probeset
      return this->ImportProbesetAnnotation(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[3]) == 0) {     //control
      return this->ImportControlAnnotation(input, option, sep, delim, split);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//Import

//______________________________________________________________________________
Int_t XDNAChip::ImportXML(ifstream &input, Option_t *option, const char *sep,
                char delim, Int_t split)
{
   // Import XML-data from input. Option consists of "option.type" with type
   // determining type of tree, e.g. "CREATE.scm" imports new scheme
   // If option contains "UPDATE" replace existing data
   if(kCS) cout << "------XDNAChip::ImportXML------" << endl;

   TString exten = Path2Name(option,".","");
   if (strcmp(exten.Data(), kExtenScheme[0]) == 0) {
      return this->ImportSchemeXML(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenScheme[1]) == 0) {
      return this->ImportLayoutXML(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenScheme[3]) == 0) {
      return this->ImportProbeInfoXML(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[0]) == 0) {
      return this->ImportTransAnnotationXML(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[1]) == 0) {
      return this->ImportExonAnnotationXML(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[2]) == 0) {
      return this->ImportProbesetAnnotationXML(input, option, sep, delim, split);
   } else if (strcmp(exten.Data(), kExtenAnnot[3]) == 0) {
      return this->ImportControlAnnotationXML(input, option, sep, delim, split);
   } else {
      return fManager->HandleError(errExtension, exten);
   }//if
}//ImportXML

//______________________________________________________________________________
Int_t XDNAChip::ExportTrees(const char *exten, const char * /*varlist*/, 
                ofstream &/*output*/, const char * /*sep*/)
{
   // Export variables from varlist for all trees in current set
   // with extension exten
   if(kCS) cout << "------XDNAChip::ExportTrees------" << endl;

   cout << "Note: Option <*." << exten << "> is not available for schemes" << endl;
   cout << "      since scheme trees have different number of entries." << endl;

   return errNoErr;
}//ExportTrees

//______________________________________________________________________________
Int_t XDNAChip::ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output
   // varlist can be at most: "fProbeLen:fMask" or "*"
   if(kCS) cout << "------XDNAChip::ExportSchemeTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasPLen = kFALSE;
   Bool_t hasMask = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasPLen = kTRUE;
      hasMask = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fProbeLen")   == 0) {hasPLen = kTRUE;}
         if (strcmp(name,"fMask")       == 0) {hasMask = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get scheme tree for this chip
   XScheme *scheme = 0;
   fTree = (TTree*)gDirectory->Get((names[0]).Data());
   if (!fTree) return errGetTree;
   fTree->SetBranchAddress("ScmBranch", &scheme);

   Int_t nentries = (Int_t)(fTree->GetEntries());
//??   Int_t size     = fNRows*fNCols;
   Int_t size     = fNProbes;
   if (nentries != size) {
      TString str = ""; str += size;
//??      return fManager->HandleError(errNumTreeEntries, fTree->GetName(), str);
   }//if

// Output header
   output << "UNIT_ID" << sep << "X" << sep << "Y";
   if (hasPLen) output << sep << "ProbeLength";
   if (hasMask) output << sep << "Mask";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << scheme->GetUnitID() << sep
             << scheme->GetX() << sep << scheme->GetY();
      if (hasPLen) output << sep << scheme->GetProbeLength();
      if (hasMask) output << sep << scheme->GetMask();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportSchemeTree

//______________________________________________________________________________
Int_t XDNAChip::ExportTransAnnotTree(Int_t n, TString *names, const char *varlist,
                ofstream &output, const char *sep)
{
   // Export data stored in annotation tree to file output
   if(kCS) cout << "------XDNAChip::ExportTransAnnotTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasUnit = kFALSE;
   Bool_t hasGene = kFALSE;
   Bool_t hasSymb = kFALSE;
   Bool_t hasChro = kFALSE;
   Bool_t hasCyto = kFALSE;
   Bool_t hasSequ = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasUnit = kTRUE;
      hasGene = kTRUE;
      hasSymb = kTRUE;
      hasChro = kTRUE;
      hasCyto = kTRUE;
      hasSequ = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName")   == 0) {hasUnit = kTRUE;}
         if (strcmp(name,"fName")       == 0) {hasGene = kTRUE;}
         if (strcmp(name,"fSymbol")     == 0) {hasSymb = kTRUE;}
         if (strcmp(name,"fChromosome") == 0) {hasChro = kTRUE;}
         if (strcmp(name,"fCytoBand")   == 0) {hasCyto = kTRUE;}
         if (strcmp(name,"fSequence")   == 0) {hasSequ = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get unit tree for this chip
   TString unitname = fName + "." + TString(kExtenScheme[2]);
   XUnit *unit = 0;
   TTree *unittree = (TTree*)gDirectory->Get(unitname); 
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch", &unit);

// Get annotation tree
   XAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get((names[0].Data())); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

   Int_t unitentries  = (Int_t)(unittree->GetEntries());
   Int_t annotentries = (Int_t)(anntree->GetEntries());
   if (annotentries != unitentries) {
      return fManager->HandleError(errEQTreeEntries, unittree->GetName(), anntree->GetName());
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit) output << sep << "UnitName";
   if (hasGene) output << sep << "GeneName";
   if (hasSymb) output << sep << "GeneSymbol";
   if (hasChro) output << sep << "Chromosome";
   if (hasCyto) output << sep << "CytoBand";
   if (hasSequ) output << sep << "TargetSequence";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<unitentries; i++) {
      unittree->GetEntry(i);
      anntree->GetEntry(i);
      output << unit->GetUnitID();
      if (hasUnit) output << sep << unit->GetUnitName();
      if (hasGene) output << sep << annot->GetName();
      if (hasSymb) output << sep << annot->GetSymbol();
      if (hasChro) output << sep << annot->GetChromosome();
      if (hasCyto) output << sep << annot->GetCytoBand();
      if (hasSequ) output << sep << annot->GetSequence();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << unitentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportTransAnnotTree

//______________________________________________________________________________
Int_t XDNAChip::ExportTransAnnotXML(Int_t n, TString *names, const char *varlist,
                ofstream &output, const char *sep)
{
   // Export data stored in annotation tree to file output.xml
   if(kCS) cout << "------XDNAChip::ExportTransAnnotXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportTransAnnotXML

//______________________________________________________________________________
Int_t XDNAChip::ImportTransAnnotation(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree "tree.ann":
   // UnitID, UnitName, GeneName, GeneSymbol, Chromosome, CytoBand, TargetSequence
   if(kCS) cout << "------XDNAChip::ImportTransAnnotation------" << endl;

   char nextline[kBuf4096];
   Int_t i, k, idx;
   Int_t err = errNoErr;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Get unit tree
   TString unitname = fName + "." + TString(kExtenScheme[2]);
   XUnit *unit = 0;
   TTree *unittree = (TTree*)gDirectory->Get(unitname);
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch",&unit);

// Create new annotation tree
   fAnnTreeName   = TString(fName) + "." + exten;
   TTree *anntree = new TTree(fAnnTreeName, "annotation for units");
   if (anntree == 0) return errCreateTree;
   XAnnotation *annot = 0;
   annot = new XAnnotation();
   anntree->Branch("AnnBranch", "XAnnotation", &annot, 64000, split);

// Check if header line begins with "UnitName"
   input.getline(nextline, kBuf4096, delim);
   if ( strncmp("UnitName", nextline, 8) != 0) return errMissingColumn;
//?? determine columns according to column names?

// Create hash table to store unit names from input
   XIdxString *idxstr = 0;
   THashTable *htable = 0;
   if (!(htable = new THashTable(2*fNUnits))) {err = errInitMemory;}

// Create local arrays to store data from input
   TString *arrUnit  = 0;
   TString *arrSymb  = 0;
   TString *arrName  = 0;
   TString *arrChrom = 0;
   TString *arrCyto  = 0;
   TString *arrSeq   = 0;

   // initialize memory for local arrays
   Int_t size = fNUnits;
   if (!(arrUnit  = new (nothrow) TString[size])) {err = errInitMemory;}
   if (!(arrSymb  = new (nothrow) TString[size])) {err = errInitMemory;}
   if (!(arrName  = new (nothrow) TString[size])) {err = errInitMemory;}
   if (!(arrChrom = new (nothrow) TString[size])) {err = errInitMemory;}
   if (!(arrCyto  = new (nothrow) TString[size])) {err = errInitMemory;}
   if (!(arrSeq   = new (nothrow) TString[size])) {err = errInitMemory;}

   if (!err) {
   // Initialize local arrays
      for (i=0; i<size; i++) {
         arrUnit[i]  = "NA";
         arrSymb[i]  = "NA";
         arrName[i]  = "NA";
         arrChrom[i] = "NA";
         arrCyto[i]  = "NA";
         arrSeq[i]   = "NA";
      }//for_i

   // Read data
      TString str, str1, str2, str3, str4, str5, str6;
      idx = 0;
      while (input.good()) {
         input.getline(nextline, kBuf4096, delim);
         if (input.fail() || (idx == size)) break;
         // read fields
         str1  = strtok(&nextline[0], sep);
         str2  = strtok(NULL, sep);
         str3  = strtok(NULL, sep);
         str4  = strtok(NULL, sep);
         str5  = strtok(NULL, sep);
         str6  = strtok(NULL, sep);

         // remove non-alphanumeric characters from ends
         arrUnit[idx]  = RemoveEnds(str1);
         arrSymb[idx]  = RemoveEnds(str2);
         arrName[idx]  = RemoveEnds(str3);
         arrChrom[idx] = RemoveEnds(str4);
         arrCyto[idx]  = RemoveEnds(str5);
         arrSeq[idx]   = RemoveEnds(str6);

         idxstr = new XIdxString(idx, arrUnit[idx]);
         htable->Add(idxstr);
         idx++;
      }//while

      if (idx != size) {
         cout << "Warning: Number of lines read <" << idx  
              << "> is not equal to number of units <" << size << ">" << endl;
      }//if

   // Fill annotation tree
      const char *unitname;  //ccc new??
      for (i=0; i<size; i++) {
         unittree->GetEntry(i);
         unitname = unit->GetUnitName();

         idxstr = (XIdxString*)(htable->FindObject(unitname));
         if (idxstr) {
            k = idxstr->GetIndex();
            annot->SetName(arrName[k]);
            annot->SetSymbol(arrSymb[k]);
            annot->SetChromosome(arrChrom[k]);
            annot->SetCytoBand(arrCyto[k]);
            annot->SetSequence(arrSeq[k]);
         } else {
            annot->SetName("NA");
            annot->SetSymbol("NA");
            annot->SetChromosome("NA");
            annot->SetCytoBand("NA");
            annot->SetSequence("NA");
         }//if
         anntree->Fill();

         if (XManager::fgVerbose && i%10000 == 0) {
            cout << "   <" << i+1 << "> records imported...\r" << flush;
         }//if
      }//for_i
      if (XManager::fgVerbose) {
         cout << "   <" << i << "> records imported...Finished" << endl;
      }//if
//TO DO:
   // Add tree info to tree
//      AddAnnotationTreeInfo(anntree, xxxxx);

      // Write annotation tree to file
      if ((err = WriteTree(anntree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(anntree->GetName(), 0);
      }//if
      //delete tree from memory
      anntree->Delete("");
      anntree = 0;
//no      SafeDelete(idxstr);  //deleted in htable!!
   }//if

// Clean up
   if (arrUnit)  {delete [] arrUnit;  arrUnit  = 0;}
   if (arrSymb)  {delete [] arrSymb;  arrSymb  = 0;}
   if (arrName)  {delete [] arrName;  arrName  = 0;}
   if (arrChrom) {delete [] arrChrom; arrChrom = 0;}
   if (arrCyto)  {delete [] arrCyto;  arrCyto  = 0;}
   if (arrSeq)   {delete [] arrSeq;   arrSeq   = 0;}
   if (htable)   {htable->Delete(); delete htable; htable = 0;}
   SafeDelete(annot);

   return err;
}//ImportTransAnnotation

//______________________________________________________________________________
Int_t XDNAChip::ImportTransAnnotationXML(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile in XML-format and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportTransAnnotationXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportTransAnnotation

//______________________________________________________________________________
Int_t XDNAChip::ImportExonAnnotation(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportExonAnnotation------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportExonAnnotation

//______________________________________________________________________________
Int_t XDNAChip::ImportExonAnnotationXML(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile in XML-format and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportExonAnnotationXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportExonAnnotationXML

//______________________________________________________________________________
Int_t XDNAChip::ImportProbesetAnnotation(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportProbesetAnnotation------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportProbesetAnnotation

//______________________________________________________________________________
Int_t XDNAChip::ImportProbesetAnnotationXML(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile in XML-format and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportProbesetAnnotationXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportProbesetAnnotationXML

//______________________________________________________________________________
Int_t XDNAChip::ImportControlAnnotation(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportControlAnnotation------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportControlAnnotation

//______________________________________________________________________________
Int_t XDNAChip::ImportControlAnnotationXML(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import annotation from infile in XML-format and store as annotation tree
   if(kCS) cout << "------XDNAChip::ImportControlAnnotationXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportControlAnnotationXML



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XMicroArray                                                          //
//                                                                      //
// Class for cDNA microarrays                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XMicroArray::XMicroArray()
            :XDNAChip()
{
   // Default MicroArray constructor
   if(kCS) cout << "---XMicroArray::XMicroArray(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XMicroArray::XMicroArray(const char *name, const char *title)
            :XDNAChip(name, title)
{
   // Normal MicroArray constructor
   if(kCS) cout << "---XMicroArray::XMicroArray------" << endl;

}//Constructor

//______________________________________________________________________________
XMicroArray::XMicroArray(const XMicroArray &chip) 
            :XDNAChip(chip)
{
   // MicroArray copy constructor
   if(kCS) cout << "---XMicroArray::XMicroArray(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XMicroArray& XMicroArray::operator=(const XMicroArray& rhs)
{
   // MicroArray assignment operator.

   if (this != &rhs) {
      XDNAChip::operator=(rhs);
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XMicroArray::~XMicroArray()
{
   // MicroArray destructor
   if(kCS) cout << "---XMicroArray::~XMicroArray------" << endl;

}//Destructor

//______________________________________________________________________________
void XMicroArray::PrintInfo()
{
   // Print microarray information
   if(kCS) cout << "------XMicroArray::PrintInfo------" << endl;

}//PrintInfo

//______________________________________________________________________________
Int_t XMicroArray::ExportProbeTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in probe information tree to file output
   if(kCS) cout << "------XMicroArray::ExportProbeTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportProbeTree

//______________________________________________________________________________
Int_t XMicroArray::ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output
   if(kCS) cout << "------XMicroArray::ExportSchemeTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportSchemeTree

//______________________________________________________________________________
Int_t XMicroArray::ExportUnitTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output
   if(kCS) cout << "------XMicroArray::ExportUnitTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportUnitTree

//______________________________________________________________________________
Int_t XMicroArray::ExportProbeXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in probe information tree to file output.xml
   if(kCS) cout << "------XMicroArray::ExportProbeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportProbeXML

//______________________________________________________________________________
Int_t XMicroArray::ExportSchemeXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output.xml
   if(kCS) cout << "------XMicroArray::ExportSchemeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportSchemeXML

//______________________________________________________________________________
Int_t XMicroArray::ExportUnitXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output.xml
   if(kCS) cout << "------XMicroArray::ExportUnitXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportUnitXML

//______________________________________________________________________________
Int_t XMicroArray::ImportLayout(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import layout information from infile and store as probe tree
   if(kCS) cout << "------XMicroArray::ImportLayout------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportLayout

//______________________________________________________________________________
Int_t XMicroArray::ImportLayoutXML(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import layout information from infile in XML-format and store as probe tree
   if(kCS) cout << "------XMicroArray::ImportLayoutXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportLayoutXML

//______________________________________________________________________________
Int_t XMicroArray::ImportProbeInfo(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import probe information from infile and store as probe tree
   if(kCS) cout << "------XMicroArray::ImportProbeInfo------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportProbeInfo

//______________________________________________________________________________
Int_t XMicroArray::ImportProbeInfoXML(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import probe information from infile in XML-format and store as probe tree
   if(kCS) cout << "------XMicroArray::ImportProbeInfoXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportProbeInfoXML

//______________________________________________________________________________
Int_t XMicroArray::ImportScheme(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import scheme data from infile and store in branch of tree
   if(kCS) cout << "------XMicroArray::ImportScheme------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportScheme

//______________________________________________________________________________
Int_t XMicroArray::ImportSchemeXML(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import scheme data from infile in XML-format and store in branch of tree
   if(kCS) cout << "------XMicroArray::ImportSchemeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportSchemeXML


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XOligoArray                                                          //
//                                                                      //
// Class for oligonucleotide arrays                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XOligoArray::XOligoArray()
            :XDNAChip()
{
   // Default OligoArray constructor
   if(kCS) cout << "---XOligoArray::XOligoArray(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XOligoArray::XOligoArray(const char *name, const char *title)
            :XDNAChip(name, title)
{
   // Normal OligoArray constructor
   if(kCS) cout << "---XOligoArray::XOligoArray------" << endl;

}//Constructor

//______________________________________________________________________________
XOligoArray::XOligoArray(const XOligoArray &chip) 
            :XDNAChip(chip)
{
   // OligoArray copy constructor
   if(kCS) cout << "---XOligoArray::XOligoArray(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XOligoArray& XOligoArray::operator=(const XOligoArray& rhs)
{
   // OligoArray assignment operator.

   if (this != &rhs) {
      XDNAChip::operator=(rhs);
   }//if
   return *this;
}//operator=

//______________________________________________________________________________
XOligoArray::~XOligoArray()
{
   // OligoArray destructor
   if(kCS) cout << "---XOligoArray::~XOligoArray------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XOligoArray::ExportLayoutTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in layout tree to file output
   if(kCS) cout << "------XOligoArray::ExportLayoutTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportLayoutTree

//______________________________________________________________________________
Int_t XOligoArray::ExportProbeTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in probe information tree to file output
   if(kCS) cout << "------XOligoArray::ExportProbeTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportProbeTree

//______________________________________________________________________________
Int_t XOligoArray::ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output
   if(kCS) cout << "------XOligoArray::ExportSchemeTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportSchemeTree

//______________________________________________________________________________
Int_t XOligoArray::ExportUnitTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output
   if(kCS) cout << "------XOligoArray::ExportUnitTree------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   return errNoErr;
}//ExportUnitTree

//______________________________________________________________________________
Int_t XOligoArray::ExportLayoutXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in layout tree to file output.xml
   if(kCS) cout << "------XOligoArray::ExportLayoutXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportLayoutXML

//______________________________________________________________________________
Int_t XOligoArray::ExportProbeXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in probe information tree to file output.xml
   if(kCS) cout << "------XOligoArray::ExportProbeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportProbeXML

//______________________________________________________________________________
Int_t XOligoArray::ExportSchemeXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output.xml
   if(kCS) cout << "------XOligoArray::ExportSchemeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportSchemeXML

//______________________________________________________________________________
Int_t XOligoArray::ExportUnitXML(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output.xml
   if(kCS) cout << "------XOligoArray::ExportUnitXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportUnitXML

//______________________________________________________________________________
Int_t XOligoArray::ImportLayout(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import layout information from infile and store as probe tree
   if(kCS) cout << "------XOligoArray::ImportLayout------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportLayout

//______________________________________________________________________________
Int_t XOligoArray::ImportLayoutXML(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import layout information from infile in XML-format and store as probe tree
   if(kCS) cout << "------XOligoArray::ImportLayoutXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportLayoutXML

//______________________________________________________________________________
Int_t XOligoArray::ImportProbeInfo(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import probe information from infile and store as probe tree
   if(kCS) cout << "------XOligoArray::ImportProbeInfo------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportProbeInfo

//______________________________________________________________________________
Int_t XOligoArray::ImportProbeInfoXML(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import probe information from infile in XML-format and store as probe tree
   if(kCS) cout << "------XOligoArray::ImportProbeInfoXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportProbeInfoXML

//______________________________________________________________________________
Int_t XOligoArray::ImportScheme(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import scheme data from infile and store in branch of tree
   if(kCS) cout << "------XOligoArray::ImportScheme------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   return errNoErr;
}//ImportScheme

//______________________________________________________________________________
Int_t XOligoArray::ImportSchemeXML(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import scheme data from infile in XML-format and store in branch of tree
   if(kCS) cout << "------XOligoArray::ImportSchemeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportSchemeXML


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGeneChip                                                            //
//                                                                      //
// Class for Affymetrix GeneChip arrays                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGeneChip::XGeneChip()
          :XOligoArray()
{
   // Default GeneChip constructor
   if(kCS) cout << "---XGeneChip::XGeneChip(default)------" << endl;

   fAnpTreeName = "";
   fNProbesets  =  0;
   fNProbes     =  0;
   fNAffx       =  0;
   fMinGC       = -1;
   fMaxGC       = -1;
   fMinTm       = -1.0;
   fMaxTm       = -1.0;
}//Constructor

//______________________________________________________________________________
XGeneChip::XGeneChip(const char *name, const char *title)
          :XOligoArray(name, title)
{
   // Normal GeneChip constructor
   if(kCS) cout << "---XGeneChip::XGeneChip------" << endl;

   fAnpTreeName = "NA";
   fNProbesets  =  0;
   fNProbes     =  0;
   fNAffx       =  0;
   fMinGC       = -1;
   fMaxGC       = -1;
   fMinTm       = -1.0;
   fMaxTm       = -1.0;
}//Constructor

//______________________________________________________________________________
XGeneChip::XGeneChip(const XGeneChip &chip) 
          :XOligoArray(chip)
{
   // GeneChip copy constructor
   if(kCS) cout << "---XGeneChip::XGeneChip(copy)------" << endl;

   fAnpTreeName = chip.fAnpTreeName;
   fNProbesets  = chip.fNProbesets;
   fNProbes     = chip.fNProbes;
   fNAffx       = chip.fNAffx;
   fMinGC       = chip.fMinGC;
   fMaxGC       = chip.fMaxGC;
   fMinTm       = chip.fMinTm;
   fMaxTm       = chip.fMaxTm;
}//CopyConstructor

//______________________________________________________________________________
XGeneChip& XGeneChip::operator=(const XGeneChip& rhs)
{
   // GeneChip assignment operator.

   if (this != &rhs) {
      XGeneChip::operator=(rhs);

      fAnpTreeName = rhs.fAnpTreeName;
      fNProbesets  = rhs.fNProbesets;
      fNProbes     = rhs.fNProbes;
      fNAffx       = rhs.fNAffx;
      fMinGC       = rhs.fMinGC;
      fMaxGC       = rhs.fMaxGC;
      fMinTm       = rhs.fMinTm;
      fMaxTm       = rhs.fMaxTm;
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XGeneChip::~XGeneChip()
{
   // GeneChip destructor
   if(kCS) cout << "---XGeneChip::~XGeneChip------" << endl;

}//Destructor

//______________________________________________________________________________
void XGeneChip::PrintInfo()
{
   // Print GeneChip information
   if(kCS) cout << "------XGeneChip::PrintInfo------" << endl;

   if (fgPrintHeader) {
      cout << "==============================================================================" << endl;
      cout << setw(14) << "ChipName" << setw(12) << "Title"
           << setw(17) << "SchemeTree" << setw(17) << "UnitTree" 
           << setw(17) << "ProbeTree" << setw(17) << "AnnotTree" 
           << setw(9) << "NUnits" << setw(9) << "NGenes" << endl;
      cout << "==============================================================================" << endl;
      fgPrintHeader = kFALSE;
   }//if

   cout << setw(14) << this->GetName()     << setw(12) << this->GetTitle()
        << setw(17) << fScmTreeName.Data() << setw(17) << fIdxTreeName.Data() 
        << setw(17) << fPrbTreeName.Data() << setw(17) << fAnnTreeName.Data() 
        << setw(9)  << fNUnits             << setw(9)  << fNGenes << endl;

   cout << "------------------------------------------------------------------------------" << endl;
}//PrintInfo

//______________________________________________________________________________
void XGeneChip::AddLayoutTreeInfo(TTree *tree, const char *name, Option_t *option)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGeneChip::AddLayoutTreeInfo------" << endl;

   XLayoutTreeInfo *info = new XLayoutTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   tree->GetUserInfo()->Add(info);
}//AddLayoutTreeInfo

//______________________________________________________________________________
void XGeneChip::AddProbeTreeInfo(TTree *tree, const char *name, Option_t *option)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGeneChip::AddProbeTreeInfo------" << endl;

   XProbeTreeInfo *info = new XProbeTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   tree->GetUserInfo()->Add(info);
}//AddProbeTreeInfo

//______________________________________________________________________________
void XGeneChip::AddSchemeTreeInfo(TTree *tree, const char *name, Option_t *option)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGeneChip::AddSchemeTreeInfo------" << endl;

   XSchemeTreeInfo *info = new XSchemeTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(this);

   tree->GetUserInfo()->Add(info);
}//AddSchemeTreeInfo

//______________________________________________________________________________
void XGeneChip::AddUnitTreeInfo(TTree *tree, const char *name, Option_t *option,
                Int_t nctrls, Int_t naffx, Int_t ngenes,
                Int_t nmin, Int_t min, Int_t nmax, Int_t max)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGeneChip::AddUnitTreeInfo------" << endl;

   XUnitTreeInfo *info = new XUnitTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(nctrls, naffx, ngenes, nmin, min, nmax, max);

   tree->GetUserInfo()->Add(info);
}//AddUnitTreeInfo

//______________________________________________________________________________
Int_t XGeneChip::ContentGC(const char *sequence, Double_t &Tm, const char *method)
{
   // Return number of G/C in sequence and melting temperature for method:
   // method = "none": Tm is set to Tm = -1.0
   // method = "empirical": Tm = 4*(G + C) + 2*(A + T) - 5
   if(kCSa) cout << "------XGeneChip::ContentGC------" << endl;

   // calculate number of G/C in sequence
   Int_t length = strlen(sequence);
   Int_t numGC  = 0;
   for (Int_t k=0; k<length; k++) {
      (sequence[k] == 'C') ? numGC++ : numGC;
      (sequence[k] == 'G') ? numGC++ : numGC;
   }//for_k

   // calculate melting temperature
   Tm = this->MeltingTemperature(numGC, length, method);

   return numGC;
}//ContentGC

//______________________________________________________________________________
Double_t XGeneChip::MeltingTemperature(Int_t numGC, Int_t length, const char *method)
{
   // Return oligo melting temperature dependent on G/C number for method:
   // method = "none": Tm is set to Tm = -1.0
   // method = "empirical": Tm = 4*(G + C) + 2*(A + T) - 5
   if(kCSa) cout << "------XGeneChip::MeltingTemperature------" << endl;

   // calculate melting temperature
   Double_t Tm = -1.0;
   if (strcmp(method, "none") == 0) {
      Tm = -1.0;
   } else if (strcmp(method, "empirical") == 0) {
      Tm = 4.0*numGC + 2.0*(length - numGC) - 5.0;
   }//if

   return Tm;
}//MeltingTemperature

//______________________________________________________________________________
Int_t XGeneChip::ExportProbeTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in probe information tree to file output
   if(kCS) cout << "------XGeneChip::ExportProbeTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasIPos = kFALSE;
   Bool_t hasSequ = kFALSE;
   Bool_t hasNrGC = kFALSE;
   Bool_t hasTmel = kFALSE;
   Bool_t hasOrie = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasIPos = kTRUE;
      hasSequ = kTRUE;
      hasNrGC = kTRUE;
      hasTmel = kTRUE;
      hasOrie = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fPosition")    == 0) {hasIPos = kTRUE;}
         if (strcmp(name,"fSequence")    == 0) {hasSequ = kTRUE;}
         if (strcmp(name,"fNumberGC")    == 0) {hasNrGC = kTRUE;}
         if (strcmp(name,"fTm")          == 0) {hasTmel = kTRUE;}
         if (strcmp(name,"fIsAntisense") == 0) {hasOrie = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get scheme tree
   TString schemename = fName + "." + TString(kExtenScheme[0]);
   XScheme *scheme = 0;
   TTree *schemetree = (TTree*)gDirectory->Get(schemename); 
   if (schemetree == 0) return errGetTree;
   schemetree->SetBranchAddress("ScmBranch", &scheme);

// Get probe tree
   XGCProbe  *probe = 0;
   TTree *probetree = (TTree*)gDirectory->Get(names[0].Data()); 
   if (probetree == 0) return errGetTree;
   probetree->SetBranchAddress("PrbBranch", &probe);

   Int_t scm_entries = (Int_t)(schemetree->GetEntries());
   Int_t prb_entries = (Int_t)(probetree->GetEntries());
   if (scm_entries != prb_entries) {
      cerr << "Error: Number of scheme tree entries <" << scm_entries
           << "> is not equal to number of probe tree entries <" << prb_entries
           << ">." << endl;
      return errGeneral;
   }//if

// Output header
   output << "ProbeSetID" << sep << "ProbeX" << sep << "ProbeY";
   if (hasIPos) output << sep << "ProbeInterrogationPosition";
   if (hasSequ) output << sep << "ProbeSequence";
   if (hasNrGC) output << sep << "ContentGC";
   if (hasTmel) output << sep << "Tm";
   if (hasOrie) output << sep << "ProbeType";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<prb_entries; i++) {
      schemetree->GetEntry(i);
      probetree->GetEntry(i);
      output << scheme->GetUnitID() << sep
             << probe->GetX() << sep << probe->GetY();
      if (hasIPos) output << sep << probe->GetPosition();
      if (hasSequ) output << sep << probe->GetSequence();
      if (hasNrGC) output << sep << probe->GetNumberGC();
      if (hasTmel) output << sep << probe->GetTMelting();
      if (hasOrie) output << sep << probe->GetProbeType();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << prb_entries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportProbeTree

//______________________________________________________________________________
Int_t XGeneChip::ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output
   // varlist can be at most: "fProbeLen:fMask:fAtomNo:fPBase:fTBase" or "*"
   if(kCS) cout << "------XGeneChip::ExportSchemeTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasPLen = kFALSE;
   Bool_t hasMask = kFALSE;
   Bool_t hasAtom = kFALSE;
   Bool_t hasPBas = kFALSE;
   Bool_t hasTBas = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasPLen = kTRUE;
      hasMask = kTRUE;
      hasAtom = kTRUE;
      hasPBas = kTRUE;
      hasTBas = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fProbeLen") == 0) {hasPLen = kTRUE;}
         if (strcmp(name,"fMask")     == 0) {hasMask = kTRUE;}
         if (strcmp(name,"fAtomNo")   == 0) {hasAtom = kTRUE;}
         if (strcmp(name,"fPBase")    == 0) {hasPBas = kTRUE;}
         if (strcmp(name,"fTBase")    == 0) {hasTBas = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get scheme tree for this chip
   XGCScheme *scheme = 0;
   fTree = (TTree*)gDirectory->Get((names[0]).Data());
   if (!fTree) return errGetTree;
   fTree->SetBranchAddress("ScmBranch", &scheme);

   Int_t nentries = (Int_t)(fTree->GetEntries());
   Int_t size     = fNRows*fNCols;
   if (nentries != size) {
      TString str = ""; str += size;
      return fManager->HandleError(errNumTreeEntries, fTree->GetName(), str);
   }//if

// Output header
   output << "UNIT_ID" << sep << "X" << sep << "Y";
   if (hasPLen) output << sep << "ProbeLength";
   if (hasMask) output << sep << "Mask";
   if (hasAtom) output << sep << "Atom";
   if (hasPBas) output << sep << "ProbeBase";
   if (hasTBas) output << sep << "TargetBase";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << scheme->GetUnitID() << sep
             << scheme->GetX() << sep << scheme->GetY();
      if (hasPLen) output << sep << scheme->GetProbeLength();
      if (hasMask) output << sep << scheme->GetMask();
      if (hasAtom) output << sep << scheme->GetAtomNumber();
      if (hasPBas) output << sep << scheme->GetProbeBase();
      if (hasTBas) output << sep << scheme->GetTargetBase();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportSchemeTree

//______________________________________________________________________________
Int_t XGeneChip::ExportUnitTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output
   if(kCS) cout << "------XGeneChip::ExportUnitTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasUnit = kFALSE;
   Bool_t hasCell = kFALSE;
   Bool_t hasAtom = kFALSE;
   Bool_t hasType = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasUnit = kTRUE;
      hasCell = kTRUE;
      hasAtom = kTRUE;
      hasType = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName") == 0) {hasUnit = kTRUE;}
         if (strcmp(name,"fNCells")   == 0) {hasCell = kTRUE;}
         if (strcmp(name,"fNAtoms")   == 0) {hasAtom = kTRUE;}
         if (strcmp(name,"fUnitType") == 0) {hasType = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get unit tree for this chip
   XGCUnit *unit = 0;
   fTree = (TTree*)gDirectory->Get((names[0]).Data());
   if (fTree == 0) return errGetTree;
   fTree->SetBranchAddress("IdxBranch", &unit);

   Int_t nentries = (Int_t)(fTree->GetEntries());
   if (nentries != fNUnits) {
      TString str = ""; str += fNUnits;
      return fManager->HandleError(errNumTreeEntries, fTree->GetName(), str);
   }//if

// Output header
   output << "UNIT_ID";
   if (hasUnit) output << sep << "UnitName";
   if (hasCell) output << sep << "NumCells";
   if (hasAtom) output << sep << "NumAtoms";
   if (hasType) output << sep << "UnitType";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << unit->GetUnitID();
      if (hasUnit) output << sep << unit->GetUnitName();
      if (hasCell) output << sep << unit->GetNumCells();
      if (hasAtom) output << sep << unit->GetNumAtoms();
      if (hasType) output << sep << unit->GetUnitType();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportUnitTree

//______________________________________________________________________________
Int_t XGeneChip::ExportTransAnnotTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in transcript annotation tree to file output
   if(kCS) cout << "------XGeneChip::ExportTransAnnotTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasTran = kFALSE;  //transcript_id
   Bool_t hasName = kFALSE;  //gene name 
   Bool_t hasSymb = kFALSE;  //gene symbol
   Bool_t hasAccn = kFALSE;  //mRNA accession
   Bool_t hasEntr = kFALSE;  //entrez ID
   Bool_t hasChro = kFALSE;  //chromosome
   Bool_t hasCyto = kFALSE;  //cytoband
   Bool_t hasStar = kFALSE;  //start position
   Bool_t hasStop = kFALSE;  //stop position
   Bool_t hasStrd = kFALSE;  //strand

   if (strcmp(varlist,"*")  == 0) {
      hasTran = kTRUE;
      hasName = kTRUE;
      hasSymb = kTRUE;
      hasAccn = kTRUE;
      hasEntr = kTRUE;
      hasChro = kTRUE;
      hasCyto = kTRUE;
      hasStar = kTRUE;
      hasStop = kTRUE;
      hasStrd = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fTranscriptID") == 0) {hasTran = kTRUE;}
         if (strcmp(name,"fName")         == 0) {hasName = kTRUE;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymb = kTRUE;}
         if (strcmp(name,"fAccession")    == 0) {hasAccn = kTRUE;}
         if (strcmp(name,"fEntrezID")     == 0) {hasEntr = kTRUE;}
         if (strcmp(name,"fChromosome")   == 0) {hasChro = kTRUE;}
         if (strcmp(name,"fCytoBand")     == 0) {hasCyto = kTRUE;}
         if (strcmp(name,"fStart")        == 0) {hasStar = kTRUE;}
         if (strcmp(name,"fStop")         == 0) {hasStop = kTRUE;}
         if (strcmp(name,"fStrand")       == 0) {hasStrd = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get transcript annotation tree
   XTransAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get((names[0].Data())); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

// Output header
   output << "UNIT_ID";
   if (hasTran) output << sep << "ProbesetID";
   if (hasName) output << sep << "GeneName";
   if (hasSymb) output << sep << "GeneSymbol";
   if (hasAccn) output << sep << "GeneAccession";
   if (hasEntr) output << sep << "EntrezID";
   if (hasChro) output << sep << "Chromosome";
   if (hasCyto) output << sep << "Cytoband";
   if (hasStar) output << sep << "Start";
   if (hasStop) output << sep << "Stop";
   if (hasStrd) output << sep << "Strand";
   output << endl;

// Export selected variables
   Int_t nentries = (Int_t)(anntree->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      anntree->GetEntry(i);
      output << annot->GetUnitID();
      if (hasTran) output << sep << annot->GetTranscriptID();
      if (hasName) output << sep << annot->GetName();
      if (hasSymb) output << sep << annot->GetSymbol();
      if (hasAccn) output << sep << annot->GetAccession();
      if (hasEntr) output << sep << annot->GetEntrezID();
      if (hasChro) output << sep << annot->GetChromosome();
      if (hasCyto) output << sep << annot->GetCytoBand();
      if (hasStar) output << sep << annot->GetStart();
      if (hasStop) output << sep << annot->GetStop();
      if (hasStrd) output << sep << annot->GetStrand();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportTransAnnotTree

//______________________________________________________________________________
Int_t XGeneChip::ExportProbeXML(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in probe information tree to file output.xml
   if(kCS) cout << "------XGeneChip::ExportProbeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportProbeXML

//______________________________________________________________________________
Int_t XGeneChip::ExportSchemeXML(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output.xml
   if(kCS) cout << "------XGeneChip::ExportSchemeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportSchemeXML

//______________________________________________________________________________
Int_t XGeneChip::ExportUnitXML(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output.xml
   if(kCS) cout << "------XGeneChip::ExportUnitXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   n = 0; names = 0; varlist = 0; output << endl; sep = 0;

   cout << "Error: Export of trees as XML-files is not yet implemented" << endl;
   return 1;
}//ExportUnitXML

//______________________________________________________________________________
Int_t XGeneChip::IsBinaryFile(ifstream &input)
{
   if(kCS) cout << "------XGeneChip::IsBinaryFile------" << endl;

   // read magic number from input.
	Int_t magic = 0;
	READ_INT(input, magic);

   // need to rewind input to start
   input.seekg(ios::beg);

	return (magic == kCDFMagicNumber);
}//IsBinaryFile

//______________________________________________________________________________
Int_t XGeneChip::ImportProbeInfo(ifstream &input, Option_t *option,
                 const char *sep, char delim, Int_t split)
{
   // Import probe information from infile and store as probe tree "tree.prb":
   // ProbeSetID, ProbeX, ProbeY, ProbeInterrogationPosition, ProbeSequence,
   // ContentGC, Tm, ProbeType
   // Note: ContentGC will be calculated from ProbeSequence
   if(kCS) cout << "------XGeneChip::ImportProbeInfo------" << endl;

   char nextline[kBufSize];
   Int_t    err   = errNoErr;
   Int_t    hc    = -1;
   Int_t    nprb  = kNPrbCols;
   Int_t    idx   = 0;
   Int_t    count = 0;
   Int_t    minGC = 9999999;
   Int_t    maxGC = -1;
   Double_t minTm = DBL_MAX;  //defined in float.h
   Double_t maxTm = -1.0;

   Int_t   order, x, y, position, id;
   TString name, sequence, orientation;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Get scheme tree
   TString schemename = fName + "." + TString(kExtenScheme[0]);
   XGCScheme *scheme = 0;
   TTree *schemetree = (TTree*)gDirectory->Get(schemename); 
   if (schemetree == 0) return errGetTree;
   schemetree->SetBranchAddress("ScmBranch",&scheme);

// Create new probe tree
   fPrbTreeName = TString(fName) + "." + exten;
   TTree *probetree = new TTree(fPrbTreeName, "probe info for scheme");
   if (probetree == 0) return errCreateTree;
   XGCProbe *probe = 0;
   probe = new XGCProbe();
   probetree->Branch("PrbBranch", "XGCProbe", &probe, 64000, split);

// Create local arrays to store data from input
   Int_t    *arrPos  = 0;
   TString  *arrSeq  = 0;
   Short_t  *arrPrb  = 0;
   Int_t    *arrGC   = 0;
   Double_t *arrTm   = 0;
   Int_t   *isColumn = 0;
   TString  *names   = 0;

   // initialize memory for local arrays
   Int_t size = fNRows*fNCols;
   if (!(arrPos = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
   if (!(arrSeq = new (nothrow) TString[size]))  {err = errInitMemory; goto cleanup;}
   if (!(arrPrb = new (nothrow) Short_t[size]))  {err = errInitMemory; goto cleanup;}
   if (!(arrGC  = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
   if (!(arrTm  = new (nothrow) Double_t[size])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<size; i++) {
      arrPos[i] = -1;
      arrSeq[i] = "N";
      arrPrb[i] = -1;
      arrGC[i]  = -1;
      arrTm[i]  = -1.0;
   }//for_i

   if (!(isColumn = new (nothrow) Int_t[nprb]))   {err = errInitMemory; goto cleanup;}
   if (!(names    = new (nothrow) TString[nprb])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<nprb; i++) isColumn[i] = 0;

// Check column headers from header line: which header column is used as first column
   input.getline(nextline, kBufSize, delim);
   if      (strncmp(kPrbHeader[0], nextline, strlen(kPrbHeader[0])) == 0) hc = 0;
   else if (strncmp(kPrbHeader[1], nextline, strlen(kPrbHeader[1])) == 0) hc = 1;
   else {
      cerr << "Error: Header column must begin with <" << kPrbHeader[0] << "> or <"
           << kPrbHeader[1] << ">" << endl;
      err = errHeaderLine; goto cleanup;
   }//if

// Check column headers from header line: isColumn = -1 if column is absent
   nprb = GetHeaderOrder(&nextline[0], kPrbHeader, kNPrbCols, isColumn, sep);
   if (nprb > 0) {
      cout << "Note: The following header columns are missing: " << endl;
      for (Int_t i=0; i<kNPrbCols/2; i++) {
         if (isColumn[i*2+hc] == -1) cout << "<" << kPrbHeader[i*2+hc] << ">" << endl;
         if (i == 2 && isColumn[i*2+hc] == -1) {err = errMissingColumn; goto cleanup;}  //Probe X
         if (i == 3 && isColumn[i*2+hc] == -1) {err = errMissingColumn; goto cleanup;}  //Probe Y
         if (i == 5 && isColumn[i*2+hc] == -1) {err = errMissingColumn; goto cleanup;}  //Probe Sequence
      }//for_i
   }//if 

// Import data
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.fail() || (idx == size)) break;

      // Import fields
      nprb = TokenizeString(&nextline[0], nprb, names, sep);
      // Probe Set Name
      name = names[0];
      // Serial Order
      if (isColumn[2+hc] >= 0) order = atoi(names[isColumn[2+hc]]);
      else                     order = -1;
      // Probe X, Probe Y
      x = atoi(names[isColumn[4+hc]]);
      y = atoi(names[isColumn[6+hc]]);
      // Probe Interrogation Position
      if (isColumn[8+hc] >= 0) position = atoi(names[isColumn[8+hc]]);
      else                     position = -1;
      // Probe Sequence
      sequence = names[isColumn[10+hc]];
      // Target Strandedness
      if (isColumn[12+hc] >= 0) orientation = names[isColumn[12+hc]];
      else                      orientation = "NA";

      idx         = XY2Index(x, y);
      arrPos[idx] = position;
      arrSeq[idx] = RemoveEnds(sequence);
      arrPrb[idx] = (Int_t)(strcmp((RemoveEnds(orientation)).Data(), "Antisense") == 0);

      // calculate number of G/C in sequence and melting temperature
      arrGC[idx] = this->ContentGC(sequence, arrTm[idx], "empirical");

      // minimum/maximum GC content
      minGC = (minGC <= arrGC[idx]) ? minGC : arrGC[idx];
      maxGC = (maxGC >= arrGC[idx]) ? maxGC : arrGC[idx];

      // minimum/maximum Tm
      minTm = (minTm <= arrTm[idx]) ? minTm : arrTm[idx];
      maxTm = (maxTm >= arrTm[idx]) ? maxTm : arrTm[idx];

      count++;
      if (XManager::fgVerbose && count%10000 == 0) {
         if(count%10000 == 0) cout << "   <" << count << "> records read...\r" << flush;
      }//if
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << count << "> records read...Finished" << endl;
   }//if

// Add GC content, Tm and ProbeType to MM probes (problem: will be set for controls, too)
   for (Int_t i=0; i<fNCols; i++) {
      for (Int_t k=0; k<fNRows-1; k++) { //must be fNRows-1 since k+1
         Int_t ik0 = XY2Index(i, k);
         Int_t ik1 = XY2Index(i, k+1);

         // if no data exist for MM[x,y+1] then set to data of PM[x,y]
         if((arrGC[ik1] == -1) && (arrGC[ik0] != -1)) {
            arrGC[ik1]  = arrGC[ik0];
            arrTm[ik1]  = arrTm[ik0];
            arrPrb[ik1] = arrPrb[ik0];
         }//if
      }//for_k
   }//for_i

// Fill probe tree
   for (Int_t i=0; i<size; i++) {
      schemetree->GetEntry(i);
      id  = scheme->GetUnitID();
      x   = scheme->GetX();
      y   = scheme->GetY();
      idx = XY2Index(x, y);

      probe->SetX(x);
      probe->SetY(y);
      probe->SetSequence(arrSeq[idx]);
      probe->SetPosition(arrPos[idx]);
      // set only for MMs, not for controls (see above)
      probe->SetNumberGC(id  >= 0 ? arrGC[idx]  : -1);
      probe->SetTMelting(id  >= 0 ? arrTm[idx]  : -1);
      probe->SetProbeType(id >= 0 ? arrPrb[idx] : -1);
      probetree->Fill();

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << size << "> records imported...Finished" << endl;
   }//if

   fMinGC = minGC;
   fMaxGC = maxGC;
   fMinTm = minTm;
   fMaxTm = maxTm;

   if (XManager::fgVerbose) {
      cout << "   probe info: " << endl;
      cout << "      GC content: minimum GC is <" << minGC << ">  maximum GC is <"
           << maxGC << ">" << endl;
      cout << "      Melting Tm: minimum Tm is <" << minTm << ">  maximum Tm is <"
           << maxTm << ">" << endl;
   }//if

// Add tree info to probetree
   this->AddProbeTreeInfo(probetree, fPrbTreeName);

// Write probe tree to file
   if ((err = WriteTree(probetree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(probetree->GetName(), 0);
   }//if
   //delete tree from memory
   probetree->Delete("");
   probetree = 0;

// Clean up
cleanup:
   if (names)    {delete [] names;    names    = 0;}
   if (isColumn) {delete [] isColumn; isColumn = 0;}
   if (arrTm)    {delete [] arrTm;    arrTm    = 0;}
   if (arrGC)    {delete [] arrGC;    arrGC    = 0;}
   if (arrPrb)   {delete [] arrPrb;   arrPrb   = 0;}
   if (arrSeq)   {delete [] arrSeq;   arrSeq   = 0;}
   if (arrPos)   {delete [] arrPos;   arrPos   = 0;}
   SafeDelete(probe);

   return err;
}//ImportProbeInfo

//______________________________________________________________________________
Int_t XGeneChip::ImportProbeInfoXML(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import probe information from infile in XML-format and store as probe tree
   if(kCS) cout << "------XGeneChip::ImportProbeInfoXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportProbeInfoXML

//______________________________________________________________________________
Int_t XGeneChip::ReadHeader(ifstream &input, const char * /*sep*/, char delim)
{
   // Read header from input CDF file. 
   if(kCS) cout << "------XGeneChip::ReadHeader------" << endl;

   char  nextline[kBufSize];
   Int_t err = errNoErr;

// Check for header line [CDF]
   input.getline(nextline, kBufSize, delim);
   if ( strncmp("[CDF]", nextline, 5) != 0) return errHeaderLine;

// Read chip information following line [CHIP]
   // get row number: fNRows
   while (strncmp("Rows=", nextline, 5) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[5]), "%d", &fNRows);

   // get column number: fNCols
   while (strncmp("Cols=", nextline, 5) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[5]), "%d", &fNCols);

   // get number of genes: fNGenes
   while (strncmp("NumberOfUnits=", nextline, 14) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[14]), "%d", &fNGenes);

   // get number of QC units
   while (strncmp("NumQCUnits=", nextline, 11) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[11]), "%d", &fNControls);

   // get number of units
   fNUnits     = fNGenes + fNControls;
   fNProbesets = fNUnits;  //??

   return err;
}//ReadHeader

//______________________________________________________________________________
Int_t XGeneChip::ReadData(ifstream &input, Option_t *option,
                 const char * /*sep*/, char delim, Int_t split)
{
   // Import scheme data from input CDF file and store in the following trees:
   // Scheme tree "tree.scm": UnitID, X, Y, PLen, Mask, Atom, PBase, TBase
   // Unit   tree "tree.idx": UnitID, UnitName, NumCells, NumAtoms, UnitType
   // Controls: sign for unitID and unittype of QC set is negative and mask=-1
   // Note: missing (x,y) pairs will be added to tree.scm with mask=-2
   if(kCS) cout << "------XGeneChip::ReadData------" << endl;

   char  nextline[kBufSize];
   Int_t i, k;
   Int_t err    = errNoErr;
   Int_t min    = 99999999;
   Int_t max    = 0;
   Int_t nummin = 0;
   Int_t nummax = 0;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Create new scheme tree
   fScmTreeName      = TString(fName) + "." + exten;
   TTree *schemetree = new TTree(fScmTreeName, "scheme information");
   if (schemetree == 0) return errCreateTree;
   XGCScheme *scheme = 0;
   scheme = new XGCScheme();
   schemetree->Branch("ScmBranch", "XGCScheme", &scheme, 64000, split);

// Create new unit tree
   fIdxTreeName    = TString(fName) + "." + TString(kExtenScheme[2]);
   TTree *unittree = new TTree(fIdxTreeName, "unit information");
   if (unittree == 0) return errCreateTree;
   XGCUnit *unit = 0;
   unit = new XGCUnit();
   unittree->Branch("IdxBranch", "XGCUnit", &unit, 64000, split);

// Create table to store flag for (x,y) containing oligo data
   Int_t   **arrXY = 0;
   if (!(arrXY = new (nothrow) Int_t*[fNCols])) return errInitMemory;
   for (i=0; i<fNCols; i++) {
      arrXY[i] = 0;
   }//for_i
   for (k=0; k<fNRows; k++) {
      if (!(arrXY[k] = new (nothrow) Int_t[fNCols])) {err = errInitMemory;}
   }//for_k

   if (err == errNoErr) {
   // Initialize arrXY
      char pbase = 'N';
      char tbase = 'N';
      for (i=0; i<fNCols; i++) {
         for (k=0; k<fNRows; k++) {
            arrXY[i][k] = 0;
         }//for_k
      }//for_i

   // Read scheme data for controls
      TString unitname;
      Int_t x, y, plen, atom;
      Int_t unitID, numcells, numatoms;
      Int_t unittype = 0;  //control
      Int_t count    = 0;
      for (i=0; i<fNControls; i++) {
         // set unitID (negative ID for controls)
         unitID = i - fNControls;

         // get name of control
         while (strncmp("[QC", nextline, 3) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         unitname = &nextline[1];
         unitname = RemoveEnds(unitname.Data()); // remove "\n"

         // get type ID of control
         while (strncmp("Type=", nextline, 5) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[5]), "%d", &unittype);
         // analgous to negative unitID for controls, use negative type ID
         unittype = (unittype > 0) ? (-1)*unittype : 0;

         // get number of cells
         while (strncmp("NumberCells=", nextline, 12) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[12]), "%d", &numcells);

         // get cell header
         while (strncmp("CellHeader=", nextline, 11) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while

         numatoms = 0;
         for (k=0; k<numcells; k++) {
            input.getline(nextline, kBufSize, delim);
            sscanf(nextline, "%*4c%*[0123456789]%*[=]%d %d %*s %d %d \n",
                   &x, &y, &plen, &atom);
///////////
//?? better?: .. atom, index, match -> since match=1 is PM and match=0 is MM
// -> ev arrMask=-1 for PM and arrMask=0 (or -2 ??) for MM
// but only few QC have match!!
///////////

            numatoms    = TMath::Max(numatoms, atom);
            arrXY[x][y] = 1;

            // fill scheme tree with control data
            scheme->SetUnitID(unitID);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(plen);
            scheme->SetMask(-1);
            scheme->SetAtomNumber(atom);
            scheme->SetProbeBase(pbase);
            scheme->SetTargetBase(tbase);
            schemetree->Fill();

            if (XManager::fgVerbose && ++count%10000 == 0) {
               cout << "   <" << count + 1 << "> records imported...\r" << flush;
            }//if
         }//for_k
         numatoms = numatoms + 1;

         // fill unit tree with control data
         unit->SetUnitName(unitname);
         unit->SetUnitID(unitID);
         unit->SetUnitType(unittype);
         unit->SetNumCells(numcells);
         unit->SetNumAtoms(numatoms);
         unittree->Fill();
      }//for_i

   // Read scheme data for genes
      plen = kProbeLength;
      for (i=0; i<fNGenes; i++) {
         // get unit ID
         while (strncmp("UnitNumber=", nextline, 11) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         unitID = i;

         // get direction
//?? ev use direction?

         // get unit type
         while (strncmp("UnitType=", nextline, 9) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[9]), "%d", &unittype);

         // get unit block name
         while (strncmp("Name=", nextline, 5) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         unitname = &nextline[5];
         unitname = RemoveEnds(unitname.Data()); // remove "\n"

         // get number of atoms
         while (strncmp("NumAtoms=", nextline, 9) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[9]), "%d", &numatoms);

         // get number of cells
         while (strncmp("NumCells=", nextline, 9) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[9]), "%d", &numcells);

         // get cell header
         while (strncmp("CellHeader=", nextline, 11) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while

         // check for presence of PM/MM pairs
         if (numcells != 2*numatoms) {
            cout << "Warning: probeset <" << unitname.Data() << "> has NumAtoms="
                 << numatoms << " and NumCells=" << numcells << endl;
         }//if

         // get data
         for (k=0; k<numcells; k++) {
            input.getline(nextline, kBufSize, delim);
            sscanf(nextline,
               "%*4c%*[0123456789]%*[=]%d %d %*s %*s %*s %*d %*d %*c %c %c %d \n",
               &x, &y, &pbase, &tbase, &atom);

            arrXY[x][y] = 1;

            // fill scheme tree with gene data
            scheme->SetUnitID(unitID);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(plen);
            scheme->SetMask((pbase==tbase) ? 0 : 1);  //MM=0 and PM=1
            scheme->SetAtomNumber(atom);
            scheme->SetProbeBase(pbase);
            scheme->SetTargetBase(tbase);
            schemetree->Fill();

            if (XManager::fgVerbose && ++count%10000 == 0) {
               cout << "   <" << count + 1 << "> records imported...\r" << flush;
            }//if
         }//for_k

         // minimal number of cells
         if      (numcells <  min) {min = numcells; nummin = 1;}
         else if (numcells == min) {nummin++;}

         // maximal number of cells
         if      (numcells >  max) {max = numcells; nummax = 1;}
         else if (numcells == max) {nummax++;}

         // fill unit tree with gene data
         unit->SetUnitName(unitname);
         unit->SetUnitID(unitID);
         unit->SetUnitType(unittype);
         unit->SetNumCells(numcells);
         unit->SetNumAtoms(numatoms);
         unittree->Fill();
      }//for_i

   // Number of probes on array
      fNProbes = count;

   // Increase max number of cells for odd max since maxpairs = max/2!!
      if (TMath::Odd(max)) {
         // necessary for Citrus.CDF array since it contains some PMs w/o MMs!!
         max = max + 1;  //avoid buffer overflow if max is odd
      }//if

   // Add tree info to tree
      this->AddUnitTreeInfo(unittree, unittree->GetName(), "",
                            fNControls, fNAffx, fNGenes,
                            nummin, min/2, nummax, max/2);

   // Write unit tree to file
      if ((err = WriteTree(unittree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(unittree->GetName(), 0);
      }//if
      //delete tree from memory
      unittree->Delete("");
      unittree = 0;

   // Fill index for missing (x,y) pairs
      unitID = -(fNControls + 1);
      pbase  = 'N';
      tbase  = 'N';
      plen = atom = 0;
      for (i=0; i<fNCols; i++) {
         for (k=0; k<fNRows; k++) {
            if(arrXY[i][k] == 0) {
               // fill scheme tree for missing (x,y)
               scheme->SetUnitID(unitID);
               scheme->SetX(i);
               scheme->SetY(k);
               scheme->SetProbeLength(plen);
               scheme->SetMask(-2);
               scheme->SetAtomNumber(atom);
               scheme->SetProbeBase(pbase);
               scheme->SetTargetBase(tbase);
               schemetree->Fill();

               if (XManager::fgVerbose && ++count%10000 == 0) {
                  cout << "   <" << count + 1 << "> records imported...\r" << flush;
               }//if
            }//if
         }//for_k
      }//for_i
      if (XManager::fgVerbose) {
         cout << "   <" << count << "> records imported...Finished" << endl;
      }//if

      if (XManager::fgVerbose) {
         cout << "   PM/MM statistics: " << endl;
         cout << "      " << nummin << " cells with minimum number of PM/MM pairs: " << min/2 << endl;
         cout << "      " << nummax << " cells with maximum number of PM/MM pairs: " << max/2 << endl;
      }//if

//?? or to unittree
   // Add tree info to tree
      this->AddSchemeTreeInfo(schemetree, schemetree->GetName());

   // Write scheme tree to file
      if ((err = WriteTree(schemetree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(schemetree->GetName(), 0);
      }//if
      //delete tree from memory
      schemetree->Delete("");
      schemetree = 0;
   }//if

// Clean up
   SafeDelete(unit);
   SafeDelete(scheme);

   for (i=0; i<fNRows; i++) {
      if (arrXY[i]) {delete [] arrXY[i]; arrXY[i] = 0;}
   }//for_i
   delete [] arrXY;

   return err;
}//ReadData

//______________________________________________________________________________
Int_t XGeneChip::ReadBinaryHeader(ifstream &input, const char * /*sep*/, char /*delim*/)
{
   // Read header from binary input CDF file. 
   if(kCS) cout << "------XGeneChip::ReadBinaryHeader------" << endl;

   Int_t err = errNoErr;

// Read magic number
   Int_t magic;
   READ_INT(input, magic);
   if (magic != kCDFMagicNumber) {
      TString str = ""; str += magic;
      return fManager->HandleError(errCDFVersion, str);
   }//if

// Read version
   Int_t version;
   READ_INT(input, version);
   if (version != kCDFVersionNumber) {
      TString str = ""; str += version;
      return fManager->HandleError(errCDFVersion, str);
   }//if

// Read array dimensions: fNCols, fNRows
   unsigned short value; 
   READ_USHORT(input, value); fNCols = (Int_t)value;
   READ_USHORT(input, value); fNRows = (Int_t)value;

   // get number of units: fNGenes, fNControls
   READ_INT(input, fNGenes);
   READ_INT(input, fNControls);
   fNUnits     = fNGenes + fNControls;
   fNProbesets = fNUnits;  //??

// Read reference sequence
   char *str = NULL;
   READ_STRING(input, str);
   delete[] str; str = NULL;

   return err;
}//ReadBinaryHeader

//______________________________________________________________________________
Int_t XGeneChip::ReadBinaryData(ifstream &input, Option_t *option,
                 const char * /*sep*/, char /*delim*/, Int_t split)
{
   // Import scheme data from input CDF file and store in the following trees:
   // Scheme tree "tree.scm": UnitID, X, Y, PLen, Mask, Atom, PBase, TBase
   // Unit   tree "tree.idx": UnitID, UnitName, NumCells, NumAtoms, UnitType
   // Controls: sign for unitID and unittype of QC set is negative and mask=-1
   // Note: missing (x,y) pairs will be added to tree.scm with mask=-2
   if(kCS) cout << "------XGeneChip::ReadBinaryData------" << endl;

   Int_t i, k;
   Int_t err = errNoErr;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Create new scheme tree
   fScmTreeName      = TString(fName) + "." + exten;
   TTree *schemetree = new TTree(fScmTreeName, "scheme information");
   if (schemetree == 0) return errCreateTree;
   XGCScheme *scheme = 0;
   scheme = new XGCScheme();
   schemetree->Branch("ScmBranch", "XGCScheme", &scheme, 64000, split);

// Create new unit tree
   fIdxTreeName    = TString(fName) + "." + TString(kExtenScheme[2]);
   TTree *unittree = new TTree(fIdxTreeName, "unit information");
   if (unittree == 0) return errCreateTree;
   XGCUnit *unit = 0;
   unit = new XGCUnit();
   unittree->Branch("IdxBranch", "XGCUnit", &unit, 64000, split);

// Init variables
   TString *arrUnit = 0;
   Int_t  **arrXY   = 0;

   TString unitname;
   Int_t x, y, plen, flag;
   Int_t unitID, numcells, blockatoms, blockcells;
   Int_t atom      = 0;  //index to group atoms
   Int_t numatoms  = 1;  //number of atoms
   Int_t numblocks = 1;  //number of blocks
   Int_t unittype  = 0;  //control=1
   Int_t count     = 0;
   Int_t min       = 99999999;
   Int_t max       = 0;
   Int_t nummin    = 0;
   Int_t nummax    = 0;
   Int_t size      = fNRows*fNCols;

   char pbase = 'N';
   char tbase = 'N';
   char name[kUnitNameLength1];
   int  qcindex;
   int  geneindex;

   unsigned short usvalue;
   unsigned char  ucvalue;
            int   ivalue;

// Create array to store unit name
   if (!(arrUnit = new (nothrow) TString[fNGenes]))  {err = errInitMemory; goto cleanup;}

// Create table to store flag for (x,y) containing oligo data
   if (!(arrXY = new (nothrow) Int_t*[fNRows]))      {err = errInitMemory; goto cleanup;}
   for (i=0; i<fNRows; i++) {
      arrXY[i] = 0;
   }//for_i
   for (k=0; k<fNRows; k++) {
      if (!(arrXY[k] = new (nothrow) Int_t[fNCols])) {err = errInitMemory; goto cleanup;}
   }//for_k

// Initialize arrXY
   for (i=0; i<fNCols; i++) {
      for (k=0; k<fNRows; k++) {
         arrXY[i][k] = 0;
      }//for_k
   }//for_i

// Read gene names
   for (Int_t i=0; i<fNGenes; i++) {
      READ_FIXED_STRING(input, name, kUnitNameLength);
      arrUnit[i] = TString(name);
   }//for_i

// Skip over control indices
   for (Int_t i=0; i<fNControls; i++) {
      READ_INT(input, qcindex);
   }//for_i

// Skip over the gene indices
   for (Int_t i=0; i<fNGenes; i++) {
      READ_INT(input, geneindex);
   }//for_i

// Read scheme data for controls
   for (Int_t i=0; i<fNControls; i++) {
      // set unitID (negative ID for controls)
      unitID = i - fNControls;

      // get name of control
      unitname = "QC";
      unitname += (i + 1);

      // get type ID of control
      READ_USHORT(input, usvalue);
      unittype = (Int_t)((usvalue > 0) ? (-1)*usvalue : 0);

      // get number of cells
      READ_INT(input, numcells);

      // read cells
      for (Int_t k=0; k<numcells; k++) {
         // get X coordinate
         READ_USHORT(input, usvalue);
         x = (Int_t)usvalue;

         // get Y coordinate
         READ_USHORT(input, usvalue);
         y = (Int_t)usvalue;

         // get probe length
         READ_UCHAR(input, ucvalue);
         plen = (Int_t)ucvalue;

         // get perfect match flag
         READ_UCHAR(input, ucvalue);
         flag = (Int_t)ucvalue;

         // get background probe flag
         READ_UCHAR(input, ucvalue);

         arrXY[x][y] = 1;

         // fill scheme tree with control data
         scheme->SetUnitID(unitID);
         scheme->SetX(x);
         scheme->SetY(y);
         scheme->SetProbeLength(plen);
         scheme->SetMask(-1);
         scheme->SetAtomNumber(atom);
         scheme->SetProbeBase(pbase);
         scheme->SetTargetBase(tbase);
         schemetree->Fill();

         if (XManager::fgVerbose && ++count%10000 == 0) {
            cout << "   <" << count + 1 << "> records imported...\r" << flush;
         }//if
      }//for_k

      // fill unit tree with control data
      unit->SetUnitName(unitname);
      unit->SetUnitID(unitID);
      unit->SetUnitType(unittype);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();
   }//for_i

// Read scheme data for genes
   plen = kProbeLength;
   for (Int_t i=0; i<fNGenes; i++) {
   // Unit information
      // get unit ID
      unitID = i;

      // get type ID
      READ_USHORT(input, usvalue);
      unittype = (Int_t)usvalue;

      // get direction
      READ_UCHAR(input, ucvalue);
//?? ev use direction?

      // get number of atoms
      READ_INT(input, numatoms);

      // get number of blocks
      READ_INT(input, numblocks);

      // get number of cells
      READ_INT(input, numcells);

      // get unit number (not used)
      READ_INT(input, ivalue);

      // get number of cells per atom
      READ_UCHAR(input, ucvalue);

   // Read units
      for (Int_t j=0; j<numblocks; j++) {
      // Unit_Block information
         // get number of atoms
         READ_INT(input, blockatoms);

         // get number of cells
         READ_INT(input, blockcells);

         // get number of cells per atom
         READ_UCHAR(input, ucvalue);

         // get direction
         READ_UCHAR(input, ucvalue);
//?? ev use direction?

         // get position of the first atom
         READ_INT(input, ivalue);
//?? need to be added??!! to get start position of oligo in sequence??

         // get position of the last atom (not used)
         READ_INT(input, ivalue);

         // get block name
         READ_FIXED_STRING(input, name, kUnitNameLength);

      // Read cells
         for (Int_t k=0; k<blockcells; k++) {
         // Cell information
            // get number of atoms
            READ_INT(input, atom);

            // get X coordinate
            READ_USHORT(input, usvalue);
            x = (Int_t)usvalue;

            // get Y coordinate
            READ_USHORT(input, usvalue);
            y = (Int_t)usvalue;

            // get index position
            READ_INT(input, ivalue);
//?? not used??

            // get base of probe at substitution position
            READ_CHAR(input, pbase);

            // get base of target at interrogation position
            READ_CHAR(input, tbase);

            arrXY[x][y] = 1;

            // fill scheme tree with gene data
            scheme->SetUnitID(unitID);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(plen);
            scheme->SetMask((pbase==tbase) ? 0 : 1);  //MM=0 and PM=1
            scheme->SetAtomNumber(atom);
            scheme->SetProbeBase(pbase);
            scheme->SetTargetBase(tbase);
            schemetree->Fill();

            if (XManager::fgVerbose && ++count%10000 == 0) {
               cout << "   <" << count + 1 << "> records imported...\r" << flush;
            }//if
         }//for_k
      }//for_j

      // minimal number of cells
      if      (numcells <  min) {min = numcells; nummin = 1;}
      else if (numcells == min) {nummin++;}

      // maximal number of cells
      if      (numcells >  max) {max = numcells; nummax = 1;}
      else if (numcells == max) {nummax++;}

   // Increase max number of cells for odd max since maxpairs = max/2!!
      if (TMath::Odd(max)) {
         // necessary for plate array CDF since it contains PMs w/o MMs!!
         max = max + 1;  //avoid buffer overflow if max is odd
      }//if

      // fill unit tree with control data
      unit->SetUnitName(arrUnit[i]);
      unit->SetUnitID(unitID);
      unit->SetUnitType(unittype);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << count << "> records imported...Finished" << endl;
   }//if

// Number of probesets on chip
   fNProbes = count;

// Add tree info to tree
   this->AddUnitTreeInfo(unittree, unittree->GetName(), "",
                         fNControls, fNAffx, fNGenes,
                         nummin, min/2, nummax, max/2);

// Write unit tree to file
   if ((err = WriteTree(unittree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(unittree->GetName(), 0);
   }//if
   //delete tree from memory
   unittree->Delete("");
   unittree = 0;

// Fill scheme for missing (x,y) pairs
   unitID = -(fNControls + 1);
   pbase  = 'N';
   tbase  = 'N';
   plen = atom = 0;
   for (i=0; i<fNCols; i++) {
      for (k=0; k<fNRows; k++) {
         if(arrXY[i][k] == 0) {
            // fill scheme tree for missing (x,y)
            scheme->SetUnitID(unitID);
            scheme->SetX(i);
            scheme->SetY(k);
            scheme->SetProbeLength(plen);
            scheme->SetMask(-2);
            scheme->SetAtomNumber(atom);
            scheme->SetProbeBase(pbase);
            scheme->SetTargetBase(tbase);
            schemetree->Fill();

            if (XManager::fgVerbose && ++count%10000 == 0) {
               cout << "   <" << count + 1 << "> records imported...\r" << flush;
            }//if
         }//if
      }//for_k
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << count << "> records imported...Finished" << endl;
   }//if

   if (XManager::fgVerbose && count != size) {
      cout << "Warning: Number of scheme tree entries <" << count 
           << "> is not <" << size << ">" << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "   PM/MM statistics: " << endl;
      cout << "      " << nummin << " cells with minimum number of PM/MM pairs: " << min/2 << endl;
      cout << "      " << nummax << " cells with maximum number of PM/MM pairs: " << max/2 << endl;
   }//if

// Add tree info to tree
   this->AddSchemeTreeInfo(schemetree, schemetree->GetName());

// Write scheme tree to file
   if ((err = WriteTree(schemetree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(schemetree->GetName(), 0);
   }//if
   //delete tree from memory
   schemetree->Delete("");
   schemetree = 0;

// Clean up
cleanup:
   if (arrUnit) {delete [] arrUnit; arrUnit = 0;}
   for (i=0; i<fNRows; i++) {
      if (arrXY[i]) {delete [] arrXY[i]; arrXY[i] = 0;}
   }//for_i
   delete [] arrXY;

   SafeDelete(unit);
   SafeDelete(scheme);

   return err;
}//ReadBinaryData

//______________________________________________________________________________
Int_t XGeneChip::ImportScheme(ifstream &input, Option_t *option,
                 const char *sep, char delim, Int_t split)
{
   // Import scheme data from binary input.CDF file and store in branch of tree
   if(kCS) cout << "------XGeneChip::ImportScheme------" << endl;

   Int_t err = errNoErr;

   if (this->IsBinaryFile(input)) {
      if (err == errNoErr) err = this->ReadBinaryHeader(input, sep, delim);
      if (err == errNoErr) err = this->ReadBinaryData(input, option, sep, delim, split);
   } else {
      if (err == errNoErr) err = this->ReadHeader(input, sep, delim);
      if (err == errNoErr) err = this->ReadData(input, option, sep, delim, split);
   }//if

   return err;
}//ImportScheme

//______________________________________________________________________________
Int_t XGeneChip::ImportTransAnnotation(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import transcript annotation file
   if(kCS) cout << "------XGeneChip::ImportTransAnnotation------" << endl;

   Int_t err  = errNoErr;
   Int_t size = 0;
   Int_t idx  = 0;
   Int_t nsub = kNumSub;
   Int_t nalg = kNumAlg;
   Int_t namb = 0;  //number of ambigous annotations
   Int_t unitID;
   Int_t k;

   std::string nextline;
   streampos position;

   TString str, dummy;
   TString name, symbol, cytoband, accession, alignment, entrez;
   TString start, stop, strand;
   TString bestalign, bestcyto;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Get unit tree
   TString unitname = fName + "." + TString(kExtenScheme[2]);
   XUnit *unit = 0;
   TTree *unittree = (TTree*)gDirectory->Get(unitname);
   if (unittree == 0) return errGetTree;
   unittree->SetBranchAddress("IdxBranch",&unit);

// Create new annotation tree
   fAnnTreeName      = TString(fName) + "." + exten;
   TTree *anntree = new TTree(fAnnTreeName, "transcript annotation");
   if (anntree == 0) return errCreateTree;
   XTransAnnotation *ann = 0;
   ann = new XTransAnnotation();
   anntree->Branch("AnnBranch", "XTransAnnotation", &ann, 64000, split);

// Check if header line begins with "Probe Set ID"
   while (nextline.compare(0, 14, "\"Probe Set ID\"") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errMissingColumn;
   }//while
   Int_t numsep = NumSeparators(nextline.c_str(), sep);

// Init variables to tokenize lines
   const char *csv = "\",\""; //necessary to tokenize lines, better?: sep = "\",\"";
   const char *tab = "\t";    //replace csv with tab

// Init hash table
   XIdxString *idxstr = 0;
   THashTable *htable = 0;

// Init local arrays to store data from input columns
   TString *names    = 0;  //array containing tokenized string
   TString *align    = 0;  //array containing alignment data
   TString *arrTrans = 0;  //Probe Set ID
   TString *arrName  = 0;  //Gene Title
   TString *arrSymb  = 0;  //Gene Symbol
   TString *arrAcces = 0;  //RefSeq Transcript ID
   Int_t   *arrEntrz = 0;  //Entrez Gene
   TString *arrCyto  = 0;  //Chromosomal Location

// Init local arrays for derived data
   TString *arrChrom = 0;  //chromosome, extract from Alignments
   Int_t   *arrStart = 0;  //start position, extract from Alignments
   Int_t   *arrStop  = 0;  //stop position, extract from Alignments
   char    *arrStran = 0;  //strand, extract from Alignments

   // remove all "\"" from header line
   str = TString(nextline);
   str.ReplaceAll("\"", "");

// Create array hasColumn 
   Int_t *hasColumn = 0;
   if (!(hasColumn = new (nothrow) Int_t[kNAnnotCols])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<kNAnnotCols; i++) {
      hasColumn[i] = 0;
   }//for_i

//TO DO: see ImportProbesetAnnotation()
//better   err = CheckHeaderOrder(str, kTranscriptHeader, kNTranscriptCols, hasColumn, sep);
// Check column headers from header line, hasColumn = 1 if column is present
   err = CheckHeader(str, kAnnotHeader, kNAnnotCols, hasColumn, sep);
   if (err > 0) {
      cout << "Note: The following header columns are missing: " << endl;
      for (Int_t i=0; i<kNAnnotCols; i++) {
         if (hasColumn[i] == 0) cout << "<" << kAnnotHeader[i] << ">" << endl;
         if (i ==  0 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //Probe Set ID
         if (i == 12 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //Alignments
         if (i == 13 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //Gene Title
         if (i == 14 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //Gene Symbol
         if (i == 15 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //Chromosomal Location
         if (i == 18 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //Entrez Gene
         if (i == 23 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //RefSeq Transcript ID
      }//for_i
   }//if 

// Get number of transcripts
   // get current streamposition for rewinding input
   position = input.tellg();
   while (1) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      size++;
   }//while
   // reset input file to position
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

   if (size > fNGenes) {
      cout << "Warning: Number of lines read <" << size  
           << "> is greater than number of genes <" << fNGenes << ">" << endl;
   } else if (size < fNGenes) {
      cout << "Warning: Number of lines read <" << size  
           << "> is less than number of genes <" << fNGenes << ">" << endl;
      size = fNGenes; //to avoid initializing not enough memory
   }//if

// Create hash table to store unit names from input
   if (!(htable = new THashTable(2*fNUnits))) {err = errInitMemory; goto cleanup;}

// Initialize memory for local arrays
   if (!(names    = new (nothrow) TString[nsub])) {err = errInitMemory; goto cleanup;}
   if (!(align    = new (nothrow) TString[nalg])) {err = errInitMemory; goto cleanup;}
   if (!(arrTrans = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrName  = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrSymb  = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrAcces = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrEntrz = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(arrCyto  = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}

// Initialize memory for sorting arrays
   if (!(arrChrom = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(arrStart = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(arrStop  = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(arrStran = new (nothrow) char[size+1]))  {err = errInitMemory; goto cleanup;}

// Initialize local arrays
   for (Int_t i=0; i<size; i++) {
      arrTrans[i] = "NA";
      arrName[i]  = "NA";
      arrSymb[i]  = "NA";
      arrAcces[i] = "NA";
      arrEntrz[i] = -1;
      arrCyto[i]  = "NA";
      arrChrom[i] = "NA";
      arrStart[i] = -1;
      arrStop[i]  = -1;
      arrStran[i] = '?';
   }//for_i

// Read data
   while (input.good()) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      if (input.fail()) {
         cout << "Error: Failed reading line:" << endl;
         cout << nextline << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // replace all "\",\"" with tab "\t" and remove "\""
      str = TString(nextline);
      // first need to replace "" with NA
      str.ReplaceAll("\"\"", "\"NA\"");
      // replace tab with space to eliminate wrong tabs in Affymetrix transcript annotation files
      str.ReplaceAll(tab, " ");
      // replace csv with tab
      str.ReplaceAll(csv, tab);
      // remove all "\"" from line
      str.ReplaceAll("\"", "");

      // check number of separators
      if (numsep != NumSeparators(str, tab)) {
         cout << "Error: Wrong number of separators in line:" << endl;
         cout << str.Data() << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // import fields
      arrTrans[idx] = RemoveEnds(strtok((char*)str.Data(), tab));
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      alignment     = strtok(0, tab);
      name          = strtok(0, tab);
      symbol        = strtok(0, tab);
      cytoband      = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      entrez        = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      dummy         = strtok(0, tab);
      accession     = strtok(0, tab);

      // get chromosome, start, stop, strand 
      Double_t maxscore = 0;
      if (strcmp(alignment.Data(),"---") == 0) {
         arrChrom[idx] = "NA";
         arrStart[idx] = -1;
         arrStop[idx]  = -1;
         arrStran[idx] = '?';
         bestcyto      = "NA";
      } else {
         Int_t index = 0;
         if (alignment.Index(kSepSl3, kNumSl3, index, TString::kExact) > 0) {
            nsub = kNumSub;
            index = TokenizeString(alignment.Data(), nsub, names, kNumSl3, kSepSl3);
         } else {
            nsub = 1;
            names[0] = alignment;
         }//if

         // find Alignments entries with best score
         for(Int_t i=0; i<nsub; i++){
            nalg  = kNumAlg;
            index = TokenizeString(names[i].Data(), nalg, align, kNumSl2, kSepSl2);
            if (nalg < 2) {
               cerr << "Error: wrong token number for <" << names[i].Data() << ">" << endl;
               err = errReadingInput; goto cleanup;
            }//if

            if (align[1].Atof() > maxscore) {;
               bestalign = align[0];
               maxscore  = align[1].Atof();
               bestcyto  = (nalg > 2) ? align[2] : "NA";
            }//if
         }//for_i

         arrChrom[idx] = RemoveEnds(strtok((char*)bestalign.Data(), ":"));
         start         = strtok(0, "-");
         stop          = strtok(0, "(");
         strand        = strtok(0, ")");
         arrStart[idx] = start.Atoi();
         arrStop[idx]  = stop.Atoi();
         arrStran[idx] = (strand.Data())[0];

         if (strcmp(bestcyto.Data(),"NA") != 0) {
            bestcyto = RemoveEnds((arrChrom[idx] + bestcyto));
         }//if
      }//if

      // ceck if name exists
      if (strcmp(name.Data(),"---") != 0) {
         Int_t index = 0;
         index = name.Index(kSepSl3, kNumSl3, index, TString::kExact);
         // index>0 only if name is a multipart entry
         name = (index > 0) ? name(0, index) : name;
         arrName[idx] = name;
      } else {
         arrName[idx] = "NA";
      }//if

      // ceck if symbol exists
      if (strcmp(symbol.Data(),"---") != 0) {
         symbol = RemoveEnds(strtok((char*)symbol.Data(), "/"));
         arrSymb[idx] = symbol;
      } else {
         arrSymb[idx] = "NA";
      }//if

      // ceck if cytoband exists
      if (strcmp(cytoband.Data(),"---") != 0) {
         Int_t index = 0;
         if (cytoband.Index(kSepSl3, kNumSl3, index, TString::kExact) > 0) {
            nsub = kNumSub;
            index = TokenizeString(cytoband.Data(), nsub, names, kNumSl3, kSepSl3);
         } else {
            nsub = 1;
            names[0] = cytoband;
         }//if

         // find cytoband equal to chromosome with best score
         for(Int_t i=0; i<nsub; i++){
            if (names[i].Contains(arrChrom[idx])) {
               cytoband = bestcyto;
            }//if
         }//for_i

         cytoband     = RemoveEnds(strtok((char*)cytoband.Data(), "/"));
         arrCyto[idx] = cytoband;

         if ((strcmp(bestcyto.Data(),"NA") != 0) &&
            !(cytoband.Contains(arrChrom[idx]))) {
//??            !(bestcyto.Contains(cytoband) || cytoband.Contains(bestcyto))) {
            if (XManager::fgVerbose == 11) {
               cout << arrTrans[idx].Data() << ": Column 'Alignments' <" << bestcyto.Data()
                    << "> with best score <" << maxscore
                    << "> is not equal to column 'Chromosomal Location' <" << cytoband.Data()
                    << ">." << endl;
            }//if
            namb++;
//??         arrCyto[idx] = cytoband + "/(" + bestcyto + " ?)";
         }//if
      } else {
         arrCyto[idx] = "NA";
      }//if

      // get Entrez ID as integer
      if (strcmp(entrez.Data(),"---") != 0) {
         arrEntrz[idx] = atoi(entrez.Data());
      } else {
         arrEntrz[idx] = -1;
      }//if

      // ceck if accession exists
      if (strcmp(accession.Data(),"---") != 0) {
         accession = RemoveEnds(strtok((char*)accession.Data(), "/"));
         arrAcces[idx] = accession;
      } else {
         arrAcces[idx] = "NA";
      }//if

      idxstr = new XIdxString(idx, arrTrans[idx]);
      htable->Add(idxstr);
      idx++;
   }//while

   if (XManager::fgVerbose) {
      if (idx == fNGenes) {
         cout << "   Number of annotated transcripts is <" << idx << ">." << endl;
      } else {
         cout << "   Note: Number of annotated transcripts <" << idx  
              << "> is not equal to number of genes <" << fNGenes << ">" << endl;
      }//if
   }//if

   if (XManager::fgVerbose) {
      if (namb >0) {
         cout << "   Note: Number of transcripts with ambigous annotation is <"
              << namb  << ">" << endl;
      }//if
   }//if

// Fill annotation tree (
   for (Int_t i=0; i<fNUnits; i++) {
      unittree->GetEntry(i);
      unitID = unit->GetUnitID();

      //loop over control units
      if (unitID < 0) continue;

      // fill tree with annotation data
      ann->SetUnitID(unitID);

      // fill in order of unit tree entries
      unitname = unit->GetUnitName();
      idxstr   = (XIdxString*)(htable->FindObject(unitname.Data()));
      if (idxstr) {
         k = idxstr->GetIndex();
         ann->SetTranscriptID(arrTrans[k]);
         ann->SetName(arrName[k]);
         ann->SetSymbol(arrSymb[k]);
         ann->SetAccession(arrAcces[k]);
         ann->SetEntrezID(arrEntrz[k]);
         ann->SetChromosome(arrChrom[k]);
         ann->SetCytoBand(arrCyto[k]);
         ann->SetStart(arrStart[k]);
         ann->SetStop(arrStop[k]);
         ann->SetStrand(arrStran[k]);
      } else {
         ann->SetTranscriptID("NA");
         ann->SetName("NA");
         ann->SetSymbol("NA");
         ann->SetAccession("NA");
         ann->SetEntrezID(-1);
         ann->SetChromosome("NA");
         ann->SetCytoBand("NA");
         ann->SetStart(-1);
         ann->SetStop(-1);
         ann->SetStrand('?');
      }//if
      anntree->Fill();

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << size << "> records imported...Finished" << endl;
   }//if

//TO DO:
// Add tree info to tree
//   AddAnnotationTreeInfo(anntree, fAnnTreeName);

// Write annotation tree to file
   if ((err = WriteTree(anntree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(anntree->GetName(), 0);
   }//if
   //delete tree from memory
   anntree->Delete("");
   anntree = 0;

// Clean up
cleanup:
   SafeDelete(ann);

   if (arrStran) {delete [] arrStran; arrStran  = 0;}
   if (arrStop)  {delete [] arrStop;  arrStop   = 0;}
   if (arrStart) {delete [] arrStart; arrStart  = 0;}
   if (arrChrom) {delete [] arrChrom; arrChrom  = 0;}
   if (arrCyto)  {delete [] arrCyto;  arrCyto   = 0;}
   if (arrEntrz) {delete [] arrEntrz; arrEntrz  = 0;}
   if (arrAcces) {delete [] arrAcces; arrAcces  = 0;}
   if (arrSymb)  {delete [] arrSymb;  arrSymb   = 0;}
   if (arrName)  {delete [] arrName;  arrName   = 0;}
   if (arrTrans) {delete [] arrTrans; arrTrans  = 0;}
   if (align)    {delete [] align;    align     = 0;}
   if (names)    {delete [] names;    names     = 0;}
   if (htable)   {htable->Delete(); delete htable; htable = 0;}

   return err;
}//ImportTransAnnotation

//______________________________________________________________________________
Int_t XGeneChip::ImportSchemeXML(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import scheme data from infile in XML-format and store in branch of tree
   if(kCS) cout << "------XGeneChip::ImportSchemeXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of data as XML-file is not yet implemented" << endl;
   return 1;
}//ImportSchemeXML

//______________________________________________________________________________
Short_t XGeneChip::ProbeType(const char *type)
{
   // Convert probe type to id EProbeType
   if(kCSa) cout << "------XGeneChip::ProbeType------" << endl;

   if (strcmp(type, "pm:st") == 0) {
      return ePMST;
   } else if (strcmp(type, "pm:at") == 0) {
      return ePMAT;
   } else if (strcmp(type, "mm:st") == 0) {
      return eMMST;
   } else if (strcmp(type, "mm:at") == 0) {
      return eMMAT;
   }//if

   return eMMAT;
}//ProbeType


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XSNPChip                                                             //
//                                                                      //
// Class for Affymetrix Mapping arrays                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XSNPChip::XSNPChip()
         :XGeneChip()
{
   // Default SNPChip constructor
   if(kCS) cout << "---XSNPChip::XSNPChip(default)------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChip::XSNPChip(const char *name, const char *title)
         :XGeneChip(name, title)
{
   // Normal SNPChip constructor
   if(kCS) cout << "---XSNPChip::XSNPChip------" << endl;

}//Constructor

//______________________________________________________________________________
XSNPChip::XSNPChip(const XSNPChip &chip) 
         :XGeneChip(chip)
{
   // SNPChip copy constructor
   if(kCS) cout << "---XSNPChip::XSNPChip(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XSNPChip& XSNPChip::operator=(const XSNPChip& rhs)
{
   // SNPChip assignment operator.

   if (this != &rhs) {
      XGeneChip::operator=(rhs);
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XSNPChip::~XSNPChip()
{
   // XNPChip destructor
   if(kCS) cout << "---XSNPChip::~XSNPChip------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XSNPChip::ImportProbeInfo(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import probe information from infile and store as probe tree
   if(kCS) cout << "------XSNPChip::ImportProbeInfo------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of probe info for SNPChip is not yet implemented" << endl;
   return 1;
}//ImportProbeInfo

//______________________________________________________________________________
Int_t XSNPChip::ImportProbeInfoXML(ifstream &input, Option_t *option, 
                const char *sep, char delim, Int_t split)
{
   // Import probe information from infile in XML-format and store as probe tree
   if(kCS) cout << "------XSNPChip::ImportProbeInfoXML------" << endl;

// NOT YET USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Error: Import of probe info as XML-file is not yet implemented" << endl;
   return 1;
}//ImportProbeInfoXML

//______________________________________________________________________________
Int_t XSNPChip::ReadData(ifstream &input, Option_t *option,
                const char * /*sep*/, char delim, Int_t split)
{
   // Import scheme data from input CDF file and store in the following trees:
   // Scheme tree "tree.scm": UnitID, X, Y, PLen, Mask, Atom, PBase, TBase
   // Unit   tree "tree.idx": UnitID, UnitName, NumCells, NumAtoms, UnitType
   // Controls: sign for unitID and unittype of QC set is negative and mask=-1
   // Note: missing (x,y) pairs will be added to tree.scm with mask=-2
   if(kCS) cout << "------XSNPChip::ReadData------" << endl;

   char  nextline[kBufSize];
   Int_t i, k, idx;
   Int_t err = errNoErr;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Create new scheme tree
   fScmTreeName      = TString(fName) + "." + exten;
   TTree *schemetree = new TTree(fScmTreeName, "scheme information");
   if (schemetree == 0) return errCreateTree;
   XGCScheme *scheme = 0;
   scheme = new XGCScheme();
   schemetree->Branch("ScmBranch", "XGCScheme", &scheme, 64000, split);

// Create new unit tree
   fIdxTreeName    = TString(fName) + "." + TString(kExtenScheme[2]);
   TTree *unittree = new TTree(fIdxTreeName, "unit information");
   if (unittree == 0) return errCreateTree;
   XGCUnit *unit = 0;
   unit = new XGCUnit();
   unittree->Branch("IdxBranch", "XGCUnit", &unit, 64000, split);

// Create local arrays to store data from input
   Int_t   *index     = 0;
   Int_t   *arrUnitID = 0;
   Int_t   *arrX      = 0;
   Int_t   *arrY      = 0;
   Int_t   *arrAtom   = 0;
   Int_t   *arrPLen   = 0;
   char    *arrPBase  = 0;
   char    *arrTBase  = 0;
   Short_t *arrMask   = 0;
   Int_t   **arrXY    = 0;
   Int_t min    = 99999999;
   Int_t max    = 0;
   Int_t nummin = 0;
   Int_t nummax = 0;

   // initialize memory for 2-dim array
   if (!(arrXY = new (nothrow) Int_t*[fNRows])) return errInitMemory;
   for (i=0; i<fNRows; i++) {
      arrXY[i] = 0;
   }//for_i
   for (k=0; k<fNRows; k++) {
      if (!(arrXY[k] = new (nothrow) Int_t[fNCols])) {err = errInitMemory;}
   }//for_k

   // initialize memory for local arrays
   Int_t size = fNRows*fNCols;
   if (!(index     = new (nothrow) Int_t[size]))  {err = errInitMemory;}
   if (!(arrUnitID = new (nothrow) Int_t[size]))  {err = errInitMemory;}
   if (!(arrX      = new (nothrow) Int_t[size]))  {err = errInitMemory;}
   if (!(arrY      = new (nothrow) Int_t[size]))  {err = errInitMemory;}
   if (!(arrAtom   = new (nothrow) Int_t[size]))  {err = errInitMemory;}
   if (!(arrPLen   = new (nothrow) Int_t[size]))  {err = errInitMemory;}
   if (!(arrPBase  = new (nothrow) char[size+1])) {err = errInitMemory;} //size+1 or size?
   if (!(arrTBase  = new (nothrow) char[size+1])) {err = errInitMemory;}
   if (!(arrMask   = new (nothrow) Short_t[size])){err = errInitMemory;}

   if (!err) {
   // Initialize local arrays
      char pbase = 'N';
      char tbase = 'N';
      for (i=0; i<fNCols; i++) {
         for (k=0; k<fNRows; k++) {
            arrXY[i][k]    = 0;
            idx            = i*fNRows + k;
            index[idx]     = -1;
            arrUnitID[idx] = -(fNControls + 1);
            arrX[idx]      = k;
            arrY[idx]      = i;
            arrAtom[idx]   = 0;
            arrPLen[idx]   = 0;
            arrPBase[idx]  = pbase;
            arrTBase[idx]  = tbase;
            arrMask[idx]   = -2;   //-2 for missing  (x,y) pairs
         }//for_k
      }//for_i

   // Read scheme data for controls
      TString unitname;
      Int_t x, y, plen, atom;
      Int_t unitID, numcells, numatoms, numblocks, blockatoms, blockcells;
      Int_t unittype = 0;  //control
      Int_t count    = 0;
      for (i=0; i<fNControls; i++) {
         // set unitID (negative ID for controls)
         unitID = i - fNControls;

         // get name of control
         while (strncmp("[QC", nextline, 3) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         unitname = &nextline[1];
         unitname = RemoveEnds(unitname.Data()); // remove "\n"

         // get type ID of control
         while (strncmp("Type=", nextline, 5) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[5]), "%d", &unittype);
         // analgous to negative unitID for controls, use negative type ID
         unittype = (unittype > 0) ? (-1)*unittype : 0;

         // get number of cells
         while (strncmp("NumberCells=", nextline, 12) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[12]), "%d", &numcells);

         // get cell header
         while (strncmp("CellHeader=", nextline, 11) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while

         numatoms = 0;
         for (k=0; k<numcells; k++) {
            input.getline(nextline, kBufSize, delim);
            sscanf(nextline, "%*4c%*[0123456789]%*[=]%d %d %*s %d %d \n",
                   &x, &y, &plen, &atom);
///////////
//?? better: .. atom, index, match -> since match=1 is PM and match=0 is MM
// -> ev arrMask=-1 for PM and arrMask=0 (or -2 ??) for MM
// but only few QC have match!!
///////////

            numatoms       = TMath::Max(numatoms, atom);
            arrXY[x][y]    = 1;
            idx            = XY2Index(x, y);
            arrUnitID[idx] = unitID;
            arrAtom[idx]   = atom;
            arrPLen[idx]   = plen;
            arrMask[idx]   = -1;   //-1 for controls
            index[count++] = idx;
         }//for_k
         numatoms = numatoms + 1;

         // minimal number of cells
         if      (numcells <  min) {min = numcells; nummin = 1;}
         else if (numcells == min) {nummin++;}

         // maximal number of cells
         if      (numcells >  max) {max = numcells; nummax = 1;}
         else if (numcells == max) {nummax++;}

         // fill unit tree with control data
         unit->SetUnitName(unitname);
         unit->SetUnitID(unitID);
         unit->SetUnitType(unittype);
         unit->SetNumCells(numcells);
         unit->SetNumAtoms(numatoms);
         unittree->Fill();
      }//for_i

      if (XManager::fgVerbose) {
         cout << "   PM/MM statistics: " << endl;
         cout << "      " << nummin << " cells with minimum number of PM/MM pairs: " << min/2 << endl;
         cout << "      " << nummax << " cells with maximum number of PM/MM pairs: " << max/2 << endl;
      }//if

   // Read scheme data for genes
      plen = kProbeLength;
      for (i=0; i<fNGenes; i++) {
         // get unit name
         while (strncmp("Name=", nextline, 5) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         unitname = &nextline[5];
         unitname = RemoveEnds(unitname.Data()); // remove "\n"

         // get direction
//?? use direction?

         // get number of atoms
         while (strncmp("NumAtoms=", nextline, 9) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[9]), "%d", &numatoms);

         // get number of cells
         while (strncmp("NumCells=", nextline, 9) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[9]), "%d", &numcells);

         // get unit ID
         while (strncmp("UnitNumber=", nextline, 11) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         unitID = i;

         // get unit type
         while (strncmp("UnitType=", nextline, 9) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[9]), "%d", &unittype);

         // get number of blocks
         while (strncmp("NumberBlocks=", nextline, 13) != 0) {
            input.getline(nextline, kBufSize, delim);
         }//while
         sscanf(&(nextline[13]), "%d", &numblocks);

      // Read units
         for (Int_t j=0; j<numblocks; j++) {
         // Unit_Block information
            // get unit name for expression units where unitname is NONE
            if ((j == 0) && (strcmp(unitname.Data(),"NONE") == 0)) {
               while (strncmp("Name=", nextline, 5) != 0) {
                  input.getline(nextline, kBufSize, delim);
               }//while
               unitname = &nextline[5];
               unitname = unitname.Remove(unitname.Length()-1); // remove "\n"
            }//if

            // get block number

            // get number of atoms
            while (strncmp("NumAtoms=", nextline, 9) != 0) {
               input.getline(nextline, kBufSize, delim);
            }//while
            sscanf(&(nextline[9]), "%d", &blockatoms);

            // get number of cells
            while (strncmp("NumCells=", nextline, 9) != 0) {
               input.getline(nextline, kBufSize, delim);
            }//while
            sscanf(&(nextline[9]), "%d", &blockcells);

            // get cell header
            while (strncmp("CellHeader=", nextline, 11) != 0) {
               input.getline(nextline, kBufSize, delim);
            }//while

            // get start position
            // get stop position
            // get direction
//?? use direction?

         // get data for current block
            for (k=0; k<blockcells; k++) {
               input.getline(nextline, kBufSize, delim);
               sscanf(nextline,
                  "%*4c%*[0123456789]%*[=]%d %d %*s %*s %*s %*d %*d %*c %c %c %d \n",
                  &x, &y, &pbase, &tbase, &atom);

               arrXY[x][y]    = 1;
               idx            = XY2Index(x, y);
               arrUnitID[idx] = unitID;
               arrAtom[idx]   = atom;
               arrPLen[idx]   = plen;
               arrPBase[idx]  = pbase;
               arrTBase[idx]  = tbase;
               arrMask[idx]   = (pbase==tbase) ? 0 : 1;  //MM=0 and PM=1
               index[count++] = idx;
            }//for_k
         }//for_j

      // fill unit tree with gene data
         unit->SetUnitName(unitname);
         unit->SetUnitID(unitID);
         unit->SetUnitType(unittype);
         unit->SetNumCells(numcells);
         unit->SetNumAtoms(numatoms);
         unittree->Fill();
      }//for_i

   // Add tree info to tree
      this->AddUnitTreeInfo(unittree, unittree->GetName(), "", fNControls,
                            fNAffx, fNGenes, nummin, min/2, nummax, max/2);

      // Write unit tree to file
      if ((err = WriteTree(unittree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(unittree->GetName(), 0);
      }//if
      //delete tree from memory
      unittree->Delete("");
      unittree = 0;

   // Fill index for missing (x,y) pairs
      for (i=0; i<fNCols; i++) {
         for (k=0; k<fNRows; k++) {
            if(arrXY[i][k] == 0) index[count++] = XY2Index(i, k);
         }//for_k
      }//for_i

      if (count != size) {
         cerr << "Error: Number of scheme tree entries is not <"
              << size << ">" << endl;
      }//if

   // Fill scheme tree
      for (i=0; i<size; i++) {
         idx = index[i];
         scheme->SetUnitID(arrUnitID[idx]);
         scheme->SetX(arrX[idx]);
         scheme->SetY(arrY[idx]);
         scheme->SetProbeLength(arrPLen[idx]);
         scheme->SetMask(arrMask[idx]);
         scheme->SetAtomNumber(arrAtom[idx]);
         scheme->SetProbeBase(arrPBase[idx]);
         scheme->SetTargetBase(arrTBase[idx]);
         schemetree->Fill();

         if (XManager::fgVerbose && i%10000 == 0) {
            cout << "   <" << i+1 << "> records imported...\r" << flush;
         }//if
      }//for_i
      if (XManager::fgVerbose) {
         cout << "   <" << i << "> records imported...Finished" << endl;
      }//if

//?? or: ??unittree
   // Add tree info to tree
      AddSchemeTreeInfo(schemetree, schemetree->GetName());

   // Write scheme tree to file
      if ((err = WriteTree(schemetree, TObject::kOverwrite)) == errNoErr) {
         // add tree header to list
         AddTreeHeader(schemetree->GetName(), 0);
      }//if
      //delete tree from memory
      schemetree->Delete("");
      schemetree = 0;
   }//if

// Clean up
   SafeDelete(unit);
   SafeDelete(scheme);

   for (i=0; i<fNRows; i++) {
      if (arrXY[i]) {delete [] arrXY[i]; arrXY[i] = 0;}
   }//for_i
   delete [] arrXY;

   if (index)     {delete [] index;     index     = 0;}
   if (arrUnitID) {delete [] arrUnitID; arrUnitID = 0;}
   if (arrX)      {delete [] arrX;      arrX      = 0;}
   if (arrY)      {delete [] arrY;      arrY      = 0;}
   if (arrAtom)   {delete [] arrAtom;   arrAtom   = 0;}
   if (arrPLen)   {delete [] arrPLen;   arrPLen   = 0;}
   if (arrPBase)  {delete [] arrPBase;  arrPBase  = 0;}
   if (arrTBase)  {delete [] arrTBase;  arrTBase  = 0;}
   if (arrMask)   {delete [] arrMask;   arrMask   = 0;}

   return err;
}//ReadData


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XGenomeChip                                                          //
//                                                                      //
// Class for Affymetrix HuGene arrays                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XGenomeChip::XGenomeChip()
            :XGeneChip()
{
   // Default GenomeChip constructor
   if(kCS) cout << "---XGenomeChip::XGenomeChip(default)------" << endl;

   fAffxNames = 0;
}//Constructor

//______________________________________________________________________________
XGenomeChip::XGenomeChip(const char *name, const char *title)
            :XGeneChip(name, title)
{
   // Normal GenomeChip constructor
   if(kCS) cout << "---XGenomeChip::XGenomeChip------" << endl;

   fAffxNames = new TList();
}//Constructor

//______________________________________________________________________________
XGenomeChip::XGenomeChip(const XGenomeChip &chip) 
            :XGeneChip(chip)
{
   // GenomeChip copy constructor
   if(kCS) cout << "---XGenomeChip::XGenomeChip(copy)------" << endl;

}//CopyConstructor

//______________________________________________________________________________
XGenomeChip& XGenomeChip::operator=(const XGenomeChip& rhs)
{
   // GenomeChip assignment operator.

   if (this != &rhs) {
      XGeneChip::operator=(rhs);
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XGenomeChip::~XGenomeChip()
{
   // GenomeChip destructor
   if(kCS) cout << "---XGenomeChip::~XGenomeChip------" << endl;

   if (fAffxNames) {fAffxNames->Delete(); delete fAffxNames; fAffxNames = 0;}
}//Destructor

//______________________________________________________________________________
void XGenomeChip::AddGenomeTreeInfo(TTree *tree, const char *name, Option_t *option,
                  Int_t nctrls, Int_t naffx, Int_t nsubs, Int_t mincells, Int_t maxcells)
{
   // Add tree info to list fUserInfo of tree
   if(kCS) cout << "------XGenomeChip::AddGenomeTreeInfo------" << endl;

   XGenomeTreeInfo *info = new XGenomeTreeInfo(name, "");

   // store class, and name and class of treeset
   info->SetTitle(info->ClassName());
   info->SetOption(option);
   info->SetTreeSetName(GetName());
   info->SetTreeSetClass(ClassName());

   // add user info (using this->AddUserInfo())
   info->AddUserInfo(nctrls, naffx, nsubs, mincells, maxcells);

   tree->GetUserInfo()->Add(info);
}//AddGenomeTreeInfo

//______________________________________________________________________________
Int_t XGenomeChip::ExportLayoutTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in scheme layout tree to file output
   if(kCS) cout << "------XGenomeChip::ExportLayoutTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasX = kFALSE;
   Bool_t hasY = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasX = kTRUE;
      hasY = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fX") == 0) {hasX = kTRUE;}
         if (strcmp(name,"fY") == 0) {hasY = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get layout tree for this chip
   XLayout *layout = 0;
   fTree = (TTree*)gDirectory->Get((names[0]).Data());
   if (!fTree) return errGetTree;
   fTree->SetBranchAddress("CxyBranch", &layout);

   Int_t nentries = (Int_t)(fTree->GetEntries());
   Int_t size     = fNRows*fNCols;
   if (nentries != size) {
      cout << "Warning: Number of entries <" << nentries << "> is not equal to rows*cols <"
           << fNRows*fNCols << ">." << endl;
   }//if

// Output header
   output << "PROBE_ID";
   if (hasX) output << sep << "X";
   if (hasY) output << sep << "Y";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << layout->GetProbeID();
      if (hasX) output << sep << layout->GetX();
      if (hasY) output << sep << layout->GetY();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportLayoutTree

//______________________________________________________________________________
Int_t XGenomeChip::ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                  ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output
   // varlist can be at most: "fProbeLen:fMask" or "*"
   if(kCS) cout << "------XGenomeChip::ExportSchemeTree------" << endl;

   Int_t err = errNoErr;

   err = XDNAChip::ExportSchemeTree(n, names, varlist, output, sep);

   return err;
}//ExportSchemeTree

//______________________________________________________________________________
Int_t XGenomeChip::ExportTransAnnotTree(Int_t n, TString *names, const char *varlist,
                   ofstream &output, const char *sep)
{
   // Export data stored in transcript annotation tree to file output
   if(kCS) cout << "------XGenomeChip::ExportTransAnnotTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasTran = kFALSE;  //transcript_cluster_id
   Bool_t hasName = kFALSE;  //gene name 
   Bool_t hasSymb = kFALSE;  //gene symbol
   Bool_t hasAccn = kFALSE;  //mRNA accession
   Bool_t hasEntr = kFALSE;  //entrez ID
   Bool_t hasChro = kFALSE;  //chromosome
   Bool_t hasCyto = kFALSE;  //cytoband
   Bool_t hasStar = kFALSE;  //start position
   Bool_t hasStop = kFALSE;  //stop position
   Bool_t hasStrd = kFALSE;  //strand
   Bool_t hasXHyb = kFALSE;  //crosshyb_type
   Bool_t hasPTyp = kFALSE;  //probeset type

   if (strcmp(varlist,"*")  == 0) {
      hasTran = kTRUE;
      hasName = kTRUE;
      hasSymb = kTRUE;
      hasAccn = kTRUE;
      hasEntr = kTRUE;
      hasChro = kTRUE;
      hasCyto = kTRUE;
      hasStar = kTRUE;
      hasStop = kTRUE;
      hasStrd = kTRUE;
      hasXHyb = kTRUE;
      hasPTyp = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fTranscriptID") == 0) {hasTran = kTRUE;}
         if (strcmp(name,"fName")         == 0) {hasName = kTRUE;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymb = kTRUE;}
         if (strcmp(name,"fAccession")    == 0) {hasAccn = kTRUE;}
         if (strcmp(name,"fEntrezID")     == 0) {hasEntr = kTRUE;}
         if (strcmp(name,"fChromosome")   == 0) {hasChro = kTRUE;}
         if (strcmp(name,"fCytoBand")     == 0) {hasCyto = kTRUE;}
         if (strcmp(name,"fStart")        == 0) {hasStar = kTRUE;}
         if (strcmp(name,"fStop")         == 0) {hasStop = kTRUE;}
         if (strcmp(name,"fStrand")       == 0) {hasStrd = kTRUE;}
         if (strcmp(name,"fCrossHybType") == 0) {hasXHyb = kTRUE;}
         if (strcmp(name,"fProbesetType") == 0) {hasPTyp = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get transcript annotation tree
   XGenomeAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get((names[0].Data())); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

// Output header
   output << "UNIT_ID";
   if (hasTran) output << sep << "TranscriptClusterID";
   if (hasName) output << sep << "GeneName";
   if (hasSymb) output << sep << "GeneSymbol";
   if (hasAccn) output << sep << "GeneAccession";
   if (hasEntr) output << sep << "EntrezID";
   if (hasChro) output << sep << "Chromosome";
   if (hasCyto) output << sep << "Cytoband";
   if (hasStar) output << sep << "Start";
   if (hasStop) output << sep << "Stop";
   if (hasStrd) output << sep << "Strand";
   if (hasXHyb) output << sep << "CrossHybridization";
   if (hasPTyp) output << sep << "ProbesetType";
   output << endl;

// Export selected variables
   Int_t nentries = (Int_t)(anntree->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      anntree->GetEntry(i);
      output << annot->GetUnitID();
      if (hasTran) output << sep << annot->GetTranscriptID();
      if (hasName) output << sep << annot->GetName();
      if (hasSymb) output << sep << annot->GetSymbol();
      if (hasAccn) output << sep << annot->GetAccession();
      if (hasEntr) output << sep << annot->GetEntrezID();
      if (hasChro) output << sep << annot->GetChromosome();
      if (hasCyto) output << sep << annot->GetCytoBand();
      if (hasStar) output << sep << annot->GetStart();
      if (hasStop) output << sep << annot->GetStop();
      if (hasStrd) output << sep << annot->GetStrand();
      if (hasXHyb) output << sep << annot->GetCrossHybType();
      if (hasPTyp) output << sep << ProbesetTypeID2Name(annot->GetProbesetType());
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportTransAnnotTree

//______________________________________________________________________________
Int_t XGenomeChip::ImportLayout(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import layout information from infile and store as layout tree "tree.cxy"
   if(kCS) cout << "------XGenomeChip::ImportLayout------" << endl;

   char nextline[kBufSize];
   Int_t err = errNoErr;

   TString clf_format_version;
   TString chip_type;
   TString lib_set_name, lib_set_version;
   TString order;

// Check for line "#%clf_format_version"
   while (strncmp("#%clf_format_version", nextline, 20) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   clf_format_version = strtok(&nextline[21], sep);

// Check for line(s) "#%chip_type"
   // reset input file to start
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%chip_type", nextline, 11) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;

      chip_type = strtok(&nextline[12], sep);
   }//while

// Check for line "#%lib_set_name"
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%lib_set_name", nextline, 14) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_name = strtok(&nextline[15], sep);
//TO DO: check if lib_set_name is equal to fName (or fName.Contains())

// Check for line "#%lib_set_version"
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%lib_set_version", nextline, 17) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_version = strtok(&nextline[18], sep);

// Check for line "#%rows": fNRows
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%rows", nextline, 6) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[7]), "%d", &fNRows);

// Check for line "#%cols": fNCols
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%cols", nextline, 6) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[7]), "%d", &fNCols);

// Check for line "#%sequential": fSequential
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%sequential", nextline, 12) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   sscanf(&(nextline[13]), "%hd", &fSequential);

// Check for line "#%order"
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%order", nextline, 7) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   order = strtok(&nextline[8], sep);
//TO DO:   if (strcmp(order,"col_major") == 0) {} else {}

// Check for line "#%header0"
   while (strncmp("#%header0", nextline, 9) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   if (strcmp(strtok(&nextline[10], sep),"probe_id") != 0) {
      cerr << "Error: Column <probe_id> is missing." << endl;
      return errReadingInput;
   }//if

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Create new layout tree
   fCxyTreeName = TString(fName) + "." + exten;
   TTree *layouttree = new TTree(fCxyTreeName, "layout information");
   if (layouttree == 0) return errCreateTree;
   XLayout *layout = 0;
   layout = new XLayout();
   layouttree->Branch("CxyBranch", "XLayout", &layout, 64000, split);

//   Long_t probeid;
   Int_t probeid, x, y;
   Int_t idx = 0;
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) break;
      if (!input.good()) {err = errPrematureEOF; break;}
      sscanf(nextline,"%d %d %d", &probeid, &x, &y);

      // fill layout tree with control data
      layout->SetProbeID(probeid);
      layout->SetX(x);
      layout->SetY(y);
      layouttree->Fill();

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx + 1 << "> records imported...\r" << flush;
      }//if
      idx++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records imported...Finished" << endl;
   }//if

   if (idx != fNRows*fNCols) {
      cout << "Warning: Number of entries <" << idx << "> is not equal to <"
           << fNRows*fNCols << ">." << endl;
   }//if

// Add tree info to layouttree
   AddLayoutTreeInfo(layouttree, fCxyTreeName);

// Write layout tree to file
   if ((err = WriteTree(layouttree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(layouttree->GetName(), 0);
   }//if
   //delete tree from memory
   layouttree->Delete("");
   layouttree = 0;

// Clean up
   SafeDelete(layout);

   return err;
}//ImportLayout

//______________________________________________________________________________
Int_t XGenomeChip::ImportProbeInfo(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import probe information from infile and store as probe tree
   // Note: Not necessary since info is available in *.pfg" file!
   if(kCS) cout << "------XGenomeChip::ImportProbeInfo------" << endl;

// NOT USED - to prevent compiler warnings:
   input.eof(); option = 0; sep = 0; delim = 0; split = 0;

   cout << "Note: Import of probe info for Genome and Exon Arrays is not necessary."
        << endl;
   return errNoErr;
}//ImportProbeInfo

//______________________________________________________________________________
Int_t XGenomeChip::ReadHeader(ifstream &input, const char *sep, char delim)
{
   // Read header from input PGF file. 
   if(kCS) cout << "------XGenomeChip::ReadHeader------" << endl;

   char  nextline[kBufSize];
   Int_t err = errNoErr;

   TString str;
   Int_t   begin, end;

   TString pgf_format_version;
   TString chip_type;
   TString lib_set_name, lib_set_version;

// Check for line "#%pgf_format_version"
   while (strncmp("#%pgf_format_version", nextline, 20) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   pgf_format_version = strtok(&nextline[21], sep);

// Check for line(s) "#%chip_type"
   // reset input file to start
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%chip_type", nextline, 11) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;

      chip_type = strtok(&nextline[12], sep);
   }//while

// Check for line "#%lib_set_name"
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%lib_set_name", nextline, 14) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_name = strtok(&nextline[15], sep);

// Check for line "#%lib_set_version"
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%lib_set_version", nextline, 17) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_version = strtok(&nextline[18], sep);
//TO DO? check lib_set_version

// Check for line "#%header0"
   input.clear(); input.seekg(0, ios::beg);
   while (strncmp("#%header0", nextline, 9) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   if (strcmp(strtok(&nextline[10], sep),"probeset_id") != 0) {
      cerr << "Error: Header for <probeset_id> is missing." << endl;
      return errReadingInput;
   }//if

// Check for line "#%header1"
   while (strncmp("#%header1", nextline, 9) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   str = RemoveEnds(&nextline[10], begin, end);
   if (strcmp(str.Data(), "atom_id") != 0) {
      cerr << "Error: Header for <atom_id> is missing." << endl;
      return errReadingInput;
   }//if

// Check for line "#%header2"
   while (strncmp("#%header2", nextline, 9) != 0) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) return errPrematureEOF;
   }//while

   str = RemoveEnds(&nextline[10], begin, end);
   str = strtok((char*)str.Data(), sep);
   if (strcmp(str.Data(), "probe_id") != 0) {
      cerr << "Error: Header for <probe_id> is missing." << endl;
      return errReadingInput;
   }//if

   return err;
}//ReadHeader

//______________________________________________________________________________
Int_t XGenomeChip::ReadData(ifstream &input, Option_t *option,
                   const char *sep, char delim, Int_t split)
{
   // Import scheme data from input PGF file and store in the following trees:
   if(kCS) cout << "------XGenomeChip::ReadData------" << endl;

   char  nextline[kBufSize];
   Int_t err = errNoErr;

// Get transcript annotation tree
   TString treename = TString(fName) + "." + kExtenAnnot[0];
   XGenomeAnnotation *ann = 0;
   TTree *anntree = (TTree*)gDirectory->Get(treename.Data());
   if (anntree == 0) {
      cerr << "Error: Missing transcript annotation tree <" << treename.Data() << ">."
           << endl;
      return errGetTree;
   }//if
   anntree->SetBranchAddress("AnnBranch", &ann);

   Int_t nentries = anntree->GetEntries();

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Create new scheme tree
   fScmTreeName      = TString(fName) + "." + exten;
   TTree *schemetree = new TTree(fScmTreeName, "scheme information");
   if (schemetree == 0) return errCreateTree;
   XScheme *scheme = 0;
   scheme = new XScheme();
   schemetree->Branch("ScmBranch", "XScheme", &scheme, 64000, split);

// Create new unit tree
   fIdxTreeName    = TString(fName) + "." + TString(kExtenScheme[2]);
   TTree *unittree = new TTree(fIdxTreeName, "unit information");
   if (unittree == 0) return errCreateTree;
   XGCUnit *unit = 0;
   unit = new XGCUnit();
   unittree->Branch("IdxBranch", "XGCUnit", &unit, 64000, split);

// Create new probe tree
   fPrbTreeName     = TString(fName) + "." + TString(kExtenScheme[3]);
   TTree *probetree = new TTree(fPrbTreeName, "probe info for scheme");
   if (probetree == 0) return errCreateTree;
   XGCProbe *probe = 0;
   probe = new XGCProbe();
   probetree->Branch("PrbBranch", "XGCProbe", &probe, 64000, split);

   TString geneaccess, genesymbol, mrnaaccess, seqname;
   Int_t   translevel;

   TString order;
   TString str;
   Int_t   size;  //probeset_count
   Int_t   begin, end;
   Int_t   probeset_id;
   TString probeset_type;

   Int_t   probe_id, gc_count, probe_len, inter_pos;
   TString probe_type, sequence;
   Int_t   x, y;

   Int_t    unitID = 0;   //???replace with ??
   Int_t    mask   = 0;   //???replace with ??
   Int_t    idx;
   Int_t    numcells, numatoms;
   Double_t Tm;
   Short_t  strd;
   Int_t    geneID;

// Init counting variables
   Int_t    scmcount = 0;       //count number of scheme data (<=fNRows*fNCols)
   Int_t    idxcount = 0;       //count number of units
   Int_t    ctlcount = 0;       //count number of negative controls (EProbesetType<0)
   Int_t    maxunits = 0;       //maximum number of unit cells
   Int_t    minunits = 9999999; //minimum number of unit cells

   Bool_t   didWhile = kFALSE;  //to avoid circular for-loop

// Init local arrays
   Long_t  *pos = 0;  //positions of probeset_ids
   Int_t   *psi = 0;  //probeset_id
   Int_t   *pst = 0;  //probeset_type as EProbesetType
   Short_t *msk = 0;  //mask: set msk=1 for probeset_ids already read
   Long_t  *arr = 0;  //temporary array to store other arrays

// Init local arrays for sorting
   Int_t    *index = 0;  //sort index
   Long64_t *value = 0;  //value array to store the two columns to sort
   Long64_t majorv, minorv;  //values to store in array value

// Init local arrays for sorted data from transcript annotation tree
   Int_t   *arrGene  = 0;  //array for transcript_cluster_id
   Short_t *arrXType = 0;  //array for crosshyb_type (unique, similar, mixed)

// Get position of first probeset_id after header lines
   streampos position = input.tellg();

// NEED to parse "category" for scheme mask: 
//  1 = unique  => eMETACORE
//  2 = similar (does not exist)
//  3 = mixed   => ePARACORE
// => use SchemeMask(Short_t xtype) analogously to exon array

//--------------------------------------------------------------------
// Read data and sort for probeset_type and position of probeset_id
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Reading data from input file..." << endl;
   }//if

// Read data to get number of probesets (header0)
   size = 0;
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) break;
      if (!input.good()) {err = errPrematureEOF; break;}
      str = RemoveEnds(&nextline[0], begin, end);

      if (begin==0) size++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   Number of probesets is <" << size << ">." << endl;
   }//if

// Check if number of annotated probesets is equal to probeset_count
   if (fNProbesets != size) {
      if (XManager::fgVerbose) {
         cout << "   Note: Number of annotated probesets <" << fNProbesets  
              << "> is not equal to number of probesets <" << size << ">." << endl;
      }//if
      // set final number of probesets
      fNProbesets = size;
   }//if

// Initialize memory for local arrays
   if (!(pos = new (nothrow) Long_t[size]))  {err = errInitMemory; goto cleanup;}
   if (!(psi = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(pst = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(msk = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize memory for sorting arrays
   if (!(index = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
   if (!(value = new (nothrow) Long64_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize array pst to eINIT since eCONTROLCHIP=0!
//not necessary???   for (Int_t i=0; i<size; i++) pst[i] = eINIT;

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Read data to get probeset_id and type (header0)
   idx = 0;
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) break;
      if (!input.good()) {err = errPrematureEOF; break;}
      str = RemoveEnds(&nextline[0], begin, end);

      // get position, and probeset_id and type (header0)
      if (begin==0) {
         pos[idx] = (Long_t)input.tellg() - (Long_t)(str.Length() + end + 1);
         psi[idx] = atoi(strtok((char*)str.Data(), sep));
         pst[idx] = ProbesetType(strtok(0, sep));
         msk[idx] = 0;

         // fill sort values
         majorv = (Long64_t)pst[idx];  //sort for probeset_type first
         minorv = (Long64_t)pos[idx];  //then sort for position of probeset_id
         value[idx]  = majorv << 31;
         value[idx] += minorv;

         if (XManager::fgVerbose && idx%10000 == 0) {
            cout << "   <" << idx + 1 << "> records read...\r" << flush;
         }//if
         idx++;
      }//if
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records read...Finished" << endl;
   }//if

   if (fNProbesets != idx) {
      cerr << "Error: Number of probeset_ids read <" << idx  
           << "> is not equal to number of probesets <" << fNProbesets << ">" << endl;
      err = errReadingInput;
      goto cleanup;
   }//if

// Sort for probeset_type first and then for position of probeset_id
   if (XManager::fgVerbose) {
      cout << "   Sorting data for probeset_type and position..." << endl;
   }//if
   TMath::Sort(size, value, index, 0);

   // delete value here to free memory
   if (value) {delete [] value; value = 0;}

// Replace data in arrays with sorted data using temporary array arr
   // create temporary array (one temporary array only to save memory)
   if (!(arr = new (nothrow) Long_t[size])) {err = errInitMemory; goto cleanup;}

   // replace probeset_id positions with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pos[index[i]];
   for (Int_t i=0; i<size; i++) pos[i] = arr[i];

   // replace probeset_ids with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = psi[index[i]];
   for (Int_t i=0; i<size; i++) psi[i] = arr[i];

   // replace probeset_types with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pst[index[i]];
   for (Int_t i=0; i<size; i++) pst[i] = arr[i];

   // delete temporary array here to free memory
   if (arr) {delete [] arr; arr = 0;}

// count number of negative controls (EProbesetType<0)
   for (Int_t i=0; i<size; i++) {
      if (pst[i] < eCONTROLCHIP) ctlcount++;
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   Total number of controls is <" << ctlcount << ">" << endl;
   }//if


//--------------------------------------------------------------------
// 1. Fill tree(s) with data for probeset type: control->chip
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
//      cout << "Filling trees with data for probeset type: control->chip..." << endl;
      cout << "   Note: no data for probeset type: control->chip..." << endl;
   }//if


//--------------------------------------------------------------------
// 2. Fill tree(s) with data for probeset types:
//    normgene->intron, normgene->exon (not main!), rescue->FLmRNA->unmapped
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: normgene, rescue..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for probe_types of normgene, rescue
   for (Int_t i=0; i<size; i++) {
      if (pst[i] > eINTRON) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break; //ev goto cleanup;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            // control->chip: blank
            if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = SchemeMask(pst[i], 0);

            // fill scheme tree
            scheme->SetUnitID(-ctlcount);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if(scmcount%10000 == 0) cout << " <" << scmcount + 1 << "> records imported...\r" << flush;
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(-ctlcount);
      unit->SetUnitType(pst[i]);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();

      ctlcount--;
      idxcount++;
   }//for_i


//--------------------------------------------------------------------
// 3. Fill tree(s) with data for probeset types: 
//    control->bgp->antigenomic, control->bgp->genomic
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: control->bgp..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for probe_types of antigenomic, genomic
   for (Int_t i=0; i<size; i++) {
      if (!((pst[i] == eANTIGENOMIC) || (pst[i] == eGENOMIC))) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            // control->chip: blank
            if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = SchemeMask(pst[i], 0);

            // fill scheme tree
            scheme->SetUnitID(-ctlcount);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if (XManager::fgVerbose && scmcount%10000 == 0) {
               cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
            }//if
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(-ctlcount);
      unit->SetUnitType(pst[i]);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();

      ctlcount--;
      idxcount++;
   }//for_i

// Set number of controls
   fNControls = idxcount;


//--------------------------------------------------------------------
// 4. Fill tree(s) with data for probeset type: control->affx
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: control->affx..." << endl;
   }//if

// Sort for probeset_id (necessary for BinarySearch below)
   TMath::Sort(size, psi, index, 0);

// Replace data in arrays with sorted data using temporary array arr
   // create temporary array (one temporary array only to save memory)
   if (!(arr = new (nothrow) Long_t[size])) {err = errInitMemory; goto cleanup;}

   // replace probeset_id positions with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pos[index[i]];
   for (Int_t i=0; i<size; i++) pos[i] = arr[i];

   // replace probeset_ids with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = psi[index[i]];
   for (Int_t i=0; i<size; i++) psi[i] = arr[i];

   // replace probeset_types with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pst[index[i]];
   for (Int_t i=0; i<size; i++) pst[i] = arr[i];

   // replace msk with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = msk[index[i]];
   for (Int_t i=0; i<size; i++) msk[i] = arr[i];

   // delete temporary array here to free memory
   if (arr) {delete [] arr; arr = 0;}

// Initialize memory for anntree arrays
   if (!(arrGene  = new (nothrow) Int_t[nentries]))   {err = errInitMemory; goto cleanup;}
   if (!(arrXType = new (nothrow) Short_t[nentries])) {err = errInitMemory; goto cleanup;}

// Fill arrays with anntree entries
   for (Int_t i=0; i<nentries; i++) {
      anntree->GetEntry(i);
      arrGene[i]  = atoi(ann->GetTranscriptID());
      arrXType[i] = ann->GetCrossHybType();
      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> probeset tree entries read...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << nentries << "> probeset tree entries read...Finished" << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for probe_types of bgp, normgene, rescue
   unitID = 0;
   for (Int_t i=0; i<fNAffx; i++) {
      geneID = arrGene[i];

      Int_t z = (Int_t)(TMath::BinarySearch((Long64_t)size, psi, arrGene[i]));
      if (psi[z] != arrGene[i]) continue;

      input.seekg(pos[z], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[z] != probeset_id) {
         cerr << "Error: Probeset_id <" << psi[z] << "> is not <" << probeset_id << ">."
              << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[z] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = SchemeMask(eCONTROLAFFX, strd);

            // fill scheme tree
            scheme->SetUnitID(unitID);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if (XManager::fgVerbose && scmcount%10000 == 0) {
               cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
            }//if
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(unitID);
      unit->SetUnitType(mask);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();

      // minimal/maximal number of unit cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      maxunits = (numcells >  maxunits) ? numcells : maxunits;

      unitID++;
      idxcount++;
   }//for_i

// Number of AFFX controls
   if (XManager::fgVerbose) {
      cout << "   Number of control->affx probesets is <" << unitID << ">." << endl;
   }//if

// Check if total number of AFFX controls is equal to fNAffx
   if ((fNAffx > 0) && (fNAffx != unitID)) {
      cerr << "Error: Number of control->affx imported <" << unitID  
           << "> is not equal to number of annotated AFFX controls <" << fNAffx << ">."
           << endl;
//??      fNAffx = unitID;
      err = errCDFVersion;
      goto cleanup;
   } else if (fNAffx == 0) {
      cout << "Warning: Missing annotation for control->affx probesets." << endl;
      fNAffx = unitID;
   }//if


//--------------------------------------------------------------------
// 5. Fill tree(s) with data for probeset type: main
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: main..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

   // control->affx are already imported
   for (Int_t i=fNAffx; i<nentries; i++) {
      geneID = arrGene[i];

      Int_t z = (Int_t)(TMath::BinarySearch((Long64_t)size, psi, arrGene[i]));
      if (psi[z] != arrGene[i] || pst[z] < eMAIN) continue;

      numcells = 0;
      numatoms = 0;
      didWhile = kFALSE;
      translevel = eNOLEVEL;
      while (geneID == arrGene[i]) {

         input.seekg(pos[z], ios::beg);
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {err = errPrematureEOF; goto cleanup;}

         probeset_id   = atoi(strtok(&nextline[0], sep));
         probeset_type = TString(strtok(NULL, sep));

         if (psi[z] != probeset_id) {
            cerr << "Error: Probeset_id <" << psi[z] << "> is not <" << probeset_id
                 << ">." << endl;
            err = errReadingInput;
            break; //ev goto cleanup;
         }//if
         msk[z] = 1;

         begin = 1;
         while (begin != 0) { //loop over header0
            input.getline(nextline, kBufSize, delim);
            if (input.eof()) {
               // last line of *.pfg: need to reset input file to continue
               input.clear();  //clear all flags
               input.seekg(position, ios::beg);
               break;
            }//if
            str = RemoveEnds(&nextline[0], begin, end);

            if (begin == 1) numatoms++; //count atoms for header1
            if (begin == 2) { //data for header2
               probe_id   = atoi(strtok((char*)str.Data(), sep));
               probe_type = TString(strtok(NULL, sep));
               gc_count   = atoi(strtok(NULL, sep));
               probe_len  = atoi(strtok(NULL, sep));
               inter_pos  = atoi(strtok(NULL, sep));
               sequence   = TString(strtok(NULL, sep));

               // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
               x = Index2X(probe_id - 1);
               y = Index2Y(probe_id - 1);
//to do               x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_i);
//to do               y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

               // control->chip: blank
               if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

               Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
               strd = ProbeType(probe_type);
               mask = SchemeMask(arrXType[i]);

               // fill scheme tree
               scheme->SetUnitID(unitID);
               scheme->SetX(x);
               scheme->SetY(y);
               scheme->SetProbeLength(probe_len);
               scheme->SetMask(mask);
               schemetree->Fill();

              // fill probe tree: (x,y) has identical order to (x,y) of schemetree
               probe->SetX(x);
               probe->SetY(y);
               probe->SetSequence(sequence);
               probe->SetPosition(inter_pos);
               probe->SetNumberGC(gc_count);
               probe->SetTMelting(Tm);
               probe->SetProbeType(strd);
               probetree->Fill();

               if (XManager::fgVerbose && scmcount%10000 == 0) {
                  cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
               }//if
               numcells++;
               scmcount++;
            }//if
         }//while

         translevel = (mask > translevel) ? mask : translevel;

         didWhile = kTRUE;
         i++;

//?? ev at beginning of while with:  if(didWhile) ?? see XExonChip
         z = (Int_t)(TMath::BinarySearch((Long64_t)size, psi, arrGene[i]));
         if (psi[z] != arrGene[i]) continue;
      }//while

      // fill transcript unit tree
      str = TString(""); str += geneID;
      unit->SetUnitName(str);
      unit->SetUnitID(unitID);
      unit->SetUnitType(translevel);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();

      // minimal/maximal number of unit cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      maxunits = (numcells >  maxunits) ? numcells : maxunits;

      unitID++;
      idxcount++;
      if (didWhile) i--;
   }//for_i

// Set number of genes
   fNGenes = unitID - fNAffx;


//--------------------------------------------------------------------
// 6. Fill tree(s) with data for non-annotated probesets
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for non-annotated probesets..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for non-annotated probeset_ids
//!! fill tree.scm and tree.prb with rest of probe_set data!
  for (Int_t i=0; i<size; i++) {
      if (msk[i] == 1) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            // control->chip: blank
            if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = SchemeMask(pst[i], strd);

            // fill scheme tree
            scheme->SetUnitID(eUNKNOWNTYPE);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if (XManager::fgVerbose && scmcount%10000 == 0) {
               cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
            }//if
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(eUNKNOWNTYPE);
      unit->SetUnitType(eNOLEVEL);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unittree->Fill();

      // minimal/maximal number of unit cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      maxunits = (numcells >  maxunits) ? numcells : maxunits;

      idxcount++;
   }//for_i

// Set number of probes
   if (XManager::fgVerbose) {
      cout << "   <" << scmcount << "> records imported...Finished" << endl;
   }//if
   fNProbes = scmcount;

// Set number of transcript units
   if (XManager::fgVerbose) {
      cout << "   <" << idxcount << "> total transcript units imported." << endl;
   }//if
   fNUnits = idxcount;

   if (XManager::fgVerbose) {
      cout << "   Genome cell statistics: " << endl;
      cout << "      Number of unit cells: minimum = " << minunits << ",  maximum = " << maxunits << endl;
   }//if

//?? or to unittree
// Add tree info to tree
   AddSchemeTreeInfo(schemetree, schemetree->GetName());

// Write scheme tree to file
   if ((err = WriteTree(schemetree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(schemetree->GetName(), 0);
   }//if
   //delete tree from memory
   schemetree->Delete("");
   schemetree = 0;

// Add tree info to unittree
   this->AddGenomeTreeInfo(unittree, unittree->GetName(), "", fNControls, fNAffx,
                           fNGenes, minunits, maxunits);

// Write unit tree to file
   if ((err = WriteTree(unittree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(unittree->GetName(), 0);
   }//if
   //delete tree from memory
   unittree->Delete("");
   unittree = 0;

// Add tree info to probetree
   AddProbeTreeInfo(probetree, fPrbTreeName);

// Write probe tree to file
   if ((err = WriteTree(probetree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(probetree->GetName(), 0);
   }//if
   //delete tree from memory
   probetree->Delete("");
   probetree = 0;

//Clean up
cleanup:
   if (arrXType) {delete [] arrXType; arrXType = 0;}
   if (arrGene)  {delete [] arrGene;  arrGene  = 0;}

   if (value) {delete [] value; value = 0;}
   if (index) {delete [] index; index = 0;}

   if (arr) {delete [] arr; arr = 0;}
   if (msk) {delete [] msk; msk = 0;}
   if (pst) {delete [] pst; pst = 0;}
   if (psi) {delete [] pst; pst = 0;}
   if (pos) {delete [] pos; pos = 0;}

   SafeDelete(probe);
   SafeDelete(unit);
   SafeDelete(scheme);

   return err;
}//ReadData

//______________________________________________________________________________
Int_t XGenomeChip::ImportTransAnnotation(ifstream &input, Option_t *option, 
                   const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree "tree.ann":
   if(kCS) cout << "------XGenomeChip::ImportTransAnnotation------" << endl;

   Int_t err  = errNoErr;
   Int_t size = 0;
   Int_t idx  = 0;
   Int_t nsub = 5;  //number of sub-fields

   std::string nextline;  //to read next line
   streampos position;

   TString lib_set_name, lib_set_version;
   TString genome_species, genome_version;
   TString annot_type;
   TString assigngene, assignmrna, category;
   TString str, dummy;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Init variables to tokenize lines
   const char *csv = "\",\""; //necessary to tokenize lines, better?: sep = "\",\"";
   const char *tab = "\t";    //replace csv with tab

// Init local arrays to store data from input columns
   Int_t   *geneid     = 0;  //transcript_cluster_id
   TString *seqname    = 0;  //seqname, i.e. chromosome
   TString *strand     = 0;  //strand of chromosome
   Int_t   *start      = 0;  //start position on chromosome
   Int_t   *stop       = 0;  //stop position on chromosome
   Short_t *crosshtype = 0;  //crosshyb_type (unique, similar, mixed)
//   Int_t   *nprobes    = 0;  //total number of probes
   Int_t    nprobes;

// Init local arrays for derived data
   TString *names      = 0;  //array containing tokenized string
   TString *geneaccess = 0;  //gene_assignment, accession number
   TString *genename   = 0;  //gene_assignment, name
   TString *genesymbol = 0;  //gene_assignment, symbol
   TString *cytoband   = 0;  //gene_assignment, cytoband
   Int_t   *entrezid   = 0;  //gene_assignment, entrez ID
   TString *mrnaaccess = 0;  //mrna_assignment, accession number
   Short_t *psettype   = 0;  //category (main, normgene->intron, control->affx, etc)

// Init local arrays for sorting
   Int_t    *index     = 0;  //sort index
   Long64_t *value     = 0;  //value array to store the two columns to sort
   Long64_t majorv, minorv;  //values to store in array value

// Create new transcript annotation tree
   fAnnTreeName = TString(fName) + "." + exten;
   TTree *anntree = new TTree(fAnnTreeName, "transcript annotation");
   if (anntree == 0) return errCreateTree;
   XGenomeAnnotation *ann = 0;
   ann = new XGenomeAnnotation();
   anntree->Branch("AnnBranch", "XGenomeAnnotation", &ann, 64000, split);

// Check for line "#%lib-set-name"
   while ((nextline.compare(0, 14, "#%lib-set-name") != 0) &&
          (nextline.compare(0, 14, "#%lib_set_name") != 0)) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_name = strtok((&((char*)nextline.c_str())[15]), sep);

// Check for line "#%lib-set-version"
   input.clear();
   input.seekg(position, ios::beg);
   while ((nextline.compare(0, 17, "#%lib-set-version") != 0) &&
          (nextline.compare(0, 17, "#%lib_set_version") != 0)) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_version = strtok((&((char*)nextline.c_str())[18]), sep);

// Check for line "#%genome-species"
   Bool_t hasSpecies = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 16, "#%genome-species") != 0) {
      std::getline(input, nextline, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasSpecies = kFALSE; break;}
   }//while
   if (hasSpecies == kTRUE) {
      genome_species = strtok((&((char*)nextline.c_str())[17]), sep);
   } else {
      cout << "   Note: The following header line is missing: %genome-species=" << endl;
      genome_species = "NA";
   }//if

// Check for line "#%genome-version"
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 16, "#%genome-version") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   genome_version = strtok((&((char*)nextline.c_str())[17]), sep);

// Check for line "#%netaffx-annotation-netaffx-build"
   Bool_t hasNetBuild = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 34, "#%netaffx-annotation-netaffx-build") != 0) {
      std::getline(input, nextline, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasNetBuild = kFALSE; break;}
   }//while
   if (hasNetBuild == kTRUE) {
      fVersionAnnot = strtok((&((char*)nextline.c_str())[35]), sep);
   } else {
      cout << "   Note: The following header line is missing: %netaffx-annotation-netaffx-build=" << endl;
      fVersionAnnot = "NA";
   }//if

// Check for line "#%netaffx-annotation-data-type"
   Bool_t hasNetType = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 30, "#%netaffx-annotation-data-type") != 0) {
      std::getline(input, nextline, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasNetType = kFALSE; break;}
   }//while
   if (hasNetType == kTRUE) {
      annot_type = strtok((&((char*)nextline.c_str())[31]), sep);
   } else {
      cout << "   Note: The following header line is missing: %netaffx-annotation-data-type=" << endl;
      annot_type = "NA";
   }//if

// Check for line ""transcript_cluster_id""
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 23, "\"transcript_cluster_id\"") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   Int_t numsep = NumSeparators(nextline.c_str(), sep);

   // remove all "\"" from header line
   str = TString(nextline);
   str.ReplaceAll("\"", "");

// Create array hasColumn 
   Int_t *hasColumn = 0;
   if (!(hasColumn = new (nothrow) Int_t[kNTranscriptCols])) return errInitMemory;
   for (Int_t i=0; i<kNTranscriptCols; i++) {
      hasColumn[i] = 0;
   }//for_i

//TO DO: see ImportProbesetAnnotation()
//better   err = CheckHeaderOrder(str, kTranscriptHeader, kNTranscriptCols, hasColumn, sep);
// Check column headers from header line, hasColumn = 1 if column is present
   err = CheckHeader(str, kTranscriptHeader, kNTranscriptCols, hasColumn, sep);
   if (err > 0) {
      cout << "Note: The following header columns are missing: " << endl;
      for (Int_t i=0; i<kNTranscriptCols; i++) {
         if (hasColumn[i] == 0) cout << "<" << kTranscriptHeader[i] << ">" << endl;
         if (i == 0  && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //transcript_cluster_id
         if (i == 2  && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //seqname
         if (i == 3  && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //strand
         if (i == 4  && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //start
         if (i == 7  && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //gene_assignment
         if (i == 16 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //crosshyb_type
         if (i == 17 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //category
      }//for_i
   }//if 

// Get number of transcripts
   // get current streamposition for rewinding input
   position = input.tellg();
   while (1) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      size++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   Number of transcripts is <" << size << ">." << endl;
   }//if

   // set number of annotated probesets (less than probeset_count in ImportScheme())
   fNProbesets = size;

   // reset input file to position
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

///////////////
//TO DO: ???
//PROBLEM: need to sort according to UnitID and UnitName of unittree tree.idx!
//use: Int_t z = (Int_t)(TMath::BinarySearch((Long64_t)size, psi, arrPSet[i]));

// Initialize memory for local arrays
   if (!(geneid     = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(seqname    = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(strand     = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(start      = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(stop       = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
//   if (!(nprobes    = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(crosshtype = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}

   if (!(names      = new (nothrow) TString[nsub])) {err = errInitMemory; goto cleanup;}
   if (!(geneaccess = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(genename   = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(genesymbol = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(cytoband   = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(entrezid   = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(mrnaaccess = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(psettype   = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize memory for sorting arrays
   if (!(index = new (nothrow) Int_t[size]))        {err = errInitMemory; goto cleanup;}
   if (!(value = new (nothrow) Long64_t[size]))     {err = errInitMemory; goto cleanup;}

   while (input.good()) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      if (input.fail()) {
         cout << "Error: Failed reading line:" << endl;
         cout << nextline << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // replace all "\",\"" with tab "\t" and remove "\""
      str = TString(nextline);
      // first need to replace "" with NA
      str.ReplaceAll("\"\"", "\"NA\"");
      // replace tab with space to eliminate wrong tabs in Affymetrix annotation files
      str.ReplaceAll(tab, " ");
      // replace csv with tab
      str.ReplaceAll(csv, tab);
      // remove all "\"" from line
      str.ReplaceAll("\"", "");

      // check number of separators
      if (numsep != NumSeparators(str, tab)) {
         cout << "Error: Wrong number of separators in line:" << endl;
         cout << str.Data() << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // import fields
      dummy           = strtok((char*)str.Data(), tab);
      geneid[idx]     = atoi((RemoveEnds(dummy)).Data());
      dummy           = strtok(NULL, tab);
      seqname[idx]    = strtok(NULL, tab);
      strand[idx]     = strtok(NULL, tab);
      start[idx]      = atoi(strtok(NULL, tab));
      stop[idx]       = atoi(strtok(NULL, tab));
      nprobes         = atoi(strtok(NULL, tab));
//      if (nprobes == 0) continue;
      assigngene      = strtok(NULL, tab);
      assignmrna      = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      crosshtype[idx] = atoi(strtok(NULL, tab));
      category        = strtok(NULL, tab);

      // check seqname and strand
      if (strcmp(seqname[idx].Data(),"---") == 0) {seqname[idx] = "NA";}
      if (strcmp(strand[idx].Data(), "---") == 0) {strand[idx]  = "?";}

      // get gene accession and symbol
      geneaccess[idx] = strtok((char*)assigngene.Data(), "\\");
      if (strcmp(geneaccess[idx].Data(),"---") == 0) {
         geneaccess[idx] = "NA";
         genesymbol[idx] = "NA";
         genename[idx]   = "NA";
         cytoband[idx]   = "NA";
         entrezid[idx]   = -1;
      } else {
         nsub = 5;  //since TokenizeString() can change nsub
         for (Int_t i=0; i<nsub; i++) names[i] = "---";

         Int_t index = 0;
         index = assigngene.Index(kSepSl3, kNumSl3, index, TString::kExact);
         // index>0 only if assigngene is a multipart entry
         assigngene = (index > 0) ? assigngene(0, index) : assigngene;

         index = TokenizeString(assigngene.Data(), nsub, names, kNumSl2, kSepSl2);

         geneaccess[idx] = (strcmp(names[0].Data(),"---") != 0) ? names[0] : "NA";
         genesymbol[idx] = (strcmp(names[1].Data(),"---") != 0) ? names[1] : "NA";
         genename[idx]   = (strcmp(names[2].Data(),"---") != 0) ? names[2] : "NA";
         cytoband[idx]   = (strcmp(names[3].Data(),"---") != 0) ? ("chr" + names[3]) : "NA";
         entrezid[idx]   = (strcmp(names[4].Data(),"---") != 0) ? names[4].Atoi() : -1;
      }//if

      // convert category to probesettype_id
      psettype[idx] = this->ProbesetType(category);

      // get mrna accession for genes or for control->affx 
      if (psettype[idx] == eCONTROLAFFX) {
         mrnaaccess[idx] = RemoveEnds(assignmrna);

         // replace gene accession and symbol for control->affx
         geneaccess[idx] = mrnaaccess[idx];
         geneaccess[idx] += TString("_");
         geneaccess[idx] += (fNAffx + 1);

         // add control->affx names to list fAffxNames
         XIdxString *idxstr = new XIdxString(geneid[idx], geneaccess[idx]);
         fAffxNames->Add(idxstr);  //idxstr will be deleted by fAffxNames

         fNAffx++;
      } else if ((psettype[idx] == eINTRON) ||
                 (psettype[idx] == eEXON)   ||
                 (psettype[idx] == eUNMAPPED)) {
         nsub = 5;  //since TokenizeString() can change nsub
         for (Int_t i=0; i<nsub; i++) names[i] = "---";

         Int_t index = 0;
         index = assignmrna.Index(kSepSl3, kNumSl3, index, TString::kExact);
         // index>0 only if assignmrna is a multipart entry
         assignmrna = (index > 0) ? assignmrna(0, index) : assignmrna;

         index = TokenizeString(assignmrna.Data(), nsub, names, kNumSl2, kSepSl2);

         geneaccess[idx] = (strcmp(names[0].Data(),"---") != 0) ? names[0] : "NA";
         genesymbol[idx] = (strcmp(names[1].Data(),"---") != 0) ? names[1] : "NA";
         genename[idx]   = (strcmp(names[2].Data(),"---") != 0) ? names[2] : "NA";
         cytoband[idx]   = (strcmp(names[3].Data(),"---") != 0) ? ("chr" + names[3]) : "NA";
         entrezid[idx]   = (strcmp(names[4].Data(),"---") != 0) ? names[4].Atoi() : -1;
      } else {
         mrnaaccess[idx] = RemoveEnds(strtok((char*)assignmrna.Data(), "/"));
         if (strcmp(mrnaaccess[idx].Data(),"") == 0) {mrnaaccess[idx] = "NA";}

         // for geneaccess = NA use mrnaaccess
         if (strcmp(geneaccess[idx].Data(),"NA") == 0) {
            geneaccess[idx] = mrnaaccess[idx];
         }//if
      }//if

      // fill sort values
      majorv = (Long64_t)geneid[idx]; //sort for transcript_cluster_id first
      minorv = (Long64_t)start[idx];  //then sort for start position
      value[idx]  = majorv << 31;
      value[idx] += minorv;

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx + 1 << "> records read...\r" << flush;
      }//if
      idx++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records read...Finished" << endl;
   }//if

   // set size to number of records read
   size = idx;

// Sort for geneid first and then for start (assuming transcipt on one chromosome)
   TMath::Sort(size, value, index, 0);

// Fill transcript annotation tree with control->affx 
   for (Int_t i=0; i<fNAffx; i++) {
      XIdxString *idxstr   = (XIdxString*)(fAffxNames->At(i));
      TString     affxname = idxstr->GetString();

      // affxid as string
      str.Form("%d", idxstr->GetIndex());

      // fill annotation tree with data
      ann->SetUnitID(i); 
      ann->SetTranscriptID(str);
      ann->SetName(affxname);
      ann->SetSymbol("NA");
      ann->SetAccession(str);
      ann->SetEntrezID(-1);
      ann->SetChromosome("NA");
      ann->SetCytoBand("NA");
      ann->SetStrand('?');
      ann->SetStart(0);
      ann->SetStop(0);
      ann->SetCrossHybType(0);
      ann->SetProbesetType(eCONTROLAFFX);
      anntree->Fill();
   }//for_i

// Fill transcript annotation tree with sorted "main" annotation data
   idx = fNAffx; //start at fNAffx
   for (Int_t i=0; i<size; i++) {
      Int_t k = index[i];
      if (psettype[k] != eMAIN) continue;

      // geneid as string
      str.Form("%d", geneid[k]);

      // fill annotation tree with data
      ann->SetUnitID(idx++);
      ann->SetTranscriptID(str);
      ann->SetName(genename[k]);
      ann->SetSymbol(genesymbol[k]);
      ann->SetAccession(geneaccess[k]);
      ann->SetEntrezID(entrezid[k]);
      ann->SetChromosome(seqname[k]);
      ann->SetCytoBand(cytoband[k]);
      ann->SetStrand((strand[k].Data())[0]);
      ann->SetStart(start[k]);
      ann->SetStop(stop[k]);
      ann->SetCrossHybType(crosshtype[k]);
      ann->SetProbesetType(psettype[k]);
      anntree->Fill();

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx << "> records imported...\r" << flush;
      }//if
   }//for_i

// Fill transcript annotation tree with sorted "rescue->FLmRNA->unmapped" annotation data
   for (Int_t i=0; i<size; i++) {
      Int_t k = index[i];
      if (!((psettype[k] == eINTRON) ||
            (psettype[k] == eEXON)   ||
            (psettype[k] == eUNMAPPED))) continue;

      // geneid as string
      str.Form("%d", geneid[k]);

      // fill annotation tree with data
      ann->SetUnitID(idx++);
      ann->SetTranscriptID(str);
      ann->SetName(genename[k]);
      ann->SetSymbol(genesymbol[k]);
      ann->SetAccession(geneaccess[k]);
      ann->SetEntrezID(entrezid[k]);
      ann->SetChromosome(seqname[k]);
      ann->SetCytoBand(cytoband[k]);
      ann->SetStrand((strand[k].Data())[0]);
      ann->SetStart(start[k]);
      ann->SetStop(stop[k]);
      ann->SetCrossHybType(crosshtype[k]);
      ann->SetProbesetType(psettype[k]);
      anntree->Fill();

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx << "> records imported...\r" << flush;
      }//if
   }//for_i


   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records imported...Finished" << endl;
   }//if

//TO DO:
// Add tree info to tree
//   AddAnnotationTreeInfo(anntree, fAnnTreeName);

// Write annotation tree to file
   if ((err = WriteTree(anntree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(anntree->GetName(), 0);
   }//if
   //delete tree from memory
   anntree->Delete("");
   anntree = 0;

// Clean up
cleanup:
   SafeDelete(ann);

   if (value)      {delete [] value;      value      = 0;}
   if (index)      {delete [] index;      index      = 0;}
   if (psettype)   {delete [] psettype;   psettype   = 0;}
   if (mrnaaccess) {delete [] mrnaaccess; mrnaaccess = 0;}
   if (entrezid)   {delete [] entrezid;   entrezid   = 0;}
   if (cytoband)   {delete [] cytoband;   cytoband   = 0;}
   if (genesymbol) {delete [] genesymbol; genesymbol = 0;}
   if (genename)   {delete [] genename;   genename   = 0;}
   if (geneaccess) {delete [] geneaccess; geneaccess = 0;}
   if (names)      {delete [] names;      names      = 0;}
   if (crosshtype) {delete [] crosshtype; crosshtype = 0;}
   if (stop)       {delete [] stop;       stop       = 0;}
   if (start)      {delete [] start;      start      = 0;}
   if (strand)     {delete [] strand;     strand     = 0;}
   if (seqname)    {delete [] seqname;    seqname    = 0;}
   if (geneid)     {delete [] geneid;     geneid     = 0;}
   if (hasColumn)  {delete [] hasColumn;  hasColumn  = 0;}

   return err;
}//ImportTransAnnotation

//______________________________________________________________________________
Short_t XGenomeChip::ProbesetType(const char *type)
{
   // Convert probeset type to id EProbesetType
   if(kCSa) cout << "------XGenomeChip::ProbesetType------" << endl;

   if (strstr(type, "exon:main") != 0) {
      return eMAINEXON;
   } else if (strstr(type, "main") != 0) {
//?   if (strstr(type, "main") != 0) {
      return eMAIN;
   } else if (strstr(type, "control->affx") != 0) {
      return eCONTROLAFFX;
   } else if (strstr(type, "control->chip") != 0) {
      return eCONTROLCHIP;
   } else if (strstr(type, "bgp->antigenomic") != 0) {
      return eANTIGENOMIC;
   } else if (strstr(type, "bgp->genomic") != 0) {
      return eGENOMIC;
   } else if (strstr(type, "intron") != 0) {
      return eINTRON;
   } else if (strstr(type, "exon") != 0) {
      return eEXON;
   } else if (strstr(type, "unmapped") != 0) {
      return eUNMAPPED;
   }//if

   return eUNKNOWNTYPE;
}//ProbesetType

//______________________________________________________________________________
TString XGenomeChip::ProbesetTypeID2Name(Short_t id)
{
   // Convert level_id to probeset level
   if(kCSa) cout << "------XGenomeChip::ProbesetTypeID2Name------" << endl;

   switch (id) {
      case eMAIN:
         return TString("main");
         break;

      case eMAINEXON:
         return TString("normgene->exon:main");
         break;

      case eCONTROLAFFX:
         return TString("control->affx");
         break;

      case eCONTROLCHIP:
         return TString("control->chip");
         break;

      case eANTIGENOMIC:
         return TString("control->bgp->antigenomic");
         break;

      case eGENOMIC:
         return TString("control->bgp->genomic");
         break;

      case eINTRON:
         return TString("normgene->intron");
         break;

      case eEXON:
         return TString("normgene->exon");
         break;

      case eUNMAPPED:
         return TString("rescue->FLmRNA->unmapped");
         break;

      default: //NA
         return TString("NA");
         break;
   }//switch

   return TString("NA");
}//ProbesetTypeID2Name

//______________________________________________________________________________
Int_t XGenomeChip::SchemeMask(Int_t type, Int_t subtype)
{
   // Convert probeset type "type" with optional probeset level "subtype" to mask:
   // fMask =  1024: type=pm:st, type=main, subtype=core
   // fMask =  512:  type=pm:st, type=main, subtype=extended
   // fMask =  256:  type=pm:st, type=main, subtype=full
   // fMask =  128:  type=pm:st, type=main, subtype=ambiguous
   // fMask =  64:   type=pm:st, type=main, subtype=free
   // fMask =  -64:  type=pm:st, type=main, subtype=unknown level
   // fMask =  32:   type=pm:st, type=control->affx
   // fMask =  8:    type=mm:st, type=control->affx
   // fMask =  16:   type=pm:at, type=control->affx
   // fMask =  4:    type=mm:at, type=control->affx
   // fMask = -1:    type=mm:st, type=control->bgp->genomic
   // fMask = -2:    type=pm:st, type=control->bgp->antigenomic ????
   // fMask = -4:    type=pm:st, type=normgene->intron
   // fMask = -8:    type=pm:st, type=normgene->exon (w/o main)
   // fMask = -16:   type=pm:st, type=rescue->FLmRNA->unmapped
   // fMask = -100:  type=pm:at, type=control->chip, subtype=jumbo-checkerboard:at
   // fMask = -101:  type=pm:at, type=control->chip, subtype=jumbo-checkerboard:st
   // fMask = -102:  type=pm:at, type=control->chip, subtype=thermo:at
   // fMask = -103:  type=pm:st, type=control->chip, subtype=thermo:st
   // fMask = -104:  type=pm:at, type=control->chip, subtype=trigrid:at
   // fMask = -105:  type=pm:st, type=control->chip, subtype=trigrid:st
   // fMask = -106:  type=pm:at, type=control->chip, subtype=generic:at
   // fMask = -107:  type=pm:st, type=control->chip, subtype=generic:st
   // fMask = -108:  type=?????, type=control->chip, subtype=blank
   // fMask = -109:  type=pm:st, type=control->chip, subtype=pm:st
   if(kCSa) cout << "------XGenomeChip::SchemeMask------" << endl;

   if (type == eMAIN) {
      if (subtype == ePARACORE) {
         return ePARACORE;
      } else if (subtype == ePARAEXTENDED) {
         return ePARAEXTENDED;
      } else if (subtype == ePARAFULL) {
         return ePARAFULL;
      } else if (subtype == eAMBIGUOUS) {
         return eAMBIGUOUS;
      } else if (subtype == eFREE) {
         return eFREE;
      } else if (subtype == eNOLEVEL) {
         return eNOLEVEL;
      }//if
   } else if (type == eCONTROLAFFX) {
      if (subtype == ePMST) {
         return ePMST;
      } else if (subtype == ePMAT) {
         return ePMAT;
      } else if (subtype == eMMST) {
         return eMMST;
      } else if (subtype == eMMAT) {
         return eMMAT;
      }//if
   } else if (type == eCONTROLCHIP) {
      if (subtype == eJUMBOAT) {
         return eJUMBOAT;
      } else if (subtype == eJUMBOST) {
         return eJUMBOST;
      } else if (subtype == eTHERMOAT) {
         return eTHERMOAT;
      } else if (subtype == eTHERMOST) {
         return eTHERMOST;
      } else if (subtype == eTRIGRIDAT) {
         return eTRIGRIDAT;
      } else if (subtype == eTRIGRIDST) {
         return eTRIGRIDST;
      } else if (subtype == eGENERICAT) {
         return eGENERICAT;
      } else if (subtype == eGENERICST) {
         return eGENERICST;
      } else if (subtype == eBLANK) {
         return eBLANK;
      } else if (subtype == eCTRLPMST) {
         return eCTRLPMST;
      } else if (subtype == eUNKNOWN) {
         return eUNKNOWN;
      }//if
   } else if (type == eANTIGENOMIC) {
      return eANTIGENOMIC;
   } else if (type == eGENOMIC) {
      return eGENOMIC;
   } else if (type == eINTRON) {
      if (subtype == 0) {
         return eINTRON;
      } else if (subtype == ePARACORE) {
         return ePARACORE;
      } else if (subtype == ePARAEXTENDED) {
         return ePARAEXTENDED;
      } else if (subtype == ePARAFULL) {
         return ePARAFULL;
      } else if (subtype == eAMBIGUOUS) {
         return eAMBIGUOUS;
      } else if (subtype == eFREE) {
         return eFREE;
      }//if
      return eINTRON;
   } else if (type == eEXON) {
      if (subtype == 0) {
         return eEXON;
      } else if (subtype == ePARACORE) {
         return ePARACORE;
      } else if (subtype == ePARAEXTENDED) {
         return ePARAEXTENDED;
      } else if (subtype == ePARAFULL) {
         return ePARAFULL;
      } else if (subtype == eAMBIGUOUS) {
         return eAMBIGUOUS;
      } else if (subtype == eFREE) {
         return eFREE;
      }//if
      return eEXON;
   } else if (type == eUNMAPPED) {
      return eUNMAPPED;
   }//if

   return eUNKNOWNTYPE;
}//SchemeMask

//______________________________________________________________________________
Int_t XGenomeChip::SchemeMask(Short_t xtype)
{
   // Convert cross-hybridization probeset_type "xtype" to mask
   // Note: since e.g. HuGene array contains only well annotated genes from
   // full-length RefSeq, Ensembl etc, we can use eMETACORE and ePARACORE
   if(kCSa) cout << "------XGenomeChip::SchemeMask(main)------" << endl;

   if (xtype == 1) { //(crosshyb_type==unique)
      return eMETACORE;
   }//if

   return ePARACORE;
}//SchemeMask


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XExonChip                                                            //
//                                                                      //
// Class for Affymetrix Exon arrays                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
XExonChip::XExonChip()
          :XGenomeChip()
{
   // Default ExonArray constructor
   if(kCS) cout << "---XExonChip::XExonChip(default)------" << endl;

   fExnTreeName = "";
   fPbsTreeName = "";
   fAnxTreeName = "";
   fNExonUnits  = 0;
   fNExons      = 0;
}//Constructor

//______________________________________________________________________________
XExonChip::XExonChip(const char *name, const char *title)
          :XGenomeChip(name, title)
{
   // Normal ExonArray constructor
   if(kCS) cout << "---XExonChip::XExonChip------" << endl;

   fExnTreeName = "NA";
   fPbsTreeName = "NA";
   fAnxTreeName = "NA";
   fNExonUnits  = 0;
   fNExons      = 0;
}//Constructor

//______________________________________________________________________________
XExonChip::XExonChip(const XExonChip &chip) 
          :XGenomeChip(chip)
{
   // ExonArray copy constructor
   if(kCS) cout << "---XExonChip::XExonChip(copy)------" << endl;

   fExnTreeName = chip.fExnTreeName.Copy();
   fPbsTreeName = chip.fPbsTreeName.Copy();
   fAnxTreeName = chip.fAnxTreeName.Copy();
   fNExonUnits  = chip.fNExonUnits;
   fNExons      = chip.fNExons;
//TO DO   fAffxNames   = chip.fAffxNames;
}//CopyConstructor

//______________________________________________________________________________
XExonChip& XExonChip::operator=(const XExonChip& rhs)
{
   // ExonArray assignment operator.

   if (this != &rhs) {
      XExonChip::operator=(rhs);

      fExnTreeName = rhs.fExnTreeName.Copy();
      fPbsTreeName = rhs.fPbsTreeName.Copy();
      fAnxTreeName = rhs.fAnxTreeName.Copy();
      fNExonUnits  = rhs.fNExonUnits;
      fNExons      = rhs.fNExons;
//TO DO      fAffxNames   = chip.fAffxNames;
   }//if

   return *this;
}//operator=

//______________________________________________________________________________
XExonChip::~XExonChip()
{
   // ExonArray destructor
   if(kCS) cout << "---XExonChip::~XExonChip------" << endl;

}//Destructor

//______________________________________________________________________________
Int_t XExonChip::ExportSchemeTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in scheme tree to file output
   // varlist can be at most: "fProbeLen:fMask:fExonID" or "*"
   if(kCS) cout << "------XExonChip::ExportSchemeTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasPLen = kFALSE;
   Bool_t hasMask = kFALSE;
   Bool_t hasExon = kFALSE;
   Bool_t hasPSet = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasPLen = kTRUE;
      hasMask = kTRUE;
      hasExon = kTRUE;
      hasPSet = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fProbeLen")   == 0) {hasPLen = kTRUE;}
         if (strcmp(name,"fMask")       == 0) {hasMask = kTRUE;}
         if (strcmp(name,"fExonID")     == 0) {hasExon = kTRUE;}
         if (strcmp(name,"fProbesetID") == 0) {hasPSet = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get scheme tree for this chip
   XExonScheme *scheme = 0;
   fTree = (TTree*)gDirectory->Get((names[0]).Data());
   if (!fTree) return errGetTree;
   fTree->SetBranchAddress("ScmBranch", &scheme);

   Int_t nentries = (Int_t)(fTree->GetEntries());
//??   Int_t size     = fNRows*fNCols;
   Int_t size     = fNProbes;
   if (nentries != size) {
      TString str = ""; str += size;
//??      return fManager->HandleError(errNumTreeEntries, fTree->GetName(), str);
   }//if

// Output header
   output << "UNIT_ID" << sep << "X" << sep << "Y";
   if (hasPLen) output << sep << "ProbeLength";
   if (hasMask) output << sep << "Mask";
   if (hasExon) output << sep << "EXON_ID";
   if (hasPSet) output << sep << "PROBESET_ID";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << scheme->GetUnitID() << sep
             << scheme->GetX() << sep << scheme->GetY();
      if (hasPLen) output << sep << scheme->GetProbeLength();
      if (hasMask) output << sep << scheme->GetMask();
      if (hasExon) output << sep << scheme->GetExonID();
      if (hasPSet) output << sep << scheme->GetProbesetID();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportSchemeTree

//______________________________________________________________________________
Int_t XExonChip::ExportUnitTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in unit tree to file output
   if(kCS) cout << "------XExonChip::ExportUnitTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasName = kFALSE;
   Bool_t hasUnit = kFALSE;
   Bool_t hasCell = kFALSE;
   Bool_t hasAtom = kFALSE;
   Bool_t hasSubu = kFALSE;
   Bool_t hasType = kFALSE;

   if (strcmp(varlist,"*")  == 0) {
      hasName = kTRUE;
      hasUnit = kTRUE;
      hasCell = kTRUE;
      hasAtom = kTRUE;
      hasSubu = kTRUE;
      hasType = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fUnitName")  == 0) {hasName = kTRUE;}
         if (strcmp(name,"fSubUnitID") == 0) {hasUnit = kTRUE;}
         if (strcmp(name,"fNCells")    == 0) {hasCell = kTRUE;}
         if (strcmp(name,"fNAtoms")    == 0) {hasAtom = kTRUE;}
         if (strcmp(name,"fNSubunits") == 0) {hasSubu = kTRUE;}
         if (strcmp(name,"fUnitType")  == 0) {hasType = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get unit tree for this chip
   XExonUnit *unit = 0;
   fTree = (TTree*)gDirectory->Get((names[0]).Data());
   if (fTree == 0) return errGetTree;
   fTree->SetBranchAddress("IdxBranch", &unit);

   TString exten  = Path2Name(fTree->GetName(),".",";");
   TString substr = "SubunitID";

   Int_t nentries = (Int_t)(fTree->GetEntries());
   Int_t numunits = 0;
   if (strcmp(exten.Data(),kExtenScheme[2]) == 0) {
      numunits = fNUnits;
      substr   = "TranscriptClusterID";
   } else if (strcmp(exten.Data(),kExtenScheme[4]) == 0) {
      // Branch UnitName is not filled for exon tree
      if (hasName) {
         hasName = kFALSE;
         hasUnit = kTRUE;
      }//if

      numunits = fNExonUnits;
      substr   = "ExonID";
   } else if (strcmp(exten.Data(),kExtenScheme[5]) == 0) {
      // Branch UnitName is not filled for probeset tree
      if (hasName) {
         hasName = kFALSE;
         hasUnit = kTRUE;
      }//if

      numunits = fNProbesets;
      substr   = "ProbesetID";
   }//if

   if (nentries != numunits) {
      TString str = ""; str += numunits;
//??      return fManager->HandleError(errNumTreeEntries, fTree->GetName(), str);
   }//if

// Output header
   output << "UNIT_ID";
   if (hasName) output << sep << "UnitName";
   if (hasUnit) output << sep << substr.Data();
   if (hasCell) output << sep << "NumCells";
   if (hasAtom) output << sep << "NumAtoms";
   if (hasSubu) output << sep << "NumSubunits";
   if (hasType) output << sep << "UnitType";
   output << endl;

// Export selected variables
   for (Int_t i=0; i<nentries; i++) {
      fTree->GetEntry(i);
      output << unit->GetUnitID();
      if (hasName) output << sep << unit->GetUnitName();
      if (hasUnit) output << sep << unit->GetSubUnitID();
      if (hasCell) output << sep << unit->GetNumCells();
      if (hasAtom) output << sep << unit->GetNumAtoms();
      if (hasSubu) output << sep << unit->GetNumSubunits();
      if (hasType) output << sep << unit->GetUnitType();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportUnitTree

//______________________________________________________________________________
Int_t XExonChip::ExportExonAnnotTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in exon annotation tree to file output
   if(kCS) cout << "------XExonChip::ExportExonAnnotTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasTran = kFALSE;  //transcript_cluster_id
   Bool_t hasName = kFALSE;  //gene accession 
   Bool_t hasSymb = kFALSE;  //gene symbol
   Bool_t hasAccn = kFALSE;  //mRNA accession
//   Bool_t hasEntr = kFALSE;  //entrez ID
   Bool_t hasChro = kFALSE;  //chromosome
//   Bool_t hasCyto = kFALSE;  //cytoband
   Bool_t hasStar = kFALSE;  //start position
   Bool_t hasStop = kFALSE;  //stop position
   Bool_t hasStrd = kFALSE;  //strand
   Bool_t hasExon = kFALSE;  //exon_id

   if (strcmp(varlist,"*")  == 0) {
      hasTran = kTRUE;
      hasName = kTRUE;
      hasSymb = kTRUE;
      hasAccn = kTRUE;
//      hasEntr = kTRUE;
      hasChro = kTRUE;
//      hasCyto = kTRUE;
      hasStar = kTRUE;
      hasStop = kTRUE;
      hasStrd = kTRUE;
      hasExon = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fTranscriptID") == 0) {hasTran = kTRUE;}
         if (strcmp(name,"fName")         == 0) {hasName = kTRUE;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymb = kTRUE;}
         if (strcmp(name,"fAccession")    == 0) {hasAccn = kTRUE;}
//         if (strcmp(name,"fEntrezID")     == 0) {hasEntr = kTRUE;}
         if (strcmp(name,"fChromosome")   == 0) {hasChro = kTRUE;}
//         if (strcmp(name,"fCytoBand")     == 0) {hasCyto = kTRUE;}
         if (strcmp(name,"fStart")        == 0) {hasStar = kTRUE;}
         if (strcmp(name,"fStop")         == 0) {hasStop = kTRUE;}
         if (strcmp(name,"fStrand")       == 0) {hasStrd = kTRUE;}
         if (strcmp(name,"fExonID")       == 0) {hasExon = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get exon annotation tree
   XExonAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get((names[0].Data())); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

// Output header
   output << "UNIT_ID";
   if (hasExon) output << sep << "ExonID";
   if (hasTran) output << sep << "TranscriptClusterID";
   if (hasName) output << sep << "GeneAccession";
   if (hasSymb) output << sep << "GeneSymbol";
   if (hasAccn) output << sep << "mRNAAccession";
//   if (hasEntr) output << sep << "EntrezID";
   if (hasChro) output << sep << "Chromosome";
//   if (hasCyto) output << sep << "Cytoband";
   if (hasStar) output << sep << "Start";
   if (hasStop) output << sep << "Stop";
   if (hasStrd) output << sep << "Strand";
   output << endl;

// Export selected variables
   Int_t nentries = (Int_t)(anntree->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      anntree->GetEntry(i);
      output << annot->GetUnitID();
      if (hasExon) output << sep << annot->GetExonID();
      if (hasTran) output << sep << annot->GetTranscriptID();
      if (hasName) output << sep << annot->GetName();
      if (hasSymb) output << sep << annot->GetSymbol();
      if (hasAccn) output << sep << annot->GetAccession();
//      if (hasEntr) output << sep << annot->GetEntrezID();
      if (hasChro) output << sep << annot->GetChromosome();
//      if (hasCyto) output << sep << annot->GetCytoBand();
      if (hasStar) output << sep << annot->GetStart();
      if (hasStop) output << sep << annot->GetStop();
      if (hasStrd) output << sep << annot->GetStrand();
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportExonAnnotTree

//______________________________________________________________________________
Int_t XExonChip::ExportProbesetAnnotTree(Int_t n, TString *names, const char *varlist,
                 ofstream &output, const char *sep)
{
   // Export data stored in probeset annotation tree to file output
   if(kCS) cout << "------XExonChip::ExportProbesetAnnotTree------" << endl;

// Crosscheck
   if (n != 1) return errGeneral;

// Decompose varlist
   Bool_t hasTran = kFALSE;  //transcript_cluster_id
   Bool_t hasName = kFALSE;  //gene accession 
   Bool_t hasSymb = kFALSE;  //gene symbol
   Bool_t hasAccn = kFALSE;  //mRNA accession
   Bool_t hasChro = kFALSE;  //chromosome
//   Bool_t hasCyto = kFALSE;  //cytoband
   Bool_t hasStar = kFALSE;  //start position
   Bool_t hasStop = kFALSE;  //stop position
   Bool_t hasStrd = kFALSE;  //strand
   Bool_t hasExon = kFALSE;  //exon_id
   Bool_t hasPSet = kFALSE;  //probeset_id
   Bool_t hasNPrb = kFALSE;  //probe_count
   Bool_t hasXTyp = kFALSE;  //crosshyb_type
   Bool_t hasLevl = kFALSE;  //level, converted to id
   Bool_t hasBoun = kFALSE;  //bounded, noBoundedEvidence
   Bool_t hasFull = kFALSE;  //flag for putative full-length mRNA
   Bool_t hasPTyp = kFALSE;  //probeset_type

   if (strcmp(varlist,"*")  == 0) {
      hasTran = kTRUE;
      hasName = kTRUE;
      hasSymb = kTRUE;
      hasAccn = kTRUE;
      hasChro = kTRUE;
//      hasCyto = kTRUE;
      hasStar = kTRUE;
      hasStop = kTRUE;
      hasStrd = kTRUE;
      hasExon = kTRUE;
      hasPSet = kTRUE;
      hasNPrb = kTRUE;
      hasXTyp = kTRUE;
      hasLevl = kTRUE;
      hasBoun = kTRUE;
      hasFull = kTRUE;
      hasPTyp = kTRUE;
   } else {
      char *name  = new char[strlen(varlist) + 1];
      char *dname = name;
      name = strtok(strcpy(name,varlist),":");
      while(name) {
         if (strcmp(name,"fTranscriptID") == 0) {hasTran = kTRUE;}
         if (strcmp(name,"fName")         == 0) {hasName = kTRUE;}
         if (strcmp(name,"fSymbol")       == 0) {hasSymb = kTRUE;}
         if (strcmp(name,"fAccession")    == 0) {hasAccn = kTRUE;}
         if (strcmp(name,"fChromosome")   == 0) {hasChro = kTRUE;}
//         if (strcmp(name,"fCytoBand")     == 0) {hasCyto = kTRUE;}
         if (strcmp(name,"fStart")        == 0) {hasStar = kTRUE;}
         if (strcmp(name,"fStop")         == 0) {hasStop = kTRUE;}
         if (strcmp(name,"fStrand")       == 0) {hasStrd = kTRUE;}
         if (strcmp(name,"fExonID")       == 0) {hasExon = kTRUE;}
         if (strcmp(name,"fProbesetID")   == 0) {hasPSet = kTRUE;}
         if (strcmp(name,"fNProbes")      == 0) {hasNPrb = kTRUE;}
         if (strcmp(name,"fCrossHybType") == 0) {hasXTyp = kTRUE;}
         if (strcmp(name,"fLevelID")      == 0) {hasLevl = kTRUE;}
         if (strcmp(name,"fBounded")      == 0) {hasBoun = kTRUE;}
         if (strcmp(name,"fFull")         == 0) {hasFull = kTRUE;}
         if (strcmp(name,"fProbesetType") == 0) {hasPTyp = kTRUE;}
         name = strtok(NULL, ":");
         if (name == 0) break;
      }//while
      delete [] dname;
   }//if

// Get probeset annotation tree
   XProbesetAnnotation *annot = 0;
   TTree *anntree = (TTree*)gDirectory->Get((names[0].Data())); 
   if (anntree == 0) return errGetTree;
   anntree->SetBranchAddress("AnnBranch", &annot);

// Output header
   output << "UNIT_ID";
   if (hasPSet) output << sep << "ProbesetID";
   if (hasExon) output << sep << "ExonID";
   if (hasTran) output << sep << "TranscriptClusterID";
   if (hasName) output << sep << "GeneAccession";
   if (hasSymb) output << sep << "GeneSymbol";
   if (hasAccn) output << sep << "mRNAAccession";
   if (hasChro) output << sep << "Chromosome";
//   if (hasCyto) output << sep << "Cytoband";
   if (hasStar) output << sep << "Start";
   if (hasStop) output << sep << "Stop";
   if (hasStrd) output << sep << "Strand";
   if (hasNPrb) output << sep << "ProbeCount";
   if (hasXTyp) output << sep << "CrossHybType";
   if (hasLevl) output << sep << "Level"; //convert ID to level
   if (hasBoun) output << sep << "Bounded_NoBoundedEvidence";
   if (hasFull) output << sep << "FullLength";
   if (hasPTyp) output << sep << "ProbesetType";
   output << endl;

// Export selected variables
   Int_t nentries = (Int_t)(anntree->GetEntries());
   for (Int_t i=0; i<nentries; i++) {
      anntree->GetEntry(i);
      output << annot->GetUnitID();
      if (hasPSet) output << sep << annot->GetProbesetID();
      if (hasExon) output << sep << annot->GetExonID();
      if (hasTran) output << sep << annot->GetTranscriptID();
      if (hasName) output << sep << annot->GetName();
      if (hasSymb) output << sep << annot->GetSymbol();
      if (hasAccn) output << sep << annot->GetAccession();
      if (hasChro) output << sep << annot->GetChromosome();
//      if (hasCyto) output << sep << annot->GetCytoBand();
      if (hasStar) output << sep << annot->GetStart();
      if (hasStop) output << sep << annot->GetStop();
      if (hasStrd) output << sep << annot->GetStrand();
      if (hasNPrb) output << sep << annot->GetNumProbes();
      if (hasXTyp) output << sep << annot->GetCrossHybType();
//TO DO      if (hasXTyp) output << sep << CrossHybID2Name(annot->GetCrossHybType());
      if (hasLevl) output << sep << LevelID2Level(annot->GetLevelID());
      if (hasBoun) output << sep << annot->GetBounded();
      if (hasFull) output << sep << annot->GetFullFlag();
      if (hasPTyp) output << sep << ProbesetTypeID2Name(annot->GetProbesetType());
      output << endl;

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "<" << i+1 << "> records exported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "<" << nentries << "> records exported...Finished" << endl;
   }//if

   return errNoErr;
}//ExportProbesetAnnotTree

//______________________________________________________________________________
Int_t XExonChip::ReadData(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import probe information from infile and store as probe tree
   if(kCS) cout << "------XExonChip::ReadData------" << endl;

   char nextline[kBufSize];
   Int_t err = errNoErr;

//TO DO
//check: if (fNCols == 0) error: 
//better: if (tree.cxy == 0) error: missing layout tree!!

// Get probeset annotation tree
   TString treename = TString(fName) + "." + kExtenAnnot[2];
   XProbesetAnnotation *psann = 0;
   TTree *psanntree = (TTree*)gDirectory->Get(treename.Data());
   if (psanntree == 0) {
      cerr << "Error: Missing probeset annotation tree <" << treename.Data() << ">."
           << endl;
      return errGetTree;
   }//if
   psanntree->SetBranchAddress("AnnBranch", &psann);

   Int_t nentries = psanntree->GetEntries();

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Create new scheme tree
   fScmTreeName      = TString(fName) + "." + exten;
   TTree *schemetree = new TTree(fScmTreeName, "scheme information");
   if (schemetree == 0) return errCreateTree;
   XExonScheme *scheme = 0;
   scheme = new XExonScheme();
   schemetree->Branch("ScmBranch", "XExonScheme", &scheme, 64000, split);

// Create new unit tree
   fIdxTreeName    = TString(fName) + "." + TString(kExtenScheme[2]);
   TTree *unittree = new TTree(fIdxTreeName, "unit information");
   if (unittree == 0) return errCreateTree;
   XExonUnit *unit = 0;
   unit = new XExonUnit();
   unittree->Branch("IdxBranch", "XExonUnit", &unit, 64000, split);

// Create new exon unit tree
   fExnTreeName    = TString(fName) + "." + TString(kExtenScheme[4]);
   TTree *exontree = new TTree(fExnTreeName, "exon unit information");
   if (exontree == 0) return errCreateTree;
   XExonUnit *exon = 0;
   exon = new XExonUnit();
   exontree->Branch("IdxBranch", "XExonUnit", &exon, 64000, split);

// Create new probeset unit tree
   fPbsTreeName    = TString(fName) + "." + TString(kExtenScheme[5]);
   TTree *psettree = new TTree(fPbsTreeName, "probeset unit information");
   if (psettree == 0) return errCreateTree;
   XExonUnit *pset = 0;
   pset = new XExonUnit();
   psettree->Branch("IdxBranch", "XExonUnit", &pset, 64000, split);

// Create new probe tree
   fPrbTreeName     = TString(fName) + "." + TString(kExtenScheme[3]);
   TTree *probetree = new TTree(fPrbTreeName, "probe info for scheme");
   if (probetree == 0) return errCreateTree;
   XGCProbe *probe = 0;
   probe = new XGCProbe();
   probetree->Branch("PrbBranch", "XGCProbe", &probe, 64000, split);

   XBitSet exonbitmsk, transbitmsk;
   Int_t   exonlevel  = eNOLEVEL;
   Int_t   translevel = eNOLEVEL;

   TString order;
   TString str;
   Int_t   size;  //probeset_count
   Int_t   begin, end;
   Int_t   probeset_id;
   TString probeset_type;
   TString geneaccess, genesymbol, mrnaaccess, seqname;

   Int_t   idx;
   Int_t   ctrlsize, numctrl;
   Bool_t  flag;

   Int_t   probe_id, gc_count, probe_len, inter_pos;
   TString probe_type, sequence;
   Int_t   x, y;

   Int_t    unitID   = 0;
   Int_t    mask     = 0;
   Int_t    mainexon = 0;
   Int_t    numcells, numatoms, numprobs;
   Double_t Tm;
   Short_t  strd;
   Short_t  unittype;
   Int_t    probesetype = eUNKNOWNTYPE;
   Int_t    geneID, exonID, psetID;

// Init counting variables
   Int_t    scmcount = 0;       //count number of scheme data (<=fNRows*fNCols)
   Int_t    idxcount = 0;       //count number of units
   Int_t    exncount = 0;       //count number of exons
   Int_t    ctlcount = 0;       //count number of negative controls (EProbesetType<0)
//?   Int_t    numpsets = -1;      //count number of probesets per transcript_cluster
   Int_t    numexons = -1;      //count number of exons per transcript_cluster
   Int_t    numcelex = -1;      //number of cells per exon
   Int_t    numperex = -1;      //number of probesets per exon
   Int_t    maxunits = 0;       //maximum number of unit cells
   Int_t    maxexons = 0;       //maximum number of exon cells
   Int_t    maxpsets = 0;       //maximum number of probeset cells
   Int_t    minunits = 9999999; //minimum number of unit cells
   Int_t    minexons = 9999999; //minimum number of exon cells
   Int_t    minpsets = 9999999; //minimum number of probeset cells

   Int_t    intexonID = -1;     //internal exon ID
   Int_t    intpsetID = -1;     //internal probeset ID
   Bool_t   didWhile  = kFALSE; //to avoid circular for-loop

// Init local arrays
   Long_t  *pos = 0;  //positions of probeset_ids
   Int_t   *psi = 0;  //probeset_id
   Int_t   *pst = 0;  //probeset_type as EProbesetType
   Short_t *msk = 0;  //mask: set msk=1 for probeset_ids already read
   Long_t  *arr = 0;  //temporary array to store other arrays

// Init local arrays for sorting
   Int_t    *index = 0;  //sort index
   Long64_t *value = 0;  //value array to store the two columns to sort
   Long64_t majorv, minorv;  //values to store in array value

// Init local arrays for control->chip
   Long_t   *ctrlpos = 0;  //positions of probe_ids for control->chip
   Int_t    *ctrlpbt = 0;  //probe_type as EControlchipType
   Int_t    *ctrlidx = 0;  //sort index
   Long64_t *ctrlval = 0;  //value array to store the two columns to sort

// Init local arrays for sorted data from probeset annotation tree
   Int_t   *arrPSet  = 0;  //array for probeset_id
   Int_t   *arrCount = 0;  //array for probe_count
   Int_t   *arrGene  = 0;  //array for transcript_cluster_id
   Int_t   *arrExon  = 0;  //array for exon_id
   Short_t *arrXType = 0;  //array for crosshyb_type (unique, similar, mixed)
   Short_t *arrLevel = 0;  //array for level id
   Short_t *arrBound = 0;  //array for bounded

// Get position of first probeset_id after header lines
   streampos position = input.tellg();

//--------------------------------------------------------------------
// Read data and sort for probeset_type and position of probeset_id
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Reading data from input file..." << endl;
   }//if

// Read data to get number of probesets (header0)
   size = 0;
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) break;
      if (!input.good()) {err = errPrematureEOF; break;}
      str = RemoveEnds(&nextline[0], begin, end);

      if (begin==0) size++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   Number of probesets is <" << size << ">." << endl;
   }//if

// Check if number of annotated probesets is equal to probeset_count
   if (fNProbesets != size) {
      if (XManager::fgVerbose) {
         cout << "   Note: Number of annotated probesets <" << fNProbesets  
              << "> is not equal to number of probesets <" << size << ">." << endl;
      }//if
      // set final number of probesets
      fNProbesets = size;
   }//if

// Initialize memory for local arrays
   if (!(pos = new (nothrow) Long_t[size]))  {err = errInitMemory; goto cleanup;}
   if (!(psi = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(pst = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(msk = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize memory for sorting arrays
   if (!(index = new (nothrow) Int_t[size]))    {err = errInitMemory; goto cleanup;}
   if (!(value = new (nothrow) Long64_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize array pst to eINIT since eCONTROLCHIP=0!
//not necessary???   for (Int_t i=0; i<size; i++) pst[i] = eINIT;

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Read data to get probeset_id and type (header0)
   idx = 0;
   while (input.good()) {
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) break;
      if (!input.good()) {err = errPrematureEOF; break;}
      str = RemoveEnds(&nextline[0], begin, end);

      // get position, and probeset_id and type (header0)
      if (begin==0) {
         pos[idx] = (Long_t)input.tellg() - (Long_t)(str.Length() + end + 1);
         psi[idx] = atoi(strtok((char*)str.Data(), sep));
         pst[idx] = ProbesetType(strtok(0, sep));
         msk[idx] = 0;

         // fill sort values
         majorv = (Long64_t)pst[idx];  //sort for probeset_type first
         minorv = (Long64_t)pos[idx];  //then sort for position of probeset_id
         value[idx]  = majorv << 31;
         value[idx] += minorv;

         if (XManager::fgVerbose && idx%10000 == 0) {
            cout << "   <" << idx + 1 << "> records read...\r" << flush;
         }//if
         idx++;
      }//if
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records read...Finished" << endl;
   }//if

   if (fNProbesets != idx) {
      cerr << "Error: Number of probeset_ids read <" << idx  
           << "> is not equal to number of probesets <" << fNProbesets << ">" << endl;
      err = errReadingInput;
      goto cleanup;
   }//if

// Sort for probeset_type first and then for position of probeset_id
   if (XManager::fgVerbose) {
      cout << "   Sorting data for probeset_type and position..." << endl;
   }//if
   TMath::Sort(size, value, index, 0);

   // delete value here to free memory
   if (value) {delete [] value; value = 0;}

// Replace data in arrays with sorted data using temporary array arr
   // create temporary array (one temporary array only to save memory)
   if (!(arr = new (nothrow) Long_t[size])) {err = errInitMemory; goto cleanup;}

   // replace probeset_id positions with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pos[index[i]];
   for (Int_t i=0; i<size; i++) pos[i] = arr[i];

   // replace probeset_ids with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = psi[index[i]];
   for (Int_t i=0; i<size; i++) psi[i] = arr[i];

   // replace probeset_types with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pst[index[i]];
   for (Int_t i=0; i<size; i++) pst[i] = arr[i];

   // delete temporary array here to free memory
   if (arr) {delete [] arr; arr = 0;}

// count number of negative controls (EProbesetType<0)
   for (Int_t i=0; i<size; i++) {
      if (pst[i] == eGENOMIC || pst[i] == eANTIGENOMIC) ctlcount++;
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   Total number of controls is <" << ctlcount << ">" << endl;
   }//if


//--------------------------------------------------------------------
// 1. Fill tree(s) with data for probeset type: control->chip
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: control->chip..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get number of probes for control->chip
   ctrlsize = numctrl = 0;
   for (Int_t i=0; i<size; i++) {
      if (pst[i] != eCONTROLCHIP) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      begin = 1;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) break;
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) continue;   //skip data for header1
         if (begin == 2) ctrlsize++; //count data for header2 (probe_ids)
      }//while

      numctrl++;  //count number of probeset_ids
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   Number of control->chip items is <" << numctrl << ">." << endl;
   }//if

   if (!(ctrlpos = new (nothrow) Long_t[ctrlsize]))   {err = errInitMemory; goto cleanup;}
   if (!(ctrlpbt = new (nothrow) Int_t[ctrlsize]))    {err = errInitMemory; goto cleanup;}
   if (!(ctrlidx = new (nothrow) Int_t[ctrlsize]))    {err = errInitMemory; goto cleanup;}
   if (!(ctrlval = new (nothrow) Long64_t[ctrlsize])) {err = errInitMemory; goto cleanup;}

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for probe_types of control->chip
   idx = 0;
   for (Int_t i=0; i<size; i++) {
      if (pst[i] != eCONTROLCHIP) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin = 1;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) continue; //skip data for header1
         if (begin == 2) { //data for header2
            ctrlpos[idx] = (Long_t)input.tellg() - (Long_t)(str.Length() + end + 1);

            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(0, sep);
            ctrlpbt[idx] = ControlchipType(probe_type);

            // fill sort values
            majorv = (Long64_t)ctrlpbt[idx];  //sort for probe_type first
            minorv = (Long64_t)ctrlpos[idx];  //then sort for position of probe_id
            ctrlval[idx]  = majorv << 31;
            ctrlval[idx] += minorv;

            idx++;
         }//if
      }//while
   }//for_i

// Sort for control->chip probe_type first and then for position of probe_id
   TMath::Sort(ctrlsize, ctrlval, ctrlidx, 0);

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get data for probe_types of control->chip
   idx = -9999;
   flag = kFALSE;
   numcells = 0;
   numatoms = 1;
   for (Int_t i=0; i<=ctrlsize; i++) {
      if (ctrlsize == 0) break;
      if (i == ctrlsize) {

         // fill transcript unit tree
         str = TString(""); str += ctrlpbt[ctrlidx[i-1]];
         unit->SetUnitName(str);
//x         unit->SetUnitName(probe_type);
         unit->SetUnitID(ctrlpbt[ctrlidx[i-1]]);
         unit->SetSubUnitID(ctrlpbt[ctrlidx[i-1]]);
         unit->SetUnitType(ctrlpbt[ctrlidx[i-1]]);
         unit->SetNumCells(numcells);
         unit->SetNumAtoms(numatoms);
         unit->SetNumSubunits(numexons);
         unittree->Fill();

         // fill exon unit tree
         exon->SetUnitName(str);
         exon->SetUnitID(ctrlpbt[ctrlidx[i-1]]);
         exon->SetSubUnitID(ctrlpbt[ctrlidx[i-1]]);
         exon->SetUnitType(ctrlpbt[ctrlidx[i-1]]);
         exon->SetNumCells(numcells);
         exon->SetNumAtoms(numatoms);
         exon->SetNumSubunits(numexons);
         exontree->Fill();

//TO DO: > 1 probesets per ctrl!!!!
         // fill probeset unit tree
         pset->SetUnitName(str);
         pset->SetUnitID(ctrlpbt[ctrlidx[i-1]]);
         pset->SetSubUnitID(ctrlpbt[ctrlidx[i-1]]); //not probe_type
//to do: correct probeset_id => subdivide e.g. thermo.at also by probeset_id
         pset->SetUnitType(ctrlpbt[ctrlidx[i-1]]);
         pset->SetNumCells(numcells);
         pset->SetNumAtoms(numatoms);
         pset->SetNumSubunits(numexons);
         psettree->Fill();

      // minimal number of cells
         minunits = (numcells <  minunits) ? numcells : minunits;
         minexons = (numcells <  minexons) ? numcells : minexons;
         minpsets = (numcells <  minpsets) ? numcells : minpsets;

      // maximal number of cells
         maxunits = (numcells >  maxunits) ? numcells : maxunits;
         maxexons = (numcells >  maxexons) ? numcells : maxexons;
         maxpsets = (numcells >  maxpsets) ? numcells : maxpsets;

         idxcount++;
         break;
      }//if

      Int_t k = ctrlidx[i];
      if (ctrlpbt[k] != idx) {
         idx = ctrlpbt[k];
         if (flag) {

            // fill transcript unit tree
            str = TString(""); str += ctrlpbt[ctrlidx[i-1]];
            unit->SetUnitName(str);
//x            unit->SetUnitName(probe_type);
            unit->SetUnitID(ctrlpbt[ctrlidx[i-1]]);
            unit->SetSubUnitID(ctrlpbt[ctrlidx[i-1]]);
//to do: correct probeset_id => subdivide e.g. thermo.at also by probeset_id
            unit->SetUnitType(ctrlpbt[ctrlidx[i-1]]);
            unit->SetNumCells(numcells);
            unit->SetNumAtoms(numatoms);
            unit->SetNumSubunits(numexons);
            unittree->Fill();

            // fill exon unit tree
            exon->SetUnitName(str);
            exon->SetUnitID(ctrlpbt[ctrlidx[i-1]]);
            exon->SetSubUnitID(ctrlpbt[ctrlidx[i-1]]);
            exon->SetUnitType(ctrlpbt[ctrlidx[i-1]]);
            exon->SetNumCells(numcells);
            exon->SetNumAtoms(numatoms);
            exon->SetNumSubunits(numexons);
            exontree->Fill();

//TO DO: > 1 probesets per ctrl!!!!
            // fill probeset unit tree
            pset->SetUnitName(str);
            pset->SetUnitID(ctrlpbt[ctrlidx[i-1]]);
            pset->SetSubUnitID(ctrlpbt[ctrlidx[i-1]]); //not probe_type
            pset->SetUnitType(ctrlpbt[ctrlidx[i-1]]);
            pset->SetNumCells(numcells);
            pset->SetNumAtoms(numatoms);
            pset->SetNumSubunits(numexons);
            psettree->Fill();

            // minimal/maximal number of unit cells
            minunits = (numcells <  minunits) ? numcells : minunits;
            maxunits = (numcells >  maxunits) ? numcells : maxunits;

            numcells = 0;
            idxcount++;
         }//if
         flag = kTRUE;
      }//if

      input.seekg(ctrlpos[k], ios::beg);
      input.getline(nextline, kBufSize, delim);

      probe_id   = atoi(strtok(&nextline[0], sep));
      probe_type = TString(strtok(NULL, sep));
      gc_count   = atoi(strtok(NULL, sep));
      probe_len  = atoi(strtok(NULL, sep));
      inter_pos  = atoi(strtok(NULL, sep));
      sequence   = TString(strtok(NULL, sep));

      // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
      x = Index2X(probe_id - 1);
      y = Index2Y(probe_id - 1);
//to do      x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do      y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

      // control->chip: blank
      if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

      Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
      strd = ProbeType(probe_type);
      mask = XGenomeChip::SchemeMask(eCONTROLCHIP, ctrlpbt[k]);

      // fill scheme tree
      scheme->SetUnitID(ctrlpbt[k]);
      scheme->SetX(x);
      scheme->SetY(y);
      scheme->SetProbeLength(probe_len);
      scheme->SetMask(mask);
      scheme->SetExonID(ctrlpbt[k]);
      scheme->SetProbesetID(ctrlpbt[k]);
      schemetree->Fill();

     // fill probe tree: (x,y) has identical order to (x,y) of schemetree
      probe->SetX(x);
      probe->SetY(y);
      probe->SetSequence(sequence);
      probe->SetPosition(inter_pos);
      probe->SetNumberGC(gc_count);
      probe->SetTMelting(Tm);
      probe->SetProbeType(strd);
      probetree->Fill();

      if (XManager::fgVerbose && scmcount%10000 == 0) {
         cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
      }//if

      numcells++;
      scmcount++;
   }//for_i

   // delete here to free memory
   if (ctrlval) {delete [] ctrlval; ctrlval = 0;}
   if (ctrlidx) {delete [] ctrlidx; ctrlidx = 0;}
   if (ctrlpbt) {delete [] ctrlpbt; ctrlpbt = 0;}
   if (ctrlpos) {delete [] ctrlpos; ctrlpos = 0;}


//--------------------------------------------------------------------
// 2. Fill tree(s) with data for probeset types: 
//    control->bgp->antigenomic, control->bgp->genomic
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: control->bgp..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for probe_types of antigenomic, genomic
   for (Int_t i=0; i<size; i++) {
      if (!((pst[i] == eANTIGENOMIC) || (pst[i] == eGENOMIC))) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            // control->chip: blank
            if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = XGenomeChip::SchemeMask(pst[i], 0);

            // fill scheme tree
            scheme->SetUnitID(-ctlcount);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            scheme->SetExonID(-ctlcount);
            scheme->SetProbesetID(-ctlcount);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if (XManager::fgVerbose && scmcount%10000 == 0) {
               cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
            }//if
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(-ctlcount);
      unit->SetSubUnitID(probeset_id);
      unit->SetUnitType(pst[i]);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unit->SetNumSubunits(numexons);
      unittree->Fill();

      // fill exon unit tree
      exon->SetUnitName(str);
      exon->SetUnitID(-ctlcount);
      exon->SetSubUnitID(probeset_id);
      exon->SetUnitType(pst[i]);
      exon->SetNumCells(numcells);
      exon->SetNumAtoms(numatoms);
      exon->SetNumSubunits(numexons);
      exontree->Fill();

      // fill probeset unit tree
      pset->SetUnitName(str);
      pset->SetUnitID(-ctlcount);
      pset->SetSubUnitID(probeset_id);
      pset->SetUnitType(pst[i]);
      pset->SetNumCells(numcells);
      pset->SetNumAtoms(numatoms);
      pset->SetNumSubunits(numexons);
      psettree->Fill();

      // minimal number of cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      minexons = (numcells <  minexons) ? numcells : minexons;
      minpsets = (numcells <  minpsets) ? numcells : minpsets;

      // maximal number of cells
      maxunits = (numcells >  maxunits) ? numcells : maxunits;
      maxexons = (numcells >  maxexons) ? numcells : maxexons;
      maxpsets = (numcells >  maxpsets) ? numcells : maxpsets;

      ctlcount--;
      idxcount++;
   }//for_i

// Set number of controls
   fNControls = idxcount;


//--------------------------------------------------------------------
// 3. Fill tree(s) with data for probeset type: control->affx
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: control->affx..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for probe_types of bgp, normgene, rescue
   unitID   = 0;
   numexons = 1;
   for (Int_t i=0; i<size; i++) {
      if (pst[i] != eCONTROLAFFX) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      unittype = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            // control->chip: blank
            if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = XGenomeChip::SchemeMask(pst[i], strd);

            // use only mask for PM not for MM
            unittype = (unittype > mask) ? unittype : mask;

            // fill scheme tree
            scheme->SetUnitID(unitID);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            scheme->SetExonID(unitID);
            scheme->SetProbesetID(unitID);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if (XManager::fgVerbose && scmcount%10000 == 0) {
               cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
            }//if
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(unitID);
      unit->SetSubUnitID(probeset_id);
      unit->SetUnitType(unittype);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unit->SetNumSubunits(numexons);
      unittree->Fill();

      // fill exon unit tree
      exon->SetUnitName(str);
      exon->SetUnitID(unitID);
      exon->SetSubUnitID(probeset_id);
      exon->SetUnitType(unittype);
      exon->SetNumCells(numcells);
      exon->SetNumAtoms(numatoms);
      exon->SetNumSubunits(numexons);
      exontree->Fill();

      // fill probeset unit tree
      pset->SetUnitName(str);
      pset->SetUnitID(unitID);
      pset->SetSubUnitID(probeset_id);
      pset->SetUnitType(unittype);
      pset->SetNumCells(numcells);
      pset->SetNumAtoms(numatoms);
      pset->SetNumSubunits(numexons);
      psettree->Fill();

      // minimal number of cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      minexons = (numcells <  minexons) ? numcells : minexons;
      minpsets = (numcells <  minpsets) ? numcells : minpsets;

      // maximal number of cells
      maxunits = (numcells >  maxunits) ? numcells : maxunits;
      maxexons = (numcells >  maxexons) ? numcells : maxexons;
      maxpsets = (numcells >  maxpsets) ? numcells : maxpsets;

      unitID++;
      idxcount++;
   }//for_i

// Number of AFFX controls
   if (XManager::fgVerbose) {
      cout << "   Number of control->affx probesets is <" << unitID << ">." << endl;
   }//if

// Check if total number of AFFX controls is equal to fNAffx
   if ((fNAffx > 0) && (fNAffx != unitID)) {
      cerr << "Error: Number of control->affx imported <" << unitID  
           << "> is not equal to number of annotated AFFX controls <" << fNAffx << ">."
           << endl;
//??      fNAffx = unitID;
      err = errCDFVersion;
      goto cleanup;
   } else if (fNAffx == 0) {
      cout << "Warning: Missing annotation for control->affx probesets." << endl;
      fNAffx = unitID;
   }//if


//--------------------------------------------------------------------
// 4. Fill tree(s) with data for probeset annotation tree
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for probeset type: main..." << endl;
   }//if

// Sort for probeset_id (necessary for BinarySearch below)
   TMath::Sort(size, psi, index, 0);

// Replace data in arrays with sorted data using temporary array arr
   // create temporary array (one temporary array only to save memory)
   if (!(arr = new (nothrow) Long_t[size])) {err = errInitMemory; goto cleanup;}

   // replace probeset_id positions with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pos[index[i]];
   for (Int_t i=0; i<size; i++) pos[i] = arr[i];

   // replace probeset_ids with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = psi[index[i]];
   for (Int_t i=0; i<size; i++) psi[i] = arr[i];

   // replace probeset_types with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = pst[index[i]];
   for (Int_t i=0; i<size; i++) pst[i] = arr[i];

   // replace msk with sorted values
   for (Int_t i=0; i<size; i++) arr[i] = msk[index[i]];
   for (Int_t i=0; i<size; i++) msk[i] = arr[i];

   // delete temporary array here to free memory
   if (arr) {delete [] arr; arr = 0;}

// Initialize memory for exontree arrays
   if (!(arrPSet  = new (nothrow) Int_t[nentries]))   {err = errInitMemory; goto cleanup;}
   if (!(arrCount = new (nothrow) Int_t[nentries]))   {err = errInitMemory; goto cleanup;}
   if (!(arrGene  = new (nothrow) Int_t[nentries]))   {err = errInitMemory; goto cleanup;}
   if (!(arrExon  = new (nothrow) Int_t[nentries]))   {err = errInitMemory; goto cleanup;}
   if (!(arrXType = new (nothrow) Short_t[nentries])) {err = errInitMemory; goto cleanup;}
   if (!(arrLevel = new (nothrow) Short_t[nentries])) {err = errInitMemory; goto cleanup;}
   if (!(arrBound = new (nothrow) Short_t[nentries])) {err = errInitMemory; goto cleanup;}

// Fill arrays with psanntree entries
   for (Int_t i=0; i<nentries; i++) {
      psanntree->GetEntry(i);
      arrPSet[i]  = psann->GetProbesetID();
      arrCount[i] = psann->GetNumProbes();
      arrGene[i]  = psann->GetTranscriptID();
      arrExon[i]  = psann->GetExonID();
      arrXType[i] = psann->GetCrossHybType();
      arrLevel[i] = psann->GetLevelID();
      arrBound[i] = psann->GetBounded();

      if (XManager::fgVerbose && i%10000 == 0) {
         cout << "   <" << i+1 << "> probeset tree entries read...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << nentries << "> probeset tree entries read...Finished" << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

   idx = 0;
   intexonID = fNAffx; //continue after affx controls
   intpsetID = fNAffx; //continue after affx controls
   numexons  = 1;
   numcelex  = 0;

   // control->affx are already imported
   for (Int_t i=fNAffx; i<nentries; i++) {
      if (arrGene[i] == 0) continue; //rescue->FLmRNA->unmapped

      geneID = arrGene[i];
      exonID = arrExon[i];
      psetID = arrPSet[i];

      Int_t z = (Int_t)(TMath::BinarySearch((Long64_t)size, psi, arrPSet[i]));
      if (psi[z] != arrPSet[i]) continue;

      numcells = 0;
      numatoms = 0;
      numexons = 0;
      numperex = 0;
      didWhile = kFALSE;
      exonbitmsk.ResetBit(XBitSet::kBitMask);
      transbitmsk.ResetBit(XBitSet::kBitMask);
      while (geneID == arrGene[i]) {
         if (didWhile) {  // loop at least once
            if (numperex == 0) intexonID++;
            intpsetID++;
            idx++;
            i++;

            z = (Int_t)(TMath::BinarySearch((Long64_t)size, psi, arrPSet[i]));
            if (psi[z] != arrPSet[i]) continue;
            if (geneID != arrGene[i]) break;
         }//if

         input.seekg(pos[z], ios::beg);
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {err = errPrematureEOF; goto cleanup;}

         probeset_id   = atoi(strtok(&nextline[0], sep));
         probeset_type = TString(strtok(NULL, sep));

         // check for exon:main
         mainexon = (strstr(probeset_type.Data(),"exon:main") == 0) ? 0 : eMAINEXON;

         if (psi[z] != probeset_id) {
            cerr << "Error: Probeset_id <" << psi[z] << "> is not <" << probeset_id
                 << ">." << endl;
            err = errReadingInput;
            break; //ev goto cleanup;
         }//if
         msk[z] = 1;

         // minimal/maximal number of probeset cells
         minpsets = (arrCount[i] <  minpsets) ? arrCount[i] : minpsets;
         maxpsets = (arrCount[i] >  maxpsets) ? arrCount[i] : maxpsets;

         numprobs = 0;
         begin    = 1;
         while (begin != 0) { //loop over header0
            input.getline(nextline, kBufSize, delim);
            if (input.eof()) {
               // last line of *.pfg: need to reset input file to continue
               input.clear();  //clear all flags
               input.seekg(position, ios::beg);
               break;
            }//if
            str = RemoveEnds(&nextline[0], begin, end);

//x            if (begin == 1) numatoms++; //count atoms for header1
            if (begin == 1) { // counts for header1
               numatoms++; //count number of atoms
               numprobs++; //count number of probes
            }//if

            if (begin == 2) { //data for header2
               probe_id   = atoi(strtok((char*)str.Data(), sep));
               probe_type = TString(strtok(NULL, sep));
               gc_count   = atoi(strtok(NULL, sep));
               probe_len  = atoi(strtok(NULL, sep));
               inter_pos  = atoi(strtok(NULL, sep));
               sequence   = TString(strtok(NULL, sep));

               // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
               x = Index2X(probe_id - 1);
               y = Index2Y(probe_id - 1);
//to do               x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_i);
//to do               y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

               // control->chip: blank
               if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

               // get unitID and mask
//??               Type2SchemeID(probeset_type, probe_type, unitID, mask);

               Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
               strd = ProbeType(probe_type);

               // main/intron/exon vs unmapped
               if (pst[z] == eMAIN || mainexon) {
                  mask = this->SchemeMask(arrXType[i], arrLevel[i], arrBound[i]);
               } else {
                  mask = XGenomeChip::SchemeMask(pst[z], arrLevel[i]);
               }//if
               // convert for XBitSet: exonbitmsk, transbitmsk
               mask = (mask >= 0) ? mask : (abs(mask) << 15);

               // fill scheme tree
               scheme->SetUnitID(unitID);
               scheme->SetX(x);
               scheme->SetY(y);
               scheme->SetProbeLength(probe_len);
               scheme->SetMask(mask);
               scheme->SetExonID(intexonID);
               scheme->SetProbesetID(intpsetID);
               schemetree->Fill();

              // fill probe tree: (x,y) has identical order to (x,y) of schemetree
               probe->SetX(x);
               probe->SetY(y);
               probe->SetSequence(sequence);
               probe->SetPosition(inter_pos);
               probe->SetNumberGC(gc_count);
               probe->SetTMelting(Tm);
               probe->SetProbeType(strd);
               probetree->Fill();

               exonbitmsk.SetBit(mask);

               if (XManager::fgVerbose && scmcount%10000 == 0) {
                  cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
               }//if
               numcells++;
               numcelex++;
               scmcount++;
            }//if
         }//while

         // fill probeset unit tree
         str = TString(""); str += probeset_id;
         pset->SetUnitName(str);
         pset->SetUnitID(intpsetID);
         pset->SetSubUnitID(probeset_id);
         pset->SetUnitType(mask);
//x         pset->SetNumCells(arrCount[i]);
//x         pset->SetNumAtoms(arrCount[i]);
         pset->SetNumCells(numprobs);
         pset->SetNumAtoms(numprobs);
         pset->SetNumSubunits(1);
         psettree->Fill();

         if (XManager::fgVerbose && arrCount[i] != numprobs) {
            cout << "Warning: number of  probes <" << numprobs << "> for probeset <"
                 << probeset_id << "> is not equal to annotated probe_count <"
                 << arrCount[i] << ">" << endl;
         }//if

         if (exonID == arrExon[i+1]) {
            numperex++;
         } else {
            exonlevel = exonbitmsk.TestBits(XBitSet::kBitMask);

            // fill exon unit tree
            str = TString(""); str += exonID;
            exon->SetUnitName(str);
            exon->SetUnitID(intexonID);
            exon->SetSubUnitID(exonID);
            exon->SetUnitType(exonlevel);
            exon->SetNumCells(numcelex);
            exon->SetNumAtoms(numcelex);
            exon->SetNumSubunits(numperex + 1);
            exontree->Fill();

            // minimal/maximal number of exon cells
            minexons = (numcelex <  minexons) ? numcelex : minexons;
            maxexons = (numcelex >  maxexons) ? numcelex : maxexons;

            numexons++;
            exonID = arrExon[i+1];

            exonbitmsk.ResetBit(XBitSet::kBitMask);

            numperex = 0;
            numcelex = 0;
         }//if

         transbitmsk.SetBit(mask);

         didWhile = kTRUE;
      }//while

      // count number of exons
      exncount += numexons;

      translevel = transbitmsk.TestBits(XBitSet::kBitMask);

      // fill transcript unit tree
      str = TString(""); str += geneID;
      unit->SetUnitName(str);
      unit->SetUnitID(unitID);
      unit->SetSubUnitID(geneID);
      unit->SetUnitType(translevel);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unit->SetNumSubunits(numexons);
      unittree->Fill();

      // minimal/maximal number of unit cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      maxunits = (numcells >  maxunits) ? numcells : maxunits;

      unitID++;
      idxcount++;
      if (didWhile) i--;
   }//for_i

// Check if total number of exons is equal to fNExons (determined in ImportExonAnnotation)
   if (exncount != fNExons) {
      if (XManager::fgVerbose) {
         cout << "   Note: Number of exons imported <" << exncount  
              << "> is not equal to number of annotated exons <" << fNExons << ">."
              << endl;
      }//if
//??      fNExons = exncount;
   }//if

// Set number of genes
   fNGenes = unitID - fNAffx;


//--------------------------------------------------------------------
// 4. Fill tree(s) with data for non-annotated probesets
//--------------------------------------------------------------------
   if (XManager::fgVerbose) {
      cout << "   Filling trees with data for non-annotated probesets..." << endl;
   }//if

   // reset input file to first probeset_id
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Get positions and probe_type IDs for non-annotated probeset_ids
//!! fill tree.scm and tree.prb with rest of probe_set data!
  numexons = 1;
  for (Int_t i=0; i<size; i++) {
      if (msk[i] == 1) continue;

      input.seekg(pos[i], ios::beg);
      input.getline(nextline, kBufSize, delim);
      if (input.eof()) {err = errPrematureEOF; goto cleanup;}

      probeset_id   = atoi(strtok(&nextline[0], sep));
      probeset_type = TString(strtok(0, sep));
      probesetype   = ProbesetType(probeset_type);

      if (psi[i] != probeset_id) {
         cerr << "Error: Probeset_id is not <" << probeset_id << ">." << endl;
         err = errReadingInput;
         break; //ev goto cleanup;
      }//if
      msk[i] = 1;

      begin    = 1;
      numcells = 0;
      numatoms = 0;
      while (begin != 0) { //loop over header0
         input.getline(nextline, kBufSize, delim);
         if (input.eof()) {
            // last line of *.pfg: need to reset input file to continue
            input.clear();  //clear all flags
            input.seekg(position, ios::beg);
            break;
         }//if
         str = RemoveEnds(&nextline[0], begin, end);

         if (begin == 1) numatoms++; //count atoms for header1
         if (begin == 2) { //data for header2
            probe_id   = atoi(strtok((char*)str.Data(), sep));
            probe_type = strtok(NULL, sep);
            gc_count   = atoi(strtok(NULL, sep));
            probe_len  = atoi(strtok(NULL, sep));
            inter_pos  = atoi(strtok(NULL, sep));
            sequence   = TString(strtok(NULL, sep));

            // get (x,y) coordinates (1st probe_id=1 at x=0, y=0 => probe_id - 1)
            x = Index2X(probe_id - 1);
            y = Index2Y(probe_id - 1);
//to do            x = fIsSequential ? Index2X(probe_id - 1) : LayoutToX(probe_id);
//to do            y = fIsSequential ? Index2Y(probe_id - 1) : LayoutToY(probe_id);

            // control->chip: blank
            if (strcmp(sequence.Data(), "!")  == 0) sequence = "N";

            Tm   = MeltingTemperature(gc_count, probe_len, "empirical");
            strd = ProbeType(probe_type);
            mask = XGenomeChip::SchemeMask(pst[i], strd);

            // fill scheme tree
            scheme->SetUnitID(unitID);
            scheme->SetX(x);
            scheme->SetY(y);
            scheme->SetProbeLength(probe_len);
            scheme->SetMask(mask);
            scheme->SetExonID(intexonID);
            scheme->SetProbesetID(intpsetID);
            schemetree->Fill();

           // fill probe tree: (x,y) has identical order to (x,y) of schemetree
            probe->SetX(x);
            probe->SetY(y);
            probe->SetSequence(sequence);
            probe->SetPosition(inter_pos);
            probe->SetNumberGC(gc_count);
            probe->SetTMelting(Tm);
            probe->SetProbeType(strd);
            probetree->Fill();

            if (XManager::fgVerbose && scmcount%10000 == 0) {
               cout << "   <" << scmcount + 1 << "> records imported...\r" << flush;
            }//if
            numcells++;
            scmcount++;
         }//if
      }//while

      // fill transcript unit tree
      str = TString(""); str += probeset_id;
      unit->SetUnitName(str);
      unit->SetUnitID(unitID);
      unit->SetSubUnitID(probeset_id);
      unit->SetUnitType(probesetype);
      unit->SetNumCells(numcells);
      unit->SetNumAtoms(numatoms);
      unit->SetNumSubunits(numexons);
      unittree->Fill();

      // fill exon unit tree
      exon->SetUnitName(str);
      exon->SetUnitID(intexonID);
      exon->SetSubUnitID(probeset_id);
      exon->SetUnitType(probesetype);
      exon->SetNumCells(numcells);
      exon->SetNumAtoms(numatoms);
      exon->SetNumSubunits(numexons);
      exontree->Fill();

      // fill probeset unit tree
      pset->SetUnitName(str);
      pset->SetUnitID(intpsetID);
      pset->SetSubUnitID(probeset_id);
      pset->SetUnitType(probesetype);
      pset->SetNumCells(numcells);
      pset->SetNumAtoms(numatoms);
      pset->SetNumSubunits(numexons);
      psettree->Fill();

      // minimal number of cells
      minunits = (numcells <  minunits) ? numcells : minunits;
      minexons = (numcells <  minexons) ? numcells : minexons;
      minpsets = (numcells <  minpsets) ? numcells : minpsets;

      // maximal number of cells
      maxunits = (numcells >  maxunits) ? numcells : maxunits;
      maxexons = (numcells >  maxexons) ? numcells : maxexons;
      maxpsets = (numcells >  maxpsets) ? numcells : maxpsets;

      exncount++;
      idxcount++;
      intexonID++;
      intpsetID++;
      unitID++;
   }//for_i

// Set number of probes
   if (XManager::fgVerbose) {
      cout << "   <" << scmcount << "> records imported...Finished" << endl;
   }//if
   fNProbes = scmcount;

// Set number of transcript units
   if (XManager::fgVerbose) {
      cout << "   <" << idxcount << "> total transcript units imported." << endl;
   }//if
   fNUnits = idxcount;

// Set number of exon units
   fNExonUnits = exncount + fNControls + fNAffx;
   if (XManager::fgVerbose) {
      cout << "   <" << fNExonUnits << "> total exon array units imported." << endl;
   }//if

   if (XManager::fgVerbose) {
      cout << "   Exon cell statistics: " << endl;
      cout << "      Number of probeset cells: minimum = " << minpsets << ",  maximum = " << maxpsets << endl;
      cout << "      Number of exon cells:     minimum = " << minexons << ",  maximum = " << maxexons << endl;
      cout << "      Number of unit cells:     minimum = " << minunits << ",  maximum = " << maxunits << endl;
   }//if

//?? or to unittree
// Add tree info to tree
   AddSchemeTreeInfo(schemetree, schemetree->GetName());

// Write scheme tree to file
   if ((err = WriteTree(schemetree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(schemetree->GetName(), 0);
   }//if
   //delete tree from memory
   schemetree->Delete("");
   schemetree = 0;

// Add tree info to unittree
   AddGenomeTreeInfo(unittree, unittree->GetName(), "", fNControls, fNAffx, fNGenes,
                     minunits, maxunits);

// Write unit tree to file
   if ((err = WriteTree(unittree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(unittree->GetName(), 0);
   }//if
   //delete tree from memory
   unittree->Delete("");
   unittree = 0;

// Add tree info to probetree
   AddProbeTreeInfo(probetree, fPrbTreeName);

// Write probe tree to file
   if ((err = WriteTree(probetree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(probetree->GetName(), 0);
   }//if
   //delete tree from memory
   probetree->Delete("");
   probetree = 0;

// Add tree info to exontree
   AddGenomeTreeInfo(exontree, fExnTreeName, "", fNControls, fNAffx, exncount,
                     minexons, maxexons);

// Write exon tree to file
   if ((err = WriteTree(exontree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(exontree->GetName(), 0);
   }//if
   //delete tree from memory
   exontree->Delete("");
   exontree = 0;

// Add tree info to psettree
   AddGenomeTreeInfo(psettree, fPbsTreeName, "", fNControls, fNAffx, idx,
                     minpsets, maxpsets);

// Write exon tree to file
   if ((err = WriteTree(psettree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(psettree->GetName(), 0);
   }//if
   //delete tree from memory
   psettree->Delete("");
   psettree = 0;

//Clean up
cleanup:
   if (arrBound) {delete [] arrBound; arrBound = 0;}
   if (arrLevel) {delete [] arrLevel; arrLevel = 0;}
   if (arrXType) {delete [] arrXType; arrXType = 0;}
   if (arrExon)  {delete [] arrExon;  arrExon  = 0;}
   if (arrGene)  {delete [] arrGene;  arrGene  = 0;}
   if (arrCount) {delete [] arrCount; arrCount = 0;}
   if (arrPSet)  {delete [] arrPSet;  arrPSet  = 0;}

   if (ctrlval) {delete [] ctrlval; ctrlval = 0;}
   if (ctrlidx) {delete [] ctrlidx; ctrlidx = 0;}
   if (ctrlpbt) {delete [] ctrlpbt; ctrlpbt = 0;}
   if (ctrlpos) {delete [] ctrlpos; ctrlpos = 0;}

   if (value) {delete [] value; value = 0;}
   if (index) {delete [] index; index = 0;}

   if (arr) {delete [] arr; arr = 0;}
   if (msk) {delete [] msk; msk = 0;}
   if (pst) {delete [] pst; pst = 0;}
   if (psi) {delete [] pst; pst = 0;}
   if (pos) {delete [] pos; pos = 0;}

   SafeDelete(probe);
   SafeDelete(pset);
   SafeDelete(exon);
   SafeDelete(unit);
   SafeDelete(scheme);

   return err;
}//ReadData

//______________________________________________________________________________
Int_t XExonChip::ImportTransAnnotation(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree "tree.ann":
   if(kCS) cout << "------XExonChip::ImportTransAnnotation------" << endl;

   Int_t err  = errNoErr;
   Int_t size = 0;
   Int_t idx  = 0;
   Int_t nsub = 5;  //number of sub-fields

   //annotation file is at least version 1.8 (human), 1.3 (mouse), 1.3 (rat)?
   Bool_t hasPST = kFALSE;
   Bool_t isNA22 = kFALSE; //is at least na22

   std::string nextline;  //to read next line
   streampos position;

   TString lib_set_name, lib_set_version;
   TString genome_species, genome_version;
   TString annot_version, annot_type, annps_version;
   TString assigngene, assignmrna, category;
   TString str, dummy;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

// Init variables to tokenize lines
   const char *csv = "\",\""; //necessary to tokenize lines, better?: sep = "\",\"";
   const char *tab = "\t";    //replace csv with tab

// Init local arrays to store data from input columns
   Int_t   *geneid  = 0;  //transcript_cluster_id
   TString *seqname = 0;  //seqname, i.e. chromosome
   TString *strand  = 0;  //strand of chromosome
   Int_t   *start   = 0;  //start position on chromosome
   Int_t   *stop    = 0;  //stop position on chromosome
//   Int_t   *nprobes = 0;  //total number of probes
   Int_t    nprobes;

// Init local arrays for derived data
   TString *names      = 0;  //array containing tokenized string
   TString *geneaccess = 0;  //gene_assignment, accession number
   TString *genename   = 0;  //gene_assignment, name
   TString *genesymbol = 0;  //gene_assignment, symbol
   TString *cytoband   = 0;  //gene_assignment, cytoband
   Int_t   *entrezid   = 0;  //gene_assignment, entrez ID
   Short_t *psettype   = 0;  //probeset_type (main, normgene->intron, control->affx, etc)

// Init local arrays for sorting
   Int_t    *index     = 0;  //sort index
   Long64_t *value     = 0;  //value array to store the two columns to sort
   Long64_t majorv, minorv;  //values to store in array value

// Create new transcript annotation tree
   fAnnTreeName = TString(fName) + "." + exten;
   TTree *anntree = new TTree(fAnnTreeName, "transcript annotation");
   if (anntree == 0) return errCreateTree;
   XGenomeAnnotation *ann = 0;
   ann = new XGenomeAnnotation();
   anntree->Branch("AnnBranch", "XGenomeAnnotation", &ann, 64000, split);

// Check for line "#%lib-set-name"
   while ((nextline.compare(0, 14, "#%lib-set-name") != 0) &&
          (nextline.compare(0, 14, "#%lib_set_name") != 0)) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_name = strtok((&((char*)nextline.c_str())[15]), sep);

// Check for line "#%lib-set-version"
   input.clear();
   input.seekg(position, ios::beg);
   while ((nextline.compare(0, 17, "#%lib-set-version") != 0) &&
          (nextline.compare(0, 17, "#%lib_set_version") != 0)) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_version = strtok((&((char*)nextline.c_str())[18]), sep);

// Check for line "#%genome-species"
   Bool_t hasSpecies = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 16, "#%genome-species") != 0) {
      std::getline(input, nextline, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasSpecies = kFALSE; break;}
   }//while
   if (hasSpecies == kTRUE) {
      genome_species = strtok((&((char*)nextline.c_str())[17]), sep);
   } else {
      cout << "   Note: The following header line is missing: %genome-species=" << endl;
      genome_species = "NA";
   }//if

// Check for line "#%genome-version"
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 16, "#%genome-version") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   genome_version = strtok((&((char*)nextline.c_str())[17]), sep);

// Check for line "#%netaffx-annotation-netaffx-build"
   Bool_t hasNetBuild = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 34, "#%netaffx-annotation-netaffx-build") != 0) {
      std::getline(input, nextline, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasNetBuild = kFALSE; break;}
   }//while
   if (hasNetBuild == kTRUE) {
      annot_version = strtok((&((char*)nextline.c_str())[35]), sep);
   } else {
      cout << "   Note: The following header line is missing: %netaffx-annotation-netaffx-build=" << endl;
      annot_version = "NA";
   }//if

   // compare annotation versions, but only main version, e.g. 27 for 27.2
   annps_version = fVersionAnnot;
   if (annot_version.Index(".") > 0) annot_version.Resize(annot_version.Index("."));
   if (annps_version.Index(".") > 0) annps_version.Resize(annps_version.Index("."));
   if (strcmp(annot_version.Data(), annps_version.Data()) != 0) {
      return fManager->HandleError(errAnnVersion, annot_version);
   }//if

   // is at least na22
   if (atoi(annot_version.Data()) >= kVersionAnnot22) {isNA22 = kTRUE;}

// Check for line "#%netaffx-annotation-data-type"
   Bool_t hasNetType = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 30, "#%netaffx-annotation-data-type") != 0) {
      std::getline(input, nextline, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasNetType = kFALSE; break;}
   }//while
   if (hasNetType == kTRUE) {
      annot_type = strtok((&((char*)nextline.c_str())[31]), sep);
   } else {
      cout << "   Note: The following header line is missing: %netaffx-annotation-data-type=" << endl;
      annot_type = "NA";
   }//if

// Check for line ""transcript_cluster_id""
   input.clear();
   input.seekg(position, ios::beg);
   while (nextline.compare(0, 23, "\"transcript_cluster_id\"") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   Int_t numsep = NumSeparators(nextline.c_str(), sep);

   // remove all "\"" from header line
   str = TString(nextline);
   str.ReplaceAll("\"", "");

// Create array hasColumn 
   Int_t *hasColumn = 0;
   if (!(hasColumn = new (nothrow) Int_t[kNTranscriptCols])) return errInitMemory;
   for (Int_t i=0; i<kNTranscriptCols; i++) {
      hasColumn[i] = 0;
   }//for_i

//TO DO: see ImportProbesetAnnotation()
//better   err = CheckHeaderOrder(str, kTranscriptHeader, kNTranscriptCols, hasColumn, sep);
// Check column headers from header line, hasColumn = 1 if column is present
   err = CheckHeader(str, kTranscriptHeader, kNTranscriptCols, hasColumn, sep);
   if (err > 0) {
      cout << "Note: The following header columns are missing: " << endl;
      for (Int_t i=0; i<kNTranscriptCols; i++) {
         if (hasColumn[i] == 0) cout << "<" << kTranscriptHeader[i] << ">" << endl;
         if (i == 0 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //transcript_cluster_id
         if (i == 2 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //seqname
         if (i == 3 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //strand
         if (i == 4 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //start
         if (i == 7 && hasColumn[i] == 0) {err = errMissingColumn; goto cleanup;} //gene_assignment
      }//for_i
   }//if 

// Check if annotation file is at least build version na22 with new column category
//no   if (hasColumn[kNTranscriptCols-1] == kNTranscriptCols-1) {hasPST = kTRUE;} //category
   if (hasColumn[kNTranscriptCols-2] == 1) {hasPST = kTRUE;} //crosshyb_type

// Get number of transcripts
   // get current streamposition for rewinding input
   position = input.tellg();
   while (1) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      size++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   Number of annotated transcripts is <" << size << ">." << endl;
   }//if

   // reset input file to position
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Initialize memory for local arrays
   if (!(geneid     = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(seqname    = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(strand     = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(start      = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(stop       = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
//   if (!(nprobes    = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}

   if (!(names      = new (nothrow) TString[nsub])) {err = errInitMemory; goto cleanup;}
   if (!(geneaccess = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(genename   = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(genesymbol = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(cytoband   = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(entrezid   = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(psettype   = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize memory for sorting arrays
   if (!(index = new (nothrow) Int_t[size]))        {err = errInitMemory; goto cleanup;}
   if (!(value = new (nothrow) Long64_t[size]))     {err = errInitMemory; goto cleanup;}

   while (input.good()) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      if (input.fail()) {
         cout << "Error: Failed reading line:" << endl;
         cout << nextline << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // replace all "\",\"" with tab "\t" and remove "\""
      str = TString(nextline);
      // first need to replace "" with NA
      str.ReplaceAll("\"\"", "\"NA\"");
      // replace tab with space to eliminate wrong tabs in Affymetrix annotation files
      str.ReplaceAll(tab, " ");
      // replace csv with tab
      str.ReplaceAll(csv, tab);
      // remove all "\"" from line
      str.ReplaceAll("\"", "");

      // check number of separators
      if (numsep != NumSeparators(str, tab)) {
         cout << "Error: Wrong number of separators in line:" << endl;
         cout << str.Data() << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // import fields
      dummy        = strtok((char*)str.Data(), tab);
      geneid[idx]  = atoi((RemoveEnds(dummy)).Data());
      dummy        = strtok(NULL, tab);
      seqname[idx] = strtok(NULL, tab);
      strand[idx]  = strtok(NULL, tab);
      start[idx]   = atoi(strtok(NULL, tab));
      stop[idx]    = atoi(strtok(NULL, tab));
      // need to test for nprobes==0 because of bug in transcript-annotation files from March 2007
      nprobes      = atoi(strtok(NULL, tab));
      if (nprobes == 0) continue;
      assigngene   = strtok(NULL, tab);
      assignmrna   = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      dummy        = strtok(NULL, tab);
      category     = hasPST ? strtok(NULL, tab) : dummy; //for >=na22

      // chromosome
      seqname[idx] = strcmp(seqname[idx].Data(),"---") ? seqname[idx] : "NA";
      strand[idx]  = strcmp(strand[idx].Data(), "---") ? strand[idx]  : "?";

      // get gene accession and symbol
      geneaccess[idx] = strtok((char*)assigngene.Data(), "\\");
      if (strcmp(geneaccess[idx].Data(),"---") == 0) {
         geneaccess[idx] = "NA";
         genesymbol[idx] = "NA";
         genename[idx]   = "NA";
         cytoband[idx]   = "NA";
         entrezid[idx]   = -1;
      } else {
         nsub = 5;  //since TokenizeString() can change nsub
         for (Int_t i=0; i<nsub; i++) names[i] = "---";

         Int_t index = 0;
         index = assigngene.Index(kSepSl3, kNumSl3, index, TString::kExact);
         // index>0 only if assigngene is a multipart entry
         assigngene = (index > 0) ? assigngene(0, index) : assigngene;

         index = TokenizeString(assigngene.Data(), nsub, names, kNumSl2, kSepSl2);

         geneaccess[idx] = (strcmp(names[0].Data(),"---") != 0) ? names[0] : "NA";
         genesymbol[idx] = (strcmp(names[1].Data(),"---") != 0) ? names[1] : "NA";
         genename[idx]   = (strcmp(names[2].Data(),"---") != 0) ? names[2] : "NA";
         cytoband[idx]   = (strcmp(names[3].Data(),"---") != 0) ? ("chr" + names[3]) : "NA";
         entrezid[idx]   = (strcmp(names[4].Data(),"---") != 0) ? names[4].Atoi() : -1;
      }//if

      // convert category to probesettype_id
      psettype[idx] = this->ProbesetType(category);

      // get mrna accession for genes or for control->affx 
      if (psettype[idx] == eCONTROLAFFX) {
         geneaccess[idx] = RemoveEnds(assignmrna);
      }//if

      // get mrna accession for genes or for control->affx 
      if ((psettype[idx] == eINTRON) ||
          (psettype[idx] == eEXON)   ||
          (psettype[idx] == eUNMAPPED)) {
         nsub = 5;  //since TokenizeString() can change nsub
         for (Int_t i=0; i<nsub; i++) names[i] = "---";

         Int_t index = 0;
         index = assignmrna.Index(kSepSl3, kNumSl3, index, TString::kExact);
         // index>0 only if assignmrna is a multipart entry
         assignmrna = (index > 0) ? assignmrna(0, index) : assignmrna;

         index = TokenizeString(assignmrna.Data(), nsub, names, kNumSl2, kSepSl2);

         geneaccess[idx] = (strcmp(names[0].Data(),"---") != 0) ? names[0] : "NA";
         genesymbol[idx] = (strcmp(names[1].Data(),"---") != 0) ? names[1] : "NA";
         genename[idx]   = (strcmp(names[2].Data(),"---") != 0) ? names[2] : "NA";
         cytoband[idx]   = (strcmp(names[3].Data(),"---") != 0) ? ("chr" + names[3]) : "NA";
         entrezid[idx]   = (strcmp(names[4].Data(),"---") != 0) ? names[4].Atoi() : -1;
      }//if

      // fill sort values
      majorv = (Long64_t)geneid[idx]; //sort for transcript_cluster_id first
      minorv = (Long64_t)start[idx];  //then sort for start position
      value[idx]  = majorv << 31;
      value[idx] += minorv;

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx + 1 << "> records read...\r" << flush;
      }//if
      idx++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records read...Finished" << endl;
   }//if

   // set size to number of records read
   size = idx;

// Sort for geneid first and then for start (assuming transcipt on one chromosome)
   TMath::Sort(size, value, index, 0);

// Fill probeset/exon annotation tree with control->affx 
   idx = 0;
   if ((hasPST == kTRUE) && (fNAffx == 0)) {
      for (Int_t i=0; i<size; i++) {
         Int_t k = index[i];
         if (psettype[k] != eCONTROLAFFX) continue;

         str.Form("%d", geneid[k]);

         ann->SetUnitID(idx++);
         ann->SetTranscriptID(str); //may interfere with transcript_cluster_id???
         ann->SetName(geneaccess[k]);
         ann->SetSymbol("NA");
         ann->SetAccession(str);
         ann->SetEntrezID(-1);
         ann->SetChromosome("NA");
         ann->SetCytoBand("NA");
         ann->SetStrand('?');
         ann->SetStart(0);
         ann->SetStop(0);
         ann->SetCrossHybType(0);
         ann->SetProbesetType(eCONTROLAFFX);
         anntree->Fill();

         fNAffx++;
      }//for_i
   } else if ((fNAffx > 0) && (fAffxNames->GetSize() == fNAffx)) {
      for (Int_t i=0; i<fNAffx; i++) {
         XIdxString *idxstr   = (XIdxString*)(fAffxNames->At(i));
         TString     affxname = idxstr->GetString();

         // affxid as string
         str.Form("%d", idxstr->GetIndex());

         // fill annotation tree with data
         ann->SetUnitID(idx++);
         ann->SetTranscriptID(str); //may interfere with transcript_cluster_id???
         ann->SetName(affxname);
         ann->SetSymbol("NA");
         ann->SetAccession(str);
         ann->SetEntrezID(-1);
         ann->SetChromosome("NA");
         ann->SetCytoBand("NA");
         ann->SetStrand('?');
         ann->SetStart(0);
         ann->SetStop(0);
         ann->SetCrossHybType(0);
         ann->SetProbesetType(eCONTROLAFFX);
         anntree->Fill();
      }//for_i
   }//if

// Fill transcript annotation tree with sorted annotation data
   for (Int_t i=0; i<size; i++) {
      Int_t k = index[i];
      if (!(isNA22 && psettype[k] == eMAIN)) continue;

      // geneid as string
      str.Form("%d", geneid[k]);

      // fill annotation tree with data
      ann->SetUnitID(idx++);
      ann->SetTranscriptID(str);
      ann->SetName(genename[k]);
      ann->SetSymbol(genesymbol[k]);
      ann->SetAccession(geneaccess[k]);
      ann->SetEntrezID(entrezid[k]);
      ann->SetChromosome(seqname[k]);
      ann->SetCytoBand(cytoband[k]);
      ann->SetStrand((strand[k].Data())[0]);
      ann->SetStart(start[k]);
      ann->SetStop(stop[k]);
      ann->SetCrossHybType(0);
      ann->SetProbesetType(psettype[k]);
      anntree->Fill();

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx << "> records imported...\r" << flush;
      }//if
   }//for_i

// Fill transcript annotation tree with sorted "rescue->FLmRNA->unmapped" annotation data
   for (Int_t i=0; i<size; i++) {
      Int_t k = index[i];
      if (isNA22 == kFALSE) break;
      if (!((psettype[k] == eINTRON) ||
            (psettype[k] == eEXON)   ||
            (psettype[k] == eUNMAPPED))) continue;

      // geneid as string
      str.Form("%d", geneid[k]);

      // fill annotation tree with data
      ann->SetUnitID(idx++);
      ann->SetTranscriptID(str);
      ann->SetName(genename[k]);
      ann->SetSymbol(genesymbol[k]);
      ann->SetAccession(geneaccess[k]);
      ann->SetEntrezID(entrezid[k]);
      ann->SetChromosome(seqname[k]);
      ann->SetCytoBand(cytoband[k]);
      ann->SetStrand((strand[k].Data())[0]);
      ann->SetStart(start[k]);
      ann->SetStop(stop[k]);
      ann->SetCrossHybType(0);
      ann->SetProbesetType(psettype[k]);
      anntree->Fill();

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx << "> records imported...\r" << flush;
      }//if
   }//for_i

   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records imported...Finished" << endl;
   }//if

//TO DO:
// Add tree info to tree
//   AddAnnotationTreeInfo(anntree, fAnnTreeName);

// Write annotation tree to file
   if ((err = WriteTree(anntree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(anntree->GetName(), 0);
   }//if
   //delete tree from memory
   anntree->Delete("");
   anntree = 0;

// Clean up
cleanup:
   SafeDelete(ann);

   if (value)      {delete [] value;      value      = 0;}
   if (index)      {delete [] index;      index      = 0;}
   if (psettype)   {delete [] psettype;   psettype   = 0;}
   if (entrezid)   {delete [] entrezid;   entrezid   = 0;}
   if (cytoband)   {delete [] cytoband;   cytoband   = 0;}
   if (genesymbol) {delete [] genesymbol; genesymbol = 0;}
   if (genename)   {delete [] genename;   genename   = 0;}
   if (geneaccess) {delete [] geneaccess; geneaccess = 0;}
   if (names)      {delete [] names;      names      = 0;}
   if (stop)       {delete [] stop;       stop       = 0;}
   if (start)      {delete [] start;      start      = 0;}
   if (strand)     {delete [] strand;     strand     = 0;}
   if (seqname)    {delete [] seqname;    seqname    = 0;}
   if (geneid)     {delete [] geneid;     geneid     = 0;}
   if (hasColumn)  {delete [] hasColumn;  hasColumn  = 0;}

   return err;
}//ImportTransAnnotation

//______________________________________________________________________________
Int_t XExonChip::ImportExonAnnotation(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree "tree.anx":
   if(kCS) cout << "------XExonChip::ImportExonAnnotation------" << endl;

// Currently, exon annotation is extracted from probeset annotation input
   return this->ImportProbesetAnnotation(input, option, sep, delim, split);
}//ImportExonAnnotation

//______________________________________________________________________________
Int_t XExonChip::ImportProbesetAnnotation(ifstream &input, Option_t *option, 
                 const char *sep, char delim, Int_t split)
{
   // Import annotation from infile and store as annotation tree "tree.anx":
   if(kCS) cout << "------XExonChip::ImportProbesetAnnotation------" << endl;

   char nextline[kAnnBuf];
   Int_t err  = errNoErr;
   Int_t size = 0;
   Int_t idx  = 0;

   Int_t exonID    = 0;
   Int_t numpsets  = 0;
   Int_t numexons  = 0;
   Int_t exonstart = 0;
   Int_t exonstop  = 0;

   TString opt   = Path2Name(option, dSEP, ".");
   TString exten = Path2Name(option, ".", "");

   //annotation file is at least version 1.8 (human), 1.3 (mouse), 1.3 (rat)?
   Bool_t hasPST = kFALSE;

   streampos position;

   TString lib_set_name, lib_set_version;
   TString genome_species, genome_version;
   TString annot_type;
   TString assigngene, assignmrna, level, probesettype;
   TString str, dummy;

   Short_t bound, nobound;

// Init variables to tokenize lines
   const char *csv = "\",\""; //necessary to tokenize lines, better?: sep = "\",\"";
   const char *tab = "\t";    //replace csv with tab

// Init local arrays to store data from input columns
   Int_t   *probesetid = 0;  //probeset_id
   TString *seqname    = 0;  //seqname, i.e. chromosome
   TString *strand     = 0;  //strand of chromosome
   Int_t   *start      = 0;  //start position on chromosome
   Int_t   *stop       = 0;  //stop position on chromosome
   Int_t   *nprobes    = 0;  //probe_count
   Int_t   *geneid     = 0;  //transcript_cluster_id
   Int_t   *exonid     = 0;  //exon_id
   Short_t *crosshtype = 0;  //crosshyb_type (unique, similar, mixed)
   Short_t *full       = 0;  //fl, i.e. flag for putative full-length mRNA

// Init local arrays for derived data
   TString *geneaccess = 0;  //gene_assignment, accession number
   TString *genesymbol = 0;  //gene_assignment, symbol
   TString *mrnaaccess = 0;  //mrna_assignment, accession number
   Short_t *levelid    = 0;  //level, converted to id
   Short_t *bounded    = 0;  //bounded and noBoundedEvidence
   Short_t *psettype   = 0;  //probeset_type (main, normgene->intron, control->affx, etc)

// Init local arrays for sorting
   Int_t    *index     = 0;  //sort index
   Long64_t *value     = 0;  //value array to store the two columns to sort
   Long64_t majorv, minorv;  //values to store in array value

// Check for line "#%lib-set-name"
   while ((strncmp("#%lib-set-name", nextline, 14) != 0) &&
          (strncmp("#%lib_set_name", nextline, 14) != 0)) {
      input.getline(nextline, kAnnBuf, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_name = strtok(&nextline[15], sep);

// Check for line "#%lib-set-version"
   input.clear();
   input.seekg(position, ios::beg);
   while ((strncmp("#%lib-set-version", nextline, 17) != 0) &&
          (strncmp("#%lib_set_version", nextline, 17) != 0)) {
      input.getline(nextline, kAnnBuf, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_version = strtok(&nextline[18], sep);

// Check for line "#%genome-species"
   Bool_t hasSpecies = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (strncmp("#%genome-species", nextline, 16) != 0) {
      input.getline(nextline, kAnnBuf, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasSpecies = kFALSE; break;}
   }//while
   if (hasSpecies == kTRUE) {
      genome_species = strtok(&nextline[17], sep);
   } else {
      cout << "   Note: The following header line is missing: %genome-species=" << endl;
      genome_species = "NA";
   }//if

// Check for line "#%genome-version"
   input.clear();
   input.seekg(position, ios::beg);
   while (strncmp("#%genome-version", nextline, 16) != 0) {
      input.getline(nextline, kAnnBuf, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   genome_version = strtok(&nextline[17], sep);

// Check for line "#%netaffx-annotation-netaffx-build"
   Bool_t hasNetBuild = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (strncmp("#%netaffx-annotation-netaffx-build", nextline, 34) != 0) {
      input.getline(nextline, kAnnBuf, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasNetBuild = kFALSE; break;}
   }//while
   if (hasNetBuild == kTRUE) {
      fVersionAnnot = strtok(&nextline[35], sep);
   } else {
      cout << "   Note: The following header line is missing: %netaffx-annotation-netaffx-build=" << endl;
      fVersionAnnot = "NA";
   }//if

// Check for line "#%netaffx-annotation-data-type"
   Bool_t hasNetType = kTRUE;
   input.clear();
   input.seekg(position, ios::beg);
   while (strncmp("#%netaffx-annotation-data-type", nextline, 30) != 0) {
      input.getline(nextline, kAnnBuf, delim);
//      if (input.eof()) return errPrematureEOF;
      if (input.eof()) {hasNetType = kFALSE; break;}
   }//while
   if (hasNetType == kTRUE) {
      annot_type = strtok(&nextline[31], sep);
      if ( strncmp(annot_type, "probe_set", 9) != 0) {
         cout << "Warning: Annotation data type is not <probe_set>." << endl;
//         return errReadingInput;
      }//if
   } else {
      cout << "   Note: The following header line is missing: %netaffx-annotation-data-type=" << endl;
      annot_type = "NA";
   }//if

// Check for line ""probeset_id""
   input.clear();
   input.seekg(position, ios::beg);
   while (strncmp("\"probeset_id\"", nextline, 13) != 0) {
      input.getline(nextline, kAnnBuf, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   Int_t numsep = NumSeparators(&nextline[0], sep);

   // remove all "\"" from line
   str = TString(&nextline[0]);
   str.ReplaceAll("\"", "");

   fAnxTreeName = TString(fName) + "." + kExtenAnnot[1];
   TTree *exanntree = new TTree(fAnxTreeName, "exon annotation");
   if (exanntree == 0) return errCreateTree;
   XExonAnnotation *exann = 0;
   exann = new XExonAnnotation();
   exanntree->Branch("AnnBranch", "XExonAnnotation", &exann, 64000, split);

// Create new probeset annotation tree (for use in ImportScheme())
   fAnpTreeName = TString(fName) + "." + kExtenAnnot[2];
   TTree *psanntree = new TTree(fAnpTreeName, "probeset annotation");
   if (psanntree == 0) return errCreateTree;
   XProbesetAnnotation *psann = 0;
   psann = new XProbesetAnnotation();
   psanntree->Branch("AnnBranch", "XProbesetAnnotation", &psann, 64000, split);

// Create array hasColumn 
   Int_t *hasColumn = 0;
   if (!(hasColumn = new (nothrow) Int_t[kNProbesetCols])) {err = errInitMemory; goto cleanup;}
   for (Int_t i=0; i<kNProbesetCols; i++) {
      hasColumn[i] = 0;
   }//for_i

// Check column headers from header line, hasColumn = 1 if column is present
   err = CheckHeaderOrder(str, kProbesetHeader, kNProbesetCols, hasColumn, sep);
   if (err > 0) {
      cout << "Note: The following header columns are missing or in wrong order: " << endl;
      for (Int_t i=0; i<kNProbesetCols; i++) {
         if (hasColumn[i] != i) cout << "   <" << kProbesetHeader[i] << ">" << endl;
         if (i == 0  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //probeset_id
         if (i == 1  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //seqname
         if (i == 2  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //strand
         if (i == 3  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //start
         if (i == 6  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //transcript_cluster_id
         if (i == 7  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //exon_id
         if (i == 9  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //gene_assignment
         if (i == 10 && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //mrna_assignment
         if (i == 15 && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //level
         if (i == 16 && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //bounded
         if (i == 17 && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //noBoundedEvidence
         if (i == 19 && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //fl
      }//for_i
   }//if 

// Check if annotation file is at least version 1.8 with new column probeset_type
   if (hasColumn[kNProbesetCols-1] == kNProbesetCols-1) {hasPST = kTRUE;} //probeset_type

// Get number of probesets
   // get current streamposition for rewinding input
   position = input.tellg();
   while (1) {
      input.getline(nextline, kAnnBuf, delim);
      if (input.eof()) break;
      size++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   Number of probesets is <" << size << ">." << endl;
   }//if
   // set number of annotated probesets (less than probeset_count in ImportScheme())
   fNProbesets = size;

   // reset input file to position
   input.clear();  //clear all flags
   input.seekg(position, ios::beg);

// Initialize memory for local arrays
   if (!(probesetid = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(seqname    = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(strand     = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(start      = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(stop       = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(nprobes    = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(geneid     = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(exonid     = new (nothrow) Int_t[size]))   {err = errInitMemory; goto cleanup;}
   if (!(crosshtype = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(full       = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(geneaccess = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(genesymbol = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(mrnaaccess = new (nothrow) TString[size])) {err = errInitMemory; goto cleanup;}
   if (!(levelid    = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(bounded    = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}
   if (!(psettype   = new (nothrow) Short_t[size])) {err = errInitMemory; goto cleanup;}

// Initialize memory for sorting arrays
   if (!(index = new (nothrow) Int_t[size]))        {err = errInitMemory; goto cleanup;}
   if (!(value = new (nothrow) Long64_t[size]))     {err = errInitMemory; goto cleanup;}

   while (input.good()) {
      input.getline(nextline, kAnnBuf, delim);
      if (input.eof()) break;
      if (input.fail()) {
         cout << "Error: Failed reading line:" << endl;
         cout << &nextline[0] << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // replace all "\",\"" with tab "\t" and remove "\""
      str = TString(&nextline[0]);
      // first need to replace "" with NA
      str.ReplaceAll("\"\"", "\"NA\"");
      // replace tab with space to eliminate wrong tabs in Affymetrix annotation files
      str.ReplaceAll(tab, " ");
      // replace csv with tab
      str.ReplaceAll(csv, tab);
      // remove all "\"" from line
      str.ReplaceAll("\"", "");

      // check number of separators
      if (numsep != NumSeparators(str, tab)) {
         cout << "Error: Wrong number of separators in line:" << endl;
         cout << &nextline[0] << endl;
//??         return errReadingInput;
      }//if

      // import fields
      probesetid[idx] = atoi(strtok((char*)str.Data(), tab));
      seqname[idx]    = strtok(NULL, tab);
      strand[idx]     = strtok(NULL, tab);
      start[idx]      = atoi(strtok(NULL, tab));
      stop[idx]       = atoi(strtok(NULL, tab));
      nprobes[idx]    = atoi(strtok(NULL, tab));
      geneid[idx]     = atoi(strtok(NULL, tab));
      exonid[idx]     = atoi(strtok(NULL, tab));
      dummy           = strtok(NULL, tab);
      assigngene      = strtok(NULL, tab);
      assignmrna      = strtok(NULL, tab);
      crosshtype[idx] = atoi(strtok(NULL, tab));
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      level           = strtok(NULL, tab);
      bound           = atoi(strtok(NULL, tab));
      nobound         = atoi(strtok(NULL, tab));
      dummy           = strtok(NULL, tab);
      full[idx]       = atoi(strtok(NULL, tab));
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      dummy           = strtok(NULL, tab);
      probesettype    = hasPST ? strtok(NULL, tab) : "NA";

      // get gene accession and symbol
      geneaccess[idx] = strtok((char*)assigngene.Data(), "/");
      if (strcmp(geneaccess[idx].Data(),"---") == 0) {
         geneaccess[idx] = "NA";
         genesymbol[idx] = "NA";
      } else {
         geneaccess[idx] = RemoveEnds(geneaccess[idx]);
         genesymbol[idx] = RemoveEnds(strtok(0, "/"));

         if (strcmp(genesymbol[idx].Data(),"") == 0) {genesymbol[idx] = "NA";}
      }//if

      // get bounded as combination of bounded and noBoundedEvidence
      bounded[idx] = nobound ? (bound ? eBND_T_NOBND_T : eBND_F_NOBND_T)
                             : (bound ? eBND_T_NOBND_F : eBND_F_NOBND_F);

      // convert probesettype to probesettype_id
      psettype[idx] = ProbesetType(probesettype);

      // convert level to level_id
      levelid[idx] = this->ProbesetLevel(level, psettype[idx]);

      // get mrna accession for genes or for control->affx 
      if (psettype[idx] != eCONTROLAFFX) {
         mrnaaccess[idx] = RemoveEnds(strtok((char*)assignmrna.Data(), "/"));
         if (strcmp(mrnaaccess[idx].Data(),"") == 0) {mrnaaccess[idx] = "NA";}
      } else {
         mrnaaccess[idx] = RemoveEnds(assignmrna);

         // replace gene accession and symbol for control->affx
         geneaccess[idx] = mrnaaccess[idx];

         // add control->affx names to list fAffxNames
         XIdxString *idxstr = new XIdxString(probesetid[idx], mrnaaccess[idx]);
         fAffxNames->Add(idxstr);  //idxstr will be deleted by fAffxNames

         fNAffx++;
      }//if

      // fill sort values
      majorv = (Long64_t)geneid[idx]; //sort for transcript_cluster_id first
      minorv = (Long64_t)start[idx];  //then sort for start position
      value[idx]  = majorv << 31;
      value[idx] += minorv;
// TO DO:
// ev. sort for "geneid" then "start" then "probesetid"
// see e.g. "HuGene-1_1-st-v1.na30.hg19.probeset.csv"

      if (XManager::fgVerbose && idx%10000 == 0) {
         cout << "   <" << idx + 1 << "> records read...\r" << flush;
      }//if
      idx++;
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << idx << "> records read...Finished" << endl;
   }//if

// Sort for geneid first and then for start (assuming transcipt on one chromosome)
   TMath::Sort(size, value, index, 0);

// Fill probeset/exon annotation tree with control->affx 
   if ((hasPST == kTRUE) && (fNAffx == 0)) {
      for (Int_t i=0; i<size; i++) {
         if (psettype[i] != eCONTROLAFFX) continue;

         // fill probeset annotation tree
         psann->SetUnitID(fNAffx);
         psann->SetProbesetID(probesetid[i]);
         psann->SetName(geneaccess[i]);
         psann->SetSymbol(genesymbol[i]);
         psann->SetAccession(mrnaaccess[i]);
         psann->SetChromosome(seqname[i]);
         psann->SetStrand((strand[i].Data())[0]);
         psann->SetStart(start[i]);
         psann->SetStop(stop[i]);
         psann->SetNumProbes(nprobes[i]);
         psann->SetTranscriptID(geneid[i]);
         psann->SetExonID(exonid[i]);
         psann->SetCrossHybType(crosshtype[i]);
         psann->SetLevelID(levelid[i]);
         psann->SetBounded(bounded[i]);
         psann->SetFullFlag(full[i]);
         psann->SetProbesetType(psettype[i]);
         psanntree->Fill();

         // fill exon annotation tree
         exann->SetUnitID(fNAffx);
         exann->SetName(geneaccess[i]);
         exann->SetSymbol(genesymbol[i]);
         exann->SetAccession(mrnaaccess[i]);
         exann->SetChromosome(seqname[i]);
         exann->SetStrand((strand[i].Data())[0]);
         exann->SetStart(start[i]);
         exann->SetStop(stop[i]);
         exann->SetTranscriptID(geneid[i]);
         exann->SetExonID(exonid[i]);
         exanntree->Fill();

         fNAffx++;
      }//for_i
   } else if ((fNAffx > 0) && (fAffxNames->GetSize() == fNAffx)) {
      for (Int_t i=0; i<fNAffx; i++) {
         XIdxString *idxstr   = (XIdxString*)(fAffxNames->At(i));
         TString     affxname = idxstr->GetString();
         Int_t       affxidx  = idxstr->GetIndex();

         // affxid as string
         str.Form("%d", idxstr->GetIndex());

         // fill probeset annotation tree
         psann->SetUnitID(i);
         psann->SetProbesetID(affxidx);
         psann->SetName(affxname);
         psann->SetSymbol("NA");
         psann->SetAccession(str);
         psann->SetChromosome("NA");
         psann->SetStrand('?');
         psann->SetStart(0);
         psann->SetStop(0);
         psann->SetNumProbes(0);
//TO DO:         psann->SetNumProbes(???);
         psann->SetTranscriptID(0);
//??         psann->SetTranscriptID(affxidx); //may interfere with transcript_cluster_id???
         psann->SetExonID(0);
         psann->SetCrossHybType(0);
         psann->SetLevelID(eNOLEVEL);
         psann->SetBounded(0);
         psann->SetFullFlag(0);
         psann->SetProbesetType(eCONTROLAFFX);
         psanntree->Fill();

         // fill exon annotation tree
         exann->SetUnitID(i);
         exann->SetName(affxname);
         exann->SetSymbol(str);
         exann->SetSymbol("NA");
         exann->SetAccession(str);
         exann->SetChromosome("NA");
         exann->SetStrand('?');
         exann->SetStart(0);
         exann->SetStop(0);
         exann->SetTranscriptID(0);
//??         exann->SetTranscriptID(affxidx); //may interfere with transcript_cluster_id???
         exann->SetExonID(0);
         exanntree->Fill();
      }//for_i
   } else {
      cerr << "Error: Problem with number of AFFX control <" << fNAffx << ">. "
           << "Aborting import of annotation file." << endl;
      err = errReadingInput;
      goto cleanup;
   }//if

// Fill probeset/exon annotation tree with sorted "main" annotation data
   idx      = 0; //use as switch for chromosome start
   numpsets = fNAffx;  //start at fNAffx
   numexons = fNAffx;  //start at fNAffx
   for (Int_t i=0; i<size; i++) {
      Int_t k = index[i];
      if (hasPST && !((psettype[k] == eMAIN) || 
                      (psettype[k] == eEXON) ||
                      (psettype[k] == eINTRON))) continue;

      // fill probeset annotation tree with annotation data
      psann->SetUnitID(numpsets++);
      psann->SetProbesetID(probesetid[k]);
      psann->SetName(geneaccess[k]);
      psann->SetSymbol(genesymbol[k]);
      psann->SetAccession(mrnaaccess[k]);
      psann->SetChromosome(seqname[k]);
      psann->SetStrand((strand[k].Data())[0]);
      psann->SetStart(start[k]);
      psann->SetStop(stop[k]);
      psann->SetNumProbes(nprobes[k]);
      psann->SetTranscriptID(geneid[k]);
      psann->SetExonID(exonid[k]);
      psann->SetCrossHybType(crosshtype[k]);
      psann->SetLevelID(levelid[k]);
      psann->SetBounded(bounded[k]);
      psann->SetFullFlag(full[k]);
      psann->SetProbesetType(psettype[k]);
      psanntree->Fill();

      // start position for current exon with exonID
      exonID = exonid[k];
      if (idx == 0) {
         exonstart = start[k];
         idx = 1;
      }//if

      // fill exon annotation tree with annotation data
      Int_t k1 = index[i+1];
      if (exonID != exonid[k1]) {
         // stop position for current exon with exonID
         exonstop = stop[k];

         // fill exon annotation tree with annotation data
         exann->SetUnitID(numexons++);
         exann->SetName(geneaccess[k]);
         exann->SetSymbol(genesymbol[k]);
         exann->SetAccession(mrnaaccess[k]);
         exann->SetChromosome(seqname[k]);
         exann->SetStrand((strand[k].Data())[0]);
         exann->SetStart(exonstart);
         exann->SetStop(exonstop);
         exann->SetTranscriptID(geneid[k]);
         exann->SetExonID(exonid[k]);
         exanntree->Fill();

         exonID = exonid[k1];
         idx    = 0; //reset idx to get start position of next exon
      }//if

      if (XManager::fgVerbose && numpsets%10000 == 0) {
         cout << "   <" << numpsets << "> records imported...\r" << flush;
      }//if
   }//for_i

// Fill probeset/exon annotation tree with sorted unmapped annotation data
   for (Int_t i=0; i<size; i++) {
      Int_t k = index[i];
      if (!hasPST) break;  //???
      if (hasPST && psettype[k] != eUNMAPPED) continue;

      // fill probeset annotation tree with annotation data
      psann->SetUnitID(numpsets++);
      psann->SetProbesetID(probesetid[k]);
      psann->SetName(geneaccess[k]);
      psann->SetSymbol(genesymbol[k]);
      psann->SetAccession(mrnaaccess[k]);
      psann->SetChromosome(seqname[k]);
      psann->SetStrand((strand[k].Data())[0]);
      psann->SetStart(start[k]);
      psann->SetStop(stop[k]);
      psann->SetNumProbes(nprobes[k]);
      psann->SetTranscriptID(geneid[k]);
      psann->SetExonID(exonid[k]);
      psann->SetCrossHybType(crosshtype[k]);
      psann->SetLevelID(levelid[k]);
      psann->SetBounded(bounded[k]);
      psann->SetFullFlag(full[k]);
      psann->SetProbesetType(psettype[k]);
      psanntree->Fill();

      // start position for current exon with exonID
      exonID = exonid[k];
      if (idx == 0) {
         exonstart = start[k];
         idx = 1;
      }//if

      // fill exon annotation tree with annotation data
      Int_t k1 = index[i+1];
      if (exonID != exonid[k1]) {
         // stop position for current exon with exonID
         exonstop = stop[k];

         // fill exon annotation tree with annotation data
         exann->SetUnitID(numexons++);
         exann->SetName(geneaccess[k]);
         exann->SetSymbol(genesymbol[k]);
         exann->SetAccession(mrnaaccess[k]);
         exann->SetChromosome(seqname[k]);
         exann->SetStrand((strand[k].Data())[0]);
         exann->SetStart(exonstart);
         exann->SetStop(exonstop);
         exann->SetTranscriptID(geneid[k]);
         exann->SetExonID(exonid[k]);
         exanntree->Fill();

         exonID = exonid[k1];
         idx    = 0; //reset idx to get start position of next exon
      }//if

      if (XManager::fgVerbose && numpsets%10000 == 0) {
         cout << "   <" << numpsets << "> records imported...\r" << flush;
      }//if
   }//for_i
   if (XManager::fgVerbose) {
      cout << "   <" << numpsets << "> records imported...Finished" << endl;
   }//if

// Set number of exons
   cout << "   <" << numexons - fNAffx << "> exon annotations imported." << endl;
   fNExons = numexons - fNAffx;

//TO DO:
// Add tree info to tree
//   AddAnnotationTreeInfo(exanntree, fAnnTreeName);

// Write exon annotation tree to file
   if ((err = WriteTree(exanntree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(exanntree->GetName(), 0);
   }//if
   //delete tree from memory
   exanntree->Delete("");
   exanntree = 0;

//TO DO:
// Add tree info to tree
//   AddAnnotationTreeInfo(psanntree, fPbsTreeName);

// Write probeset annotation tree to file
   if ((err = WriteTree(psanntree, TObject::kOverwrite)) == errNoErr) {
      // add tree header to list
      AddTreeHeader(psanntree->GetName(), 0);
   }//if
   //delete tree from memory
   psanntree->Delete("");
   psanntree = 0;

// Clean up
cleanup:
   SafeDelete(psann);
   SafeDelete(exann);

   if (value)      {delete [] value;      value      = 0;}
   if (index)      {delete [] index;      index      = 0;}
   if (psettype)   {delete [] psettype;   psettype   = 0;}
   if (bounded)    {delete [] bounded;    bounded    = 0;}
   if (levelid)    {delete [] levelid;    levelid    = 0;}
   if (mrnaaccess) {delete [] mrnaaccess; mrnaaccess = 0;}
   if (genesymbol) {delete [] genesymbol; genesymbol = 0;}
   if (geneaccess) {delete [] geneaccess; geneaccess = 0;}
   if (full)       {delete [] full;       full       = 0;}
   if (crosshtype) {delete [] crosshtype; crosshtype = 0;}
   if (exonid)     {delete [] exonid;     exonid     = 0;}
   if (geneid)     {delete [] geneid;     geneid     = 0;}
   if (nprobes)    {delete [] nprobes;    nprobes    = 0;}
   if (stop)       {delete [] stop;       stop       = 0;}
   if (start)      {delete [] start;      start      = 0;}
   if (strand)     {delete [] strand;     strand     = 0;}
   if (seqname)    {delete [] seqname;    seqname    = 0;}
   if (probesetid) {delete [] probesetid; probesetid = 0;}
   if (hasColumn)  {delete [] hasColumn;  hasColumn  = 0;}

   return err;
}//ImportProbesetAnnotation

//______________________________________________________________________________
Int_t XExonChip::ImportControlAnnotation(ifstream &input, Option_t * /*option*/, 
                 const char *sep, char delim, Int_t /*split*/)
{
   // Import control->affx annotation from infile and store in list fAffxNames
   // Note: The AFFX controls are NOT stored in a control tree, but in the
   //       list to be stored in the other annotation trees!
   if(kCS) cout << "------XExonChip::ImportControlAnnotation------" << endl;

   Int_t err = errNoErr;
   Int_t probesetid = 0;

   std::string nextline;  //to read next line

   TString lib_set_name, lib_set_version;
   TString str, dummy, affxname;

// Check for line "#%lib-set-name"
   while (nextline.compare(0, 14, "#%lib_set_name") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_name = strtok((&((char*)nextline.c_str())[15]), sep);

// Check for line "#%lib-set-version"
   while (nextline.compare(0, 17, "#%lib_set_version") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   lib_set_version = strtok((&((char*)nextline.c_str())[18]), sep);

// Check for line ""probeset_id""
   while (nextline.compare(0, 11, "probeset_id") != 0) {
      std::getline(input, nextline, delim);
      if (input.eof()) return errPrematureEOF;
   }//while
   Int_t numsep = NumSeparators(nextline.c_str(), sep);

// Create array hasColumn 
   Int_t *hasColumn = 0;
   if (!(hasColumn = new (nothrow) Int_t[kNControlCols])) return errInitMemory;
   for (Int_t i=0; i<kNControlCols; i++) {
      hasColumn[i] = 0;
   }//for_i

// Check column headers from header line, hasColumn = 1 if column is present
   str = RemoveEnds(nextline.c_str());  //need to remove <cr>
   err = CheckHeaderOrder(str, kControlHeader, kNControlCols, hasColumn, sep);
   if (err > 0) {
      cout << "Note: The following header columns are missing or in wrong order: " << endl;
      for (Int_t i=0; i<kNControlCols; i++) {
         if (hasColumn[i] != i) cout << "   <" << kControlHeader[i] << ">" << endl;
         if (i == 0  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //probeset_id
         if (i == 2  && hasColumn[i] != i) {err = errMissingColumn; goto cleanup;} //probeset_name
      }//for_i
      err = errNoErr; //reset err!
   }//if 

// check number of separators
   if (fNAffx > 0) {
      cout << "Warning: Number of AFFX controls is already <" << fNAffx << ">."
           << " Thus control file will not be imported!" << endl;
      return errNoErr;
   }//if

   while (input.good()) {
      std::getline(input, nextline, delim);
      if (input.eof()) break;
      if (input.fail()) {
         cout << "Error: Failed reading line:" << endl;
         cout << nextline << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      str = TString(nextline);

      // check number of separators
      if (numsep != NumSeparators(str, sep)) {
         cout << "Error: Wrong number of separators in line:" << endl;
         cout << str.Data() << endl;
         err = errReadingInput;
         goto cleanup;
      }//if

      // import fields
      probesetid = atoi(strtok((char*)str.Data(), sep));
      dummy      = strtok(NULL, sep);
      affxname   = RemoveEnds(strtok(NULL, sep));  //need to remove <cr>

      // add control->affx names to list fAffxNames
      if (strncmp("AFFX", affxname, 4) == 0) {
         XIdxString *idxstr = new XIdxString(probesetid, affxname);
         fAffxNames->Add(idxstr);

         fNAffx++;
      }//if
   }//while
   if (XManager::fgVerbose) {
      cout << "   <" << fNAffx << "> AFFX controls imported." << endl;
      cout << "   Note: No tree will be created but controls are saved for other annotation trees!"
           << endl;
   }//if

// Clean up
cleanup:
   if (hasColumn)  {delete [] hasColumn;  hasColumn  = 0;}

   return err;
}//ImportControlAnnotation

//______________________________________________________________________________
Int_t XExonChip::LayoutToX(Int_t index)
{
   // Convert index to X coordinate if fSequential==FALSE
   if(kCSa) cout << "------XExonChip::LayoutToX------" << endl;

//to do

   return 0;
}//LayoutToX

//______________________________________________________________________________
Int_t XExonChip::LayoutToY(Int_t index)
{
   // Convert index to Y coordinate if fSequential==FALSE
   if(kCSa) cout << "------XExonChip::LayoutToY------" << endl;

//to do

   return 0;
}//LayoutToY

//______________________________________________________________________________
Int_t XExonChip::ControlchipType(const char *type)
{
   // Convert controlchip type to id EControlchipType
   if(kCSa) cout << "------XExonChip::ControlchipType------" << endl;

   if (strstr(type, "board:at") != 0) {
      return eJUMBOAT;
   } else if (strstr(type, "board:st")   != 0) {
      return eJUMBOST;
   } else if (strstr(type, "thermo:at")  != 0) {
      return eTHERMOAT;
   } else if (strstr(type, "thermo:st")  != 0) {
      return eTHERMOST;
   } else if (strstr(type, "trigrid:at") != 0) {
      return eTRIGRIDAT;
   } else if (strstr(type, "trigrid:st") != 0) {
      return eTRIGRIDST;
   } else if (strstr(type, "generic:at") != 0) {
      return eGENERICAT;
   } else if (strstr(type, "generic:st") != 0) {
      return eGENERICST;
   } else if (strstr(type, "blank")      != 0) {
      return eBLANK;
   } else if (strstr(type, "pm:st")      != 0) {
      return eCTRLPMST;
   }//if

   return eUNKNOWN;
}//ControlchipType

//______________________________________________________________________________
TString XExonChip::LevelID2Level(Short_t id)
{
   // Convert level_id to probeset level
   if(kCSa) cout << "------XExonChip::LevelID2Level------" << endl;

   switch (id) {
      case ePARACORE:
         return TString("core");
         break;

      case ePARAEXTENDED:
         return TString("extended");
         break;

      case ePARAFULL:
         return TString("full");
         break;

      case eAMBIGUOUS:
         return TString("ambiguous");
         break;

      case eFREE:
         return TString("free");
         break;

      default: //eNOLEVEL
         return TString("nolevel");
         break;
   }//switch

   return TString("nolevel");
}//LevelID2Level

//______________________________________________________________________________
Int_t XExonChip::ProbesetLevel(const char *level, Short_t type)
{
   // Convert probeset level to level_id ELevel
   if(kCSa) cout << "------XExonChip::ProbesetLevel------" << endl;

   if (strcmp(level,"core") == 0) {
      return ePARACORE;
   } else if (strcmp(level,"extended") == 0) {
      return ePARAEXTENDED;
   } else if (strcmp(level,"full") == 0) {
      return ePARAFULL;
   } else if (strcmp(level,"ambiguous") == 0) {
      return eAMBIGUOUS;
   } else if (strcmp(level,"free") == 0) {
      return eFREE;
   } else if (strcmp(level,"NA") == 0) {
      // for whole genome arrays with annotation <= na30
      return ePARACORE;
   } else if (strcmp(level,"---") == 0 && type == eMAIN) {
      // for whole genome arrays with annotation >= na31
      return ePARACORE;
   }//if

   return eNOLEVEL;
}//ProbesetLevel

//______________________________________________________________________________
Int_t XExonChip::SchemeMask(Short_t xtype, Short_t level, Short_t bound)
{
   // Convert cross-hybridization probeset_type "xtype" with probeset level "level"
   // and bounded probeset group "bound" to mask:
   if(kCSa) cout << "------XExonChip::SchemeMask(main)------" << endl;

   if ((xtype == 1) && (bound == 0)) { //(probeset_type==unique) && (bounded==0)
      // mask for meta-probesets
      if (level == ePARACORE) {
         return eMETACORE;
      } else if (level == ePARAEXTENDED) {
         return eMETAEXTENDED;
      } else if (level == ePARAFULL) {
         return eMETAFULL;
      } else if (level == eAMBIGUOUS) {
         return eAMBIGUOUS;
      } else if (level == eFREE) {
         return eFREE;
      } else if (level == eNOLEVEL) {
         return eNOLEVEL;
      }//if
   } else {
      if (level == ePARACORE) {
         return ePARACORE;
      } else if (level == ePARAEXTENDED) {
         return ePARAEXTENDED;
      } else if (level == ePARAFULL) {
         return ePARAFULL;
      } else if (level == eAMBIGUOUS) {
         return eAMBIGUOUS;
      } else if (level == eFREE) {
         return eFREE;
      } else if (level == eNOLEVEL) {
         return eNOLEVEL;
      }//if
   }//if

   return eUNKNOWNTYPE;
}//SchemeMask

