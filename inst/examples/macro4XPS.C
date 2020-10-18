/******************************************************************************
* Macro to test XPS command line                                              *
*                                                                             *
* Author: Christian Stratowa, Vienna, Austria    .                            *
* Created: 15 Mar 2003                             Last modified: 25 Nov 2007 *
******************************************************************************/
//
///////////////////////////
//
// in new root session(s):
//   0.step: initialize
//     > .L macro4XPS.C 
//     > Init("/Volumes/CoreData/ROOT/rootcode/xps-0.4.1/src/xps.so") 
//   1.step: create new root scheme file and add schemes
//     Test3:
//     > ImportSchemeTest3("SchemeTest3","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1")
//     U133P2:
//     > ImportSchemeU133P2("SchemeU133P2","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1")
//     HuGene:
//     > ImportSchemeHuGene("HumanGeneScheme_v1","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1")
//     HuExon:
//     > ImportSchemeHuExon("HumanExonScheme_r2","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1")
//   2.step: create new root data file and add raw data (hybridizations)
//     Test3:
//     > ImportDataTest3("DataTest3","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1")
//     U133P2:
//     > ImportDataU133P2("U133P2Tissues","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1")
//     HuGene:
//     > ImportDataHuGene("HuGeneTissues","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1")
//     HuExon:
//     > ImportDataHuExon("HuExonTissues","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1")
//   3.step: create new root expression file and summarize raw data
//     Test3:
//     > RMATest3("Test3_rma","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1")
//     > MAS5Test3("Test3_mas5","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1")
//     U133P2:
//     > RMAU133P2("U133P2_rma","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1")
//     > MAS5U133P2("U133P2_mas5","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1")
//     HuGene:
//     > RMAHuGene("HuGene_rma",9216,"/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1")
//     > MAS5HuGene("HuGene_mas5",9216,"/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1")
//     HuExon:
//     > RMAHuExon("HuExon_rma_tr9216",9216,"log2")
//     > RMAHuExon("HuExon_rma_ps9216",9216,"probeset:log2")
//     > MAS5HuExon("HuExon_mas5_tr9216",9216,"log2")
//     > MAS5HuExon("HuExon_mas5_ps9216",9216,"probeset:log2")
//   4.step: create new root call file and compute detection call
//     Test3:
//     > MAS5CallTest3("Test3_dc5","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1")
//     > DABGCallTest3("Test3_dabg","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1")
//     U133P2:
//     > MAS5CallU133P2("U133P2_dc5","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1")
//     HuGene:
//     > MAS5CallHuGene("HuGene_dc5",9216,"/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1")
//     > DABGCallHuGene("HuGene_dabg",9216,"/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1")
//     HuExon:
//     > MAS5CallHuExon("HuExon_dc5_tr9216",9216,"adjusted")
//     > MAS5CallHuExon("HuExon_dc5_ps9216",9216,"probeset:adjusted")
//     > DABGCallHuExon("HuExon_dabg_tr9216",9216,"raw")
//     > DABGCallHuExon("HuExon_dabg_ps9216",9216,"probeset:raw")
//   5.step: export trees
//     Test3:
//     > ExportScheme("/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","Test3.Test3.scm","*")
//     > ExportScheme("/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","Test3.Test3.idx","*")
//     > ExportScheme("/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","Test3.Test3.ann","*")
//     > ExportScheme("/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","Test3.Test3.prb","*")
//     > ExportDataset("/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root","DataSet/*.cel","*","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.txt")
//     > ExportPreprocess("GeneChip","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_rma.root","RMASet.*.mdp","*","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_mdp.txt")
//     > ExportPreprocess("GeneChip","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_rma.root","RMASet.*.mdp","fLevel","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_mdp_level.txt")
//     > ExportPreprocess("GeneChip","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_mas5.root","MAS5Set.*.tbw","*","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_tbw.txt")
//     > ExportPreprocess("GeneChip","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_dc5.root","MAS5CallSet.*.dc5","*","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_dc5.txt")
//     > ExportPreprocess("GeneChip","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_dabg.root","DABGCallSet.*.dab","*","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/Test3_dab.txt")
//     U133P2:
//     > ExportScheme("/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root","HG-U133_Plus_2.HG-U133_Plus_2.scm","*")
//     > ExportScheme("/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root","HG-U133_Plus_2.HG-U133_Plus_2.idx","*")
//     > ExportPreprocess("GeneChip","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2_rma.root","RMASet.*.mdp","*","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2_mdp.txt")
//     > ExportPreprocess("GeneChip","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2_mas5.root","MAS5Set.*.tbw","*","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2_tbw.txt")
//     > ExportPreprocess("GeneChip","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2_dc5.root","MAS5CallSet.*.dc5","*","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2_dc5.txt")
//     HuGene:
//     > ExportScheme("/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","HuGene-1_0-st-v1.HuGene-1_0-st-v1.ann","*")
//     > ExportDataset("/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root","DataSet/*.cel","*","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.txt")
//     > ExportPreprocess("GenomeChip","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_rma.root","RMASet.*.mdp","*")
//     > ExportPreprocess("GenomeChip","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_rma.root","RMASet.*.mdp","*","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_mdp.txt")
//     > ExportPreprocess("GenomeChip","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_rma.root","RMASet.*.mdp","fLevel","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_mdp_level.txt")
//     > ExportPreprocess("GenomeChip","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_mas5.root","MAS5Set.*.tbw","*","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_tbw.txt")
//     > ExportPreprocess("GenomeChip","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_dc5.root","MAS5CallSet.*.dc5","*","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_dc5.txt")
//     > ExportPreprocess("GenomeChip","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_dabg.root","DABGCallSet.*.dab","*","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGene_dab.txt")
//     HuExon:
//     > ExportScheme("/Volumes/GigaDrive/ROOT/rootdata/huexon/HumanExonScheme_r2.root","HuEx-1_0-st-v2.HuEx-1_0-st-v2.scm","*")
//     > ExportScheme("/Volumes/GigaDrive/ROOT/rootdata/huexon/HumanExonScheme_r2.root","HuEx-1_0-st-v2.HuEx-1_0-st-v2.idx","*")
//     > ExportPreprocess("ExonChip","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_rma_tr9216.root","RMASet.*.mdp","*","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_mdp_tr9216.txt")
//     > ExportPreprocess("ExonChip","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_rma_tr9216.root","RMASet.*.mdp","fLevel","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_mdp_tr9216_level.txt")
//     > ExportPreprocess("ExonChip","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_mas5_tr9216.root","MAS5Set.*.tbw","*","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_tbw_tr9216.txt")
//     > ExportPreprocess("ExonChip","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_dc5_tr9216.root","MAS5CallSet.*.dc5","*","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_dc5_tr9216.txt")
//     > ExportPreprocess("ExonChip","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_dabg_tr9216.root","DABGCallSet.*.dab","*","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExon_dab_tr9216.txt")
//
///////////////////////////

//______________________________________________________________________________
void Init(const char *libxps = "/Volumes/CoreData/ROOT/rootcode/xps-0.4.1/src/xps.so")
{
// load libraries
   gSystem->Load("libGui.so");
   gSystem->Load(libxps);
}//Init


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Create ROOT scheme files                                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void ImportSchemeTest3(const char *filename, const char *filedir)
{
// Import Affymetrix chip definition file Test3.CDF into root file

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager");

// create new root schemes file
   manager->New(filename, filedir, "GeneChip", "Schemes");

// store chip definitions,  probe sequences and annotations
   // Test3:
   manager->NewScheme("Test3","/Volumes/CoreData/Firma/Affy/libraryfiles/Test3.CDF");
   manager->NewScheme("Test3","/Volumes/CoreData/Firma/Affy/Annotation/Test3_probe.tab","probe");
   manager->NewAnnotation("Test3","/Volumes/CoreData/Firma/Affy/Annotation/Test3.na21.annot.csv");

// cleanup
   manager->Close();
   delete manager;
}//ImportSchemeTest3

//______________________________________________________________________________
void ImportSchemeU133P2(const char *filename, const char *filedir, Int_t verbose = 1)
{
// Import Affymetrix chip definition file HG-U133_Plus_2.CDF into root file

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager", "", verbose);

// create new root schemes file
   manager->New(filename, filedir, "GeneChip", "Schemes");

// store chip definitions,  probe sequences and annotations
   // HG-U133_Plus_2:
   manager->NewScheme("HG-U133_Plus_2","/Volumes/GigaDrive/Affy/libraryfiles/HG-U133_Plus_2.cdf");
   manager->NewScheme("HG-U133_Plus_2","/Volumes/GigaDrive/Affy/libraryfiles/HG-U133-PLUS_probe.tab","probe");
   manager->NewAnnotation("HG-U133_Plus_2","/Volumes/GigaDrive/Affy/Annotation/Version07Jul/HG-U133_Plus_2.na23.annot.csv");

// cleanup
   manager->Close();
   delete manager;
}//ImportSchemeU133P2

//______________________________________________________________________________
void ImportSchemeHuGene(const char *filename, const char *filedir)
{
// Import Affymetrix chip definition files *.CLF and *.PGF files into XPS

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager");

// create new root schemes file
   manager->New(filename, filedir, "GenomeChip", "Schemes");

// store chip definitions,  probe sequences and annotations
// note: setname must be identical to CLF-name (since CEL-files contain CLF name!)
   // HuGene:
   manager->NewScheme("HuGene-1_0-st-v1","/Volumes/GigaDrive/Affy/libraryfiles/HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf","layout");
   manager->NewAnnotation("HuGene-1_0-st-v1","/Volumes/GigaDrive/Affy/Annotation/Version07Jul/HuGene-1_0-st-v1.na23.hg18.transcript.csv","transcript");
   manager->NewScheme("HuGene-1_0-st-v1","/Volumes/GigaDrive/Affy/libraryfiles/HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf","scheme");

// cleanup
   manager->Close();
   delete manager;
}//ImportSchemeHuGene

//______________________________________________________________________________
void ImportSchemeHuExon(const char *filename, const char *filedir)
{
// Import Affymetrix chip definition files *.CLF and *.PGF files into XPS

// create new scheme manager
   XSchemeManager *manager = new XSchemeManager("SchemeManager");

// create new root schemes file
   manager->New(filename, filedir, "ExonChip", "Schemes");

// store chip definitions,  probe sequences and annotations
// note: setname must be identical to PGF-name (since CEL-files contain PGF name!)
   // Exon-v2:
   manager->NewScheme("HuEx-1_0-st-v2","/Volumes/GigaDrive/Affy/libraryfiles/HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf","layout");
   manager->NewAnnotation("HuEx-1_0-st-v2","/Volumes/GigaDrive/Affy/Annotation/Version07Jul/HuEx-1_0-st-v2.na23.hg18.probeset.csv","exon");
   manager->NewScheme("HuEx-1_0-st-v2","/Volumes/GigaDrive/Affy/libraryfiles/HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf","scheme");
   manager->NewAnnotation("HuEx-1_0-st-v2","/Volumes/GigaDrive/Affy/Annotation/Version07Jul/HuEx-1_0-st-v2.na23.hg18.transcript.csv","transcript");

// cleanup
   manager->Close();
   delete manager;
}//ImportSchemeHuExon

//______________________________________________________________________________
void ExportScheme(const char *filename, const char *treename, const char *varlist="*")
{
// Export tree "treename" from file "filename"

// create new scheme manager
   XSchemeManager *manager = 0;
   manager = new XSchemeManager("SchemeManager");

// open root scheme file
   manager->Open(filename);

// export tree
   manager->Export(treename, varlist);

// cleanup
   manager->Close();
   delete manager;
}//ExportScheme


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Create ROOT raw data files                                                 //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void ImportDataTest3(const char *filename   = "DataTest3",
                     const char *filedir    = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1",
                     const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root")
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager");

// initialize chip type and variable list
   manager->Initialize("GeneChip");
   manager->InitInput("Test3","cel","MEAN/D:STDV/D:NPIXELS/I","RawData");

// create new root data file 
   manager->New(filename, filedir, "GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// store *.CEL data as tree in data file
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestA1.CEL","TestA1");
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestA2.CEL","TestA2");
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestB1.CEL","TestB1");
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestB2.CEL","TestB2");

// cleanup
   manager->Close();
   delete manager;
}//ImportDataTest3

//______________________________________________________________________________
void ImportDataTest3A(const char *filename   = "DataTest3",
                      const char *filedir    = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.3.12",
                      const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.3.12/SchemeTest3.root")
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager");

// initialize chip type and variable list
   manager->Initialize("GeneChip");
   manager->InitInput("Test3","cel","MEAN/D:STDV/D:NPIXELS/I","RawData");

// create new root data file 
   manager->New(filename, filedir, "GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// optionally add database info
   manager->BeginTransaction();
//?   manager->LoginInfo("strato", "password");
   manager->ProjectInfo("TestProject",20060126, "Project Type","use Test3 data for testing","my comment");
   manager->AuthorInfo("Stratowa","Christian","Project Leader","Company","Dept","cstrato.at.aon.at","++43-1-1234","my comment");
   manager->DatasetInfo("DataSet","MC","Tissue","Stratowa",20060126,"description","my comment");
   manager->SourceInfo("Unknown","source type","Homo sapiens","caucasian","description","my comment");
   // demo only: use only one of SampleInfo, CellLineInfo, PrimaryCellInfo, TissueInfo, BiopsyInfo!!!
   manager->SampleInfo("Yeast","sample type","NA","my pheno","my genotype","RNA extraction",0,"","",0.0,"", "my comment");
   manager->CellLineInfo("HeLa-S3","cell type","HeLa","ATCC-12.3","pCSV transfected","female","my pheno","my genotype","RNA extraction",0,"","",0.0,"", "my comment");
   manager->PrimaryCellInfo("Mel31","primary cell",20071123,"extracted from patient","male","my pheno","my genotype","RNA extraction",1,"NMRI","female",7.0,"months", "my comment");
   manager->TissueInfo("Liver","tissue type","adult","morphology","carcinoma","Grade 3",34.6,"years","dead","male","my pheno","my genotype","RNA extraction",0,"","",0.0,"", "my comment");
   manager->BiopsyInfo("Breast","needle biopsy","morphology","carcinoma","Grade 3",45.5,"years","alive","female","my pheno","my genotype","RNA extraction",0,"","",0.0,"", "my comment");
   manager->ArrayInfo("Test3","GeneChip","description","my comment");
   manager->HybridizationInfo("TestA1","hyb type","TestA1.CEL",20071117,"my prep1","standard protocol","A1",1,"my comment");
   manager->HybridizationInfo("TestA2","hyb type","TestA2.CEL",20071117,"my prep2","standard protocol","A2",1,"my comment");
   manager->TreatmentInfo("TestA1","DMSO",4.3,"mM",1.0,"hours","intravenous","my comment");
   manager->TreatmentInfo("TestA2","DMSO",4.3,"mM",8.0,"hours","intravenous","my comment");
   manager->CommitTransaction();

// store *.CEL data as tree in data file
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestA1.CEL","TestA1");
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestA2.CEL","TestA2");

// cleanup
   manager->Close();
   delete manager;
}//ImportDataTest3

//______________________________________________________________________________
void UpdateDataTest3A(const char *filename   = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.3.12/DataTest3_cel.root",
                      const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.3.12/SchemeTest3.root")
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager");

// initialize chip type and variable list
   manager->Initialize("GeneChip");
   manager->InitInput("Test3","cel","MEAN/D:STDV/D:NPIXELS/I","RawData");

// open root scheme file
   manager->OpenSchemes(schemefile);

// update root data file (filname must end with *.root)
   manager->Update(filename);

   manager->BeginTransaction();
   manager->ProjectInfo("TestProject",20070523,"Project Update","use Test3 data for updating","my comment");
   manager->AuthorInfo("Stratowa","Christian","Project Leader","Home","home","email","++43-1-1234","my comment",kTRUE);
   manager->DatasetInfo("DataSet","MC","Tissue","StratowaUp",20070106,"description up","my comment",kTRUE);
   manager->HybridizationInfo("TestB1","hyb type","TestB1.CEL",20071117,"my prep1","standard protocol","B1",2,"my comment");
   manager->HybridizationInfo("TestB2","hyb type","TestB2.CEL",20071117,"my prep2","standard protocol","B2",2,"my comment");
   manager->TreatmentInfo("TestB1","DrugA2",4.3,"mM",1.0,"hours","intravenous","my comment");
   manager->TreatmentInfo("TestB2","DrugA2",4.3,"mM",8.0,"hours","intravenous","my comment");
   manager->CommitTransaction();

// store *.CEL data for mix as tree in data file
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestB1.CEL","TestB1");
   manager->Import("DataSet","/Volumes/CoreData/ROOT/rootdata/testAB/raw/TestB2.CEL","TestB2");

// cleanup
   manager->Close();
   delete manager;
}//UpdateDataTest3

//______________________________________________________________________________
void ImportDataU133P2(const char *filename   = "U133P2Tissues",
                      const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1",
                      const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root")
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager");

// initialize chip type and variable list
   manager->Initialize("GeneChip");
   manager->InitInput("HG-U133_Plus_2","cel","MEAN/D:STDV/D:NPIXELS/I","RawData");

// create new root data file 
   manager->New(filename, filedir, "GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// store *.CEL data as tree in data file
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/u1332plus_ivt_breast_A.CEL","BreastA");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/u1332plus_ivt_breast_B.CEL","BreastB");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/u1332plus_ivt_breast_C.CEL","BreastC");

   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/u1332plus_ivt_prostate_A.CEL","ProstateA");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/u1332plus_ivt_prostate_B.CEL","ProstateB");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/u1332plus_ivt_prostate_C.CEL","ProstateC");

// cleanup
   manager->Close();
   delete manager;
}//ImportDataU133P2

//______________________________________________________________________________
void ImportDataHuGene(const char *filename   = "HuGeneTissues",
                      const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1",
                      const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root")
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager");

// initialize chip type and variable list
   manager->Initialize("GenomeChip");
   manager->InitInput("HuGene-1_0-st-v1","cel","MEAN/D:STDV/D:NPIXELS/I","RawData");

// create new root data file 
   manager->New(filename, filedir, "GenomeChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// store *.CEL data as tree in data file
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuGene/TisMap_Breast_01_v1_WTGene1.CEL","Breast01");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuGene/TisMap_Breast_02_v1_WTGene1.CEL","Breast02");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuGene/TisMap_Breast_03_v1_WTGene1.CEL","Breast03");

   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuGene/TisMap_Prostate_01_v1_WTGene1.CEL","Prostate01");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuGene/TisMap_Prostate_02_v1_WTGene1.CEL","Prostate02");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuGene/TisMap_Prostate_03_v1_WTGene1.CEL","Prostate03");

// cleanup
   manager->Close();
   delete manager;
}//ImportDataHuGene

//______________________________________________________________________________
void ImportDataHuExon(const char *filename   = "HuExonTissues",
                      const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1",
                      const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root")
{
// Import Affymetrix *.CEL files into XPS

// create new data manager
   XDataManager *manager = new XDataManager("DataManager");

// initialize chip type and variable list
   manager->Initialize("ExonChip");
   manager->InitInput("HuEx-1_0-st-v2","cel","MEAN/D:STDV/D:NPIXELS/I","RawData");

// create new root data file 
   manager->New(filename, filedir, "ExonChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// store *.CEL data as tree in data file
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/huex_wta_breast_A.CEL","BreastA");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/huex_wta_breast_B.CEL","BreastB");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/huex_wta_breast_C.CEL","BreastC");

   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/huex_wta_prostate_A.CEL","ProstateA");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/huex_wta_prostate_B.CEL","ProstateB");
   manager->Import("DataSet","/Volumes/GigaDrive/ChipData/Exon/HuMixture/huex_wta_prostate_C.CEL","ProstateC");

// cleanup
   manager->Close();
   delete manager;
}//ImportDataHuExon

//______________________________________________________________________________
void ExportDataset(const char *schemefile, const char *datafile, const char *treeset,
                   const char *varlist = "*", const char *outfile = "PivotDataset",
                   const char *sep = "\t")
{
// Export raw data trees to "outfile"

// create new data manager
   XDataManager *manager = 0;
   manager = new XDataManager("DataManager");

   manager->Initialize("GenomeChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// open root file (necessary!)
   manager->Open(datafile);

// all trees from treeset for export
   manager->AddTree("SetExport", treeset);

// export all trees from Test3Set in datafile
   manager->ExportSet("SetExport", "cel", varlist, outfile, sep);

// cleanup
   manager->Close();
   delete manager;
}//ExportDataset


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Preprocess raw data                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void RMATest3(const char *filename   = "Test3_rma",
              const char *filedir    = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1",
              const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root")
{
// Preprocess RMA

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize backgrounder
//?   manager->InitAlgorithm("selector","probe","none", 0, 0);
   manager->InitAlgorithm("selector","probe","pmonly", 0, 0);
   manager->InitAlgorithm("backgrounder","rma","pmonly:epanechnikov",0, 1,16384);

// initialize normalizer
   manager->InitAlgorithm("selector","probe","pmonly", 0);
   manager->InitAlgorithm("normalizer","quantile","together:none:0",0, 1,0.0);

// initialize expressor
   manager->InitAlgorithm("selector","probe","pmonly", 0);
   manager->InitAlgorithm("expressor","medianpolish","log2",0,3,10,0.01,1.0);

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("RMASet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA1.cel");
   manager->AddTree("RMASet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA2.cel");
   manager->AddTree("RMASet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB1.cel");
   manager->AddTree("RMASet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB2.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("RMASet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//RMATest3

//______________________________________________________________________________
void MAS5Test3(const char *filename   = "Test3_mas5",
               const char *filedir    = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1",
               const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root")
{
// Preprocess MAS5

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","both", 0, 0);
   manager->InitAlgorithm("backgrounder","weightedsector","correctbg", "", 6,0.02,4,4,0,100,0.5);

// initialize expressor
   manager->InitAlgorithm("selector","probe","none", 0);
   manager->InitAlgorithm("expressor","TukeyBiweight","log2","",7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);

// initialize call detector
//   manager->InitAlgorithm("selector","default","none", 0);
//   manager->InitAlgorithm("calldetector","dc5","raw",0, 6,0.015,0.04,0.06,1,0,0);

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("MAS5Set","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA1.cel");
   manager->AddTree("MAS5Set","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA2.cel");
   manager->AddTree("MAS5Set","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB1.cel");
   manager->AddTree("MAS5Set","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB2.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("MAS5Set", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5Test3

//______________________________________________________________________________
void MAS5CallTest3(const char *filename   = "Test3_dc5",
                   const char *filedir    = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1",
                   const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root")
{
// Preprocess MAS5 call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize call detector
   manager->InitAlgorithm("selector","probe","none", 0);
   manager->InitAlgorithm("calldetector","dc5","raw",0, 6,0.015,0.04,0.06,1,0,0);

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("MAS5CallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA1.cel");
   manager->AddTree("MAS5CallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA2.cel");
   manager->AddTree("MAS5CallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB1.cel");
   manager->AddTree("MAS5CallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB2.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("MAS5CallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5CallTest3

//______________________________________________________________________________
void DABGCallTest3(const char *filename   = "Test3_dabg",
                   const char *filedir    = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1",
                   const char *schemefile = "/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/SchemeTest3.root")
{
// Preprocess DABG call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize call detector
   manager->InitAlgorithm("selector","probe","none", 0);
   manager->InitAlgorithm("calldetector","dab","raw",0, 3, 2.0, 0.01, 0.015); //cut, alpha1, alpha2

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("DABGCallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA1.cel");
   manager->AddTree("DABGCallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestA2.cel");
   manager->AddTree("DABGCallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB1.cel");
   manager->AddTree("DABGCallSet","/Volumes/CoreData/ROOT/rootdata/testAB/xps-0.4.1/DataTest3_cel.root/DataSet/TestB2.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("DABGCallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//DABGCallTest3

//______________________________________________________________________________
void RMAU133P2(const char *filename   = "U133P2_rma",
               const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1",
               const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root")
{
// Preprocess RMA

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","pmonly", 0, 0);
   manager->InitAlgorithm("backgrounder","rma","pmonly:epanechnikov",0, 1,16384);

// initialize normalizer
   manager->InitAlgorithm("selector","probe","pmonly", 0);
   manager->InitAlgorithm("normalizer","quantile","together:none:0",0, 1,0.0);

// initialize expressor
   manager->InitAlgorithm("selector","probe","pmonly", 0);
   manager->InitAlgorithm("expressor","medianpolish","log2",0,3,10,0.01,1.0);

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateC.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("RMASet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//RMAU133P2

//______________________________________________________________________________
void MAS5U133P2(const char *filename   = "U133P2_mas5",
                const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1",
                const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root")
{
// Preprocess MAS5

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","both", 0, 0);
   manager->InitAlgorithm("backgrounder","weightedsector","correctbg", "", 6,0.02,4,4,0,100,0.5);

// initialize expressor
   manager->InitAlgorithm("selector","probe","none", 0);
   manager->InitAlgorithm("expressor","TukeyBiweight","log2","",7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);

// initialize call detector
//   manager->InitAlgorithm("selector","default","none", 0);
//   manager->InitAlgorithm("calldetector","dc5","raw",0, 6,0.015,0.04,0.06,1,0,0);

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateC.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("MAS5Set", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5U133P2

//______________________________________________________________________________
void MAS5CallU133P2(const char *filename   = "U133P2_dc5",
                    const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1",
                    const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/SchemeU133P2.root")
{
// Preprocess MAS5 call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GeneChip");

// initialize call detector
   manager->InitAlgorithm("selector","probe","none", 0);
   manager->InitAlgorithm("calldetector","dc5","raw",0, 6,0.015,0.04,0.06,1,0,0);

// create new root data file 
   manager->New(filename,filedir,"GeneChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for preprocessing
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/u133p2/xps-0.4.1/U133P2Tissues_cel.root/DataSet/ProstateC.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("MAS5CallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5CallU133P2

//______________________________________________________________________________
void RMAHuGene(const char *filename   = "HuGene_rma", Int_t level = 9216,
               const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1",
               const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root")
{
// Preprocess RMA

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GenomeChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","genome",0, 2, level, -2);
   manager->InitAlgorithm("backgrounder","rma","pmonly:epanechnikov",0, 1,16384);

// initialize normalizer
   manager->InitAlgorithm("selector","probe","genome",0, 1, level);
   manager->InitAlgorithm("normalizer","quantile","together:none:0","", 1,0.0);
//   manager->InitAlgorithm("normalizer","quantile","together:none:0","/Volumes/GigaDrive/ROOT/rootcode/xps-0.4.1/tmp_rkq.root", 1,0.0);

// initialize expressor
   manager->InitAlgorithm("selector","probe","genome",0, 1, level);
   manager->InitAlgorithm("expressor","medianpolish","log2","",3,10,0.01,1.0);
//   manager->InitAlgorithm("expressor","medianpolish","log2","/Volumes/GigaDrive/ROOT/rootcode/xps-0.4.1/tmp.root",3,10,0.01,1.0);

// create new root data file 
   manager->New(filename, filedir, "GenomeChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast01.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast02.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast03.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate01.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate02.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate03.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("RMASet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//RMAHuGene

//______________________________________________________________________________
void MAS5HuGene(const char *filename   = "HuGene_mas5", Int_t level = 9216,
                const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1",
                const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root")
{
// Preprocess MAS5

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GenomeChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","genome",0, 1, level);
   manager->InitAlgorithm("backgrounder","weightedsector","correctbg",0, 6,0.02,4,4,0,100,0.5);

// initialize expressor
   manager->InitAlgorithm("selector","probe","genome",0, 1, level);
   manager->InitAlgorithm("expressor","TukeyBiweight","log2","",7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);

// create new root data file 
   manager->New(filename, filedir, "GenomeChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast01.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast02.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast03.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate01.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate02.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate03.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("MAS5Set", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5HuGene

//______________________________________________________________________________
void MAS5CallHuGene(const char *filename   = "HuGene_dc5", Int_t level = 9216,
                    const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1",
                    const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root")
{
// Preprocess MAS5 call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GenomeChip");

// initialize backgrounder (bgrd needed to compute MM values)
   manager->InitAlgorithm("selector","probe","genome",0, 1, level); 
   manager->InitAlgorithm("backgrounder","weightedsector","correctbg",0, 6,0.02,4,4,0,100,0.5);

// initialize call detector
   manager->InitAlgorithm("selector","probe","genome",0, 2, level, -2);
   manager->InitAlgorithm("calldetector","dc5","adjusted",0, 6,0.015,0.04,0.06,1,0,0);

// create new root data file 
   manager->New(filename, filedir, "GenomeChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast01.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast02.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast03.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate01.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate02.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate03.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("MAS5CallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5CallHuGene

//______________________________________________________________________________
void DABGCallHuGene(const char *filename   = "HuGene_dabg", Int_t level = 9216,
                    const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1",
                    const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HumanGeneScheme_v1.root")
{
// Preprocess DABG call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("GenomeChip");

// initialize call detector
   manager->InitAlgorithm("selector","probe","genome",0, 2, level, -2);
   manager->InitAlgorithm("calldetector","dab","raw",0, 3, 2.0, 0.01, 0.015); //cut, alpha1, alpha2

// create new root data file 
   manager->New(filename, filedir, "GenomeChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast01.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast02.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Breast03.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate01.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate02.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/hugene/xps-0.4.1/HuGeneTissues_cel.root/DataSet/Prostate03.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("DABGCallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//DABGCallHuGene

//______________________________________________________________________________
void RMAHuExon(const char *filename   = "HuExon_rma", Int_t level = 9216,
               const char *exproption = "log2",
               const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1",
               const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root")
{
// Preprocess RMA

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("ExonChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","exon",0, 2, level, -2);
   manager->InitAlgorithm("backgrounder","rma","pmonly:epanechnikov",0, 1,16384);

// initialize normalizer
   manager->InitAlgorithm("selector","probe","exon",0, 1, level);
   manager->InitAlgorithm("normalizer","quantile","together:none:0","", 1,0.0);
//   manager->InitAlgorithm("normalizer","quantile","together:none:0","/Volumes/GigaDrive/ROOT/rootcode/xps-0.4.1/tmp_rkq.root", 1,0.0);

// initialize expressor
   manager->InitAlgorithm("selector","probe","exon",0, 1, level);
   manager->InitAlgorithm("expressor","medianpolish",exproption,"",3,10,0.01,1.0);
//   manager->InitAlgorithm("expressor","medianpolish",exproption,"/Volumes/GigaDrive/ROOT/rootcode/xps-0.4.1/tmp.root",3,10,0.01,1.0);

// create new root data file 
   manager->New(filename, filedir, "ExonChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("RMASet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateC.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("RMASet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//RMAHuExon

//______________________________________________________________________________
void MAS5HuExon(const char *filename   = "HuExon_mas5", Int_t level = 9216,
                const char *exproption = "log2",
                const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1",
                const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root")
{
// Preprocess MAS5

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("ExonChip");

// initialize backgrounder
   manager->InitAlgorithm("selector","probe","exon",0, 1, level);
   manager->InitAlgorithm("backgrounder","weightedsector","correctbg",0, 6,0.02,4,4,0,100,0.5);

// initialize expressor
   manager->InitAlgorithm("selector","probe","exon",0, 1, level);
   manager->InitAlgorithm("expressor","TukeyBiweight",exproption,"",7,0.03,10.0,2.0e-20,5.0,0.0001,1.0,0.5);

// create new root data file 
   manager->New(filename, filedir, "ExonChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("MAS5Set","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateC.cel");

// preprocess expression values and store as trees in new file
   manager->Preprocess("MAS5Set", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5HuExon

//______________________________________________________________________________
void MAS5CallHuExon(const char *filename   = "HuExon_dc5", Int_t level = 9216,
                    const char *calloption = "adjusted",
                    const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1",
                    const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root")
{
// Preprocess MAS5 call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("ExonChip");

// initialize backgrounder (bgrd needed to compute MM values)
   manager->InitAlgorithm("selector","probe","exon",0, 1, level); 
//?   manager->InitAlgorithm("backgrounder","weightedsector","correctbg",0, 6,0.02,4,4,0,100,0.5);
// the following setting calculates bg but does not subtract bg from intensity:
   manager->InitAlgorithm("backgrounder","weightedsector","none", "", 6,0.02,4,4,0,100,-1.0);

// initialize call detector
   manager->InitAlgorithm("selector","probe","exon",0, 2, level, -2);
   manager->InitAlgorithm("calldetector","dc5",calloption, 0, 6,0.015,0.04,0.06,1,0,0);

// create new root data file 
   manager->New(filename, filedir, "ExonChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("MAS5CallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateC.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("MAS5CallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//MAS5CallHuExon

//______________________________________________________________________________
void DABGCallHuExon(const char *filename   = "HuExon_dabg", Int_t level = 9216,
                    const char *calloption = "raw",
                    const char *filedir    = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1",
                    const char *schemefile = "/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HumanExonScheme_r2.root")
{
// Preprocess DABG call

// create new preprocessing manager
   XPreProcessManager *manager = new XPreProcessManager("PreProcessManager");

// initialize preprocessing algorithms
   manager->Initialize("ExonChip");

// initialize call detector
   manager->InitAlgorithm("selector","probe","exon",0, 2, level, -2); 
   manager->InitAlgorithm("calldetector","dab",calloption, 0, 3, 2.0, 0.01, 0.015); //cut, alpha1, alpha2

// create new root data file 
   manager->New(filename, filedir, "ExonChip");

// open root scheme file
   manager->OpenSchemes(schemefile);

// add trees for mas
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastA.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastB.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/BreastC.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateA.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateB.cel");
   manager->AddTree("DABGCallSet","/Volumes/GigaDrive/ROOT/rootdata/huexon/xps-0.4.1/HuExonTissues_cel.root/DataSet/ProstateC.cel");

// preprocess call values and store as trees in new file
   manager->Preprocess("DABGCallSet", "preprocess");

// cleanup
   manager->Close();
   delete manager;
}//DABGCallHuExon

//______________________________________________________________________________
void ExportPreprocess(const char *chiptype, const char *schemefile, const char *rootfile, 
                      const char *treename, const char *varlist, const char *outfile = "")
{
// macro to export tree "treename" from file "rootfile"

// create new preprocessing manager
   XPreProcessManager *manager = 0;
   manager = new XPreProcessManager("PreProcessManager");

   manager->Initialize(chiptype);

// open root scheme file
   manager->OpenSchemes(schemefile);

// open root file
   manager->Open(rootfile);

// export tree
   manager->Export(treename, varlist, outfile);

// cleanup
   manager->Close();
   delete manager;
}//ExportPreprocess

//______________________________________________________________________________
//______________________________________________________________________________
