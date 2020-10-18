#------------------------------------------------------------------------------#
# Script: step-by-step functions to create ROOT scheme files for package "xps"
#
# Note: please feel free to copy-paste the examples of interest and adapt the
#       examples to your own needs
#
# Copyright (c) 2010-2012 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# 1. step: import Affymetrix chip definition and annotation files into
#          ROOT scheme files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files need only be created when new annotation files
#          are available. In order to allow all users to access the ROOT 
#          scheme files it is recommended to create all ROOT scheme files
#          in a common system directory and use function "root.scheme()"
#          to load the desired scheme into the current R session to analyse
#          a project.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: do not separate name of ROOT files with dots, use underscores instead.
#------------------------------------------------------------------------------#

### new R session: load library xps
library(xps)

### define directories:
# directory containing Affymetrix library files
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
# directory containing Affymetrix annotation files
anndir <- "/Volumes/GigaDrive/Affy/Annotation"
# directory to store ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"


#------------------------------------------------------------------------------#
# Nov 2012: Affymetrix annotation "na33"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na33 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na33"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version12Nov", "Test3.na33.annot.csv"))

# Hu6800: 
# annotation file for na32 no longer provided by Affymetrix (see na30 annotations)

# HG-U133A:
scheme.hgu133a.na33 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na33"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version12Nov", "HG-U133A.na33.annot.csv"))

# HG-U133B:
scheme.hgu133b.na33 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na33"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version12Nov", "HG-U133B.na33.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na33 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na33"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version12Nov", "HG-U133_Plus_2.na33.annot.csv"))

# Rice:
scheme.rice.na33 <- import.expr.scheme("rice", filedir = file.path(scmdir, "na33"),
                     schemefile = file.path(libdir, "Rice.cdf"), 
                     probefile  = file.path(libdir, "Rice.probe.tab"), 
                     annotfile  = file.path(anndir, "Version12Nov", "Rice.na33.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na33 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version12Nov", "HuGene-1_0-st-v1.na33.hg19.probeset.csv"),
                            file.path(anndir, "Version12Nov", "HuGene-1_0-st-v1.na33.hg19.transcript.csv"))

# HuGene-2_0-st: used as exon array
scheme.hugene20st.na33 <- import.exon.scheme("hugene20stv1", filedir = file.path(scmdir, "na33"),
                          file.path(libdir, "HuGene-2_0-st", "HuGene-2_0-st.clf"),
                          file.path(libdir, "HuGene-2_0-st", "HuGene-2_0-st.pgf"),
                          file.path(anndir, "Version12Nov", "HuGene-2_0-st-v1.na33.hg19.probeset.csv"),
                          file.path(anndir, "Version12Nov",  "HuGene-2_0-st-v1.na33.hg19.transcript.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na33 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version12Nov", "MoGene-1_0-st-v1.na33.mm9.probeset.csv"),
                            file.path(anndir, "Version12Nov", "MoGene-1_0-st-v1.na33.mm9.transcript.csv"))

# MoGene-2_0-st: used as exon array
scheme.mogene20st.na33 <- import.exon.scheme("mogene20stv1", filedir = file.path(scmdir, "na33"),
                          file.path(libdir, "MoGene-2_0-st", "MoGene-2_0-st.clf"),
                          file.path(libdir, "MoGene-2_0-st", "MoGene-2_0-st.pgf"),
                          file.path(anndir, "Version12Nov", "MoGene-2_0-st-v1.na33.mm10.probeset.csv"),
                          file.path(anndir, "Version12Nov",  "MoGene-2_0-st-v1.na33.mm10.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na33 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version12Nov", "RaGene-1_0-st-v1.na33.rn4.probeset.csv"),
                            file.path(anndir, "Version12Nov", "RaGene-1_0-st-v1.na33.rn4.transcript.csv"))

# RaGene-2_0-st: used as exon array
scheme.ragene20st.na33 <- import.exon.scheme("ragene20stv1", filedir = file.path(scmdir, "na33"),
                          file.path(libdir, "RaGene-2_0-st", "RaGene-2_0-st.clf"),
                          file.path(libdir, "RaGene-2_0-st", "RaGene-2_0-st.pgf"),
                          file.path(anndir, "Version12Nov", "RaGene-2_0-st-v1.na33.rn4.probeset.csv"),
                          file.path(anndir, "Version12Nov",  "RaGene-2_0-st-v1.na33.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na33 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na33"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version12Nov", "HuEx-1_0-st-v2.na33.hg19.probeset.csv"),
                          file.path(anndir, "Version12Nov", "HuEx-1_0-st-v2.na33.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1.na33 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na33"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version12Nov", "MoEx-1_0-st-v1.na33.mm9.probeset.csv"),
                          file.path(anndir, "Version12Nov", "MoEx-1_0-st-v1.na33.mm9.transcript.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na33 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na33"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version12Nov", "RaEx-1_0-st-v1.na33.rn4.probeset.csv"),
                          file.path(anndir, "Version12Nov", "RaEx-1_0-st-v1.na33.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na33 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version12Nov", "HuGene-1_1-st-v1.na33.hg19.probeset.csv"),
                            file.path(anndir, "Version12Nov", "HuGene-1_1-st-v1.na33.hg19.transcript.csv"))

# HuGene-2_1-st: 
scheme.hugene21st.na33 <- import.exon.scheme("hugene21stv1", filedir = file.path(scmdir, "na33"),
                          file.path(libdir, "HuGene-2_1-st", "HuGene-2_1-st.clf"),
                          file.path(libdir, "HuGene-2_1-st", "HuGene-2_1-st.pgf"),
                          file.path(anndir, "Version12Nov",  "HuGene-2_1-st-v1.na33.hg19.probeset.csv"),
                          file.path(anndir, "Version12Nov",  "HuGene-2_1-st-v1.na33.hg19.transcript.csv"))

# MoGene-1_1-st-v1.r4:
scheme.mogene11stv1.na33 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version12Nov", "MoGene-1_1-st-v1.na33.mm9.probeset.csv"),
                            file.path(anndir, "Version12Nov", "MoGene-1_1-st-v1.na33.mm9.transcript.csv"))

# MoGene-2_1-st:
scheme.mogene21stv1.na33 <- import.exon.scheme("mogene21stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "MoGene-2_1-st", "MoGene-2_1-st.clf"),
                            file.path(libdir, "MoGene-2_1-st", "MoGene-2_1-st.pgf"),
                            file.path(anndir, "Version12Nov", "MoGene-2_1-st-v1.na33.mm10.probeset.csv"),
                            file.path(anndir, "Version12Nov", "MoGene-2_1-st-v1.na33.mm10.transcript.csv"))

# RaGene-1_1-st-v1.r4:
scheme.ragene11stv1.na33 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version12Nov", "RaGene-1_1-st-v1.na33.rn4.probeset.csv"),
                            file.path(anndir, "Version12Nov", "RaGene-1_1-st-v1.na33.rn4.transcript.csv"))

# RaGene-2_1-st:
scheme.ragene21stv1.na33 <- import.exon.scheme("ragene21stv1", filedir = file.path(scmdir, "na33"),
                            file.path(libdir, "RaGene-2_1-st", "RaGene-2_1-st.clf"),
                            file.path(libdir, "RaGene-2_1-st", "RaGene-2_1-st.pgf"),
                            file.path(anndir, "Version12Nov", "RaGene-2_1-st-v1.na33.rn4.probeset.csv"),
                            file.path(anndir, "Version12Nov", "RaGene-2_1-st-v1.na33.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for model organisms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# RUSGene-1_1-st-v1:
scheme.rusgene11stv1 <- import.exon.scheme("rusgene11stv1", filedir = file.path(scmdir, "designtime"),
                        file.path(libdir, "RUSGene-1_1-st_rev2", "RUSGene-1_1-st.clf"),
                        file.path(libdir, "RUSGene-1_1-st_rev2", "RUSGene-1_1-st.pgf"),
                        file.path(anndir, "ModelOrganisms", "RUSGene-1_1-st-v1.design-time.20130417.probeset.csv"),
                        file.path(anndir, "ModelOrganisms", "RUSGene-1_1-st-v1.design-time.20130417.transcript.csv"))



#------------------------------------------------------------------------------#
# Jul 2011: Affymetrix annotation "na32"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na32 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na32"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version11Jul", "Test3.na32.annot.csv"))

# Hu6800: 
# annotation file for na32 no longer provided by Affymetrix (see na30 annotations)

# HG_U95A:
scheme.hgu95a.na32 <- import.expr.scheme("hgu95a", filedir = file.path(scmdir, "na32"),
                      schemefile = file.path(libdir, "HG_U95A.CDF"), 
                      probefile  = file.path(libdir, "HG-U95A_probe.tab"), 
                      annotfile  = file.path(anndir, "Version11Jul", "HG_U95A.na32.annot.csv"))

# HG_U95Av2:
scheme.hgu95av2.na32 <- import.expr.scheme("hgu95av2", filedir = file.path(scmdir, "na32"),
                        schemefile = file.path(libdir, "HG_U95Av2.CDF"), 
                        probefile  = file.path(libdir, "HG-U95Av2_probe.tab"), 
                        annotfile  = file.path(anndir, "Version11Jul", "HG_U95Av2.na32.annot.csv"))

# HG-U133A:
scheme.hgu133a.na32 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na32"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version11Jul", "HG-U133A.na32.annot.csv"))

# HG-U133B:
scheme.hgu133b.na32 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na32"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version11Jul", "HG-U133B.na32.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na32 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na32"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version11Jul", "HG-U133_Plus_2.na32.annot.csv"))

# PrimeView:
scheme.primeview.na32 <- import.expr.scheme("primeview", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "PrimeView.CDF"), 
                         probefile  = file.path(libdir, "PrimeView.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "PrimeView.na32.annot.csv"))

# MOE430A:
scheme.moe430a.na32 <- import.expr.scheme("moe430a",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/MOE430A.na32.annot.csv"))

# MOE430B:
scheme.moe430b.na32 <- import.expr.scheme("moe430b",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/MOE430B.na32.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na32 <- import.expr.scheme("mouse4302",filedir=file.path(scmdir, "na32"),
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version11Jul/Mouse430_2.na32.annot.csv"))

# RG_U34A:
scheme.rgu34a.na32 <- import.expr.scheme("rgu34a",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RG_U34A.na32.annot.csv"))

# RG_U34B:
scheme.rgu34b.na32 <- import.expr.scheme("rgu34b",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RG_U34B.na32.annot.csv"))

# RAE230A:
scheme.rae230a.na32 <- import.expr.scheme("rae230a",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RAE230A.na32.annot.csv"))

# RAE230B:
scheme.rae230b.na32 <- import.expr.scheme("rae230b",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/RAE230B.na32.annot.csv"))

# Rat230_2:
scheme.rat2302.na32 <- import.expr.scheme("rat2302",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/Rat230_2.na32.annot.csv"))

# Bovine:
scheme.bovine.na32 <- import.expr.scheme("bovine",filedir=file.path(scmdir, "na32"),
                      schemefile = file.path(libdir,"Bovine.cdf"),
                      probefile  = file.path(libdir,"Bovine.probe.tab"),
                      annotfile  = file.path(anndir,"Version11Jul/Bovine.na32.annot.csv"))

# Porcine:
scheme.porcine.na32 <- import.expr.scheme("porcine",filedir=file.path(scmdir, "na32"),
                       schemefile = file.path(libdir,"Porcine.cdf"),
                       probefile  = file.path(libdir,"Porcine.probe.tab"),
                       annotfile  = file.path(anndir,"Version11Jul/Porcine.na32.annot.csv"))

# Rhesus:
scheme.rhesus.na32 <- import.expr.scheme("rhesus",filedir=file.path(scmdir, "na32"),
                      schemefile = file.path(libdir,"Rhesus.cdf"),
                      probefile  = file.path(libdir,"Rhesus.probe.tab"),
                      annotfile  = file.path(anndir,"Version11Jul/Rhesus.na32.annot.csv"))

# Zebrafish:
scheme.zebrafish.na32 <- import.expr.scheme("zebrafish",filedir=file.path(scmdir, "na32"),
                         schemefile = file.path(libdir,"Zebrafish.cdf"),
                         probefile  = file.path(libdir,"Zebrafish.probe.tab"),
                         annotfile  = file.path(anndir,"Version11Jul/Zebrafish.na32.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na32 <- import.expr.scheme("hthgu133a", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133A.na32.annot.csv"))

# HT_HG-U133B
scheme.hthgu133b.na32 <- import.expr.scheme("hthgu133b", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_HG-U133B.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133B.na32.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na32 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na32"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133_Plus_PM.na32.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na32 <- import.expr.scheme("htmg430a", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_MG-430A.na32.annot.csv"))

# HT_MG-430B
scheme.htmg430b.na32 <- import.expr.scheme("htmg430b", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_MG-430B.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_MG-430B.na32.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na32 <- import.expr.scheme("htmg430pm", filedir = file.path(scmdir, "na32"),
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version11Jul", "HT_MG-430_PM.na32.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na32 <- import.expr.scheme("htrat230pm", filedir = file.path(scmdir, "na32"),
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version11Jul", "HT_Rat230_PM.na32.annot.csv"))

# HG-U219
scheme.hgu219.na32 <- import.expr.scheme("hgu219", filedir = file.path(scmdir, "na32"),
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version11Jul", "HG-U219.na32.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na32 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.probeset.csv"),
                            file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.transcript.csv"))

# HuGene-2_0-st: used as exon array
# use perl script "HuGene20_update_AFFX.pl" to create corrected annotation files, see: 
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2012-August/047755.html
scheme.hugene20st.na32 <- import.exon.scheme("hugene20stv1", filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "HuGene-2_0-st", "HuGene-2_0-st.clf"),
                          file.path(libdir, "HuGene-2_0-st", "HuGene-2_0-st.pgf"),
                          file.path(anndir, "Version11Jul", "HuGene-2_0-st-v1.na32.hg19.probeset.csv", "HuGene-2_0-st-v1.na32.hg19.probeset.corr.csv"),
                          file.path(anndir, "Version11Jul", "HuGene-2_0-st-v1.na32.hg19.transcript.csv", "HuGene-2_0-st-v1.na32.hg19.transcript.corr.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na32 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "MoGene-1_0-st-v1.na32.mm9.probeset.csv"),
                            file.path(anndir, "Version11Jul", "MoGene-1_0-st-v1.na32.mm9.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na32 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "RaGene-1_0-st-v1.na32.rn4.probeset.csv"),
                            file.path(anndir, "Version11Jul", "RaGene-1_0-st-v1.na32.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na32 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.probeset.csv"),
                          file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1.na32 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version11Jul", "MoEx-1_0-st-v1.na32.mm9.probeset.csv"),
                          file.path(anndir, "Version11Jul", "MoEx-1_0-st-v1.na32.mm9.transcript.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na32 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na32"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version11Jul", "RaEx-1_0-st-v1.na32.rn4.probeset.csv"),
                          file.path(anndir, "Version11Jul", "RaEx-1_0-st-v1.na32.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na32 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "HuGene-1_1-st-v1.na32.hg19.probeset.csv"),
                            file.path(anndir, "Version11Jul", "HuGene-1_1-st-v1.na32.hg19.transcript.csv"))

# HuGene-2_1-st: 
# use perl script "HuGene21_update_AFFX.pl" to create corrected annotation files, see: 
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2012-August/047755.html
scheme.hugene21st.na32 <- import.exon.scheme("hugene21stv1", filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "HuGene-2_1-st", "HuGene-2_1-st.clf"),
                          file.path(libdir, "HuGene-2_1-st", "HuGene-2_1-st.pgf"),
                          file.path(anndir, "Version11Jul", "HuGene-2_1-st-v1.na32.hg19.probeset.csv", "HuGene-2_1-st-v1.na32.hg19.probeset.corr.csv"),
                          file.path(anndir, "Version11Jul", "HuGene-2_1-st-v1.na32.hg19.transcript.csv", "HuGene-2_1-st-v1.na32.hg19.transcript.corr.csv"))

# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na32 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "MoGene-1_1-st-v1.na32.mm9.probeset.csv"),
                            file.path(anndir, "Version11Jul", "MoGene-1_1-st-v1.na32.mm9.transcript.csv"))

# RaGene-1_1-st-v1.r4
# use updated Affymetrix annotation files from 09/17/10
scheme.ragene11stv1.na32 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "RaGene-1_1-st-v1.na32.rn4.probeset.csv"),
                            file.path(anndir, "Version11Jul", "RaGene-1_1-st-v1.na32.rn4.transcript.csv"))



#------------------------------------------------------------------------------#
# Sep 2010: Affymetrix annotation "na31"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na31 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na31"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version10Sep", "Test3.na31.annot.csv"))

# Hu6800: 
# annotation file for na31 no longer provided by Affymetrix (see na30 annotations)

# HG_U95A:
# annotation file for na31 no longer provided by Affymetrix (see na30 annotations)

# HG_U95Av2:
# annotation file for na31 no longer provided by Affymetrix (see na30 annotations)

# HG-U133A:
scheme.hgu133a.na31 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na31"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version10Sep", "HG-U133A.na31.annot.csv"))

# HG-U133B:
scheme.hgu133b.na31 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na31"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version10Sep", "HG-U133B.na31.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na31 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na31"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version10Sep", "HG-U133_Plus_2.na31.annot.csv"))

# MOE430A:
scheme.moe430a.na31 <- import.expr.scheme("moe430a",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/MOE430A.na31.annot.csv"))

# MOE430B:
scheme.moe430b.na31 <- import.expr.scheme("moe430b",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/MOE430B.na31.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na31 <- import.expr.scheme("mouse4302",filedir=file.path(scmdir, "na31"),
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version10Sep/Mouse430_2.na31.annot.csv"))

# RG_U34A:
scheme.rgu34a.na31 <- import.expr.scheme("rgu34a",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RG_U34A.na31.annot.csv"))

# RG_U34B:
scheme.rgu34b.na31 <- import.expr.scheme("rgu34b",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RG_U34B.na31.annot.csv"))

# RG_U34C:
scheme.rgu34c.na31 <- import.expr.scheme("rgu34c",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RG_U34C.cdf"),
                       probefile  = file.path(libdir,"RG_U34C.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RG_U34C.na31.annot.csv"))

# RAE230A:
scheme.rae230a.na31 <- import.expr.scheme("rae230a",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RAE230A.na31.annot.csv"))

# RAE230B:
scheme.rae230b.na31 <- import.expr.scheme("rae230b",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/RAE230B.na31.annot.csv"))

# Rat230_2:
scheme.rat2302.na31 <- import.expr.scheme("rat2302",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/Rat230_2.na31.annot.csv"))

# Bovine:
scheme.bovine.na31 <- import.expr.scheme("bovine",filedir=file.path(scmdir, "na31"),
                      schemefile = file.path(libdir,"Bovine.cdf"),
                      probefile  = file.path(libdir,"Bovine.probe.tab"),
                      annotfile  = file.path(anndir,"Version10Sep/Bovine.na31.annot.csv"))

# Porcine:
scheme.porcine.na31 <- import.expr.scheme("porcine",filedir=file.path(scmdir, "na31"),
                       schemefile = file.path(libdir,"Porcine.cdf"),
                       probefile  = file.path(libdir,"Porcine.probe.tab"),
                       annotfile  = file.path(anndir,"Version10Sep/Porcine.na31.annot.csv"))

# Rhesus:
scheme.rhesus.na31 <- import.expr.scheme("rhesus",filedir=file.path(scmdir, "na31"),
                      schemefile = file.path(libdir,"Rhesus.cdf"),
                      probefile  = file.path(libdir,"Rhesus.probe.tab"),
                      annotfile  = file.path(anndir,"Version10Sep/Rhesus.na31.annot.csv"))

# Zebrafish:
scheme.zebrafish.na31 <- import.expr.scheme("zebrafish",filedir=file.path(scmdir, "na31"),
                         schemefile = file.path(libdir,"Zebrafish.cdf"),
                         probefile  = file.path(libdir,"Zebrafish.probe.tab"),
                         annotfile  = file.path(anndir,"Version10Sep/Zebrafish.na31.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na31 <- import.expr.scheme("hthgu133a", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_HG-U133A.na31.annot.csv"))

# HT_HG-U133B
scheme.hthgu133b.na31 <- import.expr.scheme("hthgu133b", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_HG-U133B.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_HG-U133B.na31.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na31 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na31"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version10Sep", "HT_HG-U133_Plus_PM.na31.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na31 <- import.expr.scheme("htmg430a", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_MG-430A.na31.annot.csv"))

# HT_MG-430B
scheme.htmg430b.na31 <- import.expr.scheme("htmg430b", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_MG-430B.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430B.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_MG-430B.na31.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na31 <- import.expr.scheme("htmg430pm", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version10Sep", "HT_MG-430_PM.na31.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na31 <- import.expr.scheme("htrat230pm", filedir = file.path(scmdir, "na31"),
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version10Sep", "HT_Rat230_PM.na31.annot.csv"))

# HG-U219
scheme.hgu219.na31 <- import.expr.scheme("hgu219", filedir = file.path(scmdir, "na31"),
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version10Sep", "HG-U219.na31.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na31 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "HuGene-1_0-st-v1.na31.hg19.probeset.csv"),
                            file.path(anndir, "Version10Sep", "HuGene-1_0-st-v1.na31.hg19.transcript.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na31 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "MoGene-1_0-st-v1.na31.mm9.probeset.csv"),
                            file.path(anndir, "Version10Sep", "MoGene-1_0-st-v1.na31.mm9.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na31 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "RaGene-1_0-st-v1.na31.rn4.probeset.csv"),
                            file.path(anndir, "Version10Sep", "RaGene-1_0-st-v1.na31.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na31 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na31"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version10Sep", "HuEx-1_0-st-v2.na31.hg19.probeset.csv"),
                          file.path(anndir, "Version10Sep", "HuEx-1_0-st-v2.na31.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
# use updated Affymetrix annotation files from 09/08/10
scheme.moex10stv1.na31 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na31"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.probeset.csv"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.transcript.csv"))

# original annotation files from 08/30/10: need to delete control->affx which are neg_control
source(paste(path.package("xps"),"examples/updateAnnotation.R",sep="/"))
deleteNegControlFromAffxControl("MoEx-1_0-st-v1.na31.mm9.probeset.csv", "MoEx-1_0-st-v1.na31.mm9.probeset.fixed.csv", eol="\n")
deleteNegControlFromAffxControl("MoEx-1_0-st-v1.na31.mm9.transcript.csv", "MoEx-1_0-st-v1.na31.mm9.transcript.fixed.csv", eol="\n")
# use fixed annotation files
scheme.moex10stv1.na31 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na31"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.probeset.fixed.csv"),
                          file.path(anndir, "Version10Sep", "MoEx-1_0-st-v1.na31.mm9.transcript.fixed.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na31 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na31"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version10Sep", "RaEx-1_0-st-v1.na31.rn4.probeset.csv"),
                          file.path(anndir, "Version10Sep", "RaEx-1_0-st-v1.na31.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na31 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "HuGene-1_1-st-v1.na31.hg19.probeset.csv"),
                            file.path(anndir, "Version10Sep", "HuGene-1_1-st-v1.na31.hg19.transcript.csv"))

# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na31 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "MoGene-1_1-st-v1.na31.mm9.probeset.csv"),
                            file.path(anndir, "Version10Sep", "MoGene-1_1-st-v1.na31.mm9.transcript.csv"))

# RaGene-1_1-st-v1.r4
# use updated Affymetrix annotation files from 09/17/10
scheme.ragene11stv1.na31 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.probeset.csv"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.transcript.csv"))

# original annotation files from 08/27/10: need to add missing header lines
scheme.ragene11stv1.na31 <- import.exon.scheme("ragene11stv1", filedir = file.path(scmdir, "na31"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.probeset.fixed.csv"),
                            file.path(anndir, "Version10Sep", "RaGene-1_1-st-v1.na31.rn4.transcript.fixed.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# E_coli_2:
scheme.ecoli2.na31 <- import.expr.scheme("ecoli2", filedir = file.path(scmdir, "na31"),
                      schemefile = file.path(libdir, "E_coli_2.cdf"),
                      probefile  = file.path(libdir, "E_coli_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version10Sep", "E_coli_2.na31.annot.csv"))

# Ecoli_ASv2:
scheme.ecoliasv2.na31 <- import.expr.scheme("ecoliasv2", filedir = file.path(scmdir, "na31"),
                         schemefile = file.path(libdir, "Ecoli_ASv2.CDF"),
                         probefile  = file.path(libdir, "Ecoli_ASv2.probe.tab"),
                         annotfile  = file.path(anndir, "Version10Sep", "Ecoli_ASv2.na31.annot.csv"))

# Pae_G1a:
scheme.paeg1a.na31 <- import.expr.scheme("paeg1a", filedir = file.path(scmdir, "na31"),
                      schemefile = file.path(libdir, "Pae_G1a.CDF"),
                      probefile  = file.path(libdir, "Pae_G1a.probe.tab"), 
                      annotfile  = file.path(anndir, "Version10Sep", "Pae_G1a.na31.annot.csv"))

# Yeast_2:
scheme.yeast2.na31 <- import.expr.scheme("yeast2", filedir = file.path(scmdir, "na31"),
                      schemefile = file.path(libdir, "Yeast_2.cdf"),
                      probefile  = file.path(libdir, "Yeast_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version10Sep", "Yeast_2.na31.annot.csv"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna10 <- import.expr.scheme("mirna10", filedir = file.path(scmdir, "na31"),
                  schemefile = file.path(libdir, "miRNA-1_0.CDF"),
                  probefile  = file.path(libdir, "miRNA-1_0.probe.tab"))



#------------------------------------------------------------------------------#
# Nov 2009: Affymetrix annotation "na30"
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na30 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na30"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version09Nov", "Test3.na30.annot.csv"))

# Hu6800:
scheme.hu6800.na30 <- import.expr.scheme("hu6800", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "Hu6800.CDF"), 
                      probefile  = file.path(libdir, "HuGeneFL_probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "Hu6800.na30.annot.csv"))

# HG_U95A:
scheme.hgu95a.na30 <- import.expr.scheme("hgu95a", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "HG_U95A.CDF"), 
                      probefile  = file.path(libdir, "HG-U95A_probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "HG_U95A.na30.annot.csv"))

# HG_U95Av2:
scheme.hgu95av2.na30 <- import.expr.scheme("hgu95av2", filedir = file.path(scmdir, "na30"),
                        schemefile = file.path(libdir, "HG_U95Av2.CDF"), 
                        probefile  = file.path(libdir, "HG-U95Av2_probe.tab"), 
                        annotfile  = file.path(anndir, "Version09Nov", "HG_U95Av2.na30.annot.csv"))

# HG-U133A:
scheme.hgu133a.na30 <- import.expr.scheme("hgu133a", filedir = file.path(scmdir, "na30"),
                       schemefile = file.path(libdir, "HG-U133A.CDF"), 
                       probefile  = file.path(libdir, "HG-U133A_probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U133A.na30.annot.csv"))

# HG-U133B:
scheme.hgu133b.na30 <- import.expr.scheme("hgu133b", filedir = file.path(scmdir, "na30"),
                       schemefile = file.path(libdir, "HG-U133B.CDF"), 
                       probefile  = file.path(libdir, "HG-U133B_probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U133B.na30.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na30 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na30"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version09Nov", "HG-U133_Plus_2.na30.annot.csv"))

# MOE430A:
scheme.moe430a.na30 <- import.expr.scheme("moe430a",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"MOE430A.CDF"),
                       probefile  = file.path(libdir,"MOE430A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/MOE430A.na30.annot.csv"))

# MOE430B:
scheme.moe430b.na30 <- import.expr.scheme("moe430b",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"MOE430B.CDF"),
                       probefile  = file.path(libdir,"MOE430B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/MOE430B.na30.annot.csv"))

# Mouse430_2:
scheme.mouse4302.na30 <- import.expr.scheme("mouse4302",filedir=file.path(scmdir, "na30"),
                         schemefile = file.path(libdir,"Mouse430_2.cdf"),
                         probefile  = file.path(libdir,"Mouse430_2.probe.tab"),
                         annotfile  = file.path(anndir,"Version09Nov/Mouse430_2.na30.annot.csv"))

# RG_U34A:
scheme.rgu34a.na30 <- import.expr.scheme("rgu34a",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RG_U34A.cdf"),
                       probefile  = file.path(libdir,"RG_U34A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34A.na30.annot.csv"))

# RG_U34B:
scheme.rgu34b.na30 <- import.expr.scheme("rgu34b",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RG_U34B.cdf"),
                       probefile  = file.path(libdir,"RG_U34B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34B.na30.annot.csv"))

# RG_U34C:
scheme.rgu34c.na30 <- import.expr.scheme("rgu34c",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RG_U34C.cdf"),
                       probefile  = file.path(libdir,"RG_U34C.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RG_U34C.na30.annot.csv"))

# RAE230A:
scheme.rae230a.na30 <- import.expr.scheme("rae230a",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RAE230A.CDF"),
                       probefile  = file.path(libdir,"RAE230A.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RAE230A.na30.annot.csv"))

# RAE230B:
scheme.rae230b.na30 <- import.expr.scheme("rae230b",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"RAE230B.CDF"),
                       probefile  = file.path(libdir,"RAE230B.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/RAE230B.na30.annot.csv"))

# Rat230_2:
scheme.rat2302.na30 <- import.expr.scheme("rat2302",filedir=file.path(scmdir, "na30"),
                       schemefile = file.path(libdir,"Rat230_2.cdf"),
                       probefile  = file.path(libdir,"Rat230_2.probe.tab"),
                       annotfile  = file.path(anndir,"Version09Nov/Rat230_2.na30.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133a.na30 <- import.expr.scheme("hthgu133a", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "HT_HG-U133A.cdf"), 
                         probefile  = file.path(libdir, "HT_HG-U133A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_HG-U133A.na30.annot.csv"))

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na30 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na30"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version09Nov", "HT_HG-U133_Plus_PM.na30.annot.csv"))

# HT_MG-430A
scheme.htmg430a.na30 <- import.expr.scheme("htmg430a", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "HT_MG-430A.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430A.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_MG-430A.na30.annot.csv"))

# HT_MG-430_PM
scheme.htmg430pm.na30 <- import.expr.scheme("htmg430pm", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "HT_MG-430_PM.cdf"), 
                         probefile  = file.path(libdir, "HT_MG-430_PM.probe.tab"), 
                         annotfile  = file.path(anndir, "Version09Nov", "HT_MG-430_PM.na30.annot.csv"))

# HT_Rat230_PM
scheme.htrat230pm.na30 <- import.expr.scheme("htrat230pm", filedir = file.path(scmdir, "na30"),
                          schemefile = file.path(libdir, "HT_Rat230_PM.cdf"), 
                          probefile  = file.path(libdir, "HT_Rat230_PM.probe.tab"), 
                          annotfile  = file.path(anndir, "Version09Nov", "HT_Rat230_PM.na30.annot.csv"))

# HG-U219
scheme.hgu219.na30 <- import.expr.scheme("hgu219", filedir = file.path(scmdir, "na30"),
                       schemefile = file.path(libdir, "HG-U219.cdf"), 
                       probefile  = file.path(libdir, "HG-U219.probe.tab"), 
                       annotfile  = file.path(anndir, "Version09Nov", "HG-U219.na30.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na30 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.1.hg19.probeset.csv"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.hg19.transcript.csv"))

# MoGene-1_0-st-v1.r4: used as exon array
scheme.mogene10stv1.na30 <- import.exon.scheme("mogene10stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_0-st-v1.r4.analysis-lib-files", "MoGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.probeset.csv"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.transcript.csv"))

# RaGene-1_0-st-v1.r4: used as exon array
scheme.ragene10stv1.na30 <- import.exon.scheme("ragene10stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_0-st-v1.r4.analysis-lib-files", "RaGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.probeset.csv"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na30 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na30"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version09Nov", "HuEx-1_0-st-v2.na30.hg19.probeset.csv"),
                          file.path(anndir, "Version09Nov", "HuEx-1_0-st-v2.na30.hg19.transcript.csv"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1.na30 <- import.exon.scheme("moex10stv1",filedir = file.path(scmdir, "na30"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "MoEx_libraryfile", "MoEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version09Nov", "MoEx-1_0-st-v1.na30.mm9.probeset.csv"),
                          file.path(anndir, "Version09Nov", "MoEx-1_0-st-v1.na30.mm9.transcript.csv"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1.na30 <- import.exon.scheme("raex10stv1",filedir=file.path(scmdir, "na30"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.clf"),
                          file.path(libdir, "RaEx_libraryfile", "RaEx-1_0-st-v1.r2.pgf"),
                          file.path(anndir, "Version09Nov", "RaEx-1_0-st-v1.na30.rn4.probeset.csv"),
                          file.path(anndir, "Version09Nov", "RaEx-1_0-st-v1.na30.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na30 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.1.hg19.probeset.csv"),
                            file.path(anndir, "Version09Nov", "HuGene-1_1-st-v1.na30.hg19.transcript.csv"))

# MoGene-1_1-st-v1.r4
scheme.mogene11stv1.na30 <- import.exon.scheme("mogene11stv1", filedir = file.path(scmdir, "na30"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "MoGene-1_1-st-v1.r4.analysis-lib-files", "MoGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.probeset.csv"),
                            file.path(anndir, "Version09Dec", "MoGene-1_1-st-v1.na30.1.mm9.transcript.csv"))

# RaGene-1_1-st-v1.r4
scheme.ragene11stv1.na30 <- import.exon.scheme("ragene11stv1", filedir = scmdir,
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "RaGene-1_1-st-v1.r4.analysis-lib-files", "RaGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.probeset.csv"),
                            file.path(anndir, "Version09Dec", "RaGene-1_1-st-v1.na30.1.rn4.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# E_coli_2:
scheme.ecoli2.na30 <- import.expr.scheme("ecoli2", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "E_coli_2.cdf"),
                      probefile  = file.path(libdir, "E_coli_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version09Nov", "E_coli_2.na30.annot.csv"))

# Ecoli_ASv2:
scheme.ecoliasv2.na30 <- import.expr.scheme("ecoliasv2", filedir = file.path(scmdir, "na30"),
                         schemefile = file.path(libdir, "Ecoli_ASv2.CDF"),
                         probefile  = file.path(libdir, "Ecoli_ASv2.probe.tab"),
                         annotfile  = file.path(anndir, "Version09Nov", "Ecoli_ASv2.na30.annot.csv"))

# Pae_G1a:
scheme.paeg1a.na30 <- import.expr.scheme("paeg1a", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "Pae_G1a.CDF"),
                      probefile  = file.path(libdir, "Pae_G1a.probe.tab"), 
                      annotfile  = file.path(anndir, "Version09Nov", "Pae_G1a.na30.annot.csv"))

# Yeast_2:
scheme.yeast2.na30 <- import.expr.scheme("yeast2", filedir = file.path(scmdir, "na30"),
                      schemefile = file.path(libdir, "Yeast_2.cdf"),
                      probefile  = file.path(libdir, "Yeast_2.probe.tab"),
                      annotfile  = file.path(anndir, "Version09Nov", "Yeast_2.na30.annot.csv"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna10 <- import.expr.scheme("mirna10", filedir = file.path(scmdir, "na30"),
                  schemefile = file.path(libdir, "miRNA-1_0.CDF"),
                  probefile  = file.path(libdir, "miRNA-1_0.probe.tab"))



#------------------------------------------------------------------------------#
# schemes for older Affymetrix annotation files
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na21 <- import.expr.scheme("Scheme_Test3_na21",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version06Nov/Test3.na21.annot.csv",sep="/"))
scheme.test3.na22 <- import.expr.scheme("Scheme_Test3_na22",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07May/Test3.na22.annot.csv",sep="/"))
scheme.test3.na23 <- import.expr.scheme("Scheme_Test3_na23",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07Jul/Test3.na23.annot.csv",sep="/"))
scheme.test3.na24 <- import.expr.scheme("Scheme_Test3_na24",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version07Nov/Test3.na24.annot.csv",sep="/"))
scheme.test3.na25 <- import.expr.scheme("Scheme_Test3_na25",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version08Mar/Test3.na25.annot.csv",sep="/"))
scheme.test3.na26 <- import.expr.scheme("Scheme_Test3_na26",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version08Jul/Test3.na26.annot.csv",sep="/"))
scheme.test3.na27 <- import.expr.scheme("Scheme_Test3_na27",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version08Nov/Test3.na27.annot.csv",sep="/"))
scheme.test3.na29 <- import.expr.scheme("Scheme_Test3_na29",filedir=scmdir,paste(libdir,"Test3.CDF",sep="/"),paste(libdir,"Test3_probe.tab",sep="/"),paste(anndir,"Version09Jul/Test3.na29.annot.csv",sep="/"))

# Hu6800:
scheme.hu6800.na21 <- import.expr.scheme("Scheme_Hu6800_na21",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version06Nov/Hu6800.na21.annot.csv",sep="/"))
scheme.hu6800.na22 <- import.expr.scheme("Scheme_Hu6800_na22",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07May/Hu6800.na22.annot.csv",sep="/"))
scheme.hu6800.na23 <- import.expr.scheme("Scheme_Hu6800_na23",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07Jul/Hu6800.na23.annot.csv",sep="/"))
scheme.hu6800.na24 <- import.expr.scheme("Scheme_Hu6800_na24",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version07Nov/Hu6800.na24.annot.csv",sep="/"))
scheme.hu6800.na25 <- import.expr.scheme("Scheme_Hu6800_na25",filedir=scmdir,paste(libdir,"Hu6800.CDF",sep="/"),paste(libdir,"HuGeneFL_probe.tab",sep="/"),paste(anndir,"Version08Mar/Hu6800.na25.annot.csv",sep="/"))

# HG_U95A:
scheme.hgu95a <- import.expr.scheme("Scheme_HGU95A",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),annotfile="")
scheme.hgu95a.na23 <- import.expr.scheme("Scheme_HGU95A_na23",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG_U95A.na23.annot.csv",sep="/"))
scheme.hgu95a.na24 <- import.expr.scheme("Scheme_HGU95A_na24",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG_U95A.na24.annot.csv",sep="/"))
scheme.hgu95a.na25 <- import.expr.scheme("Scheme_HGU95A_na25",filedir=scmdir,paste(libdir,"HG_U95A.CDF",sep="/"),paste(libdir,"HG-U95A_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG_U95A.na25.annot.csv",sep="/"))

# HG_U95Av2:
scheme.hgu95av2.na21 <- import.expr.scheme("Scheme_HGU95Av2_na21",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG_U95Av2.na21.annot.csv",sep="/"))
scheme.hgu95av2.na22 <- import.expr.scheme("Scheme_HGU95Av2_na22",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07May/HG_U95Av2.na22.annot.csv",sep="/"))
scheme.hgu95av2.na23 <- import.expr.scheme("Scheme_HGU95Av2_na23",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG_U95Av2.na23.annot.csv",sep="/"))
scheme.hgu95av2.na24 <- import.expr.scheme("Scheme_HGU95Av2_na24",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG_U95Av2.na24.annot.csv",sep="/"))
scheme.hgu95av2.na25 <- import.expr.scheme("Scheme_HGU95Av2_na25",filedir=scmdir,paste(libdir,"HG_U95Av2.CDF",sep="/"),paste(libdir,"HG-U95Av2_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG_U95Av2.na25.annot.csv",sep="/"))

# HG-U133A:
scheme.hgu133a.na21 <- import.expr.scheme("Scheme_HGU133A_na21",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133A.na21.annot.csv",sep="/"))
scheme.hgu133a.na22 <- import.expr.scheme("Scheme_HGU133A_na22",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133A.na22.annot.csv",sep="/"))
scheme.hgu133a.na23 <- import.expr.scheme("Scheme_HGU133A_na23",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133A.na23.annot.csv",sep="/"))
scheme.hgu133a.na24 <- import.expr.scheme("Scheme_HGU133A_na24",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG-U133A.na24.annot.csv",sep="/"))
scheme.hgu133a.na25 <- import.expr.scheme("Scheme_HGU133A_na25",filedir=scmdir,paste(libdir,"HG-U133A.CDF",sep="/"),paste(libdir,"HG-U133A_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG-U133A.na25.annot.csv",sep="/"))

# HG-U133B:
scheme.hgu133b.na21 <- import.expr.scheme("Scheme_HGU133B_na21",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133B.na21.annot.csv",sep="/"))
scheme.hgu133b.na22 <- import.expr.scheme("Scheme_HGU133B_na22",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133B.na22.annot.csv",sep="/"))
scheme.hgu133b.na23 <- import.expr.scheme("Scheme_HGU133B_na23",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133B.na23.annot.csv",sep="/"))
scheme.hgu133b.na24 <- import.expr.scheme("Scheme_HGU133B_na24",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG-U133B.na24.annot.csv",sep="/"))
scheme.hgu133b.na25 <- import.expr.scheme("Scheme_HGU133B_na25",filedir=scmdir,paste(libdir,"HG-U133B.CDF",sep="/"),paste(libdir,"HG-U133B_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG-U133B.na25.annot.csv",sep="/"))

# HG-U133_Plus_2:
scheme.hgu133p2.na21 <- import.expr.scheme("Scheme_HGU133p2_na21",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version06Nov/HG-U133_Plus_2.na21.annot.csv",sep="/"))
scheme.hgu133p2.na22 <- import.expr.scheme("Scheme_HGU133p2_na22",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07May/HG-U133_Plus_2.na22.annot.csv",sep="/"))
scheme.hgu133p2.na23 <- import.expr.scheme("Scheme_HGU133p2_na23",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07Jul/HG-U133_Plus_2.na23.annot.csv",sep="/"))
scheme.hgu133p2.na24 <- import.expr.scheme("Scheme_HGU133p2_na24",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version07Nov/HG-U133_Plus_2.na24.annot.csv",sep="/"))
scheme.hgu133p2.na25 <- import.expr.scheme("Scheme_HGU133p2_na25",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version08Mar/HG-U133_Plus_2.na25.annot.csv",sep="/"))
scheme.hgu133p2.na26 <- import.expr.scheme("Scheme_HGU133p2_na26",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version08Jul/HG-U133_Plus_2.na26.annot.csv",sep="/"))
scheme.hgu133p2.na27 <- import.expr.scheme("Scheme_HGU133p2_na27",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version08Nov/HG-U133_Plus_2.na27.annot.csv",sep="/"))
scheme.hgu133p2.na28 <- import.expr.scheme("Scheme_HGU133p2_na28",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version09Mar/HG-U133_Plus_2.na28.annot.csv",sep="/"))
scheme.hgu133p2.na29 <- import.expr.scheme("Scheme_HGU133p2_na29",filedir=scmdir,paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),paste(anndir,"Version09Jul/HG-U133_Plus_2.na29.annot.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133A
scheme.hthgu133A.na28 <- import.expr.scheme("Scheme_HTHGU133A_na28",filedir=scmdir,paste(libdir,"HT_HG-U133A.cdf",sep="/"),paste(libdir,"HT_HG-U133A.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_HG-U133A.na28.annot.csv",sep="/"))

# HT_HG-U133_Plus_PM
scheme.hthgu133ppm.na27 <- import.expr.scheme("Scheme_HTHGU133pPM_na27",filedir=scmdir,paste(libdir,"HT_HG-U133_Plus_PM.CDF",sep="/"),paste(libdir,"HT_HG-U133_Plus_PM.probe.tab",sep="/"),paste(anndir,"Version09Feb/HT_HG-U133_Plus_PM.na27.1.annot.csv",sep="/"))
scheme.hthgu133ppm.na28 <- import.expr.scheme("Scheme_HTHGU133pPM_na28",filedir=scmdir,paste(libdir,"HT_HG-U133_Plus_PM.CDF",sep="/"),paste(libdir,"HT_HG-U133_Plus_PM.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_HG-U133_Plus_PM.na28.annot.csv",sep="/"))

# HT_MG-430_PM
scheme.htmg430pm.na28 <- import.expr.scheme("Scheme_HTMG430PM_na28",filedir=scmdir,paste(libdir,"HT_MG-430_PM.cdf",sep="/"),paste(libdir,"HT_MG-430_PM.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_MG-430_PM.na28.annot.csv",sep="/"))

# HT_Rat230_PM
scheme.htrat230.na28 <- import.expr.scheme("Scheme_HTRat230PM_na28",filedir=scmdir,paste(libdir,"HT_Rat230_PM.cdf",sep="/"),paste(libdir,"HT_Rat230_PM.probe.tab",sep="/"),paste(anndir,"Version09Mar/HT_Rat230_PM.na28.annot.csv",sep="/"))

# HG-U219
scheme.hgu219.na29 <- import.expr.scheme("Scheme_HGU219_na29",filedir=scmdir,paste(libdir,"HG-U219.cdf",sep="/"),paste(libdir,"HG-U219.probe.tab",sep="/"),paste(anndir,"Version09Oct/HG-U219.na29.annot.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for alternative CDF-files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# alternative CDF-files, e.g.:
#AffyProbeMiner
apmdir <- "/Volumes/GigaDrive/Affy/CDF_alternative/AffyProbeMiner"
scheme.hgu133p2.apm <- import.expr.scheme("Scheme_HGU133p2_apm_refseq",filedir=scmdir,paste(apmdir,"HG-U133_Plus_2_transcript_refseq/HG-U133_Plus_2_transcript_refseq.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),"","HG-U133_Plus_2_apm_refseq")
#BrainArray (UniMichigan)
umidir <- "/Volumes/GigaDrive/Affy/CDF_alternative/UniMichigan"
scheme.hgu133p2.umi <- import.expr.scheme("Scheme_HGU133p2_umi_refseq_v10",filedir=scmdir,paste(umidir,"Version10/Hs133P_Hs_REFSEQ/Hs133P_Hs_REFSEQ.cdf",sep="/"),paste(umidir,"Version10/Hs133P_Hs_REFSEQ/Hs133P_Hs_REFSEQ_probe.tab",sep="/"),"","HG-U133_Plus_2_umi_refseq")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r3: used as whole genome array
scheme.hugene10stv1r3.na23 <- import.genome.scheme("Scheme_HuGene10stv1r3_na23",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Jul/HuGene-1_0-st-v1.na23.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na24 <- import.genome.scheme("Scheme_HuGene10stv1r3_na24",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Nov/HuGene-1_0-st-v1.na24.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na25 <- import.genome.scheme("Scheme_HuGene10stv1r3_na25",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Mar/HuGene-1_0-st-v1.na25.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na26 <- import.genome.scheme("Scheme_HuGene10stv1r3_na26",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Jul/HuGene-1_0-st-v1.na26.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na27 <- import.genome.scheme("Scheme_HuGene10stv1r3_na27",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Nov/HuGene-1_0-st-v1.na27.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r3.na28 <- import.genome.scheme("Scheme_HuGene10stv1r3_na28",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r3.analysis_libraryfile/HuGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version09Mar/HuGene-1_0-st-v1.na28.hg18.transcript.csv",sep="/"))

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1r4.na27 <- import.exon.scheme("Scheme_HuGene10stv1r4_na27",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Feb/HuGene-1_0-st-v1.na27.2.hg18.probeset.csv",sep="/"),
                              paste(anndir,"Version09Feb/HuGene-1_0-st-v1.na27.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r4.na28 <- import.exon.scheme("Scheme_HuGene10stv1r4_na28",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Mar/HuGene-1_0-st-v1.na28.hg18.probeset.csv",sep="/"),
                              paste(anndir,"Version09Mar/HuGene-1_0-st-v1.na28.hg18.transcript.csv",sep="/"))

scheme.hugene10stv1r4.na29 <- import.exon.scheme("Scheme_HuGene10stv1r4_na29",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Jul/HuGene-1_0-st-v1.na29.hg18.probeset.csv",sep="/"),
                              paste(anndir,"Version09Jul/HuGene-1_0-st-v1.na29.hg18.transcript.csv",sep="/"))

# note: need to use corrected the probeset annotation file "na30.1"
scheme.hugene10stv1r4.na30 <- import.exon.scheme("Scheme_HuGene10stv1r4_na30_hg19",filedir=scmdir,
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Nov/HuGene-1_1-st-v1.na30.1.hg19.probeset.csv",sep="/"),
                              paste(anndir,"Version09Nov/HuGene-1_1-st-v1.na30.hg19.transcript.csv",sep="/"))

# MoGene-1_0-st-v1.r3:
scheme.mogene10stv1r3.na24 <- import.genome.scheme("Scheme_MoGene10stv1r3_na24",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Nov/MoGene-1_0-st-v1.na24.mm8.transcript.csv",sep="/"))

scheme.mogene10stv1r3.na25 <- import.genome.scheme("Scheme_MoGene10stv1r3_na25",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Mar/MoGene-1_0-st-v1.na25.mm9.transcript.csv",sep="/"))

scheme.mogene10stv1r3.na27 <- import.genome.scheme("Scheme_MoGene10stv1r3_na27",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r3.analysis-lib-files/MoGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Nov/MoGene-1_0-st-v1.na27.mm9.transcript.csv",sep="/"))

# MoGene-1_0-st-v1.r4: used as exon array
# need to fix problem with annotation files first, see:
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2009-August/029049.html
source(paste(path.package("xps"),"examples/updateAnnotation.R",sep="/"))
updateAnnotation("MoGene-1_0-st-v1.na29.mm9.probeset.csv", "MoGene-1_0-st-v1.na29.mm9.probeset.fixed.csv", probeset="10338063", skip=18, eol="\n")
updateAnnotation("MoGene-1_0-st-v1.na29.mm9.transcript.csv", "MoGene-1_0-st-v1.na29.mm9.transcript.fixed.csv", probeset="10338063", skip=19, eol="\n")

scheme.mogene10stv1r4.na29 <- import.exon.scheme("Scheme_MoGene10stv1r4_na29",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Jul/MoGene-1_0-st-v1.na29.mm9.probeset.fixed.csv",sep="/"),
                              paste(anndir,"Version09Jul/MoGene-1_0-st-v1.na29.mm9.transcript.fixed.csv",sep="/"))

# in Dec 2009 Affymetrix has corrected the annotation files for release "na30.1"
scheme.mogene10stv1r4.na30 <- import.exon.scheme("Scheme_MoGene10stv1r4_na30_1",filedir=scmdir,
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"MoGene-1_0-st-v1.r4.analysis-lib-files/MoGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Dec/MoGene-1_1-st-v1.na30.1.mm9.probeset.csv",sep="/"),
                              paste(anndir,"Version09Dec/MoGene-1_1-st-v1.na30.1.mm9.transcript.csv",sep="/"))


# RaGene-1_0-st-v1.r3:
scheme.ragene10stv1r3.na24 <- import.genome.scheme("Scheme_RaGene10stv1r3_na24",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version07Nov/RaGene-1_0-st-v1.na24.rn4.transcript.csv",sep="/"))

scheme.ragene10stv1r3.na25 <- import.genome.scheme("Scheme_RaGene10stv1r3_na25",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r3.analysis-lib-files/RaGene-1_0-st-v1.r3.pgf",sep="/"),
                              paste(anndir,"Version08Mar/RaGene-1_0-st-v1.na25.rn4.transcript.csv",sep="/"))

# RaGene-1_0-st-v1.r4: used as exon array
# need to fix problem with annotation files first, see:
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2009-August/029049.html
source(paste(path.package("xps"),"examples/updateAnnotation.R",sep="/"))
updateAnnotation("RaGene-1_0-st-v1.na29.rn4.probeset.csv", "RaGene-1_0-st-v1.na29.rn4.probeset.fixed.csv", probeset="10700063", skip=18, eol="\n")
updateAnnotation("RaGene-1_0-st-v1.na29.rn4.transcript.csv", "RaGene-1_0-st-v1.na29.rn4.transcript.fixed.csv", probeset="10700063", skip=19, eol="\n")

scheme.ragene10stv1r4.na29 <- import.exon.scheme("Scheme_RaGene10stv1r4_na29",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Jul/RaGene-1_0-st-v1.na29.rn4.probeset.fixed.csv",sep="/"),
                              paste(anndir,"Version09Jul/RaGene-1_0-st-v1.na29.rn4.transcript.fixed.csv",sep="/"))

# in Dec 2009 Affymetrix has corrected the annotation files for release "na30.1"
scheme.ragene10stv1r4.na30 <- import.exon.scheme("Scheme_RaGene10stv1r4_na30_1",filedir=scmdir,
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.clf",sep="/"),
                              paste(libdir,"RaGene-1_0-st-v1.r4.analysis-lib-files/RaGene-1_0-st-v1.r4.pgf",sep="/"),
                              paste(anndir,"Version09Dec/RaGene-1_1-st-v1.na30.1.rn4.probeset.csv",sep="/"),
                              paste(anndir,"Version09Dec/RaGene-1_1-st-v1.na30.1.rn4.transcript.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2r2.na21 <- import.exon.scheme("Scheme_HuEx10stv2r2_na21",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version06Nov/HuEx-1_0-st-v2.na21.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version06Nov/HuEx-1_0-st-v2.na21.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na22 <- import.exon.scheme("Scheme_HuEx10stv2r2_na22",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version07Mar/HuEx-1_0-st-v2.na22.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version07Mar/HuEx-1_0-st-v2.na22.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na23 <- import.exon.scheme("Scheme_HuEx10stv2r2_na23",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version07Jul/HuEx-1_0-st-v2.na23.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version07Jul/HuEx-1_0-st-v2.na23.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na24 <- import.exon.scheme("Scheme_HuEx10stv2r2_na24",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version07Nov/HuEx-1_0-st-v2.na24.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version07Nov/HuEx-1_0-st-v2.na24.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na25 <- import.exon.scheme("Scheme_HuEx10stv2r2_na25",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version08Mar/HuEx-1_0-st-v2.na25.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version08Mar/HuEx-1_0-st-v2.na25.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na26 <- import.exon.scheme("Scheme_HuEx10stv2r2_na26",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version08Jul/HuEx-1_0-st-v2.na26.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version08Jul/HuEx-1_0-st-v2.na26.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na27 <- import.exon.scheme("Scheme_HuEx10stv2r2_na27",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version08Nov/HuEx-1_0-st-v2.na27.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version08Nov/HuEx-1_0-st-v2.na27.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na28 <- import.exon.scheme("Scheme_HuEx10stv2r2_na28",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version09Mar/HuEx-1_0-st-v2.na28.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version09Mar/HuEx-1_0-st-v2.na28.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na29 <- import.exon.scheme("Scheme_HuEx10stv2r2_na29",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na29.hg18.probeset.csv",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na29.hg18.transcript.csv",sep="/"))

scheme.huex10stv2r2.na30 <- import.exon.scheme("Scheme_HuEx10stv2r2_na30",filedir=scmdir,
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                            paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na30.hg19.probeset.csv",sep="/"),
                            paste(anndir,"Version09Jul/HuEx-1_0-st-v2.na30.hg19.transcript.csv",sep="/"))

# HuEx-1_0-st-v2.r2 old annotation:
scheme.huex10stv2r2.old <- import.exon.scheme("Scheme_HuEx10stv2r2_old",filedir=scmdir,
                           paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
                           paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
                           paste(anndir,"Version06Jun/HuEx-1_0-st-probeset-annot.csv",sep="/"),
                           paste(anndir,"Version06Jun/HuEx-1_0-st-transcript-annot.csv",sep="/"),
                           paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.controls.ps",sep="/"))

# MoEx-1_0-st-v1.r2:
scheme.moex10stv1r2.na21 <- import.exon.scheme("Scheme_MoEx10stv2r2_na21",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version06Nov/MoEx-1_0-st-v1.na21.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version06Nov/MoEx-1_0-st-v1.na21.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na22 <- import.exon.scheme("Scheme_MoEx10stv2r2_na22",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07May/MoEx-1_0-st-v1.na22.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version07May/MoEx-1_0-st-v1.na22.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na23 <- import.exon.scheme("Scheme_MoEx10stv2r2_na23",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Jul/MoEx-1_0-st-v1.na23.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version07Jul/MoEx-1_0-st-v1.na23.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na24 <- import.exon.scheme("Scheme_MoEx10stv2r2_na24",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Nov/MoEx-1_0-st-v1.na24.mm8.probeset.csv",sep="/"),
                            paste(anndir,"Version07Nov/MoEx-1_0-st-v1.na24.mm8.transcript.csv",sep="/"))

scheme.moex10stv1r2.na25 <- import.exon.scheme("Scheme_MoEx10stv2r2_na25",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version08Mar/MoEx-1_0-st-v1.na25.mm9.probeset.csv",sep="/"),
                            paste(anndir,"Version08Mar/MoEx-1_0-st-v1.na25.mm9.transcript.csv",sep="/"))

scheme.moex10stv1r2.na30 <- import.exon.scheme("Scheme_MoEx10stv2r2_na30",filedir=scmdir,
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"MoEx_libraryfile/MoEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version09Nov/MoEx-1_0-st-v1.na30.mm9.probeset.csv",sep="/"),
                            paste(anndir,"Version09Nov/MoEx-1_0-st-v1.na30.mm9.transcript.csv",sep="/"))

# RaEx-1_0-st-v1.r2:
scheme.raex10stv1r2.na21 <- import.exon.scheme("Scheme_RaEx10stv2r2_na21",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version06Nov/RaEx-1_0-st-v1.na21.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version06Nov/RaEx-1_0-st-v1.na21.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na22 <- import.exon.scheme("Scheme_RaEx10stv2r2_na22",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07May/RaEx-1_0-st-v1.na22.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version07May/RaEx-1_0-st-v1.na22.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na23 <- import.exon.scheme("Scheme_RaEx10stv2r2_na23",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Jul/RaEx-1_0-st-v1.na23.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version07Jul/RaEx-1_0-st-v1.na23.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na24 <- import.exon.scheme("Scheme_RaEx10stv2r2_na24",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version07Nov/RaEx-1_0-st-v1.na24.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version07Nov/RaEx-1_0-st-v1.na24.rn4.transcript.csv",sep="/"))

scheme.raex10stv1r2.na25 <- import.exon.scheme("Scheme_RaEx10stv2r2_na25",filedir=scmdir,
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.clf",sep="/"),
                            paste(libdir,"RaEx_libraryfile/RaEx-1_0-st-v1.r2.pgf",sep="/"),
                            paste(anndir,"Version08Mar/RaEx-1_0-st-v1.na25.rn4.probeset.csv",sep="/"),
                            paste(anndir,"Version08Mar/RaEx-1_0-st-v1.na25.rn4.transcript.csv",sep="/"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for other Affymetrix arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Yeast_2:
scheme.yeast2.na28 <- import.expr.scheme("Scheme_Yeast2_na27",filedir=scmdir,
                      paste(libdir,"Yeast_2.cdf",sep="/"),
                      paste(libdir,"Yeast_2.probe.tab",sep="/"),
                      paste(anndir,"Version09Mar/Yeast_2.na28.annot.csv",sep="/"))

# Pae_G1a:
scheme.paeg1a.na28 <- import.expr.scheme("Scheme_PaeG1a_na28",filedir=scmdir,
                      paste(libdir,"Pae_G1a.CDF",sep="/"),
                      paste(libdir,"Pae_G1a.probe.tab",sep="/"), 
                      paste(anndir,"Version09Mar/Pae_G1a.na28.annot.csv",sep="/"))

# Mouse430_2:
scheme.mouse430.2.na28 <- import.expr.scheme("Scheme_Mouse430_2_na28",filedir=scmdir,
                          paste(libdir,"Mouse430_2.cdf",sep="/"),
                          paste(libdir,"Mouse430_2.probe.tab",sep="/"),
                          paste(anndir,"Version09Mar/Mouse430_2.na28.annot.csv",sep="/"))

# miRNA-1_0:
# note: you need to rename "miRNA-1_0.probe_list.20081203.txt" to "miRNA-1_0.probe.tab"
scheme.mirna <- import.expr.scheme("Scheme_miRNA_1_0",filedir=scmdir,
                paste(libdir,"miRNA-1_0.CDF",sep="/"),
                paste(libdir,"miRNA-1_0.probe.tab",sep="/"))



#------------------------------------------------------------------------------#
