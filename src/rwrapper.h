/*
 *******************************************************************************
 *  An R wrapper for the XPS libraries
 *
 * File: rwrapper.h
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
*/

#ifndef R_WRAPPER
#define R_WRAPPER

extern "C" {
void ImportExprSchemes(char **filename, char **dirname, char **chipname,
                       char **scheme, char **probe, char **annot,
                       int *verbose, int *err);
void ImportExonSchemes(char **filename, char **dirname, char **chipname,
                       char **layout, char **scheme, char **probeset,
                       char **transcript, char **control,
                       int *verbose, int *err);
void ImportGenomeSchemes(char **filename, char **dirname, char **chipname,
                         char **layout, char **scheme, char **transcript,
                         int *verbose, int *err);
void ExportScheme(char **filename, char **treename, char **varlist,
                  char **outfile, char **sep, int *verbose);

void ImportData(char **filename, char **dirname, char **chiptype,
                char **chipname, char **schemefile, char **treeset,
                char **celfiles, char **celnames, int *numdata,
                char **project, int *nproject, char **author, int *nauthor,
                char **dataset, int *ndataset, char **source, int *nsource,
                char **sample, int *nsample, char **cell, int *ncell,
                char **pcell, int *npcell, char **tissue, int *ntissue,
                char **biopsy, int *nbiopsy, char **array, int *narray,
                char **hyb, int *nhyb, char **treat, int *ntreat, 
                int *replace, int *update, int *verbose, char **result);
void ExportData(char **filename, char **schemefile, char **chiptype,
                char **datatype, char **treenames, int *ntrees, char **exten,
                char **varlist, char **outfile, char **sep,
                int *verbose, int *err);

void PreprocessRMA(char **filename, char **dirname, char **chipname,
                   char **chiptype, char **schemefile, char **tmpdir,
                   char **bgrdoption, char **exproption, char **treeset, char **datafile, 
                   char **treenames, int *ntrees, int *normalize, double *pars,
                   int *bgrdlevel, int *normlevel, int *exprlevel,
                   int *verbose, char **result);
void PreprocessMAS4(char **filename, char **dirname, char **chipname,
                    char **chiptype, char **schemefile, char **tmpdir,
                    char **exproption, char **treeset, char **treenames,
                    int *ntrees, int *bgrdlevel, int *exprlevel,
                    int *verbose, char **result);
void PreprocessMAS5(char **filename, char **dirname, char **chipname,
                    char **chiptype, char **schemefile, char **tmpdir,
                    char **exproption, char **treeset, char **datafile,
                    char **treenames, int *ntrees, int *bgrdlevel, int *exprlevel,
                    int *verbose, char **result);

void PreprocessFIRMA(char **filename, char **dirname, char **chipname,
                     char **chiptype, char **schemefile, char **tmpdir,
                     char **bgrdoption, char **exproption, char **treeset, char **datafile, 
                     char **treenames, int *ntrees, int *normalize, double *pars,
                     int *bgrdlevel, int *normlevel, int *exprlevel,
                     int *verbose, char **result);

void PreprocessMAS5Call(char **filename, char **dirname, char **chipname,
                        char **chiptype, char **schemefile, char **tmpdir,
                        char **calloption, char **treeset, char **datafile,
                        char **treenames, int *ntrees, double *tau,
                        double *alpha1, double *alpha2, int *ignore,
                        char **bgrdoption, int *bgrdlevel, int *callevel,
                        int *verbose, char **result);
void PreprocessDABGCall(char **filename, char **dirname, char **chipname,
                        char **chiptype, char **schemefile, char **calloption,
                        char **treeset, char **datafile, char **treenames,
                        int *ntrees, double *alpha1, double *alpha2, int *level,
                        int *verbose, char **result);
void PreprocessINICall(char **filename, char **dirname, char **chipname,
                       char **chiptype, char **schemefile, char **tmpdir,
                       char **option, char **treeset, char **datafile, char **treenames,
                       int *ntrees, int *version, double *weight, double *mu,
                       double *scale, double *tol, int *cyc, double *alpha1,
                       double *alpha2, int *normlevel, int *callevel,
                       int *verbose, char **result);

void Preprocess(char **filename, char **dirname, char **chipname, char **chiptype,
                char **schemefile, char **tmpdir, int *update,
                char **bgrdtype, char **bgrdselection, char **bgrdoption,
                int *nbgrdpar, double *bgrdpars,
                char **normtype, char **normselection, char **normoption,
                int *nnormpar, double *normpars,
                char **algorithm, char **sumtype, char **sumselection, char **sumoption,
                int *nsumpar, double *sumpars,
                char **reftree, char **refmethod, double *refparam,
                char **treeset, char **datafile, char **treenames, int *ntrees,
                int *bgrdlevel, int *normlevel, int *sumlevel,
                int *bufsize, int *verbose, char **result);
void BgCorrect(char **filename, char **dirname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption,
               char **bgrdtype, char **bgrdoption, int *npars, double *pars,
               char **treeset, char **treenames, int *ntrees,
               int *update, int *level, int *verbose, char **result);
void Normalize(char **filename, char **dirname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption,
               char **type, char **option, int *npar, double *pars,
               int *level, char **treeset, char **treenames, int *ntrees,
               char **reftree, char **refmethod, int *update,
               int *verbose, char **result);
void Summarize(char **filename, char **dirname, char **chipname, char **chiptype,
               char **schemefile, char **tmpdir, char **seloption, char **algorithm,
               char **type, char **option, int *npar, double *pars, int *level,
               char **treeset, char **treenames, int *ntrees,
               int *update, int *verbose, char **result);

void Normxpress(char **filename, char **dirname, char **chiptype,
                char **schemefile, char **tmpdir, char **seloption, double *pc,
                char **type, char **option, int *npar, double *pars,
                int *level, char **treeset, char **datafile, char **treenames, int *ntrees,
                char **reftree, char **refmethod, int *update,
                int *verbose, char **result);

void PreFilter(char **filename, char **dirname, char **chiptype, char **chipname,
               char **schemefile, char **treeset, char **treename,
               int *min, char **logbase,
               int *nvarpar, double *varpars, char **varoption,
               int *nlowpar, double *lowpars,  char **lowoption,
               int *nhighpar, double *highpars, char **highoption,
               int *nquanpar, double *quanpars, char **quanoption,
               int *ncallpar, double *callpars, char **calloption,
               char **exprtrees, int *nexpr, char **calltrees, int *ncall,
               int *verbose, char **result);
void UniFilter(char **filename, char **dirname, char **chiptype, char **chipname,
               char **schemefile, char **treeset, char **treename,
               int *min, char **logbase,
               char **unitest, int *nunipar, double *unipars, char **unioption,
               int *nfcpar, double *fcpars,  char **fcoption,
               int *nufpar, double *ufpars,  char **ufoption,
               int *ncallpar, double *callpars, char **calloption,
               char **exprtrees, int *nexpr, char **calltrees, int *ncall,
               char **group, int *grpidx, char **fltrtree, int *nfltr,
               char **varlist, int *update, int *verbose, char **result);

void ChipNameType(char **filename, char **nametype);
void GeneChipProbeInfo(char **filename, int *value);
void GetNumberOfTrees(char **filename, char **setname, char **exten, int *numtrees);
void GetNumberOfTrees4Exten(char **filename, char **exten, int *numtrees);
void GetTreeNames(char **filename, char **setname, char **exten, int *gettitle,
                  char **treenames);
void GetTreeNames4Exten(char **filename, char **exten, int *gettitle,
                        char **treenames);
void GetRawCELNames(char **datafile, int *numtrees, char **treename, char **celname);
void MetaProbesets(char **schemefile, char **infile, char **outfile,
                   int *level, int *meta, int *err);

void PlotImage(char **filename, char **canvasname, char **treename,
               char **varlist, char **logbase, char **option, int *err);
}

#endif
