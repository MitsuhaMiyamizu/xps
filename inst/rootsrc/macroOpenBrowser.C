// File created: 04/11/2007                          last modified: 07/12/2007
// Author: Christian Stratowa 04/11/2007

/*
 *******************************************************************************
 *
 * NOTE:
 * This macro is called by method "root.browser" of S4 class TreeSet.
 * It allows to open the ROOT browser from within R to inspect the ROOT files
 * and their content.
 *
 *******************************************************************************
 */

void macroOpenBrowser(const char *xps)
{
   gSystem->Load("libGui"); //necessary for ROOT version >=5.15/09
   gSystem->Load("libTreeViewer.so");  //necessary for TParallelCoord
   gSystem->Load(xps);
   gROOT->ProcessLine("TBrowser b;");
}
