// File created: 08/02/2008                          last modified: 08/22/2010
// Author: Christian Stratowa 08/02/2008

/*
 *******************************************************************************
 *
 * NOTE:
 * This macro is called by method "root.profilechart" of package "xps".
 * It allows to open a ROOT canvas to display a parallel coordinates plot
 *
 *******************************************************************************
 */

void macroDrawProfilePlot(const char *xps, const char *filename, const char *canvasname,
                          const char *setname, const char *treename, const char *exten,
                          const char *varlist, const char *savename, 
                          Double_t min, Double_t max, Int_t aslog, Int_t globalscale,
                          Int_t candle, Int_t width, Int_t height)
{
   gSystem->Load("libGui");    //necessary for ROOT version >=5.15/09
   gSystem->Load(xps);

// create new manager
   XPlot *plotter = new XPlot("Plot");

// open root file
   plotter->Open(filename);

// add trees
   TString tname = setname;
   if (strcmp(treename,"*") == 0) {
      tname += "/"; tname += treename;
      tname += "."; tname += exten;

      plotter->AddTree(tname.Data());
   } else {
      char    *name  = new char[strlen(treename) + 1];
      char    *dname = name;
      int      num   = TString(treename).CountChar(':') + 1;
      int      count = 0;
      TString *trees = new TString[num];

      name = strtok(strcpy(name, treename), ":");
      while(name) {
         tname  = setname;
         tname += "/"; tname += name;
         tname += "."; tname += exten;
         trees[count++] = tname;

         name = strtok(0, ":");
         if (name == 0) break;
      }//while

      for (int i=0; i<num; i++) {
         plotter->AddTree(trees[i].Data());
      }//for

      delete [] trees;
      delete [] dname;
   }//if

   plotter->NewCanvas(canvasname, "title", 20, 20, width, height);

// draw in separate canvases
   plotter->DrawParallelCoord("", varlist, min, max, aslog, globalscale, candle);

// save figure and close canvas
   if (strcmp(savename, "") != 0) {
      plotter->SaveAs(savename);
      plotter->CloseCanvas("");
   }//if

// cleanup
   delete plotter;
}//macroDrawProfilePlot
