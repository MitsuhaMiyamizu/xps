// File created: 10/16/2007                          last modified: 09/05/2008
// Author: Christian Stratowa 10/16/2007

/*
 *******************************************************************************
 *
 * NOTE:
 * This macro is called by methods "root.density" of package "xps".
 *
 *******************************************************************************
 */

void macroDrawDensity(const char *xps, const char *filename, const char *canvasname,
                      const char *treename, const char *varlist, const char *logbase,
                      const char *option, const char *savename, Int_t width, Int_t height)
{
   gSystem->Load("libGui"); //necessary for ROOT version >=5.15/09
   gSystem->Load(xps);

// create new manager
   XPlot *plotter = new XPlot("Plot");

// open root file
   plotter->Open(filename);

// add trees
   plotter->AddTree(treename);

   plotter->NewCanvas(canvasname, "title", 20, 20, width, height);

// set attributes
//   Int_t colors[] = {1,2,3,4};
   Int_t colors[] = {kBlack,kBlue,kRed,kOrange};
//   Int_t styles[] = {1,2,3,4};
   plotter->SetLineColor(4,colors,0);
//   plotter->SetLineStyle(4,styles,2);

// draw density
   plotter->DrawDensity("", varlist, logbase, "epanechnikov", option, 512, 0, 0);

// save figure and close canvas
   if (strcmp(savename, "") != 0) {
      plotter->SaveAs(savename);
      plotter->CloseCanvas("");
   }//if

// cleanup
   delete plotter;
}//macroDrawDensity
