// File created: 10/13/2007                          last modified: 09/05/2008
// Author: Christian Stratowa 10/13/2007

/*
 *******************************************************************************
 *
 * NOTE:
 * This macro is called by method "root.image" of package "xps".
 * It allows to open a ROOT canvas to display the DNA-chip image
 *
 *******************************************************************************
 */

void macroDrawImage(const char *xps, const char *filename, const char *canvasname,
                    const char *treename, const char *varlist, const char *logbase,
                    const char *option, const char *savename, Double_t min,
                    Double_t max, Int_t width, Int_t height)
{
   gSystem->Load("libGui"); //necessary for ROOT version >=5.15/09
   gSystem->Load(xps);

// create new manager
   XPlot *plotter = new XPlot("Plot");

// open root file
   plotter->Open(filename);

   plotter->NewCanvas(canvasname, "title", 20, 20, width, height);

// draw in separate canvases
//   plotter->DrawImage("", treename, varlist, logbase, option, 0, 0, "U");
   plotter->DrawImage("", treename, varlist, logbase, option, min, max, "U");

// save figure and close canvas
   if (strcmp(savename, "") != 0) {
      plotter->SaveAs(savename);
      plotter->CloseCanvas("");
   }//if

// cleanup
   delete plotter;
}//macroDrawImage
