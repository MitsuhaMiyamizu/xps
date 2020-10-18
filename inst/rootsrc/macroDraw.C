// File created: 10/14/2007                          last modified: 09/05/2008
// Author: Christian Stratowa 10/13/2007

/*
 *******************************************************************************
 *
 * NOTE:
 * This macro is called by method "root.draw" of package "xps".
 *
 *******************************************************************************
 */

void macroDraw(const char *xps, const char *filename, const char *canvasname,
               Int_t ntree, const char *treename1, const char *treename2,
               const char *treename3, const char *varlist, const char *logbase,
               const char *type, const char *option, const char *savename, 
               Int_t width, Int_t height)
{
   gSystem->Load("libGui"); //necessary for ROOT version >=5.15/09
   gSystem->Load(xps);

// create new manager
   XPlot *plotter = new XPlot("Plot");

// open root file
   plotter->Open(filename);

   plotter->NewCanvas(canvasname, "title", 20, 20, width, height);

// draw
   if (ntree == 1) {
      plotter->Draw("", treename1, varlist, logbase, type,option);
   } else if (ntree == 2) {
      plotter->Draw("", treename1, treename2, varlist, logbase, type,option);
   } else if (ntree == 3) {
      plotter->Draw("", treename1, treename2, treename3, varlist, logbase, type,option);
   }//if

// save figure and close canvas
   if (strcmp(savename, "") != 0) {
      plotter->SaveAs(savename);
      plotter->CloseCanvas("");
   }//if

// cleanup
   delete plotter;
}//macroDraw
