#==============================================================================#
# Constructors.R: functions serving as class contructors
# Note: suggested by Martin Morgan to save the user from calling new()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ProjectInfo:
# PreFilter:
# UniFilter:
#==============================================================================#

"ProjectInfo" <-
function(submitter      = character(),
         laboratory     = character(),
         contact        = character(),
         project        = character(),
         author         = character(),
         dataset        = character(),
         source         = character(),
         sample         = character(),
         celline        = character(),
         primarycell    = character(),
         tissue         = character(),
         biopsy         = character(),
         arraytype      = character(),
         hybridizations = character(),
         treatments     = character()) 
{
   if (debug.xps()) print("------ProjectInfo constructor------")

   info <- new("ProjectInfo",
               submitter   = submitter,
               laboratory  = laboratory,
               contact     = contact,
               project     = as.list(project),
               author      = as.list(author),
               dataset     = as.list(dataset),
               source      = as.list(source),
               sample      = as.list(sample),
               celline     = as.list(celline),
               primarycell = as.list(primarycell),
               tissue      = as.list(tissue),
               biopsy      = as.list(biopsy),
               arraytype   = as.list(arraytype));

   if (length(hybridizations) > 0) {
      hybridizInfo(info) <- hybridizations;
   }#if

   if (length(treatments) > 0) {
      treatmentInfo(info) <- treatments;
   }#if

   return(info);
}#ProjectInfo constructor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"PreFilter" <-
function(mad         = character(),
         cv          = character(),
         variance    = character(),
         difference  = character(),
         ratio       = character(),
         gap         = character(),
         lothreshold = character(),
         hithreshold = character(),
         quantile    = character(),
         prescall    = character()) 
{
   if (debug.xps()) print("------PreFilter constructor------")

   fltr <- new("PreFilter",
               mad         = as.list(mad),
               cv          = as.list(cv),
               variance    = as.list(variance),
               difference  = as.list(difference),
               ratio       = as.list(ratio),
               gap         = as.list(gap),
               lothreshold = as.list(lothreshold),
               hithreshold = as.list(hithreshold),
               quantile    = as.list(quantile),
               prescall    = as.list(prescall))

   return(fltr);
}#PreFilter constructor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"UniFilter" <-
function(unitest    = "t.test",
         foldchange = character(),
         prescall   = character(),
         unifilter  = character()) 
{
   if (debug.xps()) print("------UniFilter constructor------")

   fltr <- new("UniFilter",
               foldchange = as.list(foldchange),
               prescall   = as.list(prescall),
               unifilter  = as.list(unifilter),
               unitest    = as.list(unitest));

   return(fltr);
}#UniFilter constructor


#------------------------------------------------------------------------------#
