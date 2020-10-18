#==============================================================================#
# TreeSetClasses.R: contains all class definitions and method definitions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# methods.ProjectInfo.R:     contains methods for ProjectInfo
# methods.TreeSet.R:         contains methods for TreeSet
# methods.SchemeTreeSet.R:   contains methods for SchemeTreeSet
# methods.ProcesSet.R:       contains methods for ProcesSet
# methods.DataTreeSet.R:     contains methods for DataTreeSet
# methods.ExprTreeSet.R:     contains methods for ExprTreeSet
# methods.CallTreeSet.R:     contains methods for CallTreeSet
# methods.QualTreeSet.R:     contains methods for QualTreeSet
# methods.Filter.R:          contains methods for Filter
# methods.PreFilter.R:       contains methods for PreFilter
# methods.UniFilter.R:       contains methods for UniFilter
# methods.FilterTreeSet.R:   contains methods for FilterTreeSet
# methods.AnalysisTreeSet.R: contains methods for AnalysisTreeSet
#==============================================================================#


#------------------------------------------------------------------------------#
# ProjectInfo:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class ProjectInfo
setClass("ProjectInfo",
   representation(submitter      = "character",
                  laboratory     = "character",
                  contact        = "character",
                  project        = "list",
                  author         = "list",
                  dataset        = "list",
                  source         = "list",
                  sample         = "list",
                  celline        = "list",
                  primarycell    = "list",
                  tissue         = "list",
                  biopsy         = "list",
                  arraytype      = "list",
                  hybridizations = "data.frame",
                  treatments     = "data.frame"
   ),
   prototype(submitter      = "",
             laboratory     = "",
             contact        = "",
             project        = list(),
             author         = list(),
             dataset        = list(),
             source         = list(),
             sample         = list(),
             celline        = list(),
             primarycell    = list(),
             tissue         = list(),
             biopsy         = list(),
             arraytype      = list(),
             hybridizations = data.frame(matrix(nr=0,nc=0)),
             treatments     = data.frame(matrix(nr=0,nc=0))
   )
)#ProjectInfo

# generic methods for class ProjectInfo
setGeneric("projectInfo",     function(object)        standardGeneric("projectInfo"));
setGeneric("projectInfo<-",   function(object, value) standardGeneric("projectInfo<-"));
setGeneric("authorInfo",      function(object)        standardGeneric("authorInfo"));
setGeneric("authorInfo<-",    function(object, value) standardGeneric("authorInfo<-"));
setGeneric("datasetInfo",     function(object)        standardGeneric("datasetInfo"));
setGeneric("datasetInfo<-",   function(object, value) standardGeneric("datasetInfo<-"));
setGeneric("sourceInfo",      function(object)        standardGeneric("sourceInfo"));
setGeneric("sourceInfo<-",    function(object, value) standardGeneric("sourceInfo<-"));
setGeneric("sampleInfo",      function(object)        standardGeneric("sampleInfo"));
setGeneric("sampleInfo<-",    function(object, value) standardGeneric("sampleInfo<-"));
setGeneric("cellineInfo",     function(object)        standardGeneric("cellineInfo"));
setGeneric("cellineInfo<-",   function(object, value) standardGeneric("cellineInfo<-"));
setGeneric("primcellInfo",    function(object)        standardGeneric("primcellInfo"));
setGeneric("primcellInfo<-",  function(object, value) standardGeneric("primcellInfo<-"));
setGeneric("tissueInfo",      function(object)        standardGeneric("tissueInfo"));
setGeneric("tissueInfo<-",    function(object, value) standardGeneric("tissueInfo<-"));
setGeneric("biopsyInfo",      function(object)        standardGeneric("biopsyInfo"));
setGeneric("biopsyInfo<-",    function(object, value) standardGeneric("biopsyInfo<-"));
setGeneric("arrayInfo",       function(object)        standardGeneric("arrayInfo"));
setGeneric("arrayInfo<-",     function(object, value) standardGeneric("arrayInfo<-"));
setGeneric("hybridizInfo",    function(object)        standardGeneric("hybridizInfo"));
setGeneric("hybridizInfo<-",  function(object, value) standardGeneric("hybridizInfo<-"));
setGeneric("treatmentInfo",   function(object)        standardGeneric("treatmentInfo"));
setGeneric("treatmentInfo<-", function(object, value) standardGeneric("treatmentInfo<-"));


#------------------------------------------------------------------------------#
# TreeSet: virtual superset for 'SchemeTreeSet', 'DataTreeSet' etc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class TreeSet
setClass("TreeSet",
   representation(setname   = "character",
                  settype   = "character",
                  rootfile  = "character",
                  filedir   = "character",
                  numtrees  = "numeric",
                  treenames = "list",
                  "VIRTUAL"
   ),
   prototype(setname   = "",
             settype   = "",
             rootfile  = "ROOTFile",
             filedir   = getwd(),
             numtrees  = 0,
             treenames = list()
   )
)#TreeSet

# generic methods for class TreeSet
setGeneric("rootFile",     function(object)        standardGeneric("rootFile"));
setGeneric("rootFile<-",   function(object, value) standardGeneric("rootFile<-"));
setGeneric("fileDir",      function(object)        standardGeneric("fileDir"));
setGeneric("fileDir<-",    function(object, value) standardGeneric("fileDir<-"));
setGeneric("setName",      function(object)        standardGeneric("setName"));
setGeneric("setName<-",    function(object, value) standardGeneric("setName<-"));
setGeneric("setType",      function(object)        standardGeneric("setType"));
setGeneric("setType<-",    function(object, value) standardGeneric("setType<-"));
setGeneric("treeNames",    function(object)        standardGeneric("treeNames"));
setGeneric("treeInfo",     function(object,...)    standardGeneric("treeInfo"));
setGeneric("export",       function(object,...)    standardGeneric("export"));
setGeneric("root.browser", function(object)        standardGeneric("root.browser"));


#------------------------------------------------------------------------------#
# SchemeTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("SchemeTreeSet",
   representation(chipname  = "character",
                  chiptype  = "character",
                  probeinfo = "list",
                  probe     = "data.frame",
                  unitname  = "data.frame",
                  mask      = "data.frame"
   ),
   contains=c("TreeSet"),
   prototype(chipname  = "",
             chiptype  = "GeneChip",
             probeinfo = list(),
             probe     = data.frame(matrix(nr=0,nc=0)),
             unitname  = data.frame(matrix(nr=0,nc=0)),
             mask      = data.frame(matrix(nr=0,nc=0))
   )
)#SchemeTreeSet

# generic methods for class SchemeTreeSet
setGeneric("chipName",             function(object)        standardGeneric("chipName"));
setGeneric("chipType",             function(object)        standardGeneric("chipType"));
setGeneric("chipType<-",           function(object, value) standardGeneric("chipType<-"));
setGeneric("chipMask",             function(object)        standardGeneric("chipMask"));
setGeneric("chipMask<-",           function(object, value) standardGeneric("chipMask<-"));
setGeneric("chipProbe",            function(object)        standardGeneric("chipProbe"));
setGeneric("chipProbe<-",          function(object, value) standardGeneric("chipProbe<-"));
setGeneric("unitNames",            function(object, ...)   standardGeneric("unitNames"));
setGeneric("unitNames<-",          function(object, value) standardGeneric("unitNames<-"));
setGeneric("probeContentGC",       function(object, ...)   standardGeneric("probeContentGC"));
setGeneric("probeSequence",        function(object, ...)   standardGeneric("probeSequence"));
setGeneric("probeInfo",            function(object)        standardGeneric("probeInfo"));
setGeneric("nrows",                function(object)        standardGeneric("nrows"));
setGeneric("ncols",                function(object)        standardGeneric("ncols"));
setGeneric("attachMask",           function(object)        standardGeneric("attachMask"));
setGeneric("removeMask",           function(object)        standardGeneric("removeMask"));
setGeneric("attachProbe",          function(object, ...)   standardGeneric("attachProbe"));
setGeneric("removeProbe",          function(object)        standardGeneric("removeProbe"));
setGeneric("attachProbeContentGC", function(object)        standardGeneric("attachProbeContentGC"));
setGeneric("removeProbeContentGC", function(object)        standardGeneric("removeProbeContentGC"));
setGeneric("attachProbeSequence",  function(object)        standardGeneric("attachProbeSequence"));
setGeneric("removeProbeSequence",  function(object)        standardGeneric("removeProbeSequence"));
setGeneric("attachUnitNames",      function(object, ...)   standardGeneric("attachUnitNames"));
setGeneric("removeUnitNames",      function(object)        standardGeneric("removeUnitNames"));
setGeneric("unitID2transcriptID",  function(object, ...)   standardGeneric("unitID2transcriptID"));
setGeneric("unitID2probesetID",    function(object, ...)   standardGeneric("unitID2probesetID"));
setGeneric("transcriptID2unitID",  function(object, ...)   standardGeneric("transcriptID2unitID"));
setGeneric("probesetID2unitID",    function(object, ...)   standardGeneric("probesetID2unitID"));
setGeneric("symbol2unitID",        function(object, ...)   standardGeneric("symbol2unitID"));
setGeneric("unitID2symbol",        function(object, ...)   standardGeneric("unitID2symbol"));


#------------------------------------------------------------------------------#
# ProcesSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("ProcesSet",
   representation(scheme = "SchemeTreeSet",
                  data   = "data.frame",
                  params = "list"
   ),
   contains=c("TreeSet"),
   prototype(scheme = new("SchemeTreeSet"),
             data   = data.frame(matrix(nr=0,nc=0)),
             params = list()
   )
)#ProcesSet

# generic methods for class ProcesSet
setGeneric("schemeFile",  function(object)        standardGeneric("schemeFile"));
setGeneric("schemeFile<-",function(object, value) standardGeneric("schemeFile<-"));
setGeneric("schemeSet",   function(object)        standardGeneric("schemeSet"));
setGeneric("schemeSet<-", function(object, value) standardGeneric("schemeSet<-"));
setGeneric("getTreeData", function(object, ...)   standardGeneric("getTreeData"));
setGeneric("attachData",  function(object, ...)   standardGeneric("attachData"));
setGeneric("removeData",  function(object)        standardGeneric("removeData"));
setGeneric("treeData",    function(object)        standardGeneric("treeData"));
setGeneric("validData",   function(object, ...)   standardGeneric("validData"));
setGeneric("boxplot",     function(x, ...)        standardGeneric("boxplot"));
setGeneric("mboxplot",    function(x, ...)        standardGeneric("mboxplot"));
setGeneric("hist",        function(x, ...)        standardGeneric("hist"));
setGeneric("image",       function(x, ...)        standardGeneric("image"));


#------------------------------------------------------------------------------#
# DataTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("DataTreeSet",
   representation(bgtreenames = "list",
                  bgrd        = "data.frame",
                  projectinfo = "ProjectInfo"
   ),
   contains=c("ProcesSet"),
   prototype(bgtreenames = list(),
             bgrd        = data.frame(matrix(nr=0,nc=0)),
             projectinfo = new("ProjectInfo")
   )
)#DataTreeSet

# generic methods for class DataTreeSet
setGeneric("intensity",         function(object)             standardGeneric("intensity"));
setGeneric("intensity<-",       function(object, ..., value) standardGeneric("intensity<-"));
setGeneric("background",        function(object)             standardGeneric("background"));
setGeneric("background<-",      function(object, value)      standardGeneric("background<-"));
setGeneric("bgtreeNames",       function(object)             standardGeneric("bgtreeNames"));
setGeneric("attachInten",       function(object, ...)        standardGeneric("attachInten"));
setGeneric("removeInten",       function(object)             standardGeneric("removeInten"));
setGeneric("attachBgrd",        function(object, ...)        standardGeneric("attachBgrd"));
setGeneric("removeBgrd",        function(object)             standardGeneric("removeBgrd"));
setGeneric("attachDataXY",      function(object)             standardGeneric("attachDataXY"));
setGeneric("removeDataXY",      function(object)             standardGeneric("removeDataXY"));
setGeneric("validBgrd",         function(object, ...)        standardGeneric("validBgrd"));
setGeneric("addData",           function(object, ...)        standardGeneric("addData"));
setGeneric("indexUnits",        function(object, ...)        standardGeneric("indexUnits"));
setGeneric("pmindex",           function(object, ...)        standardGeneric("pmindex"));
setGeneric("mmindex",           function(object, ...)        standardGeneric("mmindex"));
setGeneric("pm",                function(object, ...)        standardGeneric("pm"));
setGeneric("mm",                function(object, ...)        standardGeneric("mm"));
setGeneric("rawCELName",        function(object, ...)        standardGeneric("rawCELName"));
setGeneric("xpsRMA",            function(object, ...)        standardGeneric("xpsRMA"));
setGeneric("xpsFIRMA",          function(object, ...)        standardGeneric("xpsFIRMA"));
setGeneric("xpsMAS4",           function(object, ...)        standardGeneric("xpsMAS4"));
setGeneric("xpsMAS5",           function(object, ...)        standardGeneric("xpsMAS5"));
setGeneric("xpsMAS5Call",       function(object, ...)        standardGeneric("xpsMAS5Call"));
setGeneric("xpsDABGCall",       function(object, ...)        standardGeneric("xpsDABGCall"));
setGeneric("xpsINICall",        function(object, ...)        standardGeneric("xpsINICall"));
setGeneric("xpsPreprocess",     function(object, ...)        standardGeneric("xpsPreprocess"));
setGeneric("xpsQualityControl", function(object, ...)        standardGeneric("xpsQualityControl"));
setGeneric("xpsBgCorrect",      function(object, ...)        standardGeneric("xpsBgCorrect"));
setGeneric("xpsNormalize",      function(object, ...)        standardGeneric("xpsNormalize"));
setGeneric("xpsSummarize",      function(object, ...)        standardGeneric("xpsSummarize"));
setGeneric("xpsQualify",        function(object, ...)        standardGeneric("xpsQualify"));
setGeneric("pmplot",            function(x, ...)             standardGeneric("pmplot"));
setGeneric("probesetplot",      function(x, ...)             standardGeneric("probesetplot"));
setGeneric("intensity2GCplot",  function(x, ...)             standardGeneric("intensity2GCplot"));


#------------------------------------------------------------------------------#
# ExprTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("ExprTreeSet",
   representation(exprtype = "character",
                  normtype = "character"
   ),
   contains=c("ProcesSet"),
   prototype(exprtype = "none",
             normtype = "none"
   )
)#ExprTreeSet

# generic methods for class ExprTreeSet
setGeneric("exprType",     function(object)             standardGeneric("exprType"));
setGeneric("exprType<-",   function(object, value)      standardGeneric("exprType<-"));
setGeneric("normType",     function(object)             standardGeneric("normType"));
setGeneric("normType<-",   function(object, value)      standardGeneric("normType<-"));
setGeneric("exprs",        function(object)             standardGeneric("exprs"));
setGeneric("exprs<-",      function(object, ..., value) standardGeneric("exprs<-"));
setGeneric("se.exprs",     function(object)             standardGeneric("se.exprs"));
setGeneric("validExpr",    function(object, ...)        standardGeneric("validExpr"));
setGeneric("validSE",      function(object, ...)        standardGeneric("validSE"));
setGeneric("attachExpr",   function(object, ...)        standardGeneric("attachExpr"));
setGeneric("removeExpr",   function(object)             standardGeneric("removeExpr"));
setGeneric("xpsPreFilter", function(object, ...)        standardGeneric("xpsPreFilter"));
setGeneric("xpsUniFilter", function(object, ...)        standardGeneric("xpsUniFilter"));
setGeneric("corplot",      function(x, ...)             standardGeneric("corplot"));
setGeneric("madplot",      function(x, ...)             standardGeneric("madplot"));
setGeneric("mvaplot",      function(x, ...)             standardGeneric("mvaplot"));
setGeneric("nuseplot",     function(x, ...)             standardGeneric("nuseplot"));
setGeneric("pcaplot",      function(x, ...)             standardGeneric("pcaplot"));
setGeneric("rleplot",      function(x, ...)             standardGeneric("rleplot"));


#------------------------------------------------------------------------------#
# CallTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("CallTreeSet",
   representation(calltype = "character",
                  detcall  = "data.frame"
   ),
   contains=c("ProcesSet"),
   prototype(calltype = "mas5",
             detcall  = data.frame(matrix(nr=0,nc=0))
   )
)#CallTreeSet

# generic methods for class CallTreeSet
setGeneric("pvalData",   function(object)             standardGeneric("pvalData"));
setGeneric("pvalData<-", function(object, ..., value) standardGeneric("pvalData<-"));
setGeneric("presCall",   function(object)             standardGeneric("presCall"));
setGeneric("presCall<-", function(object, ..., value) standardGeneric("presCall<-"));
setGeneric("validPVal",  function(object, ...)        standardGeneric("validPVal"));
setGeneric("validCall",  function(object, ...)        standardGeneric("validCall"));
setGeneric("attachPVal", function(object, ...)        standardGeneric("attachPVal"));
setGeneric("removePVal", function(object)             standardGeneric("removePVal"));
setGeneric("attachCall", function(object, ...)        standardGeneric("attachCall"));
setGeneric("removeCall", function(object)             standardGeneric("removeCall"));
setGeneric("callplot",   function(x, ...)             standardGeneric("callplot"));


#------------------------------------------------------------------------------#
# QualTreeSet: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("QualTreeSet",
   representation(qualtype = "character",
                  qualopt  = "character"
   ),
   contains=c("ProcesSet"),
   prototype(qualtype = "rlm",
             qualopt  = "raw"
   )
)#QualTreeSet

# generic methods for class QualTreeSet
setGeneric("qualType",     function(object)        standardGeneric("qualType"));
setGeneric("qualType<-",   function(object, value) standardGeneric("qualType<-"));
setGeneric("qualOption",   function(object)        standardGeneric("qualOption"));
setGeneric("qualOption<-", function(object, value) standardGeneric("qualOption<-"));
setGeneric("borders",      function(object, ...)   standardGeneric("borders"));
setGeneric("residuals",    function(object, ...)   standardGeneric("residuals"));
setGeneric("weights",      function(object, ...)   standardGeneric("weights"));
setGeneric("xpsRNAdeg",    function(object, ...)   standardGeneric("xpsRNAdeg"));
setGeneric("borderplot",   function(x, ...)        standardGeneric("borderplot"));
setGeneric("coiplot",      function(x, ...)        standardGeneric("coiplot"));
setGeneric("NUSE",         function(x, ...)        standardGeneric("NUSE"));
setGeneric("RLE",          function(x, ...)        standardGeneric("RLE"));


#------------------------------------------------------------------------------#
# Filter:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class Filter
setClass("Filter",
   representation(numfilters = "numeric"
   ),
   prototype(numfilters = 0
   )
)#Filter

# generic methods for class Filter
setGeneric("numberFilters", function(object) standardGeneric("numberFilters"));


#------------------------------------------------------------------------------#
# PreFilter:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class PreFilter
setClass("PreFilter",
   representation(mad         = "list",
                  cv          = "list",
                  variance    = "list",
                  difference  = "list",
                  ratio       = "list",
                  gap         = "list",
                  lothreshold = "list",
                  hithreshold = "list",
                  quantile    = "list",
                  prescall    = "list"
   ),
   contains=c("Filter"),
   prototype(mad         = list(),
             cv          = list(),
             variance    = list(),
             difference  = list(),
             ratio       = list(),
             gap         = list(),
             lothreshold = list(),
             hithreshold = list(),
             quantile    = list(),
             prescall    = list()
   )
)#PreFilter

# generic methods for class PreFilter
setGeneric("madFilter",        function(object)        standardGeneric("madFilter"));
setGeneric("madFilter<-",      function(object, value) standardGeneric("madFilter<-"));
setGeneric("cvFilter",         function(object)        standardGeneric("cvFilter"));
setGeneric("cvFilter<-",       function(object, value) standardGeneric("cvFilter<-"));
setGeneric("varFilter",        function(object)        standardGeneric("varFilter"));
setGeneric("varFilter<-",      function(object, value) standardGeneric("varFilter<-"));
setGeneric("diffFilter",       function(object)        standardGeneric("diffFilter"));
setGeneric("diffFilter<-",     function(object, value) standardGeneric("diffFilter<-"));
setGeneric("ratioFilter",      function(object)        standardGeneric("ratioFilter"));
setGeneric("ratioFilter<-",    function(object, value) standardGeneric("ratioFilter<-"));
setGeneric("gapFilter",        function(object)        standardGeneric("gapFilter"));
setGeneric("gapFilter<-",      function(object, value) standardGeneric("gapFilter<-"));
setGeneric("lowFilter",        function(object)        standardGeneric("lowFilter"));
setGeneric("lowFilter<-",      function(object, value) standardGeneric("lowFilter<-"));
setGeneric("highFilter",       function(object)        standardGeneric("highFilter"));
setGeneric("highFilter<-",     function(object, value) standardGeneric("highFilter<-"));
setGeneric("quantileFilter",   function(object)        standardGeneric("quantileFilter"));
setGeneric("quantileFilter<-", function(object, value) standardGeneric("quantileFilter<-"));
setGeneric("callFilter",       function(object)        standardGeneric("callFilter"));
setGeneric("callFilter<-",     function(object, value) standardGeneric("callFilter<-"));


#------------------------------------------------------------------------------#
# UniFilter:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# class UniFilter
setClass("UniFilter",
   representation(foldchange = "list",
                  prescall   = "list",
                  unifilter  = "list",
                  unitest    = "list"
   ),
   contains=c("Filter"),
   prototype(foldchange = list(),
             prescall   = list(),
             unifilter  = list(),
             unitest    = list()
   )
)#UniFilter

# generic methods for class UniFilter
setGeneric("fcFilter",        function(object)        standardGeneric("fcFilter"));
setGeneric("fcFilter<-",      function(object, value) standardGeneric("fcFilter<-"));
setGeneric("unitestFilter",   function(object)        standardGeneric("unitestFilter"));
setGeneric("unitestFilter<-", function(object, value) standardGeneric("unitestFilter<-"));
setGeneric("uniTest",         function(object)        standardGeneric("uniTest"));
setGeneric("uniTest<-",       function(object, value) standardGeneric("uniTest<-"));


#------------------------------------------------------------------------------#
# FilterTreeSet:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("FilterTreeSet",
   representation(filter  = "Filter",
                  exprset = "ExprTreeSet",
                  callset = "CallTreeSet"
   ),
   contains=c("ProcesSet"),
   prototype(filter  = new("Filter"),
             exprset = new("ExprTreeSet"),
             callset = new("CallTreeSet")
#??             callset = NULL
   )
)#FilterTreeSet

# generic methods for class FilterTreeSet
setGeneric("exprTreeset", function(object) standardGeneric("exprTreeset"));
setGeneric("callTreeset", function(object) standardGeneric("callTreeset"));


#------------------------------------------------------------------------------#
# AnalysisTreeSet:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("AnalysisTreeSet",
   representation(fltrset = "FilterTreeSet"
   ),
   contains=c("ProcesSet"),
   prototype(fltrset = new("FilterTreeSet")
   )
)#AnalysisTreeSet

# generic methods for class AnalysisTreeSet
setGeneric("filterTreeset", function(object) standardGeneric("filterTreeset"));
setGeneric("validFilter",   function(object) standardGeneric("validFilter"));
setGeneric("volcanoplot",   function(x, ...) standardGeneric("volcanoplot"));


#------------------------------------------------------------------------------#


