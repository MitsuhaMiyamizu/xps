#==============================================================================#
# class SchemeTreeSet: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# chipName:
# chipType:
# chipType<-:
# chipMask:
# chipMask<-:
# unitNames:
# unitNames<-:
# probeContentGC:
# probeSequence:
# probeInfo:
# nrows:
# ncols:
# attachMask:
# removeMask:
# attachProbe:
# removeProbe:
# attachProbeContentGC:
# removeProbeContentGC:
# attachProbeSequence:
# removeProbeSequence:
# attachUnitNames:
# removeUnitNames:
# unitID2transcriptID:
# unitID2probesetID:
# transcriptID2unitID:
# probesetID2unitID:
# unitID2symbol:
# symbol2unitID:
# export:
#==============================================================================#


#------------------------------------------------------------------------------#
# SchemeTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "SchemeTreeSet", 
   function(.Object,
            chipname  = character(),
            chiptype  = character(),
            probeinfo = list(),
            ...) 
   {
      if (debug.xps()) print("------initialize:SchemeTreeSet------")

      ## set default chiptype
      if (missing(chiptype) || chiptype == "") {
         chiptype <- "GeneChip";
      }#if

      .Object <- callNextMethod(.Object,
                                chipname  = chipname,
                                chiptype  = chiptype,
                                probeinfo = probeinfo,
                                ...);
      .Object@chipname  <- chipname;
      .Object@chiptype  <- chiptype;
      .Object@probeinfo <- probeinfo;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("SchemeTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:SchemeTreeSet------")

      msg <- NULL;

      ## check for correct chiptype
      if(!isChipType(object@chiptype)) {
         msg <- validMsg(msg, paste(sQuote("chiptype"), "must be",
                         sQuote("<GeneChip,GenomeChip,ExonChip>")));
      }#if

      ## check for correct settype
      if (object@settype != "scheme") {
         msg <- validMsg(msg,
                         paste(sQuote("settype"), "is not", sQuote("scheme")));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# SchemeTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("chipName", signature(object="SchemeTreeSet"),
   function(object) object@chipname
)#chipName

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipType", signature(object="SchemeTreeSet"),
   function(object) object@chiptype
)#chipType

setReplaceMethod("chipType", signature(object="SchemeTreeSet", value="character"),
   function(object, value) {
      object@chiptype <- value;
      return(object);
   }
)#chipType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipMask", signature(object="SchemeTreeSet"),
   function(object) object@mask
)#chipMask

setReplaceMethod("chipMask", signature(object="SchemeTreeSet", value="data.frame"),
   function(object, value) {
      object@mask <- value;
      return(object);
   }
)#chipMask<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipProbe", signature(object="SchemeTreeSet"),
   function(object) object@probe
)#chipProbe

setReplaceMethod("chipProbe", signature(object="SchemeTreeSet", value="data.frame"),
   function(object, value) {
      object@probe <- value;
      return(object);
   }
)#chipProbe<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("unitNames", signature(object="SchemeTreeSet"),
   function(object, as.list = FALSE) {
      if (as.list) {
         return(split(object@unitname[,2], object@unitname[,"UNIT_ID"]))
      } else {
         rownames(object@unitname) <- object@unitname[,"UNIT_ID"];
         return(object@unitname);
      }#if
   }
)#unitNames

setReplaceMethod("unitNames", signature(object="SchemeTreeSet", value="data.frame"),
   function(object, value) {
      object@unitname <- value;
      return(object);
   }
)#unitNames<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("probeContentGC", signature(object="SchemeTreeSet"),
   function(object,
            which    = "",
            unitID   = NULL,
            unittype = "transcript")
   {
      if (debug.xps()) print("------probeContentGC.SchemeTreeSet------")

      if (min(dim(object@probe)) == 0) {
         stop(paste("slot", sQuote("probe"), "is empty"));
      } else if (length(grep("ContentGC", colnames(object@probe))) == 0) {
         stop(paste("missing column", sQuote("ContentGC")));
      }#if

      if (chipType(object) == "ExonChip" && unittype == "probeset") {
         varlist <- "fMask:fProbesetID";
         column  <- "PROBESET_ID";
      } else {
         varlist <- "fMask";
         column  <- "UNIT_ID";
      }#if

      ## get scheme mask for (x,y)
      msk <- chipMask(object);
      if (min(dim(msk)) == 0) {
         cat("slot", sQuote("mask"), "is empty, importing mask from scheme.root...\n");
         msk <- export(object,
                       treetype     = "scm",
                       varlist      = varlist,
                       as.dataframe = TRUE,
                       verbose      = FALSE);
      }#if

      if (is.null(unitID) || unitID[1] == "" || unitID[1] == "*") {
         probe <- cbind(msk[,c(column,"X","Y","Mask")], object@probe[,"ContentGC",drop=FALSE]);
      } else {
         ## check for valid UNIT_IDs
         id  <- unique(msk[,column]);
         len <- length(intersect(id, unitID));
         if (length(unitID) != len) {
            stop(paste("only", len, "of", length(unitID), column, "are valid"));
         }#if

         id    <- unlist(sapply(unitID, function(x) which(msk[,column] == x)));
         probe <- cbind(msk[id,c(column,"X","Y","Mask")], object@probe[id,"ContentGC",drop=FALSE]);
      }#if

      id <- 1:nrow(probe);
      if (chipType(object) == "GeneChip") {
         if (which == "pm") {
            id <- which(probe[,"Mask"] == 1);
         } else if (which == "mm") {
            id <- which(probe[,"Mask"] == 0);
         }#if
      } else if (chipType(object) == "GenomeChip") {
         if (which[1] == "antigenomic") {
            id <- which(probe[,"Mask"] == -2);
         } else {
            level <- exonLevel(which, "GenomeChip", as.sum=FALSE);
            id    <- unlist(sapply(level, function(x) which(probe[,"Mask"] == x)));
         }#if
      } else if (chipType(object) == "ExonChip") {
         if (which[1] == "genomic") {
            id <- which(probe[,"Mask"] == -1);
         } else if (which[1] == "antigenomic") {
            id <- which(probe[,"Mask"] == -2);
         } else {
            level <- exonLevel(which, "ExonChip", as.sum=FALSE);
            id    <- unlist(sapply(level, function(x) which(probe[,"Mask"] == x)));
         }#if
      }#if

      return(probe[id,]);
   }
)#probeContentGC

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("probeSequence", signature(object="SchemeTreeSet"),
   function(object,
            unitID = NULL) 
   {
      if (debug.xps()) print("------probeSequence.SchemeTreeSet------")

      if (min(dim(object@probe)) == 0) {
         stop(paste("slot", sQuote("probe"), "is empty"));
      } else if (length(grep("ProbeSequence", colnames(object@probe))) == 0) {
         stop(paste("missing column", sQuote("ProbeSequence")));
      }#if

      if (!(is.null(unitID) || unitID[1] == "" || unitID[1] == "*")) {
         id    <- unlist(sapply(unitID, function(x) which(object@probe[,"ProbeSetID"] == x)));
         probe <- object@probe[id,c("ProbeSetID", "ProbeX", "ProbeY", "ProbeSequence")];

         return(probe[which(probe[,"ProbeSequence"] != "N"),]);
      }#if

      return(object@probe[,c("ProbeSetID", "ProbeX", "ProbeY", "ProbeSequence")]);
   }
)#probeSequence

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("probeInfo", signature(object="SchemeTreeSet"),
   function(object) object@probeinfo
)#probeInfo

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("nrows", signature(object="SchemeTreeSet"),
   function(object) object@probeinfo$nrows
)#nrows

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("ncols", signature(object="SchemeTreeSet"),
   function(object) object@probeinfo$ncols
)#ncols


#------------------------------------------------------------------------------#
# SchemeTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("attachMask", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------attachMask.SchemeTreeSet------")

      varlist <- "fMask";
      if (chipType(object) == "ExonChip") varlist <- "fMask:fProbesetID";

      chipMask(object) <- export(object,
                                 treetype     = "scm",
                                 varlist      = varlist,
                                 as.dataframe = TRUE,
                                 verbose      = FALSE);
      return(object);
   }
)#attachMask

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeMask", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeMask.SchemeTreeSet------")

      chipMask(object) <- data.frame(matrix(nrow=0,ncol=0));
      gc(); #????
      return(object);
   }
)#removeMask

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachProbe", signature(object="SchemeTreeSet"),
   function(object, varlist) {
      if (debug.xps()) print("------attachProbe.SchemeTreeSet------")

      object@probe <- export(object,
                             treetype     = "prb",
                             varlist      = varlist,
                             as.dataframe = TRUE,
                             verbose      = FALSE);
      return(object);
   }
)#attachProbe

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeProbe", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeProbe.SchemeTreeSet------")

      object@probe <- data.frame(matrix(nrow=0,ncol=0));
      gc(); #????
      return(object);
   }
)#removeProbe

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachProbeContentGC", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------attachProbeContentGC.SchemeTreeSet------")

      return(attachProbe(object, "fNumberGC"));
   }
)#attachProbeContentGC

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeProbeContentGC", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeProbeContentGC.SchemeTreeSet------")

      return(removeProbe(object));
   }
)#removeProbeContentGC

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachProbeSequence", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------attachProbeSequence.SchemeTreeSet------")

      return(attachProbe(object, "fSequence"));
   }
)#attachProbeSequence

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeProbeSequence", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeProbeSequence.SchemeTreeSet------")

      return(removeProbe(object));
   }
)#removeProbeSequence

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachUnitNames", signature(object="SchemeTreeSet"),
   function(object, treetype = "idx") {
      if (debug.xps()) print("------attachUnitNames.SchemeTreeSet------")

      unitNames(object) <- export(object,
                                  treetype     = treetype,
                                  varlist      = "fUnitName",
                                  as.dataframe = TRUE,
                                  verbose      = FALSE);

      return(object);
   }
)#attachUnitNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeUnitNames", signature(object="SchemeTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeUnitNames.SchemeTreeSet------")

      unitNames(object) <- data.frame(matrix(nrow=0,ncol=0));
      gc(); #????
      return(object);
   }
)#removeUnitNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"unitID2affyID" <-
function(object,
         unitID,
         treetype,
         unitname,
         nunits,
         as.list = FALSE)
{
   if (debug.xps()) print("------unitID2affyID.SchemeTreeSet------")

   if (nrow(object@unitname) == nunits) {
      id <- unitNames(object);
   } else {
      treename <- unlist(treeNames(object));
      cat("slot", sQuote("unitname"), "is empty, importing data from",
          sQuote(treename[grep(treetype, treename)]), "...\n");

      id <- export(object,
                   treetype     = treetype,
                   varlist      = "fUnitName",
                   as.dataframe = TRUE,
                   verbose      = FALSE);
   }#if

   ## return mapping for all unitIDs
   if (is.null(unitID)) {
      if (as.list == TRUE) {
         id <- split(id[,unitname], id[,"UNIT_ID"]);
      } else {
         rownames(id) <- id[,"UNIT_ID"];
      }#if
      return(id);
   }#if

   id  <- split(id[,unitname], id[,"UNIT_ID"]);
   id  <- lapply(unitID, function(x) eval(parse(text=paste("id$'", x, "'", sep=""))));
   len <- length(unlist(id));
   if (len != length(unitID)) {
      stop(paste("only", len, "of", length(unitID), "UNIT_IDs are valid"));
   }#if

   ## return mapping for selected unitIDs
   if (as.list == TRUE) {
      names(id) <- unitID;
      return(id);
   } else {
      return(unlist(id));
   }#if
}#unitID2affyID

setMethod("unitID2transcriptID", signature(object="SchemeTreeSet"),
   function(object, unitID = NULL, as.list = TRUE) {
      return(unitID2affyID(object, unitID, "idx", "UnitName", object@probeinfo$nunits, as.list));
   }
)#unitID2transcriptID

setMethod("unitID2probesetID", signature(object="SchemeTreeSet"),
   function(object, unitID = NULL, as.list = TRUE) {
      if (chipType(object) == "GeneChip") {
         return(unitID2affyID(object, unitID, "idx", "UnitName",   object@probeinfo$nunits,     as.list));
      } else {
         return(unitID2affyID(object, unitID, "pbs", "ProbesetID", object@probeinfo$nprobesets, as.list));
      }#if
   }
)#unitID2probesetID

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"affyID2unitID" <-
function(object,
         affyID,
         treetype,
         unitname,
         nunits,
         as.list = FALSE)
{
   if (debug.xps()) print("------affyID2unitID.SchemeTreeSet------")

   if (nrow(object@unitname) == nunits) {
      id <- unitNames(object);
   } else {
      treename <- unlist(treeNames(object));
      cat("slot", sQuote("unitname"), "is empty, importing data from",
          sQuote(treename[grep(treetype, treename)]), "...\n");

      id <- export(object,
                   treetype     = treetype,
                   varlist      = "fUnitName",
                   as.dataframe = TRUE,
                   verbose      = FALSE);
   }#if

   ## return mapping for all affyIDs
   if (is.null(affyID)) {
      if (as.list == TRUE) {
         id <- split(id[,"UNIT_ID"], id[,unitname]);
      } else {
         id <- id[,c(unitname,"UNIT_ID")];
         if (length(unique(id[,unitname])) == nrow(id)) rownames(id) <- id[,unitname];
      }#if
      return(id);
   }#if

   id  <- split(id[,unitname], id[,"UNIT_ID"]);
   id  <- lapply(affyID, function(x) names(which(id == x)));
   len <- length(unlist(id));
   if (len != length(affyID)) {
      stop(paste("only", len, "of", length(affyID), "probeset IDs are valid"));
   }#if

   ## return mapping for selected affyIDs
   if (as.list == TRUE) {
      names(id) <- affyID;
      return(id);
   } else {
      return(unlist(id));
   }#if
}#affyID2unitID

setMethod("transcriptID2unitID", signature(object="SchemeTreeSet"),
   function(object, transcriptID = NULL, as.list = TRUE) {
      return(affyID2unitID(object, transcriptID, "idx", "UnitName", object@probeinfo$nunits, as.list));
   }
)#transcriptID2unitID

setMethod("probesetID2unitID", signature(object="SchemeTreeSet"),
   function(object, probesetID = NULL, as.list = TRUE) {
      if (chipType(object) == "GeneChip") {
         return(affyID2unitID(object, probesetID, "idx", "UnitName",  object@probeinfo$nunits, as.list));
      } else {
         return(affyID2unitID(object, probesetID, "pbs", "ProbesetID", object@probeinfo$nprobesets, as.list));
      }#if
   }
)#probesetID2unitID

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("unitID2symbol", signature(object="SchemeTreeSet"),
   function(object,
            unitID,
            unittype = "transcript", 
            as.list  = TRUE) 
   {
      if (debug.xps()) print("------unitID2symbol.SchemeTreeSet------")

      if (unittype == "probeset") {
         varlist  <- "fProbesetID:fSymbol";
         treetype <- "anp";
         colname  <- "ProbesetID";
      } else {
         varlist  <- "fTranscriptID:fSymbol";
         treetype <- "ann";
         if(chipType(object) == "ExonChip") colname  <- "TranscriptClusterID"
         else                               colname  <- "ProbesetID";
      }#if

      ann <- export(object,
                    treetype     = treetype,
                    varlist      = varlist,
                    outfile      = "tmp.txt",
                    as.dataframe = TRUE,
                    verbose      = FALSE);

      if (unittype == "probeset") {
         affid <- unitID2probesetID(object, unitID=unitID, as.list=FALSE);
      } else {
         affid <- unitID2transcriptID(object, unitID=unitID, as.list=FALSE);
      }#if

      id <- split(ann[,"GeneSymbol"], ann[,colname]);
      id <- id[unlist(lapply(affid, function(x) which(names(id) == x)))];

      if (as.list) {
         id <- as.list(id);
         names(id) <- unitID;
         return(id);
      }#if
      return(unlist(id));
   }
)#unitID2symbol

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("symbol2unitID", signature(object="SchemeTreeSet"),
   function(object,
            symbol,
            unittype = "transcript", 
            as.list  = TRUE) 
   {
      if (debug.xps()) print("------symbol2unitID.SchemeTreeSet------")

      if (unittype == "probeset") {
         varlist  <- "fProbesetID:fSymbol";
         treetype <- "anp";
         colname  <- "ProbesetID";
      } else {
         varlist  <- "fTranscriptID:fSymbol";
         treetype <- "ann";
         if(chipType(object) == "ExonChip") colname  <- "TranscriptClusterID"
         else                               colname  <- "ProbesetID";
      }#if

      ann <- export(object,
                    treetype     = treetype,
                    varlist      = varlist,
                    outfile      = "tmp.txt",
                    as.dataframe = TRUE,
                    verbose      = FALSE);

      id <- split(ann[,"GeneSymbol"], ann[,colname]);
      id <- unlist(lapply(symbol, function(x) names(which(id == x))));

      if (unittype == "probeset") {
         return(probesetID2unitID(object, probesetID=id, as.list=as.list));
      } else {
         return(transcriptID2unitID(object, transcriptID=id, as.list=as.list));
      }#if
   }
)#symbol2unitID

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"exportSchemeTreeSet" <-
function(object,
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE,
         ...) 
{
   if (debug.xps()) print("------export.SchemeTreeSet------")

   callNextMethod();

   ## get scheme file
   scheme <- rootFile(object);

   ## check for presence of parameters
   if (!(is.character(treetype) && is.character(varlist))) {
      stop("arguments <treetype,varlist> must be of type character");
   }#if

   ## check for presence of valid tree type (see utils.R)
   TYPE <- SCMTYPE;
   if (chipType(object) == "ExonChip") {
      TYPE <- c(SCMTYPE, EXNTYPE);
   }#if
   type <- match(treetype, TYPE);
   if (length(type) == 0 || is.na(type)) {
      stop(paste("invalid parameter", sQuote("treetype")));
   }#if

   if (nchar(varlist) < 1) {
      varlist <- "*";
   }#if

   ## get tree name "treeset.treename.treetype"
   treename <- paste(chipName(object), chipName(object), treetype, sep=".");

   ## check for presence of outfile=/path/outname
   outname <- paste(chipName(object), treetype, sep="_");
   outfile <- validOutfile(outname, outfile);

   ## check for presence of valid separator
   validSeparator(sep);

   ## export scheme from root file
   r <- .C("ExportScheme",
           as.character(scheme),
           as.character(treename),
           as.character(varlist),
           as.character(outfile),
           as.character(sep),
           as.integer(verbose),
           PACKAGE="xps");

   ## import outfile as dataframe
   ds <- NULL;
   if (as.dataframe) {
      ## use quote = "\"", since annotation has single quotes e.g. "2'-PDE"
      ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep=sep,
                       row.names=NULL, stringsAsFactors=FALSE, quote = "\""); #"
   }#if

   return(ds);
}#exportSchemeTreeSet

setMethod("export", "SchemeTreeSet", exportSchemeTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
