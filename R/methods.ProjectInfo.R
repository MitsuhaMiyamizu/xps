#==============================================================================#
# methods.ProjectInfo.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# projectInfo:
# projectInfo<-:
# authorInfo:
# authorInfo<-:
# datasetInfo:
# datasetInfo<-:
# sourceInfo:
# sourceInfo<-:
# sampleInfo:
# sampleInfo<-:
# cellineInfo:
# cellineInfo<-:
# primcellInfo:
# primcellInfo<-:
# tissueInfo:
# tissueInfo<-:
# biopsyInfo:
# biopsyInfo<-:
# arrayInfo:
# arrayInfo<-:
# hybridizInfo:
# hybridizInfo<-:
# treatmentInfo:
# treatmentInfo<-:
# show:
#==============================================================================#


#------------------------------------------------------------------------------#
# ProjectInfo initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "ProjectInfo", 
   function(.Object, ...) {
      if (debug.xps()) print("------initialize:ProjectInfo------")

      .Object <- callNextMethod(.Object, ...);
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("ProjectInfo",
   function(object) {
      if (debug.xps()) print("------setValidity:ProjectInfo------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# ProjectInfo accessors:
#------------------------------------------------------------------------------#

setMethod("projectInfo", signature(object="ProjectInfo"),
   function(object) object@project
)#projectInfo

setReplaceMethod("projectInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 5) {
         stop(paste(sQuote("project"),
              "must have <name,date,type,description,comments>"));
      }#if

      object@project <- list(name        = value[1],
                             date        = value[2],
                             type        = value[3],
                             description = value[4],
                             comments    = value[5]);
      return(object);
   }
)#projectInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("authorInfo", signature(object="ProjectInfo"),
   function(object) object@author
)#authorInfo

setReplaceMethod("authorInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 8) {
         stop(paste(sQuote("author"),
              "must have <lastname,firstname,type,company,department,email,phone,comments>"));
      }#if

      object@author <- list(lastname   = value[1],
                            firstname  = value[2],
                            type       = value[3],
                            company    = value[4],
                            department = value[5],
                            email      = value[6],
                            phone      = value[7],
                            comments   = value[8]);
      return(object);
   }
)#authorInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("datasetInfo", signature(object="ProjectInfo"),
   function(object) object@dataset
)#datasetInfo

setReplaceMethod("datasetInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 7) {
         stop(paste(sQuote("dataset"),
              "must have <name,type,sample,submitter,date,description,comments>"));
      }#if

      TYPE <- c("UD", "TS", "DR", "MC");
      if (is.na(match(value[2], TYPE))) {
         stop(paste(sQuote("type"), "must be <UD,TS,DR,MC>"));
      }#if

      object@dataset <- list(name        = value[1],
                             type        = value[2],
                             sample      = value[3],
                             submitter   = value[4],
                             date        = value[5],
                             description = value[6],
                             comments    = value[7]);
      return(object);
   }
)#datasetInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("sourceInfo", signature(object="ProjectInfo"),
   function(object) object@source
)#sourceInfo

setReplaceMethod("sourceInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 6) {
         stop(paste(sQuote("source"),
              "must have <name,type,species,subspecies,description,comments>"));
      }#if

      object@source <- list(name        = value[1],
                            type        = value[2],
                            species     = value[3],
                            subspecies  = value[4],
                            description = value[5],
                            comments    = value[6]);

      return(object);
   }
)#sourceInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("sampleInfo", signature(object="ProjectInfo"),
   function(object) object@sample
)#sampleInfo

setReplaceMethod("sampleInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 12) {
         stop(paste(sQuote("sample"),
              "must have <name,type,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments>"));
      }#if

      object@sample <- list(name        = value[1],
                            type        = value[2],
                            sex         = value[3],
                            phenotype   = value[4],
                            genotype    = value[5],
                            extraction  = value[6],
                            isxenograft = as.logical(value[7]),
                            xenostrain  = value[8],
                            xenosex     = value[9],
                            xenoage     = as.double(value[10]),
                            xenoageunit = value[11],
                            comments    = value[12]);

      return(object);
   }
)#sampleInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("cellineInfo", signature(object="ProjectInfo"),
   function(object) c(object@celline, object@sample[3:12])
)#cellineInfo

setReplaceMethod("cellineInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 15) {
         stop(paste(sQuote("celline"),
              "must have <name,type,parent,atcc,modification,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments>"));
      }#if

      object@celline <- list(name        = value[1],
                            type         = value[2],
                            parent       = value[3],
                            atcc         = value[4],
                            modification = value[5]);
      sampleInfo(object) <- c(value[1], "celline", value[6:15]);

      return(object);
   }
)#cellineInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("primcellInfo", signature(object="ProjectInfo"),
   function(object) c(object@primarycell, object@sample[3:12])
)#primcellInfo

setReplaceMethod("primcellInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 14) {
         stop(paste(sQuote("primarycell"),
              "must have <name,type,date,description,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments>"));
      }#if

      object@primarycell <- list(name        = value[1],
                                 type        = value[2],
                                 date        = value[3],
                                 description = value[4]);
      sampleInfo(object) <- c(value[1], "primarycell", value[5:14]);

      return(object);
   }
)#primcellInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("tissueInfo", signature(object="ProjectInfo"),
   function(object) c(object@tissue, object@sample[3:12])
)#tissueInfo

setReplaceMethod("tissueInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 19) {
         stop(paste(sQuote("tissue"),
              "must have <name,type,development,morphology,disease,stage,donorage,ageunit,status,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments>"));
      }#if

      object@tissue <- list(name        = value[1],
                            type        = value[2],
                            development = value[3],
                            morphology  = value[4],
                            disease     = value[5],
                            stage       = value[6],
                            donorage    = as.double(value[7]),
                            ageunit     = value[8],
                            status      = value[9]);
      sampleInfo(object) <- c(value[1], "tissue", value[10:19]);

      return(object);
   }
)#tissueInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("biopsyInfo", signature(object="ProjectInfo"),
   function(object) c(object@biopsy, object@sample[3:12])
)#biopsyInfo

setReplaceMethod("biopsyInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 18) {
         stop(paste(sQuote("biopsy"),
              "must have <name,type,morphology,disease,stage,donorage,ageunit,status,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments>"));
      }#if

      object@biopsy <- list(name       = value[1],
                            type       = value[2],
                            morphology = value[3],
                            disease    = value[4],
                            stage      = value[5],
                            donorage   = as.double(value[6]),
                            ageunit    = value[7],
                            status     = value[8]);
      sampleInfo(object) <- c(value[1], "biopsy", value[9:18]);

      return(object);
   }
)#biopsyInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("arrayInfo", signature(object="ProjectInfo"),
   function(object) object@arraytype
)#arrayInfo

setReplaceMethod("arrayInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      if (length(value) != 4) {
         stop(paste(sQuote("arraytype"),
              "must have <chipname,chiptype,description,comments>"));
      }#if

      if (!isChipType(value[2])) {
         stop(paste(sQuote("arraytype"), "must be <GeneChip,GenomeChip,ExonChip>"));
      }#if

      object@arraytype <- list(chipname    = value[1],
                               chiptype    = value[2],
                               description = value[3],
                               comments    = value[4]);
      return(object);
   }
)#arrayInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("hybridizInfo", signature(object="ProjectInfo"),
   function(object) object@hybridizations
)#hybridizInfo

setReplaceMethod("hybridizInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      ## get number of hybridizations
      len  <- length(unlist(value));
      npar <- 9;   ##number of hybridization parameters
      nhyb <- as.integer(len/npar);

      if (len != npar*nhyb) {
         stop(paste("each", sQuote("hybridization"),
              "must have <name,type,inputname,date,preparation,protocol,repname,replica,comments>"));
      }#if

      ## convert to data.frame
      ds <- data.frame(matrix(unlist(value), nrow=nhyb, byrow=TRUE));
      rownames(ds) <- ds[,1];
      colnames(ds) <- c("name","type","inputname","date","preparation","protocol","repname","replica","comments");

      object@hybridizations <- ds;
      return(object);
   }
)#hybridizInfo<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("treatmentInfo", signature(object="ProjectInfo"),
   function(object) object@treatments
)#treatmentInfo

setReplaceMethod("treatmentInfo", signature(object="ProjectInfo", value="character"),
   function(object, value) {
      ## get number of treatments
      len  <- length(unlist(value));
      npar <- 8;   ##number of treatment parameters
      ntrt <- as.integer(len/npar);

      if (len != npar*ntrt) {
         stop(paste("each", sQuote("treatment"),
              "must have <name,type,concentration,concentrationunit,time,timeunit,administration,comments>"));
      }#if

      ## convert to data.frame
      ds <- data.frame(matrix(unlist(value), nrow=ntrt, byrow=TRUE));
      rownames(ds) <- ds[,1];
      colnames(ds) <- c("name","type","concentration","concentrationunit","time","timeunit","administration","comments");

      object@treatments <- ds;
      return(object);
   }
)#treatmentInfo<-


#------------------------------------------------------------------------------#
# ProjectInfo methods:
#------------------------------------------------------------------------------#

setMethod("show", signature(object="ProjectInfo"),
   function(object) {
      if (debug.xps()) print("------show.ProjectInfo------")

      type  <- c("project", "author", "dataset", "source", "sample", "celline",
                 "primarycell", "tissue", "biopsy", "arraytype",
                 "hybridizations", "treatments");
      index <- c(length(object@project)        > 0,
                 length(object@author)         > 0,
                 length(object@dataset)        > 0,
                 length(object@source)         > 0,
                 length(object@sample)         > 0,
                 length(object@celline)        > 0,
                 length(object@primarycell)    > 0,
                 length(object@tissue)         > 0,
                 length(object@biopsy)         > 0,
                 length(object@arraytype)      > 0,
                 length(object@hybridizations) > 0,
                 length(object@treatments)     > 0);

      cat("Project information:\n");
      cat("   Submitter: ", object@submitter,  "\n");
      cat("   Laboratory:", object@laboratory, "\n");
      cat("   Contact:   ", object@contact,    "\n");

      if (any(index)) {
         cat("   Information is available on:",
             paste(type[index],collapse=", "),"\n");
      }#if
   }
)#show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

