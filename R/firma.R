#==============================================================================#
# firma.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# firma: RMA preprocessing
# firma.data
# firma.expr
# firma.score
#==============================================================================#


#------------------------------------------------------------------------------#
# firma: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"firma" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         background = "antigenomic",
         normalize  = TRUE,
         option     = "probeset",
         exonlevel  = "metacore",
         method     = "mdp",
         params     = list(16384, 0.0, 1.0, 10, 0.01, 1.0),
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE) 
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsFIRMA(xps.data,
                      filename   = filename,
                      filedir    = filedir,
                      tmpdir     = tmpdir,
                      background = background,
                      normalize  = normalize,
                      option     = option,
                      exonlevel  = exonlevel,
                      method     = method,
                      params     = params,
                      xps.scheme = xps.scheme,
                      add.data   = add.data,
                      verbose    = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#firma

#------------------------------------------------------------------------------#
# firma.data: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"firma.data" <-
function(xps.data,
         probeset = NULL,
         option   = "probeset",
         type     = "LEVEL_PS")
{
   ## check for presence of valid transcript option
   option <- validTranscriptOption(option);

   if (is.na(match(type, c("LEVEL_TS", "LEVEL_PS", "SCORE_PS")))) {
      stop(paste(sQuote("type"), "must be one of <LEVEL_TS,LEVEL_PS,SCORE_PS>"));
   }#if

   if (!is(xps.data, "ExprTreeSet")) {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("ExprTreeSet")));
   }#if

   ## check if data were created with firma
   data <- xps.data@data;
   if (length(grep("fir", extenPart(colnames(data)[ncol(data)]))) == 0) {
      stop(paste(sQuote("xps.data"), "was not created with method", sQuote("firma")));
   }#if
   rownames(data) <- data[,"UNIT_ID"];

   ## get only columns of "type"
   id   <- grep(type, colnames(data));
   data <- data[,c("UnitName", "TranscriptID", colnames(data)[id])];
   colnames(data) <- namePart(colnames(data));

   ## subset of data
   if (is.null(probeset)) {
      if (option == "transcript") {
         tmp  <- split(data, data$TranscriptID);
         data <- t(sapply(names(tmp),function(x){tmp[[x]][1,c(2:ncol(data))]}));
      }#if
   } else if (option == "transcript") {
      data <- data[data[,"TranscriptID"] == probeset, , drop=FALSE];
      data <- data[1, , drop=FALSE];
      rownames(data) <- data[,"TranscriptID"];
   } else {
#      tmp <- subset(data, UnitName == probeset, drop=FALSE);
      tmp <- data[data[,"UnitName"] == probeset, , drop=FALSE];

      ## probeset can be UnitName or TranscriptID
      if (nrow(tmp)) {data <- tmp;}
      else           {data <-data[data[,"TranscriptID"] == probeset, , drop=FALSE];}
   }#if

   ## remove column TranscriptID
   data <- data[,-which(colnames(data) == "TranscriptID")];

   ## remove optional column UnitName
   colid <- which(colnames(data) == "UnitName");
   if (length(colid)) {
      ## necessary since option="exon" has duplicate ids
      len <- length(which(duplicated(data[,"UnitName"])==TRUE));
      if (len) {
         warning(paste("cannot use ", sQuote("UnitName"), "as row.names since it has <",
                 len, "> non-unique values.", sep=""));
      } else {
         if (option != "transcript") rownames(data) <- data[,"UnitName"];
         data <- data[,-colid];
      }#if
   }#if

   return(data);
}#firma.data

#------------------------------------------------------------------------------#
# firma.expr: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"firma.expr" <-
function(xps.data,
         probeset = NULL,
         option   = "probeset")
{
   if (option == "transcript") type = "LEVEL_TS" else type = "LEVEL_PS";
   return(firma.data(xps.data, probeset=probeset, option=option, type=type));
}#firma.expr

#------------------------------------------------------------------------------#
# firma.score: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"firma.score" <-
function(xps.data,
         probeset = NULL,
         option   = "probeset")
{
   if (option == "transcript") {
      stop(paste(sQuote("option"), "must be <probeset,exon>"));
   }#if
   return(firma.data(xps.data, probeset=probeset, option=option, type="SCORE_PS"));
}#firma.score
