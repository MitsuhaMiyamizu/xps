#===========================================================================#
# xps package initialization
#===========================================================================#

#.onLoad <- function(libname, pkgname) {
#   require(methods);
#   require(utils);
#}

#.onAttach <- function(libname, pkgname) {
#   packageStartupMessage(paste("\nWelcome to", pkgname, "version", packageDescription(pkgname, fields="Version")));
#   packageStartupMessage("    an R wrapper for XPS - eXpression Profiling System");
#   packageStartupMessage("    (c) Copyright 2001-2018 by Christian Stratowa");
#   packageStartupMessage("    ");
#}
.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.13")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}

.onUnload <- function(libpath) {
   library.dynam.unload( "xps", libpath )
}
