\name{ROOT}
\alias{ROOT}
\title{ROOT An Object-Oriented Data Analysis Framework}
\description{
  ROOT system overview
}
\details{
ROOT is a modular object-oriented framework aimed at solving the data analysis challenges 
of high-energy physics. The relevant features of ROOT are as follows:

Architecture: The ROOT architecture is a layered class hierarchy with over 500 classes 
divided into different categories. Most of the classes inherit from a common base class 
TObject, which provides the default behavior and protocol for all objects. 

ROOT Files: Object input/output is handled by class TFile, which has a UNIX-like directory 
structure and provides a hierarchical sequential and direct access persistent object store. 
ROOT files store information in a machine independent format and support on-the-fly data 
compression. Furthermore, ROOT files are self-describing: for every object stored in TFile, 
a dictionary describing the corresponding class is written to the file. A dictionary generator, 
called ROOTCINT, parses the class header files and generates a dictionary.
Note: TFile can be considered to be the ROOT analogon to an R environment.

Data Trees: Any object derived from TObject can be written to a file with an associated 
key TKey. However, each key has an overhead in the directory structure in memory. To 
reduce this overhead, a novel concept, called Trees (class TTree) has been developed. 
Trees are designed to support very large numbers of complex objects in a large number 
of files. A Tree consists of branches (TBranch) with each branch described by its 
leaves (TLeaf). Trees allow direct and random access to any entry of a selected subset 
of branches. Thus, Trees extend and replace the usual data tables. The concept of Tree 
friends allows the joining of many trees as one virtual tree. However, unlike table 
joins in an RDBMS, the processing time is independent of the number of tree friends.
Note: TTree can be considered to be the ROOT analogon to an R data.frame.

CINT: CINT is an interactive C/C++ interpreter, which is aimed at processing C/C++ scripts, 
called macros. Currently, CINT covers 99\% of ANSI C and 95\% of ANSI C++. CINT offers 
a gdb-like debugger for interpreted programs and allows the automatic compilation of 
scripts using ACLiC, the automatic compiler of libraries for CINT. Although available 
as independent program, CINT is embedded in ROOT as command line interpreter and macro 
processor, as well as dictionary generator. 

User interaction: The ROOT system can be accessed from the command line, by writing macros, 
or via a graphic user interface (e.g. RootBrowser). Furthermore, it is possible to write 
libraries and applications. The ROOT GUI classes allow the development of full-featured 
standalone applications.
Note: A macro can be considered to be the ROOT analogon of an R script. The RootBrowser 
can be opened using function \code{\link{root.browser}}

Platform independence: The ROOT system is available for most platforms and operating systems, 
including Linux, MacOS X, and the major flavors of UNIX and Windows. ROOT and ROOT-derived 
applications can be compiled for any supported platform. 
}
\references{
 ROOT User Guide \url{http://root.cern.ch/root/doc/RootDoc.html}

 ROOT publications \url{http://root.cern.ch/root/Publications.html}

 Christian Stratowa (2003), Distributed Storage and Analysis of Microarray Data in the 
 Terabyte Range: An Alternative to BioConductor
 \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/Stratowa.pdf}
}
\author{The ROOT team \url{http://root.cern.ch/root/Authors.html}}
\keyword{misc}
