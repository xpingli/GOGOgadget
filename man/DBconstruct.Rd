% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DBconstruct.R
\name{DBconstruct}
\alias{DBconstruct}
\title{DBconstruct function search for PANTHER ftp database and downloads the goterms and forms a csv file if write == TRUE}
\usage{
DBconstruct(ftp, species, short, write = F, ...)
}
\description{
DBconstruct takes 4 arguments: ftp, species, short, write and ...: "ftp" is the ftp file name in PANTHER. "species" is full latin name of your organisim. "short" is the short name for your organism. If write = TRUE, the function will write a csv file of the full ontology into your current directory, ... what you want this csv file to name.
}
\examples{
\dontrun{
DBconstruct("PTHR11.0_deinococcus", "Deinococcus radiodurans", "dra", T, "GOdatabase for dra.csv" )}


}
\author{
Xiaoping Li
}

