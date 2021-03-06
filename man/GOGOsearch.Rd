% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GOGOsearch.R
\name{GOGOsearch}
\alias{GOGOsearch}
\title{GOGOsearch function search the genes of interest in a GO types}
\usage{
GOGOsearch(scav, help = T, view = c("full", "search"), goterm,
  gotype = "PC", ...)
}
\description{
GOGOsearch function helps you find the specific genes of interest. It takes 5 arguments. "scav": input a scav dataset processed by 'GOGOscavenger funciton'. "help": logical if you are not sure what terms your data set has. "view": "full" gives all the regulation pattern in each GO type (MF, BP, CC, PC); "search" will not show the full view results. "goterm": the GO terms you are interested to look into in the dataset. "gotype": the four gotypes to choose from (MF, BP, CC, PC), default is PC)
}
\examples{
\dontrun{
GOGOsearch(scav, help = T, view = "full", "transporter", gotype = "MF")

}

}
\author{
Xiaoping Li
}

