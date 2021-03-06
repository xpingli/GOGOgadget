% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GOGOscavenger.R
\name{GOGOscavenger}
\alias{GOGOscavenger}
\title{GOGOscavenger function looks for the genes that have not been picked by GOGOpick(), but also significantly expressed. It turns to UniprotKB database, and search for the go terms.}
\usage{
GOGOscavenger(whole, pick, taxiID, write = F, ...)
}
\description{
GOGOscavenger takes 5 arguments. whole data is from GOGOpick processing where it generates a dataset containing all the significantly determined up and down genes. Screening it against the "pick", it finds out the genes that are not missed by the GOGOpick because there is not hit in PANTHER ftp. It looks for the GO terms for those genes in Uniprot instead. TaxiID is the taxonomy ID for the species, it can be found by UniProt.ws package function availableUniprotSpecies(). The function produced a whole data set of all the up and down genes ontology and the ones GOGO pick missed.
}
\examples{
\dontrun{GOGOscavenger(whole, pick, 243230, write = T, "dra24hr_ontology.csv" )}

}
\author{
Xiaoping Li
}

