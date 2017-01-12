#' GOGOpick automatically adds goTerms and other elements to the edgeR or DESeq2 processed file.
#'
#'
#' @description  GOGOpick takes 4 arguments: "dataset" contains differentially expressed genes from edgeR and DESeq2. "database" is the Ontology database file you created use "DBconstruct" function. write" defualt is TRUE would write results in a csv file in the current directory. Use ... to specify the name of the csv file.
#'
#' @author Xiaoping Li
#'
#' @examples \dontrun{GOGOpick("DE_drad_24hr_adjst.csv", [dra]OntologyDataBase, write = T, "Up and Down regulated genes by using DESeq2.csv")}
#'
#' @import dplyr
#'
#' @export
GOGOpick <- function(dataset, database, write = TRUE, ... ){


        Sig <- read.csv(dataset, stringsAsFactors = FALSE)
        names(Sig)[1] <- "refID"

        Colnamez <- colnames(Sig)
        pattern <- "log2FoldChange"

        #DESeq2
        if(pattern %in% Colnamez){
                up <-  Sig %>%
                        filter(padj <= 0.05 & log2FoldChange >= 1) %>%
                        select(refID, log2FoldChange) %>%
                        mutate(DE = "1")

                down <-  Sig %>%
                        filter(padj <= 0.05 & log2FoldChange <= -1) %>%
                        select(refID, log2FoldChange) %>%
                        mutate(DE = "-1")


        } else {
                up <-   Sig %>%
                        filter(FDR <= 0.05 & logFC >= 1) %>%
                        select(refID, logFC) %>%
                        mutate(DE = "1")

                down <- Sig %>%
                        filter(FDR <= 0.05 & logFC <= -1) %>%
                        select(refID, logFC)
                        mutate(DE = "-1")
        }

        joint_df <- as.data.frame(rbind(up, down), stringsAsFactors = FALSE)




        if(grepl(".csv", database) == TRUE){

                DB <- read.csv(database, stringsAsFactors = F)
        } else {

                DB <- database
        }

        refids <- joint_df$refID
        share <- DB$refID %in% refids

        Ont <- DB[which(share == TRUE), ]

        if(write == FALSE){
                warning("Write = F: not to produce a csv file")
        } else {

                write.csv(Ont, ...)
        }

        Ont
}

