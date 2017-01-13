#' GOGOpick automatically adds goTerms and other elements to the edgeR or DESeq2 processed file.
#'
#'
#' @description  GOGOpick takes 6 arguments: "dataset" contains differentially expressed genes from edgeR and DESeq2. "database" is the Ontology database file you created use "DBconstruct" function. write" defualt is TRUE would write results in a csv file in the current directory. Use ... to specify the name of the csv file. Species is the full name of your organism, short is the short name for it. It also returns up regulation, down regulation, and a whole csv files.
#'
#' @author Xiaoping Li
#'
#' @examples \dontrun{GOGOpick("DE_drad_24hr_adjst.csv", [dra]OntologyDataBase, "Deinococcus radiodurans", "dra",write = T, "Up and Down regulated genes by using DESeq2.csv")}
#'
#' @import dplyr
#' @import biomartr
#' @import KEGGREST
#'
#' @export
GOGOpick <- function(dataset, database, species, short, write = TRUE, ... ){


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

                colnames(up)[2] <- "logFC"



                down <-  Sig %>%
                        filter(padj <= 0.05 & log2FoldChange <= -1) %>%
                        select(refID, log2FoldChange) %>%
                        mutate(DE = "-1")

                colnames(down)[2] <- "logFC"



        } else {

                up <-   Sig %>%
                        filter(FDR <= 0.05 & logFC >= 1) %>%
                        select(refID, logFC) %>%
                        mutate(DE = "1")



                down <- Sig %>%
                        filter(FDR <= 0.05 & logFC <= -1) %>%
                        select(refID, logFC) %>%
                        mutate(DE = "-1")


        }


        #==============================================================

        #library(KEGGREST)
        #library(biomartr)

        dra.proteome <- getProteome(db = "refseq", organism = species, path = file.path("_ncbi_downloads", "proteome"))

        pro <- read_proteome(dra.proteome, format = "fasta")
        pro.anno <- strsplit(pro@ranges@NAMES, ".1 ", fixed = T)

        pro.index <- list(method = "vector")
        for(i in 1:length(pro.anno)){
                pro.index[i] <- pro.anno[[i]][1]
        }

        pro.index <- unlist(pro.index)

        pro.annotation <- list(method = "vector")
        for(i in 1:length(pro.anno)){
                pro.annotation[i] <- pro.anno[[i]][2]
        }

        pro.annotation <- unlist(pro.annotation)

        # set strings as characters
        annotation_df <- data.frame(proteinID = pro.index, annotation = pro.annotation, stringsAsFactors = F)

        ##===========================================================
        ## Construct a dataset that can bridge proteinID and refID
        kg <- keggConv(short, "ncbi-proteinid")
        kg_split <- strsplit(names(kg), ":", fixed = TRUE)


        ##protein ID
        kg_proteinID <- list(method = "vector")
        for(i in 1:length(kg)){
                kg_proteinID[i] <- kg_split[[i]][2]
        }
        kg_proteinID <- unlist(kg_proteinID)

        ##refID
        kg_refID <- list(method = "vector")
        for(i in 1:length(kg)){
                split <- strsplit(kg, ":", fixed = TRUE)
                kg_refID[i] <- split[[i]][2]
        }

        kg_refID <- unlist(kg_refID)

        # set strings as characters
        kg.index <- data.frame(refID = kg_refID, proteinID = kg_proteinID, stringsAsFactors = F)
        ##========================================================

        ## gene ID
        geneid <- keggConv(short,"ncbi-geneid")


        #x has to be a keggConv object
        get_geneID <- function(x = geneid){
                #x has to be a keggConv object
                genex <- vector()
                refx <- vector()

                for ( i in 1:length(x)){
                        genex[i] <- sapply(strsplit(names(x[i]), ":"), "[", 2)
                        refx[i] <- sapply(strsplit(x[i], ":"), "[", 2)
                }

                data.frame(refID = refx, geneID = genex, stringsAsFactors = F)
        }

        kg.geneid <- get_geneID(geneid)

###########################################################

        construct <- function(x){

                id_ensembl_x <- x$refID
                id_ensembl_gene <- kg.geneid$refID
                id_ensembl_pr <- kg.index$refID
                anno_proid <- annotation_df$proteinID


                Genid <- vector()
                Proid <- vector()
                Annos <- vector()

                for(i in 1:nrow(x)){



                        if(id_ensembl_x[i] %in% id_ensembl_gene){

                                Genid[i] <- kg.geneid[kg.geneid$refID == id_ensembl_x[i], ]$geneID

                        } else {
                                Genid[i] <- "No match"
                        }

                        if(id_ensembl_x[i] %in% id_ensembl_pr){

                                Proid[i] <- kg.index[kg.index$refID == id_ensembl_x[i], ]$proteinID

                                if(Proid[i] %in% anno_proid){

                                        Annos[i] <- annotation_df[annotation_df$proteinID == Proid[i], ]$annotation

                                } else {

                                        Annos[i] <- "No match"
                                }

                        } else {
                                Proid[i] <- "No match"
                        }


                }

                data.frame(refID = id_ensembl_x, geneID = Genid, proteinID = Proid, annotation = Annos, DE = x$DE, stringsAsFactors = F)

        }

        Up <- construct(up)
        Down <- construct(down)

        whole <- rbind(Up, Down)


        write.csv(Up, "processed/[dra]upReg.csv")
        write.csv(Down, "processed/[dra]downReg.csv")
        write.csv(whole, "processed/[dra]UpDown.csv")





        if(grepl(".csv", database)[1] == TRUE){

                DB <- read.csv(database, stringsAsFactors = F)
        } else {

                DB <- as.data.frame(database, stringAsFactors = F)
        }

        joint_df <- as.data.frame(rbind(up, down), stringAsFactors = FALSE)

        refids <- joint_df$refID
        DBref <- DB$refID

        filtered <- DB[which(DBref %in% refids == TRUE),]

        for( i in 1:nrow(filtered)){


                if(filtered$refID[i] %in% up$refID){
                        filtered$DE[i] = "1"

                } else if(filtered$refID[i] %in% down$refID){
                        filtered$DE[i] = "-1"
                } else {
                        filtered$DE[i] = "No match"
                }


        }



        if(write == FALSE){
                warning("Write = F: not to produce a csv file")
        } else {

                write.csv(filtered, ...)

        }

        filtered
}

