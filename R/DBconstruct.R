#' DBconstruct function search for PANTHER ftp database and downloads the goterms and forms a csv file if write == TRUE
#'
#' @import biomartr
#' @import KEGGREST
#'
#' @author Xiaoping Li
#'
#' @description DBconstruct takes 4 arguments: ftp, species, short, write and ...: "ftp" is the ftp file name in PANTHER. "species" is full latin name of your organisim. "short" is the short name for your organism. If write = TRUE, the function will write a csv file of the full ontology into your current directory, ... what you want this csv file to name.
#' @examples
#' \dontrun{
#' DBconstruct("PTHR11.0_deinococcus", "Deinococcus radiodurans", "dra", T, "GOdatabase for dra.csv" )}
#'
#'
#' @export
DBconstruct <- function(ftp, species, short, write = F, ...){

        # Download data base from PANTHER
        tmpfil <- tempfile()
        url <- paste("ftp://ftp.pantherdb.org/sequence_classifications/11.0/PANTHER_Sequence_Classification_files/", ftp,sep = "")

        # for deinococcus: ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR11.0_deinococcus

        #for shewanella: ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR11.0_shewanella

        download.file(url, tmpfil)

        GOdb <- read.delim(tmpfil, sep = "\t", head = F, na.strings = "",stringsAsFactors = F)

        unlink(tmpfil)
        rm(tmpfil)
        # change the column names
        colnames(GOdb) <- c("Identifier", "proteinID", "PANTHER_SF_ID", "family_name","Sub_family_name", "MF", "BP", "CC", "PC", "Pathway")



        refID <- vector()
        uniprot <- vector()
        for(i in 1:nrow(GOdb)){

                uniprot[i] <- unlist(strsplit( GOdb$Identifier[i], "[[:punct:]]"))[6]

                refID[i] <- paste(unlist(strsplit( GOdb$Identifier[i],  "[[:punct:]]"))[3], unlist(strsplit( GOdb$Identifier[i],  "[[:punct:]]"))[4], sep = "_")

        }

        GOdb$refID <- refID
        GOdb$uniprotKB <- uniprot

        GOdb <- GOdb[, c(11:12, 1:10)]
        GOdb <- GOdb[, -3] # NO NEED FOR THE IDENTIFIER
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

        ##############################################################

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

                data.frame(refID = id_ensembl_x, geneID = Genid, proteinID = Proid, annotation = Annos, SF_ID = x$PANTHER_SF_ID, family_name = x$family_name, sub_family_name = x$Sub_family_name, MF = x$MF, BP = x$BP, CC = x$CC, PC = x$PC, pathway = x$Pathway, stringsAsFactors = F)



        }



        constructed_db <- construct(GOdb)



        tidy_db <- function(x){

                RefID <- vector()
                annos <- vector()
                geneids <- vector()
                proids <- vector()
                goTerm <- vector()
                goID <- vector()
                type <- vector()
                fam <- vector()
                sub_fam <- vector()
                pathways <- vector()

                count <- 0
                # x an construct() object
                for(i in 1:nrow(x)){

                        mf0 <- unlist(strsplit(unlist(strsplit(x[i,"MF"], split= ";")), split = "#"))
                        index_mf <- grep("GO:.?", mf0)
                        mf <- mf0[-index_mf]
                        GOid_mf <- mf0[index_mf]

                        if(length(mf) == 0 & length(GOid_mf) == 0){
                                mf <- 0
                                GOid_mf <- 0
                                Type_mf <- rep("MF", length(mf))
                        } else {

                                mf <- mf0[-index_mf]
                                GOid_mf <- mf0[index_mf]
                                Type_mf <- rep("MF", length(mf))
                        }


                        bp0 <- unlist(strsplit(unlist(strsplit(x[i,"BP"], split= ";")), split = "#"))
                        index_bp <- grep("GO:.?", bp0)
                        bp <- bp0[-index_bp]
                        GOid_bp <- bp0[index_bp]

                        if(length(bp) == 0 & length(GOid_bp) == 0){

                                bp <- 0
                                GOid_bp <- 0
                                Type_bp <- rep("BP", length(bp))
                        }else{

                                bp <- bp0[-index_bp]
                                GOid_bp <- bp0[index_bp]
                                Type_bp <- rep("BP", length(bp))
                        }

                        cc0 <- unlist(strsplit(unlist(strsplit(x[i,"CC"], split= ";")), split = "#"))
                        index_cc <- grep("GO:.?", cc0)
                        cc <- cc0[-index_cc]
                        GOid_cc <- cc0[index_cc]

                        if(length(cc) == 0 & length(GOid_cc) == 0){

                                cc <- 0
                                GOid_cc <- 0
                                Type_cc <- rep("CC", length(cc))

                        } else {
                                cc <- cc0[-index_cc]
                                GOid_cc <- cc0[index_cc]
                                Type_cc <- rep("CC", length(cc))
                        }

                        pc0 <- unlist(strsplit(unlist(strsplit(x[i,"PC"], split= ";")), split = "#"))
                        index_pc <- grep("PC.?", pc0)
                        pc <- pc0[-index_pc]
                        GOid_pc <- pc0[index_pc]

                        if(length(pc) == 0 & length(GOid_pc) == 0){

                                pc <- 0
                                GOid_pc <- 0
                                Type_pc <- rep("PC", length(pc))

                        } else {

                                pc <- pc0[-index_pc]
                                GOid_pc <- pc0[index_pc]
                                Type_pc <- rep("PC", length(pc))
                        }


                        lenMax <- length(mf) + length(bp) + length(cc) + length(pc)



                        for( j in lenMax ){


                                RefID[I(count+1):I(count+j)] <- x[i,]$refID
                                geneids[I(count + 1): I(count + j)] <- x[i,]$geneID
                                proids[I(count + 1): I(count + j)] <- x[i,]$proteinID
                                annos[I(count + 1): I(count + j)] <- x[i,]$annotation
                                fam[I(count + 1): I(count + j)] <- x[i,]$family_name
                                sub_fam[I(count + 1):I(count + j)]<-x[i,]$sub_family_name
                                goTerm[I(count + 1): I(count + j)] <- c(mf, bp, cc, pc)
                                goID[I(count + 1): I(count + j)] <- c(GOid_mf, GOid_bp, GOid_cc, GOid_pc)

                                type[I(count + 1): I(count + j)] <- c(Type_mf, Type_bp, Type_cc, Type_pc)

                                pathways[I(count + 1): I(count + j)] <- x[i,]$pathway

                                count <- count + j
                        }

                }

                rbind(data.frame(refID = RefID, geneID = geneids, proteinID = proids, Family_name = fam, Sub_family_name = sub_fam, goIDs = goID, annotation = annos, goType = type, goTerms = goTerm, Pathway = pathways, stringsAsFactors = F), data.frame(refID = c("DR_t13","DR_t0253", "DR_t02",   "DR_t02",   "DR_t0253", "DR_t18",   "DR_t11" ,"DR_r04",   "DR_r03",   "DR_r08",   "DR_t0253", "DR_t01",   "DR_t11",   "DR_t20","DR_2366"), geneID = NA, proteinID = NA, Family_name = NA, Sub_family_name = NA, goIDs = NA, annotation = c("tRNA-Ala", "tRNA-Asp", "tRNA-Asp", "tRNA-Asp", "tRNA-Asp", "tRNA-Pro", "tRNA-Met", "5S ribosomal RNA", "5S ribosomal RNA", "5S ribosomal RNA", "tRNA-Asp", "tRNA-Phe", "tRNA-Met", "tRNA-Leu", "50S ribosomal protein L32"), goType = rep("PC", 15), goTerms = c(rep("tRNA", 7), rep("ribosomal RNA", 3), rep("tRNA", 4), rep("ribosomal RNA", 1)), Pathway = NA, stringsAsFactors = F ))

        }

        go <- tidy_db(constructed_db)
        if(write == TRUE){

                write.csv(go, file = paste("[GO]totalDB/", ..., sep = ""))

        } else {
                warning("Write = F: not to produce a whole PANTHER Gene Ontology .csv file")
        }
        go

}
