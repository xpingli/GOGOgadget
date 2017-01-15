GOGOscavenger <- function(updown, pick, taxiID, write = F, ...){


        updown_ids <- updown$refID

        pick <- pick[,-1]
        pick_ids <- pick$refID

        #pick the genes that is not in the GOGOpick processed but do exist as significantly expressed in UpDown.csv

        shareIDs <- updown_ids %in% pick_ids
        Un_Onted <- updown[!shareIDs,]

        Un_Onted_refID <- Un_Onted$refID

        IDmaps <- data("UniprotIDMaps")
        mf_go <- data("mf_go")
        bp_go <- data("bp_go")
        cc_go <- data("cc_go")

        ###-------------------------------------------------------
        Uniprot_id <- vector()
        for( i in 1:length(Un_Onted_refID)){


                index <- which(IDmaps$refID == Un_Onted_refID[i])
                go <- IDmaps[index,]$UniprotKB

                if(length(go) >= 1){
                        Uniprot_id[i] <- go
                } else {
                        Uniprot_id[i] <- "NA"
                }
        }



        ##-------------------------------------------------------

        refIDs <- vector()
        uniprot_go <- vector()
        geneids <- vector()
        de <- vector()
        proteinid <- vector()
        annos <- vector()
        goid <- vector()

        for( i in 1:length(Un_Onted_refID)){

                #convert to Uniprot id

                refIDs[i] <- Un_Onted_refID[i]

                geneids[i] <- Un_Onted[Un_Onted$refID == refIDs[i],]$geneID

                proteinid[i] <- Un_Onted[Un_Onted$refID == refIDs[i],]$proteinID

                annos[i] <- Un_Onted[Un_Onted$refID == refIDs[i],]$annotation

                de[i] <- Un_Onted[Un_Onted$refID == refIDs[i],]$DE


                #Use UniProt.ws function
                org <- UniProt.ws(taxiID)

                if(Uniprot_id[i] != "NA"){

                        uniprot_go[i] <- select(org,
                                                keys = Uniprot_id[i],
                                                columns = "GO",
                                                keytype = "UNIPROTKB")[,2]


                } else {

                        uniprot_go[i] <- "Not found"

                }

        }



        scav <- data.frame(refID = refIDs, geneID = geneids, proteinID = proteinid, Family_name = NA, Sub_family_name = NA, goIDs = NA, annotation = annos, goType = NA, goTerms = uniprot_go, Pathway = NA, DE = de, stringsAsFactors = F )

        id_ref <- vector()
        geneid0 <- vector()
        proids <- vector()
        annos0 <- vector()
        goterms <- vector()
        goids <- vector()
        gotypes <- vector()
        de0 <- vector()
        n <- 0
        for(i in 1:nrow(scav)){

                terms <- scav[i,]$goTerms

                len <- length(terms)


                for(k in len){

                        id_ref[(n + 1):(n + k)] <- scav$refID[i]
                        geneid0[(n + 1):(n + k)] <- scav$geneID[i]
                        proids[(n + 1):(n + k)] <- scav$proteinID[i]
                        annos0[(n + 1):(n + k)] <- scav$annotation[i]
                        de0[(n + 1):(n + k)] <- scav$DE[i]


                        if(!is.na(terms) & terms != "Not found"){

                                goid[( n + 1) : (n + k)] <- c(gsub(".* \\[(GO:[0-9]{7,})\\]", "\\1", terms))

                                goterms[( n + 1) : (n + k)] <- c(gsub(" \\[.*", "", terms))



                        } else {
                                goid[(n + 1): (n + k)] <- NA
                                goterms[(n + 1): (n + k)]<- NA

                        }

                        if(goid[(n + 1): (n + k)] %in% mf_go){

                                gotypes[(n + 1):(n + k)] <- "MF"
                        } else if(goid[(n + 1): (n + k)] %in% bp_go){

                                gotypes[(n + 1):(n + k)] <- "BP"
                        } else if(goid[(n + 1): (n + k)] %in% cc_go){

                                gotypes[(n + 1):(n + k)] <- "CC"
                        } else {

                                gotypes[(n + 1):(n + k)] <- NA
                        }



                        n <- n + k


                }
        }

        scavenger <- data.frame(refID = id_ref, geneID = geneid0, proteinID = proids, Family_name = NA, Sub_family_name = NA, goIDs = goid, annotation = annos0,goType = gotypes, goTerms = goterms, Pathway = rep("Uniprot", n), DE = de0, stringsAsFactors = F )

        vulture <- rbind(pick, scavenger)
        if(write == TRUE){
                write.csv(vulture, ...)
        } else {
                warning("Write = F: not save the file as in .csv.")
        }

}
















