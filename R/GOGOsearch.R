#' GOGOsearch function search the genes of interest in a GO types
#'
#' @description GOGOsearch function helps you find the specific genes of interest. It takes 5 arguments. "scav": input a scav dataset processed by 'GOGOscavenger funciton'. "help": logical if you are not sure what terms your data set has. "view": "full" gives all the regulation pattern in each GO type (MF, BP, CC, PC); "search" will not show the full view results. "goterm": the GO terms you are interested to look into in the dataset. "gotype": the four gotypes to choose from (MF, BP, CC, PC), default is PC)
#'
#' @author Xiaoping Li
#'
#' @examples
#' \dontrun{
#' GOGOsearch(scav, help = T, view = "full", "transporter", gotype = "MF")
#'
#' }
#'
#' @export
GOGOsearch <- function(scav, help = T, view = c("full", "search"), goterm, gotype = "PC", ...){



        MF <- scav[scav$goType == "MF" & !scav$goTerms ==0 & !is.na(scav$goTerms) ,]
        BP <- scav[scav$goType == "BP" & !scav$goTerms ==0 & !is.na(scav$goTerms),]
        CC <- scav[scav$goType == "CC" & !scav$goTerms ==0 & !is.na(scav$goTerms),]
        PC <- scav[scav$goType == "PC" & !scav$goTerms ==0 & !is.na(scav$goTerms),]


        # reponse to goterm
        pattern <- grep(paste0(".*",goterm, ".*"), scav$goTerms, value = F)

        pattern2 <- grep(paste0(".*", goterm, ".*"), scav$goTerms, value = T)


        filtered <- scav[pattern, ]

        mf <- filtered[filtered$goType == "MF", ]
        bp <- filtered[filtered$goType == "BP", ]
        cc <- filtered[filtered$goType == "CC", ]
        pc <- filtered[filtered$goType == "PC", ]

        if(help == TRUE ){

                goterms_mf <- unique(MF$goTerms)
                goterms_bp <- unique(BP$goTerms)
                goterms_cc <- unique(CC$goTerms)
                goterms_pc <- unique(BP$goTerms)

                cat("MF GO terms:\n")
                print(goterms_mf)
                cat("\n")
                cat("BP GO terms:\n")
                print(goterms_bp)
                cat("\n")
                cat("CC GO terms:\n")
                print(goterms_cc)
                cat("\n")
                cat("PC GO terms:\n")
                print(goterms_pc)
                cat("\n")
                cat("\n")
        } else {

                warning("Argument help = F: not showing all the GO therms in your scav data set. To see, set help = T.")
                cat("\n")
        }

        if(view == "full"){

                tb1 <- table(MF$goTerms, regulation = ifelse(MF$DE == -1, "down", "up"))
                tb2 <- table(BP$goTerms, regulation = ifelse(BP$DE == -1, "down", "up"))
                tb3 <- table(CC$goTerms, regulation = ifelse(CC$DE == -1, "down", "up"))
                tb4 <- table(PC$goTerms, regulation = ifelse(PC$DE == -1, "down", "up"))

                cat("MF full view:\n")
                print(tb1)
                cat("\n")
                cat("BP full view:\n")
                print(tb2)
                cat("\n")
                cat("CC full view:\n")
                print(tb3)
                cat("\n")
                cat("PC full veiw:\n")
                print(tb4)
                cat("\n")
                cat("\n")

        } else if(view == "search"){
                warning("Not showing the results for all the GO types. To see, set argument view = 'full'")
        }

        #response to gotype
        if(gotype == "MF"){

                up <- mf[mf$DE == 1,]
                down <- mf[mf$DE == -1,]


                id <- vector()
                lfc <- vector()
                annos <- vector()
                for(i in 1:nrow(up)){
                        id[i] <- up$refID[i]
                        lfc[i] <- up$logFC[i]
                        annos[i] <- up$annotation[i]
                }


                up_base <- data.frame(refID = id, logFC = lfc, annotation = annos)

                dup <- duplicated(up_base)

                df_up <- up_base[!dup, ]



                id2 <- vector()
                lfc2 <- vector()
                annos2 <- vector()
                for(i in 1:nrow(down)){
                        id2[i] <- down$refID[i]
                        lfc2[i] <- down$logFC[i]
                        annos2[i] <- down$annotation[i]
                }

                down_base <- data.frame(refID = id2, logFC = lfc2, annotation = annos2)

                dup2 <- duplicated(down_base)

                df_down <- down_base[!dup2,]


        } else if(gotype == "BP"){

                up <- bp[bp$DE == 1,]
                down <- bp[bp$DE == -1,]


                id <- vector()
                lfc <- vector()
                annos <- vector()
                for(i in 1:nrow(up)){
                        id[i] <- up$refID[i]
                        lfc[i] <- up$logFC[i]
                        annos[i] <- up$annotation[i]
                }


                up_base <- data.frame(refID = id, logFC = lfc, annotation = annos)

                dup <- duplicated(up_base)

                df_up <- up_base[!dup, ]



                id2 <- vector()
                lfc2 <- vector()
                annos2 <- vector()
                for(i in 1:nrow(down)){
                        id2[i] <- down$refID[i]
                        lfc2[i] <- down$logFC[i]
                        annos2[i] <- down$annotation[i]
                }

                down_base <- data.frame(refID = id2, logFC = lfc2, annotation = annos2)

                dup2 <- duplicated(down_base)

                df_down <- down_base[!dup2,]

        } else if(gotype == "CC"){

                up <- cc[cc$DE == 1,]
                down <- cc[cc$DE == -1,]


                id <- vector()
                lfc <- vector()
                annos <- vector()
                for(i in 1:nrow(up)){
                        id[i] <- up$refID[i]
                        lfc[i] <- up$logFC[i]
                        annos[i] <- up$annotation[i]
                }


                up_base <- data.frame(refID = id, logFC = lfc, annotation = annos)

                dup <- duplicated(up_base)

                df_up <- up_base[!dup, ]



                id2 <- vector()
                lfc2 <- vector()
                annos2 <- vector()
                for(i in 1:nrow(down)){
                        id2[i] <- down$refID[i]
                        lfc2[i] <- down$logFC[i]
                        annos2[i] <- down$annotation[i]
                }

                down_base <- data.frame(refID = id2, logFC = lfc2, annotation = annos2)

                dup2 <- duplicated(down_base)

                df_down <- down_base[!dup2,]

        } else if(gotype == "PC"){

                up <- pc[pc$DE == 1,]
                down <- pc[pc$DE == -1,]


                id <- vector()
                lfc <- vector()
                annos <- vector()
                for(i in 1:nrow(up)){
                        id[i] <- up$refID[i]
                        lfc[i] <- up$logFC[i]
                        annos[i] <- up$annotation[i]
                }


                up_base <- data.frame(refID = id, logFC = lfc, annotation = annos)

                dup <- duplicated(up_base)

                df_up <- up_base[!dup, ]



                id2 <- vector()
                lfc2 <- vector()
                annos2 <- vector()
                for(i in 1:nrow(down)){
                        id2[i] <- down$refID[i]
                        lfc2[i] <- down$logFC[i]
                        annos2[i] <- down$annotation[i]
                }

                down_base <- data.frame(refID = id2, logFC = lfc2, annotation = annos2)

                dup2 <- duplicated(down_base)

                df_down <- down_base[!dup2,]
        }

        cat(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Search Results >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        cat("\n")
        cat("\n")
        print(paste("The terms found in", gotype, "that are matched:"))
        cat("\n")
        cat(unique(pattern2),sep = ",")
        cat("\n")
        cat("\n")
        cat("\n")
        print(paste("**********",goterm, "up regulated in", gotype, "classification. ************"))
        if(is.na(df_up$refID[1])){
                cat("Up reguated genes not found")
        } else {
                print(df_up)
        }
        cat("\n")
        cat("\n")
        print(paste("*********",goterm, "down regulated in", gotype, "classification. ***********"))
        if(is.na(df_down$refID[1])){
                cat("down reguated genes not found")
        } else {
                print(df_down)
        }

}
