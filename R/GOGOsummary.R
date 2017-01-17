#' GOGOsummary function returns pie charts showing the percentage of GO terms according to the input GO types.
#'
#' @description 2 arguments: scavenger dataset processed by GOGOscavenger function and choose a GO type from "MF", "BP", "CC", "PC"
#'
#' @import dplyr
#' @import ggplot2
#' @import RColorBrewer
#'
#' @author Xiaoping Li
#'
#' @examples
#' \dontrun{
#' GOGOsummary(scav, "MF")
#' }
#'
#' @export
GOGOsummary <- function(scavenger, gotype, ...){



        tbl1 <- as.data.frame(table(scavenger$goTerms, scavenger$goType))
        colnames(tbl1) <- c("goTerms", "goType", "count")

        tbl2 <- as.data.frame(table(scavenger$goTerms, regulation = ifelse(scavenger$DE == -1, "down", "up")))
        colnames(tbl2) <- c("goTerms", "regulation", "count")

        tbl1 <- tbl1[-which(tbl1$goTerms == 0),]
        tbl1 <- tbl1[-which(tbl1$count == 0),]

        tbl2 <- tbl2[-which(tbl2$goTerms == 0),]
        tbl2 <- tbl2[-which(tbl2$count == 0),]

        tbl1 <- droplevels(tbl1)
        tbl2 <- droplevels(tbl2)




        MF_1 <- tbl1[which(tbl1$goType == "MF"),]
        BP_1 <- tbl1[which(tbl1$goType == "BP"),]
        CC_1 <- tbl1[which(tbl1$goType == "CC"),]
        PC_1 <- tbl1[which(tbl1$goType == "PC"),]



        # set up for when ... = "MF" or "BP", etc
        for(i in 1: length(tbl2$goTerms)){

                if(tbl2$goTerms[i] %in% as.character(tbl1$goTerms)){
                        tbl2$goType[i] <-as.character(tbl1[which(tbl1$goTerms == tbl2$goTerms[i]),]$goType)

                } else {
                        tbl2$goType[i] <- "Not found"
                }

        }



        MF_2 <- tbl2[which(tbl2$goType == "MF"),]
        BP_2 <- tbl2[which(tbl2$goType == "BP"),]
        CC_2 <- tbl2[which(tbl2$goType == "CC"),]
        PC_2 <- tbl2[which(tbl2$goType == "PC"),]

        library(dplyr)
        library(RColorBrewer)
        library(ggplot2)
        library(gridExtra)


        if(gotype == "MF"){

                vis <- brewer.pal(9, "Set1")
                cols <- colorRampPalette(vis)(nrow(MF_1))

                MF_1 <- MF_1 %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                # refactor so labels match each bar
                MF_1$goTerms <- factor(MF_1$goTerms, levels = unique(MF_1$goTerms))




                g1 <- ggplot(MF_1, aes(x= "", y = perc, fill = goTerms)) +
                        geom_bar(width = 1, stat = "identity", color = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols) +
                        #remove black diagonal line
                        guides(fill = guide_legend(override.aes = list(col = NA))) +
                        ggtitle(paste("Percentage of GO terms in", gotype)) +
                        theme(plot.title = element_text(vjust = -1,hjust = 0.8))
                #
                #
                #
                #
                # by regulation


                vis2 <- brewer.pal(12, "Paired")
                cols2 <- colorRampPalette(vis2)(nrow(MF_2))

                UP <- MF_2 %>%
                        filter(regulation == "up") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                UP <- droplevels(UP)

                UP$goTerms <- factor(UP$goTerms, levels = unique(UP$goTerms))

                DOWN <- MF_2 %>%
                        filter(regulation == "down") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))


                DOWN <- droplevels(DOWN)

                DOWN$goTerms <- factor(DOWN$goTerms, levels = unique(DOWN$goTerms))


                #plot up

                g2 <- ggplot(UP, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 1) +
                        theme_void() +
                        scale_fill_manual(values = cols2) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms Percentage in UP genes inM", gotype))


                ## plot down
                vis3 <- brewer.pal(8, "Dark2")
                cols3 <- colorRampPalette(vis3)(nrow(MF_2))

                g3 <- ggplot(DOWN, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.2, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols3) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms percentage in DOWN genes in", gotype))





        } else if(gotype == "BP" ){
                #********
                #       *
                #       *
                #********
                vis <- brewer.pal(9, "Set1")
                cols <- colorRampPalette(vis)(nrow(BP_1))

                BP_1 <- BP_1 %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                # refactor so labels match each bar
                BP_1$goTerms <- factor(BP_1$goTerms, levels = unique(BP_1$goTerms))




                g1 <- ggplot(BP_1, aes(x= "", y = perc, fill = goTerms)) +
                        geom_bar(width = 1, stat = "identity", color = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols) +
                        #remove black diagonal line
                        guides(fill = guide_legend(override.aes = list(col = NA))) +
                        ggtitle(paste("Percentage of GO terms in", gotype)) +
                        theme(plot.title = element_text(vjust = -1,hjust = 0.8))
                #
                #
                #
                #
                # by regulation


                vis2 <- brewer.pal(12, "Paired")
                cols2 <- colorRampPalette(vis2)(nrow(BP_2))

                UP <- BP_2 %>%
                        filter(regulation == "up") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                UP <- droplevels(UP)

                UP$goTerms <- factor(UP$goTerms, levels = unique(UP$goTerms))

                DOWN <- BP_2 %>%
                        filter(regulation == "down") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))


                DOWN <- droplevels(DOWN)

                DOWN$goTerms <- factor(DOWN$goTerms, levels = unique(DOWN$goTerms))


                #plot up

                g2 <- ggplot(UP, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 1) +
                        theme_void() +
                        scale_fill_manual(values = cols2) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms Percentage in UP genes in", gotype))


                ## plot down
                vis3 <- brewer.pal(8, "Dark2")
                cols3 <- colorRampPalette(vis3)(nrow(BP_2))

                g3 <- ggplot(DOWN, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.2, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols3) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms percentage in DOWN genes in", gotype))




        } else if(gotype == "CC"){
                #********
                #       *
                #       *
                #********
                vis <- brewer.pal(9, "Set1")
                cols <- colorRampPalette(vis)(nrow(CC_1))

                CC_1 <- CC_1 %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                # refactor so labels match each bar
                CC_1$goTerms <- factor(CC_1$goTerms, levels = unique(CC_1$goTerms))




                g1 <- ggplot(CC_1, aes(x= "", y = perc, fill = goTerms)) +
                        geom_bar(width = 1, stat = "identity", color = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols) +
                        #remove black diagonal line
                        guides(fill = guide_legend(override.aes = list(col = NA))) +
                        ggtitle(paste("Percentage of GO terms in", gotype)) +
                        theme(plot.title = element_text(vjust = -1,hjust = 0.8))
                #
                #
                #
                #
                # by regulation


                vis2 <- brewer.pal(12, "Paired")
                cols2 <- colorRampPalette(vis2)(nrow(CC_2))
                vis3 <- brewer.pal(8, "Dark2")
                cols3 <- colorRampPalette(vis3)(nrow(CC_2))

                UP <- CC_2 %>%
                        filter(regulation == "up") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                UP <- droplevels(UP)

                UP$goTerms <- factor(UP$goTerms, levels = unique(UP$goTerms))

                DOWN <- CC_2 %>%
                        filter(regulation == "down") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))


                DOWN <- droplevels(DOWN)

                DOWN$goTerms <- factor(DOWN$goTerms, levels = unique(DOWN$goTerms))


                #plot up

                g2 <- ggplot(UP, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 1) +
                        theme_void() +
                        scale_fill_manual(values = cols2) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms Percentage in UP genes in", gotype))


                ## plot down


                g3 <- ggplot(DOWN, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.2, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols3) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms percentage in DOWN genes in", gotype))



        } else if(gotype == "PC"){
                #********
                #       *
                #       *
                #********
                vis <- brewer.pal(9, "Set1")
                cols <- colorRampPalette(vis)(nrow(PC_1))

                PC_1 <- PC_1 %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                # refactor so labels match each bar
                PC_1$goTerms <- factor(PC_1$goTerms, levels = unique(PC_1$goTerms))




                g1 <- ggplot(PC_1, aes(x= "", y = perc, fill = goTerms)) +
                        geom_bar(width = 1, stat = "identity", color = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols) +
                        #remove black diagonal line
                        guides(fill = guide_legend(override.aes = list(col = NA))) +
                        ggtitle(paste("Percentage of GO terms in", gotype)) +
                        theme(plot.title = element_text(vjust = -1,hjust = 0.8))
                #
                #
                #
                #
                # by regulation


                vis2 <- brewer.pal(12, "Paired")
                cols2 <- colorRampPalette(vis2)(nrow(PC_2))
                vis3 <- brewer.pal(8, "Dark2")
                cols3 <- colorRampPalette(vis3)(nrow(PC_2))

                UP <- PC_2 %>%
                        filter(regulation == "up") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))

                UP <- droplevels(UP)

                UP$goTerms <- factor(UP$goTerms, levels = unique(UP$goTerms))

                DOWN <- PC_2 %>%
                        filter(regulation == "down") %>%
                        mutate(perc = round(100*count/sum(count), 1), y_pos = cumsum(perc) - perc/2, labels = paste0(perc, "%")) %>%
                        arrange(desc(y_pos))


                DOWN <- droplevels(DOWN)

                DOWN$goTerms <- factor(DOWN$goTerms, levels = unique(DOWN$goTerms))


                #plot up

                g2 <- ggplot(UP, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.25, y = y_pos, label = labels)) +
                        coord_polar("y", start = 1) +
                        theme_void() +
                        scale_fill_manual(values = cols2) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms Percentage in UP genes in", gotype))


                ## plot down


                g3 <- ggplot(DOWN, aes(x = "", y = perc, fill = goTerms)) +
                        geom_bar(stat = "identity", col = "black") +
                        geom_text(aes(x = 1.2, y = y_pos, label = labels)) +
                        coord_polar("y", start = 0) +
                        theme_void() +
                        scale_fill_manual(values = cols3) +
                        guides(fill = guide_legend(title = "GO terms", override.aes = list(col= NA))) +
                        ggtitle(paste("GO terms percentage in DOWN genes in", gotype))




        } else {
                warning("Define GO types: MF-Molecular Function, BP-Biological Process, CC-Cellular Component, PC-Protein Class")
        }


        if(nrow(UP) > 1 & nrow(DOWN) > 1){
                graph <- list(g1, g2, g3)
        } else if(nrow(UP) == 0 & nrow(DOWN) > 1){
                graph <- list(g1, g3)
        } else if(nrow(UP) > 1 & nrow(DOWN) == 0){
                graph <- list(g1,g2)
        }

        graph



}

