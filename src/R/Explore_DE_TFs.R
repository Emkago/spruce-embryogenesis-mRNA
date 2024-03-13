#' ---
#' title: "Explore TFs and genes important for embryogenesis in ZE and SE"
#' author: "Katja Stojkovic"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' 
#' # Aim  
#' # Setup  
#' Libraries  
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(writexl))

# save par
old.par <- par(no.readonly=T)

#' #' Set a palette
pal <- brewer.pal(12,"Paired")
pal_pastel <- brewer.pal(9, "Pastel1")

#' Import data  
#' * DEGs in ZE
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_all_one-percent-FDR-05-log2fc-cutoff.rda"))
#load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_UpAndDown_one-percent-FDR-05-log2fc-cutoff.rda"))

#' * DEGs in SE
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")
ddsDEGs_SE <- res_sig_genes
rm(res_sig_genes)

#' * gene expression, ZE
vsta_ZE <- read.csv(here("analysis/salmon/vst_normalised_aware_counts_ZEall.csv"), header = TRUE, row.names = 1)
# keep only genes with expresion > 0
vsta_ZE <- vsta_ZE[rowSums(vsta_ZE) > 0, ]

#' * gene expression, SE
vsta_SE <- read.csv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/counts/vst_normalised_aware_counts_SE.csv", header = TRUE, row.names = 1)
# keep only genes with expresion > 0
vsta_SE <- vsta_SE[rowSums(vsta_SE) > 0, ]

#' * sample info, SE & ZE
sampleInfoSEZE <- read.csv(here("doc/sampleInfo_SEZE.csv"), row.names = 1)

#' * divide sample info
sampleInfoSE <- sampleInfoSEZE[sampleInfoSEZE$Experiment == "SE", ]
sampleInfoZE <- sampleInfoSEZE[sampleInfoSEZE$Experiment == "ZE", ]

#' * sample info, ZE with duplicated S stages
sampleInfoZE_duplS <- read.csv(here("doc/sample_info_duplS_publication.csv"), row.names = 1)

#' List of TFs  
#' * from Plantgenie
TF_plantgenie <- read.table(here("doc/Pabies_TFs_Plantgenie_220130.tsv"), header = TRUE, sep = "\t")
#' * from PlantTFDB
TF_tfdb <- read.table(here("doc/Pab_TF_list_PlantTFDB_220130.txt"), header = TRUE, sep = "\t")
#' * from Nico
TF_C <- read.table(here("doc/TF_gene_names_and_families_names_picea_1841_HC_MC_LC"), header = FALSE, sep = "\t")

#' List of genes reported to be important for the process of embryogenesis
knownSEgenes <- read_excel(here("doc/MA_list_known_SE_genes_Feb_2022_UE_forImportInR.xlsx"), trim_ws = TRUE)
colnames(knownSEgenes) <- c("name", "ID", "bestDiamondAT", "process")
knownSEgenes <- data.frame(knownSEgenes)


#' # Analysis  
#' 
#' ## Which list to take?  
#' 
#' Annotate genes from Plantgenie with family information from PlantTFDB. Different references? Shouldn't be the reason?
all(TF_tfdb$TF_ID == TF_tfdb$Gene_ID)
all(TF_plantgenie$ID == TF_plantgenie$transcript_id)

# Only 621 genes are found in both databases
length(intersect(TF_plantgenie$ID, TF_tfdb$Gene_ID))

TF <- merge(TF_tfdb[ , c("Gene_ID", "Family")], TF_plantgenie[ , c("ID", "description", "atg_description")], 
            by.x = "Gene_ID", by.y = "ID", all = TRUE)

#' How do first two lists differ from the one Nico has used in his analysis?
length(intersect(TF_plantgenie$ID, TF_C$V1))
length(intersect(TF_tfdb$Gene_ID, TF_C$V1))

#' Use a TF list from PlantTFDB.  
#' 

#' ## Intersect DEGs with a TF list  
#' 
#' ### Prepare data  
#' Prepare DEGs - strip ".1" from the gene name so they can be easily compared to 
genes_S_all <- lapply(genes_S_all, function(x){
  rownames(x) <- sub("\\.1", "", rownames(x))
  return(x)
})
genes_FMG_all <- lapply(genes_FMG_all, function(x){
  rownames(x) <- sub("\\.1", "", rownames(x))
  return(x)
})
genes_ZE_all <- lapply(genes_ZE_all, function(x){
  rownames(x) <- sub("\\.1", "", rownames(x))
  return(x)
})
ddsDEGs_SE <- lapply(ddsDEGs_SE, function(x){
  rownames(x) <- sub("\\.1", "", rownames(x))
  return(x)
})

#' Prepare database, so it can be easily changed for a different one if necessary
db <- TF_tfdb[ , c("Gene_ID", "Family")]
colnames(db) <- c("gene", "family")

TF_S_all <- lapply(genes_S_all, function(x) {
  gint <- intersect(rownames(x), db$gene)
  xint <- x[gint, ]
  xfint <- cbind(xint, family = db$family[match(rownames(xint), db$gene)])
  return(xfint)
})

TF_FMG_all <- lapply(genes_FMG_all, function(x) {
  gint <- intersect(rownames(x), db$gene)
  xint <- x[gint, ]
  xfint <- cbind(xint, family = db$family[match(rownames(xint), db$gene)])
  return(xfint)
})

TF_ZE_all <- lapply(genes_ZE_all, function(x) {
  gint <- intersect(rownames(x), db$gene)
  xint <- x[gint, ]
  xfint <- cbind(xint, family = db$family[match(rownames(xint), db$gene)])
  return(xfint)
})

TF_SE_all <- lapply(ddsDEGs_SE, function(x) {
  gint <- intersect(rownames(x), db$gene)
  xint <- x[gint, ]
  xfint <- cbind(xint, family = db$family[match(rownames(xint), db$gene)])
  return(xfint)
})

#' ### Numbers of TFs found DE in different tissues  
barplot(elementNROWS(TF_S_all), ylim = c(0,250), las = 2)
barplot(elementNROWS(TF_FMG_all), ylim = c(0,250), las = 2)
barplot(elementNROWS(TF_ZE_all), ylim = c(0,250), las = 2)
barplot(elementNROWS(TF_SE_all), ylim = c(0,500), las = 2)

#' ### TF families  

#' Number of different TF families at each transition in different tissues  
#' 
# define function: count number of different TF families
nr.fam <- function(DETFs){
  sapply(DETFs, function(x){
    length(unique(x[ ,"family"]))
    })
}

barplot(nr.fam(TF_S_all), ylim = c(0,45), las = 2)
barplot(nr.fam(TF_FMG_all), ylim = c(0,45), las = 2)
barplot(nr.fam(TF_ZE_all), ylim = c(0,45), las = 2)
barplot(nr.fam(TF_SE_all), ylim = c(0,50), las = 2)


#' Number of different TF families in different tissues regardless of transition  
#' 
# define function: count number of different TF families
which.fam <- function(DETFs){
  sapply(DETFs, function(x){
    unique(x[ ,"family"])
  })
}

fam_S <- unique(unlist(which.fam(TF_S_all)))
fam_FMG <- unique(unlist(which.fam(TF_FMG_all)))
fam_ZE <- unique(unlist(which.fam(TF_ZE_all)))
fam_SE <- unique(unlist(which.fam(TF_SE_all)))

sapply(list(fam_S, fam_FMG, fam_ZE, fam_SE), length)

# Comparing tissues of ZE
pdf(here("figures/Intersect_DETFs_ZE_venn.pdf"), width = 3.15, height = 3.15)
venn::venn(list(SD = fam_S, FG = fam_FMG, ZEm = fam_ZE), 
           zcolor = c("chocolate4", "sandybrown", "gold"),
           #col = c("chocolate4", "sandybrown", "gold"),
           box = FALSE,
           opacity = 0.5,
           ilcs = 0.8, 
           sncs = 0.8)
dev.off()

# Families unique to zygotic embryo:
fam_ZE_uniq <- setdiff(setdiff(fam_ZE, fam_FMG), fam_S)
fam_ZE_uniq

# Families unique to female gametophyte:
setdiff(setdiff(fam_FMG, fam_ZE), fam_S)

# Expression of TFs in the 7 families that are unique to zygotic embryo:  
# which genes belong to these families
fam_ZE_uniq_genes <- lapply(TF_ZE_all, function(x){
  rownames(x[x[ , "family"] %in% fam_ZE_uniq, ]) 
})

fam_ZE_uniq_genes <- unique(unlist(fam_ZE_uniq_genes))
# distribution of genes in the families
db[db$gene %in% fam_ZE_uniq_genes, ]

# scale the expression data
vstaSE_scaled <- t(scale(t(vsta_SE)))
vstaZE_scaled <- t(scale(t(vsta_ZE)))

# remove extension ".1" from the gene names
rownames(vstaSE_scaled) <- sub("\\.1", "", rownames(vstaSE_scaled)) 
rownames(vstaZE_scaled) <- sub("\\.1", "", rownames(vstaZE_scaled)) 

# order stages by tissue and by time
sampleInfoZE_duplS <- sampleInfoZE_duplS[order(sampleInfoZE_duplS$Tissue, sampleInfoZE_duplS$Time), ]
vstaZE_scaled <- vstaZE_scaled[ , sampleInfoZE_duplS$NGI.ID]

# sampleInfoSE is already ordered by time, as well as vstaSE_scaled
all(colnames(vstaSE_scaled) == sampleInfoSE$NGI.ID)

# Plotting SE and ZE on the same heatmap (but separately normalised and scaled)

# use pheatmap to be able to include additional annotation
sample_annot_SEZE <- data.frame(stage = c(str_extract(sampleInfoZE_duplS$Time, "[0-9]"), str_extract(sampleInfoSE$Time, "[0-9]")),
                                tissue = c(rep("SD", 9), rep("FG", 21), rep("SD", 9), rep("Zem", 24), rep("SE", 24)))
row.names(sample_annot_SEZE) <- c(sampleInfoZE_duplS$NGI.ID, sampleInfoSE$NGI.ID)

TF_annot <- data.frame(family = c(db[match(fam_ZE_uniq_genes, db$gene), "family"]))
row.names(TF_annot) <- rownames(cbind(vstaZE_scaled[fam_ZE_uniq_genes, ], vstaSE_scaled[fam_ZE_uniq_genes, ]))

stage_colours <- pal[1:length(unique(sampleInfoZE_duplS$Time))]
names(stage_colours) <- c(as.character(unique(str_extract(sampleInfoZE_duplS$Time, "[0-9]"))))

TFfam_colours <- pal_pastel[c(4, 2, 5, 6, 9 ,7 , 8)]
names(TFfam_colours) <- unique(db[match(fam_ZE_uniq_genes, db$gene), "family"])

SEZE_colours <- list(stage = stage_colours,
                     tissue = c(SD = "chocolate4", SE = "seagreen", Zem = "gold", FG = "sandybrown"),
                     family = TFfam_colours) 


pheatmap(cbind(vstaZE_scaled[fam_ZE_uniq_genes, ], vstaSE_scaled[fam_ZE_uniq_genes, ]), 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_col = sample_annot_SEZE,
         annotation_row = TF_annot,
         labels_col = c(sub("-", ".", sampleInfoZE_duplS$TimePubl), sampleInfoSE$Time),
         labels_row = fam_ZE_uniq_genes, 
         annotation_colors = SEZE_colours,
         border_color = NA,
         treeheight_row = 100,
         cellheight = 33,
         cellwidth = 18,
         gaps_col = c(9,30,39, 63),
         #fontsize_row = 15,
         #fontsize_col = 15, 
         fontsize = 15,
         filename = here("analysis/figures/DETFs_ZEmUnique_inSEZE_pheatmap.pdf"),
         width = 28,
         height = 9)


# Including SE
pdf(here("figures/Intersect_DETFs_ZESE_venn.pdf"), width = 3.15, height = 3.0)
venn::venn(list(SD = fam_S, FG = fam_FMG, Zem = fam_ZE, SE = fam_SE), 
           zcolor = c("chocolate4", "sandybrown", "gold", "seagreen"), 
           box = FALSE,
           opacity = 0.5,
           ilcs = 0.8, 
           sncs = 0.8)
dev.off()

# Families unique to somatic embryogenesis:
setdiff(fam_SE, union(union(fam_S, fam_ZE), fam_FMG))


#' Distribution of TF families in all the transitions in different tissues  
TFfamDist <- function(DETF){
  df1 <- lapply(DETF, function(x){
    data.frame(table(x[ , "family"]))
  })
  df2 <- df1 %>%
    reduce(full_join, by = "Var1")
  colnames(df2) <- c("family", names(DETF))
  df2[is.na(df2)] <- 0
  return(df2)
}

# Make a table with number and percent of DE TFs from each family at every transition
famDist_S <- TFfamDist(TF_S_all)
#famDist_S_proc <- cbind(family = famDist_S[, "family"], data.frame(proportions(data.matrix(famDist_S[,-1]), 2)*100))

# exclude the last transition from FMG, as there were no DE TFs found
famDist_FMG <- TFfamDist(TF_FMG_all[c(1:(length(TF_FMG_all)-1))])
#famDist_FMG_proc <- cbind(family = famDist_FMG[, "family"], data.frame(proportions(data.matrix(famDist_FMG[,-1]), 2)*100))

famDist_ZE <- TFfamDist(TF_ZE_all)
#famDist_ZE_proc <- cbind(family = famDist_ZE[, "family"], data.frame(proportions(data.matrix(famDist_ZE[,-1]), 2)*100))

famDist_SE <- TFfamDist(TF_SE_all)
#famDist_SE_proc <- cbind(family = famDist_SE[, "family"], data.frame(proportions(data.matrix(famDist_SE[,-1]), 2)*100))

# combine all tissues in one df
famDist_all <- list(famDist_S, famDist_FMG, famDist_ZE, famDist_SE) %>%
  reduce(full_join, "family")
famDist_all[is.na(famDist_all)] <- 0

# calculate %
famDist_all_proc <- cbind(family = famDist_all[, "family"], data.frame(round(proportions(data.matrix(famDist_all[,-1]), 2)*100, 1)))

# plot  
# create a palette that will include all TF families (union from different tissues)
# this looks like a nice set
col_distinct21 <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#ffffff', '#000000')
# expand to 53 colours
col_distinct21to53 <- colorRampPalette(col_distinct21)(53)

# number of TFs in different families  
# change names of the stages/transitions as in the manuscript
publNamesTransitions <- colnames(famDist_all[, -1])
publNamesTransitions <- gsub("B", "Z", publNamesTransitions)
publNamesTransitions <- gsub("ZE", "Zem", publNamesTransitions)
publNamesTransitions <- gsub("FMG", "FG", publNamesTransitions)
publNamesTransitions <- gsub("S", "SD", publNamesTransitions)
publNamesTransitions <- sub("res_", "S", publNamesTransitions)
publNamesTransitions <- sub("vs", "_S", publNamesTransitions)

# rearrange stages in transitions - e.g. S1-S2 instead of S2-S1, use "-" instead of "_"
transition_part <- str_extract_all(publNamesTransitions, "[SZ][0-9]\\.?[A-Za-z]*")
publNamesTransitions <- unlist(
  lapply(transition_part, function(x){
    paste0(x[2], "-", x[1])
  })
)

pdf(here("figures/NrTFs_AllTransitions_AllTissues.pdf"), width = 10, height = 7)
par(mar=c(9,3,2,2)+0.1)
barplot(as.matrix(famDist_all[,-1]), 
        col = col_distinct21to53[as.integer(factor(famDist_all$family))], 
        names.arg = publNamesTransitions, 
        ylim = c(0, 500),
        ylab = "number of TFs",
        #xlab = "transitions in different tissues",
        las = 2
        )
dev.off()

# Make a 3 plot figure for the manuscript:

pdf(here("figures/NumberDETFs_TFfamDistribution_SEandZE.pdf"), width = 8, height = 11)
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE), widths = 4, heights = c(1,3,1.5))
# create plot with only black bars to integrate it above "PercentTFs_AllTransitions...pdf"
par(mar = c(2, 7, 2, 2))
barplot(as.matrix(famDist_all[,-1]), 
        col = "black",
        names.arg = rep("", ncol(famDist_all[,-1])),
        space = c(rep(0.2, 2), 0.7, rep(0.2, 4), 0.7, rep(0.2, 5), 0.7, rep(0.2, 6)),
        #ylim = c(0, 500),
        #ylab = "number of TFs",
        #xlab = "transitions in different tissues",
        las = 2, 
        cex.axis = 1.5, 
        #cex.lab = 1.5, 
        #horiz = T
        
)
title(ylab = "number of TFs", line =  5, cex.lab = 1.5)

# share of TFs in different families
# pdf(here("figures/PercentTFs_AllTransitions_AllTissues_distinct21to53.pdf"), width = 10, height = 7)
par(mar=c(11,7,2,2)+0.1, xpd = TRUE)
barplot(as.matrix(famDist_all_proc[,-1]),
        col = col_distinct21to53[as.integer(factor(famDist_all_proc$family))], 
        names.arg = publNamesTransitions, 
        space = c(rep(0.2, 2), 0.7, rep(0.2, 4), 0.7, rep(0.2, 5), 0.7, rep(0.2, 6)),
        ylim = c(0, 100),
        #ylab = "% of TFs",
        #xlab = "transitions in different tissues",
        las = 2,
        cex.axis = 1.5,
        cex.names = 1.5
        )
title(ylab = "% of TFs", line = 5, cex.lab = 1.5)

plot.new()
#par(mar = c(1,1,1,1))
legend("top", #x = 0, y = -70, ######### the legend hides behind the plot
       #inset = -0.1,
       ncol = 5,
       title = "TF families",
       legend = famDist_all_proc$family, 
       col = col_distinct21to53[as.integer(factor(famDist_all_proc$family))], 
       fill = col_distinct21to53[as.integer(factor(famDist_all_proc$family))], 
       bty = "n", cex = 1.3, title.font = 2)

dev.off()

# to restore par
par(old.par)


# Try to make a horizontal 3 plot figure for the manuscript:

pdf(here("figures/NumberDETFs_TFfamDistribution_SEandZE_horizontal.pdf"), width = 12, height = 8)
layout(matrix(c(1,2,3), 1, 3), widths = c(4,1,2), heights = 4)

# share of TFs in different families
par(mar=c(6,12,2,2))
barplot(as.matrix(famDist_all_proc[,-1]),
        col = col_distinct21to53[as.integer(factor(famDist_all_proc$family))], 
        names.arg = publNamesTransitions, 
        space = c(rep(0.2, 2), 0.7, rep(0.2, 4), 0.7, rep(0.2, 5), 0.7, rep(0.2, 6)),
        #xlab = "% of TFs",
        las = 1,
        cex.axis = 1.5,
        cex.names = 1.5,
        #cex.lab = 1.5,
        horiz = T
)
title(xlab = "% of TFs", line = 4, cex.lab = 1.5)

# create plot with only black bars to integrate it above "PercentTFs_AllTransitions...pdf"
par(mar = c(6, 2, 2, 2))
barplot(as.matrix(famDist_all[,-1]), 
        col = "black",
        names.arg = rep("", ncol(famDist_all[,-1])),
        space = c(rep(0.2, 2), 0.7, rep(0.2, 4), 0.7, rep(0.2, 5), 0.7, rep(0.2, 6)),
        xlim = c(0, 500),
        #xlab = "number of TFs",
        cex.axis = 1.5, 
        cex.lab = 1.5, 
        horiz = T
        
)
title(xlab = "number of TFs", line =  4, cex.lab = 1.5)

plot.new()
par(mar = c(6,4,10,1))
legend("top", #x = 0, y = -70, ######### the legend hides behind the plot
       #inset = -0.1,
       ncol = 2,
       title = "TF families",
       legend = famDist_all_proc$family, 
       col = col_distinct21to53[as.integer(factor(famDist_all_proc$family))], 
       fill = col_distinct21to53[as.integer(factor(famDist_all_proc$family))], 
       bty = "n", cex = 1.5, title.font = 2)

dev.off()

# to restore par
par(old.par)


#' # Expression of genes important for the process of embryogenesis (not neccessary TFs)  
#' Are all the genes that we are interested in present (expressed) in both datasets?  
table(knownSEgenes$ID %in% rownames(vstaZE_scaled))
table(knownSEgenes$ID %in% rownames(vstaSE_scaled))

# These are not present:
setdiff(knownSEgenes$ID, rownames(vstaZE_scaled))
setdiff(knownSEgenes$ID, rownames(vstaSE_scaled))
# expression of MA_69685g0020 in ZE does not change, there is only one sample in which the expression value of this gene is different
table(vstaZE_scaled["MA_69685g0020", ])

knownSEgenes_notExp <- union(setdiff(knownSEgenes$ID, rownames(vstaZE_scaled)),
                             setdiff(knownSEgenes$ID, rownames(vstaSE_scaled)))

knownSEgenes[ knownSEgenes$ID %in% knownSEgenes_notExp, ]

# exclude them from analysis
knownSEgenes <- knownSEgenes[!(knownSEgenes$ID %in% knownSEgenes_notExp), ]

# use pheatmap to be able to include additional annotation
knownSE_annot <- data.frame(process = knownSEgenes$process)
#row.names(knownSE_annot) <- rownames(cbind(vstaZE_scaled[knownSEgenes$ID, ], vstaSE_scaled[knownSEgenes$ID, ]))
# There seem to be duplicate row.names. Check!
length(knownSEgenes$ID)
length(unique(knownSEgenes$ID))

knownSEgenes[knownSEgenes$ID %in% names(which(table(knownSEgenes$ID) > 1)), ]

# what to do with them? remove ABA? Why is it the same as VP1 & 2? delete one arginase, they are the same.
knownSEgenes <- subset(knownSEgenes, name != "Absisic acid")
knownSEgenes <- knownSEgenes[-64, ]

# check again if names are unique
length(knownSEgenes$ID) == length(unique(knownSEgenes$ID))

# assign cleaned list to annotation
knownSE_annot <- data.frame(process = knownSEgenes$process)
row.names(knownSE_annot) <- rownames(cbind(vstaZE_scaled[knownSEgenes$ID, ], vstaSE_scaled[knownSEgenes$ID, ]))

process_colours <- pal_pastel[c(3, 6, 8)]
names(process_colours) <- unique(knownSEgenes$process)

SEZE_colours <- list(stage = stage_colours,
                     tissue = c(SD = "chocolate4", SE = "seagreen", Zem = "gold", FG = "sandybrown"),
                     process = process_colours)

# genes clustered by expression pattern
pheatmap(cbind(vstaZE_scaled[knownSEgenes$ID, ], vstaSE_scaled[knownSEgenes$ID, ]), 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_col = sample_annot_SEZE,
         annotation_row = knownSE_annot,
         labels_col = c(sampleInfoZE_duplS$TimePubl, sampleInfoSE$Time),
         labels_row = knownSEgenes$ID, 
         annotation_colors = SEZE_colours,
         border_color = NA,
         treeheight_row = 100,
         cellheight = 23,
         cellwidth = 18,
         gaps_col = c(9,30,39, 63),
         #fontsize_row = 15,
         #fontsize_col = 15, 
         fontsize = 15,
         filename = here("analysis/figures/knownSEgenes_inSEZE_pheatmap.pdf"),
         width = 31,
         height = 25)

# genes unclustered
pheatmap(cbind(vstaZE_scaled[knownSEgenes$ID, ], vstaSE_scaled[knownSEgenes$ID, ]), 
         scale = "none", 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         #clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_col = sample_annot_SEZE,
         annotation_row = knownSE_annot,
         labels_col = c(sampleInfoZE_duplS$TimePubl, sampleInfoSE$Time),
         labels_row = knownSEgenes$ID, 
         annotation_colors = SEZE_colours,
         border_color = NA,
         treeheight_row = 100,
         cellheight = 23,
         cellwidth = 18,
         gaps_col = c(9,30,39, 63),
         gaps_row = c(23, 33),
         #fontsize_row = 15,
         #fontsize_col = 15, 
         fontsize = 15,
         filename = here("analysis/figures/knownSEgenes_unclustered_inSEZE_pheatmap.pdf"),
         width = 31,
         height = 25)

#' Check if any of these genes are DE in our data  
knownSE_S_all <- sapply(genes_S_all, function(x) {
  knownSEgenes$ID %in% rownames(x)
})
knownSE_FMG_all <- sapply(genes_FMG_all, function(x) {
  knownSEgenes$ID %in% rownames(x)
})
knownSE_ZE_all <- sapply(genes_ZE_all, function(x) {
  knownSEgenes$ID %in% rownames(x)
})

knownSE_SE_all <- sapply(ddsDEGs_SE, function(x) {
  knownSEgenes$ID %in% rownames(x)
})

# add columns with information if a gene is DE in any of the SE or ZE stages of the same tissue or experiment
knownSE_anyDE <- data.frame(sapply(list(knownSE_S_all, knownSE_FMG_all, knownSE_ZE_all), function(x){
  apply(data.frame(x), 1, any)  
}))
colnames(knownSE_anyDE) <- c("anySD", "anyFG", "anyZem")

knownSE_DE_ZE <- data.frame(cbind(knownSE_S_all, knownSE_FMG_all, knownSE_ZE_all))
knownSE_anyDE$anyZygEmb <- apply(knownSE_DE_ZE, 1, any)

knownSE_DE_SE <- data.frame(knownSE_SE_all)
knownSE_anyDE$anySomEmb <- apply(knownSE_DE_SE, 1, any)

knownSE_anyDE <- knownSE_anyDE[ , c(5,4,1:3)]

# join DE information and annotation together
knownSEgenes_DE <- cbind(knownSEgenes, knownSE_anyDE, knownSE_DE_ZE, knownSE_DE_SE)

# write out
write.csv(knownSEgenes_DE, here("doc/knownSE_DE_SE.csv"), quote = FALSE, row.names = FALSE)
write_xlsx(knownSEgenes_DE, here("doc/knownSE_DE_SE.xlsx"))

