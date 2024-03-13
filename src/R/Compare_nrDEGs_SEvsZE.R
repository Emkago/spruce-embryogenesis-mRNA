#' ---
#' title: "Common DEGs in ZE, FMG and SE samples"
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
#' In somatic embryogenesis (SE) we have a possibility to control the processes in the lab 
#' and we know what is being induced at different stages, unlike in zygotic embryogenesis 
#' where we don't know when certain biological processes happen.  
#' One of the steps that we know is not optimal in the protocol of SE is the desiccation. 
#' Therefore we would like to compare it to the desiccation in ZE, but we first have to 
#' find out when does it actually happen and if it is quick or gradual process. In order to
#' compare the transcriptional changes in SE and ZE we will calculate number of shared genes
#' between the transitions in both types of embryogenesis. Later we will check how are DEGs 
#' important for desiccation in SE expressed in ZE.  
#' 
#' # Setup  
#' 
#' Load libraries
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(pvclust))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(formattable))
suppressPackageStartupMessages(library(SuperExactTest))
suppressPackageStartupMessages(library(stringr))

#' Get helper files
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R"))
source(here("Rtoolbox/src/plotEnrichedTreemap.R"))
source(here("Rtoolbox/src/infomapTools.R"))

#' Function to plot treemaps from enrichment saved in the list:
plotTreemapList <- function(treemapList, pathFromHere){
  lapply(names(treemapList), function(x){
    pdf(file = paste0(here(pathFromHere), x, "_treemap.pdf"), width = 10, height = 6)
    if(any(treemapList[[x]]$go$namespace == "BP")) {
      plotEnrichedTreemap(treemapList[[x]], enrichment = "go", namespace = "BP", title =  paste0(x, ": GO (BP)"), clusterColor = pal[2])
    }else{
      print("There is no enrichment in GO (BP)")
    }
    
    if(any(treemapList[[x]]$go$namespace == "CC")) {
      plotEnrichedTreemap(treemapList[[x]], enrichment = "go", namespace = "CC", title =  paste0(x, ": GO (CC)"), clusterColor = pal[11])
    }else{
      print("There is no enrichment in GO (CC)")
    }
    
    if(any(treemapList[[x]]$go$namespace == "MF")) {
      plotEnrichedTreemap(treemapList[[x]], enrichment = "go", namespace = "MF", title =  paste0(x, ": GO (MF)"), clusterColor = pal[4])
    }else{
      print("There is no enrichment in GO (MF)")
    }
    
    if(!is.null(treemapList[[x]]$mapman)) {
      plotEnrichedTreemap(treemapList[[x]], enrichment = "mapman", title = paste0(x, ": MapMan"), clusterColor = pal[6])
    }else{
      print("There is no enrichment in MapMan")
    }
    
    if(!is.null(treemapList[[x]]$pfam)) {
      plotEnrichedTreemap(treemapList[[x]], enrichment = "pfam", title = paste0(x, ": Pfam"), clusterColor = pal[7])
    }else{
      print("There is no enrichment in Pfam")
    }
    dev.off() 
  })
}

#' Set a palette
pal <- brewer.pal(12,"Paired")

#' Import data  
#' * DEGs in ZE
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_all_one-percent-FDR-05-log2fc-cutoff.rda"))

#' * DEGs in SE
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")
ddsDEGs_SE <- res_sig_genes
rm(res_sig_genes)

#' * Expression of genes in SE
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes+TEs.rda")
dds_SE <- dds
rm(dds)
# exclude TEs from this dataset, keep only genes (there are 291 TEs and ~66000 genes)
dds_SE <- dds_SE[grepl("MA_", rownames(dds_SE)), ]

vsta_SE <- read.csv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/counts/vst_normalised_aware_counts_SE.csv", header = TRUE, row.names = 1)
# keep only genes with expresion > 0
vsta_SE <- vsta_SE[rowSums(vsta_SE) > 0, ]

#' * Expression of genes in ZE
load(here("analysis/salmon/ZE-allStages_duplS_dds.rda"))
dds_ZE <- dds
rm(dds)
# reorder columns (first by tissue, then by time)
reorder_byTissueTime <- colData(dds_ZE)[order(colData(dds_ZE)$Tissue, colData(dds_ZE)$Time), ]
dds_ZE <- dds_ZE[ , reorder_byTissueTime$NGI.ID]
all(colnames(dds_ZE) == reorder_byTissueTime$NGI.ID)
all(colnames(dds_ZE) == colData(dds_ZE)$NGI.ID)

vsta_ZE <- read.csv(here("analysis/salmon/vst_normalised_aware_counts_ZEall.csv"), header = TRUE, row.names = 1)
# keep only genes with expresion > 0
vsta_ZE <- vsta_ZE[rowSums(vsta_ZE) > 0, ]
# reorder columns (first by tissue, then by time)
vsta_ZE <- vsta_ZE[ , reorder_byTissueTime$NGI.ID]
colnames(vsta_ZE) == reorder_byTissueTime$NGI.ID

#' * sample info, SE & ZE
sampleInfoSEZE <- read.csv(here("doc/sampleInfo_SEZE.csv"), row.names = 1)

#' * divide sample info
sampleInfoSE <- sampleInfoSEZE[sampleInfoSEZE$Experiment == "SE", ]
sampleInfoZE <- sampleInfoSEZE[sampleInfoSEZE$Experiment == "ZE", ]


#' # Analysis  
#' ## Find shared DEGs between SE and ZE in every transition  
#' 
#' rename stages, S for SE and Z for ZE
names(genes_S_all) <- gsub("B", "Z", names(genes_S_all))
names(genes_S_all) <- gsub("S", "SD", names(genes_S_all))
names(genes_ZE_all) <- gsub("B", "Z", names(genes_ZE_all))
names(genes_ZE_all) <- gsub("ZE", "Zem", names(genes_ZE_all))
names(genes_FMG_all) <- gsub("B", "Z", names(genes_FMG_all))
names(genes_FMG_all) <- gsub("FMG", "FG", names(genes_FMG_all))
names(ddsDEGs_SE) <- gsub("res_", "S", names(ddsDEGs_SE))
names(ddsDEGs_SE) <- gsub("vs", "_S", names(ddsDEGs_SE))

#' Calculate the intersections - nr of shared genes between SE and i) S, ii) ZEm iii) FG
sameDEGs <- sapply(ddsDEGs_SE, function(se){
  S <- lapply(genes_S_all, function(ze){
    length(intersect(rownames(se), rownames(ze)))
  })
  Z <- lapply(genes_ZE_all, function(ze){
    length(intersect(rownames(se), rownames(ze)))
  })
  F <- lapply(genes_FMG_all, function(ze){
    length(intersect(rownames(se), rownames(ze)))
  })
  return(c(unlist(S), unlist(Z), unlist(F)))
})

#' Conditional formatting of the number of shared DEGs
formattable(data.frame(sameDEGs), list(area(col = S2_S1:S8_S7) ~ color_tile("transparent", "purple")))

#' Should they be represented as % of nr of DEGs in SE or ZE samples  
#' Percent of DEGs in SE  
nrDEGs_SE <- elementNROWS(ddsDEGs_SE)
all(names(nrDEGs_SE) == colnames(sameDEGs))

sameDEGs_perSE <- t(apply(sameDEGs, 1, function(x){
  round((x/nrDEGs_SE)*100, 1)
}))
# conditional formatting
formattable(data.frame(sameDEGs_perSE), list(area(col = S2_S1:S8_S7) ~ color_tile("transparent", "purple")))

#' Percent of DEGs in ZE 
nrDEGs_ZE <- c(elementNROWS(genes_S_all), elementNROWS(genes_ZE_all), elementNROWS(genes_FMG_all))
all(names(nrDEGs_ZE) == rownames(sameDEGs))

sameDEGs_perZE <- apply(sameDEGs, 2, function(x){
  round((x/nrDEGs_ZE)*100, 1)
})
# conditional formatting 
formattable(data.frame(sameDEGs_perZE), list(area(col = S2_S1:S8_S7) ~ color_tile("transparent", "purple")))

# conditional formatting - colour each column separately
formattable(data.frame(sameDEGs_perZE), lapply(1:ncol(data.frame(sameDEGs_perZE)), function(colmn) {
    area(col = colmn) ~ color_tile("transparent", "purple")
  })
)
# if coloured by row, result is almost the same as when colouring the whole area, as the span
# of the data is very similar  

#' Are these overlaps significant? (they don't occur just because of the size of the groups in comparisons)  
# create background population of genes from all the expressed genes in SE and ZE
load(file = here("analysis/gopher/backgroundGenesSE.rda"))
load(file = here("analysis/gopher/backgroundGenesZE.rda"))
NameExpGeneSEZE <- union(NameExpGeneSE, NameExpGeneZE)
n <- length(NameExpGeneSEZE)

res_sxt <- supertest(c(lapply(ddsDEGs_SE, rownames), 
                      lapply(genes_S_all, rownames), 
                      lapply(genes_FMG_all, rownames), 
                      lapply(genes_ZE_all, rownames)), 
                     n = n, degree = 2)

# subset only the intersections of interest
comb_int <- lapply(names(ddsDEGs_SE), function(trans) {
  combS <- paste(trans, "&", names(genes_S_all))
  combFG <- paste(trans, "&", names(genes_FMG_all))
  combZE <- paste(trans, "&", names(genes_ZE_all))
  return(c(unlist(combS), unlist(combFG), unlist(combZE)))
})
comb_int <- unlist(comb_int)

# subset summary information
summ_int <- subset(summary(res_sxt)$Table, Intersections %in% comb_int)

# extract names of transitions as in the summ_int
rn <- unique(str_extract(summ_int$Intersections, "Z[0-9]\\.[A-Za-z]+_Z[0-9]\\.[A-Za-z]+"))
cn <- unique(str_extract(summ_int$Intersections, "S[0-9]_S[0-9]"))

# create a df with p-values from summ_int object in the same shape as sameDEGs
pv_int <- data.frame(matrix(summ_int$P.value, 
                     byrow = FALSE,
                     dimnames = list(rn, cn),
                     nrow = length(rn),
                     ncol = length(cn)))

# reorder rows and columns as in sameDEGs
pv_int <- pv_int[match(rownames(sameDEGs), rownames(pv_int)), ]
pv_int <- pv_int[ , match(colnames(sameDEGs), colnames(pv_int))]

pv_thr_formatter <- formatter("span",
                            style = x ~ style(color = ifelse(as.numeric(x) >= -log10(0.05), "orange", "dark grey")))

formattable(round(-log10(pv_int), 1), list(area(col = S2_S1:S8_S7) ~ color_tile("transparent", "blue")))
formattable(round(-log10(pv_int), 1), list(area(col = S2_S1:S8_S7) ~ pv_thr_formatter))

#' Final table:  
#' use % common DEGs in ZE, but use area formatting from p-value  

pv_int_log10 <- round(-log10(pv_int), 1)
colnames(pv_int_log10) <- paste0(colnames(pv_int_log10), "_p")
final <- cbind(round(data.frame(sameDEGs_perZE)), pv_int_log10)

# change labels for publication: rearrange stages in transitions - S1-S2 instead of S2-S1, use "-" instead of "_"
rownames_part <- str_extract_all(rownames(final), "Z[0-9]\\.[A-Za-z]+")
rownames(final) <- unlist(
  lapply(rownames_part, function(x){
    paste0(x[2], "-", x[1])
  })
)

colnames_part <- str_extract_all(colnames(final), "S[0-9]")
colnames(final) <- unlist(
  lapply(colnames_part, function(x){
    paste0(x[2], "-", x[1])
  })
)

# use a heatmap to display %, but colour cells based on -log10(p-value)
pdf(file = here("figures/PercentZEgenesFoundInSE.pdf"), 8, 6)
heatmap.2(x = as.matrix(final[ , 8:ncol(final)]), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          col = colorRampPalette(brewer.pal(9, "Blues")[1:6])(200),
          cellnote = as.matrix(final[ , 1:7]),
          #cellnote = as.matrix(final[ , 8:ncol(final)]), # to check if colouring is correct
          notecol = "black",
          labCol = colnames(final[ , 1:7]),
          labRow = rownames(final),
          srtCol=0,  
          adjCol = c(0.5,1),
          trace = "none",
          colsep = 5,
          key = TRUE,
          keysize = 0,
          # position elements in the matrix
          lmat=rbind( c(2, 3), c(0, 1), c(0,4) ),
          # width and height of celss in the matrix
          lwid = c(0.1, 4),
          lhei=c(0.1, 2.5, 0.6),
          margins = c(2.4, 12), 
          cexRow = 1.5, cexCol = 1.5, notecex = 1.5)
dev.off()
# cannot make legend of the right size, complaining about margins when changing parameters...

# export png, paste it together later
png(file = here("figures/PercentZEgenesFoundInSE.png"), 8, 6, units = "in", res = 150)
heatmap.2(x = as.matrix(final[ , 8:ncol(final)]), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          col = colorRampPalette(brewer.pal(9, "Blues")[1:6])(200),
          cellnote = as.matrix(final[ , 1:7]),
          notecol = "black",
          labCol = colnames(final[ , 1:7]),
          labRow = rownames(final),
          srtCol=0,  
          adjCol = c(0.5,1),
          trace = "none",
          colsep = 5,
          key = TRUE,
          keysize = 0,
          # position elements in the matrix
          lmat=rbind( c(2, 3), c(0, 1), c(0,4) ),
          # width and height of celss in the matrix
          lwid = c(0.1, 4),
          lhei=c(0.1, 2.5, 0.6),
          margins = c(2.4, 12), 
          cexRow = 1.5, cexCol = 1.5, notecex = 1.5)
dev.off()

png(file = here("figures/PercentZEgenesFoundInSE_getColorKey.png"), 8, 6, units = "in", res = 300)
heatmap.2(x = as.matrix(final[ , 8:ncol(final)]), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          col = colorRampPalette(brewer.pal(9, "Blues")[1:6])(200),
          cellnote = as.matrix(final[ , 1:7]),
          #cellnote = as.matrix(final[ , 8:ncol(final)]), # to check if colouring is correct
          notecol = "black",
          labCol = colnames(final[ , 1:7]),
          labRow = rownames(final),
          srtCol=0,  
          adjCol = c(0.5,1),
          trace = "none",
          colsep = 5,
          key = TRUE,
          keysize = 0,
          # position elements in the matrix
          lmat=rbind( c(2, 3), c(1, 1), c(0,4) ),
          # width and height of celss in the matrix
          lwid = c(1, 1),
          lhei=c(1, 4, 1.5),
          margins = c(4, 10), 
          cexRow = 1.5, 
          cexCol = 1.5, 
          notecex = 1.3)
dev.off()

pdf(file = here("figures/PercentZEgenesFoundInZE_getColourKey.pdf"), 8, 6)
heatmap.2(x = as.matrix(final[ , 8:ncol(final)]), 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none",
          col = colorRampPalette(brewer.pal(9, "Blues")[1:6])(200),
          cellnote = as.matrix(final[ , 1:7]),
          notecol = "black",
          labCol = colnames(final[ , 1:7]),
          labRow = rownames(final),
          srtCol=0,  
          adjCol = c(0.5,1),
          trace = "none",
          colsep = 5,
          key = TRUE,
          keysize = 0,
          # position elements in the matrix
          lmat=rbind( c(2, 3), c(1, 1), c(0,4) ),
          # width and height of celss in the matrix
          lwid = c(1, 1),
          lhei=c(1, 4, 1.5),
          margins = c(4, 10), 
          cexRow = 1.5, 
          cexCol = 1.5, 
          notecex = 1.3)
dev.off()


#' ## DEGs in desiccation  

#' ### Explore expression of genes DE between SE stages S5 and S6  
#' Heatmap of their expression in ZE (FG/ZEm)  i) all DEGs ii) DEGs unique to this transition and not others in SE
#' Heatmap  

#' Scale the data
vstaZE_scaled <- t(scale(t(vsta_ZE)))
vstaSE_scaled <- t(scale(t(vsta_SE)))

#' How are these genes expressed in ZE?
pdf(file = here("analysis/figures/DEGs_S5S6_inZE.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[rownames(vstaZE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]]), ], 
          trace = "none",
          Colv = FALSE,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA,
          main = "DEGs during desiccation (S5-S6) in ZE", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaZE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]])), "genes"))
dev.off() 

#' How are these genes expressed in SE?
pdf(file = here("analysis/figures/DEGs_S5S6_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]]), ], 
          trace = "none", 
          Colv = FALSE,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs during desiccation (S5-S6) in SE", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]])), "genes"))
dev.off() 

#' Display SE and ZE data on the same heatmap to see the expression patterns of the same genes in two experiments; Paste scaled data together  
#' 
# 2000 most variable genes in SE (by lfc)
DEGs_SE_S5S6 <- ddsDEGs_SE[["S6_S5"]]
names_SE_S5S6_2000 <- rownames(DEGs_SE_S5S6[order(abs(DEGs_SE_S5S6$log2FoldChange), decreasing = TRUE), ])[1:2000]
# check if they are all expressed in ZE
table(names_SE_S5S6_2000 %in% rownames(vstaZE_scaled))
# keep only those that are expressed in SE and ZE
names_SE_S5S6_2000_shared <- names_SE_S5S6_2000[names_SE_S5S6_2000 %in% rownames(vstaZE_scaled)]
# reorder rows (genes) and join the scaled data
vstaSEZE_sc_S5S6_2000 <- cbind(vstaSE_sc_S5S6_2000 <- vstaSE_scaled[names_SE_S5S6_2000_shared, ],
                               vstaZE_sc_S5S6_2000 <- vstaZE_scaled[names_SE_S5S6_2000_shared, ])

# plot heatmap with data from both experiments  

pdf(file = here("analysis/figures/DEGs_S5S6_2000_inSEZE.pdf"), width = 20, height = 12)
heatmap.2(vstaSEZE_sc_S5S6_2000, 
          trace = "none", 
          Colv = FALSE,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = c(pal[as.integer(factor(colData(dds_SE)$Stages))],
                            pal[as.integer(factor(colData(dds_ZE)$Time))]),
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = c(colData(dds_SE)$Stages, 
                     colData(dds_ZE)$Replicate),
          labRow = NA,
          key = TRUE,
          sepcolor = "white",
          colsep = c(24, 54),
          sepwidth = c(0.2, 0.2),
          main = "DEGs during desiccation (S5-S6) in SE and ZE", 
          xlab = "time points", 
          ylab = paste(length(rownames(vstaSEZE_sc_S5S6_2000)), "genes"))
legend("topright", legend = as.character(1:9), col = pal[1:9], fill = pal[1:9])
dev.off() 

# pheatmap to include the annotation ####### why it is so colour uniform? 
# change the min and max of colours (RdBu)?
sample_annot_SEZE <- data.frame(stage = c(as.character(colData(dds_SE)$Stages), colData(dds_ZE)$Time),
                                tissue = c(rep("SE", nrow(colData(dds_SE))), as.character(colData(dds_ZE)$Tissue)))
row.names(sample_annot_SEZE) <- c(rownames(colData(dds_SE)), rownames(colData(dds_ZE)))

SEZE_colours <- list(stages = c(pal[as.integer(factor(colData(dds_SE)$Stages))],
                                pal[as.integer(factor(colData(dds_ZE)$Time))]),
                     tissue = c(SE = "yellow", ZE = "dark green", FMG = "brown")) 
                     # it has to contain unique values to work!

pheatmap(vstaSEZE_sc_S5S6_2000, 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_col = sample_annot_SEZE, 
         labels_col = c(colData(dds_SE)$Stages, colData(dds_ZE)$Replicate),
         labels_row = "none", 
         annotation_colors = SEZE_colours,
         cutree_rows = 4,
         cutree_cols = 3,
         filename = here("analysis/figures/DEGs_S6vsS5_2000_inSEZE.pdf"),
         width = 15,
         height = 8)
 
#' Extract the genes with the highest expression in S5 and S6
S5S6.hclust <- hclust(dist(vstaSE_scaled[rownames(vstaSE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]]), ]), method="ward.D") 
plot(as.hclust(S5S6.hclust))

cutS5S6 <- cutree(as.hclust(S5S6.hclust), k = 4)
table(cutS5S6)
cutS5S6_df <- data.frame("cluster" = cutS5S6)

sample_annot <- data.frame(stage = colData(dds_SE)$Stages)
row.names(sample_annot) <- rownames(colData(dds_SE))

SE_colours <- list(stages = pal[as.integer(factor(colData(dds_SE)$Stages))],
                   cluster = brewer.pal(4, "Spectral"))
  
pheatmap(vstaSE_scaled[rownames(vstaSE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]]), ], 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_row = cutS5S6_df, 
         annotation_col = sample_annot, 
         labels_col = colData(dds_SE)$Stages, 
         labels_row = "none", 
         annotation_colors = SE_colours,
         cutree_rows = 4,
         filename = here("analysis/figures/DEGs_S6vsS5_inSE_byCluster.pdf"),
         width = 10,
         height = 8)

# take all the S5S6 genes, plot expression in ZE, but colour them by cluster
sample_annot_ZE <- data.frame(stage = colData(dds_ZE)$Time, tissue = colData(dds_ZE)$Tissue)
row.names(sample_annot_ZE) <- rownames(colData(dds_ZE))

S5S6_genesZE <- rownames(vstaZE_scaled[rownames(vstaZE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]]), ])
cutS5S6_dfZE <- data.frame(cluster = cutS5S6_df[rownames(cutS5S6_df) %in% S5S6_genesZE, ])
rownames(cutS5S6_dfZE) <- rownames(cutS5S6_df)[rownames(cutS5S6_df) %in% S5S6_genesZE]
############################## adapt for ZE
# create colour palette for tissue type; add annotation for the whole seed 
# (these samples were renamed and duplicated at the beginning of the analysis for easier DE analysis)
tissue_col <- apply(as.data.frame(colData(dds_ZE)), 1, function(row) {
  if (row["Time"] %in% c("B1", "B2", "B3")) {
    tissue <- "S"
  } else if (!(row["Time"] %in% c("B1", "B2", "B3")) & (row["Tissue"] == "FMG")) {
    tissue <- "FMG"
  } else {
    tissue <- "ZE"
  }
  return(tissue)
})
tissue_col

ZE_colours <- list(stages = pal[as.integer(factor(colData(dds_ZE)$Time))],
                   tissue2 = c("dark green", "brown", "olive green")[as.integer(factor(tissue_col))],
                   cluster = brewer.pal(4, "Spectral"))

pheatmap(vstaZE_scaled[rownames(vstaZE_scaled) %in% rownames(ddsDEGs_SE[["S6_S5"]]), ], 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D",
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_row = cutS5S6_dfZE,
         annotation_col = sample_annot_ZE, 
         labels_col = colData(dds_ZE)$Replicate, 
         labels_row = "none",
         annotation_colors = ZE_colours,
         filename = here("analysis/figures/DEGs_S6vsS5_inZE_byCluster2.pdf"),
         width = 13,
         height = 8)

###############################
# ################ rather plot SE and ZE together ######### take genes from Cluster 2 and plot their expression in ZE
sample_annot <- data.frame(stage = colData(dds_ZE)$Time)
row.names(sample_annot) <- rownames(colData(dds_ZE))

my_colours <- list(stages = pal[as.integer(factor(colData(dds_SE)$Stages))],
                   cluster = brewer.pal(4, "Spectral"))
pdf(file = here("analysis/figures/DEGs_S5S6_inZE_onlyCluster2.pdf"), width = 15, height = 10)
pheatmap(vstaZE_scaled[rownames(vstaZE_scaled) %in% rownames(cutS5S6_df)[which(cutS5S6_df$cluster == "2")], ], 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_row = cutS5S6_df, 
         annotation_col = sample_annot, 
         labels_col = colData(dds_SE)$Stages, 
         labels_row = "none", 
         annotation_colors = my_colours,
         cutree_rows = 4,
         filename = here("analysis/figures/DEGs_S6vsS5_inSE_byCluster.pdf"),
         width = 10,
         height = 8)
heatmap.2(vstaZE_scaled[rownames(vstaZE_scaled) %in% rownames(cutS5S6_df)[which(cutS5S6_df$cluster == "2")], ], 
          trace = "none",
          Colv = FALSE,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs during desiccation (S5-S6, Cluster 2) in ZE", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaZE_scaled) %in% rownames(cutS5S6_df)[which(cutS5S6_df$cluster == "cluster 2")]), "genes"))
dev.off() 

#' ### DEGs unique to Zem and common to Zem and SE  
#' What function have i) DEGs in Z6-Z7 common to SE and ZE, ii) DEGs in Z6-Z7 present in ZE but not SE  
desicZem_notSE <- setdiff(rownames(genes_ZE_all[["Z7.ZE_Z6.ZE"]]), rownames(ddsDEGs_SE$S6_S5))
desicZem_andSE <- intersect(rownames(genes_ZE_all[["Z7.ZE_Z6.ZE"]]), rownames(ddsDEGs_SE$S6_S5))
desicSE_notZem <- setdiff(rownames(ddsDEGs_SE$S6_S5), rownames(genes_ZE_all[["Z7.ZE_Z6.ZE"]]))

length(desicZem_notSE)
length(desicZem_andSE)
length(desicSE_notZem)

#' Functional enrichment  
enrDEGs_desicZem_notSE <- gopher(sub(".1$", "", desicZem_notSE),
                                 task = list("go","mapman","pfam"), 
                                 background = sub(".1$", "", NameExpGeneSEZE), 
                                 url="pabies")
save(enrDEGs_desicZem_notSE, file = here("analysis/gopher/enrDEGs_desicZem_notSE.rda"))

enrDEGs_desicZem_andSE <- gopher(sub(".1$", "", desicZem_andSE),
                                 task = list("go","mapman","pfam"), 
                                 background = sub(".1$", "", NameExpGeneSEZE), 
                                 url="pabies")
save(enrDEGs_desicZem_andSE, file = here("analysis/gopher/enrDEGs_desicZem_andSE.rda"))

#' Plot treemaps  
list_desic <- list(enrDEGs_desicZem_notSE, enrDEGs_desicZem_andSE)
names(list_desic) <- c("desicZem_notSE", "desicZem_andSE")
plotTreemapList(list_desic, "analysis/gopher/enrDEGS_")

#' Save also in tsv format
enr2tsv(list_desic, filePrefix = here("analysis/gopher/DEGs"))

#' Heatmap  
#' DEGs only in Zem, not SE  

#' * sample info, ZE with duplicated S stages
sampleInfoZE_duplS <- read.csv(here("doc/sample_info_duplS.csv"), row.names = 1)

# order stages by tissue and by time
sampleInfoZE_duplS <- sampleInfoZE_duplS[order(sampleInfoZE_duplS$Tissue, sampleInfoZE_duplS$Time), ]
# add names of the stages as in manuscript
sampleInfoZE_duplS$TimePubl <-sub("FMG", "FG", sampleInfoZE_duplS$Replicate)
sampleInfoZE_duplS$TimePubl <-sub("ZE", "Zem", sampleInfoZE_duplS$TimePubl)
sampleInfoZE_duplS$TimePubl <-sub("B", "Z", sampleInfoZE_duplS$TimePubl)
sampleInfoZE_duplS$TimePubl <-sub("Z1-(FG)?(Zem)?", "Z1-SD", sampleInfoZE_duplS$TimePubl)
sampleInfoZE_duplS$TimePubl <-sub("Z2-(FG)?(Zem)?", "Z2-SD", sampleInfoZE_duplS$TimePubl)
sampleInfoZE_duplS$TimePubl <-sub("Z3-(FG)?(Zem)?", "Z3-SD", sampleInfoZE_duplS$TimePubl)
sampleInfoZE_duplS$TimePubl <-sub("738", "Z8", sampleInfoZE_duplS$TimePubl)
sampleInfoZE_duplS$TimePubl <-sub("739", "Z9", sampleInfoZE_duplS$TimePubl)
write.csv(sampleInfoZE_duplS, here("doc/sample_info_duplS_publication.csv"))

vstaZE_scaled <- vstaZE_scaled[ , sampleInfoZE_duplS$NGI.ID]

sample_annot_SEZE <- data.frame(stage = c(str_extract(sampleInfoZE_duplS$Time, "[0-9]"), str_extract(sampleInfoSE$Time, "[0-9]")),
                                tissue = c(rep("SD", 9), rep("FG", 21), rep("SD", 9), rep("Zem", 24), rep("SE", 24)))
row.names(sample_annot_SEZE) <- c(sampleInfoZE_duplS$NGI.ID, sampleInfoSE$NGI.ID)
stage_colours <- pal[1:length(unique(sampleInfoZE_duplS$Time))]
names(stage_colours) <- c(as.character(unique(str_extract(sampleInfoZE_duplS$Time, "[0-9]"))))

SEZE_colours <- list(stage = stage_colours,
                     tissue = c(SD = "chocolate4", SE = "seagreen", Zem = "gold", FG = "sandybrown"))

#' Check if all the genes DE during desiccation in Zem are expressed in SE  
# all genes DE in Ze
length(desicZem_notSE)
# how many of them are expressed in SE?
length(intersect(desicZem_notSE, rownames(vstaSE_scaled)))
# gene(s) DE in Zem, but not expressed in SE
setdiff(desicZem_notSE, rownames(vstaSE_scaled))
# Annotation from Plantgenie: Nucleoside diphosphate kinase, no other useful information (no gene family)
# expression of this gene in ZE
vsta_ZE[setdiff(desicZem_notSE, rownames(vstaSE_scaled)), ]
# info from DE analysis; it has low baseMean
genes_ZE_all$Z7.ZE_Z6.ZE[setdiff(desicZem_notSE, rownames(vstaSE_scaled)), ]
# exclude it from the analysis
desicZem_notSE_cur <- desicZem_notSE[desicZem_notSE != setdiff(desicZem_notSE, rownames(vstaSE_scaled))]

pheatmap(cbind(vstaZE_scaled[desicZem_notSE_cur, ], vstaSE_scaled[desicZem_notSE_cur, ]), 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = colorRampPalette(c("white", brewer.pal(9, "OrRd")))(100),
         show_rownames = FALSE,
         annotation_col = sample_annot_SEZE,
         #annotation_row = knownSE_annot,
         labels_col = c(sampleInfoZE_duplS$TimePubl, sampleInfoSE$Time),
         #labels_row = knownSEgenes$ID, 
         annotation_colors = SEZE_colours,
         border_color = NA,
         treeheight_row = 250,
         #cellheight = 23,
         #cellwidth = 18,
         cutree_rows = 6,
         gaps_col = c(9,30,39, 63),
         #fontsize_row = 15,
         #fontsize_col = 15, 
         fontsize = 15,
         filename = here("analysis/figures/DEGS_desicZem_notSE_inSEZE_pheatmap_whiteOrRd.pdf"),
         width = 25,
         height = 20)

#' DEGs in Zem and SE  
# genes clustered by expression pattern
pheatmap(cbind(vstaZE_scaled[desicZem_andSE, ], vstaSE_scaled[desicZem_andSE, ]), 
         scale = "none", 
         cluster_cols = FALSE,
         clustering_method = "ward.D", 
         color = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
         annotation_col = sample_annot_SEZE,
         #annotation_row = knownSE_annot,
         labels_col = c(sampleInfoZE_duplS$TimePubl, sampleInfoSE$Time),
         #labels_row = NULL,
         #annotation_names_row = FALSE,
         show_rownames = FALSE,
         #angle_col = 45,
         treeheight_row = 250,
         annotation_colors = SEZE_colours,
         border_color = NA,
         cutree_rows = 6,
         #cellheight = 23,
         #cellwidth = 18,
         gaps_col = c(9,30,39, 63),
         fontsize = 15,
         filename = here("analysis/figures/DEGS_desicZem_andSE_inSEZE_pheatmap.pdf"),
         width = 25,
         height = 20)