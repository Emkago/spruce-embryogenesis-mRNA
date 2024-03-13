#' ---
#' title: "Spruce somatic embryogenesis: Differential Expression"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Setup  
#' ## Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(BiocParallel))

register(MulticoreParam(4))

#' ## Helper files
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/densityPlot.R"))
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/volcanoPlot.R"))

#' ## Functions
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

#' ## Data  
#' DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes+TEs.rda")

# save it without TEs
dds_genes <- dds[grepl("^MA_", rownames(dds)), ]
save(dds_genes, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes.rda")

#' Gene annotation from ConGenIE, downloaded 7th July 2020
gene_table <- read.delim("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/Pabies_gene_table_confidence-PFAM-GO-bestAt.tsv",
                         stringsAsFactors = FALSE)

#' ## Graphics
pal12 <- brewer.pal(12,"Paired")
mar <- par("mar")

#' Run commented code bellow only the first time, in the future skip and load the results  
#' 
#' #' # Differential Expression  
#' #' Biological quality check was already done by Camilla (BiologicalQA_genes+TEs.html), 
#' #' I will proceed directly to differential expression analysis
#' dds <- DESeq(dds, parallel = TRUE)
#' 
#' #' Plot estimated dispersion
#' plotDispEsts(dds)
#'   
#' #'  ## Obtain the results  
#' resultsNames(dds)
#' 
#' #' Compare expression of genes in consecutive stages  
#' #' specify alpha
#' alpha = 0.01
#' 
#' res_list <- lapply(2:8, function(i){
#'   results(dds, c("Stages", paste0("S", as.character(i)), paste0("S", as.character(i-1))), 
#'           filter = rowMedians(counts(dds)), alpha = alpha, parallel = TRUE)
#' })
#' 
#' names(res_list) <- lapply(2:8, function(i){
#'   paste0("res_", i, "vs", i-1)
#' })
#' 
#' #' Save object with unfiltered results:
#' save(res_list, 
#'      file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_unfiltered.rda")
#' # save results of the intercept
#' resInterceptSE <- results(dds, name = "Intercept", alpha = alpha, filter = rowMedians(counts(dds)))
#' # number of expressed genes:
#' nrow(resInterceptSE[!is.na(resInterceptSE$padj), ])
#' save(resInterceptSE, 
#'      file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_unfiltered_onlyIntercept.rda")
#' 
#' #' ## Filter results  
#' #' By default results are filtered by lfc and padj: p = 0.01, log2fc = 0.5
#' res_sig <- lapply(res_list, sigDeg)
#' 
#' #' Save object with significant differentially expressed genes and TEs
#' save(res_sig, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_padj001_lfc05.rda")
#' 
#' #' For the purpose of this project we only need to look at genes, not TEs. To make further 
#' #' analysis easier, exclude TEs from DE results.
#' res_sig_genes <- lapply(res_sig, function(x){
#'   x[grepl("MA_", rownames(x)), ]
#' })
#' 
#' #' Save object with only significant differentialy expressed genes
#' save(res_sig_genes, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")

# the second time the code is ran, load the results instead of computing them again
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")

#' # Number of DEGs  
#' 
#' ## Median vs Average
lapply(res_sig_genes, function(x){DESeq2::plotMA(x, 0.01)})  

#' ## Volcano plot  
#' a large number of genes show significant fold-changes

lapply(res_sig_genes, function(x){
  volcanoPlot(x,alpha=0.01,lfc = 0.5)
})

#' ## Number of DEGs in comparisons of consecutive stages  
#' Count all, up-regulated and down-regulated DE genes
nr_DEG <- sapply(res_sig_genes, function(x){
  all <- nrow(x)
  up <- nrow(sigDeg(x, genes = "up"))
  down <- nrow(sigDeg(x, genes = "down"))
  return(cbind(all, up, down))
})
rownames(nr_DEG) <- c("all", "up", "down")

#' Number of all DEGs
barplot(nr_DEG[1,],
        beside = TRUE, 
        main = "Number of all DEGs",
        xlab = "Stage comparison", 
        ylab = "Number of DEGs", 
        names.arg = sub("res_", "", colnames(nr_DEG)), 
        ylim = c(0,20000))

#' Number of up- and down-regulated genes
barplot(nr_DEG[c(2,3), ],
        beside = TRUE,
        main = "Number of up- and down-regulated genes",
        legend.text = TRUE,
        xlab = "Stage comparison",
        ylab = "Number of DEGs",
        names.arg = sub("res_", "", colnames(nr_DEG)), 
        ylim = c(0,12000),
        col = pal12[c(2,3)],
        args.legend = list(x = "top")
)

#' Adapt figure for the manuscript  

# prepare names
stages_man <- sub("res_", "", colnames(nr_DEG))
stages_man <- str_extract_all(stages_man, "[0-9]")
stages_man <- lapply(stages_man, function(s){
    paste0("S", s[2], "-S", s[1])
})

pdf("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/nrDEGS_SE_upDown_consec.pdf", height = 3.2, width = 5, pointsize = 9)
par(mar = c(4,6,3,2))

stages_man <- unlist(stages_man)

barplot(nr_DEG[c(2,3), ],
        beside = TRUE,
        legend.text = c("up-regulated", "down-regulated"),
        #xlab = "Stages",
        #ylab = "Number of DEGs",
        las = 1,
        names.arg = stages_man, 
        ylim = c(0,12000),
        col = pal12[c(3,2)],
        args.legend = list(x = "top", bty = "n")
)
title(ylab= "Number of DEGs", line = 4)
dev.off()


#' ## Number of DEGs in ZE compared to FMG: stages B4-B9
barplot(rbind(group_ZEFMG_up, group_ZEFMG_down),
        beside = TRUE,
        names.arg = c("B4", "B5", "B6", "B7", "B8", "B9"),
        col = c("green3", "darkblue"),
        xlab = "stage of development",
        ylab = "number of DE genes",
        main = "DEGs in ZE compared to FMG",
        ylim = c(0, 2000))
legend("topright",
       legend = c("up-regulated ", "down-regulated"),
       bty="n",
       fill=c("green3", "darkblue"))


#' # Add annotation and export lists of DEGs  
#' Prepare lists of filtered differentially expressed genes with annotation  
#' 
#' Add annotation from ConGenIE
res_annot <- lapply(res_sig_genes, function(x){
  ordx <- x[order(x$padj), ]
  rownames(ordx) <- gsub(".1$", "", rownames(ordx))
  cbind(ordx[ , c("baseMean", "log2FoldChange", "padj")], gene_table[match(rownames(ordx), gene_table$Gene), -1])
})
names(res_annot) <- sub("res_", "DE_genes_", names(res_annot))

#' Write in the files
dev.null <- lapply(names(res_annot), function(x){
  write.table(res_annot[[x]], 
              file = paste0("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/", x, "_padj001_lfc05_annotated.tsv"), 
              quote = FALSE,
              col.names = NA,
              sep = "\t")
})


#' ## Differential expression of genes, comparing all the stages of experiment  
#' 
#' Define pairs of stages to compare
names_stages <- unique(levels(dds$Stages))
#' Reverse the order of the stages, so you will always compare a later stage to an earlier stage in the experiment,
# when comparing expression of the DE genes
comparisons <- combn(rev(names_stages), 2)

#' Obtain DE genes
res_allStages <- apply(comparisons, 2, function(x){
  results(dds, contrast = c("Stages", x), filter = rowMedians(counts(dds)), alpha = alpha)
})

#' Filter them by log2fc and padj  
names(res_allStages) <- apply(comparisons, 2, paste, collapse="~")
res_sig_allStages <- lapply(res_allStages, sigDeg)

# save object with DE genes and TEs in all of the stage comparisons
save(res_sig_allStages, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_AllvsAllStages_padj001_lfc05.rda")

#' Keep only genes, exclude TEs, although there is not much difference as there are only ~300 TEs in the data and even less are DE
res_sig_allStages_genes <- lapply(res_sig_allStages, function(x){
  x[grepl("MA_", rownames(x)), ]
})
save(res_sig_allStages_genes, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_AllvsAllStages_padj001_lfc05.rda")

#' Construct a matrix with number of up- and down-regulated genes (comparisons of all the stages)  
de_mat <- matrix(NA, 
                 nrow = length(names_stages), 
                 ncol = length(names_stages))
# name rows and columns
dimnames(de_mat) <- list(names_stages, names_stages)
# calculate number of up- and down-regulated genes
de_mat[t(comparisons)] <- unlist(lapply(res_sig_allStages_genes, function(f){ - sum(f$log2FoldChange < 0, na.rm = TRUE) }))
de_mat[t(comparisons)[, c(2, 1)]] <- unlist(lapply(res_sig_allStages_genes, function(f){ sum(f$log2FoldChange > 0, na.rm = TRUE) }))

#' Print the heatmap with number of DEGs between all the stages

pdf("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/nrDEgenes_AllvsAllStages_salmon.pdf", width = 12, height = 9)
pheatmap(de_mat, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 18)
dev.off()

png("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/nrDEgenes_AllvsAllStages_salmon.png", width = 700, height = 540)
pheatmap(de_mat, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 18)
dev.off()

#' # Session info
  sessionInfo()