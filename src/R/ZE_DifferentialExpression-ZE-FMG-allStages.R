#' ---
#' title: "Differential expression of ZE vs FMG samples"
#' author: "Katja Stojkovic"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' 

#' # Setup

#' * Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(RColorBrewer))

register(MulticoreParam(2))

#' Function for selecting significant genes (FDR < 0.01, |log2FC| => 0.5)
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
}

# save par
old.par <- par(no.readonly=T)

# colour palette
pal <- brewer.pal(12, "Paired")

#' Import the data
load(here("analysis/salmon/ZE-allStages_duplS_dds.rda"))

#' # Number of expressed genes  
#' Check number of expressed genes per time point - if there would be too much variation that would break the assumptions of DE analysis.  
#' 
# expressed genes per time point (check average expression of all the replicates)
expressedGenes_byStage <- sapply(split.data.frame(t(assay(dds)),colData(dds)$Time), 
                                 function(x){colMeans(x)>0})
expressedGenes_byStage_nr <- apply(expressedGenes_byStage, 2, sum)

# genes (average of replicates) with vst value >1 in any of the time points
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

expressedGenes_vst1_byStage <- sapply(split.data.frame(t(vst), colData(dds)$Time), 
                                      function(x){colMeans(x)>1})
expressedGenes_vst1_byStage_nr <- apply(expressedGenes_vst1_byStage, 2, sum)

expressedGenes_vst5_byStage <- sapply(split.data.frame(t(vst), colData(dds)$Time), 
                                      function(x){colMeans(x)>5})
expressedGenes_vst5_byStage_nr <- apply(expressedGenes_vst5_byStage, 2, sum)

# Plot number of genes expressed or having vst value >1 in different time points
barplot(rbind(expressedGenes_byStage_nr, expressedGenes_vst1_byStage_nr, expressedGenes_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        main = "Genes expressed at different values in ZE", 
        sub = "average expression in all the replicates",
        col = c("black", "grey", "lightblue"))
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")

# Manuscript version of the figure
pdf(here("figures/nrExprGenesZE.pdf"), width = 8, height = 6)
par(mar=c(5, 5, 2.5, 1))
barplot(rbind(expressedGenes_byStage_nr, expressedGenes_vst1_byStage_nr, expressedGenes_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        # main = "Genes expressed at different values in ZE", 
        # sub = "average expression in all the replicates",
        col = c("black", "grey", "lightblue"), 
        names.arg = sub("B", "Z", names(expressedGenes_byStage_nr))
)
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")
dev.off()

#' # Differential expression analysis
# dds <- DESeq(dds, parallel = TRUE)
# message after running function DESeq:  
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# 1 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#
# Options for troubleshooting: https://support.bioconductor.org/p/65091/  
# I will rather delete 1 row than change the workflow to be as similar as it can be to the analysis of SE samples, to which we want to compare the DEGs after DE analysis.
 
# #' Explore these rows
# mcols(dds)[which(mcols(dds)$betaConv == FALSE), ]
# phosphoglycerate mutase family protein 
# 
# #' to remove this row:
# ddsClean <- dds[which(mcols(dds)$betaConv),]
# 
# create directory to save DE results in; the first time
DEout=here("analysis/DE/ZE-FMG-allStages_duplSsamples")
# if(dir.exists(DEout) == FALSE) {dir.create(DEout, recursive = TRUE)}
# save(ddsClean, file = paste0(DEout, "/dds_designTissueTime.rda"))

load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/dds_designTissueTime.rda"))
# ! ddsClean has 55315 elements, not 66069 like the unfiltered one before removing row that did not converge

#' specify alpha
alpha = 0.01


# Extract results -----------------------------------------------------------------  

#' ## Seed B1-B3  
#' ### Extract results   
dds <- ddsClean
resultsNames(dds)

B2.S_B1.S <- results(dds, name = "Time_B2_vs_B1", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B3.S_B2.S <- results(dds, contrast = list("Time_B3_vs_B1", "Time_B2_vs_B1"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)

#' ### Select and count significant DE genes
# All
group_S_all <- sapply(list(B2.S_B1.S, B3.S_B2.S), function(x) nrow(sigDeg(x)))
names(group_S_all) <- c("B2.S_B1.S", "B3.S_B2.S")

# Up-regulated
group_S_up <- sapply(list(B2.S_B1.S, B3.S_B2.S), function(x) nrow(sigDeg(x, genes = "up")))
names(group_S_up) <- names(group_S_all)

# Down-regulated
group_S_down <- sapply(list(B2.S_B1.S, B3.S_B2.S), function(x) nrow(sigDeg(x, genes = "down")))
names(group_S_up) <- names(group_S_all)

#' Table of number of significant DE genes
pander(rbind(group_S_all, group_S_up, group_S_down))

#' ### Write DE genes in the files
# Down-regulated
genes_S_down <-lapply(list(B2.S_B1.S, B3.S_B2.S), function(x) {
  sigx <- sigDeg(x, genes = "down")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_S_down) <- names(group_S_all)

# write in the variable, so the output (e.g. ## [[1]] NULL) is written in the variable instead of in the html document
dev.null <- lapply(names(genes_S_down), function(x) {
  write.table(genes_S_down[[x]],
              file = file.path(DEout,
                               paste0(x, "_down-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# Up-regulated
genes_S_up <-lapply(list(B2.S_B1.S, B3.S_B2.S), function(x) {
  sigx <- sigDeg(x, genes = "up")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_S_up) <- names(group_S_all)

# write in the variable, so the output (e.g. ## [[1]] NULL) is written in the variable instead of in the html document
dev.null <- lapply(names(genes_S_up), function(x) {
  write.table(genes_S_up[[x]],
              file = file.path(DEout,
                               paste0(x, "_up-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# All
genes_S_all <-lapply(list(B2.S_B1.S, B3.S_B2.S), function(x) {
  sigx <- sigDeg(x, genes = "all")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_S_all) <- names(group_S_all)


#' ## FMG B4-B9
#' ### Extract results

B4.FMG_B3.S <- results(dds, contrast = c("Time", "B4", "B3"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B5.FMG_B4.FMG <- results(dds, contrast = c("Time", "B5", "B4"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B6.FMG_B5.FMG <- results(dds, contrast = c("Time", "B6", "B5"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B7.FMG_B6.FMG <- results(dds, contrast = c("Time", "B7", "B6"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B8.FMG_B7.FMG <- results(dds, contrast = c("Time", "B8", "B7"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B9.FMG_B8.FMG <- results(dds, contrast = c("Time", "B9", "B8"), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)

#' ### Select significant DE genes
FMGcompStage <- c(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG)
FMGcompStage_names <- c("B4.FMG_B3.S", "B5.FMG_B4.FMG", "B6.FMG_B5.FMG", "B7.FMG_B6.FMG", "B8.FMG_B7.FMG", "B9.FMG_B8.FMG")

# All
group_FMG_all <- sapply(list(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG), function(x) nrow(sigDeg(x)))
names(group_FMG_all) <- FMGcompStage_names
# Up-regulated
group_FMG_up <- sapply(list(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG), function(x) nrow(sigDeg(x, genes = "up")))
names(group_FMG_up) <- FMGcompStage_names
# Down-regulated
group_FMG_down <- sapply(list(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG), function(x) nrow(sigDeg(x, genes = "down")))
names(group_FMG_up) <- FMGcompStage_names

#' Table of number of significant DE genes
pander(rbind(group_FMG_all, group_FMG_up, group_FMG_down))

#' ### Write DE genes in the files
# Down-regulated
genes_FMG_down <-lapply(list(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG), function(x) {
  sigx <- sigDeg(x, genes = "down")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_FMG_down) <- FMGcompStage_names

dev.null <- lapply(names(genes_FMG_down), function(x) {
  write.table(genes_FMG_down[[x]],
              file = file.path(DEout,
                               paste0(x, "_down-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# Up-regulated
genes_FMG_up <-lapply(list(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG), function(x) {
  sigx <- sigDeg(x, genes = "up")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_FMG_up) <- FMGcompStage_names

# write in the variable, so the output (e.g. ## [[1]] NULL) is written in the variable instead of in the html document
dev.null <- lapply(names(genes_FMG_up), function(x) {
  write.table(genes_FMG_up[[x]],
              file = file.path(DEout,
                               paste0(x, "_up-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# All
genes_FMG_all <-lapply(list(B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG), function(x) {
  sigx <- sigDeg(x, genes = "all")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_FMG_all) <- FMGcompStage_names

#' ## ZE B4-B9
#' ### Extract results
B4.ZE_B3.S <- results(dds, contrast = list(c("Time_B4_vs_B1", "TissueZE.TimeB4"), c("Time_B3_vs_B1", "TissueZE.TimeB3")), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B5.ZE_B4.ZE <- results(dds, contrast = list(c("Time_B5_vs_B1", "TissueZE.TimeB5"), c("Time_B4_vs_B1", "TissueZE.TimeB4")), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B6.ZE_B5.ZE <- results(dds, contrast = list(c("Time_B6_vs_B1", "TissueZE.TimeB6"), c("Time_B5_vs_B1", "TissueZE.TimeB5")), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B7.ZE_B6.ZE <- results(dds, contrast = list(c("Time_B7_vs_B1", "TissueZE.TimeB7"), c("Time_B6_vs_B1", "TissueZE.TimeB6")), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B8.ZE_B7.ZE <- results(dds, contrast = list(c("Time_B8_vs_B1", "TissueZE.TimeB8"), c("Time_B7_vs_B1", "TissueZE.TimeB7")), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B9.ZE_B8.ZE <- results(dds, contrast = list(c("Time_B9_vs_B1", "TissueZE.TimeB9"), c("Time_B8_vs_B1", "TissueZE.TimeB8")), alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)

#' ### Select significant DE genes
ZEcompStage <- c(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE)
ZEcompStage_names <- c("B4.ZE_B3.S", "B5.ZE_B4.ZE", "B6.ZE_B5.ZE", "B7.ZE_B6.ZE", "B8.ZE_B7.ZE", "B9.ZE_B8.ZE")
# All
group_ZE_all <- sapply(list(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE), function(x) nrow(sigDeg(x)))
names(group_ZE_all) <- ZEcompStage_names
# Up-regulated
group_ZE_up <- sapply(list(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE), function(x) nrow(sigDeg(x, genes = "up")))
names(group_ZE_up) <- ZEcompStage_names

# Down-regulated
group_ZE_down <- sapply(list(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE), function(x) nrow(sigDeg(x, genes = "down")))
names(group_ZE_up) <- ZEcompStage_names

#' Table of number of significant DE genes
pander(rbind(group_ZE_all, group_ZE_up, group_ZE_down))

#' ### Write DE genes in the files
# Down-regulated
genes_ZE_down <-lapply(list(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE), function(x) {
  sigx <- sigDeg(x, genes = "down")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_ZE_down) <- ZEcompStage_names


dev.null <- lapply(names(genes_ZE_down), function(x) {
  write.table(genes_ZE_down[[x]],
              file = file.path(DEout,
                               paste0(x, "_down-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# Up-regulated
genes_ZE_up <-lapply(list(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE), function(x) {
  sigx <- sigDeg(x, genes = "up")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_ZE_up) <- ZEcompStage_names

dev.null <- lapply(names(genes_ZE_up), function(x) {
  write.table(genes_ZE_up[[x]],
              file = file.path(DEout,
                               paste0(x, "_up-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# All
genes_ZE_all <-lapply(list(B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE), function(x) {
  sigx <- sigDeg(x, genes = "all")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_ZE_all) <- ZEcompStage_names

#' ## ZE vs FMG B4-B9
#' ### Extract results
B4.ZE_B4.FMG <- results(dds, name = "TissueZE.TimeB4", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B5.ZE_B5.FMG <- results(dds, name = "TissueZE.TimeB5", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B6.ZE_B6.FMG <- results(dds, name = "TissueZE.TimeB6", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B7.ZE_B7.FMG <- results(dds, name = "TissueZE.TimeB7", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B8.ZE_B8.FMG <- results(dds, name = "TissueZE.TimeB8", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)
B9.ZE_B9.FMG <- results(dds, name = "TissueZE.TimeB9", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE)

#' ### Select significant DE genes
ZEFMGcompStage <- c(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG) 
ZEFMGcompStage_names <- c("B4.ZE_B4.FMG", "B5.ZE_B5.FMG", "B6.ZE_B6.FMG", "B7.ZE_B7.FMG", "B8.ZE_B8.FMG", "B9.ZE_B9.FMG")

# All
group_ZEFMG_all <- sapply(list(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG), function(x) nrow(sigDeg(x)))
names(group_ZEFMG_all) <- ZEFMGcompStage_names
# Up-regulated
group_ZEFMG_up <- sapply(list(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG), function(x) nrow(sigDeg(x, genes = "up")))
names(group_ZEFMG_up) <- ZEFMGcompStage_names
# Down-regulated
group_ZEFMG_down <- sapply(list(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG), function(x) nrow(sigDeg(x, genes = "down")))
names(group_ZEFMG_up) <- ZEFMGcompStage_names

#' Table of number of significant DE genes
pander(rbind(group_ZEFMG_all, group_ZEFMG_up, group_ZEFMG_down))

#' ### Write DE genes in the files
# Down-regulated
genes_ZEFMG_down <-lapply(list(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG), function(x) {
  sigx <- sigDeg(x, genes = "down")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_ZEFMG_down) <- ZEFMGcompStage_names

dev.null <- lapply(names(genes_ZEFMG_down), function(x) {
  write.table(genes_ZEFMG_down[[x]],
              file = file.path(DEout,
                               paste0(x, "_down-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# Up-regulated
genes_ZEFMG_up <-lapply(list(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG), function(x) {
  sigx <- sigDeg(x, genes = "up")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_ZEFMG_up) <- ZEFMGcompStage_names

dev.null <- lapply(names(genes_ZEFMG_up), function(x) {
  write.table(genes_ZEFMG_up[[x]],
              file = file.path(DEout,
                               paste0(x, "_up-regulated_one-percent-FDR-05-log2fc-cutoff_significant-genes.tsv")),
              col.names = NA,
              sep = "\t")
})

# All
genes_ZEFMG_all <-lapply(list(B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG), function(x) {
  sigx <- sigDeg(x, genes = "all")
  ordx <- sigx[order(sigx$log2FoldChange), ]
  ordx[, c("baseMean", "log2FoldChange", "padj")]
})
names(genes_ZEFMG_all) <- ZEFMGcompStage_names

#' ## ZE vs FMG, regardless time  
#' ### Extract results
ZE_FMG <- results(dds, name = "Tissue_ZE_vs_FMG", alpha = alpha, filter = rowMedians(counts(dds)), parallel = TRUE) 

#' ### Select significant DE genes  
genes_ZE_FMG_all <- nrow(sigDeg(ZE_FMG, genes = "all"))
genes_ZE_FMG_all

#' # Export objects  
#' Save objects with filtered DEGs for later exploration
save(genes_S_all, genes_FMG_all, genes_ZE_all, genes_ZEFMG_all, 
     file = file.path(DEout, "filteredDEGsInZE_all_one-percent-FDR-05-log2fc-cutoff.rda"))
save(genes_S_up, genes_S_down, genes_FMG_up, genes_FMG_down, genes_ZE_up, genes_ZE_down, genes_ZEFMG_up, genes_ZEFMG_down, 
     file = file.path(DEout, "filteredDEGsInZE_UpAndDown_one-percent-FDR-05-log2fc-cutoff.rda")) 
#' Save objects with DEGs, not filtered by padj and lfc; add also the intercept  
# get the intercept
resIntercept <- results(dds, name = "Intercept", alpha = alpha, filter = rowMedians(counts(dds)))

save(resIntercept, 
     B2.S_B1.S, B3.S_B2.S, 
     B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG, 
     B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE,
     B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG, 
     file = file.path(DEout, "notfilteredDEGs.rda"))

#' # Figures  
#' ## DEGs in S comparing consecutive stages of the same tissue
# All
barplot(rbind(group_S_all),
         beside = TRUE,
         names.arg = c("B2-B1", "B3-B2"),
         ylim = c(0, 6000),
         main = "Nr of DE genes between stages in different tissues",
         col = "grey")

# Up- and down-regulated
barplot(rbind(group_S_up,group_S_down),
         beside = TRUE,
         names.arg = c("B2-B1", "B3-B2"),
         main = "Nr of up- and down-regulated genes between stages in the seed",
         ylim = c(0, 3000),
         col = c("green3", "darkblue"))
legend("top", legend = c("up-regulated genes", "down-regulated genes"), fill = c("green3", "darkblue"))


#' ## Number of DEGs in ZE and FMG: comparing consecutive stages of the same tissue
#'
# All
barplot(rbind(group_ZE_all, group_FMG_all),
         beside = TRUE,
         names.arg = c("B4-B3", "B5-B4", "B6-B5", "B7-B6", "B8-B7", "B9-B8"),
         ylim = c(0, 7000),
         main = "Nr of DE genes between stages in different tissues",
         col = c("darkorange", "gold"))
legend("top", legend = c("ZE", "FMG"), fill = c("salmon", "gold"))

# Up- and -down-regulated
barplot(rbind(group_ZE_up, group_FMG_up, group_ZE_down, group_FMG_down),
         beside = TRUE,
         names.arg = c("B4-B3", "B5-B4", "B6-B5", "B7-B6", "B8-B7", "B9-B8"),
         ylim = c(0, 3500),
         legend.text = FALSE,
         main = "Nr of DE genes in consecutive stages of ZE and FMG",
         col = c("lightgreen", "darkgreen", "lightblue2", "darkblue"))

legend("topright",
       legend = c("ZE: up-regulated ", "FMG: up-regulated", "ZE: down-regulated", "FMG: down-regulated"),
       bty="n",
       fill=c("lightgreen", "darkgreen", "lightblue2", "darkblue"))

#' version of the figure for the manuscript  
# paste all the data together (S, FMG, Zem)
bind_ZemFG <- rbind(group_ZE_up, group_FMG_up, group_ZE_down, group_FMG_down)
bind_all_ZE <- cbind(rbind(group_S_up,group_S_down, c(NA, NA), c(NA, NA), c(NA, NA), c(NA, NA)), rbind(c(NA, NA, NA, NA, NA, NA), c(NA, NA, NA, NA, NA, NA), bind_ZemFG))
rownames(bind_all_ZE) <- c("group_S_up","group_S_down", "group_ZE_up", "group_FMG_up", "group_ZE_down", "group_FMG_down")

# seed
pdf(here("figures/nrDEGs_SD_consec.pdf"), height = 2.2, width = 2.5, pointsize = 8)
par(mar = c(4, 5, 3, 2))
barplot(rbind(group_S_up,group_S_down),
        beside = TRUE,
        names.arg = c("Z1-Z2", "Z2-Z3"),
        ylim = c(0, 4000),
        ylab = "", 
        legend.text = FALSE,
        las = 1,
        space = c(0, 0, 1.3, 0),
        width = rep(0.10, length(rbind(group_S_up,group_S_down))),
        #xlim = c(0,3),
        # main = "Nr of DE genes in consecutive stages of ZE",
        col = c("gold", "chocolate4")
)
title(ylab = "number of DEGs", line = 3.5)
dev.off()

# Zem and FG
pdf(here("figures/nrDEGs_ZemFG_consec.pdf"), height = 2.2, width = 4.4, pointsize = 8)
par(mar = c(4, 2, 3, 2))
barplot(bind_ZemFG,
        beside = TRUE,
        names.arg = c("Z3-Z4", "Z4-Z5", "Z5-Z6", "Z6-Z7", "Z7-Z8", "Z8-Z9"),
        # add = TRUE,
        ylim = c(0, 4000), 
        yaxt = "n",
        space = c(0, rep(c(0, 0.4, 0, 1.3), 5), 0, 0.4, 0),
        legend.text = FALSE, 
        # width = c(rep(2, length(bind_ZemFG))),
        # main = "Nr of DE genes in consecutive stages of ZE",
        col = c(pal[4], pal[3], pal[2], pal[1])
)

legend("topright",
       legend = c("SD: up-regulated", "SD: down-regulated", "Zem: up-regulated", "Zem: down-regulated", "FG: up-regulated", "FG: down-regulated"),
       bty="n",
       fill=c("#F4FC0D", "#0F4365", pal[4], pal[2], pal[3], pal[1]))
dev.off()

# join SD, FG and Zem for publication
pdf(here("figures/nrDEGs_allZE_consec.pdf"), height = 3.3, width = 6.6, pointsize = 8)
par(mar = c(5, 5, 3, 2))
barplot(bind_all_ZE,
        beside = TRUE,
        names.arg = c("Z1-Z2", "Z2-Z3", "Z3-Z4", "Z4-Z5", "Z5-Z6", "Z6-Z7", "Z7-Z8", "Z8-Z9"),
        # add = TRUE,
        ylim = c(0, 4000), 
        space = c(0, rep(c(0, 0, 0, 0, 0, 0), 2), rep(c(0, 0, 0, 0.4, 0, 0), 5), 0, 0, 0, 0.4, 0),
        legend.text = FALSE,
        # main = "Nr of DE genes in consecutive stages of ZE",
        col = c(rep(c("gold", "chocolate4"), 6), rep(c("black", "black", pal[4], pal[3], pal[2], pal[1]), 7)),
        xlab = "stage of development",
        ylab = "number of DE genes"
)
legend("topright",
       legend = c("SD: up-regulated", "SD: down-regulated", "Zem: up-regulated", "Zem: down-regulated", "FG: up-regulated", "FG: down-regulated"),
       bty="n",
       fill=c("gold", "chocolate4", pal[4], pal[2], pal[3], pal[1]))
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

#' Further exploration of DEGs in the study, compared also to somatic embryogenesis, can be found in the script ExploringDEGs-ZE-FMG-SE.R



#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
