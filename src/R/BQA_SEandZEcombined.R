#' ---
#' title: "Comparison of SE and ZE samples"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Setup  

#' Load libraries
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(stringr))

#' Source some helper functions
source(here("UPSCb-common/src/R/plot.multidensity.R"))
source(here("UPSCb-common/src/R/featureSelection.R"))

#' Load data  

# Count tables, unnormalised
count_ZE <- read.csv(here("analysis/salmon/ZE-allStages_unnormalised-counts.csv"),
                       row.names = 1)

load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes+TEs.rda") 
dds_SE <- dds
rm(dds)
# exclude TEs from this dataset, keep only genes (there are 291 TEs and ~66000 genes)  
dds_SE <- dds_SE[grepl("MA_", rownames(dds_SE)), ]
count_SE <- counts(dds_SE)

#' Sample information
sample_ZE <- read.csv(here("doc/sample_info.csv"),
                       row.names = 1)
rownames(sample_ZE) <- sample_ZE$NGI.ID
sample_SE <- data.frame(colData(dds_SE))

#' Create palettes
pal <- brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
pal2 <- brewer.pal(8, "Set2")
pal_multi <- c(pal12, pal2, pal)
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' # Aim  
#' Check how similar are samples from somatic and zygotic embryogenesis.  
#'  
#' # Analysis  
#' ## Prepare and join the data from different projects  
#' If rownames are the same (they should be), join all count tables
if(all(rownames(count_SE) == rownames(count_ZE))) count_all <- cbind(count_SE, count_ZE)

#' Prepare and join sample information
sample_ZE_prep <- sample_ZE[ , c("Time", "Tissue", "Replicate", "NGI.ID")]
sample_ZE_prep$Experiment <- rep("ZE", nrow(sample_ZE_prep))

sample_SE_prep <- data.frame(Time = sample_SE[, "Stages"],
                             Tissue = rep("SE", nrow(sample_SE)),
                             Replicate = paste0(sample_SE$Stages, "-SE"),
                             NGI.ID = sample_SE$ID,
                             Experiment = rep("SE", nrow(sample_SE)))
rownames(sample_SE_prep) <- sample_SE_prep$NGI.ID


if(all(colnames(sample_SE_prep) == colnames(sample_ZE_prep))) sample_all <- rbind(sample_SE_prep, sample_ZE_prep)
rownames(sample_all) <- sample_all$NGI.ID

write.csv(sample_all, here("doc/sampleInfo_SEZE.csv"))

#' ## Data normalisation  
#'  
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
#' 
#' Create the dds object, without giving any prior on the design
dds <- DESeqDataSetFromMatrix(
  countData = count_all,
  colData = sample_all,
  design = ~Replicate)

#' Check the size factors (i.e. the sequencing library size effect)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count_all)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' Variance Stabilising Transformation  
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Validate the VST 
meanSdPlot(vst[rowSums(count_all)>0,])

#' ## QC on the normalised data
#' 
#' ### PCA
".pca" <- function(vst,fact,lgd="bottomright",pal=pal_multi){
  pc <- prcomp(t(vst))
  
  percent <- round(summary(pc)$importance[2,]*100)
  
  #' #### First 3 dimensions
  mar=c(5.1,4.1,4.1,2.1)
  scatterplot3d(pc$x[,1],
                pc$x[,2],
                pc$x[,3],
                xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                color=pal[as.integer(fact)],
                pch=19)
  legend(lgd,pch=19,
         col=pal_multi[1:nlevels(fact)],
         legend=levels(fact))
  par(mar=mar)
}

#' #### Time
.pca(vst,factor(sample_all$Time),lgd="topright",pal=c(1,pal_multi))

#' #### 1st and 2nd dims
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal_multi[as.integer(factor(sample_all$Time))],
     pch=c(17, 7, 19)[as.integer(factor(sample_all$Tissue))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("top",bty="n",col=pal_multi,levels(factor(sample_all$Time)),pch=19)
legend("bottomleft", bty = "n", pch=c(17, 7, 19), levels(factor(sample_all$Tissue)))

#' #### 2nd and 3rd dims

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal_multi[as.integer(factor(sample_all$Time))],
     pch=c(17, 7, 19)[as.integer(factor(sample_all$Tissue))],
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal_multi,levels(factor(sample_all$Time)),pch=19)
legend("topleft", bty = "n", pch=c(17, 7, 19), levels(factor(sample_all$Tissue)))


#' ### Heatmap
sel <- featureSelect(vst,factor(sample_all$Time),exp = 7,nrep = 2)
message(sprintf("We select %s genes for the heatmap",sum(sel)))
heatmap.2(vst[sel,],trace = "none",labRow = FALSE, 
          col=hpal, 
          ColSideColors = c(1,pal_multi)[as.integer(factor(sample_all$Time))],
          scale="row")

legend("bottomright", bty="n", col=c(1,pal_multi), levels(factor(sample_all$Time)),pch=19)

#' ### Sample Hierarchical clustering  
#' change labels for publication
sample_ZEpubl <- read.csv(here("doc/sample_info_abbrevAsInPaper.csv"), row.names = 1)
sample_ZEpubl <- sample_ZEpubl[ , c("Time", "Tissue", "Replicate", "NGI.ID", "TimeTissue")]
sample_ZEpubl <- cbind(sample_ZEpubl[ , c(1:4)], "Experiment" = rep("ZE", nrow(sample_ZEpubl)), "TimeTissue" = sample_ZEpubl$TimeTissue)
rownames(sample_ZEpubl) <- sample_ZEpubl$NGI.ID

# change SE names first
sample_SEpubl <- sample_SE_prep
sample_SEpubl$Replicate <- sub("-SE", "", sample_SE_prep$Replicate)
sample_SEpubl$TimeTissue <- sample_SEpubl$Replicate

# join SE and ZE sample info
sample_publ <- rbind(sample_ZEpubl, sample_SEpubl)
write.csv(sample_all, here("doc/sampleInfo_SEZE_publ.csv"))

# order samples as in vst
sample_publ <- sample_publ[colnames(vst), ]

pdf(here("figures/hierarchClustering_SEandZEjointNormalisation.pdf"), height = 5, width = 10)
hc <- hclust(dist(t(vst)))
plot(hc, 
     labels=sample_publ$TimeTissue,
     cex=0.7)
#rect.hclust(hc, k = 6, border = c("seagreen", "seagreen", "gold", "seagreen", "chocolate", "sandybrown"))
#legend = (c(SE, Zem)
dev.off()

png(here("figures/hierarchClustering_SEandZEjointNormalisation.png"), units = "in", res = 150, height = 5, width = 10)
hc <- hclust(dist(t(vst)))
plot(hc, 
     labels=sample_publ$TimeTissue,
     cex=0.7)
#rect.hclust(hc, k = 6, border = c("seagreen", "seagreen", "gold", "seagreen", "chocolate", "sandybrown"))
#legend = (c(SE, Zem)
dev.off()

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
