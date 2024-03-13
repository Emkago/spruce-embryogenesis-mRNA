#' ---
#' title: "Biological quality check of samples from spruce zygotic embryogenesis"
#' author: "Katja Stojkovic"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---


#' # Setup

#' * Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))

#' * Helper functions
source(here("UPSCb-common/src/R/plot.multidensity.R"))
source(here("UPSCb-common/src/R/featureSelection.R"))
  
#' Function for selecting significant genes (FDR < 0.01, |log2FC| => 0.5)
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
}
  
#' * Graphics
pal <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information 
samples <- read_csv(here("doc/ZE_Dataset_v4.csv"),
                    col_types = cols(col_character(),
                                     col_character(),
                                     col_factor(),
                                     col_character(),
                                     col_factor(),
                                     col_character(),
                                     col_double(),
                                     col_double())) %>% 
  mutate(Tissue=factor(Tissue)) %>% 
  mutate(Time=factor(Time))

# correct typo
samples$User.ID[samples$NGI.ID == "P11562_201"] <- sub("_S", "", samples$User.ID[samples$NGI.ID == "P11562_201"])

#' sample information including technical replicates (sequencing lanes)
samples_techRep <- samples[rep(seq_len(nrow(samples)), each = 2), ]
samples_techRep$Lane <- rep(c("L001","L002"), 57)
samples_techRep$TechRep <- paste(samples_techRep$NGI.ID, samples_techRep$Lane, sep = "_")


#' # Prepare data  
#' ## Raw data
filelist <- list.files(here("data/RNA-Seq/salmon"), 
                          recursive = TRUE,
                          pattern = "quant.sf",
                          full.names = TRUE)

#' Filter files from folder tmp
filelist <- str_subset(filelist, "tmp", negate = TRUE)

names(filelist) <- str_extract(filelist, "P11562_.+_L00[1,2]")
names(filelist) <- sub("_S[0-9]+", "", names(filelist))

setdiff(names(filelist), samples_techRep$TechRep)
setdiff(samples_techRep$TechRep, names(filelist))

# One technical replicate is missing for the sample 204 (738-FMG).
samples_techRep <- filter(samples_techRep, !grepl("P11562_204_L002",TechRep))

stopifnot(all(str_which(names(filelist), samples_techRep$TechRep) == 1:length(filelist)))

# Exclude sample representing time point B10 as there is only one replicate
samples_techRep <- filter(samples_techRep, !grepl("B10",Time))
filelist <- filelist[names(filelist) %in% samples_techRep$TechRep]
samples <- filter(samples, !grepl("B10", Time))

#' Read the expression at the gene level (there is one transcript per gene)
geneExpr <- suppressMessages(tximport(files = filelist, 
                                  type = "salmon",txOut=TRUE))
counts_techRep <- round(geneExpr$counts)

#' Export TPM to display them in Plantgenie
TPM_techRep <- geneExpr$abundance

# BQA -----------------------------------------------------------------  

#' ## Check technical replicates  
#' ### Prepare DESeq2 object  
dds_techRep <- DESeqDataSetFromMatrix(
countData = counts_techRep,
colData = samples_techRep,
design =  ~Replicate)

# As using design ~Tissue*Time gives an error, I will use ~Replicate, as I am interested in blind = TRUE evaluation. 
# Later, after joining technical replicates and before doing DE I will use design model ~Tissue*Time  

# Error in checkFullRank(modelMatrix) : 
#   the model matrix is not full rank, so the model cannot be fit as specified.
# Levels or combinations of levels without any samples have resulted in
# column(s) of zeros in the model matrix.
# 
# Please read the vignette section 'Model matrix not full rank':
#   
#   vignette('DESeq2')  

#' ### PCA  
#' Transformation of the data for visualisation purposes
vsd_t <- varianceStabilizingTransformation(dds_techRep, blind = TRUE)
vst_t <- assay(vsd_t)
vst_t <- vst_t - min(vst_t)

meanSdPlot(vst_t[rowSums(counts_techRep)>0,])

#' PCA
pc <- prcomp(t(vst_t))
percent <- round(summary(pc)$importance[2,]*100)

plot(cumsum(percent))

pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(vsd_t)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=TechRep)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Outlier  
#' Sample B4 is an outlier, it is in the group of samples that represent mature seeds. Exclude it from further analysis.  
# remove P11562_112 from count table 
counts_techRep <- counts_techRep[ , str_subset(colnames(counts_techRep), "P11562_112_L00*", negate = TRUE)]
# remove P11562_112 from sample list
samples_techRep <- filter(samples_techRep, !grepl("112",NGI.ID))
samples <- filter(samples, !grepl("112",NGI.ID))

stopifnot(all(str_which(colnames(counts_techRep), samples_techRep$TechRep) == 1:length(colnames(counts_techRep))))


#' ### Join technical replicates
counts <- sapply(split.data.frame(t(counts_techRep),samples_techRep$NGI.ID),colSums)
stopifnot(all(str_which(colnames(counts), samples$NGI.ID) == 1:length(colnames(counts))))

#' save counts and sample information
write.csv(counts, file = here("analysis/salmon/ZE-allStages_unnormalised-counts.csv"), quote = FALSE)
write.csv(samples, file = here("doc/sample_info.csv"), quote = FALSE)
write.csv(samples_techRep, file = here("doc/sample_info_techRep.csv"), quote = FALSE)

# sample information with abbreviations used in the paper
samples_paper <- data.frame(samples)
samples_paper$Replicate <- sub("FMG", "FG", samples_paper$Replicate)
samples_paper$Replicate <- sub("B", "Z", samples_paper$Replicate)
samples_paper$Replicate <- sub("S", "SD", samples_paper$Replicate)
samples_paper$Replicate <- sub("ZE", "Zem", samples_paper$Replicate)
samples_paper$Replicate <- sub("-", ".", samples_paper$Replicate)
samples_paper$Time <- sub("B", "Z", samples_paper$Time)
samples_paper$Tissue <- sub("FMG", "FG", samples_paper$Tissue)
samples_paper$Tissue <- sub("S", "SD", samples_paper$Tissue)
samples_paper$Tissue <- sub("ZE", "Zem", samples_paper$Tissue)
samples_paper$TimeTissue <- paste0(samples_paper$Time,".", samples_paper$Tissue)

write.csv(samples_paper, file = here("doc/sample_info_abbrevAsInPaper.csv"))

#' remove outlier and join replicates for TPM values
TPM_techRep <- TPM_techRep[ , str_subset(colnames(TPM_techRep), "P11562_112_L00*", negate = TRUE)]
stopifnot(all(str_which(colnames(TPM_techRep), samples_techRep$TechRep) == 1:length(colnames(TPM_techRep))))

TPM <- sapply(split.data.frame(t(TPM_techRep),samples_techRep$NGI.ID),colSums)
stopifnot(all(str_which(colnames(TPM), samples$NGI.ID) == 1:length(colnames(TPM))))

write.csv(TPM, file = here("doc/ZE-allStages_TPM.csv"), quote = FALSE)


#' ## Analysis of biological replicates  
#' ### Prepare DESeq object (design Replicate)  
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples_paper,
                              design =  ~Replicate)

#' #### PCA  
#' Transformation of the data for visualisation purposes
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

vsda <- varianceStabilizingTransformation(dds, blind = FALSE)
vsta <- assay(vsda)
vsta <- vsta - min(vsta)

meanSdPlot(vst[rowSums(counts)>0,])

#' PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

plot(cumsum(percent))

# PCA ggplot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(vsd)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Time,shape=Tissue,text=NGI.ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

# PCA base R
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

pdf(here("analysis/figures/PCA_blind_ZE.pdf"), height = 6, width = 8)
par(mar=c(5, 5, 2, 5), xpd = TRUE)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("PC 1 (",percent[1],"%)",sep=""),
     ylab=paste("PC 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples_paper$Time))],
     pch=c(17, 7, 19)[as.integer(factor(samples_paper$Tissue))],
     cex = 2,
     cex.axis = 1.5,
     cex.lab = 1.5
     #main="Principal Component Analysis",sub="variance stabilized counts"
     )

legend(130, 150, bty = "n", col=pal, levels(factor(samples_paper$Time)), pch=19, cex = 1.3)
legend(130, -50, bty = "n", pch=c(17, 7, 19), levels(factor(samples_paper$Tissue)), cex = 1.3)
dev.off()

# PCA base R, design aware (~Replicates)
pca <- prcomp(t(vsta))
percenta <- round(summary(pca)$importance[2,]*100)

pdf(here("analysis/figures/PCA_aware_designReplicate_ZE.pdf"), height = 6, width = 8)
par(mar=c(5, 5, 2, 5), xpd = TRUE)
plot(pca$x[,1],
     pca$x[,2],
     xlab=paste("Comp. 1 (",percenta[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percenta[2],"%)",sep=""),
     col=pal[as.integer(factor(samples_paper$Time))],
     pch=c(17, 7, 19)[as.integer(factor(samples_paper$Tissue))],
     #main="Principal Component Analysis",sub="variance stabilized counts, aware"
     )

legend(75, 100, bty = "n", col=pal, levels(factor(samples_paper$Time)), pch=19)
legend(75, 0, bty = "n", pch=c(17, 7, 19), levels(factor(samples_paper$Tissue)))
dev.off()


# dds: duplicated S, ~Tissue*Time----------------------------------------------------  

# use samples with changed names  
samples <- samples_paper

#' ### Prepare DESeq object (design ~Tissue*Time)  
#' Duplicate samples from first three stages and rename them from S to FMG and ZE  
#' Change sample information
samples_S <- samples[samples$Tissue == "SD", ]
# S to FMG
samples_StoFMG <- samples_S
samples_StoFMG[ , c("Tissue", "Replicate")] <- apply(samples_S[ , c("Tissue", "Replicate")], 2, function(x){gsub("SD", "FG", x)})
samples_StoFMG$User.ID <- paste(samples_StoFMG$User.ID, samples_StoFMG$Tissue, sep = "_")
# S to ZE
samples_StoZE <- samples_S
samples_StoZE[ , c("Tissue", "Replicate")] <- apply(samples_S[ , c("Tissue", "Replicate")], 2, function(x){gsub("SD", "Zem", x)})
samples_StoZE$User.ID <- paste(samples_StoZE$User.ID, samples_StoZE$Tissue, sep = "_")
# make up new NGI.ID for additional samples created above (as first samples start with 1 and replaced samples start with 2, start these with 3)
samples_StoZE$NGI.ID <- sub("_1", "_3", samples_StoZE$NGI.ID)
samples_StoZE$NGI.ID <- sub("_201", "_309", samples_StoZE$NGI.ID)
# join all in new sample info
samples_duplS <- rbind(samples_StoFMG, samples_StoZE, samples[samples$Tissue != "SD", ]) 

#' Change count table  
# duplicate counts and give them new NGI.IDs
counts_StoZE <- counts[ ,samples$Tissue == "SD"]
colnames(counts_StoZE) <- samples_StoZE$NGI.ID 
# join all in new count table
counts_duplS <- cbind(counts, counts_StoZE)

#' Reorder counts to match sample info
counts_duplS <- counts_duplS[ , match(samples_duplS$NGI.ID, colnames(counts_duplS))]

#' Save counts and sample info with duplicated "S" samples
# write.csv(counts_duplS, file = here("analysis/salmon/ZE-allStages_unnormalised-counts_duplS.csv"), quote = FALSE)
# write.csv(samples_duplS, file = here("doc/sample_info_duplS.csv"), quote = FALSE)
write.csv(counts_duplS, file = here("analysis/salmon/ZE-allStages_unnormalised-counts_duplS_msNames.csv"), quote = FALSE)
write.csv(samples_duplS, file = here("doc/sample_info_duplS_msNames.csv"), quote = FALSE)

#' Prepare DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = counts_duplS,
  colData = samples_duplS,
  design = ~ Tissue * Time)

# check factor levels
levels(dds$Tissue)
levels(dds$Time)

#' Save dds object
save(dds, file = here("analysis/salmon/ZE-allStages_duplS_dds.rda"))

load(here("analysis/salmon/ZE-allStages_duplS_dds.rda"))
# export normalised counts
vsd_f <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_f <- assay(vsd_f)
vst_f <- vst_f - min(vst_f)
write.csv(vst_f, here("analysis/salmon/ZE-allStages_duplS_vst_blind.csv"), quote = FALSE)

vsda_f <- varianceStabilizingTransformation(dds, blind = FALSE)
vsta_f <- assay(vsda_f)
vsta_f <- vsta_f - min(vsta_f)
write.csv(vsta_f, here("analysis/salmon/ZE-allStages_duplS_vst_aware.csv"), quote = FALSE)

#'
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

