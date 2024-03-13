#' title: "Compare analysis of genes with or without TEs"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Aim  
#' Get PCA plots for publications. 
#' 
#'  Setup  
#' ## Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))

#' ## DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes.rda")

#' ## Functions
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
    if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
    if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
    if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

#' ## Graphics
pal <- brewer.pal(12,"Paired")
mar <- par("mar")

#' # Analysis of dataset with only genes included  
#' ## VST  
vsd_genes <- varianceStabilizingTransformation(dds_genes, blind=TRUE)
vst_genes <- assay(vsd_genes)
vst_genes <- vst_genes - min(vst_genes)

vsda_genes <- varianceStabilizingTransformation(dds_genes, blind=FALSE)
vsta_genes <- assay(vsda_genes)
vsta_genes <- vsta_genes - min(vsta_genes)

write.csv(vst_genes, "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/geneOnly_vst_blind.csv", quote = FALSE)
write.csv(vsta_genes, "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/geneOnly_vst_aware.csv", quote = FALSE)
write.csv2(as.data.frame(colData(dds_genes)), "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/sample_info.csv", quote = FALSE)

#' ## PCA  
#' 
#' figure for manuscript, blind PCA  
pc_gb <- prcomp(t(vst_genes))
percent_gb <- round(summary(pc_gb)$importance[2,]*100)

barplot(percent_gb, ylim = c(0,40), main = "% of variation explained by PCs")
plot(cumsum(percent_gb), 
     xlab = "Number of PCs", 
     ylab = "Percent variation", 
     main = "Variation explained by PCs (cummulative)")

pdf("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/PCA_blind_SE.pdf", height = 6, width = 8)
par(mar=c(5, 5, 2, 5), xpd = TRUE)
plot(pc_gb$x[,1],
     pc_gb$x[,2],
     xlab=paste("PC 1 (",percent_gb[1],"%)",sep=""),
     ylab=paste("PC 2 (",percent_gb[2],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds_genes)$Stages))],
     pch=19,
     cex = 2,
     cex.axis = 1.5,
     cex.lab = 1.5
)

legend(270, 180, bty = "n", col=pal[1:nlevels(as.factor(colData(dds_genes)$Stages))], levels(as.factor(colData(dds_genes)$Stages)), pch=19, cex = 1.3)
dev.off()

# PCA aware
pc_g <- prcomp(t(vsta_genes))
percent_g <- round(summary(pc_g)$importance[2,]*100)

barplot(percent_g, ylim = c(0,40), main = "% of variation explained by PCs")
plot(cumsum(percent_g), 
     xlab = "Number of PCs", 
     ylab = "Percent variation", 
     main = "Variation explained by PCs (cummulative)")

plot(pc_g$x[,1],
     pc_g$x[,2],
     xlab=paste("Comp. 1 (",percent_g[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent_g[2],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds_genes)$Stages))],
     main="1st and 2nd PC",
     pch=19)
legend("top",
       pch=19,
       col=pal[1:nlevels(as.factor(colData(dds_genes)$Stages))],
       legend=levels(as.factor(colData(dds_genes)$Stages)))

plot(pc_g$x[,2],
     pc_g$x[,3],
     xlab=paste("Comp. 2 (",percent_g[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent_g[3],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds_genes)$Stages))],
     main="2nd and 3rd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(as.factor(colData(dds_genes)$Stages))],
       legend=levels(as.factor(colData(dds_genes)$Stages)))

#' ## DE
dds_genes <- DESeq(dds_genes)

res_list_genes <- lapply(2:8, function(i){
    results(dds_genes, c("Stages", paste0("S", as.character(i)), paste0("S", as.character(i-1))), 
            filter = rowMedians(counts(dds_genes)))
})

names(res_list_genes) <- lapply(2:8, function(i){
    paste0("res_", i, "vs", i-1)
})

#' Filter
res_sig_genes <- lapply(res_list_genes, sigDeg)

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
