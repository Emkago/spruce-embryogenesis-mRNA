#' ---
#' title: "Exploring DEGs in ZE, FMG and SE samples"
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
#' 
#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(venn))
suppressPackageStartupMessages(library(SuperExactTest))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(networkD3))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(htmlwidgets))

#' Get helper files
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R"))
source("~/Git/zygoticEmbryogenesis/Rtoolbox/src/plotEnrichedTreemap.R")
source("~/Git/zygoticEmbryogenesis/Rtoolbox/src/infomapTools.R")

#' * Graphics
pal <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' Import data  
#' * DEGs in ZE
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_all_one-percent-FDR-05-log2fc-cutoff.rda"))
load(here("analysis/DE/ZE-FMG-allStages_duplSsamples/filteredDEGsInZE_UpAndDown_one-percent-FDR-05-log2fc-cutoff.rda"))

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

#' * Expression of genes in ZE
load(here("analysis/salmon/ZE-allStages_duplS_dds.rda"))
dds_ZE <- dds
rm(dds)

#' # Analysis   

 
# ----------------------------------------Number of expressed genes -------------------------------------  


#' ## Number of expressed genes  
#' ### SE 
#' Check number of expressed genes per time point for the SE experiment- if there would be too much variation that would break the assumptions of DE analysis.  
# expressed genes per time point (check average expression of all the replicates)
expressedGenesSE_byStage <- sapply(split.data.frame(t(assay(dds_SE)),colData(dds_SE)$Stages), 
                                   function(x){colMeans(x)>0})
expressedGenesSE_byStage_nr <- apply(expressedGenesSE_byStage, 2, sum)

# genes with vst value >1 in any of the time points
vsd_SE <- varianceStabilizingTransformation(dds_SE, blind = TRUE)
vst_SE <- assay(vsd_SE)
vst_SE <- vst_SE - min(vst_SE)

expressedGenesSE_vst1_byStage <- sapply(split.data.frame(t(vst_SE), colData(dds_SE)$Stages), 
                                        function(x){colMeans(x)>1})
expressedGenesSE_vst1_byStage_nr <- apply(expressedGenesSE_vst1_byStage, 2, sum)

expressedGenesSE_vst5_byStage <- sapply(split.data.frame(t(vst_SE), colData(dds_SE)$Stages), 
                                        function(x){colMeans(x)>5})
expressedGenesSE_vst5_byStage_nr <- apply(expressedGenesSE_vst5_byStage, 2, sum)

# Plot number of genes expressed or having vst value >1 in different time points
barplot(rbind(expressedGenesSE_byStage_nr, expressedGenesSE_vst1_byStage_nr, expressedGenesSE_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        main = "Genes expressed at different values in SE", 
        sub = "average expression in all the replicates",
        col = c("black", "grey", "lightblue"))
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")

# Manuscript version of the figure
pdf(here("figures/nrExprGenesSE.pdf"), width = 8, height = 6)
par(mar=c(5, 5, 2.5, 1))
barplot(rbind(expressedGenesSE_byStage_nr, expressedGenesSE_vst1_byStage_nr, expressedGenesSE_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        col = c("black", "grey", "lightblue"), 
        #names.arg = sub("B", "Z", names(expressedGenes_byStage_nr))
)
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")
dev.off()

#' ### ZE 
#' Check number of expressed genes per time point for the ZE experiment- if there would be too much variation that would break the assumptions of DE analysis.  
# 
# expressed genes per time point (check average expression of all the replicates)
expressedGenesZE_byStage <- sapply(split.data.frame(t(assay(dds_ZE)),colData(dds_ZE)$Time), 
                                 function(x){colMeans(x)>0})
expressedGenesZE_byStage_nr <- apply(expressedGenesZE_byStage, 2, sum)

# genes (average of replicates) with vst value >1 in any of the time points
vsd_ZE <- varianceStabilizingTransformation(dds_ZE, blind = TRUE)
vst_ZE <- assay(vsd_ZE)
vst_ZE <- vst_ZE - min(vst_ZE)

expressedGenesZE_vst1_byStage <- sapply(split.data.frame(t(vst_ZE), colData(dds_ZE)$Time), 
                                      function(x){colMeans(x)>1})
expressedGenesZE_vst1_byStage_nr <- apply(expressedGenesZE_vst1_byStage, 2, sum)

expressedGenesZE_vst5_byStage <- sapply(split.data.frame(t(vst_ZE), colData(dds_ZE)$Time), 
                                      function(x){colMeans(x)>5})
expressedGenesZE_vst5_byStage_nr <- apply(expressedGenesZE_vst5_byStage, 2, sum)

# Plot number of genes expressed or having vst value >1 in different time points
barplot(rbind(expressedGenesZE_byStage_nr, expressedGenesZE_vst1_byStage_nr, expressedGenesZE_vst5_byStage_nr), 
        ylim = c(0, 60000), 
        beside = TRUE, 
        xlab = "time points",
        ylab = "number of expressed genes",
        main = "Genes expressed at different values in ZE", 
        sub = "average expression in all the replicates",
        col = c("black", "grey", "lightblue"))
legend("topleft", fill = c("black", "grey", "lightblue"),legend = c("raw count >0", "vst >1", "vst >5"), bty = "n")


#-----------------------------------------------UpSet plots------------------------------------------


#' ## Intersections of DEGs: UpSet plot  
#' Set the commonly used parameters
upset_default <- function(x){upset(fromList(lapply(x, rownames)), sets = names(x), keep.order = TRUE, sets.bar.color = "grey", nsets = length(x), nintersects = NA)}
# nsets = "nr of sets", others only 5 will be shown,  
# nintersects = NA to display all of them, otherwise max 40 will be shown  
# for keep.order to work, sets need to be specified; an order of sets provided to sets is kept  

#' ### DEGs in SE  
upset_default(ddsDEGs_SE)

#' ### DEGS in ZE and FMG (separately for each tissue), through time  

#' FMG: all DE genes
upset_default(genes_FMG_all)
#' ZE: all DE genes
upset_default(genes_ZE_all)

#' FMG: up-regulated genes
upset_default(genes_FMG_up)
#' FMF: down-regulated genes
upset_default(genes_FMG_down)

#' ZE: up-regulated genes
upset_default(genes_ZE_up)
#' ZE: down-regulated genes
upset_default(genes_ZE_down)

### DEGs in ZE and FMG, together
genes_ZE_experiment_all <- lapply(c(1:6), function(x){
  union(rownames(genes_FMG_all[[x]]), rownames(genes_ZE_all[[x]]))
})

elementNROWS(genes_ZE_experiment_all)

#' ### ZE vs FMG in each time point  

#' ZE vs FMG: all DEGs
upset_default(genes_ZEFMG_all)
#' ZE vs FMG: down-regulated genes
upset_default(genes_ZEFMG_down)
#' ZE vs FMG: up-regulated genes
upset_default(genes_ZEFMG_up)


#' How does union of the genes above (from every timepoint comparison) compares to Tissue_ZE_vs_FMG?
#' 

#' ### Comparisson of DEGs at any point in SE, ZE,FMG and S (make union of all DEGs in a tissue regardless of time)  
# Pool together DEGs from all timepoint comparissons in each tissue
allDEGs_S <- unique(unlist(lapply(genes_S_all, rownames), use.names = FALSE))
allDEGs_FMG <- unique(unlist(lapply(genes_FMG_all, rownames), use.names = FALSE))
allDEGs_ZE <- unique(unlist(lapply(genes_ZE_all, rownames), use.names = FALSE))

# DEGs from SE
allDEGs_SE <- unique(unlist(lapply(ddsDEGs_SE, rownames), use.names = FALSE))

#' Intersections of all DEGs regardless of time in different datasets/tissues
DEGsListAllGeneral <- list(seed = allDEGs_S, zygotic_embryo = allDEGs_ZE, FMG = allDEGs_FMG, somatic_embryo = allDEGs_SE)
venn::venn(DEGsListAllGeneral, zcolor = "style", box = FALSE)
upset(fromList(DEGsListAllGeneral), sets = names(DEGsListAllGeneral), keep.order = TRUE, sets.bar.color = "grey", nsets = length(DEGsListAllGeneral), nintersects = NA)

# # edit for the publication
# DEGsListAllGeneral <- list(SD = allDEGs_S, Zem = allDEGs_ZE, FG = allDEGs_FMG, SE = allDEGs_SE)
# 
# pdf(here("figures/Intersect_DEGs_SEZE_allTissues_venn.pdf"), width = 3.15, height = 3)
# venn::venn(DEGsListAllGeneral, 
#            zcolor = c("chocolate4", "sandybrown", "gold", "seagreen"), 
#            box = FALSE,
#            opacity = 0.5,
#            ilcs = 0.7, 
#            sncs = 0.7)
# dev.off()
# 
# upset(fromList(DEGsListAllGeneral), sets = names(DEGsListAllGeneral), keep.order = TRUE, sets.bar.color = "grey", nsets = length(DEGsListAllGeneral), nintersects = NA)

# all the tissues in ZE together vs SE
venn::venn(list(somatic_embryogenesis = DEGsListAllGeneral$somatic_embryo,
                zygotic_embryogenesis = union(DEGsListAllGeneral$seed, DEGsListAllGeneral$zygotic_embryo, DEGsListAllGeneral$FMG)),
           zcolor = "style", box = FALSE)
# what about intersection of DEGs only in tissues from zygotic embryogenesis
venn:venn(list(seed = DEGsListAllGeneral$seed, zygotic_embryo = DEGsListAllGeneral$zygotic_embryo, FMG = DEGsListAllGeneral$FMG), zcolor = "style", box = FALSE)

#' Can we get intersections of DEGs between all the stage comparisons of different kinds of the tissue?  
# combine lists of DEGs from all the tissues
DEGsListAllComp <- lapply(list(genes_S_all, genes_FMG_all, genes_ZE_all, ddsDEGs_SE), function(x){
  lapply(x, rownames)})
DEGsListAllComp <- unlist(DEGsListAllComp, recursive = FALSE)
names(DEGsListAllComp) <- sub(pattern = "res_", replacement = "SE_", names(DEGsListAllComp))
save(DEGsListAllComp, file=here("analysis/DE/ZE-FMG-allStages_duplSsamples/DEGsListAllComp_ZEandSE.rda"))

# plot intersections arranged by number of DEGs in the intersection, decreasing; first 80 intersections by size
upset(fromList(DEGsListAllComp),
      sets = names(DEGsListAllComp),
      keep.order = TRUE,
      order.by = "freq",
      sets.bar.color = "grey", 
      nsets = length(DEGsListAllComp), 
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison", 
      mb.ratio = c(0.5, 0.5),
      nintersects = 80,
      number.angles = 45)

#' ### Comparison of DEGs at two major transitions in SE and ZE  
#' Subset data  
DEGspeaks_SEZE <- DEGsListAllComp[c("B4.FMG_B3.S", "B7.FMG_B6.FMG",
                  "B4.ZE_B3.S","B7.ZE_B6.ZE",
                  "SE_3vs2", "SE_6vs5")]

upset(fromList(DEGspeaks_SEZE),
      sets = names(DEGspeaks_SEZE),
      keep.order = TRUE,
      sets.bar.color = "grey", 
      nsets = length(DEGspeaks_SEZE), 
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison", 
      #mb.ratio = c(0.5, 0.5),
      nintersects = NA,
      number.angles = 25)
# use set.metadata if you want to color the rows of the matrix  


#-------------------------------------------Background population-----------------------------------------  


#' ## Intersections of DEGs: SuperExactTest plots  
#' 
#' How likely is that intersection of this size would appear by chance (because of the size of the group)? 
#' Is it expected that this many genes would appear the same, comparing two groups of genes with certain background population size?  
#' Calculated p value tells us one-tail probability of observing equal to or larger than the number of intersect items.  

#' Comparison of setup & SET plots:  
#' In setup plot only elements specific for that group are represented, but in SET plots number of all the elements in the group is counted! 
#' (there is a possibility to choose from three different modes in setup plots, default mode is "distinct"; "intersect" and "union" also available))  

#' ### Choosing background population of genes  
#'  
#' How to choose a background population size (n)?  
#' From the vignette: in their case n = total = 18196; cis-eQTL gene sets were independently and randomly sampled from the population of 18,196 unique genes profiled in the eQTL study  

#' We should use a number of expressed genes in the dataset (experiment). When comparing two experiments (SE & ZE) use a union of expressed genes.  
#' Comment (Nico): It is easier to assume that the population consists of any gene expressed in that dataset (as we do for GO). It might vary slightly between tissues and time point, but not drastically as shown above. 
#' So it won’t affect the stats much.  

#' #### vst  

#' Calculate number of expressed genes in SE and ZE:
# i) SE
# all the genes with expression > 0 in at least one time point
expressedGenesSE_anyStage <- apply(expressedGenesSE_byStage, 1, any)
# how many are they?
SEbackground_0 <- sum(expressedGenesSE_anyStage)

# all the genes with vst > 1 in at least one time point
expressedGenesSE_vst1_anyStage <- apply(expressedGenesSE_vst1_byStage, 1, any)
# how many are they?
SEbackground_1 <- sum(expressedGenesSE_vst1_anyStage)

# ii) ZE
# all the genes with expression > 0 in at least one time point
expressedGenesZE_anyStage <- apply(expressedGenesZE_byStage, 1, any)
# how many are they?
ZEbackground_0 <- sum(expressedGenesZE_anyStage)

# all the genes with vst > 1 in at least one time point
expressedGenesZE_vst1_anyStage <- apply(expressedGenesZE_vst1_byStage, 1, any)
# how many are they?
ZEbackground_1 <- sum(expressedGenesZE_vst1_anyStage)

#' Number of expressed genes in SE and ZE experiments (raw expression > 0 or vst >1; average of replicates in at least one time)
pander(data.frame("background_0" = c(ZEbackground_0, SEbackground_0), "background_1" = c(ZEbackground_1, SEbackground_1), row.names = c("ZE", "SE")))

#' VST expression is relative to the dataset/experiment, as we subtract the minimium to get vst values and it is therefore not directly comparable between experiments. 
#' Could we use a different measure of expression?  
#' * Estimate noise and filter expressed genes (check the script in UPSCb-comon)  
#' * Count all the genes taken for doing DE analysis; that would be all the genes not discarded (padj is not NA) by independent filtering in DESeq2.  
#' 
#' #### DESeq2 filtering  
#'   
#' Import not filtered DE results for SE and ZE in order to get number of expressed genes from DESeq objects  
#' ZE
load(here("data/analysis/DE/ZE-FMG-allStages_duplSsamples/notfilteredDEGs.rda"))
# Number of genes (whose padj is not NA) reported in results of each of the comparison is not the same
listZEstages <- lapply(list(resIntercept, B2.S_B1.S, B3.S_B2.S, 
                     B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG,
                     B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE),
                     function(x) {x[!is.na(x$padj), ]})
names(listZEstages) <- c("Intercept", "B2.S_B1.S", "B3.S_B2.S", 
                         "B4.FMG_B3.S", "B5.FMG_B4.FMG", "B6.FMG_B5.FMG", "B7.FMG_B6.FMG", "B8.FMG_B7.FMG", "B9.FMG_B8.FMG",
                         "B4.ZE_B3.S", "B5.ZE_B4.ZE", "B6.ZE_B5.ZE", "B7.ZE_B6.ZE", "B8.ZE_B7.ZE", "B9.ZE_B8.ZE")

NrExpGeneZEstage <- sapply(listZEstages, nrow)

pander(NrExpGeneZEstage)
pander(table(NrExpGeneZEstage))
upset_default(listZEstages)

# Names of all the expressed genes in ZE (union of expressed genes in all the comparisons):
NameExpGeneZE <- unique(unlist(lapply(listZEstages, rownames), use.names = FALSE))
# Number of expressed genes in zygotic embryogenesis (after independent filtering) is therefore:
NrExpGeneZE <- length(NameExpGeneZE)
NrExpGeneZE
# This represents % of all genes (note: number of genes before filtering is the same in any comparison: intercept, B2.S_B3.S, ...)
NrExpGeneZE/66069*100

# Could we use intercept as the approximation of all the expressed genes in the experiment? Check also in other experiments!
nrow(listZEstages$Intercept)
# check the threshold for independent filtering: % of removed genes and the cutoff value for expression (baseMean)
metadata(resIntercept)$filterThreshold
# ! stats above represent percent of 55315 genes in the ddsClean, although initial number of the genes in the dds was equat to the number of all 
# the predicted genes in the spruce genome, that is 66069  
# recalculate stats to report percent of excluded genes:
nr_excl_genes_ZE <- as.numeric(str_extract(names(metadata(resIntercept)$filterThreshold), "[0-9]+.[0-9]+"))*nrow(resIntercept)/100
(nr_excl_genes_ZE + (66069-nrow(resIntercept)))/66069*100
# remove unnecessary objects
rm(resIntercept, 
   B2.S_B1.S, B3.S_B2.S, 
   B4.FMG_B3.S, B5.FMG_B4.FMG, B6.FMG_B5.FMG, B7.FMG_B6.FMG, B8.FMG_B7.FMG, B9.FMG_B8.FMG, 
   B4.ZE_B3.S, B5.ZE_B4.ZE, B6.ZE_B5.ZE, B7.ZE_B6.ZE, B8.ZE_B7.ZE, B9.ZE_B8.ZE,
   B4.ZE_B4.FMG, B5.ZE_B5.FMG, B6.ZE_B6.FMG, B7.ZE_B7.FMG, B8.ZE_B8.FMG, B9.ZE_B9.FMG)

#' Export names of the expressed genes to use them as a background in functional enrichment analysis  
save(NameExpGeneZE, file = here("analysis/gopher/backgroundGenesZE.rda"))

#' SE  
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_unfiltered.rda")
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genesTEs_unfiltered_onlyIntercept.rda")

listSEstages <- res_list
listSEstages[["Intercept"]] <- resInterceptSE

# exclude TEs
listSEstages <- lapply(listSEstages, function(x){
  x[grepl("^MA_", rownames(x)), ]
})

# filter for expressed genes
listSEstages_filt <- lapply(listSEstages, function(x){
  x[!is.na(x$padj), ]
})

NrExpGeneSEstage <- sapply(listSEstages_filt, nrow)

pander(NrExpGeneSEstage)
pander(table(NrExpGeneSEstage))
upset_default(listSEstages_filt)

# Names of all the expressed genes in SE (union of expressed genes in all the comparisons):
NameExpGeneSE <- unique(unlist(lapply(listSEstages_filt, rownames), use.names = FALSE))
# Number of expressed genes in zygotic embryogenesis (after independent filtering) is therefore:
NrExpGeneSE <- length(NameExpGeneSE)
NrExpGeneSE
# This represents % of all genes (note: number of genes before filtering is the same in any comparison: intercept, S2 vs S1, ... = 66069)
NrExpGeneSE/66069*100

# Export names of the expressed genes to use them as a background in functional enrichment analysis  
save(NameExpGeneSE, file = here("analysis/gopher/backgroundGenesSE.rda"))

# Could we use intercept as the approximation of all the expressed genes in the experiment? Check also in other experiments!
nrow(listSEstages_filt$Intercept)
# check the threshold for independent filtering: % of removed genes and the cutoff value for expression (baseMean)
metadata(listSEstages$Intercept)$filterThreshold
metadata(listSEstages_filt$Intercept)$filterThreshold

# remove unnecessary objects
rm(resInterceptSE, res_list)


# --------------------------------------SuperExactTest plots----------------------------------------------------  


#' ### Expressed genes  
#' #### SE  
setSEexpr <- supertest(lapply(listSEstages_filt, rownames), n = NrExpGeneSE)
plot(setSEexpr,
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:3),#degree = c(2:length(listSEstages_filt)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()

dev.new()
plot(setSEexpr,
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(4:5),#degree = c(2:length(listSEstages_filt)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()

dev.new()
plot(setSEexpr,
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(6:length(listSEstages_filt)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()

#' #### ZE: not useful due to too many combinations, keep the UpSet plot?  
#' 

#' ### DEGs in SE  
#' Define size of background population
n = NrExpGeneSE

#' All DE genes
setSE_all <- supertest(lapply(ddsDEGs_SE, rownames), n = n)
dev.new()
plot(setSE_all,
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(setSE_all)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()

#' ### DEGS in ZE and FMG (separately for each tissue), through time  
#' Define size of background population
n = NrExpGeneZE

#' FMG: all DE genes
setFMG_all <- supertest(lapply(genes_FMG_all, rownames), n = n)
dev.new()
plot(setFMG_all,
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(setFMG_all)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
    color.on = "black", color.off = "white")
dev.off()

# plot P values
# plot(setFMG_all$P.value)
# abline(h=0.05, col = "red")
# summary(setFMG_all)$Table[summary(setFMG_all)$P.value < 0.05, ]

#' ZE: all DE genes
setZE_all <- supertest(lapply(genes_ZE_all, rownames), n = n)
dev.new()
plot(setZE_all, 
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(setZE_all)),
     mar= c(2, 6, 2, 3), Layout = "landscape", 
     color.on = "black", color.off = "white")
dev.off()

# plot P values
# plot(setZE_all$P.value)
# abline(h=0.05, col = "red")
# summary(setZE_all)$Table[summary(setZE_all)$P.value < 0.05, ]

#' FMG: up-regulated genes
dev.new()
plot(supertest(lapply(genes_FMG_up, rownames), n = n),
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(genes_FMG_up)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()
#' FMG: down-regulated genes
dev.new()
plot(supertest(lapply(genes_FMG_down, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(genes_FMG_down)),
     Layout = "landscape", mar= c(2, 8, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()

#' ZE: up-regulated genes
dev.new()
plot(supertest(lapply(genes_ZE_up, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(genes_ZE_up)),
     Layout = "landscape", mar= c(1, 6, 2, 3), 
     color.on = "black", color.off = "white")
dev.off

#' ZE: down-regulated genes
dev.new()
plot(supertest(lapply(genes_ZE_down, rownames), n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(genes_ZE_down)),
     Layout = "landscape", mar= c(1, 6, 2, 3), 
     color.on = "black", color.off = "white")
dev.off()
     
     
#' ### ZE vs FMG in each time point  

#' ZE vs FMG: all DEGs
plot(supertest(lapply(genes_ZEFMG_all, rownames), n = n), 
     sort.by = c("size"), keep.empty.intersections = FALSE, degree = c(2:length(genes_ZEFMG_all)),
     Layout = "landscape", mar= c(1, 7, 2, 3), intersection.size.rotate=TRUE, #################### intersection.size.rotate does not work #################
     color.on = "black", color.off = "white")
#' ZE vs FMG: down-regulated genes
plot(supertest(lapply(genes_ZEFMG_down, rownames), n = n), 
     sort.by = c("size"), keep.empty.intersections = FALSE, degree = c(2:length(genes_ZEFMG_down)),
     Layout = "landscape", mar= c(1, 7, 2, 3), intersection.size.rotate=TRUE,
     color.on = "black", color.off = "white")
#' ZE vs FMG: up-regulated genes
plot(supertest(lapply(genes_ZEFMG_up, rownames), n = n), 
     sort.by = c("size"), keep.empty.intersections = FALSE, degree = c(2:length(genes_ZEFMG_up)),
     Layout = "landscape", mar= c(1, 7, 2, 3), intersection.size.rotate=TRUE,
     color.on = "black", color.off = "white")

#' #### Nr of genes after independent filtering in DESeq2  
venn::venn(list(SE = NameExpGeneSE, ZE = NameExpGeneZE),
           zcolor = "style", box = FALSE)

#' Union of genes expressed in SE and ZE
NameExpGeneSEZE <- union(NameExpGeneSE, NameExpGeneZE)
#' Number of genes in the union
length(NameExpGeneSEZE)

#' Intersections of all DEGs regardless of time in different datasets/tissues  
#' Define size of background population (union of unique expressed genes in SE and ZE)
n = length(NameExpGeneSEZE)

plot(supertest(DEGsListAllGeneral, n = n), 
     sort.by = "size", keep.empty.intersections = FALSE, degree = c(2:length(DEGsListAllGeneral)),
     Layout = "landscape", mar= c(1, 7, 2, 3),
     color.on = "black", color.off = "white")


#' ### Comparison of DEGs at two major transitions in SE and ZE  
dev.new()
plot(supertest(DEGspeaks_SEZE, n = n), 
     sort.by = c("size", "degree"), keep.empty.intersections = FALSE, degree = c(2:length(DEGspeaks_SEZE)),
     Layout = "landscape", mar= c(1, 7, 2, 3),
     color.on = "black", color.off = "white")
dev.off()

DEGsListTransSEupDown <- list("S2-S3_up" = rownames(ddsDEGs_SE$res_3vs2[ddsDEGs_SE$res_3vs2$log2FoldChange > 0, ]),
                              "S2-S3_down" = rownames(ddsDEGs_SE$res_3vs2[ddsDEGs_SE$res_3vs2$log2FoldChange < 0, ]),
                              "S5-S6_up" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange > 0, ]),
                              "S5-S6_down" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange < 0, ]))

upset(fromList(DEGsListTransSEupDown),
      keep.order = TRUE,
      sets.bar.color = "grey",
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison",
      nintersects = NA,
      number.angles = 25)

# include also S7vsS6, as it is similar, although not represented anymore in ZE
upset(fromList(list("S2-S3_up" = rownames(ddsDEGs_SE$res_3vs2[ddsDEGs_SE$res_3vs2$log2FoldChange > 0, ]),
                    "S2-S3_down" = rownames(ddsDEGs_SE$res_3vs2[ddsDEGs_SE$res_3vs2$log2FoldChange < 0, ]),
                    "S5-S6_up" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange > 0, ]),
                    "S5-S6_down" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange < 0, ]),
                    "S6-S7_up" = rownames(ddsDEGs_SE$res_7vs6[ddsDEGs_SE$res_7vs6$log2FoldChange > 0, ]),
                    "S6-S7_down" = rownames(ddsDEGs_SE$res_7vs6[ddsDEGs_SE$res_7vs6$log2FoldChange < 0, ]))),
      keep.order = TRUE,
      sets.bar.color = "grey",
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison",
      nintersects = NA,
      number.angles = 25)

lapply(list("S2-S3_up" = rownames(ddsDEGs_SE$res_3vs2[ddsDEGs_SE$res_3vs2$log2FoldChange > 0, ]),
     "S2-S3_down" = rownames(ddsDEGs_SE$res_3vs2[ddsDEGs_SE$res_3vs2$log2FoldChange < 0, ]),
     "S5-S6_up" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange > 0, ]),
     "S5-S6_down" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange < 0, ]),
     "S6-S7_up" = rownames(ddsDEGs_SE$res_7vs6[ddsDEGs_SE$res_7vs6$log2FoldChange > 0, ]),
     "S6-S7_down" = rownames(ddsDEGs_SE$res_7vs6[ddsDEGs_SE$res_7vs6$log2FoldChange < 0, ])),
     length)

DEGsListTransZEupDown <- list("Z4.FG-Z3.S3_up" = rownames(genes_FMG_up$B4.FMG_B3.S),
                              "Z4.FG-Z3.S3_down" = rownames(genes_FMG_down$B4.FMG_B3.S),
                              "Z4.ZE-Z3.S3_up" = rownames(genes_ZE_up$B4.ZE_B3.S),
                              "Z4.ZE-Z3.S3_down" = rownames(genes_ZE_down$B4.ZE_B3.S),
                              "Z7.FG-Z6.FG_up" = rownames(genes_FMG_up$B7.FMG_B6.FMG),
                              "Z7.FG-Z6.FG_down" = rownames(genes_FMG_down$B7.FMG_B6.FMG),
                              "Z7.ZE-Z6.ZE_up" = rownames(genes_ZE_up$B7.ZE_B6.ZE),
                              "Z7.ZE-Z6.ZE_down" = rownames(genes_ZE_down$B7.ZE_B6.ZE))

upset(fromList(DEGsListTransZEupDown),
      nsets = 8,
      keep.order = TRUE,
      sets.bar.color = "grey",
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison",
      nintersects = NA,
      number.angles = 25)

# Intersection of S5-S6 and Z6-Z7 in both tissues, FG and Ze; note changed names (not technical; S5-S6 means: which genes were downregulated from S5 to S6)
DEGsListDesiccation <- list("S5-S6_up" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange > 0, ]),
                            "S5-S6_down" = rownames(ddsDEGs_SE$res_6vs5[ddsDEGs_SE$res_6vs5$log2FoldChange < 0, ]),
                            "Z6.FG-Z7.FG_up" = rownames(genes_FMG_up$B7.FMG_B6.FMG),
                            "Z6.FG-Z7.FG_down" = rownames(genes_FMG_down$B7.FMG_B6.FMG),
                            "Z6.Ze-Z7.Ze_up" = rownames(genes_ZE_up$B7.ZE_B6.ZE),
                            "Z6.Ze-Z7.Ze_down" = rownames(genes_ZE_down$B7.ZE_B6.ZE))

upset(fromList(DEGsListDesiccation),
      keep.order = TRUE,
      sets.bar.color = "grey",
      mainbar.y.label = "Number of DEGs in the intersection", 
      sets.x.label = "Number of DEGs in the comparison",
      nintersects = NA,
      number.angles = 25)

# ###### to make upset plot more readable, use set.metadata (colour the rows and ) 
# set.metadata = list(data = metadata, plots = list(list(type = "hist", 
# column = "avgRottenTomatoesScore", assign = 20), list(type = "matrix_rows", 
#                                                       column = "Cities", colors = c(Boston = "green", NYC = "navy", LA = "purple"), 
#                                                       alpha = 0.5)))  

#===============================Sankey diagram: 2 major transition points==============================  

#' Represent changes in differential expression of genes that are DE in the first and second major transition point with Sankey diagram  
# get the data
upset_trans_ZE_upDown <- upset(fromList(DEGsListTransZEupDown), 
                           nsets = 8,
                           keep.order = TRUE,
                           nintersects = NA) 
  
# DEGsListTransZEupDown_unlisted <- unlist(DEGsListTransZEupDown, use.names = FALSE)
# DEGsListTransZEupDown_unlisted <- DEGsListTransZEupDown_unlisted[ !duplicated(DEGsListTransZEupDown_unlisted) ]
# # Now we know what does the rownumbers from "New_data" refer to in our list, if we want to extract lists of genes in the comparisons.  

# Subset data to contain only combinations of 2 sets
upsetData_trans_ZE_upDown_degree2 <- upset_trans_ZE_upDown$New_data[rowSums(upset_trans_ZE_upDown$New_data) == 2, ]
# Subset data to contain only combinations where DEGs are present in one of the sets from the first stage (Z4-Z3)
transZ4 <- colnames(upsetData_trans_ZE_upDown_degree2[grepl("Z4", colnames(upsetData_trans_ZE_upDown_degree2))])
upsetData_trans_ZE_upDown_degree2 <- upsetData_trans_ZE_upDown_degree2[apply(upsetData_trans_ZE_upDown_degree2, 1, function(x) {any(x[transZ4] == 1)}), ]

# Reformat the data to include "source", "target", "value" (early stage, late stage, nr DEGs)
transZ7 <- colnames(upsetData_trans_ZE_upDown_degree2[grepl("Z7", colnames(upsetData_trans_ZE_upDown_degree2))])

upsetDataReformat_trans_ZE_upDown_degree2 <- lapply(transZ4, function(x){
  comb <- apply(upsetData_trans_ZE_upDown_degree2, 1, function(y){
    y[x] == 1 & any(y[transZ7] == 1)
  })
subData <- upsetData_trans_ZE_upDown_degree2[comb, ]
freqData <- as.data.frame(table(subData))
freqData <- freqData[freqData$Freq != 0, ]
})

upsetDataReformat_trans_ZE_upDown_degree2 <- as.data.frame(bind_rows(upsetDataReformat_trans_ZE_upDown_degree2))

# Build a connection data frame (a list of flows with intensity for each flow)
links <- data.frame(source = apply(upsetDataReformat_trans_ZE_upDown_degree2, 1, function(x){colnames(upsetDataReformat_trans_ZE_upDown_degree2)[1:4][x[1:4] == 1]}),
                    target = apply(upsetDataReformat_trans_ZE_upDown_degree2, 1, function(x){colnames(upsetDataReformat_trans_ZE_upDown_degree2)[5:8][x[5:8] == 1]}),
                    value = upsetDataReformat_trans_ZE_upDown_degree2$Freq)

# Create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 20)
p

# save the widget
saveWidget(p, file=here("figures/IntersectDEGs_2majorTransitions_ZE_upDown_FlowOfChangingDEGs.html"))

#' As genes do not "flow" from one tissue to another (though some of the RNA molecules and proteins could be transported), 
#' separate diagrams for ZEm and FG
links_ZEm <- links[grepl("ZE", links$source) & grepl("ZE", links$target), ]
links_FG <- links[grepl("FG", links$source) & grepl("FG", links$target), ]

nodes_ZEm <- data.frame(
  name=c(as.character(links_ZEm$source), 
         as.character(links_ZEm$target)) %>% unique())
nodes_FG <- data.frame(
  name=c(as.character(links_FG$source), 
         as.character(links_FG$target)) %>% unique())

links_ZEm$IDsource <- match(links_ZEm$source, nodes_ZEm$name)-1 
links_ZEm$IDtarget <- match(links_ZEm$target, nodes_ZEm$name)-1

links_FG$IDsource <- match(links_FG$source, nodes_FG$name)-1 
links_FG$IDtarget <- match(links_FG$target, nodes_FG$name)-1

p_ZEm <- sankeyNetwork(Links = links_ZEm, Nodes = nodes_ZEm,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 20)
p_FG <- sankeyNetwork(Links = links_FG, Nodes = nodes_FG,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name", 
                       sinksRight=FALSE, fontSize = 20)

p_ZEm
p_FG

#' Sankey diagram for DEGs in SE  
# get the data
upset_trans_SE_upDown <- upset(fromList(DEGsListTransSEupDown), 
                               nsets = length(DEGsListTransSEupDown),
                               keep.order = TRUE,
                               nintersects = NA) 

# Subset data to contain only combinations of 2 sets
upsetData_trans_SE_upDown_degree2 <- upset_trans_SE_upDown$New_data[rowSums(upset_trans_SE_upDown$New_data) == 2, ]
# Subset data to contain only combinations where DEGs are present in one of the sets from the first stage (S2-S3)
transS3 <- colnames(upsetData_trans_SE_upDown_degree2[grepl("S3", colnames(upsetData_trans_SE_upDown_degree2))])
upsetData_trans_SE_upDown_degree2 <- upsetData_trans_SE_upDown_degree2[apply(upsetData_trans_SE_upDown_degree2, 1, function(x) {any(x[transS3] == 1)}), ]

# Reformat the data to include "source", "target", "value" (early stage, late stage, nr DEGs)
transS6 <- colnames(upsetData_trans_SE_upDown_degree2[grepl("S6", colnames(upsetData_trans_SE_upDown_degree2))])

upsetDataReformat_trans_SE_upDown_degree2 <- lapply(transS3, function(x){
  comb <- apply(upsetData_trans_SE_upDown_degree2, 1, function(y){
    y[x] == 1 & any(y[transS6] == 1)
  })
  subData <- upsetData_trans_SE_upDown_degree2[comb, ]
  freqData <- as.data.frame(table(subData))
  freqData <- freqData[freqData$Freq != 0, ]
})

upsetDataReformat_trans_SE_upDown_degree2 <- as.data.frame(bind_rows(upsetDataReformat_trans_SE_upDown_degree2))

# Build a connection data frame (a list of flows with intensity for each flow)
links_SE <- data.frame(source = apply(upsetDataReformat_trans_SE_upDown_degree2, 1, function(x){colnames(upsetDataReformat_trans_SE_upDown_degree2)[1:2][x[1:2] == 1]}),
                    target = apply(upsetDataReformat_trans_SE_upDown_degree2, 1, function(x){colnames(upsetDataReformat_trans_SE_upDown_degree2)[3:4][x[3:4] == 1]}),
                    value = upsetDataReformat_trans_SE_upDown_degree2$Freq)

# Create a node data frame: it lists every entities involved in the flow
nodes_SE <- data.frame(
  name=c(as.character(links_SE$source), 
         as.character(links_SE$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe. So we need to reformat it.
links_SE$IDsource <- match(links_SE$source, nodes_SE$name)-1 
links_SE$IDtarget <- match(links_SE$target, nodes_SE$name)-1

# Make the Network
p_SE <- sankeyNetwork(Links = links_SE, Nodes = nodes_SE,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 20)
p_SE


# --------------------------------------------Heatmaps and enrichment------------------------------------------------  


#' # Expression of common genes in different tissues  
#' ## Design aware vst  
#' SE  
# vsda_SE <- varianceStabilizingTransformation(dds_SE, blind = FALSE)
# vsta_SE <- assay(vsda_SE)
# vsta_SE <- vsta_SE - min(vsta_SE)
# write.csv(vsta_SE, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/counts/vst_normalised_aware_counts_SE.csv", quote = FALSE)
vsta_SE <- read.csv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/counts/vst_normalised_aware_counts_SE.csv", header = TRUE, row.names = 1)
# keep only genes with expresion > 0
vsta_SE <- vsta_SE[rowSums(vsta_SE) > 0, ]

#' ZE  
# vsda_ZE <- varianceStabilizingTransformation(dds_ZE, blind = FALSE)
# vsta_ZE <- assay(vsda_ZE)
# vsta_ZE <- vsta_ZE - min(vsta_ZE)
# write.csv(vsta_ZE, file=here("analysis/salmon/vst_normalised_aware_counts_ZEall.csv"), quote = FALSE)
vsta_ZE <- read.csv(here("analysis/salmon/vst_normalised_aware_counts_ZEall.csv"), header = TRUE, row.names = 1)
# keep only genes with expresion > 0
vsta_ZE <- vsta_ZE[rowSums(vsta_ZE) > 0, ]

#' Scale the data
vstaSE_scaled <- t(scale(t(vsta_SE)))
vstaZE_scaled <- t(scale(t(vsta_ZE)))

#' ## Expression of DEGs  
#' 
#' ### DEGs in ZE, but not in SE  
DEGs_ZE_notSE <- setdiff(union(allDEGs_S, allDEGs_ZE, allDEGs_FMG), allDEGs_SE)

# without scaling
# heatmap.2(vsta_ZE[DEGs_ZE_notSE, ], labCol = colData(dds_ZE)$Replicate, labRow = FALSE, trace = "none", main = "ZE vs SE")
# scaled data
# heatmap.2(vsta_ZE[DEGs_ZE_notSE, ], scale = "row", labCol = colData(dds_ZE)$Replicate, labRow = FALSE, trace = "none", main = "ZE vs SE, scaled")

pdf(file = here("analysis/figures/DEGs_ZE_notSE_Replicate.pdf"), width = 15, height = 10)
#png(file = here("analysis/figures/DEGs_ZE_notSE.png"), width = 1000, height = 800, pointsize = 23)
heatmap.2(vstaZE_scaled[DEGs_ZE_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs unique to zygotic embryogenesis", 
          xlab = "time points", 
          ylab = paste(length(DEGs_ZE_notSE), "genes"))
dev.off()

#' How are these genes expressed in SE?
pdf(file = here("analysis/figures/DEGs_ZE_notSE_Replicate_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_ZE_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs unique to zygotic embryogenesis", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% DEGs_ZE_notSE), "genes"))
dev.off() 

#' Functional enrichment  
enrDEGs_ZE_notSE <- gopher(sub(".1$", "", DEGs_ZE_notSE),
                           task = list("go","mapman","pfam"), 
                           background = sub(".1$", "", NameExpGeneSEZE), 
                           url="pabies")
save(enrDEGs_ZE_notSE, file = here("analysis/gopher/enrDEGs_ZE_notSE.rda"))

#' ### DEGs in ZEm (B4-B9), but not in SE
DEGs_ZEm_notSE <- setdiff(allDEGs_ZE, allDEGs_SE)

pdf(file = here("analysis/figures/DEGs_ZEm_notSE_Replicate.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_ZEm_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs in zygotic embryo, but not in SE", 
          xlab = "time points", 
          ylab = paste(length(DEGs_ZEm_notSE), "genes"))
dev.off()

#' How are these genes expressed in SE?
pdf(file = here("analysis/figures/DEGs_ZEm_notSE_Replicate_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_ZEm_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs in zygotic embryo, but not in SE", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% DEGs_ZEm_notSE), "genes"))
dev.off() 

#' Functional enrichment  
enrDEGs_ZEm_notSE <- gopher(sub(".1$", "", DEGs_ZEm_notSE),
                           task = list("go","mapman","pfam"), 
                           background = sub(".1$", "", NameExpGeneSEZE), 
                           url="pabies")
save(enrDEGs_ZEm_notSE, file = here("analysis/gopher/enrDEGs_ZEm_notSE.rda"))

#' ### DEGs in FMG (B4-B9), but not in SE
DEGs_FMG_notSE <- setdiff(allDEGs_FMG, allDEGs_SE)

pdf(file = here("analysis/figures/DEGs_FMG_notSE_Replicate.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_FMG_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs in zygotic embryo, but not in SE", 
          xlab = "time points", 
          ylab = paste(length(DEGs_FMG_notSE), "genes"))
dev.off()

#' How are these genes expressed in SE?
pdf(file = here("analysis/figures/DEGs_FMG_notSE_Replicate_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_FMG_notSE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs in zygotic embryo, but not in SE", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% DEGs_FMG_notSE), "genes"))
dev.off() 

#' Functional enrichment  
enrDEGs_FMG_notSE <- gopher(sub(".1$", "", DEGs_FMG_notSE),
                            task = list("go","mapman","pfam"), 
                            background = sub(".1$", "", NameExpGeneSEZE), 
                            url="pabies")
save(enrDEGs_FMG_notSE, file = here("analysis/gopher/enrDEGs_FMG_notSE.rda"))

#' ### DEGs in SE, but not in ZE
DEGs_SE_notZE <- setdiff(allDEGs_SE, union(allDEGs_S, allDEGs_ZE, allDEGs_FMG))

pdf(file = here("analysis/figures/DEGs_SE_notZE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[DEGs_SE_notZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs unique to somatic embryogenesis", 
          xlab = "time points",
          ylab = paste(length(DEGs_SE_notZE), "genes"))
dev.off()

#' How are these genes expressed in ZE?
pdf(file = here("analysis/figures/DEGs_SE_notZE_inZE.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[rownames(vstaZE_scaled) %in% DEGs_SE_notZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs unique to somatic embryogenesis", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaZE_scaled) %in% DEGs_SE_notZE), "genes"))
dev.off() 

#' Functional enrichment  
enrDEGs_SE_notZE <- gopher(sub(".1$", "", DEGs_SE_notZE),
                          task = list("go","mapman","pfam"), 
                          background = sub(".1$", "", NameExpGeneSEZE), 
                          url="pabies")
save(enrDEGs_SE_notZE, file = here("analysis/gopher/enrDEGs_SE_notZE.rda"))

#' ### DEGs in SE and ZE
DEGs_SEandZE <- intersect(allDEGs_SE, union(allDEGs_S, allDEGs_ZE, allDEGs_FMG))

#' Plot expression of the same genes in SE and ZE dataset.  
# check if genes are in the same order in both datasets
all(rownames(vstaZE_scaled[DEGs_SEandZE, ]) == rownames(vstaSE_scaled[DEGs_SEandZE, ]))

# make two plots and compare them side by side  
pdf(file = here("analysis/figures/DEGs_SEandZE_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[DEGs_SEandZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages,
          labRow = NA, 
          main = "DEGs common to SE and ZE (expression in SE)", 
          xlab = "time points",
          ylab = paste(length(DEGs_SEandZE), "genes"))
dev.off()

pdf(file = here("analysis/figures/DEGs_SEandZE_inZE.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_SEandZE, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "DEGs common to SE and ZE (expresion in ZE)", 
          xlab = "time points", 
          ylab = paste(length(DEGs_SEandZE), "genes"))
dev.off()

# do we want to fix order of the samples in this case to compare expression profile of the genes by the timeline in two datasets?  

#' Functional enrichment  
enrDEGs_SEandZE <- gopher(sub(".1$", "", DEGs_SEandZE),
                           task = list("go","mapman","pfam"), 
                           background = sub(".1$", "", NameExpGeneSEZE), 
                           url="pabies")
save(enrDEGs_SEandZE, file = here("analysis/gopher/enrDEGs_SEandZE.rda"))


#' ### DEGs in ZE compared to FMG at each time point - all  
#' make union of DEGs expressed at each time point
all_DEGs_ZEFMG <- unique(unlist(lapply(genes_ZEFMG_all, rownames), use.names = FALSE))

pdf(file = here("analysis/figures/DEGs_ZEvsFMG.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[all_DEGs_ZEFMG, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "Genes DE in ZE compared to FMG at any time point", 
          xlab = "time points", 
          ylab = paste(length(all_DEGs_ZEFMG), "genes"))
dev.off()

#' How are these genes expressed in SE?
pdf(file = here("analysis/figures/DEGs_ZEvsFMG_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% all_DEGs_ZEFMG, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "Genes DE in ZE compared to FMG at any time point", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% all_DEGs_ZEFMG), "genes"))
dev.off() 

#' Functional enrichment  
enrDEGs_allZEFMG <- gopher(sub(".1$", "", all_DEGs_ZEFMG),
                        task = list("go","mapman","pfam"), 
                        background = sub(".1$", "", NameExpGeneZE), 
                        url="pabies")
save(enrDEGs_allZEFMG, file = here("analysis/gopher/enrDEGs_allZEFMG.rda"))


############### change
#' ### DEGs in ZE compared to FMG at each time point - up & down
#' make union of up-regulated genes
DEGs_ZEFMG_up <- unique(unlist(lapply(genes_ZEFMG_up, rownames), use.names = FALSE))

pdf(file = here("analysis/figures/DEGs_ZEvsFMG_up.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_ZEFMG_up, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "Genes up-regulated in ZE compared to FMG at any time point", 
          xlab = "time points", 
          ylab = paste(length(DEGs_ZEFMG_up), "genes"))
dev.off()

#' Functional enrichment  
enrDEGs_ZEFMG_up <- gopher(sub(".1$", "", DEGs_ZEFMG_up),
                           task = list("go","mapman","pfam"), 
                           background = sub(".1$", "", NameExpGeneZE), 
                           url="pabies")
save(enrDEGs_ZEFMG_up, file = here("analysis/gopher/enrDEGs_ZEFMG_up.rda"))
##################  



#' ### Genes DE in B7 compared to B6  
# make a union of genes DEGs between B7 and B6 in ZE or FMG
DEGs_B7vsB6 <- unique(c(rownames(genes_FMG_all$B7.FMG_B6.FMG), rownames(genes_ZE_all$B7.ZE_B6.ZE)))
#' How many DEGs do ZE and FMG have in common in this time point comparison?
dev.new()
venn::venn(list(FMG = rownames(genes_FMG_all$B7.FMG_B6.FMG), ZE = rownames(genes_ZE_all$B7.ZE_B6.ZE)), zcolor = "style", box = FALSE)

# divide genes in three groups  
B7vB6_FMGuniq <- setdiff(rownames(genes_FMG_all$B7.FMG_B6.FMG), rownames(genes_ZE_all$B7.ZE_B6.ZE))
B7vB6_ZEuniq <- setdiff(rownames(genes_ZE_all$B7.ZE_B6.ZE), rownames(genes_FMG_all$B7.FMG_B6.FMG))
B7vB6_FMGZEint <- intersect(rownames(genes_ZE_all$B7.ZE_B6.ZE), rownames(genes_FMG_all$B7.FMG_B6.FMG))

# prepare vector with names of the groups, reorder it
B7vB6_groupByTissue <- rep(c("FMG", "ZE", "FMG&ZE"), c(length(B7vB6_FMGuniq), length(B7vB6_ZEuniq), length(B7vB6_FMGZEint)))
names(B7vB6_groupByTissue) <- c(B7vB6_FMGuniq, B7vB6_ZEuniq, B7vB6_FMGZEint)
B7vB6_groupByTissue <- B7vB6_groupByTissue[DEGs_B7vsB6]

pdf(file = here("analysis/figures/DEGs_B7vsB6.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_B7vsB6, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))), # Pearson correlation; pearson.dist does not work
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          RowSideColors = c("blue", "green", "yellow")[as.integer(factor(B7vB6_groupByTissue))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate, 
          labRow = NA, 
          main = "Genes DE in B7 compared to B6 in ZE or FMG", 
          xlab = "time points", 
          ylab = paste(length(DEGs_B7vsB6), "genes"))
# legend("bottomleft",      
#        legend = levels(factor(B7vB6_groupByTissue)),
#        col = c("blue", "green", "yellow")[1:length(levels(factor(B7vB6_groupByTissue)))], 
#        fill = c("blue", "green", "yellow")[1:length(levels(factor(B7vB6_groupByTissue)))])
dev.off()

#' How are these genes expressed in SE?  
#' How many of these genes is expressed in SE?  
sum(DEGs_B7vsB6 %in% NameExpGeneSE)
# In % of all genes DE between B7 and B6:
sum(DEGs_B7vsB6 %in% NameExpGeneSE)/length(DEGs_B7vsB6)*100
#' How many of these genes is DE in SE? 
sum(DEGs_B7vsB6 %in% allDEGs_SE) 
# In % of all genes DE between B7 and B6:
sum(DEGs_B7vsB6 %in% allDEGs_SE)/length(DEGs_B7vsB6)*100

#' In which stages are these genes DE in SE?
DEGs_B7vsB6_DEinSE <- lapply(ddsDEGs_SE, function(x){
  x[DEGs_B7vsB6[DEGs_B7vsB6 %in% rownames(x)], ]
  })
elementNROWS(DEGs_B7vsB6_DEinSE)

#' When are these genes also DE in ZE?
DEGs_B7vsB6_DEinZE <- lapply(genes_ZE_all, function(x){
  x[DEGs_B7vsB6[DEGs_B7vsB6 %in% rownames(x)], ]
})
elementNROWS(DEGs_B7vsB6_DEinZE)

DEGs_B7vsB6_DEinFMG <- lapply(genes_FMG_all, function(x){
  x[DEGs_B7vsB6[DEGs_B7vsB6 %in% rownames(x)], ]
})
elementNROWS(DEGs_B7vsB6_DEinFMG)

DEGs_B7vsB6_DEinS <- lapply(genes_S_all, function(x){
  x[DEGs_B7vsB6[DEGs_B7vsB6 %in% rownames(x)], ]
})
elementNROWS(DEGs_B7vsB6_DEinS)

# plot
par(mfrow = c(3,1))
barplot(elementNROWS(DEGs_B7vsB6_DEinSE), ylim = c(0,4000), main = "Somatic embryo")
barplot(c(elementNROWS(DEGs_B7vsB6_DEinS), elementNROWS(DEGs_B7vsB6_DEinZE)), ylim = c(0,4000), main = "Seed & Zygotic embryo")
barplot(c(elementNROWS(DEGs_B7vsB6_DEinS), elementNROWS(DEGs_B7vsB6_DEinFMG)), ylim = c(0,4000), main = "Seed & FMG")

# What about vst?
sum(DEGs_B7vsB6 %in% rownames(vstaSE_scaled))
# In % of all genes DE between B7 and B6:
sum(DEGs_B7vsB6 %in% rownames(vstaSE_scaled))/length(DEGs_B7vsB6)*100

#' Plot expression of genes DE between B7 and B6 in SE - not the same order of genes!
B7vsB6_inSE <- vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_B7vsB6, ]
B7vB6_groupByTissue_inSE <- B7vB6_groupByTissue[rownames(B7vsB6_inSE)]

pdf(file = here("analysis/figures/DEGs_B7vsB6_inSE.pdf"), width = 15, height = 10)
heatmap.2(B7vsB6_inSE,
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          RowSideColors = c("blue", "green", "yellow")[as.integer(factor(B7vB6_groupByTissue[rownames(B7vsB6_inSE)]))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages,
          labRow = NA,
          main = "Genes DE in B7 compared to B6 in ZE or FMG (expression in SE)",
          xlab = "time points",
          ylab = paste(sum(DEGs_B7vsB6 %in% rownames(vstaSE_scaled)), "genes"))
dev.off()

#' Functional enrichment  
enrDEGs_B7vB6 <- gopher(sub(".1$", "", DEGs_B7vsB6),
                        task = list("go","mapman","pfam"), 
                        background = sub(".1$", "", NameExpGeneZE), 
                        url="pabies")
save(enrDEGs_B7vB6, file = here("analysis/gopher/enrDEGs_B7vB6.rda"))

#' ### DEGs unique to each tissue  
#' UpSetR does not have yet a way to collect the intersection results, so we extract it ourselves
upset_AllGeneral <- upset(fromList(DEGsListAllGeneral), sets = names(DEGsListAllGeneral), keep.order = TRUE, sets.bar.color = "grey", nsets = length(DEGsListAllGeneral), nintersects = NA)

DEGsListAllGeneral_unlisted <- unlist(DEGsListAllGeneral, use.names = FALSE)
DEGsListAllGeneral_unlisted <- DEGsListAllGeneral_unlisted[ !duplicated(DEGsListAllGeneral_unlisted) ]
# Now we know what does the rownumbers from "New_data" refer to in our list.  

#' SE - this group we have already analysed above  
#' 
#' FMG  
DEGs_FMG_unique <- DEGsListAllGeneral_unlisted[upset_AllGeneral$New_data$FMG == 1 & 
                                               upset_AllGeneral$New_data$seed == 0 & 
                                               upset_AllGeneral$New_data$zygotic_embryo == 0 & 
                                               upset_AllGeneral$New_data$somatic_embryo == 0]
#' ZE  
DEGs_ZEm_unique <- DEGsListAllGeneral_unlisted[upset_AllGeneral$New_data$FMG == 0 & 
                                                 upset_AllGeneral$New_data$seed == 0 & 
                                                 upset_AllGeneral$New_data$zygotic_embryo == 1 & 
                                                 upset_AllGeneral$New_data$somatic_embryo == 0]
#' S  
DEGs_S_unique <- DEGsListAllGeneral_unlisted[upset_AllGeneral$New_data$FMG == 0 & 
                                                 upset_AllGeneral$New_data$seed == 1 & 
                                                 upset_AllGeneral$New_data$zygotic_embryo == 0 & 
                                                 upset_AllGeneral$New_data$somatic_embryo == 0]
#' Expression profiles  
#' 
#' ZEm
pdf(file = here("analysis/figures/DEGs_ZEm_unique.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_ZEm_unique, ],
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate,
          labRow = NA,
          main = "DEGs unique to zygotic embryo",
          xlab = "time points",
          ylab = paste(length(DEGs_ZEm_unique), "genes"))
dev.off()

# expression of ZEm unique DEGs in SE  
pdf(file = here("analysis/figures/DEGs_ZEm_unique_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_ZEm_unique, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs unique to zygotic embryo (in SE)", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% DEGs_ZEm_unique), "genes"))
dev.off() 

#' FMG
pdf(file = here("analysis/figures/DEGs_FMG_unique.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_FMG_unique, ],
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate,
          labRow = NA,
          main = "DEGs unique to FMG",
          xlab = "time points",
          ylab = paste(length(DEGs_FMG_unique), "genes"))
dev.off()

# expression of FMG unique DEGs in SE  
pdf(file = here("analysis/figures/DEGs_FMG_unique_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_FMG_unique, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs unique to FMG (in SE)", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% DEGs_FMG_unique), "genes"))
dev.off() 

#' S
pdf(file = here("analysis/figures/DEGs_S_unique.pdf"), width = 15, height = 10)
heatmap.2(vstaZE_scaled[DEGs_S_unique, ],
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_ZE)$Time))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_ZE)$Replicate,
          labRow = NA,
          main = "DEGs unique to the whole seed (B1-B3)",
          xlab = "time points",
          ylab = paste(length(DEGs_S_unique), "genes"))
dev.off()

# expression of S unique DEGs in SE  
pdf(file = here("analysis/figures/DEGs_S_unique_inSE.pdf"), width = 15, height = 10)
heatmap.2(vstaSE_scaled[rownames(vstaSE_scaled) %in% DEGs_S_unique, ], 
          trace = "none",
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(X){hclust(X,method="ward.D")},
          ColSideColors = pal[as.integer(factor(colData(dds_SE)$Stages))],
          col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(100)),
          labCol = colData(dds_SE)$Stages, 
          labRow = NA, 
          main = "DEGs unique to the whole seed, B1-B3 (in SE)", 
          xlab = "time points", 
          ylab = paste(sum(rownames(vstaSE_scaled) %in% DEGs_S_unique), "genes"))
dev.off() 

#' #### Number of unique DEGs in each stage  
#' ZEm  
DEGs_ZEm_unique_perStage <- lapply(genes_ZE_all, function(x){
  x[DEGs_ZEm_unique[DEGs_ZEm_unique %in% rownames(x)], ]
})
elementNROWS(DEGs_ZEm_unique_perStage)
#' FMG  
DEGs_FMG_unique_perStage <- lapply(genes_FMG_all, function(x){
  x[DEGs_FMG_unique[DEGs_FMG_unique %in% rownames(x)], ]
})
elementNROWS(DEGs_FMG_unique_perStage)
#' S  
DEGs_S_unique_perStage <- lapply(genes_S_all, function(x){
  x[DEGs_S_unique[DEGs_S_unique %in% rownames(x)], ]
})
elementNROWS(DEGs_S_unique_perStage)
#' SE  
DEGs_SE_unique_perStage <- lapply(ddsDEGs_SE, function(x){
  x[DEGs_SE_notZE[DEGs_SE_notZE %in% rownames(x)], ]
})
elementNROWS(DEGs_SE_unique_perStage)

# plot
par(mfrow = c(3,1))
barplot(elementNROWS(DEGs_SE_unique_perStage), ylim = c(0,12000), main = "Somatic embryo", ylab = " number of genes")
barplot(c(elementNROWS(DEGs_S_unique_perStage), elementNROWS(DEGs_ZEm_unique_perStage)), ylim = c(0,1000), main = "Seed & Zygotic embryo", ylab = " number of genes", col = rep(c("black", "grey"), c(2,6)))
barplot(c(elementNROWS(DEGs_S_unique_perStage), elementNROWS(DEGs_FMG_unique_perStage)), ylim = c(0,1000), main = "Seed & FMG", ylab = " number of genes", col = rep(c("black", "grey"), c(2,6)))


#' Functional enrichment  
enrTissueUnique <- lapply(list(DEGs_FMG_unique, DEGs_ZEm_unique, DEGs_S_unique), function(x){
  gopher(sub(".1$", "", x),
         task = list("go","mapman","pfam"), 
         background = sub(".1$", "", NameExpGeneZE), 
         url="pabies")
})
names(enrTissueUnique) <- c("enrDEGs_FMG_unique", "enrDEGs_ZEm_unique", "enrDEGs_S_unique")
save(enrTissueUnique, file = here("analysis/gopher/enrDEGs_TissueUnique.rda"))

#' ### Functional enrichment of DEGs in ZE  
#' ZEm  
enrDEGs_ZEm_allStages <- lapply(lapply(genes_ZE_all, rownames), function(x){
  gopher(sub(".1$", "", x),
         task = list("go","mapman","pfam"), 
         background = sub(".1$", "", NameExpGeneZE), 
         url="pabies")
})
save(enrDEGs_ZEm_allStages, file = here("analysis/gopher/enrDEGs_ZEm_allStages.rda"))

#' FG  
enrDEGs_FG_allStages <- lapply(lapply(genes_FMG_all, rownames), function(x){
  gopher(sub(".1$", "", x),
         task = list("go","mapman","pfam"), 
         background = sub(".1$", "", NameExpGeneZE), 
         url="pabies")
})
save(enrDEGs_FG_allStages, file = here("analysis/gopher/enrDEGs_FG_allStages.rda"))

#' S
enrDEGs_S_allStages <- lapply(lapply(genes_S_all, rownames), function(x){
  gopher(sub(".1$", "", x),
         task = list("go","mapman","pfam"), 
         background = sub(".1$", "", NameExpGeneZE), 
         url="pabies")
})
save(enrDEGs_S_allStages, file = here("analysis/gopher/enrDEGs_S_allStages.rda"))

#' DEGs in ZEm vs FMG in each stage
enrDEGs_ZEFMG_allStages <- lapply(lapply(genes_ZEFMG_all, rownames), function(x){
  gopher(sub(".1$", "", x),
         task = list("go","mapman","pfam"), 
         background = sub(".1$", "", NameExpGeneZE), 
         url="pabies")
})
save(enrDEGs_ZEFMG_allStages, file = here("analysis/gopher/enrDEGs_ZEFMG_allStages.rda"))

#--------------------------------------------Treemaps-------------------------------------------------------  


#' ### Plot treemaps  
#' 
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

# plot at once DEGs from different groups
treemapList <- list(enrDEGs_allZEFMG, enrDEGs_B7vB6, enrDEGs_SEandZE, enrDEGs_ZE_notSE, enrDEGs_ZEm_notSE, enrDEGs_FMG_notSE, enrDEGs_SE_notZE)
names(treemapList) <- c("enrDEGs_allZEFMG", "enrDEGs_B7vB6", "enrDEGs_SEandZE", "enrDEGs_ZE_notSE", "enrDEGs_ZEm_notSE", "enrDEGs_FMG_notSE", "enrDEGs_SE_notZE")
plotTreemapList(treemapList, "analysis/gopher/")

# DEGs unique to each tissue
plotTreemapList(enrTissueUnique, "analysis/gopher/")

# DEGs in ZE in all the stages ### continue
plotTreemapList(enrDEGs_ZEm_allStages, "analysis/gopher/")
plotTreemapList(enrDEGs_FG_allStages, "analysis/gopher/")
plotTreemapList(enrDEGs_S_allStages, "analysis/gopher/")
plotTreemapList(enrDEGs_ZEFMG_allStages, "analysis/gopher/")

#' ### Save enrichment as tsv files  
enrList <- dir(here("analysis/gopher"), ".rda", full.names = T)
lapply(enrList, load, .GlobalEnv) 

# New list with objects, not names ## add new objects, if they will be created
enrListObj <- list(enrDEGs_allZEFMG, enrDEGs_B7vB6, enrDEGs_FMG_notSE, enrDEGs_SEandZE, enrDEGs_SEandZE, enrDEGs_SE_notZE, enrDEGs_ZEm_notSE, enrDEGs_ZE_notSE)
names(enrListObj) <- c("enrDEGs_allZEFMG", "enrDEGs_B7vB6", "enrDEGs_FMG_notSE", "enrDEGs_SEandZE", "enrDEGs_SEandZE", "enrDEGs_SE_notZE", "enrDEGs_ZEm_notSE", "enrDEGs_ZE_notSE")
enr2tsv(enrListObj, filePrefix = here("analysis/gopher/Intersection"))

# for DEGs unique in each tissue
enr2tsv(enrTissueUnique, filePrefix = here("analysis/gopher/Intersection"))

# for DEGs in ZE in all the stages
enr2tsv(enrDEGs_ZEm_allStages, filePrefix = here("analysis/gopher/DEGs"))
enr2tsv(enrDEGs_FG_allStages, filePrefix = here("analysis/gopher/DEGs"))
enr2tsv(enrDEGs_S_allStages, filePrefix = here("analysis/gopher/DEGs"))
enr2tsv(enrDEGs_ZEFMG_allStages, filePrefix = here("analysis/gopher/DEGs"))


#---------------------------------------------% DEGs in expressed genes---------------------------------------  


#' # Percentage of DEGs in expressed genes  
#' ## SE  
#' Number of DEGs in the experiment  
length(allDEGs_SE)
#' How many percent of expressed genes do DEGs represent?
length(allDEGs_SE)/NrExpGeneSE*100
#' Number of DEGs in each stage
sapply(ddsDEGs_SE, nrow)
#' How many percent of all the expressed genes do DEGs represent in each stage?
sapply(ddsDEGs_SE, nrow)/NrExpGeneSE*100
#' How many percent of genes expressed in each stage do DEGs in each stage represent?
sapply(ddsDEGs_SE, nrow)/NrExpGeneSEstage[-8]*100

barplot(sapply(ddsDEGs_SE, nrow)/NrExpGeneSE*100, ylim = c(0,50))

#' ## ZE  
#' ### seed  
#' Number of DEGs in the tissue  
length(allDEGs_S)
#' How many percent of expressed genes do DEGs represent?
length(allDEGs_S)/NrExpGeneZE*100
#' Number of DEGs in each stage
sapply(genes_S_all, nrow)
#' How many percent of all the expressed genes do DEGs represent in each stage?
sapply(genes_S_all, nrow)/NrExpGeneZE*100
#' How many percent of genes expressed in each stage do DEGs in each stage represent?
sapply(genes_S_all, nrow)/NrExpGeneZEstage[grepl("B2.S", names(NrExpGeneZEstage))]*100

#' ### zygotic embryo  
#' Number of DEGs in the tissue  
length(allDEGs_ZE)
#' How many percent of expressed genes do DEGs represent?
length(allDEGs_ZE)/NrExpGeneZE*100
#' Number of DEGs in each stage
sapply(genes_ZE_all, nrow)
#' How many percent of all the expressed genes do DEGs represent in each stage?
sapply(genes_ZE_all, nrow)/NrExpGeneZE*100
#' How many percent of genes expressed in each stage do DEGs in each stage represent?
sapply(genes_ZE_all, nrow)/NrExpGeneZEstage[grepl("ZE", names(NrExpGeneZEstage))]*100

#' ### FMG  
#' Number of DEGs in the tissue  
length(allDEGs_FMG)
#' How many percent of expressed genes do DEGs represent?
length(allDEGs_FMG)/NrExpGeneZE*100
#' Number of DEGs in each stage
sapply(genes_FMG_all, nrow)
#' How many percent of all the expressed genes do DEGs represent in each stage?
sapply(genes_FMG_all, nrow)/NrExpGeneZE*100
#' How many percent of genes expressed in each stage do DEGs in each stage represent?
sapply(genes_FMG_all, nrow)/NrExpGeneZEstage[grepl("FMG", names(NrExpGeneZEstage))]*100

#' ### all tissues together  
# before merging ZE and FMG, elements of the lists have to have the same names  
geneNames_FMG_all <- lapply(genes_FMG_all, rownames)
names(geneNames_FMG_all) <- gsub(".FMG", "", names(geneNames_FMG_all))
geneNames_ZE_all <- lapply(genes_ZE_all, rownames)
names(geneNames_FMG_all) <- gsub(".ZE", "", names(geneNames_FMG_all))
# merge two lists
geneNames_mergedZEandFMG <- Map(c, geneNames_FMG_all, geneNames_ZE_all)
geneNames_mergedZEandFMG <- lapply(geneNames_mergedZEandFMG, unique)
#' Number of DEGs in each stage (in at least one of the tissues: zygotic embryo or FMG)
NrDEGs_mergedZEandFMG <- sapply(geneNames_mergedZEandFMG, length)
NrDEGs_mergedZEandFMG

#' Number of DEGs in the experiment
NrDEGs_ZEtotal <- length(unique(c(allDEGs_S, allDEGs_FMG, allDEGs_ZE)))
NrDEGs_ZEtotal
#' How many percent of expressed genes do DEGs represent?
NrDEGs_ZEtotal/NrExpGeneZE*100
#' Number of DEGs in each stage
c(sapply(genes_S_all, nrow),NrDEGs_mergedZEandFMG)
#' How many percent of all the expressed genes do DEGs represent in each stage?
c(sapply(genes_S_all, nrow),NrDEGs_mergedZEandFMG)/NrExpGeneZE*100
#' #' How many percent of genes expressed in each stage do DEGs in each stage represent?
#' sapply(genes_FMG_all, nrow)/NrExpGeneZEstage[grepl("FMG", names(NrExpGeneZEstage))]*100
barplot(c(sapply(genes_S_all, nrow),NrDEGs_mergedZEandFMG)/NrExpGeneZE*100, ylim = c(0,20))

#' #' 
#' #' ```{r empty,eval=FALSE,echo=FALSE}
#' #' ```
#' #' 
#' #' # Session Info
#' #' ```{r session info, echo=FALSE}
#' #' sessionInfo()
#' #' ```
#' 
