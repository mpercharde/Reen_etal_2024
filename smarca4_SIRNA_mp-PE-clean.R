

#### SMARCA4 script - Aug 2024 -
# Michelle Percharde, MRC LMS #

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# version.string R version 4.2.3 (2023-03-15)

library(tidyverse)
library(dplyr)
library(furrr)
library(limma)
library(gdata)
library(ggplot2)
library(gplots)
library("Rsamtools")
library("DESeq2")
library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(pheatmap)
library(heatmap3)
library(here)
library(glue)

setwd("/Volumes/chrom_dev$/chrom_dev/+ANALYSIS/Michelle/367_250819_4OHT_siSMARCA4/")
here::i_am("scripts/smarca4_SIRNA_mp-PE.R")
plot_dir <- "/Volumes/chrom_dev$/chrom_dev/+ANALYSIS/Michelle/367_250819_4OHT_siSMARCA4/plots/"
file_dir <- "/Volumes/chrom_dev$/chrom_dev/+ANALYSIS/Michelle/367_250819_4OHT_siSMARCA4/counts/"
count_dir <- "/Volumes/chrom_dev$/chrom_dev/+ANALYSIS/Michelle/367_250819_4OHT_siSMARCA4/telocal/"

# ## read data using loop - first do the 1st replicate manually since full join can't take a null object
counts <- read.table((paste0(count_dir,"4OHT_siNT_2_1VM.cntTable.gz")), header=T, sep="\t", na.strings="")[1:62754,] #ENSG292373 - rest of numbers is TEs

#list for the rest of the files
files <- list.files(count_dir)[2:29] #note that the length is 2:max number of files - change each time
#for loop to collapse and then merge
for (f in files) {
  dat <- read.table((paste0(count_dir,f)), header=T, sep="\t", na.strings=c("", " "))[1:62754,]
  counts <- full_join(counts, dat, by="gene.TE", relationship="one-to-one")
}

#output the very basic raw file
   write.table(counts, paste0(file_dir, today(), "_counts_genes_raw_unlabeled.txt"), sep="\t")
#
# # read gtf used for counting in as GRanges object
library("rtracklayer")
GRCh38_112_gtf <- import("../362_240801_4OHT_timecourse/ref/Homo_sapiens.GRCh38.112.gtf")
  GRCh38_112_df <- as(GRCh38_112_gtf, "data.frame") %>% # make it a dataframe
    filter(type == "gene") # get only the genes (i.e not all the isoforms, exons etc)
  rm(GRCh38_112_gtf) # tidy to save space

# make version of results with gene names - mp manual
  #----------------------------
  counts2 <- counts %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")

colnames(counts)[1] <- "gene_id"
counts2 <- left_join(counts,
                          select(GRCh38_112_df,gene_name, gene_id),
                          by = "gene_id")

##load metadata and change the colnames
metadata <-read.table("metadata2.txt", header=T)
colnames(counts2)[2:30] <-metadata$SampleID
  write.table(counts2, paste0(file_dir, today(), "_gene_counts_raw_renamed_labeled.txt"), sep="\t")

#use GY code -import in both genes and TE loci together in one file
  counts_full =
    metadata %>%
    dplyr::pull(OrigName) %>%
    furrr::future_map_dfc(function(f){
      readr::read_tsv(
        here("telocal", glue("{f}.cntTable.gz")),
        col_names=c("feature", f), col_types="ci", skip=1) %>%
        tibble::column_to_rownames("feature")
    })
  rownames(counts_full) = str_remove(rownames(counts_full), ":.*")

colnames(counts_full)[1:29] <-metadata$SampleID

  write.table(counts_full, paste0(file_dir, today(), "_gene+te_counts_raw_labeled.txt"), sep="\t")
# nrow(counts_full)
  # [1] 4810996

#george_TE loci
counts_teloci =
    DESeq2::DESeqDataSetFromMatrix(
      counts_full[!stringr::str_starts(rownames(counts_full), "ENS"), ],
      metadata,
      ~0)
counts_teloci = counts_teloci[rowSums(counts(counts_teloci)>0) >= 6, ]
counts_teloci = DESeq2::estimateSizeFactors(counts_teloci)
nrow(assay(counts_teloci)) #295198 - this is reduced from many more - 4.8 million (above) - provides more power

head(assay(counts_teloci))
  write.table(assay(counts_teloci), paste0(file_dir, today(), "_teloci_counts_raw_filtered_renamed.txt"), sep="\t")

###
counts_teloci <- read.table("counts/2024-08-20_teloci_counts_raw_filtered_renamed.txt", header=T)
#find out the rough normalisation factor
colSums(counts_teloci[,1:29])
sums <- data.frame(counts = colSums(counts_teloci[,1:29]))
  write.table(sums, paste0(file_dir, today(), "_teloci_colsums.txt"), sep="\t")

### summarise the TEs into subfamilies - GY code

  counts_tefam =
    DESeq2::DESeqDataSetFromMatrix(
      counts_full %>%
        tibble::rownames_to_column(var="feature") %>%
        dplyr::filter(!stringr::str_starts(feature, "ENS")) %>%
        tidyr::separate_wider_delim(feature, names=c("family", "element", "chrom", "start", "stop"), delim="|") %>%
        dplyr::mutate(te_type=glue("{family}_{element}")) %>%
        dplyr::select(-c("family", "element", "chrom", "start", "stop")) %>%
        dplyr::group_by(te_type) %>%
        dplyr::summarise(across(everything(), sum)) %>%
        tibble::column_to_rownames(var="te_type"),
      metadata,
      ~0)
 nrow(assay(counts_tefam))
# [1] 1261
counts_tefam = counts_tefam[rowSums(assay(counts_tefam)) > 0, ]
 nrow(assay(counts_tefam))
# [1] 1257
counts_tefam = DESeq2::estimateSizeFactors(counts_tefam)
  write.table(assay(counts_tefam), paste0(file_dir, today(), "_tesubfam_counts_raw_renamed.txt"), sep="\t")

#################################################################################################################################################

#normalise the gene-counts table - use vst
counts_genes = counts2[rowSums((counts2[,2:30])>0) >= 5, ]
# write.table(counts_genes, paste0(file_dir, today(), "_genes_counts_raw_filtered_renamed.txt"), sep="\t")
#drop nas
counts_genes <- counts_genes %>% drop_na() #26417
# rownames(counts_genes) <- counts_genes$gene_name #still has duplicate rownames
counts_genes$gene_name <- make.unique(counts_genes$gene_name, sep = "_")
# write.table(counts_genes, paste0(file_dir, today(), "_genes_counts_raw_filtered_renamed_unique.txt"), sep="\t") #think it works
#reset rownames to the gene_name
rownames(counts_genes) <- counts_genes$gene_name
counts_genes <-counts_genes[,2:30]
  write.table(counts_genes, paste0(file_dir, today(), "_genes_counts_raw_filtered_renamed_unique.txt"), sep="\t") #

#put this into a proper deseq object
genes <- DESeqDataSetFromMatrix(counts_genes, metadata, ~0)
genes = DESeq2::estimateSizeFactors(genes) #done manually here but also done under the hood by deseq2()
vsd_g <- vst(g)
g <- assay(vsd_g)

  write.table(g, paste0(file_dir, today(), "_genes_counts_vst_normalised.txt"), sep="\t") #


#PCA to see clustering of samples
PCA1<-plotPCA(vsd_g, intgroup = "Condition", returnData=T) #another way to do a PCA
percentVar <- round(100 * attr(PCA1, "percentVar"))
  #

pdf(paste0(plot_dir, today(), "_",  "PCA_genes_all_byReplicate+Condition.pdf"))
ggplot(PCA1, aes(x = PC1, y = PC2)) +
    ggtitle("PCA plot: genes - all samples") +
    theme(plot.title = element_text(hjust = 0.5, size=10)) +
    geom_point(size=3, aes(color = metadata$Condition, shape = as.factor(metadata$Replicate))) +
    theme(text = element_text(size=10), axis.text = element_text(vjust=1,size=10)) +
    theme(legend.title=element_blank()) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

#do same for TEsubfams
TEs <- counts_tefam

vsd_t <- vst(TEs)
t <- assay(vsd_t)
  write.table(t, paste0(file_dir, today(), "_tesubfam_counts_vsd_normalised.txt"), sep="\t") #

PCA2<-plotPCA(vsd_t, intgroup = "Condition", returnData=T) #another way to do a PCA
percentVar <- round(100 * attr(PCA2, "percentVar"))

pdf(paste0(plot_dir, today(), "_",  "PCA_TEfam_all_bySample+Condition.pdf"))
ggplot(PCA2, aes(x = PC1, y = PC2)) +
    ggtitle("PCA plot: TE subfam all samples") +
    theme(plot.title = element_text(hjust = 0.5, size=10)) +
    geom_point(size=3, aes(color = metadata$SampleID, shape = metadata$Condition)) +
    theme(text = element_text(size=10), axis.text = element_text(vjust=1,size=10)) +
    theme(legend.title=element_blank()) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  dev.off()

##### reorder everything

metadata2 <-read.table("metadata2.txt", header=T)
g2<-g[,order(metadata2$Order)]
t2<-t[,order(metadata2$Order)]
metadata2<-metadata2[order(metadata2$Order),]


pcorr<-function(x) as.dist(1-cor(t(x),method="pearson")) #this is a function that allows pearson correlation to drive heatmap (Better than default euclidian distfun)
my_palette2 <- colorRampPalette(c("blue4", "white", "red4"))(n=299)

keep <-rainbow(length(unique(metadata2$Vee)))[as.numeric(as.factor(metadata2$Vee))]
con <- rainbow(length(unique(metadata2$Condition)))[as.numeric(as.factor(metadata2$Condition))]

##############################################################################################################################
###DE analysis Vee samples
##############################################################################################################################
  # t2<-t[,order(metadata2$Order)]
  # metadata2<-metadata2[order(metadata2$Order),]

metadata2 <-read.table("metadata2.txt", header=T)
subf <- read.table("counts/2024-08-20_tesubfam_counts_raw_renamed.txt", sep="\t")

metadata3 <- metadata2[(metadata2$Vee == "Y"),] #keep metadata rows only that passed seq criteria

subf2 <- subf %>%
  select(all_of(metadata3$SampleID)) #worked - now we have a filtered raw counts and metadata2 to work with

#reorder both of them now for ease
subf2 <- subf2[,order(metadata3$Order)]
metadata3 <- metadata3[order(metadata3$Order),]

TEs <- DESeqDataSetFromMatrix(countData = subf2,  colData = metadata3, design = ~ Condition) #DESeq object made
TEs <- DESeq2::estimateSizeFactors(TEs)

#normalise both
vsd_t <-vst(TEs)

#PCA with the TEs own size factors first

PCA3<-plotPCA(vsd_t, intgroup = "Condition", returnData=T)
percentVars <- round(100 * attr(PCA1, "percentVar"))

pdf(paste0(plot_dir, today(), "_",  "PCA_TEs_SEL_colourbyCondition_shapebyRep.pdf"))
ggplot(PCA3, aes(x = PC1, y = PC2)) +
  ggtitle("PCA plot: TEs") +
  theme(plot.title = element_text(hjust = 0.5, size=10)) +
  geom_point(size=3, aes(color = metadata3$Condition, shape = as.factor(metadata3$Replicate))) +
  theme(text = element_text(size=10), axis.text = element_text(vjust=1,size=10)) +
  theme(legend.title=element_blank()) +
  xlab(paste0("PC1: ",percentVars[1],"% variance")) +
  ylab(paste0("PC2: ",percentVars[2],"% variance"))
dev.off()


##############################################################################################################################
### toptable analysis
##############################################################################################################################

#compare all
DEmatrix <-DESeq(TEs) #this runs without errors
#
resSvsP <- results(DEmatrix, contrast=c("Condition","siNT_4OHT","siNT_DMSO"))
resSSvsSN <- results(DEmatrix, contrast=c("Condition","siSMARCA4_4OHT","siNT_4OHT"))
resSSvsP <- results(DEmatrix, contrast=c("Condition","siSMARCA4_4OHT","siNT_DMSO"))

par(mfrow=c(1,3))
plotMA(resSvsP, main="ctl sen vs prolif")
plotMA(resSSvsSN, main="sen siS vs sen siN")
plotMA(resSSvsP, main="sen siS vs sen prolif")

#shrink the log2fc for everything - use ashr
library(ashr)

resSvsPs <- lfcShrink(DEmatrix, contrast=c("Condition","siNT_4OHT","siNT_DMSO"),type='ashr')
resSSvsSNs <- lfcShrink(DEmatrix, contrast=c("Condition","siSMARCA4_4OHT","siNT_4OHT"),type='ashr')
resSSvsPs <- lfcShrink(DEmatrix, contrast=c("Condition","siSMARCA4_4OHT","siNT_DMSO"),type='ashr')


par(mfrow=c(2,3))
plotMA(resSvsP, main="ctl sen vs prolif")
plotMA(resSSvsSN, main="sen siS vs sen siN")
plotMA(resSSvsP, main="sen siS vs sen prolif")
plotMA(resSvsPs, main="ctl sen vs prolif, shrunk")
plotMA(resSSvsSNs, main="sen siS vs sen siN, shrunk")
plotMA(resSSvsPs, main="sen siS vs sen prolif, shrunk")
#

resSvsPs<-resSvsPs[order((resSvsPs$log2FoldChange)),] #order by increasing FC
    write.table(resSvsPs, paste0(file_dir, today(), "_TEsubf_SEL_toptable_SENvPROLIF.txt"), sep="\t")
resSSvsSNs<-resSSvsSNs[order((resSSvsSNs$log2FoldChange)),] #order by increasing FC
    write.table(resSSvsSNs, paste0(file_dir, today(), "_TEsubf_SEL_toptable_SEN-SISvSEN-SIN.txt"), sep="\t")
resSSvsPs<-resSSvsPs[order((resSSvsPs$log2FoldChange)),] #order by increasing FC
    write.table(resSSvsPs, paste0(file_dir, today(), "_TEsubf_SEL_toptable_SEN-SISvPROLIF.txt"), sep="\t")

# #################### VOLCANO PLOTS #################
#
# #adapted from https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/
#

resSvsPs <- read.table("counts/2024-08-20_TEsubf_SEL_toptable_SENvPROLIF.txt", header=T)
resSSvsSNs <- read.table("counts/2024-08-20_TEsubf_SEL_toptable_SEN-SISvSEN-SIN.txt", header=T)

BiocManager::install("ggrepel")
library(ggrepel)

volSP <- data.frame(resSvsPs[,1:5]) #make it a dataframe so it can be worked on
  volSP <- volSP %>%
    mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                 log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                 TRUE ~ "ns"))

cols <- c("up" = "goldenrod", "down" = "blue3", "ns" = "grey")
sizes <- c("up" = 2, "down" = 2, "ns" = 1)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

pdf(paste0(plot_dir, today(), "_",  "Volcano_SENvsProlif.pdf"))
ggplot(data = volSP,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  geom_point(aes(colour = gene_type),
             alpha = 0.8,
             shape = 16,
             size = 2) +
  geom_text_repel(aes(label = ifelse(gene_type %in% c("up", "down"), as.character(rownames(volSP)), "")),
                  nudge_y = 2,
                  color = "black",
                  max.overlaps = 40,
                  size = 2) +
  # geom_point(data = volSP["LINE/L1_L1HS",],
  #            shape = 21,
  #            size = 2,
  #            fill = "firebrick",
  #            colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  # geom_label_repel(data = volSP["LINE/L1_L1HS",],
  #                  aes(label = rownames(volSP["LINE/L1_L1HS",])),
  #                  force = 2,
  #                  nudge_x = -3) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-8, 8, 2)),
                     limits = c(-8, 8)) +
  labs(title = "TE subfamilies - SEN vs Ctl",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression Change")
dev.off()


volSS <- data.frame(resSSvsSNs[,1:5]) #make it a dataframe so it can be worked on
volSS <- volSS %>%
  mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                               log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns"))

pdf(paste0(plot_dir, today(), "_",  "Volcano_SEN-siSvsSEN-siN.pdf"))
ggplot(data = volSS,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  geom_point(aes(colour = gene_type),
             alpha = 0.8,
             shape = 16,
             size = 2) +
  geom_text_repel(aes(label = ifelse(gene_type %in% c("up", "down"), as.character(rownames(volSS)), "")),
                  nudge_y = 2,
                  color = "black",
                  max.overlaps = 40,
                  size = 2) +
  # geom_point(data = volSP["LINE/L1_L1HS",],
  #            shape = 21,
  #            size = 2,
  #            fill = "firebrick",
  #            colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  # geom_label_repel(data = volSP["LINE/L1_L1HS",],
  #                  aes(label = rownames(volSP["LINE/L1_L1HS",])),
  #                  force = 2,
  #                  nudge_x = -3) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-8, 8, 2)),
                     limits = c(-8, 8)) +
  labs(title = "TE subfamilies - SEN SiS vs SEN SiN",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression Change")
dev.off()

#do again but relax thresholds just to do fdr

volSS <- volSS %>%
  mutate(gene_type = case_when(log2FoldChange >= 0 & padj <= 0.05 ~ "up",
                               log2FoldChange <= 0 & padj <= 0.05 ~ "down",
                               TRUE ~ "ns"))

pdf(paste0(plot_dir, today(), "_",  "Volcano_FDRonly_top10_SEN-siSvsSEN-siN.pdf"))
ggplot(data = volSS,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  geom_point(aes(colour = gene_type),
             alpha = 0.8,
             shape = 16,
             size = 2) +
  # geom_text_repel(aes(label = ifelse(gene_type %in% c("up"), as.character(rownames(volSS)), "")),
  #                 nudge_y = 2,
  #                 color = "black",
  #                 max.overlaps = 40,
  #                 size = 2) +
  geom_point(data = volSS[1247:1257,],
             shape = 21,
             size = 2,
             fill = "firebrick",
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  # geom_vline(xintercept = c(log2(0.5), log2(2)),
  #            linetype = "dashed") +
  geom_label_repel(data = volSS[1247:1257,],
                   aes(label = rownames(volSS[1247:1257,])),
                   force = 2,
                   nudge_x = -3) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-8, 8, 2)),
                     limits = c(-6, 6)) +
  labs(title = "TE subfamilies - SEN SiS vs SEN SiN",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression Change")
dev.off()


#search for satellite

sat <-volSS[grep("Satellite*", rownames(volSS)),]

pdf(paste0(plot_dir, today(), "_",  "Volcano_SEN-siSvsSEN-siN-satellites.pdf"))
ggplot(data = volSS,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  geom_point(aes(colour = gene_type),
             alpha = 0.8,
             shape = 16,
             size = 2) +
  # geom_text_repel(aes(label = ifelse(gene_type == "up", as.character(rownames(volSP)), "")),
  #                 nudge_y = 2,
  #                 color = "black",
  #                 max.overlaps = 40,
  #                 size = 2) +
  geom_point(data = sat,
             shape = 21,
             size = 2,
             fill = "firebrick",
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  geom_label_repel(data = sat, # Add labels last to appear as the top layer
                   aes(label = as.character(rownames(sat))),
                   force = 2,
                   nudge_x = -2) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-8, 8, 2)),
                     limits = c(-6, 6)) +
  labs(title = "TE subfamilies - Satellites - SEN SiS vs SEN SiN",
       x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       colour = "Expression Change")
dev.off()


########################################################################################################################
