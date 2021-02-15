
# Load in libraries
library(GEOquery)
library(BiocManager)
library(dplyr)
library(tximport)
library(biomaRt)
library(vsn)
library(ggplot2)

library(DESeq2)
library(edgeR)
library(limma)


# Code adapted from: 
# DESeq2: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#removing-hidden-batch-effects

# Create file names to get abundance files from Kallisto 
dir <- "F:/Transcriptomics/kallisto/outputs/Cardiotoxicity"

# Check for all files in the directory to ensure they are present 
list.files(dir)

# Read in meta data table
metadata <- xlsx::read.xlsx(file = "data/RNA-seq/metadata.xlsx", header = TRUE, sheetIndex = 1)

# Read in each file based off the format, and label sample numbers
files <- file.path(dir,metadata$Folder.for.kallisto,"abundance.tsv")
names(files) <- paste0("sample",1:dim(metadata)[1])

# Use biomArt to gather transcript and gene IDs
# In this case, map to Entrez IDs
# Later can map back to symbols for plotting purposes
mart  <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl")
t2g6 <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name", "entrezgene_id", "transcript_version","description"), 
                       mart = mart)
t2g6 <- within(t2g6,target_id <- paste(ensembl_transcript_id, transcript_version,sep = "."))
t2g6 <- dplyr::rename(t2g6, transcript_id = ensembl_transcript_id, 
                      ens_gene = ensembl_gene_id, 
                      ext_gene = external_gene_name, 
                      etz_gene = entrezgene_id, 
                      descript = description)

# Use Entrez IDs for mapping to transcripts
tx2gene <- t2g6[,c(7,4)] %>% filter(!is.na(etz_gene)) %>% unique()

###### Use TxImport to summarize transcript changes to the gene level #####################################
# Notes on kallisto implementation
# Used the Ensemble v96 transcriptome
# default is not corrected for gene length or library depth
txi <- tximport::tximport(files,type = "kallisto", txOut = FALSE, tx2gene = tx2gene)

# Label all conditions in the matrix in order of appearance in the each of the files variables
condition = c("5FU_24","5FU_24","5FU_24","5FU_24","5FU_24","5FU_24",
              "5FU_6", "5FU_6", "5FU_6", "5FU_6", "5FU_6", "5FU_6",
              "Ace_24", "Ace_24", "Ace_24", "Ace_24", "Ace_24", "Ace_24", "Ace_24",
              "Ace_6","Ace_6","Ace_6","Ace_6","Ace_6","Ace_6",
              "DMSO1_24","DMSO1_24","DMSO1_24","DMSO1_24","DMSO1_24","DMSO1_24",
              "DMSO1_6","DMSO1_6","DMSO1_6","DMSO1_6","DMSO1_6","DMSO1_6",
              "DMSO2_24","DMSO2_24","DMSO2_24","DMSO2_24","DMSO2_24","DMSO2_24","DMSO2_24",
              "DMSO2_6","DMSO2_6","DMSO2_6","DMSO2_6","DMSO2_6","DMSO2_6","DMSO2_6",
              "Dox_24","Dox_24","Dox_24","Dox_24","Dox_24","Dox_24","Dox_24",
              "Dox_6","Dox_6","Dox_6","Dox_6","Dox_6","Dox_6","Dox_6"
)

# Change condition vector from character vector to factor vector
condition = factor(condition)

# Convert conditions to data frame for manipulation
sampleTable <- data.frame(condition, Drug = metadata$Drug, Timepoint = factor(metadata$Timepoint))

# Create row names based on the genes in the same order
rownames(sampleTable) <- colnames(txi$counts)

############################################ DESeq2 for DEG analysis ###################################################
# Generate a DESeq2 dataset
# Import TxImport Results into DESeq2 structure
dds <- DESeq2::DESeqDataSetFromTximport(txi,sampleTable,~condition)

# Pre-filter data to remove low counts
nrow(dds)
# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep,]
nrow(dds)

# vst transformation is recommended for n > 30 (n = 65)
# blind = FALSE means that the experimental design is not used directly in the transformation
# but is used in estimating the global amount of variability in the counts
vsd <- vst(dds, blind = FALSE)

# Making comparisons between different normalizations for heteroskedasticity
# account for sequencing depth normalization, already done with vst and rld
dds <- estimateSizeFactors(dds)

# Calculate the distance between different samples
sampleDists <- dist(t(assay(vsd)))

library("pheatmap")
library("RColorBrewer")

# For this data set, we see good clustering based on different sample groups
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
# Interesting conclusions: DSMO samples cluster together based on timepoint
# Dox24 hour really stands out on it's own
# Main difference is by timepoint

# PCA analysis
# Use the vsd data rather than original counts data
PCAdata <- plotPCA(vsd, intgroup = c("Drug", "Timepoint"), returnData = TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))

# The results seen in the heatmap are also shown here
ggplot(PCAdata, aes(x = PC1, y = PC2, color = Drug, shape = Timepoint)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

# Map ENSEMBL IDs from kallisto mappings to symbols and entrez IDs
# Load rat annotation information
library("AnnotationDbi")
library("org.Rn.eg.db")

ens.str <- sapply(strsplit(rownames(assay(vsd)), split = ".", fixed = TRUE), '[',1)
convertIDs <- data.frame(kallisto = rownames(assay(vsd)))
convertIDs$symbol <- mapIds(org.Rn.eg.db,
                            keys=ens.str,
                            column="SYMBOL",
                            keytype="ENTREZID",
                            multiVals="first")
convertIDs$entrez <- mapIds(org.Rn.eg.db,
                            keys=ens.str,
                            column="ENTREZID",
                            keytype="ENTREZID",
                            multiVals="first")

#### Running differential expression analysis with DESeq2
dds <- DESeq(dds)

# Save the raw counts to a file for GEO upload
data.to.save = counts(dds)
colnames(data.to.save) = sampleTable$condition
write.csv(data.to.save, 
         file = "data/RNA-seq/raw_gene_counts.csv")

# results() function will print out comparison of last level to first level
# columns of the results data frame
# baseMean - average, normalized, and divided by size factor counts across all samples
# log2FoldChange
# lfcSE - standard error estimate for the log2FC
# stat - Wald statistic
# pvalue - p-value before correction
# padj - BH/FDR corrected p-value



######################## Save results between the three methods ###################################

# For inputs to TIMBR, each condition/timepoint combination should be saved as a separate file
# gene_id: Entrez gene identifier (integer)
# logfc: log2 fold changes represent average differences 
# in log2-transformed gene expression values between treatment and control samples
# fdr: false discovery-rate adjusted q-values (FDR) represent 
# statistical significance of gene expression changes. 
# For gene expression changes, fdr < 0.1 was considered
# significantly differentially expressed. 



# necessary column names: etz_gene,log2FoldChange,pvalue,padj,fdr,baseMean
# need to be named: gene_id, logfc, pval, fdr, ave

###### Save DESeq2 results for TIDEs analysis #########
path.directory <- 'data/RNA-seq/'

data.frame(DESeq2::results(dds,
                           contrast = c("condition","Ace_6","DMSO2_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  #mutate(EntrezID = as.character(EntrezID)) %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_deseq2.csv"), sep = "\t", quote = F, row.names = F)
# Ace 24hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Ace_24","DMSO2_24"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_ace_gene_deseq2.csv"), sep = "\t", quote = F, row.names = F)
# Dox 6hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Dox_6","DMSO2_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_dox_gene_deseq2.csv"), sep = "\t", quote = F, row.names = F)
# Dox 24hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Dox_24","DMSO2_24"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_dox_gene_deseq2.csv"), sep = "\t", quote = F, row.names = F)
# 5FU 6hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","5FU_6","DMSO1_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_5fu_gene_deseq2.csv"), sep = "\t", quote = F, row.names = F)
# 5FU 24hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","5FU_24","DMSO1_24"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_5fu_gene_deseq2.csv"), sep = "\t", quote = F, row.names = F)


# Save differential expression calls for GEO database
# file: differential_expression.xlsx
# colnames: condition, treatment, timepoint, logFC, pvalue, padj
all.results = rbind(read.csv(file = "data/RNA-seq/dougherty_rno_cardio_t24_5fu_gene_deseq2.csv", sep = "\t") %>% 
                      mutate(condition = "5FU_24", treatment = "5FU", timepoint = "24h"), 
                    read.csv(file = "data/RNA-seq/dougherty_rno_cardio_t24_ace_gene_deseq2.csv", sep = "\t") %>% 
                      mutate(condition = "Ace_24", treatment = "Ace", timepoint = "24h")) %>% 
  rbind(read.csv(file = "data/RNA-seq/dougherty_rno_cardio_t24_dox_gene_deseq2.csv", sep = "\t") %>% 
          mutate(condition = "Dox_24", treatment = "Dox", timepoint = "24h")) %>% 
  rbind(read.csv(file = "data/RNA-seq/dougherty_rno_cardio_t6_5fu_gene_deseq2.csv", sep = "\t") %>% 
          mutate(condition = "5FU_6", treatment = "5FU", timepoint = "6h")) %>% 
  rbind(read.csv(file = "data/RNA-seq/dougherty_rno_cardio_t6_ace_gene_deseq2.csv", sep = "\t") %>% 
          mutate(condition = "Ace_6", treatment = "Ace", timepoint = "6h")) %>% 
  rbind(read.csv(file = "data/RNA-seq/dougherty_rno_cardio_t6_dox_gene_deseq2.csv", sep = "\t") %>% 
          mutate(condition = "Dox_6", treatment = "Dox", timepoint = "6h")) %>% 
  dplyr::select(condition, treatment, timepoint, EntrezID, logfc, pval, fdr)
write.csv(all.results, 
           file = "data/RNA-seq/differential_expression.xlsx", 
          row.names = FALSE)

