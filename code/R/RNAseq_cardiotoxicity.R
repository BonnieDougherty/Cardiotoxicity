
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


# Code taken from: 
# DESeq2: https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#removing-hidden-batch-effects
# edge R: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# limma/voom: https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html


# Methods expect inputs as genes x samples
# DESeq2 accounts for library size internally, provide unnormalized counts
# edgeR accounts for library size internally, provide unnormalized counts

# Create file names to get abundance files from Kallisto 
dir <- "F:/Transcriptomics/kallisto/outputs/Cardiotoxicity"

# Check for all files in the directory to ensure they are present 
list.files(dir)

# Read in meta data table
metadata <- xlsx::read.xlsx(file = "F:/Transcriptomics/kallisto/outputs/Cardiotoxicity/metadata.xlsx", header = TRUE, sheetIndex =1)

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

# Use Ensembl identifiers to map to transcripts
#tx2gene <- t2g6[,c(7,2)] %>% unique()

# Use Entrez IDs for mapping to transcripts
tx2gene <- t2g6[,c(7,4)] %>% filter(!is.na(etz_gene)) %>% unique()

# Note: when mapping to Entrez IDs, lose 4037 transcripts from the original dataset
# Note: when mapping to Ensemble IDs, multiple entries for one Entrez ID. 

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

# Use the tximport function to account for differences in transcript length (allows comparison between genes in a sample)
# and differences in library size (allows comparison between samples) 
lengthScaledTPM <- tximport::tximport(files, 
                                      type = "kallisto", 
                                      txOut = FALSE,
                                      countsFromAbundance = "scaledTPM",
                                      tx2gene = tx2gene)
lengthScaledTPM <- DESeq2::DESeqDataSetFromTximport(lengthScaledTPM,sampleTable,~condition)

# Create Gene level analysis in DESeq2 format
lengthScaledTPM <- DESeq2::DESeq(lengthScaledTPM, minReplicatesForReplace = Inf)

# Pull out normalized counts for GSEA
data.to.save <- data.frame(NAME = rownames(DESeq2::counts(lengthScaledTPM, normalized = TRUE)),
                           DESCRIPTION = "na",
                           DESeq2::counts(lengthScaledTPM, normalized=TRUE))
write.table(x = data.to.save, 
            file = "C:/Users/bvd5nq/Documents/GSEA/cardiotoxicity/GSEA_counts_scaledTPMs.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Create the cls file for GSEA analysis
fileConn <- file("C:/Users/bvd5nq/Documents/GSEA/GEOdata/cardiotoxicity_phenotype.cls")
writeLines(c("65 10 1",
             "# 5FU_24 5FU_6 Ace_24 Ace_6 DMSO1_24 DMSO1_6 DMSO2_24 DMSO2_6 Dox_24 Dox_6",
             paste(unlist(condition), collapse = " ")), 
           fileConn)
close(fileConn)

# txi variable is used for all approaches
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

# PCA on non-log transformed counts means that genes with the highest counts will determine the spread of the data
# DESeq2 available functions to account for heteroskedasticity: vst, rlog

# Demonstrate the hetereoskedasticity of the data
meanSdPlot(DESeq2::counts(dds), ranks = FALSE)

# Log transformation of counts
meanSdPlot(log2(DESeq2::counts(dds)+1), ranks = FALSE)

# vst transformation is recommended for n > 30 (n = 65)
# blind = FALSE means that the experimental design is not used directly in the transformation
# but is used in estimating the global amount of variability in the counts
vsd <- vst(dds, blind = FALSE)


# Making comparisons between different normalizations for heteroskedasticity
# account for sequencing depth normalization, already done with vst and rld
dds <- estimateSizeFactors(dds)

# Only taking the first two samples
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

# These plots demonstrate that the vst accounts for changes at low counts - not as much deviation
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

# Calculate the distance between different samples
sampleDists <- dist(t(assay(vsd)))
sampleDists

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

# results() function will print out comparison of last level to first level
# columns of the results data frame
# baseMean - average, normalized, and divided by size factor counts across all samples
# log2FoldChange
# lfcSE - standard error estimate for the log2FC
# stat - Wald statistic
# pvalue - p-value before correction
# padj - BH/FDR corrected p-value

# Pulling out results from the dds object
res <- results(dds, contrast = c("condition","Dox_24","DMSO2_24"))
summary(res)

### Analyzing specific gene counts
# pull out the most significantly changes gene, determined by p-value
topGene <- rownames(res)[which.min(res$padj)]

# gene is the top gene to plot
plotCounts(dds, gene = topGene, intgroup=c("Drug"))

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition"),
                         returnData = TRUE)
geneCounts$condition <- factor(geneCounts$condition, levels = c('DMSO1_6','5FU_6','DMSO1_24','5FU_24', 
                                                                'DMSO2_6', 'Ace_6', 'Dox_6', 'DMSO2_24', 'Ace_24', 'Dox_24'))

gene.symbol <- convertIDs %>% filter(kallisto == topGene)
ggplot(geneCounts, aes(x = condition, y = count)) +
  scale_y_log10() + 
  geom_jitter(cex = 3) + 
  ggtitle(gene.symbol$symbol)

### Gene clustering, most variable genes from the transcript abundance data
library("genefilter")
library(pheatmap)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)

# Need to rename the genes from ENSEMBL
# Heat map of normalized transcript abundances, not DEGs
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

gene.symbols <- data.frame(kallisto = rownames(mat)) %>% 
  left_join(convertIDs)
rownames(mat) <- gene.symbols$symbol

anno <- as.data.frame(colData(vsd)[, c("Drug","Timepoint")])
pheatmap(mat, annotation_col = anno)

resOrdered <- res[order(res$pvalue),]
head(resOrdered)


#### Saving results from dds studies
# Save individual results files to excel files for future analysis
# independentFiltering = TRUE - removes genes with low counts from the analysis
# Shouldn't affect things too much since already filtered low gene counts before analysis
# Create a results table for the DEGs for 5FU_24hrs
Ace6hrs_DEGs <- data.frame(DESeq2::results(dds,
                                           contrast = c("condition","Dox_6","DMSO2_6"), 
                                           alpha = 0.1, 
                                           pAdjustMethod = "BH", 
                                           independentFiltering = TRUE)) %>% 
  rownames_to_column("EntrezID")


######################## edgeR for DEG analysis ##########################################################
# Create DGEList object for use with edgeR
y <- DGEList(counts=txi$counts, samples = sampleTable, group = condition)

# Filtering out lowly expressed genes
nrow(y)
# at least 3 samples with a count of 10 or higher
keep <- rowSums(y$counts >= 10) >= 6
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y)

# Recalculate library sizes after filtering for lowly expressed genes
# Calculate and apply library normalization
y <- calcNormFactors(y)

# Calculate TMM values for GSEA
#TMM Normalisation
lib.size <- colSums(y$counts)
scale.factors <- calcNormFactors(y, method = "TMM")
TMM <- y$counts
TMM <- sweep(TMM,MARGIN=2,FUN="/",STATS=colSums(TMM)*y$samples$norm.factors)*10e6

# Write file for downstream GSEA analysis
write.table(x = TMM, 
            file = "C:/Users/bvd5nq/Documents/GSEA/cardiotoxicity/GSEA_TMMs.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# MDS plot of the raw counts
# See the same separation in the MDS plot as we saw in PCA from the DESeq2 data
plotMDS(y, top = 1000, col = as.numeric(y$samples$Drug), 
        pch = as.numeric(y$samples$Timepoint))


# Create a design matrix
# 0 indicates intercept
design <- model.matrix(~0 + condition, y$samples)


# Estimate dispersion factors
# Specify the design matrix for designs other than one factor
y <- estimateDisp(y, design)

# make contrasts for comparisons between different sample groups
my.contrasts <- makeContrasts(Dox6 = conditionDox_6 - conditionDMSO2_6, 
                              Dox24 = conditionDox_24 - conditionDMSO2_24,
                              Ace6 = conditionAce_6 - conditionDMSO2_6,
                              Ace24 = conditionAce_24 - conditionDMSO2_24,
                              FiveFU6 = condition5FU_6 - conditionDMSO1_6, 
                              FiveFU24= condition5FU_24 - conditionDMSO1_24,
                              levels = design)

# Both of the following approaches use GLMs to estimate DEGs
# In general, glmLRT is more liberal than glmQLFit

# Testing for differential expression
fit <- glmFit(y, design)
# specify the coef to test for
lrt <- glmLRT(fit, contrast = my.contrasts[,'Ace6'])
tt <- topTags(lrt, n=nrow(y), adjust.method = "BH")
tt10 <- topTags(lrt) # just the top 10 by default
tt10

# plot the average log CPM vs. logFC for Ace 6hrs vs. DMSO2 6 hrs
plotMD(lrt)


# goanna() for GO analysis
# kegga() for KEGG pathway analysis
go <- goana(lrt,
            species="Rn")
topGO(go, sort="up")
keg <- kegga(lrt, 
             species="Rn")
topKEGG(keg, sort="up")


############################### limma/voom for DEG analysis #######################################
y.limma <- DGEList(counts=txi$counts, samples = sampleTable, group = condition)

# Convert raw counts to CPMs and logCPMs
cpm <- cpm(y.limma)
lcpm <- cpm(y.limma, log=TRUE)

# Filtering out lowly expressed genes
nrow(y.limma)
# at least 3 samples with a count of 10 or higher
keep <- rowSums(y.limma$counts >= 10) >= 6
y.limma <- y.limma[keep, , keep.lib.sizes=FALSE]
nrow(y.limma)

# Accounting for differences in library size
y.limma <- calcNormFactors(y.limma, method = "TMM")
y.limma$samples$norm.factors

# MDS plots of data
col.group <- sampleTable$Drug
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
# color by drug
plotMDS(y.limma, labels = condition, col = col.group)


# Create a design matrix
# 0 indicates intercept
design <- model.matrix(~0 + condition, y.limma$samples)

# make contrasts for comparisons between different sample groups
my.contrasts <- makeContrasts(Dox6 = conditionDox_6 - conditionDMSO2_6, 
                              Dox24 = conditionDox_24 - conditionDMSO2_24,
                              Ace6 = conditionAce_6 - conditionDMSO2_6,
                              Ace24 = conditionAce_24 - conditionDMSO2_24,
                              FiveFU6 = condition5FU_6 - conditionDMSO1_6, 
                              FiveFU24= condition5FU_24 - conditionDMSO1_24,
                              levels = design)


# Removing heteroskedasticity from the data
v <- voom(y.limma, design, plot=TRUE)

# determining DEGs based on the contrasts
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=my.contrasts)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

dt <- decideTests(efit, adjust.method = "BH", p.value = 0.1)
summary(dt)

# Plot the relationship between average log-expression and log fold change 
plotMD(efit, column=3, status=dt[,3])

# extracting results using topTable() 
results.limma <- topTable(efit, number = 20)

# Heatmap of most significant genes across all conditions
library(gplots)
topGenes <- data.frame(entrez = rownames(results.limma)) %>% 
  inner_join(convertIDs)
i <- which(rownames(y.limma$counts) %in% topGenes$entrez)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=topGenes$symbol, labCol=condition, 
          col=mycol, trace="none", density.info="none", 
          #margin=c(8,6), lhei=c(2,10), 
          dendrogram="column")


################### Compare results between the three methods #####################################
# DEG Venn Diagram for shared genes between the three different approaches
# Rank the fold change of genes to see how the methods rank genes
# Plot fold changes against each other
# volcano plots for each individual approach
# Transcript clustering after normalization - what are the most differential transcripts across the different normalization approaches?

# All comparisons were made using Ace 6hrs vs. DMSO2 6hrs data

# Extracting results from the DESeq2 data
# columns: baseMean, log2FoldChange lfcSE, state, pvalue, padj
DESeq2_results <- data.frame(DESeq2::results(dds,
                                           contrast = c("condition","Dox_24","DMSO2_24"), 
                                           alpha = 0.1, 
                                           pAdjustMethod = "BH", 
                                           independentFiltering = TRUE))
DESeq2_results$entrez <- rownames(DESeq2_results)
DESeq2_results <- DESeq2_results %>% dplyr::select(entrez, log2FoldChange, padj)
colnames(DESeq2_results) <- c('entrez','DESeq2_logFC','DESeq2_FDR')

# Extracting results from the edgeR data
lrt <- glmLRT(fit, contrast = my.contrasts[,'Dox24'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH"))
edgeR_results$entrez <- rownames(edgeR_results)
edgeR_results <- edgeR_results %>% 
  dplyr::select(entrez, logFC, FDR)
colnames(edgeR_results) <- c('entrez','edgeR_logFC','edgeR_FDR')


# Extractign results from the limma/voom data
# variable is efit
# colnames(efit) for coefficients
limma_results <- data.frame(topTable(efit, 
                                     coef = 2, 
                                     number = dim(v$E)[1], 
                                     adjust.method = "BH"))
limma_results$entrez <- rownames(limma_results)
limma_results <- limma_results %>% 
  dplyr::select(entrez, logFC, adj.P.Val)
colnames(limma_results) <- c('entrez','limma_logFC','limma_FDR')

# Different number of genes between the three approaches
DESeq2_results <- DESeq2_results %>% inner_join(edgeR_results %>% dplyr::select(entrez))

all.results <- inner_join(DESeq2_results, edgeR_results) %>% 
  inner_join(limma_results)

# Compare results with DESeq2 results 
table(DESeq2 = all.results$DESeq2_FDR < 0.1, edgeR = all.results$edgeR_FDR < 0.1)
table(DESeq2 = all.results$DESeq2_FDR < 0.1, limma = all.results$limma_FDR < 0.1)
table(edgeR = all.results$edgeR_FDR < 0.1, limma = all.results$limma_FDR < 0.1)


# Generating a Venn Diagram of shared DEGs
library(VennDiagram)
n123 = (all.results %>% filter(DESeq2_FDR < 0.1) %>% filter(edgeR_FDR < 0.1) %>% filter(limma_FDR < 0.1) %>% dplyr::count())$n
n13 = (all.results %>% filter(DESeq2_FDR < 0.1) %>% filter(limma_FDR < 0.1) %>% dplyr::count())$n
n23 = (all.results %>%  filter(edgeR_FDR < 0.1) %>% filter(limma_FDR < 0.1) %>% dplyr::count())$n
n12 = (all.results %>% filter(DESeq2_FDR < 0.1) %>% filter(edgeR_FDR < 0.1) %>% dplyr::count())$n
area1 = (all.results %>% filter(DESeq2_FDR < 0.1) %>% dplyr::count())$n
area2 = (all.results %>% filter(edgeR_FDR < 0.1) %>% dplyr::count())$n
area3 = (all.results %>% filter(limma_FDR < 0.1) %>% dplyr::count())$n
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123,
                 category = c("DESeq2","edgeR","limma/voom"), 
                 fill = c('#404788FF','#73D055FF','#FDE725FF'), 
                 cat.fontfamily = c('bold','bold','bold'), 
                 cat.cex = c(1.5, 1.5, 1.5), 
                 cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5))



# Comparing the ranks of genes between different conditions
# Want to have a facet_Wrap plot comparing all of the different combinations (DESeq)
# Need to make these combinations individually and then merge based on entrez ID, add new variable for facet_Wrap
rank.results <- rbind(data.frame(x = rank(all.results$DESeq2_FDR), y = rank(all.results$edgeR_FDR), data = "DESeq2 vs. edgeR"), 
                      data.frame(x = rank(all.results$DESeq2_FDR), y = rank(all.results$limma_FDR), data = "DESeq2 vs. limma"), 
                      data.frame(x = rank(all.results$edgeR_FDR), y = rank(all.results$limma_FDR), data = "edgeR vs. limma"))
Dox24hrs_rank <- ggplot(rank.results, aes(x = x, y = y)) + 
  geom_point(cex = 0.1) + 
  facet_wrap(~data) + 
  theme_bw()

ggsave(filename = 'C:/Users/bvd5nq/Documents/R scripts/Analyzing RNA-seq data/figures/Dox24hrs_ranks.png', 
       plot = Dox24hrs_rank, 
       dpi = 600, width = 10, height = 3.3, units = 'in')


# facet_wrapped volcano plots
# logFC vs. -log10(p-value), doesn't matter if p-value or FDR
volcano.data <- rbind(data.frame(logFC = all.results$DESeq2_logFC, p.val = all.results$DESeq2_FDR, approach = "DESeq2"), 
                      data.frame(logFC = all.results$edgeR_logFC, p.val = all.results$edgeR_FDR, approach = "edgeR"), 
                      data.frame(logFC = all.results$limma_logFC, p.val = all.results$limma_FDR, approach = "limma")) %>% 
  mutate(significance = p.val < 0.1)
Dox24hrs_volcano <- ggplot(volcano.data, aes(x = logFC, y = -log10(p.val), color = significance)) + 
  geom_point(cex = 0.1) +
  facet_wrap(~approach, scales = "free") + 
  theme_bw() + 
  theme(legend.position = "none")

ggsave(filename = 'C:/Users/bvd5nq/Documents/R scripts/Analyzing RNA-seq data/figures/Dox24hrs_volcano.png', 
       plot = Dox24hrs_volcano, 
       dpi = 600, width = 10, height = 3.3, units = 'in')


# Compare the logFCs for all genes between the different conditions
FC.data <- rbind(data.frame(x = all.results$DESeq2_logFC, y = all.results$edgeR_logFC, approach = "DESeq2 vs. edgeR"), 
                 data.frame(x = all.results$DESeq2_logFC, y = all.results$limma_logFC, approach = "DESeq2 vs. limma"), 
                 data.frame(x = all.results$edgeR_logFC, y = all.results$limma_logFC, approach = "edgeR vs. limma"))
Dox24hrs_FC <- ggplot(FC.data, aes(x = x, y = y)) + 
  geom_point(cex = 0.1) +
  facet_wrap(~approach) + 
  theme_bw() 

ggsave(filename = 'C:/Users/bvd5nq/Documents/R scripts/Analyzing RNA-seq data/figures/Dox24hrs_FC.png', 
       plot = Dox24hrs_FC, 
       dpi = 600, width = 10, height = 3.3, units = 'in')



# Determines the top variable genes based on the transcript data, after library normalization
# DESeq2: assay(vsd) (accounts for heteroskedasticity)
# edgeR: y$counts
# limma: v$E (accounts for heteroskedasticity)

# DESeq2
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat <- assay(vsd)[topVarGenes, ]
mat  <- mat - rowMeans(mat)
gene.symbols <- data.frame(kallisto = rownames(mat)) %>% 
  left_join(convertIDs)
rownames(mat) <- gene.symbols$symbol

anno <- as.data.frame(colData(vsd)[, c("Drug","Timepoint")])
pheatmap(mat, annotation_col = anno)


# edgeR
topVarGenes <- head(order(rowVars(y$counts), decreasing = TRUE), 20)
mat <- y$counts[topVarGenes, ]
mat <- mat - rowMeans(mat)
gene.symbols <- data.frame(kallisto = rownames(mat)) %>% 
  left_join(convertIDs)
rownames(mat) <- gene.symbols$symbol

anno <- y$samples %>% dplyr::select(Drug, Timepoint)
pheatmap(mat, annotation_col = anno)


# limma/voom
topVarGenes <- head(order(rowVars(v$E), decreasing = TRUE), 20)
mat <- v$E[topVarGenes, ]
mat <- mat - rowMeans(mat)
gene.symbols <- data.frame(kallisto = rownames(mat)) %>% 
  left_join(convertIDs)
rownames(mat) <- gene.symbols$symbol

anno <- v$targets %>% dplyr::select(Drug, Timepoint)
pheatmap(mat, annotation_col = anno)


# Notes:
# Differences in fold changes would affect predictions made by TIMBR and TIDEs
# Support for including normalized transcripts counts rather than log2FCs
# Transcript counts can be included using RIPTIDE -> place the objective function as constraints based on metabolic tasks
# -> Is there a measure of how well RIPTIDE fits the data? It will always find a solution so....
# -> Disagreements between data and reactions could indicate control points!

# RIPTIDE could give an indication of how many different pathways are available to complete tasks.


######################## Save results between the three methods ###################################

# For inputs to TIMBR, each condition/timepoint combination should be saved as a separate file
# gene_id: Entrez gene identifier (integer)
# logfc: log2 fold changes represent average differences 
# in log2-transformed gene expression values between treatment and control samples
# fdr: false discovery-rate adjusted q-values (FDR) represent 
# statistical significance of gene expression changes. 
# For gene expression changes, fdr < 0.1 was considered
# significantly differentially expressed. 


# The filename is based on the following information separated by "_": 
#   organism (hsa = human; rno = rat)
#   cell type (cardio = cardiomyocytes)
#   treatment duration (t6 = 6 hours; t24 = 24 hours)
#   treatment compound (e.g. caffeine; acetaminophen)
#   feature type (gene)
#   DEG method (deseq2; edgeR; limma)

# Save as .txt files

path.directory <- 'C:/Users/bvd5nq/Documents/R scripts/Analyzing RNA-seq data/RNA-seq DEG data/'



# necessary column names: etz_gene,log2FoldChange,pvalue,padj,fdr,baseMean
# need to be named: gene_id, logfc, pval, fdr, ave

# DESeq2
# Ace 6hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Ace_6","DMSO2_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  #mutate(EntrezID = as.character(EntrezID)) %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_deseq2.csv"), sep = ",", quote = F, row.names = F)
# Ace 24hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Ace_24","DMSO2_24"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_ace_gene_deseq2.txt"), sep = "\t", quote = F, row.names = F)
# Dox 6hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Dox_6","DMSO2_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_dox_gene_deseq2.txt"), sep = "\t", quote = F, row.names = F)
# Dox 24hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","Dox_24","DMSO2_24"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_dox_gene_deseq2.txt"), sep = "\t", quote = F, row.names = F)
# 5FU 6hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","5FU_6","DMSO1_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_5fu_gene_deseq2.txt"), sep = "\t", quote = F, row.names = F)
# 5FU 24hrs
data.frame(DESeq2::results(dds,
                           contrast = c("condition","5FU_24","DMSO1_24"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_5fu_gene_deseq2.txt"), sep = "\t", quote = F, row.names = F)








# edgeR
# Ace 6hrs
lrt <- glmLRT(fit, contrast = my.contrasts[,'Ace6'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, FDR, PValue) %>% 
  dplyr::rename("logfc" = logFC,"pval" = PValue, "fdr" = FDR) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_edgeR.txt"), sep = "\t", quote = F, row.names = F)
# Ace 24hrs
lrt <- glmLRT(fit, contrast = my.contrasts[,'Ace24'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, FDR, PValue) %>% 
  dplyr::rename("logfc" = logFC,"pval" = PValue, "fdr" = FDR) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_ace_gene_edgeR.txt"), sep = "\t", quote = F, row.names = F)


# Dox 6hrs
lrt <- glmLRT(fit, contrast = my.contrasts[,'Dox6'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, FDR, PValue) %>% 
  dplyr::rename("logfc" = logFC,"pval" = PValue, "fdr" = FDR) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_dox_gene_edgeR.txt"), sep = "\t", quote = F, row.names = F)
# Dox 24hrs
lrt <- glmLRT(fit, contrast = my.contrasts[,'Dox24'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, FDR, PValue) %>% 
  dplyr::rename("logfc" = logFC,"pval" = PValue, "fdr" = FDR) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_dox_gene_edgeR.txt"), sep = "\t", quote = F, row.names = F)

# 5FU 6hrs
lrt <- glmLRT(fit, contrast = my.contrasts[,'FiveFU6'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, FDR, PValue) %>% 
  dplyr::rename("logfc" = logFC,"pval" = PValue, "fdr" = FDR) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_5fu_gene_edgeR.txt"), sep = "\t", quote = F, row.names = F)
# 5FU 24hrs
lrt <- glmLRT(fit, contrast = my.contrasts[,'FiveFU24'])
edgeR_results <- data.frame(topTags(lrt, 
                                    n = nrow(y), 
                                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, FDR, PValue) %>% 
  dplyr::rename("logfc" = logFC,"pval" = PValue, "fdr" = FDR) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_5fu_gene_edgeR.txt"), sep = "\t", quote = F, row.names = F)




# limma results
# colnames(efit) for coefficients
# Ace 6hrs
data.frame(topTable(efit,
                    coef = 3, 
                    number = dim(v$E)[1],
                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, AveExpr, adj.P.Val, P.Value) %>% 
  dplyr::rename("logfc" = logFC, "fdr" = adj.P.Val, "pval" = P.Value, "ave" = AveExpr) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_limma.csv"), sep = ",", quote = F, row.names = F)
# Ace 24hrs
data.frame(topTable(efit,
                    coef = 4, 
                    number = dim(v$E)[1],
                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, AveExpr, adj.P.Val, P.Value) %>% 
  dplyr::rename("logfc" = logFC, "fdr" = adj.P.Val, "pval" = P.Value, "ave" = AveExpr) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_ace_gene_limma.txt"), sep = "\t", quote = F, row.names = F)


# Dox 6hrs
data.frame(topTable(efit,
                    coef = 1, 
                    number = dim(v$E)[1],
                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, AveExpr, adj.P.Val, P.Value) %>% 
  dplyr::rename("logfc" = logFC, "fdr" = adj.P.Val, "pval" = P.Value, "ave" = AveExpr) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_dox_gene_limma.txt"), sep = "\t", quote = F, row.names = F)
# Dox 24hrs
data.frame(topTable(efit,
                    coef = 2, 
                    number = dim(v$E)[1],
                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, AveExpr, adj.P.Val, P.Value) %>% 
  dplyr::rename("logfc" = logFC, "fdr" = adj.P.Val, "pval" = P.Value, "ave" = AveExpr) %>%  
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_dox_gene_limma.txt"), sep = "\t", quote = F, row.names = F)


# 5FU 6hrs
data.frame(topTable(efit,
                    coef = 5, 
                    number = dim(v$E)[1],
                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, AveExpr, adj.P.Val, P.Value) %>% 
  dplyr::rename("logfc" = logFC, "fdr" = adj.P.Val, "pval" = P.Value, "ave" = AveExpr) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_5fu_gene_limma.txt"), sep = "\t", quote = F, row.names = F)
# 5FU 24hrs
data.frame(topTable(efit,
                    coef = 6, 
                    number = dim(v$E)[1],
                    adjust.method = "BH")) %>% 
  tibble::rownames_to_column(var = "EntrezID") %>% 
  dplyr::select(EntrezID, logFC, AveExpr, adj.P.Val, P.Value) %>% 
  dplyr::rename("logfc" = logFC, "fdr" = adj.P.Val, "pval" = P.Value, "ave" = AveExpr) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t24_5fu_gene_limma.txt"), sep = "\t", quote = F, row.names = F)

### Formating DESeq2 results for TIDEs analysis
# Load transcriptomics data for DESeq2 comparison
path.directory <- 'C:/Users/bvd5nq/Documents/R scripts/Analyzing RNA-seq data/RNA-seq DEG data/'

save.directory <- 'C:/Users/bvd5nq/Documents/R scripts/Analyzing RNA-seq data/TIDEs/DEGs/'

# information needed: gene name
# Ace_6h
read.table(file = paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_deseq2.txt"), 
                                   header = TRUE) %>% 
  mutate(logfc = ifelse(fdr < 0.1, logfc, 0), 
         EntrezID = as.character(EntrezID)) %>% 
  write.table(paste0(save.directory, "dougherty_rno_cardio_t6_ace_gene_deseq2.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)


# Ace_24h
read.table(file = paste0(path.directory, "dougherty_rno_cardio_t24_ace_gene_deseq2.txt"), 
           header = TRUE) %>% 
  mutate(logfc = ifelse(fdr < 0.1, logfc, 0), 
         EntrezID = as.character(EntrezID)) %>% 
  write.table(paste0(save.directory, "dougherty_rno_cardio_t24_ace_gene_deseq2.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)


# Dox_6h
read.table(file = paste0(path.directory, "dougherty_rno_cardio_t6_dox_gene_deseq2.txt"), 
           header = TRUE) %>% 
  mutate(logfc = ifelse(fdr < 0.1, logfc, 0),
         EntrezID = as.character(EntrezID)) %>% 
  write.table(paste0(save.directory, "dougherty_rno_cardio_t6_dox_gene_deseq2.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)


# Dox_24h
read.table(file = paste0(path.directory, "dougherty_rno_cardio_t24_dox_gene_deseq2.txt"), 
           header = TRUE) %>% 
  mutate(logfc = ifelse(fdr < 0.1, logfc, 0), 
         EntrezID = as.character(EntrezID)) %>% 
  write.table(paste0(save.directory, "dougherty_rno_cardio_t24_dox_gene_deseq2.csv"),
              sep = ",", col.names = TRUE, row.names = FALSE)


# 5FU_6h
read.table(file = paste0(path.directory, "dougherty_rno_cardio_t6_5fu_gene_deseq2.txt"), 
           header = TRUE) %>% 
  mutate(logfc = ifelse(fdr < 0.1, logfc, 0),
         EntrezID = as.character(EntrezID)) %>% 
  write.table(paste0(save.directory, "dougherty_rno_cardio_t6_5fu_gene_deseq2.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)


# 5FU_24h
read.table(file = paste0(path.directory, "dougherty_rno_cardio_t24_5fu_gene_deseq2.txt"), 
           header = TRUE) %>% 
  mutate(logfc = ifelse(fdr < 0.1, logfc, 0),
         EntrezID = as.character(EntrezID)) %>% 
  write.table(paste0(save.directory, "dougherty_rno_cardio_t24_5fu_gene_deseq2.csv"), 
              sep = ",", col.names = TRUE, row.names = FALSE)



###### Save DESeq2 results for TIDEs analysis #########
path.directory = "C:/Users/bvd5nq/Documents/Cardiotoxicity/data/"

data.frame(DESeq2::results(dds,
                           contrast = c("condition","Ace_6","DMSO2_6"), 
                           alpha = 0.1, 
                           pAdjustMethod = "BH",
                           independentFiltering = TRUE)) %>%
  tibble::rownames_to_column(var = "EntrezID") %>% 
  #mutate(EntrezID = as.character(EntrezID)) %>% 
  dplyr::select(EntrezID, log2FoldChange, pvalue, padj, baseMean) %>% 
  dplyr::rename("logfc" = log2FoldChange, "pval" = pvalue, "fdr" = padj, "ave" = baseMean) %>% 
  write.table(paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_deseq2.csv"), sep = ",", quote = F, row.names = F)
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
