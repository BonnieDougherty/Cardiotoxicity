# Generating figures for the cardiotoxicity manuscript

library(xlsx)
library(ggplot2)
library(viridis)
library(forcats)
library(tidyr)
library(tibble)

library(readtext)

library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(AUCRF)

library(clusterProfiler)
library(DOSE)
library(biomaRt)
library(msigdbr)

library(forcats)

library(dplyr)
library(vegan)
library(DESeq2)
library(stringr)


#### Figure 2 ####
# Load in raw transcript count data for PCA plot
load(file = "data/RNA-seq/PCA_transcript.R")

# blind = FALSE means that the experimental design is not used directly in the transformation
# but is used in estimating the global amount of variability in the counts
vsd <- vst(dds, blind = FALSE)

transcript.counts = read.table("data/RNA-seq/scaledTPMs.txt", header = TRUE) %>% 
  select(-DESCRIPTION) %>% 
  pivot_longer(-NAME, names_to = "sample", values_to = "counts")

sample.metadata = data.frame(sample = c('sample1','sample2','sample3','sample4','sample5','sample6',
                                           'sample7','sample8','sample9','sample10','sample11','sample12',
                                           'sample13','sample14','sample15','sample16','sample17','sample18',
                                           'sample19','sample20','sample21','sample22','sample23','sample24',
                                           'sample25','sample26','sample27','sample28','sample29','sample30',
                                           'sample31','sample32','sample33','sample34','sample35','sample36',
                                           'sample37','sample38','sample39','sample40','sample41','sample42',
                                           'sample43','sample44','sample45','sample46','sample47','sample48',
                                           'sample49','sample50','sample51','sample52','sample53','sample54',
                                           'sample55','sample56','sample57','sample58','sample59','sample60',
                                           'sample61','sample62','sample63','sample64','sample65'), 
                                group = c("5FU_24h","5FU_24h","5FU_24h","5FU_24h","5FU_24h","5FU_24h",
                                            "5FU_6h", "5FU_6h", "5FU_6h", "5FU_6h", "5FU_6h", "5FU_6h",
                                            "Ace_24h", "Ace_24h", "Ace_24h", "Ace_24h", "Ace_24h", "Ace_24h", "Ace_24h",
                                            "Ace_6h","Ace_6h","Ace_6h","Ace_6h","Ace_6h","Ace_6h",
                                            "DMSO1_24h","DMSO1_24h","DMSO1_24h","DMSO1_24h","DMSO1_24h","DMSO1_24h",
                                            "DMSO1_6h","DMSO1_6h","DMSO1_6h","DMSO1_6h","DMSO1_6h","DMSO1_6h",
                                            "DMSO2_24h","DMSO2_24h","DMSO2_24h","DMSO2_24h","DMSO2_24h","DMSO2_24h","DMSO2_24h",
                                            "DMSO2_6h","DMSO2_6h","DMSO2_6h","DMSO2_6h","DMSO2_6h","DMSO2_6h","DMSO2_6h",
                                            "Dox_24h","Dox_24h","Dox_24h","Dox_24h","Dox_24h","Dox_24h","Dox_24h",
                                            "Dox_6h","Dox_6h","Dox_6h","Dox_6h","Dox_6h","Dox_6h","Dox_6h"))

condition.metadata <- data.frame(group = c("Ace_6h", "Ace_24h", "DMSO1_6h", "DMSO1_24h", "DMSO2_6h", "DMSO2_24h", "Dox_6h","Dox_24h", "5FU_6h", "5FU_24h"), 
                                 Timepoint = c('6h','24h','6h','24h','6h','24h','6h','24h','6h','24h'), 
                                 Drug = c("Ace","Ace","DMSO1","DMSO1","DMSO2","DMSO2","Dox","Dox","5FU","5FU"))

transcript.counts = transcript.counts %>% 
  left_join(sample.metadata, by = "sample") %>% 
  left_join(condition.metadata, by = "group")

# Using prcomp to plot PCA
PCA.results = prcomp(assay(vsd) %>% t(), scale = TRUE, center = TRUE)
PCA.data <- data.frame(colData(vsd), 
                       PC1 = PCA.results$x[,1], 
                       PC2 = PCA.results$x[,2]) %>% 
  mutate(Drug = as.character(Drug)) %>% 
  mutate(Drug = ifelse(Drug == "DMSO1", "DMSO1 (high)", Drug), 
         Drug = ifelse(Drug == "DMSO2", "DMSO2 (low)", Drug)) %>% 
  mutate(Timepoint = ifelse(Timepoint == 6, "6h", "24h"))

# order for PCA data: PC1, PC2, group, Timepoint, Drug, name, comparison, x, y
PCA.data$Drug <- factor(PCA.data$Drug, levels = c('5FU','DMSO1 (high)','Ace','Dox','DMSO2 (low)'))
PCA.data$Timepoint <- factor(PCA.data$Timepoint, levels = c('6h','24h'))

summary(PCA.results)

colnames(PCA.data)[3] = "Time point"

# # Colors for conditions
v_colors = viridis(5, option = "C")

# Plot the PCA for the transcript data
Figure2A <- ggplot(PCA.data, aes(x = PC1, y = PC2, fill = Drug, shape = `Time point`)) +
  geom_point(cex = 3) + 
  xlab(paste0("PC1: 35% variance")) +
  ylab(paste0("PC2: 19% variance")) + 
  scale_fill_manual(values = v_colors) +
  scale_shape_manual(values=c(21, 22)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  ggtitle("Transcriptomics") + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") 
Figure2A

# Generate a statistical test for differences between groups
gene.dist <- vegdist(assay(vsd) %>% t, method='bray')
temp.data = assay(vsd) %>% t() %>% data.frame() %>% cbind(colData(vsd))
gene.div <- adonis2(gene.dist~condition, data = temp.data, permutations = 999, method="bray")
gene.div

# Run individual comparisons for PERMANOVA testing
# treated groups (5FU and Ace) separate from control groups (DMSO1 and DMSO2) at both time points
temp.data = assay(vsd) %>% t() %>% cbind(colData(vsd)) %>% data.frame() %>% 
  mutate(treatment = ifelse(Drug == "DMSO1" | Drug == "DMSO2", "control", "treatment")) %>% 
  rownames_to_column("sample") %>% 
  filter(Drug != "Dox")

# Test difference between treatment (5FU and Ace) at 6 hours
gene.dist = vegdist(temp.data %>% filter(Timepoint == 6) %>% select(-sample, -condition, -Drug, -Timepoint, -treatment) %>% as.matrix(), 
                    method = "bray")

gene.div = adonis2(gene.dist~treatment, 
                   data = temp.data %>% filter(Timepoint == 6), 
                   permutations = 999, method = "bray")
gene.div

# Test difference between treatment (5FU and Ace) at 24 hours
gene.dist = vegdist(temp.data %>% filter(Timepoint == 24) %>% select(-sample, -condition, -Drug, -Timepoint, -treatment) %>% as.matrix(), 
                    method = "bray")

gene.div = adonis2(gene.dist~treatment, 
                   data = temp.data %>% filter(Timepoint == 24), 
                   permutations = 999, method = "bray")
gene.div

SupplementalFigure3A = ggplot(PCA.data %>% filter(Drug != "Dox"), 
                              aes(x = PC1, y = PC2, fill = Drug, shape = `Time point`)) +
  geom_point(cex = 3) + 
  xlab(paste0("PC1: 35% variance")) +
  ylab(paste0("PC2: 19% variance")) + 
  scale_fill_manual(values = v_colors[c(1,2,3,5)]) +
  scale_shape_manual(values=c(21, 22)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  ggtitle("Transcriptomics") + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") 
SupplementalFigure3A

# Generate loadings plot for PCA transcription results
top_genes <- data.frame(PCA.results$rotation) %>% 
  rownames_to_column("gene") %>% 
  select(gene, PC1, PC2) %>% 
  #mutate(weight = sqrt(PC1^2 + PC2^2)) %>% 
  pivot_longer(-gene, names_to = "PC", values_to = "value") %>% 
  group_by(PC) %>% 
  arrange(desc(abs(value))) %>% 
  slice_head(n = 100)

# Download pathways for analysis
# H - Hallmark gene sets
# C5 - GO gene sets
m_df = msigdbr(species = "Rattus norvegicus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

gene.enrichment.PC1 <- enricher((top_genes %>% filter(PC == "PC1"))$gene,
                            TERM2GENE = m_df,
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.1,
                            minGSSize = 10) %>% 
  data.frame()

gene.enrichment.PC2 <- enricher((top_genes %>% filter(PC == "PC2"))$gene,
                                TERM2GENE = m_df,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.1,
                                minGSSize = 10) %>% 
  data.frame()


# Supplemental Figure 5 - PCA of only metabolic genes included in the model
metabolic.genes = read.xlsx(file = "data/iCardio_rat_rxns.xlsx", sheetIndex = 1, header = 0)
colnames(metabolic.genes) = 'RXNID'
GPR.rules = read.xlsx(file = "data/ncomms14250-s13, iRno RAVEN.xlsx", sheetIndex = 1, header = 1)

metabolic.genes = left_join(metabolic.genes, GPR.rules, by = c('RXNID')) %>% 
  dplyr::select(RXNID, GENE.ASSOCIATION) %>% 
  mutate(GENE = strsplit(as.character(GENE.ASSOCIATION), ";")) %>%
  unnest(GENE) %>%
  select(-GENE.ASSOCIATION) %>% 
  mutate(GENE = strsplit(as.character(GENE), ":")) %>% 
  unnest(GENE) %>% 
  select(GENE) %>% 
  unique()

# Join with the assay(vsd) table and plot PCA only on these genes
genes = rownames(assay(vsd)) %>% data.frame() %>% 
  rownames_to_column("row.ID") %>% 
  inner_join(metabolic.genes, by = c("."="GENE")) %>% 
  mutate(row.ID = as.numeric(row.ID))
  
# Using prcomp to plot PCA
PCA.results = prcomp(assay(vsd)[genes$row.ID,] %>% t(), scale = TRUE, center = TRUE)
PCA.data <- data.frame(colData(vsd), 
                       PC1 = PCA.results$x[,1], 
                       PC2 = PCA.results$x[,2]) %>% 
  mutate(Drug = as.character(Drug)) %>% 
  mutate(Drug = ifelse(Drug == "DMSO1", "DMSO1 (high)", Drug), 
                                                        Drug = ifelse(Drug == "DMSO2", "DMSO2 (low)", Drug)) %>% 
  mutate(Timepoint = ifelse(Timepoint == 6, "6h", "24h"))

# order for PCA data: PC1, PC2, group, Timepoint, Drug, name, comparison, x, y
PCA.data$Drug <- factor(PCA.data$Drug, levels = c('5FU','DMSO1 (high)','Ace','Dox','DMSO2 (low)'))
PCA.data$Timepoint <- factor(PCA.data$Timepoint, levels = c('6h','24h'))

colnames(PCA.data)[3] = "Time point"

# # Colors for conditions
v_colors = viridis(5, option = "C")
# Plot the PCA for the transcript data
SupplementalFigure5 <- ggplot(PCA.data, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Drug, shape = `Time point`), color = "black", cex = 3) + 
  xlab(paste0("PC1: 35% variance")) +
  ylab(paste0("PC2: 21% variance")) + 
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values = v_colors) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme_bw() +
  ggtitle("Transcriptomics") + 
  theme(plot.title = element_text(hjust = 0.5))
SupplementalFigure5

ggsave('results/figures/SupplementalFigure5.png', SupplementalFigure5, 
       dpi = 600, width = 4.5, height = 3, units = 'in')

# enrichment on the top 100 genes separating in the first and second principal components
# Generate loadings plot for PCA transcription results
top_genes <- data.frame(PCA.results$rotation) %>% 
  rownames_to_column("gene") %>% 
  select(gene, PC1, PC2) %>% 
  #mutate(weight = sqrt(PC1^2 + PC2^2)) %>% 
  pivot_longer(-gene, names_to = "PC", values_to = "value") %>% 
  group_by(PC) %>% 
  arrange(desc(abs(value))) %>% 
  slice_head(n = 100)

# Download pathways for analysis
# H - Hallmark gene sets
# C5 - GO gene sets
m_df = msigdbr(species = "Rattus norvegicus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

gene.enrichment.PC1 <- enricher((top_genes %>% filter(PC == "PC1"))$gene,
                                TERM2GENE = m_df,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.1,
                                minGSSize = 10) %>% 
  data.frame()

gene.enrichment.PC2 <- enricher((top_genes %>% filter(PC == "PC2"))$gene,
                                TERM2GENE = m_df,
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.1,
                                minGSSize = 10) %>% 
  data.frame()

  

# Load in raw metabolomics data for PCA plot
scaled.data = xlsx::read.xlsx(file = "results/supplement/Supplemental File 3.xlsx", 
                        sheetName = "scaled", header = TRUE)

# Filter out unnamed metabolites
metabolites <- data.frame(BIOCHEMICAL = scaled.data$BIOCHEMICAL) %>% 
  separate(BIOCHEMICAL, sep = "X - ", into = c("BIOCHEMICAL","junk"),  remove = FALSE) %>% 
  filter(is.na(junk)) %>% 
  select(-junk)

scaled.data <- inner_join(scaled.data, metabolites)

# PCA of the metabolomics data
# remove the blank samples from the PCA
PCA.results <- prcomp(t(scaled.data[,-c(1,67:69)]), scale = FALSE, center = FALSE)

# Format a new data frame for PCA data
# sample name, drug, timepoint, PC1, PC2
metadata.master <- xlsx::read.xlsx(file = "data/metabolomics/metabolomics_metadata.xlsx", sheetIndex = 1) %>% as.data.frame()
metadata.master$sample.ID <- as.character(metadata.master$sample.ID)
metadata.master$Drug = as.character(metadata.master$Drug)
PCA.data <- data.frame(sample.ID = colnames(scaled.data)[-c(1,67:69)], 
                       PC1 = PCA.results$x[,1], 
                       PC2 = PCA.results$x[,2]) %>% 
  inner_join(metadata.master) %>% 
  filter(!(group == "Blank")) %>% 
  mutate(Drug = ifelse(Drug == "DMSO1", "DMSO1 (high)", Drug), 
         Drug = ifelse(Drug == "DMSO2", "DMSO2 (low)", Drug))

# order for PCA data: PC1, PC2, group, Timepoint, Drug, name, comparison, x, y
PCA.data$Drug <- factor(PCA.data$Drug, levels = c('5FU','DMSO1 (high)','Ace','Dox','DMSO2 (low)'))
PCA.data$Timepoint <- factor(PCA.data$Timepoint, levels = c('6h','24h'))
PCA.data <- PCA.data %>% mutate(comparison = "Metabolomics") %>% 
  mutate(x = 45, y = 13) %>% 
  select(PC1, PC2, Drug, Timepoint, comparison, x, y)

colnames(PCA.data)[4] = "Time point"

# Colors for conditions
v_colors = viridis(5, option = "C")

Figure2B <- ggplot(PCA.data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Drug, shape = `Time point`), cex = 3) + 
  theme_bw() +
  xlab("PC1: 45% variance") + 
  ylab("PC2: 13% variance") + 
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values = v_colors) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  ggtitle("Metabolomics") + 
  theme(plot.title = element_text(hjust = 0.5))
Figure2B

# Combine Figure2A and B into one figure
gridlayout <- data.matrix(c(1,NA,2)) %>% t()
Figure2 <- arrangeGrob(Figure2A, Figure2B,
                       ncol = 3, 
                       nrow = 1, 
                       #heights = relative_heights,
                       widths = c(0.36+0.035,0.05,0.54+0.015),
                       layout_matrix = gridlayout) 

ggsave('results/figures/Figure2.png',Figure2, 
       dpi = 600, width = 8, height = 3, units = 'in')

# Loadings plot for metabolomics data
top_metabolites <- data.frame(metabolite = metabolites, PCA.results$rotation) %>% 
  select(BIOCHEMICAL, PC1, PC2) %>% 
  mutate(weight = sqrt(PC1^2 + PC2^2)) %>% 
  # pivot_longer(-BIOCHEMICAL, names_to = "PC", values_to = "value") %>% 
  # group_by(PC) %>% 
  arrange(desc(abs(weight))) %>% 
  slice_head(n = 10)

loadings = data.frame(metabolite = metabolites, PCA.results$rotation) %>%
  select(BIOCHEMICAL, PC1, PC2) %>% 
  inner_join(top_metabolites)

# For ppt figure
loadings = loadings %>% 
  filter(BIOCHEMICAL == "uracil" | BIOCHEMICAL == "phosphate" | BIOCHEMICAL == "2'-deoxyuridine" | BIOCHEMICAL == "erythritol" | BIOCHEMICAL == "ethylmalonate")

label.locations = data.frame(BIOCHEMICAL = c('uracil','ethylmalonate','erythritol',"2'-deoxyuridine",'phosphate'),
                             PC1.label = c(0.25, 3.5, 0.2*11.5 + 0.5, 1.5, -0.15), 
                             PC2.label = c(2.25+0.1, -0.25, -0.06*11.5-0.1, -3, -3.5))

loadings = loadings %>% left_join(label.locations)

scale_factor = 11.5
SupplementalFigure3B = ggplot(PCA.data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Drug, shape = `Time point`), cex = 2.5) + 
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values = v_colors) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  theme_bw() + 
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*scale_factor, yend = PC2*scale_factor), 
               arrow = arrow(length = unit(0.03, "npc"))) + 
  #geom_text_repel(data = loadings, aes(x = PC1*scale_factor + 0.25, y = PC2*scale_factor + 0.25, label = BIOCHEMICAL)) + 
  geom_text(data = loadings, aes(x = PC1.label, y = PC2.label, label = BIOCHEMICAL)) + 
  ggtitle("Metabolomics") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("PC1: 45% variance") + 
  ylab("PC2: 13% variance") + 
  xlim(-3,4.5) + ylim(-3.5,2.5)
SupplementalFigure3B

# Combine supplemental figure3A and 3B into one figure
gridlayout <- data.matrix(c(1,NA,2)) %>% t()
SupplementalFigure3 <- arrangeGrob(SupplementalFigure3A, SupplementalFigure3B,
                       ncol = 3, 
                       nrow = 1, 
                       #heights = relative_heights,
                       widths = c(0.36+0.035,0.05,0.54+0.015),
                       layout_matrix = gridlayout) 

ggsave('results/figures/SupplementalFigure3.png', SupplementalFigure3, 
       dpi = 600, width = 8*1.25, height = 3*1.25, units = 'in')


## Getting table information for Figure 2B
# Load transcriptomics data for DESeq2 comparison
path.directory <- 'data/RNA-seq/'

# information needed: gene name
gene.expression.data <- read.csv(file = paste0(path.directory, "dougherty_rno_cardio_t6_ace_gene_deseq2.csv"), 
                                   header = TRUE) %>% 
  mutate(group = "Ace_6h", Drug = "Ace", Timepoint = '6h') %>% 
  rbind(read.csv(file = paste0(path.directory, "dougherty_rno_cardio_t24_ace_gene_deseq2.csv"), 
                   header = TRUE) %>% 
          mutate(group = "Ace_24h", Drug = "Ace", Timepoint = '24h')) %>% 
  rbind(read.csv(file = paste0(path.directory, "dougherty_rno_cardio_t6_dox_gene_deseq2.csv"), 
                   header = TRUE) %>% 
          mutate(group = "Dox_6h", Drug = "Dox", Timepoint = '6h')) %>% 
  rbind(read.csv(file = paste0(path.directory, "dougherty_rno_cardio_t24_dox_gene_deseq2.csv"), 
                   header = TRUE) %>% 
          mutate(group = "Dox_24h", Drug = "Dox", Timepoint = '24h')) %>% 
  rbind(read.csv(file = paste0(path.directory, "dougherty_rno_cardio_t6_5fu_gene_deseq2.csv"), 
                   header = TRUE) %>% 
          mutate(group = "5FU_6h", Drug = "5FU", Timepoint = '6h')) %>% 
  rbind(read.csv(file = paste0(path.directory, "dougherty_rno_cardio_t24_5fu_gene_deseq2.csv"), 
                   header = TRUE) %>% 
          mutate(group = "5FU_24h", Drug = "5FU", Timepoint = '24h'))


# Summarize results for Figure 2B
gene.expression.data %>% filter(fdr < 0.1) %>% group_by(group) %>% dplyr::count()


# Read in metabolomics data
metabolomics.data <- read.xlsx(file = "results/supplement/Supplemental File 3.xlsx", 
                            sheetName = "differential")
metabolomics.data %>% filter(p.adj < 0.1) %>% group_by(group) %>% dplyr::count()




#### Figure 3 ####
# Pathway and gene set enrichment analysis for transcriptomics data
# Identifying shared metabolites in the metabolomics data
# Include supplementary figure for metabolomics data showing all changes, broken up by consumed, ambiguous, produced


# Pathway enrichment
# Approach adapted from the clusterProfiler vignette
# Use BioMart to convert between gene symbols and gene IDs
mart  <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                          dataset = "rnorvegicus_gene_ensembl")
geneDict <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                          "external_gene_name", 
                                          "entrezgene_id",
                                          "description"), 
                           mart = mart)

# generate list of genes that are below a FDR threshold for enrichment analysis
geneList <- gene.expression.data %>%
  dplyr::select(EntrezID, fdr, logfc, group)
geneList$group = factor(geneList$group)

# Download pathways for analysis
# H - Hallmark gene sets
# C5 - GO gene sets
m_df = msigdbr(species = "Rattus norvegicus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

enrichment.results <- data.frame()
KEGG.enrichment.results <- data.frame()
set.seed(123)
for(condition in levels(geneList$group)){

  # Enrichment using Hallmark gene sets - GOOD, 47 significant results
  # Enrichment using C5 - GOOD, 3462 results that were significant
  data = geneList %>%
    filter(fdr < 0.01) %>%
    filter(group == condition)
  gene.enrichment <- enricher(data$EntrezID,
                              TERM2GENE = m_df,
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.1,
                              minGSSize = 10)
  
  enrichment.results <- enrichment.results %>% 
    rbind(data.frame(gene.enrichment) %>% mutate(condition = condition))
  
  gene.enrichment = enrichKEGG(data$EntrezID,
                               organism = "rno",
                               keyType = "kegg",
                               pvalueCutoff = 0.1,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 500,
                               qvalueCutoff = 0.2,
                               use_internal_data = FALSE
                               )
  
  KEGG.enrichment.results <- KEGG.enrichment.results %>% 
    rbind(data.frame(gene.enrichment) %>% mutate(condition = condition)) 
}


# Generate a figure for the Hallmark gene sets for this data
# facet_wrap by timepoint
plot.data <- enrichment.results %>% 
  dplyr::select(ID, pvalue,  condition) %>% 
  pivot_wider(names_from = "condition", values_from = "pvalue", values_fill = 0) %>% 
  pivot_longer(-ID, names_to = "group", values_to = "pvalue") %>% 
  separate(ID, sep = "HALLMARK_", into = c(NA,"ID"), remove = TRUE) %>% 
  mutate(ID = gsub(ID, pattern = "_", replacement = " ")) %>% 
  mutate(pvalue = ifelse(pvalue == 0, 0, 1))

condition.metadata <- data.frame(group = c("Ace_6h", "Ace_24h", "DMSO1_6h", "DMSO1_24h", "DMSO2_6h", "DMSO2_24h", "Dox_6h","Dox_24h", "5FU_6h", "5FU_24h"), 
                                 Timepoint = c('6h','24h','6h','24h','6h','24h','6h','24h','6h','24h'), 
                                 Drug = c("Ace","Ace","DMSO1","DMSO1","DMSO2","DMSO2","Dox","Dox","5FU","5FU"))

plot.data = plot.data %>% 
  left_join(condition.metadata) 

plot.data$Timepoint <- factor(plot.data$Timepoint, levels = c("6h","24h"))

# Figure 3A 
Figure3A = ggplot(plot.data, aes(x = Drug, y = fct_rev(ID))) + 
  geom_tile(aes(fill = pvalue)) + 
  facet_wrap(~Timepoint) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + 
  xlab("") + 
  ylab("") + 
  theme(legend.position = "none")
Figure3A


# Identify metabolites that were consistently produced with respect to BLANK and increased with respect to control conditions
load('data/metabolomics/blank_metabolomics.Rdata')
direction <- blank.comparisons %>% 
  filter(p.adj< 0.1) %>% 
  group_by(metabolite) %>% 
  summarize(consumed = length(which(diff < 0)), 
            produced = length(which(diff > 0))) %>% 
  mutate(direction = ifelse(produced == 0 & consumed != 0, "Consumed", 
                            ifelse(consumed == 0 & produced != 0, "Produced", "Ambiguous"))) %>% 
  right_join(data.frame(metabolite = levels(blank.comparisons$metabolite))) %>% 
  mutate(direction = ifelse(is.na(direction), "No change", direction)) %>% 
  select(metabolite, direction)

# Identify me
metabolomics.data <- read.xlsx(file = "results/supplement/Supplemental File 3.xlsx", 
                               sheetIndex = "differential")

shared.metabolites = metabolomics.data %>% 
  filter(direction != "Consumed" & direction != "Ambiguous") %>% 
  group_by(metabolite, Timepoint) %>% 
  summarize(sig = length(which(p.adj < 0.1))) %>% 
  filter(sig == 3) %>% 
  dplyr::select(metabolite) %>% 
  filter(metabolite != "HEPES") %>% 
  unique()

plot.data = metabolomics.data %>% 
  inner_join(shared.metabolites, by = c("metabolite")) %>% 
  mutate(`Difference in\nnormalized values` = ifelse(p.adj < 0.1, diff, 0))

plot.data$Timepoint <- factor(plot.data$Timepoint, levels = c("6h","24h"))

Figure3B = ggplot(plot.data, aes(x = Drug, y = fct_rev(metabolite))) + 
  geom_tile(aes(fill = `Difference in\nnormalized values`)) + 
  facet_wrap(~Timepoint) + 
  theme_bw() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", breaks = c(1.0, 0.75, 0.5, 0.25, 0, -0.25, -0.5)) + 
  xlab("") + ylab("") + 
  theme(legend.title.align = 0.5, 
        legend.direction = "vertical", 
        legend.box.just = "center")
Figure3B

# Generate Figure3
gridlayout <- matrix(c(1,NA,2), nrow = 1)
Figure3 <- arrangeGrob(Figure3A, Figure3B,
                       ncol = 3, nrow = 1,
                       #heights = c(0.475, 0.05, 0.475),
                       widths = c(0.475, 0.05, 0.475),
                       layout_matrix = gridlayout) %>% 
  as_ggplot() + 
  draw_text(text = c("A","B"), 
            size = 24,
            hjust = 0,
            x = c(0.01, 0.5), 
            y = c(0.98, 0.98))

ggsave('results/figures/Figure3.png', Figure3, 
       dpi = 800, width = 16*0.75, height = 6, units = 'in')


############# Supplemental Figure 4 for metabolomics data ####################
# Make comparisons to blank media for production/consumption of metabolites
# but it needs to be a two-level comparison

# Panel A - metabolites that were measured to be consumed

# Panel B - metabolites that were measured to be ambiguous

# Panel C - all metabolites split up by produced, consumed, no change

# Comparisons with metabolomics data:
# production/consumption with respect to blank media
# production/consumption with respect to control group
control.table = data.frame(Drug = c('5FU','Ace','Dox'), 
                           Control = c('DMSO1','DMSO2','DMSO2'))

load("data/metabolomics/blank_metabolomics.Rdata")
blank.comparisons = blank.comparisons %>% 
  inner_join(direction) %>% 
  mutate(direction = ifelse(p.adj > 0.1, "No change", 
                            ifelse(diff > 0, "Produced", "Consumed")))

plot.data = metabolomics.data %>% 
  filter(direction == "Consumed" | direction == "No") %>% 
  select(-direction) %>% 
  left_join(control.table) %>% 
  mutate(trt_group = paste0(Drug, '_', Timepoint), 
         ctl_group = paste0(Control, '_', Timepoint)) %>% 
  left_join(blank.comparisons %>% select(metabolite, direction, comparison), 
            by = c('trt_group' = 'comparison', "metabolite")) %>% 
  left_join(blank.comparisons %>% select(metabolite, direction, comparison), 
            by = c("ctl_group"="comparison", "metabolite")) %>% 
  mutate(`Difference in\n normalized values` = NA) %>% 
  mutate(`Difference in\n normalized values` = ifelse(p.adj < 0.1, diff, 
                                                     ifelse(p.adj > 0.1 & (direction.x != "No change" | direction.y != "No change"), 0, NA)))

plot.data$Timepoint <- factor(plot.data$Timepoint, levels = c("6h","24h"))

# metabolites that are blank across all of the conditions -> same consumption across conditions

# function to align legend, taken from: https://stackoverflow.com/questions/48000292/center-align-legend-title-and-legend-keys-in-ggplot2-for-long-legend-titles/48011882#48011882
align_legend <- function(p, hjust = 0.5){
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

SupplementalFigure4A = ggplot(plot.data, aes(x = Timepoint, y = fct_rev(metabolite))) + 
  geom_tile(aes(fill = -`Difference in\n normalized values`, color = "")) + 
  facet_wrap(~Drug) + 
  theme_bw() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                        na.value = "grey") +
  scale_color_manual(values = NA) + 
  xlab("") + ylab("") + 
  ggtitle("Consumption of metabolites across conditions") + 
  labs(fill = "Consumption of metabolite between\ntreated and control condition") +
  theme(legend.title.align = 0.5, 
        legend.direction = "vertical", 
        legend.box.just = "center", 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "left") + 
  guides(fill = guide_colorbar(order = 1), 
         color = guide_legend("Not consumed in either\ntreated or control condition", 
                             override.aes=list(fill="grey"), 
                             order = 2)) + 
  scale_y_discrete(position = "right")
ggdraw(align_legend(SupplementalFigure4A))

ggsave('results/figures/SupplementalFigure4A.png', ggdraw(align_legend(SupplementalFigure4A)), 
       dpi = 800, width = 8, height = 6)

# Insert subpanels for examples of metabolites that were (a) consumed across all conditions or (b) only consumed across one condition
metadata.master <- read.xlsx(file = "data/metabolomics/metabolomics_metadata.xlsx", sheetIndex = 1) %>% as.data.frame()
metadata.master$sample.ID <- as.character(metadata.master$sample.ID)
metadata.master$Drug = as.character(metadata.master$Drug)

# Load in raw metabolomics data for PCA plot
scaled.data = xlsx::read.xlsx(file = "results/supplement/Supplemental File 3.xlsx", 
                              sheetName = "scaled", header = TRUE) %>% 
  pivot_longer(-BIOCHEMICAL, names_to = "sample.ID", values_to = "value") %>% 
  left_join(metadata.master, by = c("sample.ID")) %>% 
  select(-sample.ID)

# (a) consumed across all conditions - adenine
# (b) no change in plot but was consumed in treated condition - heptanoate (7:0)
# (c) differential consumption - kynurenine
data.metadata = data.frame(condition = c("Blank", "DMSO1_6h", "5FU_6h", "DMSO2_6h", "Ace_6h", "Dox_6h", 
                                         "DMSO1_24h","5FU_24h","DMSO2_24h","Ace_24h","Dox_24h"), 
                           plot.name = c("Blank", "DMSO1, 6h", "5FU, 6h", "DMSO2, 6h", "Ace, 6h", "Dox, 6h", 
                                         "DMSO1, 24h","5FU, 24h","DMSO2, 24h","Ace, 24h","Dox, 24h"))
data.metadata$plot.name = factor(data.metadata$plot.name, 
                                 levels = c("Blank", "DMSO1, 6h", "5FU, 6h", "DMSO2, 6h", "Ace, 6h", "Dox, 6h", 
                                            "DMSO1, 24h","5FU, 24h","DMSO2, 24h","Ace, 24h","Dox, 24h"))
plot.data = scaled.data %>% 
  filter(BIOCHEMICAL %in% c("adenine","heptanoate (7:0)","kynurenine","glucose")) %>% 
  filter(!(group == "Blank" & Timepoint == "24h") & !(group == "Blank" & Timepoint == "6h")) %>% 
  dplyr::left_join(blank.comparisons, by = c("BIOCHEMICAL"="metabolite", "group"="comparison")) %>% 
  mutate(significance = ifelse(is.na(p.adj), 10, 
                               ifelse(p.adj < 0.1, 1, 0))) %>% 
  inner_join(data.metadata, by = c("group"="condition"))
plot.data$significance = factor(plot.data$significance, levels = c(0,1,10))
plot.data$group = factor(plot.data$group, 
                         levels = c("Blank","DMSO1_6h","5FU_6h",
                                    "DMSO1_24h","5FU_24h",
                                    "Ace_6h","Dox_6h","DMSO2_6h", 
                                    "Ace_24h","Dox_24h","DMSO2_24h")) 
  

blank = plot.data %>% filter(group == "Blank") %>% group_by(BIOCHEMICAL) %>% summarize(mean(value))

# Color dots by blank or significant change from blank
SupplementalFigure4A_a = ggplot(plot.data %>% filter(BIOCHEMICAL == "adenine"), 
                              aes(x = plot.name, y = value, color = significance)) + 
  geom_jitter(width = 0.1) + 
  stat_compare_means(comparisons = list(c('Dox, 24h','DMSO2, 24h'), 
                                        c('Ace, 24h','DMSO2, 24h'), 
                                        c('5FU, 24h','DMSO1, 24h'), 
                                        c('Ace, 6h','DMSO2, 6h')), 
                     label.y = c(0.6, 0.4, 0.55, 0.55)) + 
  scale_color_manual(values = c("grey","black","red")) +
  xlab("") + ylab("") +  
  theme_bw() + 
  geom_hline(aes(yintercept = 0.307), color = "red", size = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.9), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none") + 
  ggtitle("Consumption of adenine") + 
  ylim(-0.7,0.85)
SupplementalFigure4A_a

# Color dots by blank or significant change from blank
SupplementalFigure4A_b = ggplot(plot.data %>% filter(BIOCHEMICAL == "heptanoate (7:0)"), 
                              aes(x = plot.name, y = value, color = significance)) + 
  geom_jitter(width = 0.1) + 
  scale_color_manual(values = c("grey","black","red")) +
  xlab("") + ylab("") +  
  theme_bw() + 
  geom_hline(aes(yintercept = 0.186), color = "red", size = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.9), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none") + 
  ggtitle("Consumption of heptanoate (7:0)") + 
  ylim(-0.5,0.75)
SupplementalFigure4A_b

# Color dots by blank or significant change from blank
SupplementalFigure4A_c = ggplot(plot.data %>% filter(BIOCHEMICAL == "kynurenine"), 
                              aes(x = plot.name, y = value, color = significance)) + 
  geom_jitter(width = 0.1) + 
  stat_compare_means(comparisons = list(c('Ace, 6h','DMSO2, 6h'), 
                                        c('Ace, 24h','DMSO2, 24h'))) + 
  scale_color_manual(values = c("grey","black","red")) +
  xlab("") + ylab("") +  
  theme_bw() + 
  geom_hline(aes(yintercept = 0.0632), color = "red", size = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.9), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none") + 
  ggtitle("Consumption of kynurenine") + 
  ylim(-0.7,0.85)
SupplementalFigure4A_c

ggsave('results/figures/SupplementalFigure4A_a.png', SupplementalFigure4A_a, 
       dpi = 800, width = 3*1.1, height = 2.5*1.1)

ggsave('results/figures/SupplementalFigure4A_b.png', SupplementalFigure4A_b, 
       dpi = 800, width = 3*1.1, height = 2.5*1.1)

ggsave('results/figures/SupplementalFigure4A_c.png', SupplementalFigure4A_c, 
       dpi = 800, width = 3*1.1, height = 2.5*1.1)

# Supplemental figure for NO CHANGE in glucose
SupplementalFigure4C_a = ggplot(plot.data %>% filter(BIOCHEMICAL == "glucose"), 
                                aes(x = plot.name, y = value, color = significance)) + 
  geom_jitter(width = 0.1) + 
  stat_compare_means(comparisons = list(c('5FU, 6h','DMSO1, 6h'), 
                                        c('Dox, 24h','DMSO2, 24h'))) + 
  scale_color_manual(values = c("grey","red")) +
  xlab("") + ylab("") +  
  theme_bw() + 
  geom_hline(aes(yintercept = -0.00723), color = "red", size = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.9), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none") + 
  ggtitle("Changes in glucose") + 
  ylim(-0.7,0.85)
SupplementalFigure4C_a

ggsave('results/figures/SupplementalFigure4C_a.png', SupplementalFigure4C_a, 
       dpi = 800, width = 3*1.1, height = 2.5*1.1)

# Ambiguous metabolites - go back to raw data and plot with blank and control media
ambiguous.metabolites = metabolomics.data %>% 
  filter(direction == "Ambiguous") %>% 
  dplyr::select(metabolite) %>% 
  unique() %>% 
  filter(metabolite != "penicillin G")

load("data/metabolomics/blank_metabolomics.Rdata")
blank.comparisons = blank.comparisons %>% 
  inner_join(ambiguous.metabolites, by = c("metabolite")) %>% 
  dplyr::select(-group)
colnames(blank.comparisons)[1] = "group"
blank.comparisons$metabolite = factor(blank.comparisons$metabolite)
scaled.data$BIOCHEMICAL = factor(scaled.data$BIOCHEMICAL)
ambiguous.metabolites$metabolite = factor(ambiguous.metabolites$metabolite)

data <- scaled.data %>% 
  inner_join(ambiguous.metabolites, by = c("BIOCHEMICAL"="metabolite")) %>% 
  left_join(blank.comparisons %>% dplyr::select(group, metabolite, p.adj), by = c("group", "BIOCHEMICAL"="metabolite")) %>% 
  mutate(significance = ifelse(p.adj < 0.1, 1, 0)) %>% 
  mutate(significance = ifelse(Drug == "Blank", 10, significance)) %>% 
  filter(!(group == "Blank" & Timepoint == "6h")) %>% 
  filter(!(group == "Blank" & Timepoint == "24h")) 

colnames(data)[3] = "condition"
data$condition <- factor(data$condition)

blank.means = scaled.data %>% 
  inner_join(ambiguous.metabolites, by = c("BIOCHEMICAL"="metabolite")) %>% 
  filter(Drug == "Blank") %>% 
  filter(Timepoint == "0h") %>%  
  group_by(BIOCHEMICAL) %>% 
  summarize(blank.mean = mean(value))

data = data %>% left_join(blank.means)

data$condition = factor(data$condition, 
                        levels = c("Blank", "DMSO1_6h", "5FU_6h", "DMSO2_6h", "Ace_6h", "Dox_6h", 
                                   "DMSO1_24h","5FU_24h","DMSO2_24h","Ace_24h","Dox_24h"))
data.metadata = data.frame(condition = c("Blank", "DMSO1_6h", "5FU_6h", "DMSO2_6h", "Ace_6h", "Dox_6h", 
                                         "DMSO1_24h","5FU_24h","DMSO2_24h","Ace_24h","Dox_24h"), 
                           plot.name = c("Blank", "DMSO1, 6h", "5FU, 6h", "DMSO2, 6h", "Ace, 6h", "Dox, 6h", 
                                         "DMSO1, 24h","5FU, 24h","DMSO2, 24h","Ace, 24h","Dox, 24h"))
data.metadata$plot.name = factor(data.metadata$plot.name, 
                                 levels = c("Blank", "DMSO1, 6h", "5FU, 6h", "DMSO2, 6h", "Ace, 6h", "Dox, 6h", 
                                            "DMSO1, 24h","5FU, 24h","DMSO2, 24h","Ace, 24h","Dox, 24h"))

data = data %>% left_join(data.metadata)
data$significance = factor(data$significance, levels = c(0,1,10))

# Ambiguous metabolites
SupplementalFigure4B = ggplot(data) + 
  geom_jitter(aes(x = plot.name, y = value, color = significance)) + 
  geom_hline(aes(yintercept = blank.mean), color = "red", size = 1) +
  scale_color_manual(values=c("grey","black","red")) +
  facet_wrap(~BIOCHEMICAL) + 
  theme_bw() + 
  xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.9), 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        panel.spacing = unit(1, "lines")) + 
  ggtitle("Metabolites that were both produced and consumed across conditions")
SupplementalFigure4B

ggsave("results/figures/SupplementalFigure4B.png", SupplementalFigure4B, 
       dpi = 800, height = 6, width = 6)


# All metabolomics data
# read in differential metabolomics data
plot.data = metabolomics.data %>% 
  mutate(`Difference in\nnormalized values` = ifelse(p.adj < 0.1, diff, 0))

plot.data$Timepoint <- factor(plot.data$Timepoint, levels = c("6h","24h"))
plot.data$group = factor(plot.data$group, 
                         levels = c('5FU_6h','5FU_24h','Ace_6h','Ace_24h','Dox_6h','Dox_24h'))

# metabolites that are blank across all of the conditions -> same consumption across conditions
SupplementalFigure4C_produced = ggplot(plot.data %>% filter(direction == "Produced"), 
                              aes(x = Timepoint, y = fct_rev(metabolite))) + 
  geom_tile(aes(fill = `Difference in\nnormalized values`)) + 
  facet_wrap(~Drug, strip.position = "left", nrow = 3) + 
  theme_bw() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  xlab("") + ylab("") + 
  coord_flip() + 
  labs(fill = "Production of metabolite between\ntreated and control condition") +
  theme(legend.title.align = 0.5, 
        legend.direction = "vertical", 
        legend.box.just = "center", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
  ggtitle("Metabolites produced with respect to blank media")
ggdraw(align_legend(SupplementalFigure4C_produced))

ggsave('results/figures/SupplementalFigure4C_produced.png', 
       ggdraw(align_legend(SupplementalFigure4C_produced)), 
       dpi = 800, width = 14, height = 5)

SupplementalFigure4C_noChange = ggplot(plot.data %>% filter(direction == "No change"), 
                                       aes(x = Timepoint, y = (metabolite))) + 
  geom_tile(aes(fill = `Difference in\nnormalized values`)) + 
  facet_wrap(~Drug, strip.position = "left", nrow = 3) + 
  theme_bw() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") + 
  xlab("") + ylab("") + 
  labs(fill = "Difference between\ntreated and control condition") +
  coord_flip() + 
  theme(legend.title.align = 0.5, 
        legend.direction = "vertical", 
        legend.box.just = "center", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "left") +
  #scale_y_discrete(position = "right") + 
  ggtitle("Metabolites with no change\ncompared to blank media")
ggdraw(align_legend(SupplementalFigure4C_noChange))

ggsave('results/figures/SupplementalFigure4C_noChange.png', 
       ggdraw(align_legend(SupplementalFigure4C_noChange)), 
       dpi = 800, width = 12*0.75, height = 5*0.9)

################## Figure 4 - TIDEs analysis #######################
# generate list of genes that are below a FDR threshold for enrichment analysis

# Read in actual TIDEs results 
path.directory = "data/TIDEs/"
# Compare to using an FDR < 0.01 and FDR < 0.1
TIDEs = read.xlsx(file = paste0(path.directory, "Ace6hrs.xlsx"), 
                        header = FALSE, sheetIndex = 1) %>% 
  dplyr::mutate(group = "Ace_6h", Drug = "Ace", Timepoint = '6h') %>% 
  rbind(read.xlsx(file = paste0(path.directory, "Ace24hrs.xlsx"), 
                  header = FALSE, sheetIndex = 1) %>% 
          mutate(group = "Ace_24h", Drug = "Ace", Timepoint = '24h')) %>% 
  rbind(read.xlsx(file = paste0(path.directory, "Dox6hrs.xlsx"), 
                  header = FALSE, sheetIndex = 1) %>% 
          mutate(group = "Dox_6h", Drug = "Dox", Timepoint = '6h')) %>% 
  rbind(read.xlsx(file = paste0(path.directory, "Dox24hrs.xlsx"), 
                  header = FALSE, sheetIndex = 1) %>% 
          mutate(group = "Dox_24h", Drug = "Dox", Timepoint = '24h')) %>% 
  rbind(read.xlsx(file = paste0(path.directory, "FiveFU6hrs.xlsx"), 
                  header = FALSE, sheetIndex = 1) %>% 
          mutate(group = "5FU_6h", Drug = "5FU", Timepoint = '6h')) %>% 
  rbind(read.xlsx(file = paste0(path.directory, "FiveFU24hrs.xlsx"), 
                  header = FALSE, sheetIndex = 1) %>% 
          mutate(group = "5FU_24h", Drug = "5FU", Timepoint = '24h'))
colnames(TIDEs)[1:4] = c('ID','description','score','significance')

# Load in the number of reactions per task, filter for tasks that have > 3 reactions
TIDEs_df = read.csv(file = "data/TIDEs/TIDEs_rxns.csv") %>% 
  group_by(ID) %>% dplyr::count() %>% 
  filter(n > 2) %>% 
  select(-n)


# The number of significant TIDEs for each group
TIDEs %>% inner_join(TIDEs_df) %>% filter(abs(significance) <= 0.05) %>% group_by(group) %>% dplyr::count()

# The task that is most commonly different
TIDEs %>% 
  inner_join(TIDEs_df) %>% 
  filter(abs(significance) <= 0.05) %>% 
  group_by(ID, description) %>% 
  dplyr::count() %>% 
  View()

TIDEs %>% 
  mutate(significance = ifelse(abs(significance) < 0.05, 
                               ifelse(significance > 0, 1, -1), 0)) %>% 
  select(ID, description, significance, group) %>% 
  pivot_wider(names_from = "group", values_from = "significance") %>% 
  View()

# Identify the tasks that are commonly different between Dox and 5FU but not Ace
TIDEs %>% 
  filter(Drug != "Ace") %>% 
  filter(abs(significance) < 0.025) %>% 
  group_by(ID, description) %>% 
  dplyr::count() %>% 
  View()

# Plot the TIDEs categories, ordered by type of task
# Load reaction subcategories
task.categories <- read.xlsx(file = "data/TIDEs/TIDEs_categories.xlsx", 
                             sheetIndex = 2, 
                             header = 1)
task.categories$category <- factor(task.categories$category, 
                                   levels = levels(factor(task.categories$category)))
task.categories$plot.name = factor(task.categories$plot.name, 
                                   levels = task.categories$plot.name)

# Common tides
common.TIDEs = TIDEs %>% 
  filter(abs(significance) < 0.05) %>% 
  group_by(ID, description) %>% dplyr::count() %>% 
  filter(n > 0) 

# Interesting tasks to plot:
# A: changes in DNA synthesis
# B: changes in central carbon metabolism (include discussion of glucose)
# C: changes in lipid synthesis for Ace treatment
# D: shared changes across all conditions

# All the tasks are included in a task category
Figure4.data = TIDEs %>% 
  filter(!str_detect(ID, "S")) %>% 
  inner_join(task.categories) %>% 
  inner_join(common.TIDEs %>% dplyr::select(ID)) %>% 
  mutate(sig.label = ifelse(abs(significance) < 0.05, 
                            ifelse(abs(significance*2) < 0.01, '**', '*'), '')) %>% 
  mutate(significance = ifelse(abs(significance) < 0.05, 
                               ifelse(significance < 0, -1, 1), 0)) %>% 
  mutate(label.color = ifelse(significance < 0, "white", "black"))

Figure4.data$Timepoint = factor(Figure4.data$Timepoint, levels = c("6h","24h"))


# RNA/DNA synthesis and nucleotide synthesis
Figure4A = ggplot(Figure4.data %>% filter(category == "DNA/RNA synthesis" | category == "Nucleotide metabolism"), 
                  aes(x = (plot.name), y = Timepoint, fill = significance)) + 
  geom_tile(width = 0.9, height = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  #geom_text(aes(label = sig.label, color = label.color), size = 5, fontface = "bold") + 
  #scale_colour_manual(values=c("black", "white")) +
  facet_wrap(~Drug) + 
  theme_bw() + 
  ggtitle("DNA/RNA synthesis and Nucleotide metabolism") + 
  coord_flip() + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())
Figure4A

# oxidative phosphorylation
Figure4B = ggplot(Figure4.data %>% filter(category == "Central carbon metabolism"), 
                  aes(x = fct_rev(plot.name), y = Timepoint, fill = significance)) + 
  geom_tile(width = 0.9, height = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  #geom_text(aes(label = sig.label, color = label.color), size = 5, fontface = "bold") + 
  #scale_colour_manual(values=c("black", "white")) +
  facet_wrap(~Drug) + 
  theme_bw() + 
  ggtitle("Central carbon metabolism") + 
  coord_flip() + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())
Figure4B

# Phospholipid metabolism
Figure4C = ggplot(Figure4.data %>% filter(category == "Lipid synthesis"), 
                  aes(x = fct_rev(plot.name), y = Timepoint, fill = significance)) + 
  geom_tile(width = 0.9, height = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  #geom_text(aes(label = sig.label, color = label.color), size = 5, fontface = "bold") + 
  #scale_colour_manual(values=c("black", "white")) +
  facet_wrap(~Drug) + 
  theme_bw() + 
  ggtitle("Lipid synthesis") + 
  coord_flip() + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())
Figure4C

# ROS and NO synthesis
Figure4D = ggplot(Figure4.data %>% filter(category == "Signaling metabolism") %>% filter(ID != "C51"), 
                  aes(x = fct_rev(plot.name), y = Timepoint, fill = significance)) + 
  geom_tile(width = 0.9, height = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  #geom_text(aes(label = sig.label, color = label.color), size = 5, fontface = "bold") + 
  #scale_colour_manual(values=c("black", "white")) +
  facet_wrap(~Drug) + 
  theme_bw() + 
  ggtitle("Signaling metabolism") + 
  coord_flip() + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())
Figure4D

# Produce the final figure
gridlayout <- matrix(c(1,1,1,1,1,
                       NA, NA, NA, NA, NA,
                       2, NA, 3, NA, NA,
                       2, NA, 3, NA, 4), nrow = 5)
gridlayout
a = 0.098
Figure4 <- arrangeGrob(Figure4A, Figure4B, Figure4C, Figure4D,
                       ncol = 4, nrow = 5,
                       heights = c(0.32+0.045,0.01,0.37+0.045,0.01,0.2),
                       widths = c(0.45, 0.05, a, 0.5-a),
                       layout_matrix = gridlayout) %>% 
  as_ggplot() + 
  draw_text(text = c("A","B","C", "D"), 
            size = 24,
            hjust = 0,
            x = c(0.01, 0.55, 0.55, 0.55), 
            y = c(0.98, 0.98, 0.6, 0.175))

ggsave('results/figures/Figure4.png', Figure4, 
       dpi = 800, width = 12, height = 7, units = 'in')


# Supplemental figure
# Plot random data distributions for two tasks: ROS detoxification and NOS production
# IDs: C50 and C28
# Compare to using an FDR < 0.01 and FDR < 0.1
# Read in actual TIDEs results 
IDs = xlsx::read.xlsx(file = paste0(path.directory, 'Ace6hrs.xlsx'), sheetIndex = 1, header = FALSE) %>%
  dplyr::rename("ID" = "X1", "description" = "X2") %>% 
  dplyr::select(ID, description)

random.data = cbind(IDs, readxl::read_xlsx(path = paste0(path.directory, "Ace6hrs_random.xlsx"), 
                        col_names = FALSE, sheet = 1)) %>% 
  dplyr::mutate(group = "Ace_6h", Drug = "Ace", Timepoint = '6h') %>% 
  filter(ID == "C50" | ID == "C28") %>% 
  rbind(cbind(IDs, readxl::read_xlsx(path = paste0(path.directory, "Ace24hrs_random.xlsx"), 
                                     col_names = FALSE, sheet = 1)) %>% 
          dplyr::mutate(group = "Ace_24h", Drug = "Ace", Timepoint = '24h') %>% 
          filter(ID == "C50" | ID == "C28")) %>% 
  rbind(cbind(IDs, readxl::read_xlsx(path = paste0(path.directory, "Dox6hrs_random.xlsx"), 
                                     col_names = FALSE, sheet = 1)) %>% 
          dplyr::mutate(group = "Dox_6h", Drug = "Dox", Timepoint = '6h') %>% 
          filter(ID == "C50" | ID == "C28")) %>% 
  rbind(cbind(IDs, readxl::read_xlsx(path = paste0(path.directory, "Dox24hrs_random.xlsx"), 
                                     col_names = FALSE, sheet = 1)) %>% 
          dplyr::mutate(group = "Dox_24h", Drug = "Dox", Timepoint = '24h') %>% 
          filter(ID == "C50" | ID == "C28")) %>% 
  rbind(cbind(IDs, readxl::read_xlsx(path = paste0(path.directory, "5FU6hrs_random.xlsx"), 
                                     col_names = FALSE, sheet = 1)) %>% 
          dplyr::mutate(group = "5FU_6h", Drug = "5FU", Timepoint = '6h') %>% 
          filter(ID == "C50" | ID == "C28")) %>% 
  rbind(cbind(IDs, readxl::read_xlsx(path = paste0(path.directory, "5FU24hrs_random.xlsx"), 
                                     col_names = FALSE, sheet = 1)) %>% 
          dplyr::mutate(group = "5FU_24h", Drug = "5FU", Timepoint = '24h') %>% 
          filter(ID == "C50" | ID == "C28")) %>% 
  dplyr::select(-Drug, -Timepoint)
colnames(random.data)[3:1002] <- paste0('X',c(1:1000))

random.data = pivot_longer(random.data, cols = X1:X1000, names_to = "sample", values_to = "random.score")

random.data.metadata = data.frame(group = c('5FU_6h','5FU_24h','Ace_6h','Ace_24h','Dox_6h','Dox_24h'),
                                  plot.name = c('5FU, 6h','5FU, 24h','Ace, 6h','Ace, 24h','Dox, 6h','Dox, 24h'))

# Left join with the actual data
random.data = random.data %>% 
  left_join(TIDEs, by = c("ID","description","group")) %>% 
  mutate(significance = ifelse(abs(significance) < 0.05, 
                               ifelse(significance < 0, -1, 1), 0)) %>% 
  left_join(random.data.metadata, by = "group")

random.data$group = factor(random.data$plot.name, levels = c('5FU, 6h','5FU, 24h','Ace, 6h','Ace, 24h','Dox, 6h','Dox, 24h'))

annotations <- data.frame(random.data %>% select(description, group, score) %>% unique()) %>% 
  cbind(xpos = c(1.75,1.55,1.85,1.75,10,14,
                 -10,-12,-17,-12,-10,-6),
        ypos =  c(185,115,175,100,260,315,
                  185,215,575,325,325,215),
        hjustvar = c(0.95,0.95,1,1,1,1,
                     0,0,0,0,0,0),
        vjustvar = c(0.5,0.5,0.5,0.5,0.5,0.5,
                     0.5,0.5,0.5,0.5,0.5,0.5))


# random.data = random.data %>% left_join(random.data.metadata)
# random.data$plot.name = factor(random.data$plot.name, 
#                                levels = c('5FU, 6h','5FU, 24h','Ace, 6h','Ace, 24h','Dox, 6h','Dox, 24h'))

random.data$description = factor(random.data$description, 
                                 levels = c("Arginine to nitric oxide","ROS detoxification"))

SupplementalFigure6 = random.data %>% 
  ggplot() +
  geom_rect(aes(fill = significance), 
            xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, show.legend = FALSE) +
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  geom_histogram(aes(x = random.score), color = "black", fill = "black") + 
  geom_vline(aes(xintercept = score), color = "red", linetype = "dashed", size = 0.9) + 
  #geom_text(aes(x = -Inf, y = 0, label = round(score,2)), hjust = 1) + 
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=round(score,2))) +
  facet_wrap(~description + group, scales = "free", ncol = 6) +
  theme_bw() + ylab('Frequency') + xlab('Task score')
SupplementalFigure6

ggsave('results/figures/SupplementalFigure6.png', SupplementalFigure6, 
       dpi = 800, width = 10, height = 4, units = 'in')


############ Figure 5 - RIPTIDE analysis of condition-specific models #####################
# Read in the metabolomics blank data comparisons
load(file = "data/metabolomics/blank_metabolomics.Rdata")

# Load in mapping between model and metabolomics data
mapped.metabolites <- xlsx::read.xlsx(file = "data/metabolomics/mapped_metabolomics.xlsx", 
                                sheetIndex = 1, header = TRUE)

model.media = blank.comparisons %>% 
  inner_join(mapped.metabolites, by = c("metabolite"="Metabolon")) %>% 
  dplyr::select(comparison, model, reaction_id, diff, p.adj) %>% 
  mutate(direction = ifelse(p.adj < 0.1, 
                            ifelse(diff > 0, "Produced", "Consumed"), "No change")) %>% 
  dplyr::select(comparison, model, reaction_id, direction)
model.media$comparison = factor(model.media$comparison)


# Add in metabolites that were said to be in the media but not measured in metabolomics data
media = xlsx::read.xlsx(file = "data/metabolomics/Media Formulation.xlsx", 
                  sheetIndex = 1, startRow = 2, header = TRUE) %>% 
  dplyr::select(metabolite.name, transport.) %>% 
  dplyr::rename(reaction_id = "transport.") %>% 
  filter(!is.na(reaction_id)) %>% 
  filter(reaction_id != "not in model") %>% 
  filter(reaction_id != "-")

metabolites.not.detected = anti_join(media, model.media)

for(metabolite in metabolites.not.detected$reaction_id){
  model.media = model.media %>% 
    rbind(data.frame(comparison = levels(model.media$comparison), 
                     model = (metabolites.not.detected %>% filter(reaction_id == metabolite))$metabolite.name, 
                     reaction_id = metabolite, 
                     direction = "Able to consume"))
  
}


# Add in appropriate bounds based on direction of production/consumption
model.media %>% filter(comparison == "5FU_24h") %>% View()






# Adapted from https://github.com/csbl/Jenior_RIPTiDe_2020/blob/master/code/R/figure_5B.R
path.directory = "data/RIPTIDE/"
condition = c('dmso1_24h','fivefu_24h','dmso2_24h','ace_24h','dox_24h', 
              'dmso1_6h','fivefu_6h','dmso2_6h','ace_6h','dox_6h')


correlations = data.frame()
# For each fraction, read in the r2 value and p-value
for(current.condition in condition){
  current.file = paste0(path.directory, current.condition,'_media','/parameters.txt')
  correlations = readtext(current.file) %>% 
    separate_rows(text, sep = "\n") %>% 
    separate(text, sep = ":", into = c("text", "value")) %>% 
    filter(text == "Correlation between activity and transcriptome") %>% 
    separate(value, sep = ",", into = c("R", "pval")) %>% 
    separate(R, sep = "=", into = c("junk","R")) %>% dplyr::select(-junk) %>% 
    separate(pval, sep = "=", into = c("junk","pval")) %>% dplyr::select(-junk) %>% 
    mutate(condition = current.condition) %>% 
    dplyr::select(condition, R, pval) %>% 
    rbind(correlations)
}

data = data.frame()
set.seed(123)
sub.sample = sample(1:75, 75, replace = FALSE)
# Read in chosen data
for(current.condition in condition){
  current.file = paste0(path.directory, current.condition,'_media','/flux_samples.tsv')
  current.data = read.table(current.file, header = TRUE)
  
  data = rbind(data, 
               current.data[sub.sample,] %>% mutate(sample = 1:75) %>% 
                 pivot_longer(-sample, names_to = "reaction_id", values_to = "flux") %>%
                 mutate(condition = current.condition))
}

metadata.condition = data.frame(metadata = c('dmso1_24h','fivefu_24h','dmso2_24h','ace_24h','dox_24h',
                                             'dmso1_6h','fivefu_6h','dmso2_6h','ace_6h','dox_6h'),
                                comparison = c('DMSO1_24h','5FU_24h','DMSO2_24h','Ace_24h','Dox_24h',
                                               'DMSO1_6h','5FU_6h','DMSO2_6h','Ace_6h','Dox_6h'),
                                Drug = c('DMSO1','5FU','DMSO2','Ace','Dox',
                                         'DMSO1','5FU','DMSO2','Ace','Dox'), 
                                timepoint = c('24h','24h','24h','24h','24h',
                                              '6h','6h','6h','6h','6h'))

data = data %>% inner_join(metadata.condition, by = c("condition"="metadata"))

data$Drug = factor(data$Drug, levels = c("5FU","DMSO1","Ace","Dox","DMSO2"))
data$timepoint = factor(data$timepoint, levels = c("6h","24h"))

# Get number of reactions in each model
data %>% select(condition, reaction_id) %>% unique() %>% group_by(condition) %>% dplyr::count()
 
reaction.counts.6h = data %>% 
  filter(timepoint == "6h") %>% 
  select(Drug, reaction_id) %>% 
  unique() %>% 
  group_by(reaction_id) %>% 
  dplyr::count()

reaction.counts.24h = data %>% 
  filter(timepoint == "24h") %>% 
  select(Drug, reaction_id) %>% 
  unique() %>% 
  group_by(reaction_id) %>% 
  dplyr::count()

reactions.24h = data %>% filter(timepoint == "24h") %>% select(Drug, reaction_id) %>% unique()
reactions.6h = data %>% filter(timepoint == "6h") %>% select(Drug, reaction_id) %>% unique()

# Read in reaction names from model
reaction.names = xlsx::read.xlsx(file = "data/ncomms14250-s13, iRno RAVEN.xlsx", 
                           header = TRUE, sheetIndex = 1) %>% 
  dplyr::select(RXNID, NAME, EQUATION, SUBSYSTEM, GENE.ASSOCIATION, EC.NUMBER)



# Identify reactions that are unique between each pair of treatment and controls
unique.5FU.24h = reactions.24h %>% 
  filter(Drug == "5FU" | Drug == "DMSO1") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.24h) %>% 
  filter(Drug == "5FU") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.Ace.24h = reactions.24h %>% 
  filter(Drug == "Ace" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.24h) %>% 
  filter(Drug == "Ace") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.Dox.24h = reactions.24h %>% 
  filter(Drug == "Dox" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.24h) %>% 
  filter(Drug == "Dox") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.DMSO1.24h = reactions.24h %>% 
  filter(Drug == "5FU" | Drug == "DMSO1") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.24h) %>% 
  filter(Drug == "DMSO1") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.DMSO2.Ace.24h = reactions.24h %>% 
  filter(Drug == "Ace" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.24h) %>% 
  filter(Drug == "DMSO2") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.DMSO2.Dox.24h = reactions.24h %>% 
  filter(Drug == "Dox" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.24h) %>% 
  filter(Drug == "DMSO2") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.drug = inner_join(unique.Ace.24h, unique.5FU.24h, by = "reaction_id") %>% inner_join(unique.Dox.24h, by = "reaction_id")
unique.control = inner_join(unique.DMSO1.24h, unique.DMSO2.Ace.24h, by = "reaction_id") %>% inner_join(unique.DMSO2.Dox.24h, by = "reaction_id")

# glutamate is unique to the drug treatment conditions - RCR20594
plot_data = data %>% 
  filter(reaction_id %in% c("RCR11194","RCR11165"))
plot_data$Drug = factor(plot_data$Drug, levels = c('5FU','DMSO1','Ace','Dox','DMSO2'))
ggplot(plot_data, aes(x = Drug, y = flux, color = reaction_id)) + 
  geom_jitter() + 
  theme_bw() +
  facet_wrap(~timepoint, scales = "free") + 
  stat_compare_means(comparisons = list(c('5FU','DMSO1'),
                                        c('Ace','DMSO2'),
                                        c('Dox','DMSO2')), 
                     # method.args = list(alternative = "greater"), 
                     method = "wilcox.test", 
                     paired = FALSE)

# 6 hour unique reactions
unique.5FU.6h = reactions.6h %>% 
  filter(Drug == "5FU" | Drug == "DMSO1") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.6h) %>% 
  filter(Drug == "5FU") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.Ace.6h = reactions.6h %>% 
  filter(Drug == "Ace" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.6h) %>% 
  filter(Drug == "Ace") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.Dox.6h = reactions.6h %>% 
  filter(Drug == "Dox" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.6h) %>% 
  filter(Drug == "Dox") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.DMSO1.6h = reactions.6h %>% 
  filter(Drug == "5FU" | Drug == "DMSO1") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.6h) %>% 
  filter(Drug == "DMSO1") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.DMSO2.Ace.6h = reactions.6h %>% 
  filter(Drug == "Ace" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.6h) %>% 
  filter(Drug == "DMSO2") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.DMSO2.Dox.6h = reactions.6h %>% 
  filter(Drug == "Dox" | Drug == "DMSO2") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>%   
  filter(n == 1) %>% 
  left_join(reactions.6h) %>% 
  filter(Drug == "DMSO2") %>% 
  left_join(reaction.names, by = c("reaction_id"="RXNID"))

unique.drug.6h = inner_join(unique.Ace.6h, unique.5FU.6h, by = "reaction_id") %>% inner_join(unique.Dox.6h, by = "reaction_id")
unique.control.6h = inner_join(unique.DMSO1.6h, unique.DMSO2.Ace.6h, by = "reaction_id") %>% inner_join(unique.DMSO2.Dox.6h, by = "reaction_id")

# Code from original RIPTIDE implementation of NMDS: 
library(vegan)
library(ape)

shared.reactions = data %>% 
  #filter(timepoint == "6h") %>% 
  group_by(reaction_id) %>% 
  dplyr::count() %>% 
  filter(n == 750) %>% 
  select(-n)

NMDS.data = inner_join(data , shared.reactions) %>% 
  pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
  select(-sample) %>% 
  select(-condition, -Drug, -timepoint, -comparison)
NMDS.metadata = inner_join(data , shared.reactions) %>% 
  pivot_wider(names_from = "reaction_id", values_from = "flux") %>%
  select(sample, condition)

# Account for negative flux values:
NMDS.data <- as.matrix(NMDS.data) + abs(min(NMDS.data))
flux_bray_dist <- vegdist(NMDS.data, method='bray') # Bray-Curtis
library(ape)

# Try NMDS to reduce to two dimensions
# Don't worry about not converging
nmds = metaMDS(NMDS.data, distance = "bray")

plot.data = data.frame(nmds$points) %>% 
  cbind(NMDS.metadata) %>% 
  inner_join(metadata.condition, by = c("condition"="metadata"))

plot.data$timepoint = factor(plot.data$timepoint, levels = c("6h","24h"))

v_colors = viridis(5, option = "C")
plot.data$Drug = factor(plot.data$Drug, levels = c("5FU", "DMSO1", "Ace", "Dox", "DMSO2"))
SupplementalFigure7A = ggplot(data = plot.data %>% filter(timepoint == "6h"), 
                              aes(x = MDS1, y = MDS2, fill = Drug, shape = timepoint)) +
  geom_point(cex = 1.75) +
  theme_bw() + 
  xlab("NMDS1") + ylab("NMDS2") + 
  scale_fill_manual(values = v_colors) + 
  scale_shape_manual(values=c(21)) +
  theme(legend.position = "none") + xlim(-0.015,0.011) + ylim(-0.005,0.006)
SupplementalFigure7A

SupplementalFigure7B = ggplot(data = plot.data %>% filter(timepoint == "24h"), 
                              aes(x = MDS1, y = MDS2, fill = Drug, shape = timepoint)) +
  geom_point(cex = 1.75) +
  theme_bw() + 
  xlab("NMDS1") + ylab("NMDS2") + 
  scale_fill_manual(values = v_colors) + 
  scale_shape_manual(values = c(21)) + 
  guides(shape = FALSE) + 
  guides(fill=guide_legend(override.aes=list(shape=21))) + 
  xlim(-0.015,0.011) + ylim(-0.005,0.006)
SupplementalFigure7B

# Combine figures into one figure
gridlayout <- data.matrix(c(1,NA,2)) %>% t()
SupplementalFigure7 <- arrangeGrob(SupplementalFigure7A, SupplementalFigure7B,
                       ncol = 3, 
                       nrow = 1, 
                       #heights = relative_heights,
                       widths = c(0.36+0.05,0.05,0.54),
                       layout_matrix = gridlayout) 

ggsave('results/figures/SupplementalFigure7.png', SupplementalFigure7, 
       dpi = 600, width = 8, height = 3, units = 'in')

# Implementing random forest for predictors between different conditions
# Adapted from: https://github.com/csbl/Jenior_RIPTiDe_2020/blob/master/code/R/figure_S4.R

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(1234)

# Format the data correctly for RF
comparisons = data.frame(Comparison = c("5FU_24h - DMSO1_24h","Ace_24h - DMSO2_24h","Dox_24h - DMSO2_24h",
                                        "5FU_6h - DMSO1_6h","Ace_6h - DMSO2_6h","Dox_6h - DMSO2_6h"), 
                         x = c("fivefu_24h","ace_24h","dox_24h", 
                               "fivefu_6h","ace_6h","dox_6h"), 
                         y = c("dmso1_24h","dmso2_24h","dmso2_24h", 
                               "dmso1_6h","dmso2_6h","dmso2_6h"))
data$condition = factor(data$condition)
comparisons$x = factor(comparisons$x, levels = levels(data$condition))
comparisons$y = factor(comparisons$y, levels = levels(data$condition))

rf.results = data.frame()
# Filter out reactions that are correlation == 1 in both groups
library(caret)

for(row in 1:nrow(comparisons)){
  x = data %>% 
    filter(condition == (comparisons[row,]$x)) %>% 
    pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
    select(-sample, -condition, -comparison, -Drug, -timepoint)
  
  # corr_x = cor((x))
  # highCorr = findCorrelation(corr_x, cutoff = 0.9)
  # x = x[-highCorr]
  
  
  y = data %>% 
    filter(condition == (comparisons[row,]$y)) %>% 
    pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
    select(-sample, - condition, -comparison, -Drug, -timepoint)
  
  # corr_y = cor((y))
  # highCorr = findCorrelation(corr_y, cutoff = 0.9)
  # y = y[-highCorr]

  # Identify shared reactions between the two data sets
  shared <- intersect(colnames(x), colnames(y))
  
  current.data = rbind(cbind(x[,shared], condition = rep(1,times = 75)),
                       cbind(y[,shared], condition = rep(0, times = 75))) 
  current.data$condition = factor(current.data$condition, levels = c(0,1))
  # k0 - number of variables to identify and then stop
  # pdel = 0 - remove only one variable with each iteration
  current.aucrf = AUCRF(condition ~ ., data = current.data, pdel = 0, k0 = 10, ranking = "MDA")
  
  # Create a feature table
  rf.results = rf.results %>% 
    rbind(data.frame(condition = comparisons[row,]$Comparison, 
                                   current.aucrf$ranking[1:10]) %>% 
            rownames_to_column("reaction_id"))
}

# All treatment at 24 or 6 hours vs control conditions
x = data %>% 
  filter(timepoint == "24h" & Drug %in% c('5FU','Ace','Dox')) %>% 
  pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
  select(-sample, -condition, -comparison, -Drug, -timepoint) %>% 
  select_if(~ !any(is.na(.)))

# corr_x = cor((x))
# highCorr = findCorrelation(corr_x, cutoff = 0.9)
# x = x[-highCorr]


y = data %>% 
  filter(timepoint == "24h" & Drug %in% c('DMSO1','DMSO2')) %>%  
  pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
  select(-sample, - condition, -comparison, -Drug, -timepoint) %>% 
  select_if(~ !any(is.na(.)))

# corr_y = cor((y))
# highCorr = findCorrelation(corr_y, cutoff = 0.9)
# y = y[-highCorr]

# Identify shared reactions between the two data sets
shared <- intersect(colnames(x), colnames(y))

current.data = rbind(cbind(x[,shared], condition = rep(1,times = 75)),
                     cbind(y[,shared], condition = rep(0, times = 75))) 
current.data$condition = factor(current.data$condition, levels = c(0,1))
# k0 - number of variables to identify and then stop
# pdel = 0 - remove only one variable with each iteration
current.aucrf = AUCRF(condition ~ ., data = current.data, pdel = 0, k0 = 10, ranking = "MDA")

# Create a feature table
rf.results = rf.results %>% 
  rbind(data.frame(condition = "treatment - control, 24 hrs", 
                   current.aucrf$ranking[1:10]) %>% 
          rownames_to_column("reaction_id"))

# All treatment at 24 or 6 hours vs control conditions
x = data %>% 
  filter(timepoint == "6h" & Drug %in% c('5FU','Ace','Dox')) %>% 
  pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
  select(-sample, -condition, -comparison, -Drug, -timepoint) %>% 
  select_if(~ !any(is.na(.)))

# corr_x = cor((x))
# highCorr = findCorrelation(corr_x, cutoff = 0.9)
# x = x[-highCorr]


y = data %>% 
  filter(timepoint == "6h" & Drug %in% c('DMSO1','DMSO2')) %>%  
  pivot_wider(names_from = "reaction_id", values_from = "flux") %>% 
  select(-sample, - condition, -comparison, -Drug, -timepoint) %>% 
  select_if(~ !any(is.na(.)))

# corr_y = cor((y))
# highCorr = findCorrelation(corr_y, cutoff = 0.9)
# y = y[-highCorr]

# Identify shared reactions between the two data sets
shared <- intersect(colnames(x), colnames(y))

current.data = rbind(cbind(x[,shared], condition = rep(1,times = 75)),
                     cbind(y[,shared], condition = rep(0, times = 75))) 
current.data$condition = factor(current.data$condition, levels = c(0,1))
# k0 - number of variables to identify and then stop
# pdel = 0 - remove only one variable with each iteration
current.aucrf = AUCRF(condition ~ ., data = current.data, pdel = 0, k0 = 18, ranking = "MDA")

# Create a feature table
rf.results = rf.results %>% 
  rbind(data.frame(condition = "treatment - control, 6 hrs", 
                   current.aucrf$ranking[1:10]) %>% 
          rownames_to_column("reaction_id"))

rf.unique.reactions = data %>% select(reaction_id, condition) %>% 
  unique() %>% 
  left_join(reaction.names, by = c('reaction_id'='RXNID')) %>% 
  filter(SUBSYSTEM == "Transport") %>% 
  select(reaction_id) %>% 
  unique() %>% pull(reaction_id)



# RCR10294 (connect back to transcript expression differences; 25630;24651;100364062)
# RCR11159 (connects with above; same GPR rules)
# for Ace, RCR10294 -> RCR10443 -> RCR10444 (no change) -> RCR11299 (no change) -> 
# RCR11300 (no change) -> RCR11301 (no change) -> RCR14581 (293620/no change) ->
# RCR13197 (no change) -> RCR11105 (no change) -> RCR10471 (no change) -> RCR11165 (some change)
# glutamate is unique to the drug treatment conditions - RCR20594
plot_data = data %>% 
  filter(reaction_id %in% test_rxns) %>% 
  filter(timepoint == "24h") %>% 
  left_join(reaction.names %>% select(RXNID, EQUATION), by = c('reaction_id'='RXNID'))
#mutate(Drug = ifelse(Drug %in% c('5FU','Ace','Dox'), 'treatment', 'control'))
plot_data$Drug = factor(plot_data$Drug, levels = c('5FU','DMSO1','Ace','Dox','DMSO2'))
ggplot(plot_data, aes(x = Drug, y = flux, fill = Drug)) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_manual(values = v_colors) + 
  facet_wrap(~reaction_id, scales = "free", ncol = 8) + 
  # stat_compare_means(comparisons = list(c('5FU','DMSO1'),
  #                                       c('Ace','DMSO2'),
  #                                       c('Dox','DMSO2')), 
  #                    # method.args = list(alternative = "greater"), 
  #                    method = "wilcox.test", 
  #                    paired = FALSE) + 
  xlab("") + 
  ylab("Predicted flux") + 
  theme(strip.text.x = element_text(size = 7), 
        legend.position = "none")

# figure out the next reaction
test_rxns = reaction.names %>% filter(str_detect(EQUATION, "xanthosine")) %>% 
  inner_join(data %>% select(reaction_id, condition) %>% unique() %>% filter(condition == "ace_24h"), by = c('RXNID'='reaction_id')) %>% 
  pull(RXNID)


# glutamate is unique to the drug treatment conditions - RCR20594
plot_data = data %>% 
  filter(reaction_id %in% c('RCR10294','RCR11159')) %>% 
  filter(timepoint == "24h") %>% 
  left_join(reaction.names %>% select(RXNID, EQUATION), by = c('reaction_id'='RXNID'))
#mutate(Drug = ifelse(Drug %in% c('5FU','Ace','Dox'), 'treatment', 'control'))
plot_data$Drug = factor(plot_data$Drug, levels = c('5FU','DMSO1','Ace','Dox','DMSO2'))
Figure5C = ggplot(plot_data, aes(x = Drug, y = flux, fill = Drug)) + 
  geom_boxplot() + 
  theme_bw() +
  scale_fill_manual(values = v_colors) + 
  facet_wrap(~EQUATION, scales = "free", ncol = 8) + 
  # stat_compare_means(comparisons = list(c('5FU','DMSO1'),
  #                                       c('Ace','DMSO2'),
  #                                       c('Dox','DMSO2')), 
  #                    # method.args = list(alternative = "greater"), 
  #                    method = "wilcox.test", 
  #                    paired = FALSE) + 
  xlab("") + 
  ylab("Predicted flux") + 
  theme(strip.text.x = element_text(size = 7), 
        legend.position = "none")
Figure5C

# RCR10294: 25630;24651;100364062
# RCR11159: 25630;24651;100364062	

# Plot the gene expression values
plot_data = transcript.counts %>% filter(NAME %in% c('100364062')) %>% filter(Timepoint == "24h")
plot_data$Drug = factor(plot_data$Drug, levels = c('5FU','DMSO1','Ace','Dox','DMSO2'))
plot_data$Timepoint = factor(plot_data$Timepoint, levels= c('6h','24h'))
Figure5D = ggplot(plot_data, aes(x = Drug, y = counts)) + 
  geom_jitter(aes(fill = Drug), width = 0.25, color = "black", pch = 21) + 
  theme_bw() +
  scale_fill_manual(values = v_colors) + 
  # facet_wrap(~NAME + Timepoint, scales = "free") + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5)) + 
  xlab("") + ylab(expression(paste("Scaled TPMs for ", italic("PKM2")))) 
Figure5D

# Combine figures into one figure
gridlayout <- data.matrix(c(1,NA,2)) %>% t()
Figure5 <- arrangeGrob(Figure5C, Figure5D,
                                   ncol = 3, 
                                   nrow = 1, 
                                   #heights = relative_heights,
                                   widths = c(0.60+0.025,0.025,0.35),
                                   layout_matrix = gridlayout) 

ggsave('results/figures/Figure5.png',Figure5, 
       dpi = 600, width = 8, height = 3, units = 'in')
