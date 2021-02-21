# Written by Bonnie Dougherty, 2020-06-01
# Analyzing metabolomics data
# Adapted from: https://www.intechopen.com/books/metabolomics-fundamentals-and-applications/processing-and-visualization-of-metabolomics-data-using-r

library(ggplot2)
library(xlsx)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(viridis)
library(ggpubr)


######## Impute, log transform, and mean center data #################3
# Read in the original raw intensity values
raw.data <- read.xlsx(file = "results/supplement/Supplemental File 3.xlsx", 
                      sheetName = "raw", startRow = 6, 
                      header = TRUE)
colnames(raw.data)[79:81] <- c('BLANK.1','BLANK.2','BLANK.3')
# Keep PUBCHEM and KEGG for mapping back to the iRno model
raw.data <- raw.data %>% dplyr::select(-PATHWAY_SORTORDER, -SUPER_PATHWAY, -SUB_PATHWAY, 
                                       -COMP_ID, -PLATFORM, -CHEMICAL_ID, -RI, -MASS, -CAS, 
                                       -X......................................Group.HMDB, 
                                       -PUBCHEM, - KEGG)

# Determine the percent of samples that are NA for each metabolite
# Remove metabolites that are > 60% NAs (only 40% data)
metabolites <- raw.data %>% 
  separate(BIOCHEMICAL, c('X','number'), sep = "X - ", remove = FALSE) %>% 
  filter(is.na(number)) %>% 
  select(-X, -number) %>% 
  pivot_longer(-BIOCHEMICAL,
               names_to = "condition", 
               values_to = "value",
               values_drop_na = FALSE) %>% 
  group_by(BIOCHEMICAL) %>% 
  summarize(n = sum(is.na(value))) %>% 
  filter(n < 0.6*68)

# Impute values as the minimum of the row/2
# Should the BLANK values be included for imputing? or just comparisons?
imputed.data <- raw.data %>% 
  inner_join(metabolites %>% select(BIOCHEMICAL)) %>% 
  pivot_longer(-BIOCHEMICAL, 
               names_to = "condition",
               values_to = "value", 
               values_drop_na = FALSE)
minimum.values <- imputed.data %>% 
  filter(!is.na(value)) %>% 
  group_by(BIOCHEMICAL) %>% 
  summarize(min = min(value)) %>% 
  ungroup()
imputed.data <- imputed.data %>% left_join(minimum.values) %>% 
  mutate(value = ifelse(is.na(value), min/2, value)) %>% 
  select(-min) %>% 
  pivot_wider(names_from = "condition", 
              values_from = "value")

# Examine the heteroscedasticity of the data
# plot the mean of each metabolite vs variance
summary.data <- imputed.data %>% 
  pivot_longer(-BIOCHEMICAL, 
               names_to = "condition",
               values_to = "value") %>% 
  group_by(BIOCHEMICAL) %>% 
  summarize(met.mean = mean(value), 
            met.sd = sd(value))
ggplot(summary.data, aes(x = met.mean, y = met.sd)) + 
  geom_point() + 
  theme_bw()

# log2 scale intensities to remove heteroscedasticity
log.data <- imputed.data %>% 
  pivot_longer(-BIOCHEMICAL, 
               names_to = "condition",
               values_to = "value") %>% 
  mutate(value = log2(value)) %>% 
  pivot_wider(names_from = "condition", 
              values_from = "value")

# # Pareto scale values within each metabolite
# paretoscale <- function(z) {
#   rowmean <- apply(z, 1, mean) # row means
#   rowsd <- apply(z, 1, sd) # row standard deviation
#   rowsqrtsd <- sqrt(rowsd) #sqrt of sd
#   rv <- sweep(z, 1, rowmean, "-") # mean center
#   rv <- sweep(rv, 1, rowsqrtsd, "/") # divide by sqrtsd
#   return(rv)
# }

# Mean scale
meanscale <- function(z){
  rowmean <- apply(z, 1, mean)
  rowmin <- apply(z, 1, min)
  rowmax <- apply(z, 1, max)
  rv <- sweep(z, 1, rowmean, "-")
  rv <- sweep(rv, 1, (rowmax - rowmin), "/")
}

scaled.data <- log.data
scaled.data[,-1] <- meanscale(scaled.data[,-1])

write.xlsx(scaled.data %>% as.data.frame(),
           file = "results/supplement/Supplemental File 3.xlsx", 
           sheetName = "scaled", 
           row.names = FALSE, append = TRUE)

######################## Perform Mann-Whitney for time and condition ##################################

metadata.master = read.xlsx(file = "data/metabolomics/metabolomics_metadata.xlsx", 
                                sheetIndex = 1)

# Start from the scaled data
data <- scaled.data %>% 
  pivot_longer(-BIOCHEMICAL, names_to = "dataset", values_to = "value", values_drop_na = FALSE) %>% 
  left_join(metadata.master, by = c("dataset" = "sample.ID")) %>% 
  select(BIOCHEMICAL, group, Drug, Timepoint, value) %>% 
  filter(!(group == "Blank" & Timepoint == "6h")) %>% 
  filter(!(group == "Blank" & Timepoint == "24h"))

# Assign factor levels
data$Drug <- factor(data$Drug, levels = c("5FU","DMSO1","Ace","Dox","DMSO2", "Blank"))
data$Timepoint <- factor(data$Timepoint, levels = c("0h","6h","24h"))
data$BIOCHEMICAL <- factor(data$BIOCHEMICAL)

# Mann-Whitney test for differences
all.results <- data.frame()
comparisons <- data.frame(Comparison = c('5FU_24h - DMSO1_24h', '5FU_6h - DMSO1_6h',
                                         'Ace_24h - DMSO2_24h','Ace_6h - DMSO2_6h',
                                         'Dox_24h - DMSO2_24h','Dox_6h - DMSO2-6h', 
                                         '5FU_24h - Blank', '5FU_6h - Blank',
                                         'Ace_24h - Blank', 'Ace_6h - Blank',
                                         'DMSO1_24h - Blank','DMSO1_6h - Blank', 
                                         'DMSO2_24h - Blank','DMSO2_6h - Blank',
                                         'Dox_24h - Blank','Dox_6h - Blank'),
                          x = c('5FU_24h','5FU_6h','Ace_24h','Ace_6h','Dox_24h','Dox_6h',
                                '5FU_24h','5FU_6h','Ace_24h','Ace_6h','DMSO1_24h','DMSO1_6h','DMSO2_24h','DMSO2_6h','Dox_24h','Dox_6h'),
                          y = c('DMSO1_24h','DMSO1_6h','DMSO2_24h','DMSO2_6h','DMSO2_24h','DMSO2_6h',
                                'Blank','Blank','Blank','Blank','Blank','Blank','Blank','Blank','Blank','Blank'),
                          group = c('5FU_24h','5FU_6h','Ace_24h','Ace_6h','Dox_24h','Dox_6h', 
                                'Blank','Blank','Blank','Blank','Blank','Blank','Blank','Blank','Blank','Blank'))

for (metabolite in levels(data$BIOCHEMICAL)){
  for (comparison in levels(factor(comparisons$Comparison))){
      x <- comparisons %>%
        filter(Comparison == comparison) %>%
        inner_join(data %>% filter(BIOCHEMICAL == metabolite), by = c('x'='group'))
      y <- comparisons %>%
        filter(Comparison == comparison) %>%
        inner_join(data %>% filter(BIOCHEMICAL == metabolite), by = c('y'='group'))
      temp <- wilcox.test(x = x$value,
                  y = y$value)

      all.results <- all.results %>%
        rbind(data.frame(group = (comparisons %>% filter(Comparison == comparison))$group,
                         comparison = (comparisons %>% filter(Comparison == comparison))$Comparison,
                         metabolite = metabolite,
                         diff = mean(x$value) - mean(y$value),
                         p.unadj = temp$p.value))
    
  }
}

# Correct for multiple hypotheses
all.results$p.adj <- p.adjust(all.results$p.unadj, method = "BH")
all.results <- all.results %>% 
  mutate(p.adj = ifelse(is.na(p.adj), 1, p.adj))

# Plot glucose results
data$group = factor(data$group, 
                    level = c('Blank','5FU_6h','DMSO1_6h','5FU_24h','DMSO1_24h', 
                              'Ace_6h','Dox_6h','DMSO2_6h','Ace_24h','Dox_24h','DMSO2_24h'))
ggplot(data %>% filter(BIOCHEMICAL == "glucose"), 
       aes(x = group, y = value)) +
  geom_jitter(width = 0.15) + 
  stat_compare_means(comparisons = list(c('Dox_24h','DMSO2_24h'), 
                                        c('5FU_24h','DMSO1_24h'), 
                                        c('5FU_6h','DMSO1_6h')))

# Confirm direction of change in production/consumption
ggplot(data %>% filter(BIOCHEMICAL == "1-methyladenosine"), 
       aes(x = group, y = value)) +
  geom_jitter(width = 0.15) + 
  stat_compare_means(comparisons = list(c('Dox_24h','Blank'), 
                                        c('Ace_24h','Blank'), 
                                        c('DMSO1_24h','Blank'), 
                                        c('DMSO2_24h','Blank')))


blank.comparisons <- all.results %>% 
  filter(group == "Blank") %>% 
  inner_join(metadata.master %>% 
               select(group, Drug, Timepoint) %>% 
               unique() %>% 
               filter(!(group == "Blank" & Timepoint == "6h")) %>% 
               filter(!(group == "Blank" & Timepoint == "24h"))) %>% 
  separate(comparison, sep = " - Blank", into = c("comparison","junk"), remove = TRUE) %>% 
  select(-junk)

save(blank.comparisons, file = "data/blank_metabolomics.Rdata")

# Inter join with meta data
all.results <- all.results %>% 
  select(-comparison) %>% 
  inner_join(metadata.master %>% 
               select(group, Drug, Timepoint) %>% 
               unique() %>% 
               filter(!(group == "Blank" & Timepoint == "6h")) %>% 
               filter(!(group == "Blank" & Timepoint == "24h")))

# Figure out number of metabolites changing per condition
all.results %>% filter(p.adj < 0.1) %>% group_by(group) %>% dplyr::count()

# Assign produced/consumed/no change for metabolites
direction <- all.results %>% 
  filter(group == "Blank") %>% 
  filter(p.adj< 0.1) %>% 
  group_by(metabolite) %>% 
  summarize(consumed = length(which(diff < 0)), 
            produced = length(which(diff > 0))) %>% 
  mutate(direction = ifelse(produced == 0 & consumed != 0, "Consumed", 
                            ifelse(consumed == 0 & produced != 0, "Produced", "Ambiguous"))) %>% 
  right_join(data.frame(metabolite = levels(data$BIOCHEMICAL))) %>% 
  mutate(direction = ifelse(is.na(direction), "No change", direction)) %>% 
  select(metabolite, direction)


final.data <- all.results %>% 
  filter(group != "Blank") %>% 
  right_join(direction)

final.data %>% filter(p.adj < 0.1) %>% group_by(group) %>% dplyr::count()


write.xlsx(final.data,
           file = "results/supplement/Supplemental File 3.xlsx", 
           sheetName = "differential", 
           row.names = FALSE, append = TRUE)


