
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

## Analyzing Seahorse data
data <- read.csv(file = "data/experimental/20181206 -- Seahorse 24 hour drug treatment.csv")

# Add in the drugs that correspond to each of the time measurements:
inhibitor <- c("Baseline","Baseline","Baseline","Baseline","Oligomycin","Oligomycin","Oligomycin","Oligomycin","FCCP","FCCP","FCCP","FCCP","Antimycin","Antimycin","Antimycin")

data <- cbind(data, inhibitor) %>% 
  gather(key = "Condition", value = "OCR", -Time, -inhibitor) 
data$Condition <- factor(data$Condition)
data$Replicate <- lapply(data$Condition, function(x) as.numeric(last(unlist(strsplit( as.character(x), split = '[.]')))))
data$Condition <- lapply(data$Condition, function(x) unlist(strsplit( as.character(x), split = '[.]'))[1])
data <- data %>% 
  mutate(Replicate = ifelse(Replicate == '1' | Replicate == '2' | Replicate == '3' | Replicate == '4' | Replicate == '5' | Replicate == '6', Replicate, '0')) %>% 
  mutate(Condition = ifelse(Condition == "X5FU", "5FU", Condition))
data$Condition <- as.character(data$Condition)
data$Replicate <- as.character(data$Replicate)

data$Replicate <- as.character(data$Replicate)

# Filter out bad data points based on non-responding cells
data <- data %>% 
  filter(!(Condition == "5FU" & Replicate == "0"))

# Plot the raw data
data.to.plot <- data %>% 
  group_by(Time, inhibitor, Condition) %>% 
  summarize(OCR_std = sd(OCR), 
            OCR = mean(OCR))
raw.plot <- ggplot(data, aes(x = Time, y = OCR)) + 
  geom_line(aes(color = Condition, shape = Replicate))
raw.plot
raw.plot <- ggplot(data.to.plot, aes(x = Time, y = OCR)) + 
  geom_point(aes(color = Condition)) + 
  geom_line(aes(color = Condition)) + 
  geom_errorbar(aes(ymin = OCR - OCR_std, ymax = OCR + OCR_std), 
                width = 2.5)
raw.plot

Baseline <- data %>% 
  filter(inhibitor == "Baseline") %>%  
  filter(Time == 0.00) %>% 
  dplyr::select(Condition, Replicate, OCR) %>% 
  group_by(Condition, Replicate) %>% 
  summarize(Baseline_var = var(OCR),
    Baseline = mean(OCR))

Oligomycin <- data %>% 
  filter(inhibitor == "Oligomycin") %>% 
  filter(Time == 31.00) %>% 
  dplyr::select(Condition, Replicate, OCR) %>% 
  group_by(Condition, Replicate) %>% 
  summarize(Oligomycin_var = 0,
            Oligomycin = OCR)

FCCP <- data %>% 
  filter(inhibitor == "FCCP") %>%
  filter(Time == 62.00) %>% 
  dplyr::select(Condition, Replicate, OCR) %>% 
  group_by(Condition, Replicate) %>% 
  summarize(FCCP_var = var(OCR), 
            FCCP = mean(OCR))

Antimycin <- data %>% 
  filter(inhibitor == "Antimycin") %>% 
  filter(Time == 93.00) %>% 
  dplyr::select(Condition, Replicate, OCR) %>% 
  group_by(Condition, Replicate) %>% 
  summarize(Antimycin_var = var(OCR), 
            Antimycin = mean(OCR))



### Calculate measures of respiration
# Definitions from Seahorse website: 
# non-mitochondrial respiration: minimum after Antimycin
# basal respiration: last rate before injection - non-mitochondrial
# max respiration: max after FCCP - non-mitochondrial
# proton leak: min after oligomycin - non-mitochondrial
# ATP production: last rate before oligo - min after oligomycin
# spare respiratory capacity: max - basal
# coupling efficient: ATP production/basal respiration * 100

all.data <- left_join(Baseline, Oligomycin, by = c("Condition","Replicate")) %>% 
  left_join(FCCP) %>% 
  left_join(Antimycin) %>% 
  mutate(non_mito = Antimycin,
         basal_respiration = Baseline - Antimycin,
         max_respiration = FCCP - Antimycin,
         proton_leak = Oligomycin - Antimycin,
         ATP_production = Baseline - Oligomycin,
         spare_respiratory = max_respiration - basal_respiration) %>% 
  dplyr::select(Condition, Replicate, non_mito:spare_respiratory) %>% 
  gather(key = "measure", value = "OCR", -Condition, -Replicate)

all.data$Condition <- factor(all.data$Condition, levels = c("DMSO","5FU","Ace","Dox"))

### INTEGRATING CELL COUNTS WITH RAW SEAHORSE DATA
# Load well_IDs
well_IDs <- read.csv(file = "data/experimental/20181206 -- Seahorse 24 hour drug treatment, cell counts.csv") %>% 
  dplyr::select(well_ID, drug) %>% 
  filter(drug != "blank" & drug != "Blank") %>%
  unique() %>% 
  cbind(data.frame(Replicate = rep(c("0","1","2","3","4"), times = 4))) 
colnames(well_IDs)[2] <- 'Condition'

# Load cell counts 
cell.data <- read.csv(file = "data/experimental/20181206 -- Seahorse 24 hour drug treatment, cell counts.csv") %>% 
  filter(timepoint == "24h") %>% 
  dplyr::select(well_ID, drug, dead_cells_DOX, total_cells) %>% 
  filter(drug != "Blank") %>% 
  group_by(well_ID, drug) %>% 
  summarize(dead_cells = sum(dead_cells_DOX), 
            total_cells = sum(total_cells), 
            percent_dead = dead_cells/total_cells, 
            live_cells = total_cells-dead_cells) %>% 
  ungroup()

all.data.normalized <- all.data %>% 
  left_join(well_IDs, by = c("Condition","Replicate")) %>% 
  left_join(cell.data, by = c("well_ID")) %>% 
  mutate(OCR_normalized = OCR/live_cells)

all.data.normalized$Condition <- factor(all.data.normalized$Condition, 
                                        levels = c("DMSO","Low","Med","High","Ace","Dox","5FU"))

# all.data.normalized.plot <- ggplot(all.data.normalized, aes(x = Condition, y = OCR_normalized)) + 
#   geom_point() + 
#   facet_wrap(~measure, scales = "free")
# all.data.normalized.plot

# Testing out statistical tests
# statistically different: 
# ATP_production (DMSO-Dox, 0.003; 5FU-Dox, 0.04)
# basal_respiration (Dox-DMSO, 0.02; 5FU-Dox, 0.04) 
# coupling_efficiency (Ace-DMSO, 0.03; Dox-DMSO, 0.0007; 5FU-DMSO, 0.07)
# spare_respiratory: (5FU-DMSO: 0.08)

stat.test.data = all.data.normalized
#   dplyr::select(-Replicate) %>%
#   pivot_longer(cols = ATP_production:spare_respiratory, values_to = "OCR", names_to = "measure")
stat.test.data$measure = factor(stat.test.data$measure)

results = data.frame()
for(current.measure in levels(stat.test.data$measure)){
  current.data = stat.test.data %>% filter(measure == current.measure)

  fit = aov(OCR ~ Condition, data = current.data)
  results = results %>% rbind(data.frame(TukeyHSD(fit)$Condition,
                                         measure = current.measure))
}

# units of OCR are pmol/min/cell
# Plot final data
library(ggpubr)

final.data = all.data.normalized %>% select(Condition, measure, OCR_normalized) %>% 
  rename(OCR = "OCR_normalized")
final.data$Condition = factor(final.data$Condition, levels = c("DMSO","5FU","Ace","Dox"))

measures.metadata = data.frame(measure = c("ATP_production","basal_respiration"), 
                               plot.name = c("ATP production","Basal respiration"))


# Visualize the data
FigureS1A = ggplot(data = final.data %>%  
                     filter(measure == "ATP_production" | measure == "basal_respiration") %>% 
                     left_join(measures.metadata), 
                aes(x = Condition, y = OCR)) + 
  geom_jitter(width = 0.1) +
  theme_bw() + 
  stat_compare_means(method = "anova", label.y = c(0.255), size = 5) + 
  stat_compare_means(method = "wilcox.test", ref.group = "DMSO", label = "p.signif",
                     label.y = c(0.225), size = 5) +
  facet_wrap(~plot.name, scales = "free") + 
  ylim(c(-0.01, 0.275)) + 
  xlab("") + ylab("OCR normalized to live cell counts\n(pmol/min/cell)") + 
  theme(plot.title = element_text(hjust = 0.5), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 12))
FigureS1A

ggsave('results/figures/SupplementalFigure1.png', FigureS1A, 
       dpi = 600, width = 6, height = 2.75, units = 'in')


