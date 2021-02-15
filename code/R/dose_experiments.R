# Import all data and process to plot only 24h data

require(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
select <- dplyr::select
library(multcomp)

library(grid)
library(gridExtra)
library(ggpubr)
library(cowplot)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

data.dead.files <- c("data/experimental/20190103 -- 20181005 DOX correction.csv",
                "data/experimental/20190103 -- 20181102 DOX correction.csv",
                "data/experimental/20190103 -- 20181116 DOX correction.csv")

data.viability.files <- c("data/experimental/20190108 -- 20181005 Cell Viability all data.csv",
                          "data/experimental/20190108 -- 20181102 Cell Viability all data.csv",
                          "data/experimental/20190108 -- 20181116 Cell Viability all data.csv")

data.dead <- data.frame()
for(k in 1:3){
  data <- read.csv(file = data.dead.files[k]) %>% 
    filter(drug != "Empty")
  
  # Filter out data based on the number of cells identified in each FOV
  ranges <- quantile(data$total_cells)
  lower <- ranges[2] - 1.5*IQR(data$total_cells)
  upper <- ranges[4] + 1.5*IQR(data$total_cells)
  
  data <- data %>% 
    filter(total_cells > lower[[1]]) %>% 
    filter(total_cells < upper[[1]])
  
  # Group data from each FOV for a cell death measure per well ID
  data <- data %>%
    group_by(well_ID, timepoint, drug, conc) %>% 
    summarize(dead_cells = sum(dead_cells), 
              dead_cells_DOX = sum(dead_cells_DOX), 
              total_cells = sum(total_cells)) %>% 
    ungroup() %>% 
    mutate(percent_dead = dead_cells_DOX/total_cells*100) %>% 
    filter(drug != "blank" & drug != "ctl" & drug != "Blank" & drug != "Empty" & drug != "stauro") %>% 
    select(timepoint, drug, conc, percent_dead) %>% 
    mutate(replicate = k)
  
  
  # Code adapted from: https://stackoverflow.com/questions/23760401/sequential-use-of-the-williams-test-to-determine-the-minimum-effective-dose-afte
  # Incorporate statistical tests for all concentration and time points
  all_results = data.frame()
  for(drug_ID in c('Ace','Dox','5FU')) 
  {
    for(time_step in c("6h","24h","48h"))
    {
      stat_test <- data %>% 
        filter(drug == drug_ID & timepoint == time_step) %>% 
        dplyr::select(conc, percent_dead)
      
      stat_test <- data %>%
        filter(drug == "DMSO" & timepoint == time_step) %>%
        select(conc, percent_dead) %>% 
        rbind(stat_test)
      
      stat_test$conc <- factor(stat_test$conc)
      
      M <- lm(percent_dead ~ conc, data = stat_test)
      trend <- glht(M, linfct = mcp(conc = "Dunnett"), alternative = "greater")
      results <- summary(trend)
      concs <- unlist(results$model$xlevels)[2:length(unlist(results$model$xlevels))]
      temp <- data.frame(drug = rep(c(drug_ID), each = length(concs)), timepoint = rep(c(time_step), each = length(concs)),conc = concs, pvalue = results$test$pvalues)
      
      all_results <- rbind(all_results, temp)
    }
  }
  all_results <- all_results %>% 
    mutate(significance = pvalue < 0.05)
  all_results$conc <- as.numeric.factor(all_results$conc)
  
  data <- left_join(data, all_results, by = c("drug","timepoint","conc")) %>% 
    select(drug, conc, timepoint, percent_dead, significance, replicate) 
  
  data.dead <- rbind(data.dead, data)
}

data.dead <- data.dead %>% 
  filter(drug != "DMSO") %>% 
  mutate(measure = "percent_dead") %>% 
  mutate(value = percent_dead) %>% 
  select(drug, conc, timepoint, value, significance, measure, replicate)


data.viability <- data.frame()
for(k in 1:3){
  data_viability <- read.csv(file = data.viability.files[k]) %>% 
    mutate(value = raw.value - background) %>% 
    select(drug, conc, timepoint, value)
  
  DMSO_control <- data_viability %>% 
    filter(drug == "DMSO") %>% 
    group_by(timepoint) %>% 
    summarize(DMSO_value = mean(value)) %>% 
    select(timepoint, DMSO_value)
  
  data_viability <- left_join(data_viability, DMSO_control) %>% 
    mutate(value = value/DMSO_value) %>% 
    select(drug, conc, timepoint, value) %>% 
    mutate(replicate = k)
  
  # Running statistical test on all data
  all_results = data.frame()
  data <- data_viability
  for(i in c('Ace','Dox','5FU')) 
  {
    for(time_step in c("6h","24h","48h"))
    {
      stat_test <- data %>% 
        filter(drug == i & timepoint == time_step) %>% 
        dplyr::select(conc, value)
      
      stat_test <- data %>%
        filter(drug == "DMSO" & timepoint == time_step) %>% 
        select(conc, value) %>% 
        rbind(stat_test)
      
      stat_test$conc <- factor(stat_test$conc)
      
      M <- lm(value ~ conc, data = stat_test)
      trend <- glht(M, linfct = mcp(conc = "Dunnett"), alternative = "less")
      results <- summary(trend)
      concs <- unlist(results$model$xlevels)[2:length(unlist(results$model$xlevels))]
      temp <- data.frame(drug = rep(c(i), each = length(concs)), timepoint = rep(c(time_step), each = length(concs)),conc = concs, pvalue = results$test$pvalues)
      
      all_results <- rbind(all_results, temp)
    }
  }
  all_results <- all_results %>% 
    mutate(significance = (pvalue < 0.05))
  
  all_results$conc <- as.numeric.factor(all_results$conc)
  
  data_viability <- left_join(data_viability, all_results, by = c("drug","conc","timepoint")) %>% 
    filter(drug != "DMSO") %>% 
    mutate(measure = "viability") %>% 
    select(drug, conc, timepoint, value, significance, measure, replicate)
  
  data.viability = rbind(data.viability, data_viability)
}

data.viability <- data.viability %>% 
  mutate(value = value*100)

# Plot all data together
# percent dead: data.dead
# viability: data.viability

all.data <- rbind(data.dead, data.viability) %>% 
  filter(timepoint == "24h")
all.data$drug <- factor(all.data$drug, levels = c("Ace","Dox","5FU"))
all.data$replicate <- factor(all.data$replicate)

# Changing margins for plotting data
# Plot the 24 hour data (Figure 1)
Figure1A1 = ggplot(all.data %>% 
                     filter(measure == "viability") %>% 
                     filter(drug == "5FU"), 
                   aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate),  #, shape = replicate
              width = 0.05) +  
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  ylim(0, 150) +
  xlab("Concentration (mM)") + 
  ggtitle("5FU") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("Percent cell reducing\npotential relative to control") + 
  theme_bw() + 
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")
Figure1A1

Figure1A2 = ggplot(all.data %>% 
                     filter(measure == "viability") %>% 
                     filter(drug == "Ace"), 
                   aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate),  #, shape = replicate
              width = 0.05) +  
  ylim(0,150) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (mM)") + 
  ggtitle("Ace") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")
Figure1A2

Figure1A3 = ggplot(all.data %>% 
                     filter(measure == "viability") %>% 
                     filter(drug == "Dox"), 
                   aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate),  #, shape = replicate
              width = 0.05) + 
  ylim(0,150) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (uM)") + 
  ggtitle("Dox") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")
Figure1A3

Figure1B1 <- ggplot(all.data %>% 
                      filter(measure == "percent_dead") %>% 
                      filter(drug == "5FU"), 
                    aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate), #, shape = replicate
              width = 0.05) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (mM)") + 
  ggtitle("5FU") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("Percent cell death \n") + 
  theme_bw() + 
  ylim(0,100) +
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")
Figure1B1

Figure1B2 <- ggplot(all.data %>% 
                      filter(measure == "percent_dead") %>% 
                      filter(drug == "Ace"), 
                    aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate), #, shape = replicate
              width = 0.05) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (mM)") + 
  ggtitle("Ace") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  ylim(0,100) +
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")
Figure1B2

Figure1B3 <- ggplot(all.data %>% 
                      filter(measure == "percent_dead") %>% 
                      filter(drug == "Dox"), 
                    aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate), #, shape = replicate
              width = 0.05) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (uM)") + 
  ggtitle("Dox") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  ylim(0,100) +
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")
Figure1B3

# Combine figures into one figure
gridlayout <- matrix(c(1,2, 3, 
                       NA, NA, NA, 
                       4, 5, 6), nrow = 3) %>% t()
Figure1 <- arrangeGrob(Figure1A1, Figure1A2, Figure1A3, 
                       Figure1B1, Figure1B2, Figure1B3,
                       ncol = 3, 
                       nrow = 3,
                       widths = c(1/3, 1/3, 1/3), 
                       heights = c(0.475, 0.05, 0.475),
                       layout_matrix = gridlayout) %>% 
  as_ggplot() + 
  draw_text(text = c("A","B"), 
            size = 24, 
            hjust = 0, 
            x = c(0.01, 0.01),
            y = c(0.98, 0.48))
# Looks good with an 8:6 width:height ratio
ggsave('results/figures/Figure1.png',Figure1, 
       dpi = 600, width = 7*1.25, height = 5*1.1, units = 'in')

### Generate Supplemental Figure 2
all.data <- rbind(data.dead, data.viability) %>% 
  filter(timepoint == "6h")
all.data$drug <- factor(all.data$drug, levels = c("Ace","Dox","5FU"))
all.data$replicate <- factor(all.data$replicate)

# Changing margins for plotting data
# Plot the 24 hour data (Figure 1)
SuppFigure2A1 = ggplot(all.data %>% 
                     filter(measure == "viability") %>% 
                     filter(drug == "5FU"), 
                   aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate),  #, shape = replicate
              width = 0.05) +  
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  ylim(0, 150) +
  xlab("Concentration (mM)") + 
  ggtitle("5FU") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("Percent cell reducing\npotential relative to control") + 
  theme_bw() + 
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")

SuppFigure2A2 = ggplot(all.data %>% 
                     filter(measure == "viability") %>% 
                     filter(drug == "Ace"), 
                   aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate),  #, shape = replicate
              width = 0.05) +  
  ylim(0,150) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (mM)") + 
  ggtitle("Ace") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")

SuppFigure2A3 = ggplot(all.data %>% 
                     filter(measure == "viability") %>% 
                     filter(drug == "Dox"), 
                   aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate),  #, shape = replicate
              width = 0.05) + 
  ylim(0,150) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (uM)") + 
  ggtitle("Dox") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")

SuppFigure2B1 <- ggplot(all.data %>% 
                      filter(measure == "percent_dead") %>% 
                      filter(drug == "5FU"), 
                    aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate), #, shape = replicate
              width = 0.05) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (mM)") + 
  ggtitle("5FU") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("Percent cell death \n") + 
  theme_bw() + 
  ylim(0,100) +
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")

SuppFigure2B2 <- ggplot(all.data %>% 
                      filter(measure == "percent_dead") %>% 
                      filter(drug == "Ace"), 
                    aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate), #, shape = replicate
              width = 0.05) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (mM)") + 
  ggtitle("Ace") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  ylim(0,100) +
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")

SuppFigure2B3 <- ggplot(all.data %>% 
                      filter(measure == "percent_dead") %>% 
                      filter(drug == "Dox"), 
                    aes(x = conc, y = value, shape = timepoint)) + 
  geom_jitter(aes(color = significance, shape = replicate), #, shape = replicate
              width = 0.05) + 
  scale_x_log10() + 
  scale_color_manual(values=c("grey75", "grey0")) + 
  xlab("Concentration (uM)") + 
  ggtitle("Dox") + 
  #ylab(expression(atop("Reducing potential", paste("as percent of control")))) +
  ylab("") + 
  theme_bw() + 
  ylim(0,100) +
  theme(
    #strip.background = element_blank(),
    strip.placement = "outside", 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "none")

# Combine figures into one figure
gridlayout <- matrix(c(1,2, 3, 
                       NA, NA, NA, 
                       4, 5, 6), nrow = 3) %>% t()
SuppFigure2 <- arrangeGrob(SuppFigure2A1, SuppFigure2A2, SuppFigure2A3, 
                       SuppFigure2B1, SuppFigure2B2, SuppFigure2B3,
                       ncol = 3, 
                       nrow = 3,
                       widths = c(1/3, 1/3, 1/3), 
                       heights = c(0.475, 0.05, 0.475),
                       layout_matrix = gridlayout) %>% 
  as_ggplot() + 
  draw_text(text = c("A","B"), 
            size = 24, 
            hjust = 0, 
            x = c(0.01, 0.01),
            y = c(0.98, 0.48))
# Looks good with an 8:6 width:height ratio
ggsave('results/figures/SupplementalFigure2.png',SuppFigure2, 
       dpi = 600, width = 7*1.25, height = 5*1.1, units = 'in')

