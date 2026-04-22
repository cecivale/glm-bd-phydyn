# Bethany Allen  May 13th 2025
# modified by Ceci April 2026
# Code to plot and summarise proxy data

library(tidyverse)
library(deeptime)
library(gridExtra)
library(ggpubr)
library(ggthemes)
library(corrplot)
source("R/plot_opts.R")

# Specify ages and midpoints
# ages <- c(66.0, 89.8, 121.4, 145.0, 161.5, 174.7, 201.4, 237.0, 251.9)
ages <- c(66.0, 89.8, 121.4, 145.0, 161.5, 174.7, 201.4, 237.0)

# Read proxies table
proxies <- read_csv("data/predictors.csv")

proxies_long <- proxies %>%
  pivot_longer(-var, names_to = "period", values_to = "value") %>%
  mutate(age = rep(ages, 3 * 4),
         period = str_replace(period, "-", " -")) %>%
  mutate( predictor_pretty = str_replace(var, "collectionCounts", "Collections"),
          predictor_pretty = str_replace(predictor_pretty, "formationCounts", "Formations"),
          predictor_pretty = str_replace(predictor_pretty, "gridCellCounts", "Spatial Grid Cells"),
          predictor_pretty = str_replace(predictor_pretty, "_r1", " (perm1)"),
          predictor_pretty = str_replace(predictor_pretty, "_r2", " (perm2)"),
          predictor_pretty = str_replace(predictor_pretty, "_rev", " (rev)"),
          predictor_class = str_replace(predictor_pretty, " \\(perm1\\)| \\(perm2\\)| \\(rev\\)", ""))
  

cor <- cor(t(proxies[-1]))
corrplot(cor)

# Test for correlations
# Collections versus formations
cor.test(as.numeric(proxies[proxies$var == "formationCounts", -1]), 
         as.numeric(proxies[proxies$var == "collectionCounts", -1]))
# Formations versus area
cor.test(as.numeric(proxies[proxies$var == "formationCounts", -1]), 
                    as.numeric(proxies[proxies$var == "gridCellCounts", -1]))
# Collections versus area
cor.test(as.numeric(proxies[proxies$var == "collectionCounts", -1]), 
         as.numeric(proxies[proxies$var == "gridCellCounts", -1]))

proxy_plot <- ggplot(proxies_long) +
  geom_step(aes(age, value), direction = "vh") +
  facet_wrap(~predictor_pretty, scales = "free_y", ncol = 3, dir = "v") +
  theme_dino2() +
  xlab("Time before present (Ma)") +
  scale_x_reverse()

ggsave(file = "figures/proxies_step.pdf", plot = proxy_plot, width = 14.7,
       height = 11, units = "cm")
