# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: Early SARS-CoV-2 Europe Phylodynamics
#        V\ Y /V    ggplot theme and color palettes
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)
library(ggpattern)
library(scales)
library(cowplot)


theme_dino2 <- function() {
  theme_minimal() +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          legend.box="vertical",
          legend.text=element_text(size = 7),
          legend.title=element_text(size = 7),
          strip.background = element_rect(fill = "white", colour = "grey80"),
          strip.text = element_text(size = 8, color = "grey30"),
          axis.title.y = element_text(size = 8, #face = "bold", 
                                      margin = margin(r = 10)),
          axis.title.x = element_text(size = 8, #face = "bold", 
                                      margin = margin(t = 10)),
          panel.grid.major.x =  element_line(colour = "grey80", size = 0.2, linetype = 2),
          panel.grid.major.y =  element_line(colour = "grey80", size = 0.2, linetype = 2),
          #panel.grid.minor.y =  element_line(colour = "grey80", size = 0.2,  linetype = 3),
          #panel.background =  element_blank(),
          axis.text.x = element_text(size = 7, angle = 0, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size=8, face = "bold"),
          axis.ticks.length.x=unit(0.1, "cm"),
          axis.ticks.length.y=unit(0.1, "cm"),
          axis.ticks = element_line(colour = "grey50", size = 0.2))
  #plot.margin = margin(3, 7, 3, 1.5))
}

scales::show_col(ggthemes::colorblind_pal()(8))
pal_dino <- c("#000000", "#E69F00", "#56B4E9")
pal_dino2 <- c("non-GLM analysis" = "#000000", 
               "sampling-GLM" = "#467057", 
               "sampling-GLM 12p" = "#276316", 
               "sampling-GLM + error" = "#A0BC91")
ecolors <- c("#1e3d59", "#f5f0e1", "#ff6e40", "#ffc13b")

pal = c("#00504D", "#1F8E8A", "#2FB9B3")
grey_pal = c("#6C6C6C", "#383838", "#CFCFCF")


pal_performance <- c("#1F8E8A", "#000000", "#A0BC91", "#276316", "#467057", "#000000")
