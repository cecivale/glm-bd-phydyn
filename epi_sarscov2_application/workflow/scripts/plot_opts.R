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

theme_eubdmm <- function() {
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
               axis.text.x = element_text(size = 7, angle = 0, hjust = 1, vjust = 1),
               axis.text.y = element_text(size = 7),
               plot.title = element_text(size=8, face = "bold"),
               axis.ticks.length.x=unit(0.1, "cm"),
               axis.ticks.length.y=unit(0.1, "cm"),
               axis.ticks = element_line(colour = "grey50", size = 0.2))
}


pal_eubdmm <- c("#ed2f15", "#1fb9de", "#00a669", "#2e4a87", "#ffaa4fFF")
pal_eubdmm2 <- c("#000000", "#1a5f75", "#809e00", "#9b0821", "#eacab5")
demes <- c(CN = "China", FR = "France", DE = "Germany", IT = "Italy", OE = "OtherEU")
names(pal_eubdmm2) <- sort(demes)
ecolors <- c(B = "#1e3d59", OM = "#f5f0e1", IM = "#ff6e40", M = "#ffc13b")
dataset_pal = c("#00504D", "#1F8E8A", "#2FB9B3")
grey_pal = c("grey30", "grey50", "grey70")

ansys_pal <- c("#1F8E8A", "black")
