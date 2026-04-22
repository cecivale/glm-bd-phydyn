# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: Macro Dinosaurs GLM Phylodynamics
#        V\ Y /V    Figures
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(treedataverse)
library(ggridges)

source("R/utils.R")
source("R/plot_opts.R")
source("R/priors.R")


ansys_names <- c("no_GLM", "GLM_sumprior_3p", 
                 "GLM_sumprior_12p", "GLM_sumprior_errordistr_3p")


# 1. Analysis specs ------------------------------------------------------------
intervals <- tibble(epoch_age = 8:1, epoch_time = 1:8, 
                    interval = c("Early-Mid Triassic", "Late Triassic",
                              "Early Jurassic", "Middle Jurassic",
                              "Late Jurassic", "Berriasian- Barremian",
                              "Aptian- Turonian", "Coniacian- Maastrichtian"))

# 2. Trace files ---------------------------------------------------------------
traces_l <- lapply(ansys_names, function(ansys) {
  trace_file <-  paste0("results/Dinosaurs_", ansys, "_combined.log")
  
  trace_wide <- read_trace(trace_file, burninFrac = 0.1) %>%
    select(-starts_with("X"))
  
  var_pattern <- ifelse(ansys == "no_GLM", "\\.", "SVi[0-9]+$")
  var_sep     <- ifelse(ansys == "no_GLM", ".", "SVi")
  
  birth_cols <- grep(paste0("^birthRate", var_pattern), names(trace_wide))
  death_cols <- grep(paste0("^deathRate", var_pattern), names(trace_wide))
  
  idx <- str_extract(names(trace_wide)[birth_cols], "\\d+$")
  
  trace_wide[paste0("divRate", var_sep, idx)] <-
    trace_wide[birth_cols] - trace_wide[death_cols]
  
  trace_long <- trace_wide %>%
    pivot_longer(-1, names_to = "parameter", values_to = "value") %>%
    filter(!grepl("endtime", parameter)) %>%
    mutate(param = str_split(parameter, "SVi|\\.", simplify = TRUE)[, 1],
           param = ifelse(param == "samplingRateGLM", "samplingRateGLM", param),
           param = ifelse(grepl("samplingRateSV", param), "samplingRate", param),
           epoch_time = case_when(
             str_detect(parameter, "SVi") ~ 
               as.numeric(str_extract(parameter, "(?<=SVi)\\d+")) + 1,
             str_detect(parameter, "\\.") ~ 
               9 - as.numeric(str_extract(parameter, "(?<=\\.)\\d+"))),
           ansys = ansys) %>%
    left_join(intervals) %>%
    mutate(interval = factor(interval, levels = intervals$interval))
  
  return(trace_long)
})


traces <- bind_rows(traces_l) %>%
  mutate(ansys_pretty = factor(case_when(
    ansys == "no_GLM" ~ "non-GLM analysis",
    ansys == "GLM_sumprior_3p" ~ "sampling-GLM",
    ansys == "GLM_sumprior_12p" ~ "sampling-GLM 12p",
    ansys == "GLM_sumprior_errordistr_3p" ~ "sampling-GLM + error"),
    levels = c("non-GLM analysis", "sampling-GLM", "sampling-GLM 12p", "sampling-GLM + error")))

trace_summary <- traces %>%
  group_by(ansys, ansys_pretty, parameter, param, interval, epoch_age, epoch_time) %>%
  summarise(median = median(value),
            l95 = HDInterval::hdi(value)[1],
            h95 = HDInterval::hdi(value)[2], .groups = "drop") 

## 2.1 Birth - Violinplots -----------------------------------------------------

bp_birth <- ggplot() +
  geom_violin(
    data = traces %>% filter(param == "birthRate"),
    aes(x = epoch_time, y = value, fill = ansys_pretty, group = interaction(epoch_time, ansys_pretty)),
    alpha = 0.5, color = "transparent",
    position = position_dodge(0.9),
    scale = "width"
  ) +
  geom_errorbar(
    data = trace_summary %>% filter(param == "birthRate"),
    aes(x = epoch_time, ymin = l95, ymax = h95, color = ansys_pretty),
    position = position_dodge(0.9), linewidth = 0.5, width = 0.8
  ) +
  geom_point(
    data = trace_summary %>% filter(param == "birthRate"),
    aes(x = epoch_time, y = median, color = ansys_pretty),
    position = position_dodge(0.9), size = 1
  ) +
  annotate(
    "rect",
    xmin = c(1.5, 3.5, 5.5, 7.5),
    xmax = c(2.5, 4.5, 6.5, 8.5),
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.1, colour = NA
  ) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2, name = "Analysis:") +
  scale_fill_manual(values = pal_dino2, name = "Analysis:") +
  labs(x = "Interval", y = "Speciation rate") +
  coord_cartesian(ylim = c(0, 8)) +
  scale_x_continuous(
    breaks = seq(1.5, 8.5, 1),
    labels = \(x) stringr::str_wrap(intervals$interval[x], width = 10)
  )

bp_birth
ggsave(plot = bp_birth, filename = paste0("figures/", "bp_birth.pdf"), 
       width = 174, height = 80, units = "mm")


## 2.2 Death - Violinplots -----------------------------------------------------

bp_death <- ggplot() +
  geom_violin(
    data = traces %>% filter(param == "deathRate"),
    aes(x = epoch_time, y = value, fill = ansys_pretty, group = interaction(epoch_time, ansys_pretty)),
    alpha = 0.5, color = "transparent",
    position = position_dodge(0.9),
    scale = "width"
  ) +
  geom_errorbar(
    data = trace_summary %>% filter(param == "deathRate"),
    aes(x = epoch_time, ymin = l95, ymax = h95, color = ansys_pretty),
    position = position_dodge(0.9), linewidth = 0.5, width = 0.8
  ) +
  geom_point(
    data = trace_summary %>% filter(param == "deathRate"),
    aes(x = epoch_time, y = median, color = ansys_pretty),
    position = position_dodge(0.9), size = 1
  ) +
  annotate(
    "rect",
    xmin = c(1.5, 3.5, 5.5, 7.5),
    xmax = c(2.5, 4.5, 6.5, 8.5),
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.1, colour = NA
  ) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2, name = "Analysis:") +
  scale_fill_manual(values = pal_dino2, name = "Analysis:") +
  labs(x = "Interval", y = "Extinction rate") +
  coord_cartesian(ylim = c(0, 8)) +
  scale_x_continuous(
    breaks = seq(1.5, 8.5, 1),
    labels = \(x) stringr::str_wrap(intervals$interval[x], width = 10)
  )

bp_death

ggsave(plot = bp_death, filename = paste0("figures/", "bp_death.pdf"), 
       width = 174, height = 80, units = "mm")


## 2.3 Diversification - Violinplots -------------------------------------------

bp_div <- ggplot() +
  geom_violin(
    data = traces %>% filter(param == "divRate"),
    aes(x = epoch_time, y = value, fill = ansys_pretty, group = interaction(epoch_time, ansys_pretty)),
    alpha = 0.5, color = "transparent",
    position = position_dodge(0.9),
    scale = "width"
  ) +
  geom_errorbar(
    data = trace_summary %>% filter(param == "divRate"),
    aes(x = epoch_time, ymin = l95, ymax = h95, color = ansys_pretty),
    position = position_dodge(0.9), linewidth = 0.5, width = 0.8
  ) +
  geom_point(
    data = trace_summary %>% filter(param == "divRate"),
    aes(x = epoch_time, y = median, color = ansys_pretty),
    position = position_dodge(0.9), size = 1
  ) +
  annotate(
    "rect",
    xmin = c(1.5, 3.5, 5.5, 7.5),
    xmax = c(2.5, 4.5, 6.5, 8.5),
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.1, colour = NA
  ) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2, name = "Analysis:") +
  scale_fill_manual(values = pal_dino2, name = "Analysis:") +
  labs(x = "Interval", y = "Diversification rate") +
  coord_cartesian(ylim= c(-0.10, 0.7))  +
  scale_x_continuous(
    breaks = seq(1.5, 8.5, 1),
    labels = \(x) stringr::str_wrap(intervals$interval[x], width = 10)
  )

bp_div

ggsave(plot = bp_div, filename = paste0("figures/", "bp_div.pdf"), 
       width = 174, height = 80, units = "mm")

## 2.4 Sampling - Violinplots -------------------------------------------
 
bp_samp <- ggplot() +
  geom_violin(
    data = traces %>% filter(param == "samplingRate"),
    aes(x = epoch_time, y = value, fill = ansys_pretty, group = interaction(epoch_time, ansys_pretty)),
    alpha = 0.5, color = "transparent",
    position = position_dodge(0.9),
    scale = "width"
  ) +
  geom_errorbar(
    data = trace_summary %>% filter(param == "samplingRate"),
    aes(x = epoch_time, ymin = l95, ymax = h95, color = ansys_pretty),
    position = position_dodge(0.9), linewidth = 0.5, width = 0.8
  ) +
  geom_point(
    data = trace_summary %>% filter(param == "samplingRate"),
    aes(x = epoch_time, y = median, color = ansys_pretty),
    position = position_dodge(0.9), size = 1
  ) +
  annotate(
    "rect",
    xmin = c(1.5, 3.5, 5.5, 7.5),
    xmax = c(2.5, 4.5, 6.5, 8.5),
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.1, colour = NA
  ) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2, name = "Analysis:") +
  scale_fill_manual(values = pal_dino2, name = "Analysis:") +
  scale_x_continuous(
    breaks = seq(1.5, 8.5, 1),
    labels = \(x) stringr::str_wrap(intervals$interval[x], width = 10)
  ) +
  labs(x = "Interval", y = "Sampling rate") +
  coord_cartesian(ylim = c(0, 0.01)) 

bp_samp
ggsave(plot = bp_samp, filename = paste0("figures/", "bp_samp.pdf"), 
       width = 174, height = 80, units = "mm")


## 2.5 sampling GLM parameters ------------------------------------------------------

### baseline

bp_glmbaseline <- traces %>% 
  filter(param %in% c("samplingRateGLM_baseline", "scaleFactorGLM")) %>%
  # bind_rows(priors %>% filter(parameter == "GLM Baseline") %>% unnest(samples) %>%
  #             mutate(dataset = NA, analysis = NA, value = samples)) %>%
  ggplot(aes(x = ansys_pretty, y = value)) +
  geom_violin(aes(fill = ansys_pretty), alpha = 0.3, color = "transparent", position = position_dodge(0.9)) +
  geom_boxplot(aes(color = ansys_pretty, fill = ansys_pretty), position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2) +
  scale_fill_manual(values = pal_dino2) 

bp_glmbaseline

ggsave(plot = bp_glmbaseline, filename = paste0("figures/", "bp_glmbaseline.pdf"), 
       width = 85, height = 65, units = "mm")


### coefficients
# predictor_pretty_t <- tibble(predictor = unique(str_split(traces %>% 
#                                                             filter(param == "samplingRateGLM_coefficientON") %>% 
#                                                             pull(parameter), "\\.", simplify = T)[,2]),
#                              predictor_pretty = c("Collections", "Collections (perm1)", "Collections (perm2)", "Collections (rev)",
#                                                   "Formations", "Formations (perm1)", "Formations (perm2)", "Formations (rev)",
#                                                   "Grid cells", "Grid cells (perm1)", "Grid cells (perm2)", "Grid cells (rev)")) %>%
#   bind_rows(tibble(predictor = unique(str_split(traces %>% 
#                                         filter(param == "samplingRateGLM_coefficientON") %>% 
#                                         pull(parameter), "\\.", simplify = T)[,2]),
#          predictor_pretty = c("Collections",
#                               "Formations",
#                               "Grid cells"))) %>% distinct() %>%
#   mutate(predictor_pretty = forcats::fct_rev(predictor_pretty))

predictor_pretty_t <- tibble(predictor = unique(str_split(traces %>% 
                                                            filter(param == "samplingRateGLM_coefficientON") %>% 
                                                            pull(parameter), "\\.", simplify = T)[,2])) %>%
  mutate(
                             predictor_pretty = str_replace(predictor, "collectionCounts", "Collections"),
                             predictor_pretty = str_replace(predictor_pretty, "formationCounts", "Formations"),
                             predictor_pretty = str_replace(predictor_pretty, "gridCellCounts", "Spatial Grid Cells"),
                             predictor_pretty = str_replace(predictor_pretty, "_r1", " (perm1)"),
                             predictor_pretty = str_replace(predictor_pretty, "_r2", " (perm2)"),
                             predictor_pretty = str_replace(predictor_pretty, "_rev", " (rev)"),
                             predictor_class = str_replace(predictor_pretty, " \\(perm1\\)| \\(perm2\\)| \\(rev\\)", ""))
                            

traces_GLM <- traces %>% 
  filter(grepl("samplingRateGLM_coefficientON|samplingRateGLM_indicator", param)) %>%
  mutate(predictor = str_split(parameter, "\\.", simplify = T)[,2]) %>%
  select(Sample, value, param, predictor, ansys, ansys_pretty) %>% 
  left_join(predictor_pretty_t) %>%
  pivot_wider(values_from = value, names_from = param) 


included <- traces_GLM %>% 
  group_by(predictor, ansys) %>%
  summarise(p = sum(samplingRateGLM_indicator) / n(), .groups = "drop") %>% 
  mutate(inc = p > 0.55) 

trace_GLMsummary <- traces_GLM %>%
  filter(samplingRateGLM_indicator == 1) %>%
  mutate(value = samplingRateGLM_coefficientON) %>%
  group_by(predictor, predictor_pretty, ansys, ansys_pretty) %>%
  summarise(
    median = median(value),
    l95 = HDInterval::hdi(value)[1],
    h95 = HDInterval::hdi(value)[2],
    .groups = "drop")

ansys_levels <- c("non-GLM analysis", "sampling-GLM", "sampling-GLM 12p", "sampling-GLM + error")
pred_levels <- rev(levels(predictor_pretty_t$predictor_pretty))

bp_glmcoef <- traces_GLM %>% 
  filter(samplingRateGLM_indicator == 1) %>%
  mutate(value = samplingRateGLM_coefficientON,
         fac = ifelse(grepl("\\(", predictor_pretty), paste(predictor_class, "Control Predictors"), predictor_pretty),
         var = as.factor(case_when(
           fac == predictor_pretty & ansys == "GLM_sumprior_3p" ~ 1,
           fac == predictor_pretty & ansys == "GLM_sumprior_12p" ~ 2,
           fac == predictor_pretty & ansys == "GLM_sumprior_errordistr_3p" ~ 3,
           fac != predictor_pretty & str_split(predictor, "_", simplify = T)[,2] == "r1" ~ 1,
           fac != predictor_pretty & str_split(predictor, "_", simplify = T)[,2] == "r2" ~ 2,
           fac != predictor_pretty & str_split(predictor, "_", simplify = T)[,2] == "rev" ~ 3,
         ))
         ) %>%
  left_join(included) %>%
  ggplot(aes(
    x = value,
    # y = var
    y = forcats::fct_relevel(var, c("2", "3", "1"))
  )) +
  geom_violin(
    aes(fill = forcats::fct_relevel(ansys_pretty, rev(ansys_levels))),
    alpha = 0.5, color = "transparent",
    position = position_dodge(0.9), scale = "width"
  ) +
  geom_boxplot(
    aes(
      color = forcats::fct_relevel(ansys_pretty, rev(ansys_levels)),
      fill = forcats::fct_relevel(ansys_pretty, rev(ansys_levels))
    ),
    position = position_dodge(0.9),
    alpha = 0.5, width = 0.1, outlier.shape = NA
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2, breaks = ansys_levels) +
  scale_fill_manual(values = pal_dino2, breaks = ansys_levels) +
  # facet_wrap(
  #   ~ forcats::fct_relevel(predictor_pretty, pred_levels),
  #   ncol = 1
  # ) +
  facet_wrap(~ fac, ncol = 1) +
  coord_cartesian(xlim = c(-1.5, 1.5))

bp_glmcoef

prior_coef <- priors %>%
      filter(parameter == "GLM Coefficients") %>%
      unnest(samples) %>%
  ggplot(aes(samples, y = "Prior")) +
  geom_violin(fill = "grey60", alpha = 0.5, color = "transparent",
              position = position_dodge(0.9), scale = "width") + 
  geom_boxplot(color = "grey60", alpha = 0.5, width = 0.1, outlier.shape = NA) +
  theme_dino2()
  
prior_coef

### indicators

post_from_BF <- function(BF, lambda = 0.693, K) {
  p0 <- lambda / K
  (BF * p0) / (1 - p0 + BF * p0)
}

bp_glmind <- traces_GLM %>% 
  group_by(predictor_pretty, ansys_pretty, ansys, predictor, predictor_class) %>%
  summarise(p = sum(samplingRateGLM_indicator) / n(), .groups = "drop") %>%
  # mutate(predictor_pretty = factor(predictor_pretty, levels = pred_levels)) %>%
  # bind_rows(tibble(ansys_pretty = "GLM fix 3p", predictor_pretty = "Prior")) %>%
  mutate(fac = ifelse(grepl("\\(", predictor_pretty), paste(predictor_class, "Control Predictors"), predictor_pretty),
         var = as.factor(case_when(
           fac == predictor_pretty & ansys == "GLM_sumprior_3p" ~ 1,
           fac == predictor_pretty & ansys == "GLM_sumprior_12p" ~ 2,
           fac == predictor_pretty & ansys == "GLM_sumprior_errordistr_3p" ~ 3,
           fac != predictor_pretty & str_split(predictor, "_", simplify = T)[,2] == "r1" ~ 1,
           fac != predictor_pretty & str_split(predictor, "_", simplify = T)[,2] == "r2" ~ 2,
           fac != predictor_pretty & str_split(predictor, "_", simplify = T)[,2] == "rev" ~ 3,
         ))
  ) %>% 
  ggplot(aes(x = p, y = forcats::fct_relevel(var, c("2", "3", "1")))) +
  geom_bar(aes(fill = ansys_pretty), color = "transparent",
           position = position_dodge(1), stat = "identity", width = 0.6) +
  geom_vline(xintercept = post_from_BF(10, K = 7), linetype = 2) +
  geom_vline(xintercept = post_from_BF(100, K = 7), linetype = 2) +
  theme_dino2() + 
  scale_color_manual(values = pal_dino2) +
  scale_fill_manual(values = pal_dino2) +
  facet_wrap(~fac, ncol = 1) +
  scale_x_continuous(breaks = c(0,0.5,1))

bp_glmind

library(patchwork)
cp_glm <- (bp_glmcoef + theme(legend.position = "left") + 
  bp_glmind + theme(legend.position = "none",
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank())) +   plot_layout(widths = c(2.5, 1))
cp_glm
ggsave(plot = cp_glm, filename = paste0("figures/", "cp_glm_ggplot.pdf"), 
       width = 80, height = 140, units = "mm")
ggsave(plot = prior_coef, filename = paste0("figures/", "prior_cglm.pdf"), 
       width = 40, height = 30, units = "mm")

## 2.6 Tree likelihood parameter -----------------------------------------------
# 
# dp_treelikelihood <- traces %>% 
#   filter(param == "likelihood") %>%
#   ggplot(aes(x = value)) +
#   geom_density(aes(fill = ansys_pretty), alpha = 0.5, color = "transparent") + 
#   theme_dino2() +
#   scale_color_manual(values = pal_dino2) +
#   scale_fill_manual(values = pal_dino2) 
# 
# dp_treelikelihood

bp_treelikelihood <-  traces %>% 
  filter(param == "likelihood" | param == "FBD") %>%
  ggplot(aes(x = ansys_pretty, y = value)) +
  geom_violin(aes(fill = ansys_pretty), alpha = 0.5, color = "transparent", position = position_dodge(0.9)) +
  geom_boxplot(aes(color = ansys_pretty, fill = ansys_pretty), position = position_dodge(0.9), 
               alpha = 0.7, width = 0.1, outlier.shape = NA) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2) +
  scale_fill_manual(values = pal_dino2) + coord_flip() +
  theme(legend.position = "none")

bp_treelikelihood
ggsave(plot = bp_treelikelihood, filename = paste0("figures/", "bp_treelikelihood.pdf"), 
       width = 85, height = 75, units = "mm")


## 2.6 Origin parameter -----------------------------------------------

bp_origin <-  traces %>% 
  filter(param == "originFBD") %>%
  # bind_rows(priors %>% filter(parameter == "Origin") %>% unnest(samples) %>%
  #             mutate(ansys_pretty = "Prior", value = samples)) %>%
  ggplot(aes(x = ansys_pretty, y = value)) +
  geom_violin(aes(fill = ansys_pretty), alpha = 0.5, color = "transparent", position = position_dodge(0.9)) +
  geom_boxplot(aes(color = ansys_pretty, fill = ansys_pretty), position = position_dodge(0.9), 
               alpha = 0.7, width = 0.1, outlier.shape = NA) +
  theme_dino2() +
  scale_color_manual(values = pal_dino2) +
  scale_fill_manual(values = pal_dino2) + coord_flip() +
  theme(legend.position = "none")

bp_origin

ggsave(plot = bp_origin, filename = paste0("figures/", "bp_origin.pdf"), 
       width = 85, height = 75, units = "mm")

# 3. Tree files -----------------------------------------------------------

trees_l <- lapply(ansys_names, function(ansys) {
  trees_file <-  paste0("results/Dinosaurs_", ansys, "_combined.trees")
  trees <- read.beast(trees_file) })

get_node_age <- function(tree, node_id) {
  node_times <- node.depth.edgelength(tree)
  root_age <- max(node_times)
  return(root_age - node_times[node_id])
}

# Node ages for Dinosauria and Avialae MRCA nodes
df_nodeages_l <- lapply(trees_l, function(trees) {
  
  df_nodeages_trees <- lapply(1:length(trees), function(i) {
    
    tree <- trees[[i]]@phylo
    df <- tibble(node_id_dinosauria_i = getMRCA(tree, 
                                                c("Eoraptor_lunensis", "Austroraptor_cabazai")),
                 node_id_dinosauria_ii = getMRCA(tree, 
                                                 c("Nyasasaurus_parringtoni", "Austroraptor_cabazai")),

                 node_id_avialae_i = getMRCA(tree, 
                                             c("Archaeopteryx_lithographica", "Austroraptor_cabazai")),
                 node_id_avialae_ii = getMRCA(tree, 
                                              c("Aurornis_xui", "Austroraptor_cabazai")),
                 node_id_root = getMRCA(tree, 
                                   c("Nyasasaurus_parringtoni", "Eocursor_parvus")),
                 age_dinosauria_i = get_node_age(tree, node_id_dinosauria_i),
                 age_dinosauria_ii =  get_node_age(tree, node_id_dinosauria_ii),
                 age_avialae_i =  get_node_age(tree, node_id_avialae_i),
                 age_avialae_ii =  get_node_age(tree, node_id_avialae_ii),
                 age_root =  get_node_age(tree, node_id_root),
                 tree = i,
                 file = trees[[i]]@file)
    return(df)

  })
  df_ansys <- bind_rows(df_nodeages_trees)
})

df_nodeages <- bind_rows(df_nodeages_l) 

df_nodeages_long <- df_nodeages %>% 
  pivot_longer(
    cols = matches("^(node_id|age)_"),
    names_to = c(".value", "node"),
    names_pattern = "(node_id|age)_(.*)"
  ) %>%
  mutate(age_ma = age + 66,
         ansys_pretty = factor(case_when(
           file == "results/Dinosaurs_no_GLM_combined.trees"  ~ "non-GLM analysis",
           file == "results/Dinosaurs_GLM_sumprior_3p_combined.trees" ~ "sampling-GLM",
           file == "results/Dinosaurs_GLM_sumprior_12p_combined.trees"    ~ "sampling-GLM 12p",
           file ==  "results/Dinosaurs_GLM_sumprior_errordistr_3p_combined.trees" ~ "sampling-GLM + error"),
           levels = c("non-GLM analysis", "sampling-GLM", "sampling-GLM + error", "sampling-GLM 12p")),
         node = factor(node, levels = c("dinosauria_i" , "dinosauria_ii" ,"avialae_i" ,    "avialae_ii" ,   "root")))


df_nodeages_summary <- df_nodeages_long %>%
  group_by(node, ansys_pretty, file) %>%
  summarise(median_age = median(age_ma),
            l95 = HDInterval::hdi(age_ma)[1],
            h95 = HDInterval::hdi(age_ma)[2], .groups = "drop") 

# Histogram plots
hist_nodeages <- ggplot(df_nodeages_long %>% filter(node != "root", grepl("ii", node))) +
  geom_histogram(aes(age_ma, fill = ansys_pretty, y=..density..), position = "identity",
                 alpha = 0.5, color = "white", binwidth=1) +
 geom_errorbar(data = df_nodeages_summary %>% filter(node != "root", grepl("ii", node)), aes(y = 0.05, xmin = l95, xmax = h95,
                                               color = ansys_pretty),
               position = position_dodge(1), linewidth = 0.5, width = 0.03) +
 geom_point(data = df_nodeages_summary %>% filter(node != "root", grepl("ii", node)), aes(y = 0.05, x = median_age,
                                            color = ansys_pretty),
 size = 1) +
  scale_x_reverse() +
  theme_dino2() +
  scale_fill_manual(values = pal_dino2) +
  scale_color_manual(values = pal_dino2) +
  #facet_grid(ansys_pretty~node, scales = "free_x") +
  facet_wrap(node ~ ansys_pretty, scales = "free_x", ncol = 1) +
  xlab("Time before present (Ma)") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

hist_nodeages
ggsave(plot = hist_nodeages, filename = paste0("figures/", "hist_nodeages.pdf"), 
       width = 70, height = 70, units = "mm")

hist_nodeages <- ggplot(df_nodeages_long %>% filter(node != "root", grepl("ii", node))) +
  geom_histogram(aes(age_ma, fill = ansys_pretty, y=..density..), position = "identity",
                 alpha = 0.5, color = "white", binwidth=1) +
 geom_errorbar(data = df_nodeages_summary %>% filter(node != "root", grepl("ii", node)), aes(y = 0.05, xmin = l95, xmax = h95,
                                               color = ansys_pretty),
               position = position_dodge(1), linewidth = 0.5, width = 0.03) +
 geom_point(data = df_nodeages_summary %>% filter(node != "root", grepl("ii", node)), aes(y = 0.05, x = median_age,
                                            color = ansys_pretty),
 size = 1) +
  scale_x_reverse() +
  theme_dino2() +
  scale_fill_manual(values = pal_dino2) +
  scale_color_manual(values = pal_dino2) +
  #facet_grid(ansys_pretty~node, scales = "free_x") +
  facet_wrap(node ~ ansys_pretty, scales = "free_x", ncol = 1) +
  xlab("Time before present (Ma)") +
  theme(legend.position = "none",
        axis.text.y = element_blank())

cp_main <- (((bp_div + ylab("") + xlab("") + theme(legend.position = "none")) / 
               (bp_samp + ylab("") + xlab("") +  theme(legend.position = "none")) ) | 
              (hist_nodeages+ ylab("") + theme(
                strip.background = element_blank(),
                strip.text = element_blank(),
                panel.spacing = unit(0, "lines")
              ))) + plot_layout(widths = c(3.5, 1))
ggsave(plot = cp_main, filename = paste0("figures/", "cp_main_ggplot.pdf"), 
       width = 174, height = 130, units = "mm")


hist_nodeages_all <- ggplot(df_nodeages_long) +
  geom_histogram(aes(age_ma, fill = ansys_pretty, y=..density..), position = "identity",
                 alpha = 0.5, color = "white", binwidth=1) +
  geom_errorbar(data = df_nodeages_summary, aes(y = 0.05, xmin = l95, xmax = h95,color = ansys_pretty),
                position = position_dodge(1), linewidth = 0.5, width = 0.03) +
  geom_point(data = df_nodeages_summary, aes(y = 0.05, x = median_age,olor = ansys_pretty),
             size = 1) +
  scale_x_reverse() +
  theme_dino2() +
  scale_fill_manual(values = pal_dino2) +
  scale_color_manual(values = pal_dino2) +
  facet_grid(ansys_pretty~node, scales = "free_x") +
  xlab("Time before present (Ma)") +
  theme(legend.position = "none")

hist_nodeages_all

ggsave(plot = hist_nodeages_all, filename = paste0("figures/", "hist_nodeages_all.pdf"), 
       width = 174, height = 100, units = "mm")

# Tip ages ----------------------------------------------------------------

df_tipages_l <- lapply(trees_l, function(trees) {
  
  df_tipages_trees <- lapply(1:length(trees), function(i) {
    
    tree <- trees[[i]]@phylo

    tibble(
      tip_id = 1:Ntip(tree),
      tip = tree$tip.label,
      age = sapply(1:Ntip(tree), function(x) get_node_age(tree, x)),
      tree = i,
      file = trees[[i]]@file
    )
    
  })
  
  bind_rows(df_tipages_trees)
})

df_tipages <- bind_rows(df_tipages_l)  %>%
  mutate(
    age_ma = age + 66,
    ansys_pretty = factor(case_when(
      file == "results/Dinosaurs_no_GLM_combined.trees" ~ "non-GLM analysis",
      file == "results/Dinosaurs_GLM_fixedon_combined.trees" ~ "GLM fix 3p",
      file == "results/Dinosaurs_GLMall_sumprior_combined.trees" ~ "GLM 12p",
      file == "results/Dinosaurs_GLMall_sumprior_error_combined.trees" ~ "GLM 12p + error"
    ),
    levels = c("non-GLM analysis", "GLM fix 3p", "GLM 12p", "GLM 12p + error"))
  )

df_tipages_summary <- df_tipages %>%
  group_by(tip, ansys_pretty, file) %>%
  summarise(
    median_age = median(age_ma),
    l95 = HDInterval::hdi(age_ma)[1],
    h95 = HDInterval::hdi(age_ma)[2],
    width = h95 - l95,
    .groups = "drop"
  )

df_width_wide <- df_tipages_summary %>%
  select(tip, ansys_pretty, width) %>%
  pivot_wider(names_from = ansys_pretty, values_from = width)

wilcox.test(df_width_wide$`non-GLM analysis`,
            df_width_wide$`GLM fix 3p`,
            paired = TRUE)

wilcox.test(df_width_wide$`non-GLM analysis`,
            df_width_wide$`GLM 12p`,
            paired = TRUE)

wilcox.test(df_width_wide$`non-GLM analysis`,
            df_width_wide$`GLM 12p + error`,
            paired = TRUE)

df_width_wide <- df_width_wide %>%
  mutate(
    reduction_GLM12 = (`non-GLM analysis` - `GLM 12p`) / `non-GLM analysis`,
    reduction_GLM12error = (`non-GLM analysis` - `GLM 12p + error`) / `non-GLM analysis`,
    reduction_fix3 = (`non-GLM analysis` - `GLM fix 3p`) / `non-GLM analysis`
  )
df_width_wide %>%
  summarise(
    median_reduction12 = median(reduction_GLM12, na.rm = TRUE),
    mean_reduction12 = mean(reduction_GLM12, na.rm = TRUE),
    median_reduction12error = median(reduction_GLM12error, na.rm = TRUE),
    mean_reduction12error = mean(reduction_GLM12error, na.rm = TRUE),
    median_reduction3 = median(reduction_fix3, na.rm = TRUE),
    mean_reduction3 = mean(reduction_fix3, na.rm = TRUE)
  )

ggplot(df_tipages_summary, aes(x = ansys_pretty, y = width)) +
  geom_violin() +
  theme_bw() +
  labs(y = "95% HPD width (Ma)", x = "")

ggplot(df_tipages_summary, aes(x = ansys_pretty, y = width, group = tip)) +
  geom_line(alpha = 0.2) +
  stat_summary(fun = median, geom = "line", color = "red", size = 1.2, group = 1) +
  theme_bw()


# library(ggridges)
# ggplot(df_nodeages_long,
#        aes(x = age_ma,
#            y = node,
#            fill = file)) +
#   geom_density_ridges(
#     alpha = 0.6,
#     scale = 1.2
#   ) +
#   scale_x_reverse() +
#   labs(
#     x = "Time before present (Ma)",
#     y = NULL
#   ) +
#   theme_dino() +
#   scale_fill_manual(values = pal_dino2)
# 
# ggplot(df_nodeages_long, 
#        # %>% filter(node %in% c("avialae_i", "dinosauria_i", "root")), 
#        aes(x = age_ma,
#            y = interaction(file,node),
#            group = interaction(file,node),
#            fill = file)) +
#   geom_density_ridges(
#     alpha = 0.6,
#     scale = 3
#   ) +
#   scale_x_reverse() +
#   labs(
#     x = "Time before present (Ma)",
#     y = NULL
#   ) +
#   theme_dino() +
#   scale_fill_manual(values = pal_dino2) #+
#   #facet_wrap(~node, ncol = 1)
# # All
# bp_nodeages_all <- ggplot(df_nodeages_long) +
#   geom_violin(aes(age_ma, y = ansys_pretty, fill = ansys_pretty), alpha = 0.3, color = "transparent") +
#   geom_errorbar(data = df_nodeages_summary, aes(y = ansys_pretty, xmin = l95, xmax = h95, 
#                                                      color = ansys_pretty),
#                 position = position_dodge(1), linewidth = 0.5, width = 0.4) +
#   geom_point(data = df_nodeages_summary, aes(y = ansys_pretty, x = median_age, 
#                                                   color = ansys_pretty),
#              position = position_dodge(1), size = 1) +
#   facet_wrap(~node, ncol = 2, scales = "free") +
#   scale_x_reverse() +
#   theme_dino2() +
#   scale_fill_manual(values = pal_dino2) +
#   scale_color_manual(values = pal_dino2) + 
#   annotate(geom = "rect", xmin = c(160, 200, 240),
#            xmax = c(180, 220, 260),
#            ymin = -Inf, ymax = Inf,
#            colour = "grey", alpha = 0.1, linewidth = 0)  
# 
# bp_nodeages_all


