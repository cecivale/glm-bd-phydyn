# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: Early SARS-CoV-2 Europe Phylodynamics
#        V\ Y /V    Figures
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

# 0. Libraries -----------------------------------------------------------------
library(tidyverse)
library(data.table)
library(treedataverse)
library(ggridges)
library(ggimage)
library(HDInterval)
library(ggsankey)
library(patchwork)


source("workflow/scripts/utils.R")
source("workflow/scripts/priors.R")
source("workflow/scripts/plot_opts.R")


# ------------------------------------------------------------------------------
# To run outside of the snakemake workflow, run this code first 
# ------------------------------------------------------------------------------

# setClass(
#   "snakemake_object",
#    contains= "tbl_df",
#    slots = c(input = "list", output = "character", params = "character")
#  )

# # Chains from same analysis, different data replicates
# ansys_data = c("mGLM_6p7e_sumprior_lowerboundFR0/eubdmm_open.","no_GLM/eubdmm_open.")
# data_seed = c(0,1,2)
# snakemake <- new("snakemake_object", tibble(),
#                  input = list(cases = "resources/ext_ecdccases.tsv",
#                            trace_files = unlist(lapply(ansys_data, function(a) paste0("results/analysis/", a, data_seed, ".log"))),
#                            typed_node_ccd0_tree_files =  unlist(lapply(ansys_data, function(a) paste0("results/analysis/", a, data_seed, ".CCD0.typed.node.tree"))),
#                            typed_node_mcc_tree_files =  unlist(lapply(ansys_data, function(a) paste0("results/analysis/", a, data_seed, ".MCC.typed.node.tree"))),
#                            traj_files =   unlist(lapply(ansys_data, function(a) paste0("results/analysis/", a, data_seed, ".traj")))),
#                  output = c(out_folder = "results/figures/"))

# ------------------------------------------------------------------------------

# 1. Analysis specs ------------------------------------------------------------
demes <- c(CN = "China", FR = "France", DE = "Germany", IT = "Italy", OE = "OtherEU")

analysis_from <- "2019-12-01"
analysis_to <- "2020-03-08"

cases <- get_ECDCcases(snakemake@input[["cases"]], demes)


# 2. Trace files ---------------------------------------------------------------
traces_l <- lapply(snakemake@input[["trace_files"]], function(trace_file) {
  trace_wide <- read_trace(trace_file, burninFrac = 0)  %>%
  select(-starts_with("X"))
  
  trace_long <- trace_wide %>%
    pivot_longer(-1, names_to = "parameter", values_to = "value") %>%
    mutate(param = str_split(parameter, "_", n = 2, simplify = T)[, 1],
           deme = str_split(parameter, "_", n = 4, simplify = T)[, 2],
           deme_to = str_split(parameter, "_", n = 4, simplify = T)[, 4],
           epoch = case_when(
             grepl("SV", param) ~ str_split(param, "SV", n = 2, simplify = T)[, 2],
             grepl("SM", param) ~ str_split(param, "SM", n = 2, simplify = T)[, 2]),
           param = case_when(
             grepl("SV", param) ~ str_split(param, "SV", n = 2, simplify = T)[, 1],
             grepl("SM", param) ~ str_split(param, "SM", n = 2, simplify = T)[, 1],
             T ~ parameter)) %>%
    filter(deme != "endtime") %>% mutate(dataset = str_replace(str_split(trace_file, "/", simplify = T)[,4], ".log", ""),
                                         analysis = str_split(trace_file, "/", simplify = T)[,3])
})

traces <- bind_rows(traces_l) %>%
  mutate(analysis_pretty = factor(case_when(
    analysis == "no_GLM" ~ "non-GLM analysis",
    analysis == "mGLM_6p7e_sumprior_lowerboundFR0" ~ "migration-GLM analysis",
    analysis == "mGLM_6p7e_sumprior_spGLMid" ~ "migration and sampling-GLM analysis"),
    levels = c("non-GLM analysis", "migration-GLM analysis", "migration and sampling-GLM analysis" )),
    dataset_pretty = paste0("S", str_split(dataset, "\\.", simplify = T)[,2])
  )

trace_summary <- traces %>%
  group_by(dataset, analysis, dataset_pretty, analysis_pretty,
           parameter, param, deme, deme_to, epoch) %>%
  summarise(median = median(value),
            l95 = HDInterval::hdi(value)[1],
            h95 = HDInterval::hdi(value)[2], .groups = "drop")

## 2.1 Re - Boxplots -----------------------------------------------------------

bp_re <- traces %>% 
  filter(param == "Re", epoch == "i1" | (epoch == "i0" & deme == "China")) %>%
  mutate(var = paste(deme, epoch, dataset_pretty)) %>%
  bind_rows(priors %>% filter(parameter == "Re") %>% unnest(samples) %>%
              mutate(var = "Prior", value = samples, deme = "All") %>%
              cross_join(select(traces, analysis_pretty) %>% distinct())) %>%
  ggplot(aes(x = var, y = value)) +
  geom_violin(aes(fill = deme), alpha = 0.4, color = "transparent", 
              position = position_dodge(0.9), scale = "width") +
  geom_boxplot(aes(color = deme, fill = deme), position = position_dodge(0.9), 
               alpha = 0.6, width = 0.1, outlier.shape = NA) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2) +
  scale_fill_manual(values = pal_eubdmm2) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_cartesian(ylim = c(0, 6)) +
  facet_wrap(~analysis_pretty, ncol = 1) +
  scale_x_discrete(labels = function(x) sub(".* ", "", x))

bp_re

cowplot::ggsave2(bp_re, file = paste0(snakemake@output, "bp_re_ggplot.pdf"),
                 width = 178, height = 80, units = "mm")

## 2.2 Origin -------------------------------------------------------------------

dp_origin <- traces %>% 
  filter(param == "originBDMMPrime") %>%
  mutate(value_date = to_date(value, mrs = analysis_to)) %>%
  bind_rows(priors %>% filter(parameter == "Origin") %>% unnest(samples) %>%
              mutate(dataset_pretty = "Prior", value = samples,
                     value_date = to_date(value, mrs = analysis_to)) %>% 
              cross_join(select(traces, analysis_pretty) %>% distinct())) %>%
  ggplot(aes(x = value_date)) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  facet_wrap(~analysis_pretty, ncol = 1) +
  scale_fill_manual(values = c( "grey50", dataset_pal)) +
  scale_color_manual(values = c( "transparent", dataset_pal)) +
  scale_x_date(date_breaks = "month", date_labels = "%b %Y", 
               limits = c(ymd("2019-07-01"),ymd("2020-01-05")))

dp_origin 

cowplot::ggsave2(dp_origin, file = paste0(snakemake@output, "dp_origin_ggplot.pdf"),
                 width = 150, height = 70, units = "mm")

 
## 2.3 Sampling proportions -----------------------------------------------------

bp_sampling <- traces %>% 
  filter(param == "samplingProportion", epoch == "i2", analysis != "mGLM_6p7e_sumprior_spGLMid") %>%
  mutate(var = paste(epoch, deme, dataset_pretty)) %>%
  #bind_rows(priors %>% filter(parameter == "Sampling proportion") %>% unnest(samples) %>%
  #            mutate(var = "Prior", value = samples, deme = "All")) %>%
  ggplot(aes(x = var, y = value)) +
  geom_violin(aes(fill = deme), alpha = 0.3, color = "transparent", 
              position = position_dodge(0.9), scale = "width") +
  geom_boxplot(aes(color = deme, fill = deme), position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2) +
  scale_fill_manual(values = pal_eubdmm2) +
  facet_wrap(~analysis_pretty, ncol = 1)  +
  scale_x_discrete(labels = function(x) sub(".* ", "", x))

bp_sampling

cowplot::ggsave2(bp_sampling, file = paste0(snakemake@output, "bp_sampling_ggplot.pdf"),
                 width = 178, height = 80, units = "mm")

bp_sampling_spglm <- traces %>% 
  filter(param == "samplingProportion", analysis == "mGLM_6p7e_sumprior_spGLMid") %>%
  mutate(var = paste(dataset, epoch)) %>%
  #bind_rows(priors %>% filter(parameter == "Sampling proportion") %>% unnest(samples) %>%
  #            mutate(var = "Prior", value = samples, deme = "All")) %>%
  ggplot(aes(x = var, y = value)) +
  geom_violin(aes(fill = deme), alpha = 0.3, color = "transparent", 
              position = position_dodge(0.9), scale = "width") +
  geom_boxplot(aes(color = deme, fill = deme), position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2) +
  scale_fill_manual(values = pal_eubdmm2) +
  facet_wrap(~analysis_pretty + deme + dataset_pretty, ncol = 3, scales = "free")  +
  scale_x_discrete(labels = function(x) sub(".* ", "", x)) 

bp_sampling_spglm

cowplot::ggsave2(bp_sampling_spglm, file = paste0(snakemake@output, "bp_sampling_spglm_ggplot.pdf"))

## 2.4 Sampling GLM parameters -------------------------------------------------

bp_spgglmcoef <- traces %>% 
  filter(analysis == "mGLM_6p7e_sumprior_spGLMid", 
         grepl("samplingProportionValues_coefficientON", param)) %>%
  mutate(deme = demes[str_split(param, "_", simplify = T)[,3]]) %>%
  # bind_rows(priors %>% filter(parameter == "GLM Coefficients") %>% unnest(samples) %>%
  #             mutate(dataset = "", analysis = "",  value = samples, predictor = "Prior") %>%
  #             cross_join(select(traces, analysis_pretty) %>% distinct())) %>%
  ggplot(aes(x = interaction(dataset_pretty, deme), y = value)) +
  geom_violin(aes(fill = deme), alpha = 0.5, color = "transparent", 
              position = position_dodge(0.9), scale="width") +
  geom_boxplot(aes(color = deme, fill = deme), position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm) +
  scale_fill_manual(values = pal_eubdmm) 

bp_spgglmcoef

cowplot::ggsave2(bp_spgglmcoef, file = paste0(snakemake@output, "bp_spgglmcoef_ggplot.pdf"))

## 2.5 Migration rates ----------------------------------------------------------

# trace_summary_step <- trace_summary %>%
#   group_by(dataset, analysis, parameter, param, deme, deme_to) %>%
#   arrange(date) %>%
#   group_modify(~ bind_rows(
#     .x,
#     slice_tail(.x, n = 1) %>% mutate(date = analysis_from)
#   )) %>% ungroup()

epoch_dates <- tibble(
  epoch = c("i00","i0","i1","i2","i3","i4","i5", "i6"),
  date = as.Date(c(
    "2019-09-01",
    "2019-10-01",
    "2019-11-01",
    "2019-12-01",
    "2020-01-01",
    "2020-02-01",
    "2020-03-01",
    "2020-03-08"
  ))
)

lp_mig_mglm <- trace_summary %>%   
  bind_rows(trace_summary %>% ungroup() %>%
              filter(param == "migrationRate",analysis == "mGLM_6p7e_sumprior_lowerboundFR0", 
                     deme %in% demes, epoch == "i0") %>% mutate(epoch = "i00")) %>%
  left_join(epoch_dates, by = "epoch") %>%
  filter(param == "migrationRate",analysis == "mGLM_6p7e_sumprior_lowerboundFR0", 
         deme %in% demes) %>%
  mutate(var = paste(epoch, dataset, analysis),
         e = as.numeric(str_replace(epoch, "i", ""))) %>%
  ggplot(aes(x = date, y = median, color = dataset_pretty, group = dataset_pretty)) +
  pammtools::geom_stepribbon(aes(ymin = l95, ymax = h95, fill = dataset_pretty, color = dataset_pretty),
                             alpha = 0.2, linetype = 2, linewidth = 0,
                             direction = "vh") +
  geom_step(direction = "vh") +
  theme_eubdmm() +
  scale_fill_manual(values = dataset_pal, name = "Datatet") +
  scale_color_manual(values = dataset_pal, name = "Datatet") +
  facet_grid(deme ~ deme_to) +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") +
  scale_y_log10() 

lp_mig_mglm

cowplot::ggsave2(lp_mig_mglm, file = paste0(snakemake@output, "lp_mig_mglm_ggplot.pdf"),
                 width = 178, height = 150, units = "mm")

lp_mig_spglm <- trace_summary %>% 
  bind_rows(trace_summary %>% ungroup() %>%
              filter(param == "migrationRate",analysis == "mGLM_6p7e_sumprior_spGLMid", 
                     deme %in% demes, epoch == "i0") %>% mutate(epoch = "i00")) %>%
  left_join(epoch_dates, by = "epoch") %>%
  filter(param == "migrationRate",analysis == "mGLM_6p7e_sumprior_spGLMid", 
         deme %in% demes) %>%
  mutate(var = paste(epoch, dataset, analysis),
         e = as.numeric(str_replace(epoch, "i", ""))) %>%
  ggplot(aes(x = date, y = median, color = dataset_pretty, group = dataset_pretty)) +
  pammtools::geom_stepribbon(aes(ymin = l95, ymax = h95, fill = dataset_pretty),
                             alpha = 0.2, linetype = 2, linewidth = 0, 
                             direction = "vh") +
  geom_step(direction = "vh") +
  theme_eubdmm() +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  facet_grid(deme ~ deme_to) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_log10() 

lp_mig_spglm

cowplot::ggsave2(lp_mig_spglm, file = paste0(snakemake@output, "lp_mig_spglm_ggplot.pdf"))

epoch_dates_noglm <- tibble(
  epoch = c("i00", "i0","i1"),
  date = as.Date(c(
    "2019-12-01",
    "2020-01-23",
    "2020-03-08"
  ))
)

lp_mig_noglm <- trace_summary %>% ungroup() %>%
  filter(param == "migrationRate",analysis == "no_GLM", 
         deme %in% demes) %>%
  #for first step in the plot
  bind_rows(trace_summary %>% ungroup() %>%
  filter(param == "migrationRate",analysis == "no_GLM", 
         deme %in% demes, epoch == "i0") %>% mutate(epoch = "i00")) %>%
  mutate(var = paste(epoch, dataset, analysis),
         e = as.numeric(str_replace(epoch, "i", ""))) %>% 
  left_join(epoch_dates_noglm, by = "epoch") %>%
  ggplot(aes(x = date, y = median, color = dataset_pretty, group = dataset_pretty)) +
  pammtools::geom_stepribbon(aes(ymin = l95, ymax = h95, fill = dataset_pretty),
                             alpha = 0.2, linetype = 2, linewidth = 0, 
                             direction = "vh") +
  geom_step(direction = "vh") +
  geom_hline(aes(yintercept = 1), color = "grey60", linetype = 2) +
  geom_ribbon(aes(ymin = 0.141, ymax = 7.099), alpha = 0.1, fill ="grey70", color = "transparent") +
  theme_eubdmm() +
  scale_fill_manual(values = dataset_pal, name = "Datatet") +
  scale_color_manual(values = dataset_pal, name = "Datatet") +
  facet_grid(deme ~ deme_to) +
  scale_x_date(date_breaks = "1 months", date_labels = "%b") +
  scale_y_log10() 

lp_mig_noglm

cowplot::ggsave2(lp_mig_noglm, file = paste0(snakemake@output, "lp_mig_noglm_ggplot.pdf"),
                 width = 178, height = 150, units = "mm")


## 2.6 mig GLM parameters ------------------------------------------------------

### baseline

bp_migglmbaseline <- traces %>% 
  filter(param == "migrationRateValues_baseline") %>%
  bind_rows(priors %>% filter(parameter == "GLM Baseline") %>% unnest(samples) %>%
              mutate(dataset_pretty = "Prior", analysis = NA, value = samples)) %>%
  ggplot(aes(x = dataset_pretty, y = value)) +
  geom_violin(aes(fill = dataset_pretty), alpha = 0.3, color = "transparent",
              position = position_dodge(0.9), scale = "width") +
  geom_boxplot(aes(color = dataset_pretty, fill = dataset_pretty), 
               position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  theme_eubdmm() +
  scale_y_continuous(limits = c(0, 3)) +
  xlab("") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("grey50", dataset_pal), name = "Datatet") +
  scale_color_manual(values = c("grey50", dataset_pal), name = "Datatet") +
  scale_x_discrete(labels = function(x) sub(".* ", "", x))

bp_migglmbaseline

cowplot::ggsave2(bp_migglmbaseline, file = paste0(snakemake@output, "bp_migglmbaseline_ggplot.pdf"),
                 width = 60, height = 60, units = "mm")

### coefficients

traces_migGLM <- traces %>% filter(grepl("migrationRateValues_coefficientON|migrationRateValues_indicator", param)) %>%
  mutate(param = case_when(
    grepl("migrationRateValues_coefficientON", parameter) ~ "migrationRateValues_coefficientON",
    grepl("migrationRateValues_indicator", parameter) ~ "migrationRateValues_indicator"),
         predictor = str_replace(str_split(parameter, "/", simplify = T)[,3], ".tsv", "")) %>%
  select(Sample, value, param, predictor, dataset, analysis) %>% 
  pivot_wider(values_from = value, names_from = param)

included <- traces_migGLM %>% 
  group_by(predictor, dataset, analysis) %>%
  summarise(p = sum(migrationRateValues_indicator) / n(), .groups = "drop") %>% 
  mutate(inc = p > 0.55) 

bp_migglmcoef_mglm <- traces_migGLM %>% 
  filter(migrationRateValues_indicator == 1, analysis == "mGLM_6p7e_sumprior_lowerboundFR0") %>%
  mutate(value = migrationRateValues_coefficientON) %>%
  left_join(included) %>%
  bind_rows(priors %>% filter(parameter == "GLM Coefficients") %>% unnest(samples) %>%
              mutate(dataset = "", analysis = "",  value = samples, predictor = "Prior") %>%
              cross_join(select(traces, analysis_pretty) %>% distinct())) %>%
  ggplot(aes(x = value, y = paste0(predictor, dataset), group = interaction(predictor, dataset, analysis))) +
  geom_violin(aes(fill = interaction(dataset, inc)), alpha = 0.5, color = "transparent", 
              position = position_dodge(0.9), scale="width") +
  geom_boxplot(aes(color = interaction(dataset, inc), fill = interaction(dataset, inc)), position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_eubdmm() +
  scale_color_manual(values = c(grey_pal, dataset_pal)) +
  scale_fill_manual(values = c(grey_pal, dataset_pal)) 

bp_migglmcoef_mglm


bp_migglmcoef_spglm <- traces_migGLM %>% 
  filter(migrationRateValues_indicator == 1, analysis == "mGLM_6p7e_sumprior_spGLMid") %>%
  mutate(value = migrationRateValues_coefficientON) %>%
  left_join(included) %>%
  bind_rows(priors %>% filter(parameter == "GLM Coefficients") %>% unnest(samples) %>%
              mutate(dataset = "", analysis = "",  value = samples, predictor = "Prior") %>%
              cross_join(select(traces, analysis_pretty) %>% distinct())) %>%
  ggplot(aes(x = value, y = paste0(predictor, dataset), group = interaction(predictor, dataset, analysis))) +
  geom_violin(aes(fill = interaction(dataset, inc)), alpha = 0.5, color = "transparent", 
              position = position_dodge(0.9), scale="width") +
  geom_boxplot(aes(color = interaction(dataset, inc), fill = interaction(dataset, inc)), position = position_dodge(0.9), 
               alpha = 0.5, width = 0.1, outlier.shape = NA) +
  geom_vline(aes(xintercept = 0), linetype = 2) +
  theme_eubdmm() +
  scale_color_manual(values = c(grey_pal, pal)) +
  scale_fill_manual(values = c(grey_pal, pal)) 

bp_migglmcoef_spglm


### indicators

post_from_BF <- function(BF, lambda = 0.693, K) {
  p0 <- lambda / K
  (BF * p0) / (1 - p0 + BF * p0)
}


bp_migglmind_mglm <- traces_migGLM %>% 
  filter(analysis == "mGLM_6p7e_sumprior_lowerboundFR0") %>%
  group_by(predictor, dataset, analysis) %>%
  summarise(p = sum(migrationRateValues_indicator) / n()) %>%
  mutate(inc = p > 0.55) %>%
  bind_rows(priors %>% filter(parameter == "GLM Indicator") %>% unnest(samples) %>%
              mutate(dataset = "", analysis = "", value = samples)) %>%
  ggplot(aes(x = p, y = paste0(predictor, dataset), group = interaction(predictor, dataset, analysis))) +
  geom_bar(aes(fill = interaction(dataset, inc)), color = "transparent", 
           position = position_dodge(1), stat = "identity", width = 0.7, alpha = 0.8) +
  geom_vline(aes(xintercept = post_from_BF(10,  K = 7)), linetype = 2) +
  geom_vline(aes(xintercept = post_from_BF(100,  K = 7))) +
  theme_eubdmm() +
  scale_color_manual(values = c(grey_pal, dataset_pal)) +
  scale_fill_manual(values = c(grey_pal, dataset_pal)) 

bp_migglmind_mglm


bp_migglmind_spglm <- traces_migGLM %>% 
  filter(analysis == "mGLM_6p7e_sumprior_spGLMid") %>%
  group_by(predictor, dataset, analysis) %>%
  summarise(p = sum(migrationRateValues_indicator) / n()) %>%
  mutate(inc = p > 0.55) %>%
  bind_rows(priors %>% filter(parameter == "GLM Indicator") %>% unnest(samples) %>%
              mutate(dataset = "", analysis = "", value = samples)) %>%
  ggplot(aes(x = p, y = paste0(predictor, dataset), group = interaction(predictor, dataset, analysis))) +
  geom_bar(aes(fill = interaction(dataset, inc)), color = "transparent", 
           position = position_dodge(1), stat = "identity", width = 0.7, alpha = 0.8) +
  geom_vline(aes(xintercept = post_from_BF(10,  K = 7)), linetype = 2) +
  geom_vline(aes(xintercept = post_from_BF(100,  K = 7))) +
  theme_eubdmm() +
  scale_color_manual(values = c(grey_pal, pal)) +
  scale_fill_manual(values = c(grey_pal, pal)) 

bp_migglmind_spglm


p_epiglm_mglm <- cowplot::plot_grid(bp_migglmcoef_mglm + theme(legend.position = "none"), 
                   bp_migglmind_mglm + theme(legend.position = "none",
                                        axis.title.y = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank()), rel_widths=c(2,1))
cowplot::ggsave2(p_epiglm_mglm, file = paste0(snakemake@output, "p_epiglm_mglm_ggplot.pdf"))


p_epiglm_spglm <- cowplot::plot_grid(bp_migglmcoef_spglm + theme(legend.position = "none"), 
                               bp_migglmind_spglm + theme(legend.position = "none",
                                                    axis.title.y = element_blank(),
                                                    axis.text.y = element_blank(),
                                                    axis.ticks.y = element_blank()), rel_widths=c(2,1))
cowplot::ggsave2(p_epiglm_spglm, file = paste0(snakemake@output, "p_epiglm_spglm_ggplot.pdf"))

## 2.7 Tree likelihood parameters -----------------------------------------------
dp_posterior <- traces %>% 
  filter(param == "posterior") %>%
  ggplot(aes(x = value)) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  facet_wrap(~analysis_pretty, ncol=1) +
  scale_fill_manual(values = dataset_pal, name = "Dataset") +
  scale_color_manual(values = dataset_pal, name = "Dataset") +
  scale_x_continuous(limits = c(-42050, -41700))

dp_posterior
cowplot::ggsave2(dp_posterior, file = paste0(snakemake@output, "dp_posterior_ggplot.pdf"),
                 width = 80, height = 70, units = "mm")


dp_treelikelihood <- traces %>% 
  filter(param == "treeLikelihood") %>%
  ggplot(aes(x = value)) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  scale_fill_manual(values = dataset_pal, name = "Dataset") +
  scale_color_manual(values = dataset_pal, name = "Dataset") +
  scale_x_continuous(limits = c(-42100, -41800)) +
  facet_wrap(~analysis_pretty, ncol=1)

dp_treelikelihood

cowplot::ggsave2(dp_treelikelihood, file = paste0(snakemake@output, "dp_treelikelihood_ggplot.pdf"),
                 width = 80, height = 70, units = "mm") 


dp_gamma <- traces %>% 
  filter(param == "gammaShape") %>%
  ggplot(aes(x = value)) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  facet_wrap(~analysis_pretty, ncol=1) +
  scale_fill_manual(values = dataset_pal, name = "Dataset") +
  scale_color_manual(values = dataset_pal, name = "Dataset") +
  scale_x_continuous(limits = c(0, 0.5))

dp_gamma

cowplot::ggsave2(dp_gamma, file = paste0(snakemake@output, "dp_gamma_ggplot.pdf"),
                 width = 80, height = 70, units = "mm")


dp_kappa <- traces %>% 
  filter(param == "kappa") %>%
  ggplot(aes(x = value)) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  facet_wrap(~analysis_pretty, ncol=1) +
  scale_fill_manual(values = dataset_pal, name = "Dataset") +
  scale_color_manual(values = dataset_pal, name = "Dataset") 

dp_kappa
cowplot::ggsave2(dp_kappa, file = paste0(snakemake@output, "dp_kappa_ggplot.pdf"),
                 width = 80, height = 70, units = "mm")

dp_treeheight <- traces %>% 
  filter(param == "treeBDMMPrime.height") %>%
  ggplot(aes(x = to_date(value, mrs = analysis_to))) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  facet_wrap(~analysis_pretty, ncol=1) +
  scale_fill_manual(values = dataset_pal, name = "Dataset") +
  scale_color_manual(values = dataset_pal, name = "Dataset") +
  scale_x_date(date_labels = "%d %b %Y")

dp_treeheight
cowplot::ggsave2(dp_treeheight, file = paste0(snakemake@output, "dp_treeheight_ggplot.pdf"),
                 width = 80, height = 70, units = "mm")

dp_treelength <- traces %>% 
  filter(param == "treeBDMMPrime.treeLength") %>%
  ggplot(aes(x = value)) +
  geom_density(aes(fill = dataset_pretty, color = dataset_pretty), alpha = 0.5) + 
  theme_eubdmm() +
  facet_wrap(~analysis_pretty, ncol=1) +
  scale_fill_manual(values = dataset_pal, name = "Dataset") +
  scale_color_manual(values = dataset_pal, name = "Dataset") +
  coord_cartesian(xlim = c(6, 11))

dp_treelength
cowplot::ggsave2(dp_treelength, file = paste0(snakemake@output, "dp_treelength_ggplot.pdf"),
                 width = 80, height = 70, units = "mm")



# 3. Summary tree files -----------------------------------------------------------


trees_l <- lapply(list(#snakemake@input[["typed_node_ccd0_tree_files"]],
                       snakemake@input[["typed_node_mcc_tree_files"]]), function (tree_files_l) {
  lapply(tree_files_l, function(tree_file) {
  tree <- read.beast(tree_file)
  
  p <- ggtree(tree, mrsd = analysis_to) +
    geom_range(range = 'height_0.95_HPD',
               color = "grey30", alpha = 0.2, size = 2) +
    geom_nodelab(aes(label = round(as.numeric(posterior), 3)),
                 hjust = 1.3, vjust = -0.5, size = 2, alpha = 0.7) +
    theme_tree2() +
    ggtree::scale_x_ggtree(breaks = decimal_date(ymd(c("2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01"))), 
                           labels = ymd(c("2019-12-01", "2020-01-01", "2020-02-01", "2020-03-01")))
  
  types_p <- tree@data %>% 
    mutate(type.set = purrr::map(type.set, ~ str_trim(as.character(.x))),
           type.set.prob = purrr::map(type.set.prob, ~ as.numeric(.x))
    ) %>%
    select(type.set, type.set.prob, node)  %>%
    unnest(c(type.set, type.set.prob)) %>% #unnest(type.set.prob) %>%
    pivot_wider(names_from = type.set, values_from = type.set.prob, values_fill = 0)
  
  pies <- nodepie(types_p, cols = 2:6)
  pies <- lapply(pies, function(g) g+scale_fill_manual(values = pal_eubdmm2))
  p2 <- p + geom_inset(pies, width = .015, height = .015) 
  
  cowplot::ggsave2(p2, 
                   file = paste0(snakemake@output, str_replace(str_replace(tree_file, "results/analysis/", ""), "/", "_"), "_ggplot.pdf"), 
                   width = 5.9, height = 5.7)
  
})})

# Add genbank ID tip lab and posterior prob in each node

# 4. Read traj files -----------------------------------------------------------

trajs_l <- lapply(snakemake@input[["traj_files"]], function(traj_file) {
  trajs <- read_trajectories_dt(traj_file, burnin = 0, subsample_n = NULL, demes, analysis_to) %>%
    mutate(dataset = str_replace(str_split(traj_file, "/", simplify = T)[,4], ".traj", ""),
           analysis = str_split(traj_file, "/", simplify = T)[,3])
  return(trajs)
})


trajs <-  rbindlist(trajs_l)  
rm(trajs_l)
trajs[, uniqueN(Sample), by = analysis]
trajs[, min(Sample), by = analysis]
trajs[, max(Sample), by = analysis]


trajs_d <- bingrid_trajectories_dt(trajs[variable == "D"], bin = "day")
traj_d_summary <- summarise_bintrajectories_dt(trajs_d)

#trajs_b <- bingrid_trajectories_dt(trajs[variable == "B" & deme != "China"], bin = "day")
trajs_m <- bingrid_trajectories_dt(trajs[variable == "M"], bin = "day")


# Plot 1: Recovered/isolated/dead cases per deme -------------------------------
lptrajs_cases_deme <- ggplot(traj_d_summary %>% 
                               mutate(analysis_pretty = factor(case_when(
  analysis == "no_GLM" ~ "non-GLM",
  analysis == "mGLM_6p7e_sumprior_lowerboundFR0" ~ "migration-GLM"),
  levels = c("non-GLM", "migration-GLM")),
  dataset_pretty = paste0("S", str_split(dataset, "\\.", simplify = T)[,2]))) +
  geom_line(aes(x = date, y = cumvalue_med, color = deme)) +
  geom_ribbon(aes(date, ymin = cumvalue_hpd_low, ymax = cumvalue_hpd_high, 
                  group = interaction(deme, analysis, dataset), fill =deme), 
              alpha = 0.2) +
  geom_line(data = cases,
            aes(date,total_cases, group = 1), linetype = 2, linewidth = 0.5, color = "grey30") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_x_date(limits = as.Date(c("2019-10-01", "2020-03-07"))) +
  facet_grid(analysis_pretty + dataset_pretty ~deme) +
  theme_eubdmm() +
  # scale_color_manual(values = c(pal[1], pal[2], "#DC9000")) +
  # scale_fill_manual(values = c(pal[1], pal[2], "#DC9000"))
  scale_color_manual(values = pal_eubdmm2) +
  scale_fill_manual(values = pal_eubdmm2) +
  theme(legend.position = "none") +
  ylab("Infections at reporting date") +
  xlab("")

lptrajs_cases_deme

cowplot::ggsave2(lptrajs_cases_deme, file = paste0(snakemake@output, "lptrajs_cases_deme_ggplot.pdf"),
                 width = 17.5, height = 15, units = "cm")


# Plot 2: Reported cases and inference all Europe ------------------------------
trajs_d_eu_summary <- trajs_d %>% filter(variable %in% c("D"), deme != "China") %>%
  arrange(analysis, dataset, Sample, date) %>% # sum all european cases
  group_by(date, dataset, analysis, Sample) %>% summarise(cumvalue = sum(cumvalue), .groups = "drop") %>% 
  group_by(date, dataset, analysis) %>%
  summarise(mean = mean(cumvalue),
            median = median(cumvalue),
            l95 = HDInterval::hdi(cumvalue)[1],
            h95 = HDInterval::hdi(cumvalue)[2]) %>%
  ungroup() %>% make_pretty()

trajs_d_eu_summary %>% group_by(analysis, dataset) %>% arrange(desc(date)) %>% slice(1)
cases %>% filter(deme!='China') %>% group_by(date) %>% summarise(total_cases = sum(total_cases)) %>%
  arrange(desc(date)) %>% slice(1)
                                                                 
lptrajs_cases_allansys <- ggplot(trajs_d_eu_summary) +
  geom_line(aes(x = date, y =  median, group = interaction(analysis,dataset), color = analysis_pretty)) +
  geom_ribbon(aes(date, ymin =  l95, ymax =  h95, group = interaction(analysis,dataset), fill = analysis_pretty), 
              alpha = 0.3, linetype = 2) +
  geom_line(data = cases %>% filter(deme!='China') %>% group_by(date) %>% summarise(total_cases = sum(total_cases)),
            aes(date, total_cases, group = 1), linetype = 2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_date(limits = as.Date(c("2019-12-01", "2020-03-07"))) +
  theme_eubdmm() +
  facet_wrap(~dataset_pretty, ncol = 1) +
  scale_color_manual(values = c(ansys_pal[2], ansys_pal[1]), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[2], ansys_pal[1]), name = "Analysis") +
  theme(legend.position = "top") +
  ylab("Infections at reporting date") +
  xlab("")
  
  

lptrajs_cases_allansys

# cowplot::ggsave2(lptrajs_cases_allansys, file = paste0(snakemake@output, "lptrajs_cases_allansys_ggplot.pdf"),
#                  width = 17.5, height = 7, units = "cm")

lptrajs_cases <- ggplot(trajs_d %>% filter(variable %in% c("D"), deme != "China",
                                     analysis != "mGLM_6p7e_sumprior_spGLMid",
                                     dataset == "eubdmm_open.2") %>%
                    arrange(analysis, dataset, Sample, date) %>% # sum all european cases
                    group_by(date, dataset, analysis, Sample) %>% summarise(cumvalue = sum(cumvalue), .groups = "drop") %>% 
                    group_by(date, dataset, analysis) %>%
                    summarise(mean = mean(cumvalue),
                              median = median(cumvalue),
                              l95 = HDInterval::hdi(cumvalue)[1],
                              h95 = HDInterval::hdi(cumvalue)[2]) %>%
                    ungroup()) +
  geom_line(aes(x = date, y =  median, group = interaction(analysis,dataset), color = analysis)) +
  geom_ribbon(aes(date, ymin =  l95, ymax =  h95, group = interaction(analysis,dataset), fill = analysis), 
              alpha = 0.3, linetype = 2) +
  geom_line(data = cases %>% filter(deme!='China') %>% group_by(date) %>% summarise(total_cases = sum(total_cases)),
            aes(date, total_cases, group = 1), linetype = 2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_x_date(limits = as.Date(c("2019-12-01", "2020-03-07"))) +
  theme_eubdmm() +
  scale_color_manual(values = c(ansys_pal[1], ansys_pal[2]), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[1], ansys_pal[2]), name = "Analysis") 

lptrajs_cases

# Plot 3: First events ---------------------------------------------------------
events <- trajs[variable %in% c("B", "M", "D", "S") & value != 0]
setorder(events, Sample, dataset, analysis, deme, variable, date)

first_events <- events[, .SD[1], 
                    by = .(Sample, dataset, analysis, deme, variable)]

im <- events[variable == "M"]
setorder(im, Sample, dataset, analysis, deme_to, variable, date)

first_im <- im[, .SD[1], by = .(Sample, dataset, analysis, deme_to, variable)]
first_im[, "variable"] = "IM"
setnames(first_im, "deme_to", "deme")
cols_keep <- c("Sample", "dataset", "analysis", "deme", "variable", "date", "value")
first_events <- first_events[, ..cols_keep]
first_im     <- first_im[, ..cols_keep]

first_events <- rbindlist(list(first_events, first_im),  use.names=TRUE)

dptrajs_first_events <- ggplot(first_events %>% filter(variable %in% c("B", "M", "D", "IM")) %>% make_pretty()) +
  geom_density_ridges2(aes(ymd(date), y = variable, 
                           group = interaction(variable, dataset), fill = variable,), 
                       alpha = 0.7, linewidth = 0.3) +
  geom_vline(data = cases %>% group_by(deme) %>% filter(cases != 0) %>% arrange(date) %>% slice(1),
             aes(xintercept = date)) +
  geom_vline(data = first_events %>% group_by(deme) %>% filter(variable == "S") %>% arrange(date) %>% slice(1),
             aes(xintercept = date), linetype = 2) +
  facet_grid(analysis_pretty + dataset_pretty ~deme) +
  theme_eubdmm() +
  scale_fill_manual(values = as.character(ecolors)) +
  scale_color_manual(values = as.character(ecolors)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b")

dptrajs_first_events

cowplot::ggsave2(dptrajs_first_events, file = paste0(snakemake@output, "dptrajs_first_events_ggplot.pdf"),
                 width = 17.5, height = 15, units = "cm")

# Plot 4 Cumulative probability of first immigration ---------------------------

first_m_europe <- trajs_m %>%
  filter(value != 0) %>%
  group_by(Sample, analysis, dataset) %>%
  summarise(first_date = min(date), .groups = "drop")

all_dates <- seq(min(trajs_m$date), max(trajs_m$date), by = "day")

start_posterior_prob <- trajs_m %>%
  select(analysis, dataset, date) %>%
  distinct() %>%
  left_join(first_m_europe, by = c("analysis", "dataset")) %>%
  group_by(analysis, date, dataset) %>%
  summarise(
    p_europe = mean(first_date <= date, na.rm = TRUE),
    .groups = "drop"
  )

start_posterior_prob %>% group_by(analysis, dataset) %>% filter(p_europe >=5) %>% arrange(date) %>% slice(1)
start_posterior_prob %>% group_by(analysis, dataset) %>% filter(p_europe >=1) %>% arrange(date) %>% slice(1)

lp_cumeustart_all <- ggplot(start_posterior_prob  %>% make_pretty(), 
                            aes(x = date, y = p_europe, color = analysis_pretty)) +
  geom_line(size = 1) +
  labs(
    x = "",
    y = "Posterior probability Europe outbreak started"
  ) +
  theme_eubdmm() +
  facet_wrap(~dataset_pretty, ncol = 1) +
  geom_vline(xintercept = ymd("2020-01-25")) +
  scale_color_manual(values = c(ansys_pal[2], ansys_pal[1], "#DC9000"), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[2], ansys_pal[1], "#DC9000"), name = "Analysis") 

lp_cumeustart_all

#cowplot::ggsave2(lp_cumeustart_all, file = paste0(snakemake@output, "lp_cumeustart_all_ggplot.pdf"))
#width = 17.5, height = 15, units = "cm")

lp_cumeustart <- ggplot(start_posterior_prob  %>% 
                          filter(analysis != "mGLM_6p7e_sumprior_spGLMid",
                                 dataset == "eubdmm_open.2") %>% make_pretty(), 
                            aes(x = date, y = p_europe, color = analysis_pretty)) +
  geom_line(size = 1) +
  labs(
    x = "Date",
    y = "Posterior probability Europe outbreak started"
  ) +
  theme_minimal() +
  geom_vline(xintercept = ymd("2020-01-25")) +
  scale_color_manual(values = c(ansys_pal[1], ansys_pal[2], "#DC9000"), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[1], ansys_pal[2], "#DC9000"), name = "Analysis") 

lp_cumeustart

# Plot 5: Probability of local transmission larger than imported cases ---------

trajs_b <- bingrid_trajectories_dt(trajs[variable == "B" & deme != "China"], bin = "day")
setnames(trajs_b, "value", "births")

trajs_im <- bingrid_trajectories_dt(trajs[variable == "M" & deme_to != "China"], bin = "day")
setnames(trajs_im, c("deme", "deme_to", "value"), c("deme_origin", "deme", "mig"))
trajs_im_sum <- trajs_im[
  , .(
    mig = sum(mig)
  ),
  by = .(Sample, analysis, dataset, date, deme)
]

localmig_dt <- merge(
  trajs_b[, .(Sample, deme, analysis, dataset, date, births)],
  trajs_im_sum[, .(Sample, deme, analysis, dataset, date, mig)],
  by = c("Sample", "deme", "analysis", "dataset", "date"),
  all = TRUE
)
localmig_dt[is.na(births), births := 0]
localmig_dt[is.na(mig), mig := 0]

localmig_dt[, more_local := births > mig]

# Posterior probability
local_posterior_prob <- localmig_dt[
  , .(p_morelocal = mean(more_local)),
  by = .(analysis, dataset, date, deme)
]


lp_localprob_deme <- ggplot(local_posterior_prob  %>% filter(deme != "China") %>% make_pretty(), 
                            aes(x = date, y = p_morelocal, color = analysis_pretty)) +
  geom_line(size = 1) +
  labs(
    x = "Date",
    y = "P(local transmission > importations)"
  ) +
  theme_minimal() +
  geom_vline(data = cases %>% group_by(deme) %>% filter(cases != 0, deme != "China") %>% 
               arrange(date) %>% slice(1),
             aes(xintercept = date), linetype = 2) +
  #geom_vline(xintercept = ymd("2020-01-25")) +
  facet_wrap(dataset_pretty~deme, ncol = 4) +
  scale_color_manual(values = c(ansys_pal[1], ansys_pal[2], "#DC9000"), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[1], ansys_pal[2], "#DC9000"), name = "Analysis") 

lp_localprob_deme


# cowplot::ggsave2(lp_localprob_all, file = paste0(snakemake@output, "lp_localprob_all_ggplot.pdf"))
#width = 17.5, height = 15, units = "cm")

# Plot All Europe
localmig_dt_eu <- localmig_dt[
  , .(
    births = sum(births),
    mig = sum(mig)
  ),
  by = .(Sample, analysis, dataset, date)
]

localmig_dt_eu[, more_local := births > mig]
local_posterior_prob_eu <- localmig_dt_eu[
  , .(p_morelocal = mean(more_local)),
  by = .(analysis, dataset, date)
]


local_posterior_prob_eu %>% group_by(analysis, dataset) %>% filter(p_morelocal >= 0.5) %>% arrange(date) %>% slice(1)
local_posterior_prob_eu %>% group_by(analysis, dataset) %>% filter(p_morelocal >= 1) %>% arrange(date) %>% slice(1)


lp_localprob_eu <- ggplot(local_posterior_prob_eu  %>% 
                         filter(analysis != "mGLM_6p7e_sumprior_spGLMid",
                                dataset == "eubdmm_open.2") %>% make_pretty(), 
                           aes(x = date, y = p_morelocal, color = analysis_pretty)) +
  geom_line(size = 1) +
  labs(x = "Date",
       y = "P(local transmission > importations)") +
  theme_minimal() +
  geom_vline(data = cases %>% ungroup %>% filter(cases != 0, deme != "China") %>% 
               arrange(date) %>% slice(1),
             aes(xintercept = date), linetype = 2) +
  #geom_vline(xintercept = ymd("2020-01-25")) +
  scale_color_manual(values = c(ansys_pal[2], ansys_pal[1], "#DC9000"), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[2], ansys_pal[1], "#DC9000"), name = "Analysis") 

lp_localprob_eu

lp_localprob_eu_all <- ggplot(local_posterior_prob_eu  %>% 
                            filter(analysis != "mGLM_6p7e_sumprior_spGLMid") %>% make_pretty(), 
                          aes(x = date, y = p_morelocal, color = analysis_pretty)) +
  geom_line(size = 1) +
  labs(x = "",
       y = "P(local transmission > importations)") +
  theme_eubdmm() +
  geom_vline(data = cases %>% ungroup %>% filter(cases != 0, deme != "China") %>% 
               arrange(date) %>% slice(1),
             aes(xintercept = date), linetype = 2) +
  #geom_vline(xintercept = ymd("2020-01-25")) +
  scale_color_manual(values = c(ansys_pal[2], ansys_pal[1], "#DC9000"), name = "Analysis") +
  scale_fill_manual(values = c(ansys_pal[2], ansys_pal[1], "#DC9000"), name = "Analysis") +
  facet_wrap(~dataset_pretty, ncol = 1)
lp_localprob_eu_all


trajs_eu_all <- lptrajs_cases_allansys + lp_cumeustart_all + lp_localprob_eu_all +
  plot_layout(guides = "collect", widths = c(1.2, 1, 1)) &
  theme(legend.position = 'top')
trajs_eu_all

cowplot::ggsave2(trajs_eu_all, file = paste0(snakemake@output, "trajs_eu_all_ggplot.pdf"),
  width = 17.5, height = 12, units = "cm")

# Plot 6: Exportation per deme -------------------------------------------------
# Filter out last dat because creates a week value with only one day
trajs_om <- bingrid_trajectories_dt(trajs[variable == "M" & date != "2020-03-08"], bin = "week")

deme_prob <- trajs_om[
  , {
    value_h  <- as.numeric(HDInterval::hdi(value))
    
    .(
      value_median  = as.numeric(median(value)),
      value_mean    = as.numeric(mean(value)),
      value_l95     = value_h[1],
      value_h95     = value_h[2]
    )
  },
  by = .(deme, deme_to, analysis, dataset, date)
]

trajs_om_sum <- trajs_om[
  , .(
    value = sum(value)
  ),
  by = .(Sample, analysis, dataset, date, deme)
]

trajs_om_sum[, p_deme := value / sum(value),
            by = .(Sample, analysis, dataset, date)]
trajs_om_sum[is.na(p_deme), p_deme := 0]
deme_prob_summary <- trajs_om_sum[
  , {
    p_deme_h <- as.numeric(HDInterval::hdi(p_deme))
    value_h  <- as.numeric(HDInterval::hdi(value))
    
    .(
      p_deme_median = as.numeric(median(p_deme)),
      p_deme_mean   = as.numeric(mean(p_deme)),
      p_deme_l95    = p_deme_h[1],
      p_deme_h95    = p_deme_h[2],
      
      value_median  = as.numeric(median(value)),
      value_mean    = as.numeric(mean(value)),
      value_l95     = value_h[1],
      value_h95     = value_h[2]
    )
  },
  by = .(deme, analysis, dataset, date)
]

lp_trajs_probmig <- ggplot(deme_prob_summary %>% make_pretty()) +
  geom_line(aes(date, p_deme_median, color = deme)) +
  geom_ribbon(aes(date, ymin = p_deme_l95, ymax = p_deme_h95, fill = deme), alpha= 0.2) +
  facet_grid(analysis_pretty ~ dataset_pretty) +
  scale_fill_manual(values = pal_eubdmm2) +
  scale_color_manual(values = pal_eubdmm2) +
  theme_eubdmm()

lp_trajs_probmig

#cowplot::ggsave2(lp_trajs_probmig, file = paste0(snakemake@output, "lp_trajs_probmig_ggplot.pdf"))

lp_trajs_valuemig <- ggplot(deme_prob  %>% make_pretty()) +
  geom_line(aes(date, value_median, color = deme)) +
  geom_ribbon(aes(date, ymin = value_l95, ymax = value_h95, fill = deme), alpha= 0.2) +
  facet_grid(analysis_pretty + deme_to ~ dataset_pretty, scales = "free_y") +
  scale_fill_manual(values = pal_eubdmm2, name = "") +
  scale_color_manual(values = pal_eubdmm2, name = "") +
  xlab("") +
  ylab("Number of imported infections") +
  theme_eubdmm() +
  scale_y_log10() +
  #coord_cartesian(ylim = c(0, 50)) +
  theme(legend.position = "top")

lp_trajs_valuemig

cowplot::ggsave2(lp_trajs_valuemig, file = paste0(snakemake@output, "lp_trajs_valuemig_ggplot.pdf"),
                 width = 17.5, height = 17, units = "cm")

lp_trajs_valuemig <- ggplot(deme_prob_summary %>% make_pretty()) +
  geom_line(aes(date, value_median, color = deme)) +
  geom_ribbon(aes(date, ymin = value_l95, ymax = value_h95, fill = deme), alpha= 0.2) +
  facet_grid(analysis_pretty ~ dataset_pretty) +
  scale_fill_manual(values = pal_eubdmm2, name = "") +
  scale_color_manual(values = pal_eubdmm2, name = "") +
  xlab("") +
  ylab("Number of exported infections") +
  theme_eubdmm() +
  coord_cartesian(ylim = c(0, 50)) +
  theme(legend.position = "top")

lp_trajs_valuemig
cowplot::ggsave2(lp_trajs_valuemig, file = paste0(snakemake@output, "lp_trajs_valuemig_ggplot.pdf"),
                 width = 17.5, height = 12, units = "cm")


sb_trajs_mig_all <- ggplot(deme_prob_summary %>% make_pretty, 
       aes(x = date,node = deme, fill = deme, value = value_median)) +
  geom_sankey_bump(space = 0, type = "alluvial", color = "transparent", smooth = 8, alpha = 0.9) +
  scale_fill_manual(values = pal_eubdmm2, name = "Origin country") +
  facet_grid(analysis_pretty ~ dataset_pretty) +
  theme_eubdmm() +
  theme(legend.position = "top") +
  xlab("") +
  ylab("Exported infections") +
  coord_cartesian(ylim = c(0,1500))
sb_trajs_mig_all
cowplot::ggsave2(sb_trajs_mig_all, file = paste0(snakemake@output, "sb_trajs_mig_all_ggplot.pdf"),
                 width = 17.5, height = 12, units = "cm")




sb_trajs_mig_s2 <- ggplot(deme_prob_summary %>% make_pretty %>% 
                            filter(analysis != "mGLM_6p7e_sumprior_spGLMid",
                                   dataset == "eubdmm_open.2"), 
                       aes(x = date,node = deme, fill = deme, value = value_median)) +
  geom_sankey_bump(space = 0, type = "alluvial", color = "transparent", smooth = 8, alpha = 0.9) +
  scale_fill_manual(values = pal_eubdmm2) +
  facet_wrap(~analysis_pretty, ncol = 2) +
  theme_eubdmm()

sb_trajs_mig_s2


fig_epi_main <- (
  (
    lptrajs_cases + theme(legend.position = "none") +
      (
        (lp_cumeustart + theme(legend.position = "top")) /
          (lp_localprob_eu + theme(legend.position = "none"))
      )
  ) + 
    plot_layout(widths = c(2, 1))
) /
  (
    sb_trajs_mig_s2
  ) +
  plot_layout(heights = c(2, 1.5))

fig_epi_main
cowplot::ggsave2(fig_epi_main, file = paste0(snakemake@output, "fig_epi_main_ggplot.pdf"))


sankey_df <- summarise_bintrajectories_dt(trajs_om) %>%
  group_by(analysis, dataset) %>%
  mutate(bin = date, #floor_date(date, "month"),
         next_bin = bin  %m+% weeks(1)) %>%
  uncount(weights = round(value_med))  %>%
         #bin = if_else(bin == ymd("2019-12-01"), ymd("2019-12-08"), bin),
         #next_bin = if_else(next_bin == ymd("2020-04-01"), ymd("2020-03-08"), next_bin)) %>%
  select(deme, deme_to, bin, next_bin, analysis, dataset) 

sk_trajsmig <- ggplot(sankey_df %>% make_pretty(),
       aes(x = bin, 
           next_x = next_bin, 
           node = deme, 
           next_node = deme_to,
           fill = factor(deme),
           #color = factor(deme),
           label = deme)) +
  geom_alluvial(flow.alpha = .6, width = 1) +
  scale_fill_manual(values = pal_eubdmm) +
  theme_alluvial(base_size = 18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))  +
  facet_grid(analysis_pretty ~ dataset_pretty, scales = "free")

sk_trajsmig
cowplot::ggsave2(sk_trajsmig, file = paste0(snakemake@output, "sk_trajsmig_ggplot.pdf"))
