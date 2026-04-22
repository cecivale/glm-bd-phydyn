# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: Macro Dinosaurs GLM Phylodynamics
#        V\ Y /V    Get ESS/hour
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)
library(coda)
library(patchwork)

source("R/utils.R")
source("R/plot_opts.R")

# Add directories with the chains from BEAST2 results
macro_filedir = ""
epi_filedir = ""

macro_ansys_names <- c("no_GLM", "GLM_sumprior_3p", 
                 "GLM_sumprior_12p", "GLM_sumprior_errordistr_3p")
macro_ansys_names_pretty <- c("non-GLM", "GLM", 
                       "GLM 12p", "stochastic GLM")

epi_ansys_names = c("mGLM_6p7e_sumprior_lowerboundFR0","no_GLM")
epi_ansys_names_pretty = c("GLM","non-GLM")

macro_trace_info <- tibble(application = "Macro",
                          ansys_name = unlist(lapply(macro_ansys_names, function(n) rep(n, 5))),
                          ansys_name_pretty = unlist(lapply(macro_ansys_names_pretty, function(n) rep(n, 5))),
                          data_seed = 1,
                          chain_seed = rep(c("A", "B", "C", "D", "E"), 4),
                          file_dir = macro_filedir,
                          trace_file = paste0(file_dir, ansys_name, "/Dinosaurs_", ansys_name, "_", chain_seed, ".log")) 

epi_trace_info <- tibble(application = "Epi",
                         ansys_name = unlist(lapply(epi_ansys_names, function(n) rep(n, 9))),
                         ansys_name_pretty = unlist(lapply(epi_ansys_names_pretty, function(n) rep(n, 9))),
                         data_seed = rep(unlist(lapply(0:2, function(n) rep(n, 3))),2),
                         chain_seed = as.character(rep(1:3, 6)),
                         file_dir =  epi_filedir,
                         trace_file = paste0(file_dir, ansys_name, "/chains/eubdmm_open.", data_seed, ".", chain_seed, ".log"))

trace_info <- bind_rows(macro_trace_info, epi_trace_info)
                     
                     
# Read times                     
timeresults_macro <- tibble(file = list.files(path = "performance_results/macro/", pattern = "\\.log$", full.names = TRUE)) %>%
  mutate(
    filename = basename(file),
    stem = str_remove(filename, "\\.txt$"),
    stem = str_replace(stem, "3p_errordistr", "3perrordistr"),
    stem = str_replace(stem, "noGLM", "no_GLM"),
    parts = str_split(stem, "_"),
    total_calc_time = map_dbl(file, \(f) {
      read_lines(f) %>%
        tail(10) %>% 
        keep(~ str_detect(.x, "^Total calculation time:")) %>%
        last() %>%
        str_extract("[0-9.]+(?= seconds)") %>%
        as.numeric()
    }),
  ) %>%
  mutate(application = "Macro",
         ansys_name_pretty = case_when(
           grepl("noGLM", filename) ~ "non-GLM",
           grepl("GLM_3p_errordistr", filename) ~ "stochastic GLM",
           grepl("GLM_3p", filename) ~ "GLM",
           grepl("GLM_12p", filename) ~ "GLM 12p"),
         chain_seed = map_chr(parts, \(x) x[3]),
         data_seed = "1"
         ) %>% 
  mutate(total_calc_time = ifelse(filename == "GLM_12p_E_61396843.log", 5*24*3600, total_calc_time)) %>%
  filter(!is.na(total_calc_time)) %>% # Fill NA because of cancel job macro 12p E over time (5 days) and drop NA from failed jobs
  group_by(application, ansys_name_pretty, chain_seed, data_seed) %>%
  summarise(total_calc_time = sum(total_calc_time)) %>%
  select(total_calc_time, application, ansys_name_pretty, chain_seed, data_seed)

timeresults_epi <- tibble(file = list.files(path = "performance_results/epi/", pattern = "\\.txt$", full.names = TRUE)) %>%
  mutate(
    filename = basename(file),
    stem = str_remove(filename, "\\.txt$"),
    parts = str_split(stem, "_"),
    total_calc_time = map_dbl(file, \(f) {
      read_lines(f) %>%
        tail(10) %>%
        keep(~ str_detect(.x, "^Total calculation time:")) %>%
        last() %>%
        str_extract("[0-9.]+(?= seconds)") %>%
        as.numeric()
    }),
    application = "Epi",
    ansys_name_pretty = if_else(str_detect(filename, "mGLM"), "GLM", "non-GLM"),
    data_seed = str_replace(map_chr(parts, \(x) x[length(x) - 1]), "\\.", ""),
    chain_seed = str_replace(map_chr(parts, \(x) x[length(x)]), ".r0", "")
  ) %>%
  select(total_calc_time, application, ansys_name_pretty, data_seed, chain_seed)


# Read trace files
df_ess <- trace_info$trace_file %>%
  set_names() %>%
  map_dfr(function(trace_file) {
    
   trace_wide <- read_trace(trace_file, burninFrac = 0.1) %>%
      select(-starts_with("X"))
    
   steps <- max(trace_wide$Sample, na.rm = TRUE)
    
    enframe(
      effectiveSize(as.mcmc(trace_wide %>% select(-Sample))),
      name = "parameter",
      value = "ess"
    ) %>%
      mutate(
        steps = steps,
        trace_file = trace_file
      )
  })

# Take average value across chains and datasets for same analysis
performance_df <- bind_rows(timeresults_macro, timeresults_epi) %>%
  left_join(df_ess %>% left_join(trace_info, relationship = "many-to-many") %>% 
              mutate(data_seed = as.character(data_seed)), relationship = "many-to-many") %>%
  group_by(application, ansys_name_pretty, parameter) %>%
  summarise(ess_mean = mean(ess),
            ess_sd = sd(ess),
            steps = mean(steps),
            time_mean = mean(total_calc_time),
            time_sd = sd(total_calc_time)) %>%
  mutate(ess_hour = ess_mean / (time_mean / 3600),
            ess_Mstep = ess_mean / steps * 1000000,
            param = str_split(parameter, "\\.|SM", simplify = T)[, 1],
            param = ifelse(param == "samplingRateGLM_value", "samplingRate", param),
            epoch_age =  as.numeric(str_split(parameter, "\\.", simplify = T)[, 2]),) %>%
  filter(ess_mean != 0,
         !grepl("_endtime", param))

performance_df_glmparams <- performance_df %>% filter(grepl("samplingRateSV", parameter) & application == "Macro" |
                   param == "migrationRate") %>% 
  mutate(type = "GLM parameter")

performance_df_distr <- performance_df %>% filter(param %in% c("posterior", "FBD", "BDMMPrime", "treeLikelihood", "fso")) %>%
  mutate(type = case_when(param %in% c("FBD", "BDMMPrime") ~ "BD likelihood",
                          T ~ param))

performance_df_other <- performance_df %>% filter(!param %in% c(performance_df_glmparams$param, performance_df_distr$param, "prior"),
                                  !grepl("GLM|migrationRateValues_", parameter),
                                  param != "error_sigma") %>%
  mutate(type = ifelse(param %in% c("kappa", "gamma_shape"), "Nucleotide substitution parameters", "Other BD parameters"))


min_ess_all <- bind_rows(performance_df_glmparams, performance_df_distr, performance_df_other) %>%
  group_by(application, ansys_name_pretty) %>%
  summarise(min_ess_hour = min(ess_hour), min_ess_Mstep = min(ess_Mstep)) %>%
  mutate(type = "All")

av_ess_all <- bind_rows(performance_df_glmparams, performance_df_distr, performance_df_other) %>%
  group_by(application, ansys_name_pretty) %>%
  summarise(av_ess_hour = mean(ess_hour), av_ess_Mstep = mean(ess_Mstep)) %>%
  mutate(type = "All")

median_ess_all <- bind_rows(performance_df_glmparams, performance_df_distr, performance_df_other) %>%
  group_by(application, ansys_name_pretty) %>%
  summarise(med_ess_hour = median(ess_hour), med_ess_Mstep = median(ess_Mstep)) %>%
  mutate(type = "All")

bp_av_ess <- ggplot(bind_rows(performance_df_glmparams, performance_df_distr, performance_df_other) %>%
                   group_by(type, application, ansys_name_pretty) %>%
                   summarise(av_ess_hour = mean(ess_hour), av_ess_Mstep = mean(ess_Mstep)) %>%
                   ungroup %>%
                   mutate(ansys_name_pretty = factor(ansys_name_pretty, 
                                                     levels = rev(c("non-GLM", "GLM", "GLM 12p", 
                                                                "stochastic GLM")))) %>%
         bind_rows(av_ess_all)) +
  geom_bar(aes(av_ess_hour, ansys_name_pretty, 
               fill = interaction(ansys_name_pretty, application)), color = "white", 
           stat = "identity", position = "dodge") +
  facet_grid(application~type, scales = "free_x") +
  theme_dino2() +
  scale_fill_manual(values = pal_performance)

bp_median_ess <- ggplot(bind_rows(performance_df_glmparams, performance_df_distr, performance_df_other) %>%
                      group_by(type, application, ansys_name_pretty) %>%
                      summarise(med_ess_hour = median(ess_hour), med_ess_Mstep = median(ess_Mstep)) %>%
                      bind_rows(median_ess_all) %>% ungroup %>%
                      mutate(ansys_name_pretty = factor(ansys_name_pretty, 
                                                          levels = rev(c("non-GLM", "GLM", "GLM 12p", 
                                                                     "stochastic GLM"))))
                      ) +
  geom_bar(aes(med_ess_hour, ansys_name_pretty, 
               fill = interaction(ansys_name_pretty, application)), color = "white", 
           stat = "identity", position = "dodge") +
  facet_grid(application~type, scales = "free_x") +
  theme_dino2() +
  scale_fill_manual(values = pal_performance)

bp_min_ess <- ggplot(bind_rows(performance_df_glmparams, performance_df_distr, performance_df_other) %>%
         group_by(type, application, ansys_name_pretty) %>%
         summarise(min_ess_hour = min(ess_hour), min_ess_Mstep = min(ess_Mstep)) %>%
         #filter(type %in% c("GLM parameter", "Other BD parameters")) %>%
         bind_rows(min_ess_all) %>% ungroup() %>%
        mutate(ansys_name_pretty = factor(ansys_name_pretty, 
                                             levels = rev(c("non-GLM", "GLM", "GLM 12p", 
                                                        "stochastic GLM"))))) +
  geom_bar(aes(min_ess_hour, ansys_name_pretty, fill =  interaction(ansys_name_pretty, application)), color = "white", 
           stat = "identity", position = "dodge") +
  facet_grid(application~type, scales = "free_x") +
  theme_dino2() +
  scale_fill_manual(values = pal_performance) 


ggsave(plot = bp_median_ess / (bp_min_ess + theme(legend.position = "none")), filename = paste0("figures/", "ess_ggplot.pdf"), 
       width = 300, height = 110, units = "mm")

  


