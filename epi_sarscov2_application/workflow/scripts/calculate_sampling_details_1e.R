# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  epiGLM
#        V\ Y /V    Scrip to compute subsampling specs
#    (\   / - \     September 2024
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)
library(zoo)


# ------------------------------------------------------------------------------
# To run outside of the snakemake workflow, run this code first 
# ------------------------------------------------------------------------------

# setClass(
#   "snakemake_object",
#    contains= "tbl_df",
#    slots = c(input = "list", output = "character", params = "character")
#  )

# sufix = ".2"
# snakemake <- new("snakemake_object", tibble(),
#                  input = list(cases = "resources/ext_ecdccases.tsv",
#                  ids = paste0("results/data/eubdmm_open/ids_combined", sufix, ".tsv"),
#                  metadata_otherEU = "results/data/eubdmm_open/OtherEU/metadata.tsv",
#                  weights = "results/data/otherEU_weights.tsv"),
#                  output = list(sampling_details =  paste0("results/data/eubdmm_open/sampling_details", sufix, ".tsv")))
# ------------------------------------------------------------------------------


# Read files
cases <- read_tsv(snakemake@input[["cases"]])

demes <- c(China = "CN", France = "FR", Germany = "DE", Italy = "IT", OtherEU = "OE")
e1_from <- "2019-12-01"
e1_to <- "2020-03-08"

cases_ansys <- cases %>%
  mutate(country = countriesAndTerritories) %>%
  mutate(date = ymd(paste0(year, "-", month, "-", day))) %>%
  filter(date <= e1_to, date >= e1_from,
         country %in% names(demes) | continentExp == "Europe") %>%
  left_join(tibble(deme = names(demes), country = names(demes))) %>%
  replace_na(list(deme = "OtherEU")) 


# For France, we  set to 0 the sampling prop before first sequence

cases_total_epoch <- cases_ansys %>%
  mutate(epoch = case_when(
    deme == "France" & date > "2020-03-02" ~ "e2",
    TRUE ~ "e1"
  )) %>%
  group_by(deme, country, epoch) %>%
  summarise(cases = sum(cases))


sequences_ansys <- read_tsv(snakemake@input[["ids"]]) 
metadata_otherEU <- read_tsv(snakemake@input[["metadata_otherEU"]]) 
weights <- read_tsv(snakemake@input[["weights"]]) 

sequences_total_epoch <- sequences_ansys %>%
  left_join(metadata_otherEU %>% select(strain, country)) %>%
  mutate(country = ifelse(is.na(country) & !is.na(deme), deme, country),
         epoch = case_when(
        deme == "France" & date > "2020-03-02" ~ "e2",
        TRUE ~ "e1"
  )) %>%
  count(deme, country, epoch) %>%
  replace_na(list(n = 0))

# Compute sampling probabilities
max_UR = 10

sampling_props <- left_join(cases_total_epoch, sequences_total_epoch) %>%
  left_join(weights) %>%
  replace_na(list(w = 1, n = 0)) %>%
  mutate(sp_upper = n / cases,
         #sp_double_upper = sp_upper * 2,
         sp_lower = n / (cases * max_UR)) %>%
  ungroup %>%
  group_by(deme, epoch) %>%
  summarise(sp_upper = sum(w * sp_upper),
            sp_lower = sum(w * sp_lower)) %>%
  pivot_longer(c(sp_upper, sp_lower), values_to = "value", names_to = "param") %>%
  mutate(param = paste0(param, "_", deme, "_", epoch)) %>%
  ungroup() %>%
  select(param, value)

write_tsv(sampling_props, snakemake@output[["sampling_details"]])


