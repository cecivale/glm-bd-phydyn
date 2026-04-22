# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  epiGLM
#        V\ Y /V    Scrip to calculate OtherEu weights
#    (\   / - \     September 2024
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)

# ------------------------------------------------------------------------------
# To run outside of the snakemake workflow, run this code first 
# ------------------------------------------------------------------------------

# setClass(
#   "snakemake_object",
#    contains= "tbl_df",
#    slots = c(input = "list", output = "character", params = "character")
#  )

# snakemake <- new("snakemake_object", tibble(),
#                  input = list(cases = "resources/ext_ecdccases.tsv"),
#                  output = list(weights = "results/data/otherEU_weights.tsv"))
# ------------------------------------------------------------------------------



cases <- read_tsv(snakemake@input[["cases"]]) %>% mutate(date = dmy(dateRep))

# Which countries are included in the OtherEU deme
eu_countries <- eurostat::get_eurostat("avia_paocc", select_time = "M") %>%
  filter(!grepl("EU", geo)) %>%
  select(geo) %>% unique() %>%
  mutate(country = countrycode(geo, "eurostat", "country.name"))

deaths_otherEU <- cases %>% 
  filter(!countriesAndTerritories %in% c("China", "France", "Italy", "Germany")) %>%
  filter(geoId %in% eu_countries$geo) %>%
  filter(date <= "2020-03-08") %>%
  group_by(countriesAndTerritories) %>%
  summarise(totaldeaths = sum(deaths)) %>% arrange(desc(totaldeaths))

w_cases_otherEU <- cases %>% 
  filter(!countriesAndTerritories %in% c("China", "France", "Italy", "Germany")) %>%
  filter(geoId %in% eu_countries$geo) %>%
  filter(date <= "2020-03-08") %>%
  group_by(country = countriesAndTerritories, geo = geoId) %>%
  summarise(totalcases = sum(cases), .groups = "drop") %>% arrange(desc(totalcases)) %>%
  mutate(w = totalcases / sum(totalcases),
         country = str_replace_all(country, "_", " ")) %>%
  select(country, geo, w)

write_tsv(w_cases_otherEU, snakemake@output[["weights"]])



