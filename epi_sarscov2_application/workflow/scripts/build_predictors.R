
# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project:  epiGLM
#        V\ Y /V    Code for GLM predictors Migration Matrix Epi GLM
#    (\   / - \     February 2025
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(eurostat)
library(countrycode)
library(wbstats)
library(sf)
library(geosphere)
library(maps)
source("workflow/scripts/flight_data.R")
source("workflow/scripts/plot_opts.R")


# ------------------------------------------------------------------------------
# To run outside of the snakemake workflow, run this code first 
# ------------------------------------------------------------------------------
# 
# setClass(
#   "snakemake_object",
#    contains= "tbl_df",
#    slots = c(input = "character", output = "character", params = "character")
#  )
# 
# snakemake <- new("snakemake_object", tibble(),
#                  params = c(osiindex = "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/refs/heads/master/data/OxCGRT_nat_latest.csv"),
#                  input = c(weights_file = "results/data/otherEU_weights.tsv"),
#                  output = c(out_folder = "results/predictors/"))
# ------------------------------------------------------------------------------


dir.create(snakemake@output[["out_folder"]], recursive = TRUE, showWarnings = FALSE)

get_predictor_matrix <- function(df_long, deme_list, logtransform, scale, change_times, out_file = NA) {
  df <- df_long %>% 
    filter(deme.x != deme.y) %>% rename(date = time) %>%
    mutate(values = value) 
  
  if (logtransform) {
    df <- df %>%
      mutate(values = log(values + 1))
  }
  if (scale) {
    df <- df %>% ungroup %>%
      mutate(values = scale(values)[,1])
  }
  # Add all time changes in the analysis
  base <- change_times %>% cross_join(tibble(deme.x = names(deme_list)))  %>% 
    cross_join(tibble(deme.y = names(deme_list))) %>%
    filter(deme.x != deme.y) %>% mutate(date = ymd(date))
  
  df_full <- df %>% right_join(base) %>%
    group_by(deme.x, deme.y) %>%
    arrange(time, 
            factor(deme.x, levels = names(deme_list)), 
            factor(deme.y, levels = names(deme_list))) %>%
    fill(values, .direction = "up") %>%
    fill(values, .direction = "down")  %>%
    ungroup()
  
  predictor_m <- as.data.frame(matrix(df_full %>% 
                                        mutate_if(is.numeric, ~round(., 10)) %>% 
                                        pull(values), 
                                      ncol = length(demes) - 1, byrow = TRUE)) 
  if (!is.na(out_file)) write_tsv(predictor_m, out_file, col_names = F, na = "")
  
  return(df_full)
}

to_num <- function(x, mrs) {
  decimal_date(ymd(mrs)) - decimal_date(ymd(x))
}


# 0. Analysis specs ------------------------------------------------------------

demes <- c(China = "CN", France = "FR", Germany = "DE", Italy = "IT", OtherEU = "OE")
analysis_from <- "2019-12-01"
analysis_to <- "2020-03-08"
# data_from <- "2019-12-26"
# data_to <- "2020-03-08"
weights <- read_tsv(snakemake@input[["weights_file"]]) 

# Which countries are included in the OtherEU deme
eu_countries <- eurostat::get_eurostat("avia_paocc", select_time = "M") %>%
  filter(!grepl("EU", geo)) %>%
  select(geo) %>% unique() %>%
  mutate(country = countrycode(geo, "eurostat", "country.name"))


# 1. Daily flight passengers EUROSTATS -----------------------------------------
# Time and type-type-dependent info

# flight_df <- get_flightData(countries = demes, from = analysis_from, to = analysis_to) %>%
#   rowwise %>% mutate(deme.x = names(demes[demes == geo]),
#                      deme.y = names(demes[demes == partner]),
#                      value = daily_value) %>% ungroup
# 

flight_df_all <- get_flightData(countries = c(eu_countries$geo, "CN"), from = analysis_from, to = analysis_to) %>%
  rowwise %>% mutate(deme.x = ifelse(geo %in% demes, names(demes[demes == geo]), "OtherEU"),
                     deme.y = ifelse(partner %in% demes, names(demes[demes == partner]), "OtherEU")) %>% ungroup

flight_df <- flight_df_all %>%
  left_join(weights) %>%
  mutate(w = case_when(deme.x != "OtherEU" ~ 1,
                       !is.na(w) ~ w,
                       deme.x == "OtherEU" ~ 0)) %>%
  group_by(deme.x, deme.y, time) %>%
  summarise(value = sum(w * daily_value), .groups = "drop") %>%
  filter(deme.x != deme.y)


# 2a. Population origin --------------------------------------------------------
# Type-dependent info
# EU Countries considered for flight data

pop_data <- wb_data("SP.POP.TOTL", start_date = 2020, end_date = 2020) %>%
  mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country),
         iso2c = ifelse(country == "Greece", "EL", iso2c)) %>%
  filter(country %in% c(eu_countries$country, "China")) %>% 
  mutate(deme.x = case_when(iso2c %in% demes ~ country,
                            TRUE ~ "OtherEU")) %>%
  left_join(weights) %>%
  mutate(w = case_when(!is.na(w) ~ w,
                       deme.x == "OtherEU" ~ 0,
                       TRUE ~ 1)) %>%
  group_by(deme.x) %>%
  summarise(value = sum(w * SP.POP.TOTL))

pop_data_all <- wb_data("SP.POP.TOTL", start_date = 2020, end_date = 2020) %>% 
  mutate(country = ifelse(country == "Slovak Republic", "Slovakia", country),
         iso2c = ifelse(country == "Greece", "EL", iso2c),
         iso2c = ifelse(country == "United Kingdom", "UK", iso2c)) %>%
  filter(country %in% c(eu_countries$country, "China")) %>% 
  mutate(deme.x = case_when(iso2c %in% demes ~ country,
                            TRUE ~ "OtherEU"),
         value = SP.POP.TOTL)

pop_x_df <- pop_data %>% 
  cross_join(pop_data %>% select(deme.y = deme.x)) %>%
  mutate(time = ymd(analysis_from)) %>%
  filter(deme.x != deme.y)

pop_all_x_df <- pop_data_all %>% 
  cross_join(pop_data %>% select(deme.y = deme.x)) %>%
  mutate(time = ymd(analysis_from))


# 2b. Population destination ---------------------------------------------------
# Type-dependent info
pop_y_df <- pop_x_df %>% 
  mutate(temp = deme.x,
         deme.x = deme.y,
         deme.y = temp)

# 3. Daily flight passengers EUROSTATS/inhabitant origin --------------
# Time and type-type-dependent info

flight_pop_x_df <- flight_df_all %>% 
  left_join(pop_data_all %>% select(iso2c, deme.x, country, pop = value), 
            by = c("geo" = "iso2c", "deme.x")) %>%
  left_join(weights) %>%
  mutate(w = case_when(deme.x != "OtherEU" ~ 1,
                       !is.na(w) ~ w,
                       deme.x == "OtherEU" ~ 0)) %>%
  mutate(value_flightpop = daily_value * 100000 / pop) %>%
  group_by(deme.x, deme.y, time) %>%
  summarise(
    value = sum(w * value_flightpop, na.rm = T)
  ) %>%
  filter(deme.x != deme.y) %>% ungroup

# 4. Distance between countries ------------------------------------------------
# Type-type-dependent info

sf::sf_use_s2(FALSE)
map <- st_as_sf(map("world", plot = FALSE, fill = TRUE)) %>%
  mutate(country = as.character(ID),
         continent = countrycode(sourcevar = country,
                                 origin = "country.name",
                                 destination = "continent")) %>%
  filter(country %in% c(eu_countries$country, "China")) %>%
  mutate(deme = ifelse(country %in% names(demes), country, "OtherEU"))

map <- cbind(map, st_coordinates(st_centroid(map))) %>%
  rename(longitude = X, latitude = Y) #%>%

# Calculate country centroids
centroids <- data.frame(map) %>%
  select(deme, country, longitude, latitude)

centroids_join <- centroids %>%
  cross_join(centroids)

# Compute distance between centroids
dist_df <- centroids_join %>%
  mutate(dist = distHaversine(select(centroids_join, longitude.x, latitude.x),
                              select(centroids_join, longitude.y, latitude.y))) %>%
  left_join(weights, by = c("country.x" = "country")) %>%
  left_join(weights, by = c("country.y" = "country")) %>%
  mutate(w.x = case_when(deme.x != "OtherEU" ~ 1,
                       !is.na(w.x) ~ w.x,
                       #T ~ 0
                       deme.x == "OtherEU" ~ 0
                       ),
         w.y = case_when(deme.y != "OtherEU" ~ 1,
                        !is.na(w.y) ~ w.y,
                        #T ~ 0
                        deme.y == "OtherEU" ~ 0
         ),
         w = case_when(
           !is.na(w.x) & !is.na(w.y) ~ w.x * w.y,
           !is.na(w.x) & is.na(w.y) ~ w.x,
           is.na(w.x) & !is.na(w.y) ~ w.y,
           T ~ 0) ) %>%
  filter(deme.x != deme.y) %>%
  group_by(deme.x, deme.y) %>%
  summarise(value = sum(w * dist),
            time = ymd(analysis_from), .groups = "drop")

# Taking the weighted average value for OtherEU deme by cases

# 5. Shared borders ------------------------------------------------------------
# Type-type-dependent info
sharedborders_df <- tibble(deme.x = names(demes)) %>% cross_join(tibble(deme.y = names(demes))) %>%
  mutate(value = c(NA, 0, 0, 0, 0,
                   0, NA, 1, 1, 0,
                   0, 1, NA, 0, 0,
                   0, 1, 0, NA, 0,
                   0, 0, 0, 0, NA),
         time = ymd(analysis_from)) %>%
  filter(deme.x != deme.y)
# Assuming no shared border for OtherEU deme

# 6a. OSI/KOF index origin -----------------------------------------------------
# Time and type dependent info

oxfordGRT <- read_csv(snakemake@params[["osiindex"]]) #url
osi_x_df <- oxfordGRT %>%
  filter(is.na(RegionCode)) %>%
  mutate(Date = parse_date(as.character(Date), format = "%Y%m%d")) %>%
  select(c(country = CountryName, date = Date, 
           osi = StringencyIndex_Average)) %>%
  mutate(date_15 = floor_date(date, unit = period(15, "days")),
         time = case_when(day(date_15) == 31 ~ ymd(date_15) + 1,
                          TRUE ~ ymd(date_15)),
         #deme.x = ifelse(country %in% names(demes), country, "OtherEU")) %>% # It si including non europe
          deme.x = case_when(country %in% names(demes) ~ country,
                             country %in% eu_countries$country ~ "OtherEU",
                             TRUE ~ NA)) %>%
  filter(deme.x %in% c(eu_countries$country, names(demes)),
         time < analysis_to) %>%
  left_join(weights) %>%
  mutate(w = case_when(deme.x != "OtherEU" ~ 1,
                       !is.na(w) ~ w,
                       deme.x == "OtherEU" ~ 0)) %>%
  group_by(deme.x, time) %>%
  summarise(value = sum(w * osi), .groups = "drop") %>%
  #filter(time != "2020-03-31") %>%
  cross_join(tibble(deme.y = names(demes)))

# 6b. OSI/KOF index destination ------------------------------------------------
# Time and type dependent info
osi_y_df <- osi_x_df %>%
  mutate(temp = deme.x,
         deme.x = deme.y,
         deme.y = temp)

# 7. Border closure ------------------------------------------------------------
# Time and type-type dependent info
closing_dates <- tibble(country = names(demes), 
                        date = c("2020-01-23", "2020-03-16", "2020-03-16",  "2020-03-14",  "2020-03-16"))

closedborders_df <- tibble(deme.x = names(demes)) %>% cross_join(tibble(deme.y = names(demes))) %>%
  cross_join(closing_dates %>% select(time = date) %>% unique %>% bind_rows(tibble(time = c("2020-01-01"))) ) %>%
  rowwise %>%
  mutate(value = case_when(
    (time < filter(closing_dates, country == deme.x)[["date"]]) & (time < filter(closing_dates, country == deme.y)[["date"]]) ~ 0,
    TRUE ~ 1),
    time = ymd(time)) %>%
  arrange(time, deme.x, deme.y) %>%
  filter(time <= analysis_to)


# 8. Change dates --------------------------------------------------------------
get_changetimes <- function(dfs, mrs, out_file = NA) {
  
  change_dates <- bind_rows(dfs) %>% select(time) %>% unique %>% 
    mutate(date = ymd(time),
           time = to_num(ymd(date), mrs)) %>%
    unique %>% arrange(desc(date)) 
  
  change_times <- change_dates %>% select(time) %>%
    slice(1:(n() - 1))
  
  if (!is.na(out_file)) write_tsv(change_times %>% select(time), out_file, col_names = F, na = "")
  
  return(change_dates)
}

change_dates_flights <- get_changetimes(flight_df, analysis_to, 
                                        out_file = paste0(snakemake@output[["out_folder"]], 
                                                          "/changetimes_flights_4e.tsv"))
change_dates_all <-  get_changetimes(list(flight_df, osi_x_df, closedborders_df), analysis_to, 
                                     out_file = paste0(snakemake@output[["out_folder"]], 
                                                       "/changetimes_all_7e.tsv"))


# 9. Transform predictors to matrix shape and transform and save ---------------
change_dates_list <- list(change_dates_all, change_dates_flights)
pred_list <- tibble::lst(flight_df, flight_pop_x_df, pop_x_df, pop_y_df, dist_df, 
                  sharedborders_df, osi_x_df, osi_y_df, closedborders_df)

df_preds <-lapply(change_dates_list, function(change_dates) {
  
  lapply(1:length(pred_list), function(i) {
    
    lapply(c(0,1), function(logscale) {
      pred <- pred_list[[i]]
      
      file_name <- paste0(str_replace(names(pred_list[i]), "df", ""),
                          nrow(change_dates), "e",
                          ifelse(logscale, "_ls.tsv", "_nt.tsv"))
      print(file_name)
      get_predictor_matrix(pred, deme_list = demes, logtransform = logscale, scale = logscale,  
                           change_times = change_dates,
                           #out_file = paste0(snakemake@output[["out_folder"]], "/", file_name)
                           )
    })
  })
})


# Visualize predictors
# 
lp_flights <- ggplot(df_preds[[1]][[1]][[1]]) +
  geom_step(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y), direction="vh", linewidth = 0.7) +
  #geom_point(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y)) +
  facet_grid(~deme.x) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2, name = "") +
  xlab("")

lp_nflights <- ggplot(df_preds[[1]][[2]][[1]]) +
  geom_step(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y), direction="vh", linewidth = 0.7) +
  #geom_point(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y)) +
  facet_grid(~deme.x) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2, name = "") +
  xlab("")

lp_dist <- ggplot(df_preds[[1]][[5]][[1]]) +
  geom_step(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y), direction="vh", linewidth = 0.7) +
  #geom_point(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y)) +
  facet_grid(~deme.x) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2, name = "") +
  xlab("")

lp_osidest <- ggplot(df_preds[[1]][[7]][[1]]) +
  geom_step(aes(date, values, group = interaction(deme.x, deme.y), color = deme.x), direction="vh", linewidth = 0.7) +
  #geom_point(aes(date, values, group = interaction(deme.x, deme.y), color = deme.y)) +
  facet_grid(~deme.x) +
  theme_eubdmm() +
  scale_color_manual(values = pal_eubdmm2, name = "") +
  xlab("")

plot_preds <- (lp_flights / lp_nflights / lp_dist / lp_osidest) +
  plot_layout(guides = "collect") & theme(legend.position = "top")

ggsave(plot_preds, filename = "results/figures/plot_preds.pdf",
       width = 17.4, height = 17, units = "cm")


# Check pred corr
corr_df = left_join(df_preds[[1]][[1]][[1]] %>% select(deme.x, deme.y, time, flights=values),
          left_join(df_preds[[1]][[2]][[1]] %>% select(deme.x, deme.y, time, flights_pop_x=values),
          #left_join(df_preds[[1]][[3]][[1]] %>% select(deme.x, deme.y, time, pop_x=values),
          #left_join(df_preds[[1]][[4]][[1]] %>% select(deme.x, deme.y, time, pop_y=values),
          left_join(df_preds[[1]][[5]][[1]] %>% select(deme.x, deme.y, time, dist=values),
          left_join(df_preds[[1]][[6]][[1]] %>% select(deme.x, deme.y, time, shared_border=values),
          #left_join(df_preds[[1]][[7]][[1]] %>% select(deme.x, deme.y, time, osi_x=values),
          left_join(df_preds[[1]][[8]][[1]] %>% select(deme.x, deme.y, time, osi_y=values),
                    df_preds[[1]][[9]][[1]] %>% select(deme.x, deme.y, time, closed_border=values)
                    ))))) %>%
  select(flights,flights_pop_x,dist,shared_border,osi_y,closed_border)

pdf("results/figures/corplot_preds.pdf", width = 5, height = 5)

corrplot::corrplot(
  cor(corr_df),
  type = "upper",
  addCoef.col = "white",
  number.cex = 0.6,
  tl.col = "black"
)

dev.off()
