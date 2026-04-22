# ------------------------------------------------------------------------------
#          ---        
#        / o o \    Project: Early SARS-CoV-2 Europe Phylodynamics
#        V\ Y /V    Util functions for analysis workflow
#    (\   / - \     
#     )) /    |     
#     ((/__) ||     Code by Ceci VA 
# ------------------------------------------------------------------------------


library(countrycode)
library(ggpubr)
library(lubridate)
library(maps)
#library(sf)
library(wpp2019)
library(data.table)
library(boa) 

to_date <- function(x, mrs) {
  #date(date_decimal(decimal_date(ymd(mrs)) - x))
  as.Date(ymd(mrs) - x * 365.25)
}


to_num <- function(x, mrs) {
  decimal_date(ymd(mrs)) - decimal_date(ymd(x))
}

read_trace <- function(traceFile, burninFrac){
  df <- read_table2(traceFile, comment = "#")
  
  if (burninFrac > 0) {
    n <- dim(df)[1]
    df <- df[-(1:ceiling(burninFrac * n)), ]
  }
  
  return(df)
}


plot_map <- function(demes) {
  # https://www.r-bloggers.com/2019/04/zooming-in-on-maps-with-sf-and-ggplot2/
  worldmap <- st_as_sf(map("world", plot = FALSE, fill = TRUE)) %>%
    mutate(country = ifelse(ID == "UK", "United Kingdom", as.character(ID)),
           continent = countrycode(sourcevar = country,
                                   origin = "country.name",
                                   destination = "continent")) %>%
    #filter(country %in% demes$country | continent == "Europe") %>%
    mutate(deme = case_when(
      country %in% demes$country ~ country,
      continent == "Europe" ~ "Other European"))
  
  zoom_to <- c(20, 50) 
  zoom_level <- 2.5
  target_crs <- sprintf('+proj=aeqd +lon_0=%f +lat_0=%f',
                        zoom_to[1], zoom_to[2])
  C <- 40075016.686   # ~ circumference of Earth in meters
  x_span <- C / 2^zoom_level
  y_span <- C / 2^(zoom_level+1)
  
  zoom_to_xy <- st_transform(st_sfc(st_point(zoom_to), crs = 4326),
                             crs = target_crs)
  disp_window <- st_sfc(
    st_point(st_coordinates(zoom_to_xy - c(x_span/2, y_span/2))),
    st_point(st_coordinates(zoom_to_xy + c(x_span *1.2, y_span *1.25))),
    crs = target_crs
  )
  
  map <- ggplot() + geom_sf(data = worldmap, aes(fill = deme), alpha = 0.9, size = 0) +
    scale_x_continuous(breaks = seq(-100, 200, by = 10)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    #geom_sf(data = zoom_to_xy, color = 'red') +
    coord_sf(xlim = st_coordinates(disp_window)[,'X'],
             ylim = st_coordinates(disp_window)[,'Y'],
             crs = target_crs, expand = FALSE) +
    scale_fill_manual(values = dcolors, na.value = "grey90") +
    #scale_fill_manual(values = c("#E64B35FF", "grey50", "grey50", "grey50", "grey50", "grey50"), na.value = "grey90") +
    #ggsci::scale_fill_tron()
    theme_void() +
    theme(panel.grid.major = element_line(colour = "grey", size = 0.1),
          legend.position = "none")
  # panel.ontop = TRUE,
  # panel.background = element_blank())
  
  return(map)
}

read_trajectories_dt <- function(traj_file, burnin = 0.1, subsample_n = NULL,
                              type_names, mrs) {
  
  dt <- fread(
    traj_file,
    fill = T,
    select = c("Sample", "type", "type2", "variable", "value", "age"),
    showProgress = TRUE
  )
  
  if (burnin > 0) {
  burnin_from <- max(dt$Sample) * burnin
  dt <- dt[Sample >= burnin_from]
  }
  
  if (!is.null(subsample_n)) {
    s <- sample(unique(dt$Sample), subsample_n)
    dt <- dt[Sample %in% s]
  }
  
  var_map <- c(
    B="infection",
    D="death/recovery/isolation",
    M="migration",
    O="origin",
    S="sampling"
  )
  
  dt[, var_name := var_map[variable]]
  
  dt[, deme := factor(type, levels=0:4, labels=type_names)]
  dt[, deme_to := factor(type2, levels=0:4, labels=type_names)]
  
  dt[, date := to_date(age, mrs)]
  
  setkey(dt, Sample, type, type2, variable, date)
  
  dt[!is.na(variable)]
  
}


bingrid_trajectories_dt <- function(dt, bin = "day") {
  
  dt <- as.data.table(dt)
  dt[, date := lubridate::ymd(date)]
  
  dt[, bin := lubridate::floor_date(date, unit = bin)]
  
  dt <- dt[, .(value = sum(value)),
           by = .(Sample, deme, deme_to, variable, analysis, dataset, bin)]
  
  traj_ids <- unique(dt[, .(Sample, deme, deme_to, variable, analysis, dataset)])
  
  all_bins <- seq(min(dt$bin), max(dt$bin), by = bin)
  
  dt_expanded <- traj_ids[, .(bin = all_bins), 
                          by = .(Sample, deme, deme_to, variable, analysis, dataset)]
  
  dt_expanded <- dt[dt_expanded, on = .(Sample, deme, deme_to, variable, analysis, dataset, bin)]
  
  dt_expanded[is.na(value), value := 0]
  
  setorder(dt_expanded, Sample, analysis, dataset, deme, deme_to, variable, bin)
  dt_expanded[, cumvalue := cumsum(value),
              by = .(Sample, deme, deme_to, variable, analysis, dataset)]
  
  setnames(dt_expanded, "bin", "date")
  
  dt_expanded
  
}


# Function to summarize trajectories per date/deme combination
summarise_bintrajectories_dt <- function(dt) {
  
  dt <- as.data.table(dt)
  
  dt_summary <- dt[, .(
    # Value summaries
    value_low   = quantile(value, 0.025, na.rm = TRUE),
    value_med   = quantile(value, 0.5, na.rm = TRUE),
    value_mean  = mean(value, na.rm = TRUE),
    value_high  = quantile(value, 0.975, na.rm = TRUE),
    value_hpd_low  = boa::boa.hpd(value, alpha = 0.05)["Lower Bound"],
    value_hpd_high = boa::boa.hpd(value, alpha = 0.05)["Upper Bound"],
    
    # Cumulative value summaries
    cumvalue_low   = quantile(cumvalue, 0.025, na.rm = TRUE),
    cumvalue_med   = quantile(cumvalue, 0.5, na.rm = TRUE),
    cumvalue_mean  = mean(cumvalue, na.rm = TRUE),
    cumvalue_high  = quantile(cumvalue, 0.975, na.rm = TRUE),
    cumvalue_hpd_low  = boa::boa.hpd(cumvalue, alpha = 0.05)["Lower Bound"],
    cumvalue_hpd_high = boa::boa.hpd(cumvalue, alpha = 0.05)["Upper Bound"]
    
  ), by = .(date, deme, deme_to, variable, analysis, dataset)]
  
  return(dt_summary)
}


get_ECDCcases <- function(file, demes) {
  cases <- read_tsv(file)
  
  cases_ansys <- cases %>%
    mutate(country = countriesAndTerritories) %>%
    mutate(date = ymd(paste0(year, "-", month, "-", day))) %>%
    filter(date <= analysis_to, date >= analysis_from,
           country %in% demes | continentExp == "Europe") %>%
    left_join(tibble(deme = demes, country = demes)) %>%
    replace_na(list(deme = "OtherEU"))  %>%
    select(date, deme, country, cases, deaths) %>%
    group_by(deme, date) %>%
    summarise(cases = sum(cases),
              deaths = sum(deaths), .groups = "drop_last") %>%
    arrange(date) %>%
    mutate(total_cases = cumsum(cases),
           total_deaths = cumsum(deaths))
  
  return(cases_ansys)
}


make_pretty <- function(df) {
  df %>% mutate(analysis_pretty = factor(case_when(
    analysis == "no_GLM" ~ "non-GLM analysis",
    analysis == "mGLM_6p7e_sumprior_lowerboundFR0" ~ "migration-GLM analysis",
    analysis == "mGLM_6p7e_sumprior_spGLMid" ~ "migration and sampling-GLM analysis"),
    levels = c("non-GLM analysis", "migration-GLM analysis", "migration and sampling-GLM analysis" )),
    dataset_pretty = paste0("S", str_split(dataset, "\\.", simplify = T)[,2]))
}

