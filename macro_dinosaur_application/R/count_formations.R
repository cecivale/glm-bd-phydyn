#Bethany Allen  28th October 2024
#Code to quantify number of collections through time from PBDB data

#Load packages
library(tidyverse)
library(palaeoverse)

#Read in data
formations <- read_csv("data/strata.csv", skip = 17)

#Remove data with no formation, formation in quote marks, or a slash included
formations <- filter(formations, !is.na(formation))
formations <- formations[-c(1:3),]
formations <- formations[-which(str_detect(formations$formation, "/")),]

#Use tax_check to test for similar formation names
tax_check(formations, name = "formation", dis = 0.1)

#Correct errors
formations$formation[which(formations$formation == "Ashfield Shales")] <-
  "Ashfield Shale"
formations$formation[which(formations$formation == "Broome Sandstrone")] <-
  "Broome Sandstone"
formations$formation[which(formations$formation == "Bukobay")] <-
  "Bukobai"
formations$formation[which(formations$formation == "Couches du Chailley")] <-
  "Couches de Chailley"
formations$formation[which(formations$formation == "Calcaires de Rognac")] <-
  "Calcaire de Rognac"
formations$formation[which(formations$formation == "Eschuca")] <-
  "Escucha"
formations$formation[which(formations$formation == "Fix Hills")] <-
  "Fox Hills"
formations$formation[which(formations$formation == "Gailard")] <-
  "Gaillard"
formations$formation[which(formations$formation == "Grès Supérieurs")] <-
  "Grès supérieurs"
formations$formation[which(formations$formation == "Horsehoe Canyon")] <-
  "Horseshoe Canyon"
formations$formation[which(formations$formation == "Intertrappean Beds")] <-
  "Intertrappean"
formations$formation[which(formations$formation == "Ischigalasto")] <-
  "Ischigualasto"
formations$formation[which(formations$formation == "Jingchuan")] <-
  "Jianchuan"
formations$formation[which(formations$formation == "Jiagdihe")] <-
  "Jiangdihe"
formations$formation[which(formations$formation == "Kalazha")] <-
  "Kalaza"
formations$formation[which(formations$formation == "Kopanaskya")] <-
  "Kopanskaya"
formations$formation[which(formations$formation == "Lefipan")] <-
  "Lefipán"
formations$formation[which(formations$formation == "Lower Garumnian")] <-
  "Lower Red Garumnian"
formations$formation[which(formations$formation == "Nelsen")] <-
  "Neslen"
formations$formation[which(formations$formation == "Petropavlovka")] <-
  "Petropavlovskaya"
formations$formation[which(formations$formation == "Purple Claystones")] <-
  "Purple Claystone"
formations$formation[which(formations$formation == "Quarziti de Monte Serra")] <-
  "Quarziti di Monte Serra"
formations$formation[which(formations$formation == "Shinkhuduk")] <-
  "Shin-Khuduk"
formations$formation[which(formations$formation == "Steieresdorf")] <-
  "Steierdorf"
formations$formation[which(formations$formation == "Sânpetru")] <-
  "Sînpetru"
formations$formation[which(formations$formation == "Talynzhansk")] <-
  "Talynzhan"
formations$formation[which(formations$formation == "Timezgadiwine")] <-
  "Timezgadiouine"
formations$formation[which(formations$formation == "Tootebuc")] <-
  "Toolebuc"
formations$formation[which(formations$formation == "Ulan-Ereg")] <-
  "Ulaan-Ereg"
formations$formation[which(formations$formation == "Ulan Argalant")] <-
  "Ulan-Argalant"
formations$formation[which(formations$formation == "Vokhmian")] <-
  "Vokhma"
formations$formation[which(formations$formation == "Yanan")] <-
  "Yan'an"
formations$formation[which(formations$formation == "Yenchang")] <-
  "Yanchang"

#Specify ages of change times
ages <- c(66.0, 89.8, 121.4, 145.0, 161.5, 174.7, 201.4, 237.0, 251.9)
bin_names <- c("Late_Cretaceous", "Mid_Cretaceous",
               "Early_Cretaceous", "Late_Jurassic",
               "Mid_Jurassic", "Early_Jurassic",
               "Late_Triassic", "EarlyMid_Triassic")

#Loop through and determine which time bins each formation crosses
for (j in 1:(length(ages) - 1)) {
  vector <- c()
  for (i in 1:nrow(formations)) {
    if (formations$max_ma[i] >= ages[j] && formations$min_ma[i] <= ages[j+1])
      {vector[i] <- 1} else {vector[i] <- 0}
  }
  formations[,j+10] <- vector
}

#Label columns with time bin names
colnames(formations)[11:18] <- bin_names

#Filter and count unique formations (coded inefficiently to allow checking)
bin1 <- filter(formations, Late_Cretaceous == 1)
sort(unique(bin1$formation))
count1 <- length(unique(bin1$formation))

bin2 <- filter(formations, Mid_Cretaceous == 1)
sort(unique(bin2$formation))
count2 <- length(unique(bin2$formation))

bin3 <- filter(formations, Early_Cretaceous == 1)
sort(unique(bin3$formation))
count3 <- length(unique(bin3$formation))

bin4 <- filter(formations, Late_Jurassic == 1)
sort(unique(bin4$formation))
count4 <- length(unique(bin4$formation))

bin5 <- filter(formations, Mid_Jurassic == 1)
sort(unique(bin5$formation))
count5 <- length(unique(bin5$formation))

bin6 <- filter(formations, Early_Jurassic == 1)
sort(unique(bin6$formation))
count6 <- length(unique(bin6$formation))

bin7 <- filter(formations, Late_Triassic == 1)
sort(unique(bin7$formation))
count7 <- length(unique(bin7$formation))

bin8 <- filter(formations, EarlyMid_Triassic == 1)
sort(unique(bin8$formation))
count8 <- length(unique(bin8$formation))

#String together
print(c(count1, count2, count3, count4, count5, count6, count7, count8))
