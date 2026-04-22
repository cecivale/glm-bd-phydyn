#Bethany Allen  28th October 2024
#Code to quantify the amount of geographic area sampled through time

#Load packages
library(tidyverse)
library(dggridR)

#Specify ages of change times
ages <- c(66.0, 89.8, 121.4, 145.0, 161.5, 174.7, 201.4, 237.0, 251.9)
bin_names <- c("Late_Cretaceous", "Mid_Cretaceous",
               "Early_Cretaceous", "Late_Jurassic",
               "Mid_Jurassic", "Early_Jurassic",
               "Late_Triassic", "EarlyMid_Triassic")

#Construct equal area hex grid with cells approximately 1000km in area
dggs <- dgconstruct(area = 1000, metric = TRUE, resround = "down")

#Read in data
collections <- read_csv("data/collections.csv", skip = 19)

#Loop through and determine which time bins each collection crosses
for (j in 1:(length(ages) - 1)) {
  vector <- c()
  for (i in 1:nrow(collections)) {
    if (collections$max_ma[i] >= ages[j] && collections$min_ma[i] <= ages[j+1])
    {vector[i] <- 1} else {vector[i] <- 0}
  }
  collections[,j+16] <- vector
}

#Label columns with time bin names
colnames(collections)[17:24] <- bin_names

#Get the corresponding grid cells for each collection
collections$cell <- dgGEO_to_SEQNUM(dggs,collections$lng,collections$lat)$seqnum

#Filter and count unique formations (coded inefficiently to allow checking)
bin1 <- filter(collections, Late_Cretaceous == 1)
count1 <- length(unique(bin1$cell))

bin2 <- filter(collections, Mid_Cretaceous == 1)
count2 <- length(unique(bin2$cell))

bin3 <- filter(collections, Early_Cretaceous == 1)
count3 <- length(unique(bin3$cell))

bin4 <- filter(collections, Late_Jurassic == 1)
count4 <- length(unique(bin4$cell))

bin5 <- filter(collections, Mid_Jurassic == 1)
count5 <- length(unique(bin5$cell))

bin6 <- filter(collections, Early_Jurassic == 1)
count6 <- length(unique(bin6$cell))

bin7 <- filter(collections, Late_Triassic == 1)
count7 <- length(unique(bin7$cell))

bin8 <- filter(collections, EarlyMid_Triassic == 1)
count8 <- length(unique(bin8$cell))

#String together
print(c(count1, count2, count3, count4, count5, count6, count7, count8))
