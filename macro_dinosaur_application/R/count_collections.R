#Bethany Allen  28th October 2024
#Code to quantify number of collections through time from PBDB data

#Load packages
library(tidyverse)

#Specify ages of change times
ages <- c(66.0, 89.8, 121.4, 145.0, 161.5, 174.7, 201.4, 237.0, 251.9)
bin_names <- c("Late_Cretaceous", "Mid_Cretaceous",
               "Early_Cretaceous", "Late_Jurassic",
               "Mid_Jurassic", "Early_Jurassic",
               "Late_Triassic", "EarlyMid_Triassic")

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

#Filter and count unique formations (coded inefficiently to allow checking)
bin1 <- filter(collections, Late_Cretaceous == 1)
count1 <- nrow(bin1)

bin2 <- filter(collections, Mid_Cretaceous == 1)
count2 <- nrow(bin2)

bin3 <- filter(collections, Early_Cretaceous == 1)
count3 <- nrow(bin3)

bin4 <- filter(collections, Late_Jurassic == 1)
count4 <- nrow(bin4)

bin5 <- filter(collections, Mid_Jurassic == 1)
count5 <- nrow(bin5)

bin6 <- filter(collections, Early_Jurassic == 1)
count6 <- nrow(bin6)

bin7 <- filter(collections, Late_Triassic == 1)
count7 <- nrow(bin7)

bin8 <- filter(collections, EarlyMid_Triassic == 1)
count8 <- nrow(bin8)

#String together
print(c(count1, count2, count3, count4, count5, count6, count7, count8))
