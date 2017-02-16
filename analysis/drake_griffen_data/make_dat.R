# run this in the data-and-code directory of the zip file downloaded
# from http://datadryad.org/resource/doi:10.5061/dryad.q3p64

#Preprocess raw time series data for analysis
source('preprocess.R')

library(lubridate)

# data munging
timeseries$date <- mdy(as.character(timeseries$Date.))
populations$Nhat <- populations$x
populations$x <- NULL

# table of which groups were deteriorating or not
groups <- unique(timeseries[,c("ID", "Deteriorating")])

populations$deteriorating <- populations$ID %in% groups$ID[groups$Deteriorating==1]

# smooth of population size over time
pop_happy <- subset(populations, !deteriorating)
pop_unhappy <- subset(populations, deteriorating)


save.image(file="drake_griffen.RData")

