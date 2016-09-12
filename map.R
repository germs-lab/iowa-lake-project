library(ggplot2)

library(sp)
library(rgdal)
setwd("~/Box Sync/Github/germs-lab/iowa_lake_project/")
gps14 <- read.csv("lake_survey_coords_utm_zone_14.csv")
gps15 <- read.csv("lake_survey_coords_utm_zone_15.csv")

coordinates(gps14) <- c("utm.east","utm.north")
coordinates(gps15) <- c("utm.east","utm.north")

sputm14 <- SpatialPoints(gps14, proj4string=CRS("+proj=utm +zone=14 +datum=WGS84"))
sputm15<- SpatialPoints(gps15, proj4string=CRS("+proj=utm +zone=15 +datum=WGS84"))

spgeo14<- spTransform(sputm14, CRS("+proj=longlat +datum=WGS84"))
spgeo15<- spTransform(sputm15, CRS("+proj=longlat +datum=WGS84"))

gps <- data.frame(rbind(spgeo14,spgeo15))

#install.packages("ggmap")
library(ggmap)
mapImageData <- get_googlemap(center = c(lon = median(gps$utm.east), lat = median(gps$utm.north)), zoom = 6, maptype = c("terrain"))
ggmap(mapImageData, extent = "device") +   geom_point(aes(x = utm.east,  y = utm.north), data = gps, colour = "red", size = 1, pch = 20)