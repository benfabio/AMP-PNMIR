
##### 11/04/2017: R Script to examine the content and structure of the ZooScan data from the AMP PMMI project © Fabio Benedetti, OOV-UMS, AFB
##### Aims to:

#	- Map the SST and Chl-a from satellite data (AQUAMODIS, L2 product) with pnmir transects for all dates
#	- Products can be downloaded here: https://oceancolor.gsfc.nasa.gov/cgi/browse.pl (L3 products, re-gridded)
#	- or here: https://oceancolor.gsfc.nasa.gov/cgi/browse.pl?sub=level3&per=CU&day=17267&set=10&ndx=0&mon=17226&sen=am&rad=0&frc=0&dnm=D@M&prm=SST
#	- Explore 8 days composite (considering the cloud coverage of your zone of interest)
#	- Explore salinity products (might be better to locate the front in July for instance)
#	- Explore sea surface anomalies with products from AVISO (M SLA and M ADT)

### Latest update: 06/06/2017

library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("Hmisc")
library("ncdf4")
library("raster")
library("oce")
require("plyr")
library("fields")
library("akima")
library("reshape2")

# MATLAB jet color gradient
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Load coastline:
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data")
cl <- read.csv("gshhg_h_2017-04-04_17-50-04.csv")
coast <- list(
 	 	# the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	  	geom_polygon(aes(x=lon, y=lat), data= cl, fill = "grey45"),
  	  	geom_path(aes(x=lon, y=lat), data= cl, colour = "black", linetype = 1),
  	  	# appropriate projection
  	  	coord_map(),
  	  	# remove extra space around the coast
  	  	scale_x_continuous(name = "Longitude"), #breaks = c(0,10,20,30), labels= c("0°","10°E","20°E","30°E"), expand=c(0,0)), 
  		scale_y_continuous(name = "Latitude"), #breaks = c(30,35,40,45), labels= c("30°N","35°N","40°N","45°N"), expand=c(0,0)),
  		# dark gray background for the panel and legend
 	   	theme(
	   	 	panel.background = element_rect(fill="white"),  # background
	    		legend.key = element_rect(fill="black"),
	    		panel.grid.major = element_line(colour="white")
  		)
		
) # eo coast

# PNMIR stations coordinates
coords_stations <- get(load("coords_stations_pnmir.Rdata"))

quartz()
ggplot() + geom_point(aes(x=x, y=y), data = coords_stations, pch = 21, fill = "#d53e4f", size = 3) + coast + theme_linedraw()
# OKAY.


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### Example, plot the SST and Chl-a (4µm, L2, nightime) for the 04/06/2013
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
dir()

### Load netcdf
nc <- nc_open("A2013155025000.L2_LAC_SST4.nc")
# print(nc)
names(nc$var) # Names of the variables contain some attributes
### Get coords of the satellite pixels
lon <- ncvar_get(nc, "navigation_data/longitude")
# dim(lon) # 1354 2030
lon <- melt(lon)
lat <- ncvar_get(nc, "navigation_data/latitude")
# dim(lat) # 1354 2030
lat <- melt(lat)

### And sea surface temperature ?
#str(nc$var)
temp <- ncvar_get(nc, "geophysical_data/sst4")
# dim(temp); class(temp)
temp <- melt(temp)
head(temp) ; dim(temp)

# Restrict to pixels of interest:
t <- na.omit(data.frame(lon = lon[,3], lat = lat[,3], sst = temp[,3]))
colnames(t)[1:3] <- c("x","y","SST")
# summary(t)
t <- t[which(t$x >= -5.5),]
t <- t[which(t$x <= -4.0),]
t <- t[which(t$y >= 47.9),]
t <- t[which(t$y <= 48.8),]
# Map sst + stations
quartz()
ggplot() +  geom_raster(aes(x = x, y = y, fill = SST), data = t) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_fill_distiller(name = "SST (°C)", palette = "RdYlBu", limits = c(10,15)) + 
			#scale_colour_gradientn(name = "SST (°C)", colours = jet.colors(7)) + 
			coast + coord_quickmap() + theme_linedraw()

### To make a nicer map (interpolation)
require("plyr")
require("fields")
require("akima")
require("reshape2")
# Define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(t$x) + prec/2, to = max(t$x) - prec/2, by = 0.01),
   	 	y <- seq(from = min(t$y) + prec/2, to = max(t$y) - prec/2, by = 0.01)
) # eo list

# Linear interpolation
i <- interp(x = t$x, y = t$y, z = t$SST, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# -> less artifacts near the coast
# Smoothing for better map
i <- image.smooth(i, theta = 0.015)
# Convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "sst"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]
### Nice map ? Examine various palettes...beware of the sst scale
quartz()
ggplot() +  geom_raster(aes(x = lon, y = lat, fill = sst), data = out) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7), limits = c(10,15)) + 
			coast + coord_quickmap() + theme_linedraw()
			
gc()
rm(out, i, grid, prec, t, temp)


### And for Chl-a now:
nc <- nc_open("A2013155140000.L2_LAC_OC.nc")
names(nc$var) # "geophysical_data/chlor_a"
# float geophysical_data/chlor_a[pixels_per_line,number_of_lines]   (Chunking: [452,21])  (Compression: level 4)
# long_name: Chlorophyll Concentration, OCI Algorithm
# units: mg m^-3
# standard_name: mass_concentration_chlorophyll_concentration_in_sea_water
# _FillValue: -32767
# valid_min: 0.00100000004749745
# valid_max: 100
# reference: Hu, C., Lee Z., and Franz, B.A. (2012). Chlorophyll-a algorithms for oligotrophic oceans: A novel approach based on three-band reflectance difference, 
#                J. Geophys. Res., 117, C01011, doi:10.1029/2011JC007395.
			
lon <- ncvar_get(nc, "navigation_data/longitude")
# dim(lon) # 1354 2030
lon <- melt(lon)
lat <- ncvar_get(nc, "navigation_data/latitude")
# dim(lat) # 1354 2030
lat <- melt(lat)

chla <- ncvar_get(nc, "geophysical_data/chlor_a")
# dim(chla); class(chla)
chla <- melt(chla)
head(chla) ; dim(chla)

# Restrict to pixels of interest:
c <- na.omit(data.frame(lon = lon[,3], lat = lat[,3], sst = chla[,3]))
colnames(c)[1:3] <- c("x","y","Chla")
# summary(t)
c <- c[which(c$x >= -5.5),]
c <- c[which(c$x <= -4.0),]
c <- c[which(c$y >= 47.9),]
c <- c[which(c$y <= 48.8),]
# Map sst + stations
quartz()
ggplot() +  geom_point(aes(x = x, y = y, colour = log(Chla)), data = c) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_colour_distiller(name = "Chlorophyll concentration\nOCI Algorithm (mg m^-3; log)", palette = "YlGnBu") + 
			#scale_colour_gradientn(name = "SST (°C)", colours = jet.colors(7)) + 
			coast + coord_quickmap() + theme_linedraw()

### To make a nicer map (interpolation)
require("plyr")
library("fields")
library("akima")
library("reshape2")
# Define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(c$x) + prec/2, to = max(c$x) - prec/2, by = 0.01),
		y <- seq(from = min(c$y) + prec/2, to = max(c$y) - prec/2, by = 0.01)
) # eo list

# Linear interpolation
i <- interp(x = c$x, y = c$y, z = c$Chla, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# Smoothing for better map
i <- image.smooth(i, theta = 0.015)
# Convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "chla"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]
### Nice map ? Examine various palettes...beware of the sst scale
quartz()
ggplot() +  geom_raster(aes(x = lon, y = lat, fill = log(chla)), data = out) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_fill_distiller(name = "Chlorophyll a concentration\n(mg m^-3; log)", palette = "YlGnBu" ) + 
			#scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7), limits = c(11,14)) + 
			coast + coord_quickmap() + theme_linedraw()
			
gc()
rm(out, i, grid, prec, t, temp)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			
### Do it for all the dates with pnmir cruises: 
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/satellite/")
WD <- getwd()
dates <- c("02_05_2011","04_07_2011","12_10_2011","02_05_2012","04_06_2013","09_07_2015","09_10_2015")

d <- "09_10_2015"

for(d in dates) {
	
		# Go to date directory
		paste("Doing ",d, sep = "")
		setwd(paste(WD,"/",d,"/", sep = ""))
		
		##### Get SST ncdf and map
		ncname <- dir()[grep("SST", dir())]
		nc <- nc_open(ncname)
		# Extract data from nc
		lon <- ncvar_get(nc, "navigation_data/longitude")
		# dim(lon) # 1354 2030
		lon <- melt(lon)
		lat <- ncvar_get(nc, "navigation_data/latitude")
		# dim(lat) # 1354 2030
		lat <- melt(lat)
		temp <- ncvar_get(nc, "geophysical_data/sst4")
		# dim(temp); class(temp)
		temp <- melt(temp)
		# head(temp) ; dim(temp)

		# Restrict to pixels of interest:
		t <- na.omit(data.frame(lon = lon[,3], lat = lat[,3], sst = temp[,3]))
		colnames(t)[1:3] <- c("x","y","SST")
		# summary(t)
		t <- t[which(t$x >= -5.5),]
		t <- t[which(t$x <= -4.0),]
		t <- t[which(t$y >= 47.9),]
		t <- t[which(t$y <= 48.8),]
		# dim(t) ; summary(t)
		#### Check basic map
		quartz()
		ggplot() +  geom_point(aes(x = x, y = y, colour = SST), data = t) + 
					geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
					scale_colour_distiller(name = "SST (°C)", palette = "RdYlBu" ) + 
					#scale_colour_gradientn(name = "SST (°C)", colours = jet.colors(7)) + 
					coast + coord_quickmap() + theme_linedraw()
		
		# Define resulting interpolation grid
		prec <- 0.01
		grid <- list( 	
				x <- seq(from = min(t$x) + prec/2, to = max(t$x) - prec/2, by = 0.01),
		   	 	y <- seq(from = min(t$y) + prec/2, to = max(t$y) - prec/2, by = 0.01)
		) # eo list

		# Linear interpolation
		i <- interp(x = t$x, y = t$y, z = t$SST, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
		# -> less artifacts near the coast
		# Smoothing for better map
		i <- image.smooth(i, theta = 0.015)
		# Convert to data.frame
		out <- melt(i$z, varnames = c("start_lon", "start_lat") )
		names(out)[3] <- "sst"
		out$lon <- i$x[out$start_lon]
		out$lat <- i$y[out$start_lat]
		### Nice map ? Examine various palettes...beware of the sst scale
		quartz()
		map <- ggplot() +  geom_raster(aes(x = lon, y = lat, fill = sst), data = out) + 
					geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
					#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
					scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7), limits = c(11,18) ) + #, limits = c(11,17) ) + 
					coast + coord_quickmap() + theme_linedraw()
		
		# Go back to dir() and save
		setwd(WD)
		ggsave(plot = map, file = ".pdf", dpi = 300, height = 6, width = 6)
		rm(nc, out, i, grid, prec, t, temp)
		
		##### Get Chl-a ncdf and map
		setwd(paste(WD,"/",d,"/", sep = ""))
		ncname <- dir()[grep("OC", dir())]
		nc <- nc_open(ncname)
		lon <- ncvar_get(nc, "navigation_data/longitude")
		# dim(lon) # 1354 2030
		lon <- melt(lon)
		lat <- ncvar_get(nc, "navigation_data/latitude")
		# dim(lat) # 1354 2030
		lat <- melt(lat)
		chla <- ncvar_get(nc, "geophysical_data/chlor_a")
		# dim(chla); class(chla)
		chla <- melt(chla)
		# head(chla) ; dim(chla)

		# Restrict to pixels of interest:
		c <- na.omit(data.frame(lon = lon[,3], lat = lat[,3], sst = chla[,3]))
		colnames(c)[1:3] <- c("x","y","Chla")
		# summary(c)
		c <- c[which(c$x >= -5.5),]
		c <- c[which(c$x <= -4.0),]
		c <- c[which(c$y >= 47.9),]
		c <- c[which(c$y <= 48.8),]
		# dim(c) ; summary(c)
		# Map sst + stations
		quartz()
		ggplot() +  geom_point(aes(x = x, y = y, colour = log(Chla)), data = c) + 
					geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
					scale_colour_distiller(name = "Chlorophyll concentration\nOCI Algorithm (mg m^-3)", palette = "YlGnBu") + 
					#scale_colour_gradientn(name = "SST (°C)", colours = jet.colors(7)) + 
					coast + coord_quickmap() + theme_linedraw()

		### To make a nicer map (interpolation)
		require("plyr")
		library("fields")
		library("akima")
		library("reshape2")
		# Define resulting interpolation grid
		prec <- 0.01
		grid <- list( 	
				x <- seq(from = min(c$x) + prec/2, to = max(c$x) - prec/2, by = 0.01),
				y <- seq(from = min(c$y) + prec/2, to = max(c$y) - prec/2, by = 0.01)
		) # eo list

		# Linear interpolation
		i <- interp(x = c$x, y = c$y, z = c$Chla, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
		# Smoothing for better map
		i <- image.smooth(i, theta = 0.015)
		# Convert to data.frame
		out <- melt(i$z, varnames = c("start_lon", "start_lat") )
		names(out)[3] <- "chla"
		out$lon <- i$x[out$start_lon]
		out$lat <- i$y[out$start_lat]
		### Nice map ? Examine various palettes...beware of the sst scale
		quartz()
		ggplot() +  geom_raster(aes(x = lon, y = lat, fill = log(chla)), data = out) + 
					geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
					scale_fill_distiller(name = "Chlorophyll a concentration\n(mg m^-3; log)", palette = "YlGnBu" ) + 
					#scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7), limits = c(11,14)) + 
					coast + coord_quickmap() + theme_linedraw()
		
		# Go back to dir() and save
		setwd(WD)
		ggsave(plot = map, file = ".pdf", dpi = 300, height = 6, width = 6)
		
		# Clean
		gc()
		rm(nc)
	
} # eo for d in dates loop



# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##### 11/04/2017: Explore a MARS-3D re-analysis output, we might need to use those...

d <- "09_07_2015"
setwd(paste(WD,"/",d,"/", sep = ""))
nc <- nc_open("PREVIMER_F1-MARS3D-MANGAE2500_20150709T1200Z_FILTRE.nc")
names(nc$var) # Has many data :-)
# Extract data from nc
lon <- ncvar_get(nc, "longitude")
dim(lon) # 822 624
lon <- melt(lon)
lat <- ncvar_get(nc, "latitude")
dim(lat) # 822 624
lat <- melt(lat)
temp <- ncvar_get(nc, "TEMP")
# dim(temp); class(temp) ### Has 40 vertical levels, keep one of the first
temp <- temp[,,1]
temp <- melt(temp)
# head(temp) ; dim(temp)
# Filter
t <- na.omit(data.frame(lon = lon[,3], lat = lat[,3], sst = temp[,3]))
colnames(t)[1:3] <- c("x","y","SST")
t <- t[which(t$x >= -5.5),]
t <- t[which(t$x <= -4.0),]
t <- t[which(t$y >= 47.9),]
t <- t[which(t$y <= 48.8),]

quartz()
ggplot() +  geom_raster(aes(x = x, y = y, fill = SST), data = t) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_fill_distiller(name = "SST (°C)", palette = "RdYlBu", limits = c(12,20)) + 
			#scale_colour_gradientn(name = "SST (°C)", colours = jet.colors(7)) + 
			coast + coord_quickmap() + theme_linedraw()

### To make a nicer map (interpolation)
require("plyr")
require("fields")
require("akima")
require("reshape2")
# Define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(t$x) + prec/2, to = max(t$x) - prec/2, by = 0.01),
   	 	y <- seq(from = min(t$y) + prec/2, to = max(t$y) - prec/2, by = 0.01)
) # eo list

# Linear interpolation
i <- interp(x = t$x, y = t$y, z = t$SST, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# -> less artifacts near the coast
# Smoothing for better map
i <- image.smooth(i, theta = 0.015)
# Convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "sst"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]
### Nice map ? Examine various palettes...beware of the sst scale
quartz()
ggplot() +  geom_raster(aes(x = lon, y = lat, fill = sst), data = out) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7), limits = c(10,15)) + 
			coast + coord_quickmap() + theme_linedraw()


##### 12/04/17: Explore 8-day composites from the NOAA.
# http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1sstd8day.html

setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/satellite/")
nc <- nc_open("erdMH1sstd8day_2bd6_5536_e7ca.nc")
print(nc)
names(nc$var) # sst only 
sst <- ncvar_get(nc, "sst")
class(sst) ; dim(sst)
str(sst)
sst[1,,] # longitude, latitude, time

m_sst <- melt(sst)
colnames(m_sst)[1:4] <- c("x","y","time","SST")

# Extract lon, lat & time
lon <- ncvar_get(nc, "longitude") # 49 levels
lat <- ncvar_get(nc, "latitude")  # 49 levels
time <- ncvar_get(nc, "time")     # 204 levels --> need to convert seconds since 1970-01-01T00:00:00Z to actual dates !
### Use as.POSIXct()
val <- time[1]
as.POSIXct(val, origin = "1970-01-01")
# "2012-11-04 22:32:00 CST"
rm(val)
dates <- as.Date(as.POSIXct(time, origin = "1970-01-01")) ### Great

### Now, cbind in a single data frame ("m_sst")
# Provide longitudes by repeating 'lon' the adequate number of times
nrow(m_sst)/ length(lon) # 9996
m_sst$x <- rep(lon, 9996)

# Provide latitudes, like that for instance:
# m_sst[which(m_sst$y == 1), "y"] <- lat[1]
unique(m_sst$y)

for(i in unique(m_sst$y)) {
	m_sst[which(m_sst$y == i), "y"] <- lat[i]
} # eo for loop

head(m_sst) # Looks ok...time now

# Try a map first
m_sst <- na.omit( m_sst[which(m_sst$time == 204),] )
m_sst <- m_sst[which(m_sst$x >= -5.5),]
m_sst <- m_sst[which(m_sst$x <= -4.0),]
m_sst <- m_sst[which(m_sst$y >= 47.9),]
m_sst <- m_sst[which(m_sst$y <= 48.8),]

quartz()
ggplot() +  geom_raster(aes(x = x, y = y, fill = SST), data = m_sst) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7) ) + 
			coast + coord_quickmap() + theme_linedraw()
			
# Define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(m_sst$x) + prec/2, to = max(m_sst$x) - prec/2, by = 0.01),
   	 	y <- seq(from = min(m_sst$y) + prec/2, to = max(m_sst$y) - prec/2, by = 0.01)
) # eo list

# Linear interpolation
i <- interp(x = m_sst$x, y = m_sst$y, z = m_sst$SST, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# Smoothing for better map
i <- image.smooth(i, theta = 0.015)
# Convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "sst"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]

quartz()
ggplot() +  geom_raster(aes(x = lon, y = lat, fill = sst), data = out) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7) ) + 
			coast + coord_quickmap() + theme_linedraw()
			
### OK. Find the integers corresponding to the dates of interest, and map SST and Chla for each of them...	
# pnmir cruise dates = c("02_05_2011","04_07_2011","12_10_2011","02_05_2012","04_06_2013","09_07_2015","09_10_2015")		
# 02_05_2011 --> time[1]
### OK

# 04_07_2011 --> "2011-06-30" "2011-07-08" --> time[9]
### OK

# 12_10_2011 --> time[21]
### No data (strange because daily product is OK)

# 02_05_2012 --> time[47]
### OK

# 04_06_2013 --> time[97]
### D1 & D2 not covered

# 09_07_2015 --> time[193]
### OK

# 09_10_2015 --> time[204]
### OK.




##### Chl-a now
nc <- nc_open("erdMH1chla8day_9b13_590d_232d.nc")
print(nc)
names(nc$var) # sst only 
chla <- ncvar_get(nc, "chlorophyll")
class(chla) ; dim(chla)
str(chla)
chla[1,,] # longitude, latitude, time

m_chla <- melt(chla)
colnames(m_chla)[1:4] <- c("x","y","time","Chla")

# Extract lon, lat & time
lon <- ncvar_get(nc, "longitude") # 49 levels
lat <- ncvar_get(nc, "latitude")  # 49 levels
time <- ncvar_get(nc, "time")     # 204 levels --> need to convert seconds since 1970-01-01T00:00:00Z to actual dates !
### Use as.POSIXct()
val <- time[1]
as.POSIXct(val, origin = "1970-01-01")
# "2012-11-04 22:32:00 CST"
rm(val)
dates <- as.Date(as.POSIXct(time, origin = "1970-01-01")) ### Great

### Now, cbind in a single data frame ("m_sst")
# Provide longitudes by repeating 'lon' the adequate number of times
nrow(m_chla)/ length(lon) # 9996
m_chla$x <- rep(lon, 9996)

# Provide latitudes, like that for instance:
# m_sst[which(m_sst$y == 1), "y"] <- lat[1]
unique(m_chla$y)

for(i in unique(m_chla$y)) {
	m_chla[which(m_chla$y == i), "y"] <- lat[i]
} # eo for loop


# Try a map first
m_chla <- na.omit( m_chla[which(m_chla$time == 204),] )
m_chla <- m_chla[which(m_chla$x >= -5.5),]
m_chla <- m_chla[which(m_chla$x <= -4.0),]
m_chla <- m_chla[which(m_chla$y >= 47.9),]
m_chla <- m_chla[which(m_chla$y <= 48.8),]

quartz()
ggplot() +  geom_raster(aes(x = x, y = y, fill = Chla), data = m_chla) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_fill_distiller(name = "Chlorophyll a\nconcentration (mg m^-3)", palette = "YlGnBu" ) + 
			#scale_fill_gradientn(name = "SST (°C)", colours = jet.colors(7) ) + 
			coast + coord_quickmap() + theme_linedraw()
			
# Define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(m_chla$x) + prec/2, to = max(m_chla$x) - prec/2, by = 0.01),
   	 	y <- seq(from = min(m_chla$y) + prec/2, to = max(m_chla$y) - prec/2, by = 0.01)
) # eo list

# Linear interpolation
i <- interp(x = m_chla$x, y = m_chla$y, z = m_chla$Chla, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# Smoothing for better map
i <- image.smooth(i, theta = 0.015)
# Convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "chla"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]

quartz()
ggplot() +  geom_raster(aes(x = lon, y = lat, fill = chla), data = out) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_distiller(name = "Chlorophyll a\nconcentration\n(mg m^-3)", palette = "YlGnBu" ) + 
			coast + coord_quickmap() + theme_linedraw()
			
# 02_05_2011 --> time[1]
### Doesn't cover coasts

# 04_07_2011 --> time[9]
### Not OK, few data

# 12_10_2011 --> time[21]
### Nope, lot of missing data

# 02_05_2012 --> time[47]
### decent...

# 04_06_2013 --> time[97]
### OK.

# 09_07_2015 --> time[193]
### OK.

# 09_10_2015 --> time[204]
### OK.



##### 12/04/17: Explore netcdf from AVISO, downloaded @:
# http://www.aviso.altimetry.fr/en/data/data-access/gridded-data-extraction-tool.html
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/satellite/")
nc <- nc_open("dataset-duacs-dt-global-allsat-madt-h_1492012872169.nc")
print(nc)
names(nc$var) # adt = Absolute Dynamic Topography (sea_surface_height_above_geoid, meters) 
adt <- ncvar_get(nc, "adt")
class(adt) ; dim(adt)
str(adt)
adt[1,1,1000] # lon, lat, time

madt <- melt(adt)
colnames(madt)[1:4] <- c("x","y","time","adt")

# Extract lon, lat & time
lon <- ncvar_get(nc, "lon") # 49 levels
lat <- ncvar_get(nc, "lat")  # 49 levels
time <- ncvar_get(nc, "time")  # 1622 levels --> need to convert DAYS since 1950-01-01T00:00:00Z to actual dates ! SO need to multiply days by 24 and 3600 to obtain hours since 1950...
### Use as.POSIXct()
dates <- as.Date(as.POSIXct((time*24*3600), origin = "1950-01-01 00:00:00"))
# OK.
# Provide longitudes by repeating 'lon' the adequate number of times
nrow(madt)/ length(lon) # 8110
madt$x <- rep(lon, 8110)
madt$x <- (madt$x - 360)
# Provide latitudes, like that for instance:
# m_sst[which(m_sst$y == 1), "y"] <- lat[1]
unique(madt$y)

for(i in unique(madt$y)) {
	madt[which(madt$y == i), "y"] <- lat[i]
} # eo for loop


# Try a map first, choose a proper date
mmadt <- na.omit( madt[which(madt$time == 1556),] )
quartz()
ggplot() +  geom_raster(aes(x = x, y = y, fill = adt), data = mmadt) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_distiller(name = "Absolute\nDynamic\nTopography", palette = "YlGnBu" ) + 
			coast + coord_quickmap() + theme_linedraw()


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### 25/04/2017: Explore satellite product given by M. Sourisseau (IFREMER, Brest) - should be the ones available @ PREVIMER - MARC website (only as maps though)

# setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/satellite/")
nc <- nc_open("20110704-EUR-L4-SPIM-ATL-v01-fv01-OI.nc")
print(nc)
names(nc$var)

### Get data
chla <- ncvar_get(nc, "analysed_chl_a")
dim(chla)
# Extract lon, lat & time
lon <- ncvar_get(nc, "lon") # 1667
lat <- ncvar_get(nc, "lat")  # 2401
time <- ncvar_get(nc, "time")  # 1 ; seconds since 1998-01-01 00:00:00
# Melt chla data
mchla <- melt(chla)
dim(mchla)
head(mchla)
colnames(mchla)[1:3] <- c("x", "y", "Chla")

# Provide 
# nrow(mchla)/ length(lon) # 2401
mchla$x <- rep(lon, nrow(mchla)/ length(lon))

# Provide latitudes
unique(mchla$y)
# i <- unique(mchla$y)[2]

for(i in unique(mchla$y) ) {
		mchla[which(mchla$y == i), "y"] <- lat[i]	
} # eo for loop

# Check
summary(mchla) # OK.

# Restrict to zone of interest
mchla <- mchla[which(mchla$x >= -5.5),]
mchla <- mchla[which(mchla$x <= -4.0),]
mchla <- mchla[which(mchla$y >= 47.9),]
mchla <- mchla[which(mchla$y <= 48.8),]

quartz()
ggplot() +  geom_raster(aes(x = x, y = y, fill = log(Chla)), data = mchla) + 
			#geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_fill_distiller(name = "Chlorophyll_a\n(log(mg.m-3))", palette = "YlGnBu") + 
			coast + coord_quickmap() + theme_linedraw()


### Looks good !

##### 06/06/2017: From the link given by Marc Sourisseau, explore directories to get SST data
# ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic/modis/2011/
# ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/ocean-color/atlantic/viirs/2013/
nc <- nc_open("20130515.nc") # 20130515 --> VIIRS
print(nc)
names(nc$var)






