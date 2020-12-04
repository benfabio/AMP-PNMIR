
##### 18/04/2017: R Script to examine the data acquired during the PNMIR cruises © Fabio Benedetti, OOV-UMS, AFB
##### Aims to:

#	- Analyse time series of hydrobiological data (T, S, Chla etc...) - time series, Principal Component Analysis
#	- Analyse time series of phytoplankton data (counts and percentages; large groups and species-level) - Correspondance Analysis
#	- Analyse the link between both (multivariate analyses)
#	- Plot the environmental conditions of each station for the years: 2011? 2012? 2013, and 2015

### Latest update: 26/04/2017

library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("Hmisc")
library("ncdf4")
library("raster")
library("oce")
require("dplyr")
library("fields")
library("akima")
library("reshape2")
library("FactoMineR")
library("vegan")

# MATLAB jet color gradient
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Load coastline:
# setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data")
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

### Load the data as .csv:
ddf <- read.csv("PNMIR_cruise_data.csv", h = TRUE, dec = ",", sep = ";")
dim(ddf)
str(ddf)
colnames(ddf)
# Make some adjustments:
ddf$date <- paste(ddf$day, ddf$month, ddf$year, sep = "/")
ddf$date <- as.Date(ddf$date, format = "%d/%m/%Y") # should be okay
rownames(ddf) # Same as ID, ok

### 24/04/2017: Correct station coordibates plz...using 'coords_stations'
ddf[which(ddf$station == "b1"), "x"] <- -4.616667 ; ddf[which(ddf$station == "b1"), "y"] <- 48.31667
ddf[which(ddf$station == "b2"), "x"] <- -4.783333 ; ddf[which(ddf$station == "b2"), "y"] <- 48.31667
ddf[which(ddf$station == "b3"), "x"] <- -4.950000 ; ddf[which(ddf$station == "b3"), "y"] <- 48.31667
ddf[which(ddf$station == "b4"), "x"] <- -5.050000 ; ddf[which(ddf$station == "b4"), "y"] <- 48.33333
ddf[which(ddf$station == "b5"), "x"] <- -5.166667 ; ddf[which(ddf$station == "b5"), "y"] <- 48.33333
ddf[which(ddf$station == "b6"), "x"] <- -5.250000 ; ddf[which(ddf$station == "b6"), "y"] <- 48.33333
ddf[which(ddf$station == "b7"), "x"] <- -5.416667 ; ddf[which(ddf$station == "b7"), "y"] <- 48.33333
ddf[which(ddf$station == "d1"), "x"] <- -4.416667 ; ddf[which(ddf$station == "d1"), "y"] <- 48.16667
ddf[which(ddf$station == "d2"), "x"] <- -4.616667 ; ddf[which(ddf$station == "d2"), "y"] <- 48.16667
ddf[which(ddf$station == "d3"), "x"] <- -4.783333 ; ddf[which(ddf$station == "d3"), "y"] <- 48.16667
ddf[which(ddf$station == "d4"), "x"] <- -4.950000 ; ddf[which(ddf$station == "d4"), "y"] <- 48.16667
ddf[which(ddf$station == "d5"), "x"] <- -5.166667 ; ddf[which(ddf$station == "d5"), "y"] <- 48.16667
ddf[which(ddf$station == "d6"), "x"] <- -5.250000 ; ddf[which(ddf$station == "d6"), "y"] <- 48.16667

### Some PCA on environmental data:
ddf$Tdiff <- ddf$Ts - ddf$Tf
res.PCA <- PCA(X = ddf[which(ddf$year %in% c(2011,2012,2013,2015)),c(12,44,15,17,19,20,23,24)], scale.unit = T)
summary(res.PCA) # 71.66% of total variance on the 4 first PCs
### What's station 45, 46 and 106
ddf[45,] # 2 May 2012, b1, not stratified, peak of Guinardia spp. with 6.1 µgL of Chl-a ! 
ddf[46,] # Same as 45 !
ddf[106,] # 9 of July 2015, b7, stratified conditions, 21.3 µg/L of Chla- !! Problem here...

### Re-perform PCA, without 2012 and sample 106
# res.PCA <- PCA(X = ddf[which(ddf$year %in% c(2011,2013,2015) & ddf$Chla < 20 & ddf$PO4 < 4),c(12,14,15,17:20,23,24)], scale.unit = TRUE)
# With just: T, S, main nutrients, and Chla 
res.PCA <- PCA(X = ddf[which(ddf$year %in% c(2011,2012,2013,2015) & ddf$Chla <= 20 & ddf$PO4 <= 4 & ddf$NO3 <= 10 & ddf$SiOH4 <= 10),c("max_depth","Ts","Tdiff","Ss","NO3","SiOH4","PO4","Chla")], scale.unit = T)
				
summary(res.PCA) # 73.85% of total variance on first 4 PCs
plot(res.PCA)
quartz() ; plot.PCA(res.PCA, choix = "var", axes = c(3,4))

# 03/05/2017: 74.88% on 4 PCs
#                        Dim.1   Dim.2   Dim.3   Dim.4   Dim.5   Dim.6   Dim.7
# Variance               2.010   1.605   1.343   0.950   0.847   0.586   0.369
# % of var.             25.126  20.067  16.782  11.878  10.583   7.319   4.612
# Cumulative % of var.  25.126  45.193  61.975  73.853  84.436  91.755  96.367

# Better plot:
library("ggfortify")
library("ggrepel")
RES <- augment(res.PCA, dimensions= c(1,2,3,4), which = "col")
str(RES)
RES
colnames(RES)[1:8] <- c("vars", "type", "PC1", "PC2", "PC3", "PC4", "cos2", "contrib")
circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}
dat <- circleFun(c(0,0), diameter = 2, npoints = 100)
quartz()
ggplot() + 
	geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
	geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2, label= vars), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "firebrick3", data= RES) +
	geom_path(data = dat, aes(x,y) ) + xlab("PC1 (26.50%)") + ylab("PC2 (23.80%)") + theme_bw()
# Same graph but for PC3 and PC4	
quartz()
ggplot() + 
	geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
	geom_segment(aes(x = 0, y = 0, xend = PC3, yend = PC4, label= vars), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "firebrick3", data= RES) +
	geom_path(data = dat, aes(x,y) ) + xlab("PC3 (13.80%)") + ylab("PC4 (10.77%)") + theme_bw()	

### get species coordinates in niche space
# res.PCA$ind$coord
obj <- data.frame(id = rownames(res.PCA$ind$coord), PC1 = res.PCA$ind$coord[,1], 
				  PC2 = res.PCA$ind$coord[,2], PC3 = res.PCA$ind$coord[,3], 
				  PC4 = res.PCA$ind$coord[,4], month = NA, year = NA )
# Provide month based on ids:
for(i in obj$id) {
	obj[which(obj$id == i), "month"] <- ddf[which(ddf$ID == i), "month"]
	obj[which(obj$id == i), "year"] <- ddf[which(ddf$ID == i), "year"]
} # eo for loop			  

### get variables coordinates in another table
# res.PCA$var$coord
vars <- data.frame(vars = rownames(res.PCA$var$coord),
				PC1 = res.PCA$var$coord[,1],
				PC2 = res.PCA$var$coord[,2],
				PC3 = res.PCA$var$coord[,3],
				PC4 = res.PCA$var$coord[,4] )

quartz()
plot <- ggplot() + geom_point(aes(x= PC1, y= PC2, fill = factor(year)), colour = "black", size = 3, data = obj, shape = 21) + 
	 		geom_segment(aes(x=0, y=0, xend= PC1*5, yend= PC2*5), arrow = arrow(length = unit(0.4, "cm"),type = "closed"), colour = "black", data = vars) + 
			scale_fill_manual("Month", values = c("#66c2a5","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
			#scale_shape_manual(name = "Niche traits\ngroups", values = c(21,22,23,24)) + 
	 		geom_text_repel( aes(x= PC1, y= PC2, label= rownames(obj)), size = 2.75, data = obj) +
			#geom_text(aes(x= PC1*5, y= PC2*5, label= vars), hjust = -0.5, vjust = 1, size = 2.5, data = vars) + 
 		   	geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
			xlab(paste("PC1 (",round(res.PCA$eig$per[1], 2),"%)", sep = "")) + ylab(paste("PC2 (",round(res.PCA$eig$per[2], 2),"%)", sep = "")) + theme_bw()

### Very clear seasonal cycle
ggsave(plot = plot, filename = "PCA_pnmir_1_2.pdf", dpi = 300, width = 8, height = 7)

#quartz()
plot <- ggplot() + geom_point(aes(x= PC3, y= PC4, fill = factor(year)), colour = "black", size = 3, data = obj, shape = 21) + 
	 		geom_segment(aes(x=0, y=0, xend= PC3*5, yend= PC4*5), arrow = arrow(length = unit(0.4, "cm"),type = "closed"), colour = "black", data = vars) + 
			scale_fill_manual("Month", values = c("#66c2a5","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
			#scale_shape_manual(name = "Niche traits\ngroups", values = c(21,22,23,24)) + 
	 		geom_text_repel( aes(x= PC3, y= PC4, label= rownames(obj)), size = 2.75, data = obj) +
			#geom_text(aes(x= PC3*5, y= PC4*5, label= vars), hjust = -0.5, vjust = 1, size = 2.5, data = vars) + 
 			geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
			xlab(paste("PC3 (",round(res.PCA$eig$per[3],2),"%)", sep = "")) + ylab(paste("PC4 (",round(res.PCA$eig$per[4], 2),"%)", sep = "")) + theme_bw()

ggsave(plot = plot, filename = "PCA_pnmir_3_4.pdf", dpi = 300, width = 8, height = 7)

# Check stations 125, 74 and 42 & 103
ddf[which(ddf$ID == 45),]
ddf[which(ddf$ID == 125),]
ddf[which(ddf$ID == 103),]
### Both due to a PO4 peak...
ddf[which(ddf$ID == 74),] # Salinity peak 
		
		
##### Now, analyze the variations in phytoplankton community strcuture form counts		
### CA on phytoplankton counts:
res.ca <- CA(X = na.omit(ddf[which(ddf$year %in% c(2011,2012,2013)), c(26:28,40:43)]) )
summary(res.ca)
# 83.77% of total variance wuth 2 CA axes only !

### For ggplotting: need to prepare a few dataframes.
# For descriptors
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
                    Cos2 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2])

head(AFCsp)
# For objects : jellyfishes
data_for_months <- na.omit(ddf[which(ddf$year %in% c(2011,2012,2013)), c(1,4,9,26:28,40:43)])
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
                    Cos2 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2],
					month = data_for_months$month, id = data_for_months$ID, depth = data_for_months$max_depth)

head(AFCst)

quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = depth, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Localisation", values = c(21,22)) + 
	scale_fill_manual("Month", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 1 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()

	
### CA on counts looks WAY better than CA on percentages (same pattern but % of explained var much higher). Keep coordinates on CA for regression with environment data
AFCst
# Supply environmental data: 
AFCst$Temp <- NA
AFCst$Sal <- NA
AFCst$pH <- NA
AFCst$NO3 <- NA
AFCst$NH4 <- NA
AFCst$SiOH4 <- NA
AFCst$PO4 <- NA
AFCst$Chla <- NA
AFCst$Phaeo_a <- NA

# Fill with for loop:
for(i in AFCst$id) {
	
	AFCst[which(AFCst$id == i),"Temp"] <- ddf[which(ddf$ID == i),"Ts"]
	AFCst[which(AFCst$id == i),"Sal"] <- ddf[which(ddf$ID == i),"Ss"]
	AFCst[which(AFCst$id == i),"pH"] <- ddf[which(ddf$ID == i),"pH"]
	AFCst[which(AFCst$id == i),"NO3"] <- ddf[which(ddf$ID == i),"NO3"]
	AFCst[which(AFCst$id == i),"NH4"] <- ddf[which(ddf$ID == i),"NH4"]
	AFCst[which(AFCst$id == i),"SiOH4"] <- ddf[which(ddf$ID == i),"SiOH4"]
	AFCst[which(AFCst$id == i),"PO4"] <- ddf[which(ddf$ID == i),"PO4"]
	AFCst[which(AFCst$id == i),"Phaeo_a"] <- ddf[which(ddf$ID == i),"Phaeo_a"]
	AFCst[which(AFCst$id == i),"Chla"] <- ddf[which(ddf$ID == i),"Chla"]
	
} # eo for loop
	
### Check linear models between CA coordinates and environment	
colnames(AFCst)
summary(lm(Ax1 ~ Temp, data = AFCst)) # Nope
summary(lm(Ax1 ~ Sal, data = AFCst)) # Nope
summary(lm(Ax1 ~ depth, data = AFCst)) # Nope
summary(lm(Ax1 ~ Chla, data = AFCst)) # Yes ; R-squared: 0.1493 ; p-value: 0.000198
summary(lm(Ax1 ~ NO3, data = AFCst)) # Nope
summary(lm(Ax1 ~ SiOH4, data = AFCst)) # Yes, even more than Chla R-squared:  0.1657, p-value: 8.718e-05

### OK. Chl-a seasonality is driven by SIlicates utilization by large colonial diatoms --> peaks in April/ May, decrease in July (though still high) and drop in Fall, when Nanophyto take over.
### Counts of Dinoflagellates make up the second CA axis --> probably due to interannual variability + inter-station variability (time to check the time series ;-))


##### Plot time series of T, S, Chla, SiOH4 etc. facet_grid per station

# Temperature
quartz()
ggplot(data = ddf) + 
	geom_point(aes(x = date, y = Ts, fill = factor(station)), data = ddf, pch = 21, colour = "black") + 
	geom_line(aes(x = date, y = Ts), data = ddf, colour = "black", linetype = "dashed") + 
	facet_wrap(~ station, nrow = 2, ncol = NULL) + xlab("Date") + ylab("Surface temperature (°C)") + 
	theme_linedraw()

# Chl-a
quartz()
ggplot(data = ddf) + 
		geom_point(aes(x = date, y = log(Chla), fill = factor(station)), data = ddf, pch = 21, colour = "black") + 
		geom_line(aes(x = date, y = log(Chla)), data = ddf, colour = "black", linetype = "dashed") + 
		facet_wrap(~ station, nrow = 2, ncol = NULL) + xlab("Date") + ylab("Surface chlorophyll a concentration\n(µg/L, logged)") + 
		theme_linedraw()

# total phytoplankton counts
quartz()
ggplot(data = ddf) + 
	geom_point(aes(x = date, y = log(n_phyto_total), fill = factor(station)), data = ddf, pch = 21, colour = "black") + 
	geom_line(aes(x = date, y = log(n_phyto_total)), data = ddf, colour = "black", linetype = "dashed") + 
	facet_wrap(~ station, nrow = 2, ncol = NULL) + xlab("Date") + ylab("Total number of counted\nphytoplankton cells (Log)") + 
	theme_linedraw()

### Link between Chla concentration and diatoms ?
summary(lm( log(Chla) ~ pourc_diatom, data = ddf))
quartz()
ggplot(data = ddf) + 
	geom_point(aes(x = log(Chla), y = log(n_diato)), data = ddf[which(ddf$Chla < 10),], pch = 21, colour = "black", fill = "grey50") + 
	xlab("Surface chlorophyll a concentration\n(µg/L, logged)") + ylab("Percentage of diatoms") + theme_linedraw()
# Not obvious...
	
# Percentage of diatoms
quartz()
ggplot(data = ddf[which(ddf$year %in% c(2011:2013)),]) + 
		geom_point(aes(x = date, y = pourc_diatom, fill = factor(station)), data = ddf, pch = 21, colour = "black") + 
		geom_line(aes(x = date, y = pourc_diatom), data = ddf, colour = "black", linetype = "dashed") + 
		facet_wrap(~ station, nrow = 2, ncol = NULL) + xlab("Date") + ylab("% of Diatoms") + 
		theme_linedraw()

# Total dinoflagellates counts
quartz()
ggplot(data = ddf[which(ddf$year %in% c(2011:2013)),]) + 
		geom_point(aes(x = date, y = log(n_dino), fill = factor(station)), data = ddf, pch = 21, colour = "black") + 
		geom_line(aes(x = date, y = log(n_dino)), data = ddf, colour = "black", linetype = "dashed") + 
		facet_wrap(~ station, nrow = 2, ncol = NULL) + xlab("Date") + ylab("Total number of counted\nDinoflagelattes cells (Log)") + 
		theme_linedraw()
				
				
### Now, explore the difference between surface and bottom temperatures (in absolute), in order to assess which stations belong to stratified environments.
ddf$Tdiff <- abs(ddf$Ts - ddf$Tf)
quartz()
ggplot(data = ddf) + 
	geom_point(aes(x = date, y = Tdiff, fill = factor(station)), data = ddf, pch = 21, colour = "black") + 
	geom_line(aes(x = date, y = Tdiff), data = ddf, colour = "black", linetype = "dashed") + 
	facet_wrap(~ station, nrow = 2, ncol = NULL) + xlab("Date") + ylab("∆Temperature (°C, surface - bottom)") + 
	theme_linedraw() + geom_hline(yintercept = 0, "black")
	
### Do the same but by separating the grid as a function of time ; for ransect b and d separately	
quartz()
ggplot(data = ddf[which(ddf$transect == "b" & ddf$year > 2010),]) + 
	geom_point(aes(x = x, y = Ts), pch = 21, colour = "black", fill = "#2b83ba") + 
	geom_line(aes(x = x, y = Ts), colour = "black", linetype = "dashed") + 
	facet_wrap(~ date, nrow = 1, ncol = NULL) + xlab("Longitude") + ylab("Surface temperature (°C)\nTransect B") + 
	theme_linedraw()
	
quartz()
ggplot(data = ddf[which(ddf$transect == "d" & ddf$year > 2010),]) + 
		geom_point(aes(x = x, y = Ts), pch = 21, colour = "black", fill = "#d7191c") + 
		geom_line(aes(x = x, y = Ts), colour = "black", linetype = "dashed") + 
		facet_wrap(~ date, nrow = 1, ncol = NULL) + xlab("Longitude") + ylab("Surface temperature (°C)\nTransect D") + 
		theme_linedraw()	

		
### Map environmental conditions for each pnmir transect
dates <- unique(ddf$date)
# Get rid of: 2010 and July 2012 and 3rd and 7th October 2013
dates <- dates[c(2:5,7:10,12,14:16)]
# Okay. Now, for each of the 14 dates, map the env conditions + phyto counts on a map and save them in a proper directory, with a for loop.
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/cartes_env_phyto_cruises/")
map_dir <- getwd()

# Also map difference between surface and bottom SST (stratif index...)
ddf$Tdiff <- ddf$Ts - ddf$Tf

# For testing:
# d <- dates[1]

for(d in as.character(dates)[10:length(dates)] ) {
		
		# Useless message
		message(paste("Doing ",d, sep = "")) 
		t <- ddf[which(ddf$date == d),]	
		
		# Create dir() and go to it
		dir.create(path = paste(map_dir,"/",d, sep = ""))
		setwd( paste(map_dir,"/",d,"/", sep = "") )
		
		### Maps
		# Ts
		mapT <- ggplot() + geom_point(aes(x = x, y = y, fill = Ts), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("SST (°C)", palette = "YlOrRd", direction = 1) + coast + theme_linedraw()
		# save	
		ggsave(plot = mapT, filename = "map_Ts.pdf", dpi = 300, width = 6, height = 6)
			
		# Tdiff	
		mapTdiff <- ggplot() + geom_point(aes(x = x, y = y, fill = Tdiff), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_gradient2("∆SST (°C)", low = "blue", high = "red") + coast + theme_linedraw()
		
		ggsave(plot = mapTdiff, filename = "map_Tdiff.pdf", dpi = 300, width = 6, height = 6)		
		
		# Ss
		mapS <- ggplot() + geom_point(aes(x = x, y = y, fill = Ts), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("SSS", palette = "YlOrRd", direction = 1) + coast + theme_linedraw()
		
		ggsave(plot = mapS, filename = "map_Ss.pdf", dpi = 300, width = 6, height = 6)		
		
		# NO3
		mapNO3 <- ggplot() + geom_point(aes(x = x, y = y, fill = NO3), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("Nitrates\nconcentration\n(µmol/L)", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		ggsave(plot = mapNO3, filename = "map_NO3.pdf", dpi = 300, width = 6, height = 6)		
		
		# PO4
		mapPO4 <- ggplot() + geom_point(aes(x = x, y = y, fill = PO4), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("Phosphates\nconcentration\n(µmol/L)", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		ggsave(plot = mapPO4, filename = "map_PO4.pdf", dpi = 300, width = 6, height = 6)		
		
		# SiOH4
		mapSi <- ggplot() + geom_point(aes(x = x, y = y, fill = SiOH4), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("Silicates\nconcentration\n(µmol/L)", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		ggsave(plot = mapSi, filename = "map_SiOH4.pdf", dpi = 300, width = 6, height = 6)		
		
		# Chla
		mapChla <- ggplot() + geom_point(aes(x = x, y = y, fill = Chla), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("Chlorophyll a\nconcentration\n(µg/L)", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		ggsave(plot = mapChla, filename = "map_Chla.pdf", dpi = 300, width = 6, height = 6)	
		
		# Phaeo_a
		mapPhaeo <- ggplot() + geom_point(aes(x = x, y = y, fill = Phaeo_a), data = t, pch = 21, colour = "black", size = 3) +
				scale_fill_distiller("Phaeophytin a\nconcentration\n(µg/L)", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		ggsave(plot = mapPhaeo, filename = "map_Phaeo.pdf", dpi = 300, width = 6, height = 6)			
			
		# n_phyto_total
		#mapphyto <- ggplot() + geom_point(aes(x = x, y = y, fill = n_phyto_total), data = t, pch = 21, colour = "black", size = 3) +
					#scale_fill_distiller("Total\nphytoplankton\ncells", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		#ggsave(plot = mapphyto, filename = "map_totalphyto.pdf", dpi = 300, width = 6, height = 6)			
		
		# n_diato
		#mapdiatom <- ggplot() + geom_point(aes(x = x, y = y, fill = n_diato), data = t, pch = 21, colour = "black", size = 3) +
					#scale_fill_distiller("Total\ndiatoms\ncells", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		#ggsave(plot = mapdiatom, filename = "map_totaldiatom.pdf", dpi = 300, width = 6, height = 6)			
		
		# n_dino
		#mapdino <- ggplot() + geom_point(aes(x = x, y = y, fill = n_dino), data = t, pch = 21, colour = "black", size = 3) +
					#scale_fill_distiller("Total\ndinoflagellates\ncells", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		#ggsave(plot = mapdino, filename = "map_totaldino.pdf", dpi = 300, width = 6, height = 6)			
		
		# n_nano
		#mapnano <- ggplot() + geom_point(aes(x = x, y = y, fill = n_nano), data = t, pch = 21, colour = "black", size = 3) +
					#scale_fill_distiller("Total\nnanoflagellates\ncells", palette = "YlGnBu", direction = 1) + coast + theme_linedraw()
		
		#ggsave(plot = mapnano, filename = "map_totalnano.pdf", dpi = 300, width = 6, height = 6)			
				
		# Go back to main directory
		gc()
		rm(t)
		setwd(map_dir)
	
} # eo for looping


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### 25/04/2017: Perform CA on species-level phytoplankton counts

### Load the data as .csv:
ddf <- read.csv("PNMIR_cruise_phyto_counts.csv", h = TRUE, dec = ",", sep = ";")
dim(ddf)
str(ddf)
colnames(ddf)
# Make some adjustments:
ddf$date <- paste(ddf$day, ddf$month, ddf$year, sep = "/")
ddf$date <- as.Date(ddf$date, format = "%d/%m/%Y") # should be okay
rownames(ddf) # Same as ID, ok

### Correct station coordibates plz...using 'coords_stations'
ddf[which(ddf$station == "b1"), "x"] <- -4.616667 ; ddf[which(ddf$station == "b1"), "y"] <- 48.31667
ddf[which(ddf$station == "b2"), "x"] <- -4.783333 ; ddf[which(ddf$station == "b2"), "y"] <- 48.31667
ddf[which(ddf$station == "b3"), "x"] <- -4.950000 ; ddf[which(ddf$station == "b3"), "y"] <- 48.31667
ddf[which(ddf$station == "b4"), "x"] <- -5.050000 ; ddf[which(ddf$station == "b4"), "y"] <- 48.33333
ddf[which(ddf$station == "b5"), "x"] <- -5.166667 ; ddf[which(ddf$station == "b5"), "y"] <- 48.33333
ddf[which(ddf$station == "b6"), "x"] <- -5.250000 ; ddf[which(ddf$station == "b6"), "y"] <- 48.33333
ddf[which(ddf$station == "b7"), "x"] <- -5.416667 ; ddf[which(ddf$station == "b7"), "y"] <- 48.33333
ddf[which(ddf$station == "d1"), "x"] <- -4.416667 ; ddf[which(ddf$station == "d1"), "y"] <- 48.16667
ddf[which(ddf$station == "d2"), "x"] <- -4.616667 ; ddf[which(ddf$station == "d2"), "y"] <- 48.16667
ddf[which(ddf$station == "d3"), "x"] <- -4.783333 ; ddf[which(ddf$station == "d3"), "y"] <- 48.16667
ddf[which(ddf$station == "d4"), "x"] <- -4.950000 ; ddf[which(ddf$station == "d4"), "y"] <- 48.16667
ddf[which(ddf$station == "d5"), "x"] <- -5.166667 ; ddf[which(ddf$station == "d5"), "y"] <- 48.16667
ddf[which(ddf$station == "d6"), "x"] <- -5.250000 ; ddf[which(ddf$station == "d6"), "y"] <- 48.16667

# Add stratif index (though poor)
ddf$Tdiff <- ddf$Ts - ddf$Tf

### re-perform CA on seleted species (non zero)
colnames(ddf)
res.ca <- CA(X = ddf[-c(75),c(25,63,97:99)] )  ### Get rid of station 75 because of dinoflagellate peak
summary(res.ca)
# nearly 100% of total variance along first two CA axes !
### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2])

# For objects
data_for_months <- ddf[-c(75), c("ID","year","month","max_depth")]
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					month = data_for_months$month, year = data_for_months$year, id = data_for_months$ID, depth = data_for_months$max_depth )

# Add seasons:
AFCst$season <- NA
AFCst[AFCst$month %in% c(4:6) ,"season"] <- "spring"
AFCst[AFCst$month == 7 ,"season"] <- "summer"
AFCst[AFCst$month == 10 ,"season"] <- "fall"

# quartz()
plot <- ggplot() +
  			geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  			geom_point(aes(x = Ax1, y = Ax2, shape = factor(year), fill = factor(season) ), data = AFCst, alpha = 0.7, colour = "black", size = 3) +
			#scale_shape_manual("Localisation", values = c(21,22)) + 
			scale_fill_manual("Season", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
			scale_shape_manual("Year", values = c(21,22,23,24)) + 
			scale_x_continuous(paste("CA 1 (", round(res.ca$eig[,2][1], 2),"%)", sep = ""), limits = c(-1.2,1.5) ) + 
			scale_y_continuous(paste("CA 2 (", round(res.ca$eig[,2][2], 2),"%)", sep = "")) + 
			#geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  			#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) + 
			theme_bw()

ggsave(plot = plot, filename = "CA_plot_large_groups.pdf", dpi = 300, width = 7, heigh = 5)

			
### Check correlation between the stations' CA scores and env variables
AFCst$Temp <- NA
AFCst$Sal <- NA
AFCst$pH <- NA
AFCst$NO3 <- NA
AFCst$NH4 <- NA
AFCst$SiOH4 <- NA
AFCst$PO4 <- NA
AFCst$Chla <- NA
AFCst$Phaeo_a <- NA
AFCst$Tdiff

# Fill with for loop:
for(i in AFCst$id) {
	
		AFCst[which(AFCst$id == i),"Temp"] <- ddf[which(ddf$ID == i),"Ts"]
		AFCst[which(AFCst$id == i),"Sal"] <- ddf[which(ddf$ID == i),"Ss"]
		AFCst[which(AFCst$id == i),"pH"] <- ddf[which(ddf$ID == i),"pH"]
		AFCst[which(AFCst$id == i),"NO3"] <- ddf[which(ddf$ID == i),"NO3"]
		AFCst[which(AFCst$id == i),"NH4"] <- ddf[which(ddf$ID == i),"NH4"]
		AFCst[which(AFCst$id == i),"SiOH4"] <- ddf[which(ddf$ID == i),"SiOH4"]
		AFCst[which(AFCst$id == i),"PO4"] <- ddf[which(ddf$ID == i),"PO4"]
		AFCst[which(AFCst$id == i),"Phaeo_a"] <- ddf[which(ddf$ID == i),"Phaeo_a"]
		AFCst[which(AFCst$id == i),"Chla"] <- ddf[which(ddf$ID == i),"Chla"]
		AFCst[which(AFCst$id == i),"Tdiff"] <- ddf[which(ddf$ID == i),"Tdiff"]
	
} # eo for loop
	
### Check linear models between CA coordinates and environment	
colnames(AFCst)
summary(lm(Ax1 ~ Temp, data = AFCst)) # Very weakly significant...higher temperatures seem to negatively affect Guinardia sp.
summary(lm(Ax1 ~ Sal, data = AFCst)) # Most significant correlation
summary(lm(Ax1 ~ depth, data = AFCst)) # 
summary(lm(Ax1 ~ Chla, data = AFCst)) # 
summary(lm(Ax1 ~ NO3, data = AFCst)) # 
summary(lm(Ax1 ~ SiOH4, data = AFCst)) # 
summary(lm(Ax1 ~ Tdiff, data = AFCst)) # 
 
summary(lm(Ax2 ~ Temp, data = AFCst)) # yes, dinoflagelattes more abundant at warmer temp
summary(lm(Ax2 ~ Sal, data = AFCst)) # nope
summary(lm(Ax2 ~ depth, data = AFCst)) # nope
summary(lm(Ax2 ~ Chla, data = AFCst)) # yes, dinoflagellates more abundant in less productive conditions
summary(lm(Ax2 ~ NO3, data = AFCst)) #  nope
summary(lm(Ax2 ~ SiOH4, data = AFCst)) # very signif --> more diatoms in silicate depleted conditions (logical)
summary(lm(Ax2 ~ Tdiff, data = AFCst)) # nope

### Nice correlogram:
library("corrgram")
cormat <- round( cor(na.omit(AFCst[,c(1,6,7,10,12,13,14,15,16)]), method = "spearman"), 2)
head(cormat)

# Necessary FUNs
get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
quartz()
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
	 scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
  	 theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()	  
reorder_cormat <- function(cormat) {
	# Utiliser la corrélation entre les variables
	# comme mésure de distance
	dd <- as.dist((1-cormat) / 2)
	hc <- hclust(dd)
	cormat <-cormat[hc$order, hc$order]
}	
# Reordonner la matrice de corrélation
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
			theme_minimal() + # minimal theme
			theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

### Okray, let's add some text
heatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
			theme( axis.title.x = element_blank(),
  		  			axis.title.y = element_blank(),
  				  	panel.grid.major = element_blank(),
  				  	panel.border = element_blank(),
 				   	panel.background = element_blank(),
  				  	axis.ticks = element_blank(),
  				  	legend.justification = c(1, 0),
  				  	legend.position = c(0.6, 0.7),
  				  	legend.direction = "horizontal") +
  			guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
quartz()
heatmap		




### Perform correspondence analysis on species-level counts 
# Colsums to find which species have enough presences
# colnames(ddf)
sums <- colSums(as.matrix(ddf[,c(26:62, 64:93, 97:99)]))
# Keep the names which present more than X obs
names <- names( sums[which( sums > 100000)] )
names
# Perform CA
res.ca <- CA(X = ddf[,c(names)] )
summary(res.ca)
# 53.76% of total variance along the four first axes when names > 1000
# 53.93% of total variance along the four first axes when names > 5000
# 55.50% of total variance along the four first axes when names > 10000
# 55.83% of total variance along the four first axes when names > 20000
# 58.90% of total variance along the four first axes when names > 50000
# 60.81% of total variance along the four first axes when names > 100000

### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
					Ax4 = res.ca$col$coord[,4],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2], 
					Cos2_2 = res.ca$col$cos2[,3] + res.ca$col$cos2[,4])

# For objects
data_for_months <- ddf[,c("ID","month","max_depth",names)]
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
					Ax4 = res.ca$row$coord[,4],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					Cos2_2 = res.ca$row$cos2[,3] + res.ca$row$cos2[,4],
					month = data_for_months$month, id = data_for_months$ID, depth = data_for_months$max_depth )

head(AFCst)

quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Localisation", values = c(21,22)) + 
	scale_fill_manual("Month", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()


quartz()
ggplot() +
	geom_point(aes(x = Ax3, y = Ax4), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
	geom_point(aes(x = Ax3, y = Ax4, size = Cos2_2, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Localisation", values = c(21,22)) + 
	scale_fill_manual("Month", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 3 (",round(res.ca$eig$per[3], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 4 (",round(res.ca$eig$per[4], 2),"%)", sep = "")) + 
	geom_text(aes(x = Ax3, y = Ax4+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()


### And without Guinardia sp. and Skeletonema ? 
sums <- colSums(as.matrix(ddf[,c(26:62, 64:93, 97:99)]))
# Keep the names which present more than X obs
names <- names( sums[which( sums > 10000)] )
names
# Delete Guinardia sp. and Skeletonema sp.
names <- names[-c(6,14)]
names
# Perform CA
res.ca <- CA(X = ddf[,c(names)] )
summary(res.ca)
# 65.69% of total variance along the four first axes when names > 100000
# 62.79% of total variance along the four first axes when names > 50000
# 57.87% of total variance along the four first axes when names > 10000

### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
					Ax4 = res.ca$col$coord[,4],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2], 
					Cos2_2 = res.ca$col$cos2[,3] + res.ca$col$cos2[,4])

# For objects
data_for_months <- ddf[,c("ID","year","month","max_depth",names)]
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
					Ax4 = res.ca$row$coord[,4],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					Cos2_2 = res.ca$row$cos2[,3] + res.ca$row$cos2[,4],
					month = data_for_months$month, year = data_for_months$year, id = data_for_months$ID, depth = data_for_months$max_depth )

head(AFCst)

quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month), shape = factor(year) ), data = AFCst, alpha = 0.7, colour = "black" ) +
	#scale_shape_manual("Localisation", values = c(21,22)) + 
	scale_fill_manual("Month", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_shape_manual("Year", values = c(21, 22, 23, 24)) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()


quartz()
ggplot() +
	geom_point(aes(x = Ax3, y = Ax4), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
	geom_point(aes(x = Ax3, y = Ax4, size = Cos2_2, fill = factor(month), shape = factor(year) ), data = AFCst, alpha = 0.7, colour = "black" ) +
	#scale_shape_manual("Localisation", values = c(21,22)) + 
	scale_fill_manual("Month", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_shape_manual("Year", values = c(21, 22, 23, 24)) + 
	scale_x_continuous(paste("CA 3 (",round(res.ca$eig$per[3], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 4 (",round(res.ca$eig$per[4], 2),"%)", sep = "")) + 
	geom_text(aes(x = Ax3, y = Ax4+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()

	
###	Differences appear in the phytoplanktpon community between July and May.
# 1) Main differences appear between diatom-dominated communities (CA 1 when large groups only, 90% of variance) and nanoflagellate dominated comminities.
#    Diatoms show dominance in May (whatever the year) and in July 2011 (but not in July 2013 which was dominated by Nanoflag.).
#    Nanoflagellates mostly comprise: Cryptophyta, Pyramimonas (Prasinophyta) and Coccolithophores (Primnesiophyta, Emiliana-like).

# 2) The diatoms that dominated the community in July 2011 are different than the ones dominating in May (Ca on species-level counts).
#    Indeed, in May 2012, diatoms are mostly composed opf Guinardia sp. (but also true fpr 2 stations in Oct 2013). When putting Guinardia aside, we can see that
#    May 2011 and 2013 are dominated by the following diatoms: Chaetoceros sp., Cerataulina pelagica and Dactyliosolen fragilissimus. All of these are large centric colonial diatoms.
#    Meanwhile, the diatoms that dominated in July 2011 are Pseudonitzschia and Leptocylindrus sp. Together with some dinoflagellates: Lepidodinium and larger Gymnodinium. 
#    In both July 2011 and May 2011-2013, Rhizosolenia sp., larger Thalassiosira sp., Amphidinium crassum and Gyrodinium sp. are also present, to a lesser extent.

# 3) What about Octobers and July 2013 ? They are dominated by larger counts of Alexandrium pseudogoniaulax (toxic bloom-forming species), Cryptophyta, Prorocentrum sp., and to a lesser extent:
#    Primnesiophyta (Emiliana-like), Prasinophta, smaller Gymnodinium, Heterocapsa sp. (toxic dinoflagellate) and Scrippsiella sp. (dinoflagellate).
#    October 2013 presents higher counts of A. pseudogoniaulax.


##### Now, how can we analyze the environmental data with the phytoplankton counts together in an ordination technique ?
##### --> Canonical Correspondence Analysis (CCA) ! 

### CCA is the canonical form of CA. Any datatable that could be subject to a CA is suitable to be a response matrix Y for CCA, especially species P/A or abundance data tables.
### CCA combines multiple regression and CA (instead of PCA for RDA) in 2 steps : 
#   1.  First step: weighted multiple regression
#   2.  Second step: correspondence analysis (CA) on fitted values

# Note: Instead of using Y for the regression step, CCA uses matrix Q, which represents the contributions to chi-2 statistics (as in CA)
# Besides, for the regression the matrix X is standardized using as weight the diagonal matrix of row totals D(fi+)
require("vegan")
?cca

# Perform
cctable <- na.omit(ddf[,c(1:9,12,17,19,20,23,24,26:62,64:93,97:99,100)])
dim(cctable)
cctable
colnames(cctable)

# Like above, only the use the species that present at least 10000 counts or 100000
sums <- colSums(as.matrix(cctable[,c(16:85)]))
# Keep the names which present more than X obs
names <- names( sums[which( sums > 50000)] )
names
#
cca <- vegan::cca( X = as.matrix(cctable[,names]), Y = as.matrix(cctable[,c(9:15,86)]) )
summary(cca, scaling = 1, axes = 3) 

# test significance: 
anova(cca, step=1000, perm.max=1000)
#Model: cca(X = as.matrix(cctable[, names]), Y = as.matrix(cctable[, c(9:15, 86)]))
#         Df ChiSquare      F Pr(>F)    
#Model     8   0.94104 5.1299  0.001 ***
#Residual 98   2.24719 

# Quite significant ! 
# Accumulated constrained eigenvalues
#Importance of components:
#                        CCA1   CCA2   CCA3   CCA4    CCA5    CCA6     CCA7
#Eigenvalue            0.3529 0.3104 0.1260 0.0879 0.03135 0.01954 0.009747
#Proportion Explained  0.3751 0.3298 0.1339 0.0934 0.03332 0.02076 0.010360
#Cumulative Proportion 0.3751 0.7049 0.8388 0.9322 0.96556 0.98632 0.996680

# CCA triplot
# Scaling 1
plot(cca, display=c("sp","lc","cn"), scaling = 1, main = "Triplot CCA fish ~ env2 - scaling 1") 

# Scaling 2 (default scaling)
plot(cca, display=c("sp","lc","cn"), main = "Triplot CCA fish ~ env2 - scaling 2") 

# The scaling 1 triplot focuses on the distance relationships among sites, but the presence of species with extreme scores renders the plot difficult to interpret beyond trivialities. 
# Therefore, it may be useful to redraw it without the species

# CCA scaling 1 biplot without species (using lc site scores)
plot(cca, display=c("lc","cn"), scaling=1, main="Biplot CCA fish ~ env2 - scaling 1") 
# CCA scaling 2 biplot without sites 
quartz()
plot(cca, display=c("sp","cn"), main="Biplot CCA fish ~ env2 - scaling 2") 

### Nicer plots please 
str(cca$CCA$biplot)

species_scores <- data.frame(species = rownames(cca$CCA$v), CCA1 = cca$CCA$v[,"CCA1"], CCA2 = cca$CCA$v[,"CCA2"], CCA3 = cca$CCA$v[,"CCA3"])

stations_scores <- data.frame(CCA1 = cca$CCA$u[,"CCA1"], CCA2 = cca$CCA$u[,"CCA2"], CCA3 = cca$CCA$u[,"CCA3"])
stations_scores$date <- cctable$date
stations_scores$month <- cctable$month
stations_scores$year <- cctable$year
stations_scores$station <- cctable$station
stations_scores$season <- NA
stations_scores[stations_scores$month %in% c(4:6) ,"season"] <- "spring"
stations_scores[stations_scores$month == 7 ,"season"] <- "summer"
stations_scores[stations_scores$month == 10 ,"season"] <- "fall"

vars_scores <- data.frame(var = rownames(cca$CCA$biplot), CCA1 = cca$CCA$biplot[,"CCA1"], CCA2 = cca$CCA$biplot[,"CCA2"], CCA3 = cca$CCA$biplot[,"CCA3"] )


### Here we go...
quartz()
ccaplot <- ggplot() + geom_point(aes(x = CCA1*1.5, y = CCA2*1.5), data = species_scores, pch = 21, fill = "grey70", colour = "black", size = 4) + 
		   	geom_point(aes(x = CCA1, y = CCA2, fill = factor(season)), data = stations_scores, pch = 23, colour = "black", size = 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = CCA1*3, yend = CCA2*3), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour = "black", data = vars_scores) + 
		   	scale_fill_manual("Season", values = c("#fdae61","#abdda4","#d7191c"), labels = c("Spring","Summer","Fall")) + 
		   	#geom_text(aes(x = CCA1, y = CCA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
			xlab("CCA 1 (37.5%)") + ylab("CCA 2 (33.0%)") + scale_x_continuous(limits = c(-3,4)) + 
			scale_y_continuous(limits = c(-4,3)) + theme_bw()

ggsave(filename = "CCA1_nonames_14_06_18.pdf", plot = ccaplot, dpi = 300, height = 6, width = 9)	
	
	
### Great. Re-perform without Guinardia sp. which has a LOT of weight on CCA1
names <- names[-c(4)]
cca <- vegan::cca( X = as.matrix(cctable[,names]), Y = as.matrix(cctable[,c(9:15,86)]) )
summary(cca, scaling = 1, axes = 3) 
# test significance: 
anova(cca, step=1000, perm.max=1000)

# CCA scaling 1 biplot without species (using lc site scores)
plot(cca, display=c("lc","cn"), scaling=1, main="Biplot CCA fish ~ env2 - scaling 1") 
# CCA scaling 2 biplot without sites 
plot(cca, display=c("sp","cn"), main="Biplot CCA fish ~ env2 - scaling 2") 

### Nicer plots please ?
str(cca$CCA$biplot)

species_scores <- data.frame(species = rownames(cca$CCA$v), CCA1 = cca$CCA$v[,"CCA1"], CCA2 = cca$CCA$v[,"CCA2"], CCA3 = cca$CCA$v[,"CCA3"])

stations_scores <- data.frame(CCA1 = cca$CCA$u[,"CCA1"], CCA2 = cca$CCA$u[,"CCA2"], CCA3 = cca$CCA$u[,"CCA3"])
stations_scores$date <- cctable$date
stations_scores$month <- cctable$month
stations_scores$year <- cctable$year
stations_scores$station <- cctable$station
stations_scores$season <- NA
stations_scores[stations_scores$month %in% c(4:6) ,"season"] <- "spring"
stations_scores[stations_scores$month == 7 ,"season"] <- "summer"
stations_scores[stations_scores$month == 10 ,"season"] <- "fall"

vars_scores <- data.frame(var = rownames(cca$CCA$biplot), CCA1 = cca$CCA$biplot[,"CCA1"], CCA2 = cca$CCA$biplot[,"CCA2"], CCA3 = cca$CCA$biplot[,"CCA3"] )


### Here we go...
#quartz()
ccaplot <- ggplot() + geom_point(aes(x = CCA1, y = CCA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size = 4) + 
		   	geom_point(aes(x = CCA1, y = CCA2, fill = factor(month)), data = stations_scores, pch = 23, colour = "black", size = 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = CCA1*3, yend = CCA2*3), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour = "black", data = vars_scores) + 
		   	scale_fill_manual("Month", values = c("#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		   	#geom_text(aes(x = CCA1, y = CCA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
		    xlab("CCA 1 (55.1%)") + ylab("CCA 2 (28.8%)") + scale_x_continuous(limits = c(-3, 4.2)) + theme_bw()


ggsave(filename = "CCA2_nonames.pdf", plot = ccaplot, dpi = 300, height = 8, width = 12)	







