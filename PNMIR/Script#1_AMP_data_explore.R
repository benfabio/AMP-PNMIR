
##### 04/04/2017: R Script to examine the content and structure of the ZooScan data from the AMP PMMI project © Fabio Benedetti, OOV-UMS, AFB
##### Aims to:

#	- Load all .tsv data from the 'Tests_04_04_17' directory ; rbind
#	- Compute the % of validated vignettes for each sampling point: 4 months per 5 years --> 20 maps (all stations from plan A should be 100% validated)
#	- Map said % on a map with a nice coastline
#	- Compute a percentage of multiple objects within living objetcs, per sample_id (for Amanda)
#	- Map the position of all the stations (considering Platresse)
#	- Learn how to compute biovolumes, and compute the % of biovolumes that is due to multiple living objects 
#	- Explore Sardines data

### Latest update: 09/06/2017

library("dplyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("Hmisc")

# MATLAB jet color gradient
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Load coastline:
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data")
cl <- read.csv("gshhg_h_2017-04-04_17-50-04.csv")
coast <- list(
 	 	# the coast polygon itself, a bit lighter than usual to avoid taking too much attention out of the data itself
  	  	geom_polygon(aes(x=lon, y=lat), data= cl, fill = "grey30"),
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

quartz()
ggplot() + coast + theme_linedraw()
# OKAY.

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Load ZooScan data:
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/Data_08_06_17/")
proj_dir <- getwd()
projects <- dir()
# p <- "export_435_20170608_1546"
# p <- "export_377_20170921_1131"
projs <- lapply(projects, function(p) {
	
			# Go to project dir
			message(paste("Doing project #", p, sep = ""))
			setwd(paste(proj_dir,"/",p,"/", sep = ""))
			all_samples <- dir()
			
			# Read and rbind the samples from each project
			samples <- lapply(all_samples, function(s) {
					
					# Load and return
					ss <- read.table(file = s, sep = '\t', header = T)
					return(ss)
					
			}) # eo 2nd lapply 
			
			all_samples <- do.call(rbind, samples)
			# dim(all_samples) ; colnames(all_samples); summary(all_samples)
			rm(samples)
			
			# Return, by restricting to columns of interest
			return(all_samples[,c("object_lon","object_lat","object_date","sample_stationid","object_annotation_status","object_annotation_status","object_id","sample_id")])
}
) # eo lapply

all_projs_data <- do.call(rbind, projs)

str(all_projs_data)
dim(all_projs_data)
colnames(all_projs_data)
summary(all_projs_data)

summary(all_projs_data$object_annotation_status)
# predicted validated ### FOR THE 21-04-2017
# 137154     521170
#  (21%)       (79%)

### For the 08/06/17: 
# predicted validated 
#   130969    539540 
(539540 / (130969+539540))*100

### 04/04/17: Restrict to vertical tows
all_projs_data <- all_projs_data[grep(pattern='v', all_projs_data$object_id),]
dim(all_projs_data)
str(all_projs_data)
summary(all_projs_data)
### 20% non validated

### Now, per year and sampling station, compute % of predicted vs validated data and map !
# But first, find the year ^^
unique(all_projs_data$object_date) # 116
unique(all_projs_data$sample_stationid) # Name of the station: 15 names (b1-b7; d1-d6; douarnenez; molene)
str(all_projs_data$object_date) # It's an integer, need to change to factor or character...

library("lubridate") # Test lubridate package
all_projs_data$object_date <- as.character(all_projs_data$object_date)
head(year(ymd(all_projs_data$object_date)))
head(month(ymd(all_projs_data$object_date)))
head(day(ymd(all_projs_data$object_date)))
 
all_projs_data$year <- year(ymd(all_projs_data$object_date))
all_projs_data$month <- month(ymd(all_projs_data$object_date))
all_projs_data$day <- day(ymd(all_projs_data$object_date))
 
# unique(all_projs_data[which(all_projs_data$sample_stationid == "b3"), "month_object"])
# unique(all_projs_data[which(all_projs_data$sample_stationid == "b5"), "month_object"])
# unique(all_projs_data[which(all_projs_data$sample_stationid == "d3"), "month_object"])
# unique(all_projs_data[which(all_projs_data$sample_stationid == "d6"), "month_object"])
# unique(all_projs_data[which(all_projs_data$sample_stationid == "d2"), "month_object"])

# Fo tryin'
# y <- 2011
# m <- 10
# s <- "b3"

##### !!!!! 08/06/2017: because of issues between EcoTaxa and Zooprocess, some scans have corrupted information !!!!!!!!
# Get rid of them before going any further (using sample_id): 
# wp2_b1_20130603_v_tot_1
# wp2_b1_20131003_v_tot_1
# wp2_b2_20130603_v_tot_1
# wp2_b3_20130603_v_tot_1
# wp2_b4_20130604_v_tot_1
# wp2_b4_20131001_v_tot_1
# wp2_b5_20130604_v_tot_1
# wp2_b5_20131001_v_tot_1
# wp2_b6_20130604_v_tot_1
# wp2_b6_20130702_v_tot_1
# wp2_b6_20131001_v_tot_1
# wp2_b7_20130604_v_tot_1
# wp2_D1_20130515_v_tot_1
# wp2_D2_20130515_v_tot_1
# wp2_D2_20131007_v_tot_1
# wp2_d4_20130515_v_tot_1
# wp2_d4_20130704_v_tot_1
# wp2_d4_20131007_v_tot_1
# wp2_d6_20130704_v_tot_1

#corrupt <- c("wp2_b1_20130603_v_tot_1","wp2_b1_20131003_v_tot_1","wp2_b2_20130603_v_tot_1","wp2_b3_20130603_v_tot_1","wp2_b4_20130604_v_tot_1","wp2_b4_20131001_v_tot_1",
#			"wp2_b5_20130604_v_tot_1","wp2_b5_20131001_v_tot_1","wp2_b6_20130604_v_tot_1","wp2_b6_20130702_v_tot_1","wp2_b6_20131001_v_tot_1","wp2_b7_20130604_v_tot_1",
#			"wp2_D1_20130515_v_tot_1","wp2_D2_20130515_v_tot_1","wp2_D2_20131007_v_tot_1","wp2_d4_20130515_v_tot_1","wp2_d4_20130704_v_tot_1", "wp2_d4_20131007_v_tot_1","wp2_d6_20130704_v_tot_1")

head(all_projs_data)
### Add a column which will contain the sample_id + "_tot_1"
#all_projs_data$id <- paste(all_projs_data$sample_id, "_tot_1", sep = "")
### Get rid of the samples id that cause problems:
#dim(all_projs_data)
#all_projs_data <- all_projs_data[- which(all_projs_data$id %in% corrupt ),]
dim(all_projs_data)

years <- unique(all_projs_data$year)

for(y in years) {
	
		# Useless message
		paste("Doing ",y, sep = "")
	
		# Restrict
		yy <- all_projs_data[which(all_projs_data$year == y),]
	
		for(m in unique(yy$month) ) {
			
			# Useless message again
			paste("Doing ",m, sep = "")
			
			# Restrict again
			mm <- yy[which(yy$month == m),]
			
			# Watchout, sometimes, 'd1' and 'douarnenez' are separated...
			if( "douarnenez" %in% unique(mm$sample_stationid) ) {
					mm[which(mm$sample_stationid == "douarnenez"),"sample_stationid"] <- "d1"	
			} #
			
			# Compute % of predicted or validated data per station
			res <- lapply(unique(mm$sample_stationid), function(s) {
				
					station <- mm[which(mm$sample_stationid == s),]
					message(paste("Doing station ",s, sep = ""))
					
					# station$object_annotation_status
					per <- ( nrow(station[station$object_annotation_status == "validated",])/nrow(station) )*100
					
					# Return
					table <- data.frame(station = s, x = unique(station$object_lon), y = unique(station$object_lat), year = y, month = m, percentage_validated = per)
					return(table)
					
			}) # eo lapply
			
			tables <- do.call(rbind, res)
			#rm(res)
			# dim(tables)
			map <- ggplot() + geom_point(aes(x = x, y = y, fill = percentage_validated), data = tables, pch = 21, colour = "black", size = 3) + 
						scale_fill_distiller(name = "Percentage of\nvalidated vignettes", palette = "RdYlGn", direction = 1, limits = c(0,100)) + 
						#scale_size_manual(name = "Percentage of\nvalidated vignettes", limits = c(0,100)) + 
						geom_text(data = tables, aes(x = x, y = y - 0.02, label = as.character(station), hjust = 0), size = 2) + 
						ggtitle(paste("Pourcentage de vignettes validées\npour les stations échantillonnées en\n",m," ",y, sep = "")) + 
						coast + coord_quickmap() + theme_linedraw()
					
			# quartz()
			# map	
			setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
			ggsave(plot = map, filename = paste("Map_perc_validated_",m,"_",y,".pdf", sep = ""), dpi = 300, width = 6, height = 6)
			
			# rm(map, tables, res)
			setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/Data_08_06_17/")

		} # eo month for loop	
	
} # eo year for loop



##### 21/04/2017: Tell Amanda whether all samples from Plan A have been validated
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
samples_id <- read.csv("Plan_A.csv", h = TRUE, sep = ";")
### For each element of samples_id, compute % of validated 
# in 'all_projs_data', these values are stored in all_projs_data$object_id. But need to some manipulation first:
ids <- data.frame( do.call(rbind, strsplit(x = as.character(all_projs_data$object_id), split = "_")) )
all_projs_data$ID <- paste(ids$X1, ids$X2, ids$X3, ids$X4, ids$X5, ids$X6, sep = "_")

# match(all_projs_data$ID, samples_id$project_id) ### OKAY

# Create an empty column
samples_id$percentage_validated <- NA

# Fill with for loop or lapply...
# s <- samples_id$project_id[1] ### for testing

for(s in samples_id$project_id) {
		
		# Useless message
		message( paste("Doing ", s, sep = "") )
		d <- all_projs_data[which(all_projs_data$ID == s),]
		# compute percentage of validated vignettes
		# station$object_annotation_status
		per <- ( nrow(d[d$object_annotation_status == "validated",])/nrow(d) )*100
		# return
		samples_id[which(samples_id$project_id == s), "percentage_validated"] <- per
		
		gc(); rm(d, per)
	
} # eo for loop

head(samples_id$percentage_validated)
samples_id$percentage_validated
# Worked ? Yes. Plan A is complete.

##### Now, for each ID belonging to Plan B (so all not comproised in Plan A), compute number of vignettes to be validated so we can have an estimate of the quantity (hours) of work remaining
PlanB <- all_projs_data[which( !(all_projs_data$ID %in% samples_id$project_id)),]
years <- unique(PlanB$year_object)

for(y in years) {
	
		# Useless message
		paste("Doing ",y, sep = "")
	
		# Restrict
		yy <- PlanB[which(PlanB$year == y),]
	
		for(m in unique(yy$month_object) ) {
			
			# Useless message again
			paste("Doing ",m, sep = "")
			
			# Restrict again
			mm <- yy[which(yy$month == m),]
			
			# Watchout, sometimes, 'd1' and 'douarnenez' are separated...
			if( "douarnenez" %in% unique(mm$sample_stationid) ) {
					mm[which(mm$sample_stationid == "douarnenez"),"sample_stationid"] <- "d1"	
			} #
			
			# Compute % of predicted or validated data per station
			res <- lapply(unique(mm$sample_stationid), function(s) {
				
					station <- mm[which(mm$sample_stationid == s),]
					message(paste("Doing station ",s, sep = ""))
					
					# station$object_annotation_status
					n <- nrow( station[station$object_annotation_status != "validated",] )
					
					# Return
					table <- data.frame(station = s, x = unique(station$object_lon), y = unique(station$object_lat), year = y, month = m, n_objects_to_validate = n)
					return(table)
					
			}) # eo lapply
			
			tables <- do.call(rbind, res)
			#rm(res)
			# dim(tables)
			map <- ggplot() + geom_point(aes(x = x, y = y, fill = n_objects_to_validate), data = tables, pch = 21, colour = "black", size = 3) + 
						scale_fill_distiller(name = "Number of vignettes\nto validate", palette = "Spectral", direction = -1 ) + 
						#scale_size_manual(name = "Percentage of\nvalidated vignettes", limits = c(0,100)) + 
						geom_text(data = tables, aes(x = x, y = y - 0.02, label = as.character(station), hjust = 0), size = 2) + 
						ggtitle(paste("Nombre de vignettes à valider\npour les stations échantillonnées en\n",m," ",y, sep = "")) + 
						coast + coord_quickmap() + theme_linedraw()
					
			# quartz()
			# map	
			setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
			ggsave(plot = map, filename = paste("Map_n_tovalidate_",m,"_",y,".pdf", sep = ""), dpi = 300, width = 6, height = 6)
			
			# rm(map, tables, res)
			setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/Tests_21_04_17/")

		} # eo month for loop
		
	
} # eo year for loop


##### Save the information in a .csv file
to_validate <- data.frame(id = unique(PlanB$ID), n_to_validate = NA)
# fill with for loop
for(s in to_validate$id) {
		
		# Useless message
		message( paste("Doing ", s, sep = "") )
		d <- PlanB[which(PlanB$ID == s),]
		# compute percentage of validated vignettes
		# station$object_annotation_status
		n <- nrow( d[d$object_annotation_status != "validated",] )
		# return
		to_validate[which(to_validate$id == s), "n_to_validate"] <- n
		
		gc(); rm(d, n)
	
} # eo for loop

dim(to_validate)
head(to_validate)
sum(to_validate$n_to_validate)
# 162 496
sum( to_validate[c(1:40,116:159), "n_to_validate"] )
# 66 302

write.csv(to_validate, "PlanB.csv", sep = ";")


##### 06/04/2017: For Amanda Elineau, compute a percentage of multiple objects within living objetcs, per sample_id.

# For EACH acq_id, use 'all_projs_data' object to: 
#  - compute number of living objects ('object_annotation_hierarchy'), and the percentage of multiple objects within the living (from 'object_annotation_category')
#  - return the identity of of the samples with % of multiples > 10% and 15%

# summary(all_projs_data$object_annotation_hierarchy)
# class(all_projs_data$object_annotation_hierarchy) ; unique(all_projs_data$object_annotation_hierarchy)

# unique(all_projs_data$object_annotation_category)

samples <- unique(all_projs_data$acq_id)
length(samples)
# s <- samples[100]

res <- lapply(samples, function(s) {
	
				# Useless message
				message(paste("Doing ",s, sep = ""))
				ss <- all_projs_data[which(all_projs_data$acq_id == s),]
				
				# Number of living objects
				# strsplit() first
				split_status <- data.frame( do.call(rbind, strsplit(x = as.character(ss$object_annotation_hierarchy), split = ">")) )
				# head(split_status) ; tail(ss$object_annotation_hierarchy)
				colnames(split_status)[1] <- c("status")
				nlivings <- nrow( split_status[which(split_status$status == "living"),] )
				
				# And now, the number of multiples of multiples ('object_annotation_category')
				nmultiple <- nrow( ss[which(ss$object_annotation_category == "multiple"),] )
				
				# Compute percentage of multiples within libings (knowing that multiples are ALWAYS living objects)
				per <- (nmultiple/nlivings)*100
				
				# If % > 15%, then return (with %)
				#if( per >= 15 ) {
						
						d <- data.frame(sample = s, percentage_multiples = per)
						return(d)
						gc()
						rm(per, nmultiple, nlivings)
						
				#} else {
					
						#gc()
						#rm(per, nmultiple, nlivings)
					
				#} # eo if else loop
				
}) # eo lapply

percentages <- do.call(rbind, res)
rm(res)
dim(percentages)
summary(percentages)
percentages[which(percentages$percentage_multiples > 10),]
### No sample with % > 10%, max is like 8.89%
rm(percentages)


##### More interesting now: the % of living biovolume that is due to multiple objects ;-)
samples <- unique(all_projs_data$acq_id)
length(samples)
# s <- samples[10]

res <- lapply(samples, function(s) {
	
				# Useless message
				message(paste("Doing ",s, sep = ""))
				ss <- all_projs_data[which(all_projs_data$acq_id == s),]
				
				# Number of living objects
				# strsplit() first
				split_status <- data.frame( do.call(rbind, strsplit(x = as.character(ss$object_annotation_hierarchy), split = ">")) )
				# head(split_status) ; tail(ss$object_annotation_hierarchy)
				colnames(split_status)[1] <- c("status")
				nlivings <- nrow( split_status[which(split_status$status == "living"),] )
				ss$status <- split_status$status
				# And now, the number of multiples of multiples ('object_annotation_category')
				nmultiple <- nrow( ss[which(ss$object_annotation_category == "multiple"),] )
				# Compute percentage of multiples within libings (knowing that multiples are ALWAYS living objects)
				per <- (nmultiple/nlivings)*100
				
				### Now, compute biovolume ! (formulae in Vandromme et al. 2012)
				# First, need to compute the area from the 'minor' and the 'major' axes: 
				# summary(ss$object_major) ; summary(ss$object_minor) 
				areas <- (ss$object_major/2)*(ss$object_minor/2)*pi
				# gut gut
				# Now, compute spherical and elliptical biovolumes
				spher_BV <- (4/3)*pi*((areas/pi)^(3/2))
				ellipt_BV <- ( (4/3)*pi*(ss$object_major/2)*((ss$object_minor/2)^2) )
				
				### Supply to data:
				ss$areas <- areas
				ss$spher_BV <- spher_BV
				ss$ellipt_BV <- ellipt_BV
				
				# summary(lm(ellipt_BV ~ spher_BV, data = ss)) #  R-squared: 0.9807
				
				### Compute spher_BV for livings and then for multiples
				ellipt_BV_livings <- sum(ss[which(ss$status == "living"),"ellipt_BV"])
				ellipt_BV_multiple <- sum(ss[which(ss$object_annotation_category == "multiple"),"ellipt_BV"])
				per2 <- (ellipt_BV_multiple/ellipt_BV_livings)*100
				
				# Return
				d <- data.frame(sample = s, percentage_multiples = per, percentage_multiples_BV = per2)
				return(d)
				gc()
				rm(per, per2, nmultiple, nlivings, ellipt_BV_multiple, ellipt_BV_livings, ellipt_BV, spher_BV, areas)
				
}) # eo lapply

percentages <- do.call(rbind, res)

write.table(x = percentages[which(percentages$percentage_multiples_BV > 10),], file = "tableau_pourcentages.txt", sep = ";")


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### Now, map the position of the platresse station, together with the others :-)
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
dir()

unique(all_projs_data$sample_stationid)
unique(all_projs_data$object_lon)
unique(all_projs_data$object_lat)

coords_stations <- data.frame(stations = unique(all_projs_data$sample_stationid), x = NA, y = NA)

for(s in coords_stations$stations) {
	
		x <- mean(all_projs_data[which(all_projs_data$sample_stationid == s), "object_lon"])
		y <- mean(all_projs_data[which(all_projs_data$sample_stationid == s), "object_lat"])
	
		coords_stations[which(coords_stations$stations == s), "x"] <- x
		coords_stations[which(coords_stations$stations == s), "y"] <- y
		
} # eo for loop

###
# lat = 48.46666666666667
# lon =	-4.841666666666667

quartz()
ggplot() + geom_point(aes(x = -4.84166, y = 48.4666), pch = 21, colour = "black", fill = "blue", size = 4) + 
	geom_text(aes(x = -4.84166 - 0.05, y = 48.4666 - 0.02, label = "Platresse", hjust = 0), size = 2) + 
	geom_point(aes(x = x, y = y), pch = 21, colour = "black", fill = "red", size = 4, data = coords_stations) + 
	geom_text(data = coords_stations, aes(x = x - 0.05, y = y - 0.02, label = as.character(stations), hjust = 0), size = 2) + 
	coast + theme_linedraw()


### 24/07/07: Add isobaths and stations
coords_stations <- get(load( "coords_stations_pnmir.Rdata"))
library("marmap")
iso <- getNOAA.bathy(lon1 = -5.7, lon2 = -4, lat1 = 47.9, lat2 = 48.77, res = 1, keep = F)
class(iso)
str(iso)
# Transform to data.frame
str(iso)
rownames(iso) # longitudes
colnames(iso) # latitudes
iso <- fortify(iso) # use fortify to convert to data frame

quartz()
ggplot() + geom_point(aes(x = x, y = y, colour = z), data = iso[which(iso$z <= 0),]) + 
		geom_contour(aes(x = x, y = y, z = z), data = iso, colour = "black") + 
		geom_point(aes(x = x, y = y), pch = 21, colour = "black", fill = "white", size = 3, data = coords_stations[c(1:13),]) + 
		coast + theme_light()

   
### To make a nicer map (interpolation)
library("plyr")
library("fields")
library("akima")
library("reshape2")
# Define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(iso$x) + prec/2, to = max(iso$x) - prec/2, by = 0.01),
		y <- seq(from = min(iso$y) + prec/2, to = max(iso$y) - prec/2, by = 0.01)
) # eo list

# Linear interpolation
i <- interp(x = iso$x, y = iso$y, z = iso$z, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# Smoothing for better map
i <- image.smooth(i, theta = 0.015)
# Convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "z"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]
	
# Map:		
quartz()
ggplot() + geom_point(aes(x = lon, y = lat, colour = z), data = out[which(out$z <= 0),]) + 
		geom_contour(aes(x = x, y = y, z = z), data = iso, colour = "black") + 
		geom_point(aes(x = x, y = y), pch = 21, colour = "black", fill = "white", size = 3, data = coords_stations[c(1:13),]) + 
		scale_colour_distiller(name = "Depth (m)", palette = "Blues") + coast + theme_light()		

### With different isobaths levels: use breaks		
# quartz()
plot <- ggplot() + geom_point(aes(x = lon, y = lat, colour = z), data = out[which(out$z <= 0),]) + 
		geom_contour(aes(x = x, y = y, z = z), data = iso, colour = "grey30", breaks = c(-5,-25,-50,-75,-100)) + 
		geom_point(aes(x = x, y = y), pch = 21, colour = "black", fill = "white", size = 3, data = coords_stations[c(1:13),]) + 
		scale_colour_distiller(name = "Depth (m)", palette = "Blues") + coast + theme_light() 

ggsave(plot = plot, filename = "isobaths.pdf", dpi = 300, width = 6, heigh = 6)		

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##### 07/04/2017: You will need to extract the satellite data of SST and [Chl-a] to set your data into their context --> retrieve the front !
### To do so, you can first identify the dates for which you will need to extract the satellite data..
dates <- unique(all_projs_data$object_date)
true_dates <- paste(day(ymd(dates)), month(ymd(dates)), year(ymd(dates)), sep = "_")
true_dates


##### 10/04/2017: Explore Sardines data, not for sharing. 
sar <- read.csv("prod_mens_sardine_OP_bolinche.csv", h = T, sep = ";", dec = ",")
head(sar) # ?
sar <- sar[,c(1:5)]

### Per month/year, compute fishing effort (number of boats since only info available), average capture weight (kg)

sar$date <- paste(sar$month, sar$year, sep = "_")

require("dplyr")
sarsar <- sar %>%
		 group_by(date) %>%
	 	 summarise(avg_weight = mean(weight_capt), sd_weight = sd(weight_capt), effort = length(unique(name))) 

# Check bias due to sampling effort
summary(lm(avg_weight ~ effort, data = sarsar))
# R-squared:  -0.0123 ; p-value: 0.9014 ### Same with log scale
quartz()
ggplot() + geom_point(aes(x = effort, y = log(avg_weight)), data = sarsar, pch = 21, colour = "black", fill = "blue", size = 3) + 
		   xlab("Number of fishing boats") + ylab("Average capture weight (Log(kg))") + theme_bw()
		   
quartz() 
ggplot(sar, aes(x = factor(month), y = weight_capt)) + theme_bw() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Month") + ylab("Capture weight (kg)")
### Clearly an increase in autumn...but sampling effort ?
summary(lm(avg_weight ~ effort, data = sarsar)) ### Nope, no effect

# Yearly:
sarsar <- sar %>%
		 group_by(factor(year)) %>%
	 	 summarise(avg_weight = mean(weight_capt), sd_weight = sd(weight_capt), effort = length(unique(name))) 
# Or through boxplots 
quartz() 
ggplot(sar, aes(x = factor(year), y = weight_capt)) + theme_bw() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Year") + ylab("Capture weight (kg)")
### Higher production in 2009 and 2010
# Sampling effort ?
summary(lm(avg_weight ~ effort, data = sarsar)) ### Still no effect

# Seasonal: 
sar$season <- NA
sar[which(sar$month %in% c(3:5)),"season"] <- "spring"
sar[which(sar$month %in% c(6:8)),"season"] <- "summer"
sar[which(sar$month %in% c(9:11)),"season"] <- "fall"
sar[which(sar$month %in% c(1,2,12)),"season"] <- "winter"

sarsar <- sar %>%
		 group_by(factor(season)) %>%
	 	 summarise(avg_weight = mean(weight_capt), sd_weight = sd(weight_capt), effort = length(unique(name))) 
# Or through boxplots 
quartz() 
ggplot(sar, aes(x = factor(season), y = weight_capt)) + theme_bw() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Season") + ylab("Capture weight (kg)")
# Sampling effort ?
summary(lm(avg_weight ~ effort, data = sarsar)) ### Still no effect :-)

### AND now, a time series: 
library("dplyr")
head( arrange(sar, month, year) ) # ok
sar_ts <- arrange(sar, year, month)
sar_ts$date <- as.character(paste(1,sar_ts$month, sar_ts$year, sep = "-"))
require("lubridate")
sar_ts$date <- lubridate::dmy(x = sar_ts$date)

quartz()
ggplot() + geom_point(aes(x = date, y = weight_capt), data = sar_ts) + theme_bw()

### Explore with 'pastecs'
library("pastecs")

# First test, fit a polynomial of degree 6:
pol6 <- lm(weight_capt ~ poly(date, degree = 6), data = sar_ts)
str(pol6)
POL6 <- cbind(sar_ts, values = pol6$fitted.values, residuals = pol6$residuals)
quartz()
ggplot() + 
  geom_point(aes(x = date, y = weight_capt), sar_ts, colour = "black") +
  geom_line(aes(x = date, y = values), data = POL6 , colour = "blue") + 
  geom_line(aes(x = date, y = residuals), data = POL6 , colour = "darkgreen") + 
  scale_y_continuous(name = "Capture weight (kg)") + theme_bw()
  
# Works, but looks horrible

# To deduce the optimal window size, we will compute the turnogram of the data, with the function turnogram. Use a variable of marbio that you would have retained using the Escoufier criterion.
quartz()
ggplot(data = sar_ts) + geom_point(aes(x = date, y = weight_capt)) + geom_smooth(aes(x = date, y = weight_capt), colour = "red") + theme_linedraw()
?turnogram
t <- turnogram(ts(sar_ts$weight_capt))
str(t)
T <- data.frame(Lag = t$interval, Info = t$info)
quartz()
ggplot(T) + 
  geom_area(aes(x = Lag, y = Info), fill = "black", colour = "black" , alpha = 0.3) +
  geom_hline(yintercept = 0) + geom_hline(yintercept = -log(0.05, base = 2), colour = "red") + 
  geom_hline(yintercept = -4.3, colour = "red") + geom_vline(xintercept = 5, colour = "blue") + theme_linedraw()

# AT what lag value do we reach the maximal information level ?
T[which(T$Info == max(T$Info)),"Lag"]
# 102
# The functions for series decomposition all start with dec. The moving average is computed with decaverage().
?decaverage
#  Decompose a single regular time series with a moving average filtering. Return a 'tsd' object. To decompose several time series at once, use ‘tsd()’ with the argument ‘method="average"’
plot(decaverage(sar_ts$weight_capt, order = 2))
plot(decaverage(sar_ts$weight_capt, order = 5))
plot(decaverage(sar_ts$weight_capt, order = 2, times = 2))
plot(decaverage(sar_ts$weight_capt, order = 2, times = 3))
plot(decaverage(sar_ts$weight_capt, order = 2, times = 10))

## Let's weight it : 
plot(decaverage(sar_ts$weight_capt, order = 2, times = 2, weights = c(1,2,4,2,1)))
RA <- decaverage(sar_ts$weight_capt, order = 2, times = 2, weights = c(1,2,4,2,1))
plot(RA)
str(RA)
head(RA)
## values and residuals are stored as components ; let's extract those with : 
# extract(RA, component = "filtered")
# extract(RA, component = "residuals")
fil <- extract(RA, component = "filtered")
res <- extract(RA, component = "residuals")
# Joining in a data.frame for ggploting : 
X <- data.frame(date = sar_ts$date, weight_capt = sar_ts$weight_capt, values = fil, residuals = res)

quartz()
ggplot(data = X) +
  geom_point(aes(x= date, y= weight_capt), colour = "black") +
  geom_line(aes(x = date, y = values), colour = "red") +
  geom_linerange(aes(x = date, ymin = values, ymax = residuals), colour = "red", alpha = 1/2) + theme_linedraw()

trend.test(X$values)     # no trend in the values computed by the weighted running average
trend.test(X$residuals)  # no trend in the residuals

### And try some eigen vectors filtering
# The function is decevf() which has two arguments you will be interested in:
# lag : the number of repetitions of the series
# axes : which PCA axes are used to predict the smoothed series
# How do you choose the lag argument? Choose a sensible lag for weight_capt
## --> autocorrelogram to see the variation scale 
acf(ts(sar_ts$weight_capt))
## --> a sensible lag for temperature would be ???

# Compute the EVF of sar_ts$weight_capt using axes 1 and 2. Inspect the resulting object using
# plot and str.
#EVF2 <- decevf(sar_ts$weight_capt, lag = 5, axes = c(1,2)) # to keep the second PC too
EVF3 <- decevf(sar_ts$weight_capt , lag = 5, axes = c(1)) # PC1 only is enough
#str(EVF)
#plot(EVF2)
plot(EVF3)
## The predictions are stored in the "filtered" component
# Extract successively the predictions using axes 1, 2, and 3. Plot all of them and the original data on the same graph. How would you interpret this graph ?
evf3 <- extract(EVF3, component = "filtered")
E <- data.frame(Date = sar_ts$date, weight_capt = sar_ts$weight_capt, EVF = evf3)
quartz()
ggplot(data = E) + geom_point(aes(x= Date, y= weight_capt), colour = "black") + geom_line(aes(x = Date, y = EVF), colour = "blue") + theme_linedraw()


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### 10/04/2017: Examine a satellite product of SST + Chl-a form the ncdf files you can download @ ocean.color...
# L3-level product from:
# https://oceancolor.gsfc.nasa.gov/cgi/l3?per=DAY&prd=SST4_sst4.nc&sen=A&res=4km&num=100&ctg=Standard&date=4Jun2013

# But first, need to get the exact dates 
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/Tests_04_04_17/")
proj_dir <- getwd()
projects <- dir()
#p <- "export_198_20170404_1847"

projs <- lapply(projects, function(p) {
	
			# Go to project dir
			message(paste("Doing project #", p, sep = ""))
			setwd(paste(proj_dir,"/",p,"/", sep = ""))
			all_samples <- dir()
			
			# Read and rbind the samples from each project
			samples <- lapply(all_samples, function(s) {
					
					# Load and return
					ss <- read.table(file = s, sep = '\t', header = T)
					return(ss)
					
			}) # eo 2nd lapply 
			
			all_samples <- do.call(rbind, samples)
			# dim(all_samples) ; summary(all_samples)
			rm(samples)
			
			return(all_samples)
}
) # eo lapply
all_projs_data <- do.call(rbind, projs)

### Restrict to dates where pnmir stations are present
library("lubridate") # Test lubridate package
all_projs_data$object_date <- as.character(all_projs_data$object_date)
all_projs_data$year_object <- year(ymd(all_projs_data$object_date))
all_projs_data$month_object <- month(ymd(all_projs_data$object_date))
all_projs_data$day_object <- day(ymd(all_projs_data$object_date))
# summary(all_projs_data[,c(160:162)])
dates <- unique(all_projs_data[which(all_projs_data$sample_stationid == "b7"),"object_date"])
dates # Y-M-D
# 2013 06 04
# 2015 07 09
# 2015 10 09
# 2011 05 02
# 2011 07 04
# 2011 10 12
# 2012 05 02

setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
dir()
library("ncdf4")
library("raster")
library("oce")

nc <- nc_open("A2013155.L3m_DAY_SST4_sst4_4km.nc")
print(nc)

### Get coords
lon <- ncvar_get(nc, "lon")
nlon <- dim(lon) # 8640
head(lon)
lat <- ncvar_get(nc, "lat")
nlat <- dim(lat) # 4320
head(lat)

### And sea surface temperature ?
str(nc$var)
temp <- ncvar_get(nc, "sst4")
dim(temp); class(temp)
rownames(temp) <- lon
colnames(temp) <- lat

# Identify the colnames that range between 47.9 and 48.6
lats <- as.character( lat[lat > 47.9 & lat < 48.6] )
lons <- as.character( lon[lon > -6 & lon < -4.25] )
# filter latitudes
temp <- temp[,lats]
# filter longitudes
temp <- temp[lons,]

#temp <- as.data.frame(temp)
#temp$id <- rownames(temp)

t <- data.frame(melt(temp))
colnames(t)[1:3] <- c("x","y","SST")

# Change land sst values to NaN
t$SST[t$SST >= 45] <- NA

# Map sst + stations
# coords_stations
quartz()
ggplot() + geom_tile(aes(x = x, y = y, fill = SST), data = t) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_fill_distiller(name = "SST (°C)", palette = "Spectral") + 
			#scale_fill_gradientn(name = "SST (°C)", colours = rev(jet.colors(7) )) + 
			coast + coord_quickmap() + theme_linedraw()

### Test with the L2 regional product from:
# https://oceancolor.gsfc.nasa.gov/cgi/browse.pl

nc <- nc_open("A2013155025000.L2_LAC_SST4.nc")
# print(nc)
names(nc$var) # Names of the variables contain some attributes

### Get coords
lon <- ncvar_get(nc, "navigation_data/longitude")
dim(lon) # 1354 2030
lon <- melt(lon)
lat <- ncvar_get(nc, "navigation_data/latitude")
dim(lat) # 1354 2030
lat <- melt(lat)

### And sea surface temperature ?
#str(nc$var)
temp <- ncvar_get(nc, "geophysical_data/sst4")
dim(temp); class(temp)
temp <- melt(temp)
head(temp) ; dim(temp)


t <- na.omit(data.frame(lon = lon[,3], lat = lat[,3], sst = temp[,3]))
colnames(t)[1:3] <- c("x","y","SST")
# summary(t)
t <- t[which(t$x >= -5.5),]
t <- t[which(t$x <= -4.0),]
t <- t[which(t$y >= 47.9),]
t <- t[which(t$y <= 48.8),]

# Map sst + stations
# coords_stations
quartz()
ggplot() +  geom_point(aes(x = x, y = y, colour = SST), data = t) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			scale_colour_distiller(name = "SST (°C)", palette = "RdYlBu", limits = c(10,16)) + 
			#scale_colour_gradientn(name = "SST (°C)", colours = jet.colors(7)) + 
			coast + coord_quickmap() + theme_linedraw()

### To make a nicer map
require("plyr")
library("fields")
library("akima")
library("reshape2")
# define resulting interpolation grid
prec <- 0.01
grid <- list( 	
		x <- seq(from = min(t$x) + prec/2, to = max(t$x) - prec/2, by = 0.01),
   	 	y <- seq(from = min(t$y) + prec/2, to = max(t$y) - prec/2, by = 0.01)
) # eo list

# bivariate interpolation
i <- interp(x = t$x, y = t$y, z = t$SST, xo = grid[[1]], yo = grid[[2]], linear = TRUE)
# -> less artifacts near the coast
# smoothing for better map
i <- image.smooth(i, theta = 0.015)
# convert to data.frame
out <- melt(i$z, varnames = c("start_lon", "start_lat") )
names(out)[3] <- "sst"
out$lon <- i$x[out$start_lon]
out$lat <- i$y[out$start_lat]

quartz()
ggplot() +  geom_raster(aes(x = lon, y = lat, fill = sst), data = out) + 
			geom_point(aes(x = x, y = y), data = coords_stations, pch = 21, fill = "white", size = 2, colour = "black") + 
			#scale_fill_distiller(name = "SST (°C)", palette = "Spectral", limits = c(10,15)) + 
			scale_fill_gradientn(name = "SST (°C)", colours = rev(jet.colors(7)), limits = c(11,14)) + 
			coast + coord_quickmap() + theme_linedraw()
