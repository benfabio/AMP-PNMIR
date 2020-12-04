
##### 28/04/2017: R Script to examine the content and structure of the ZooScan data from the AMP PMMI project © Fabio Benedetti, OOV-UMS, AFB
##### Aims to:

#	- Load the abundances and biovolumes you extracted form ecotaxa (Script#4_AMP_pnmir_zoo.R)
#	- Select the objects/stations of interest: V versus H, PNMIR transects etc.
#	- Select the variables (biological taxa) of interest
#	- cbind() with metadata and environmental variables
#	- Community analysis (CA, RDA, CCA) with both abund and biovolumes
#	- Check if there are any big differences between H and V plankton tows
#	- Consider different taxonomic groups

### Latest update: 12/06/2017

library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
require("dplyr")
library("fields")
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
ggplot() + geom_point(aes(x=x, y=y), data = coords_stations, pch = 21, fill = "lightskyblue", size = 3) + coast + theme_linedraw()
# OKAY.


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##### 28/04/2017: Load the ZooScan data
zoo <- read.csv("wp2_allecotaxa_abund_biovol_28_04_17.csv", h = TRUE, sep = ",", dec = ".")
# Examine
head(zoo)
colnames(zoo)
str(zoo)
rm(zoo)

### If you want abundances : 
zoo <- read.csv("wp2_all_ecotaxa_abund_biovol_21_06_17.csv", h = TRUE, sep = ",", dec = ".")
# or biovolumes (but let's concentrate on abundances for now): 
# zoo <- read.csv("wp2_allecotaxa_biovol_28_04_17.csv", h = TRUE, sep = ",", dec = ".")

### Filter the objects
ids <- unique(zoo$orig_id)
molene <- ids[grep("molene",ids)] # 136 stations
douarnenez <- ids[grep("douarnenez",ids)] # 68 stations
length(ids) - (length(douarnenez) + length(molene))
### pnmir --> 198 stations ; get rid of what belongs to molene and douarnenez
zoo <- zoo[!(zoo$orig_id %in% douarnenez),]
zoo <- zoo[!(zoo$orig_id %in% molene),]
# Check:
unique(zoo$orig_id)
dim(zoo) 
### Good. Now, get rid of Year 2010 and horizontal plankton hauls
wp2010 <- ids[grep("2010",ids)]
horizontal <- ids[grep("_h",ids)]
zoo <- zoo[!(zoo$orig_id %in% wp2010),]
zoo <- zoo[!(zoo$orig_id %in% horizontal),]
# unique(zoo$orig_id)
dim(zoo) # kk


# Split orig_id in order to get the date and station of sampling
unique(zoo$orig_id)
ID <- data.frame( do.call(rbind, strsplit(x = as.character(zoo$orig_id), split = "_")) )
colnames(ID)[1:4] <- c("net", "station", "date", "type")
zoo$station <- ID$station
ID$date <- lubridate::ymd(ID$date)
zoo$date <- ID$date
zoo$month <- lubridate::month(zoo$date)
zoo$year <- lubridate::year(zoo$date)


### 02/05/2017: Examine the taxa's absolure and relative contributions to total abundance and biovolume (bar charts), according to horizontal and vertical net tows.
colnames(zoo)
colnames(zoo)[4:163] <- stringr::str_replace_all(colnames(zoo)[4:163], "concentration_", "")
zoo$other_Calanoida <- zoo$Calanoida - ( zoo$Acartiidae + zoo$Calanidae + zoo$Candaciidae + zoo$Centropagidae + zoo$Temoridae + zoo$Euchaetidae)
# Same but for biovolumes:
# colnames(zoo)[158:312] <- stringr::str_replace_all(colnames(zoo)[158:312], "biovolume_", "")
# First, melt:
mzoo <- melt(zoo[2:length(zoo)], id.vars = c("station","date","month","year","orig_id"))
head(mzoo)
colnames(mzoo)[6:7] <- c("group","abund")
str(mzoo)
require("dplyr")
groups <- unique(mzoo$group)[(grep("biovolume", unique(mzoo$group)))]
contribs <- data.frame(mzoo[!(mzoo$group %in% groups),] %>% group_by(group) %>% summarise(total_abund = round(sum(abund, na.rm = T),3), mean_abund = round(mean(abund, na.rm = T),3) ) ) 
contribs[order(contribs$mean_abund, decreasing = T),]

write.table(contribs[order(contribs$mean_abund, decreasing = T),], file = "table_abund_all_groups_18_02_19.txt", sep = ";")

#d <- data.frame(d$year)
mzoo$abund <- round(mzoo$abund, 3)

# Sort by decreasing order of total abund
mzoo <- mzoo[order(mzoo$abund, decreasing = T),]
# force factor levels to keep this order : factor() !
mzoo$abund <- factor(mzoo$abund, levels = mzoo$abund)
quartz()
plot <- ggplot(d) + geom_bar(aes(x = variable, y = total_biovol), stat = "identity") + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("Classes") + ylab("Total biovolume (mm3)")
 
#ggsave(plot = plot, filename = "histogram_biovol_all_classes_v.pdf", dpi = 300, height = 7, width = 19) 
 
 
 
### Filter the variables according to histograms above ?
colnames(zoo)
# 28/04/05: keep 1,2,6:15,17:18,20:23,25:51,53:68,70:73,75:85,87:89,91:108,110:122,126:134,136:140,142,144:148,150:158,160
rownames(zoo) <- c(1:nrow(zoo)) 
sums <- rowSums(as.matrix(zoo[,c(165:324)]))
ids <- as.numeric(names(sums[which(sums == 0)]))
zoo[ids,]

### 03/05/2017: Draw the abudnance histograms for the selected taxa to assess their absolute/ relative contributions.	
### 19/05/2017: Consider supplementary taxa that may help you to better represent the coast-offshare gradient: Bryozoa, eggs, Bivalvia, Echinodermata...
zoo <- zoo[,c("month","year","station",
			"Calanoida","Cladocera","Oikopleuridae","Hydrozoa","Cyclopoida","nauplii","Cirripedia", "Thecosomata","Poecilostomatoida","Bivalvia",
			"Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida","Acartiidae","Echinodermata","zoea",
			"Calanidae","Candaciidae","Centropagidae","Corycaeidae","Oithonidae","Euchaetidae","Oncaeidae","Temoridae","Sapphirinidae","Decapoda","egg","Bryozoa")]

### !!! Echinodermata have too much weight on the CA --> ignore it...
			
# First, melt:
zoo$other_Calanoida <- zoo$Calanoida - ( zoo$Acartiidae + zoo$Calanidae + zoo$Candaciidae + zoo$Centropagidae + zoo$Temoridae + zoo$Euchaetidae)

mzoo <- melt(zoo, id.vars = c("station","month","year"))
head(mzoo)

d <- data.frame(mzoo[-which(mzoo$variable %in% c("Poecilostomatoida","Calanoida","Cyclopoida")),] %>% 
	 group_by(variable) %>% 
	 summarise( total_abund = sum(value) ))
	 
#n <- length(unique(mzoo[which(mzoo$year == 2013),"station"]))
	 
### !!! Watch out: need to correct the possible differences	in number of sampling stations per year !!! 
d$total_abund <- round(d$total_abund, 3) 
# Sort by decreasing order of total abund
d <- d[order(d$total_abund, decreasing = T), ]
# force factor levels to keep this order : factor() !
d$variable <- factor(d$variable, levels = d$variable)
quartz()
ggplot(d) + theme_light() + geom_bar(aes(x = variable, y = total_abund), stat = "identity") + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Classes") + ylab("Total biovolume (mm3) - all years and stations")
		

### Preliminary analyses : CA on important zooplankton groups
library("FactoMineR")
colnames(zoo)
# summary( zoo[,c(5:7,9:11,13:21,23:length(zoo))] )
res.ca <- CA( zoo[,c(5:7,9:11,13:21,23:length(zoo))] )
summary(res.ca) 
### 70.12% of total variance along 4 Axes when vertical hauls
### 76.97% of total variance along 4 Axes when horizontal hauls
# 02/05/2017 : 65.34% of total variance along 4 Axes when vertical hauls; abundances
# 19/05/2017: 56%

### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
					Ax4 = res.ca$col$coord[,4],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2], 
					Cos2_2 = res.ca$col$cos2[,3] + res.ca$col$cos2[,4])

# For objects
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
					Ax4 = res.ca$row$coord[,4],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					Cos2_2 = res.ca$row$cos2[,3] + res.ca$row$cos2[,4],
					month = zoo$month, station = zoo$station, year = zoo$year )

head(AFCst)

quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year), shape = factor(year)), data = AFCst, alpha = 0.7, colour = "black") +
	scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()


quartz()
ggplot() +
	geom_point(aes(x = Ax3, y = Ax4), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
	geom_point(aes(x = Ax3, y = Ax4, size = Cos2_2, fill = factor(month), shape = factor(year)), data = AFCst, alpha = 0.7, colour = "black") +
	scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 3 (",round(res.ca$eig$per[3], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 4 (",round(res.ca$eig$per[4], 2),"%)", sep = "")) + 
	geom_text(aes(x = Ax3, y = Ax4+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()
	
	
### cbind() with other PNMIR data (env and then phyto counts)
write.table(zoo, "wp2_biovol_v_21_06_17.txt", sep = ";")

### 03/05/2017:
env <- read.csv("PNMIR_cruise_data.csv", h = T, sep = ";", dec = ",")
dim(env)
colnames(env)

#zoo_abund <- read.table("wp2_abund_v_12_06_17.txt", sep = ";")
#dim(zoo_abund)
#colnames(zoo_abund)

### Use metadata to combine both, according to date and station
zoo$id <- paste(zoo$station, zoo$month, zoo$year, sep = "_")
env$id <- paste(env$station, env$month, env$year, sep = "_")
# Add empty columms for abundances
env[,colnames(zoo)[4:36]] <- NA

matching <- na.omit(zoo$id[match(env$id, zoo$id)]) ; length(matching)

for(i in matching) {
	
			# Doing : 
			message(paste(i, sep = ""))
			# Filter biovolumes according to id
			abundances <- zoo[which(zoo$id == i),c(4:36)]
			# cbind with env
			env[which(env$id == i),colnames(zoo)[4:36] ] <- abundances
		
} # eo for loop

summary(env)
# Drop env$id
env <- env[,-c(44)]

write.table(env, "PNMIR_cruise_data_biovol_v_21_06_17.csv", sep = ";")
### OKAY


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##### 02/05/2017: Explore relationships between zooplankton abundances/ biovolumes and the environment + phytoplankton community structure...
ddf <- read.table("PNMIR_cruise_data_biovol_v_21_06_17.csv", h = T, sep = ";", dec = ".")
colnames(ddf) ; rownames(ddf)
str(ddf)
summary(ddf)

# sums <- rowSums(as.matrix(ddf[,c(44:76)]))
# ids <- as.numeric(names(sums[which(sums == 0)]))
# ddf[ids,]

# Get rid of 2010:
ddf <- ddf[which(ddf$year != 2010),]
# Difference between surface and bottom temperatures
ddf$Tdiff <- ddf$Ts - ddf$Tf

# CCA ?
colnames(ddf)
cctable <- ddf[,c(4:9,44:76)]
						
dim(cctable)
cctable
colnames(cctable)

# summary(lm(max_depth ~ x, data = cctable)) # R-squared: 0.7567

##### 05/05/2017: Perform CA and/or PCA Beware of not mixing imbricated classes (like Calanoida and Calanidae)
library("FactoMineR")
# But first, create a columns for undetermined Calanoida
# cctable$other_Calanoida <- cctable$Calanoida - ( cctable$Acartiidae + cctable$Calanidae + cctable$Candaciidae + cctable$Centropagidae + cctable$Temoridae + cctable$Euchaetidae)
# Check out its percentage
summary((cctable$other_Calanoida / cctable$Calanoida)*100) 
#    Min.  1st Qu.  Median   Mean   3rd Qu.   Max.    NA's 
#   11.05   51.64    64.82   63.07   80.27   100.00    31
### Are quite important !
summary((cctable$other_Calanoida / cctable$Copepoda)*100) 
#    Min. 1st Qu.  Median  Mean   3rd Qu.   Max.   NA's 
#   5.00   37.41   51.97   50.79   63.82   88.73    31 
### Often half of all the copepoda !

# CA: 
catable <- na.omit(cctable)
colnames(catable)
summary(catable)
# Perform CA
res.ca1 <- CA( catable[,c(8:10,12:14,16:24,26:39)], ncp = 10) # Without Echinodermata

##### ISSUES sometimes, check if it's rowSums related
# library("matrixStats")
# rownames(catable) <- c(1:nrow(catable))
# sums <- rowSums(as.matrix(catable[,c(8:10,12:14,16:24,26:39)]))
# ids <- as.numeric(names(sums[which(sums == 0)]))
# catable[ids,]

summary(res.ca1)
# Kaiser-Guttman criterion: 
# Looking at the eigenvalues
eig <- data.frame(prop = res.ca1$eig$per, nb = c(1:31))
quartz()
ggplot(eig) + geom_bar(aes(x=nb, y=prop), stat="identity") + geom_line(aes(x=nb, y=mean(prop) )) + xlab("Eigenvalue Number") +  ylab("Value") + theme_linedraw()

# 6 CA axes --> 83.72% of total variance
AFCsp <- data.frame(Ax1 = res.ca1$col$coord[,1],
                    Ax2 = res.ca1$col$coord[,2],
					Ax3 = res.ca1$col$coord[,3],
					Ax4 = res.ca1$col$coord[,4],
					Ax5 = res.ca1$col$coord[,5],
					Ax6 = res.ca1$col$coord[,6],
                    Cos2_1 = res.ca1$col$cos2[,1] + res.ca1$col$cos2[,2], 
					Cos2_2 = res.ca1$col$cos2[,3] + res.ca1$col$cos2[,4],
					Cos2_3 = res.ca1$col$cos2[,5] + res.ca1$col$cos2[,6])

# For objects
AFCst <- data.frame(Ax1 = res.ca1$row$coord[,1],
                    Ax2 = res.ca1$row$coord[,2],
					Ax3 = res.ca1$row$coord[,3],
					Ax4 = res.ca1$row$coord[,4],
					Ax5 = res.ca1$row$coord[,5],
					Ax6 = res.ca1$row$coord[,6],
                    Cos2_1 = res.ca1$row$cos2[,1] + res.ca1$row$cos2[,2], 
					Cos2_2 = res.ca1$row$cos2[,3] + res.ca1$row$cos2[,4],
					Cos2_3 = res.ca1$row$cos2[,5] + res.ca1$row$cos2[,6],
					month = catable$month, station = catable$station, year = catable$year )

head(AFCst)

# quartz()
ca1 <- ggplot() +
  		geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  		geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year)), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
		#scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		scale_x_continuous(paste("CA 1 (",round(res.ca1$eig$per[1], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 2 (",round(res.ca1$eig$per[2],2),"%)", sep = "")) + 
		geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_light()


#quartz()
ca2 <- ggplot() +
		geom_point(aes(x = Ax3, y = Ax4), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
		geom_point(aes(x = Ax3, y = Ax4, size = Cos2_2, fill = factor(year)), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
		#scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		scale_x_continuous(paste("CA 3 (",round(res.ca1$eig$per[3], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 4 (",round(res.ca1$eig$per[4], 2),"%)", sep = "")) + 
		geom_text(aes(x = Ax3, y = Ax4+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_bw()


#quartz()
ca3 <- ggplot() +
		geom_point(aes(x = Ax5, y = Ax6), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
		geom_point(aes(x = Ax5, y = Ax6, size = Cos2_3, fill = factor(month)), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
		scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		scale_x_continuous(paste("CA 5 (",round(res.ca1$eig$per[5], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 6 (",round(res.ca1$eig$per[6], 2),"%)", sep = "")) + 
		geom_text(aes(x = Ax5, y = Ax6+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_bw()

ggsave(plot = ca1, filename = "CA_biovol_zoo_V_12.pdf", dpi = 300, width = 9, height = 7)
ggsave(plot = ca2, filename = "CA_biovol_zoo_V_34.pdf", dpi = 300, width = 9, height = 7)
ggsave(plot = ca3, filename = "CA_biovol_zoo_V_56.pdf", dpi = 300, width = 9, height = 7)


### Same as above, but by considering the Copepoda families
#catable <- na.omit(cctable[,c(1:32)])
colnames(catable)
res.ca2 <- CA( catable[,c(24,27:35)], ncp = 6)
# 6 CA axes --> 75.17% of total variance
AFCsp <- data.frame(Ax1 = res.ca2$col$coord[,1],
                    Ax2 = res.ca2$col$coord[,2],
					Ax3 = res.ca2$col$coord[,3],
					Ax4 = res.ca2$col$coord[,4],
					Ax5 = res.ca2$col$coord[,5],
					Ax6 = res.ca2$col$coord[,6],
                    Cos2_1 = res.ca2$col$cos2[,1] + res.ca2$col$cos2[,2], 
					Cos2_2 = res.ca2$col$cos2[,3] + res.ca2$col$cos2[,4],
					Cos2_3 = res.ca2$col$cos2[,5] + res.ca2$col$cos2[,6])

# For objects
AFCst <- data.frame(Ax1 = res.ca2$row$coord[,1],
                    Ax2 = res.ca2$row$coord[,2],
					Ax3 = res.ca2$row$coord[,3],
					Ax4 = res.ca2$row$coord[,4],
					Ax5 = res.ca2$row$coord[,5],
					Ax6 = res.ca2$row$coord[,6],
                    Cos2_1 = res.ca2$row$cos2[,1] + res.ca2$row$cos2[,2], 
					Cos2_2 = res.ca2$row$cos2[,3] + res.ca2$row$cos2[,4],
					Cos2_3 = res.ca2$row$cos2[,5] + res.ca2$row$cos2[,6],
					month = catable$month, station = catable$station, year = catable$year, depth = catable$max_depth)

head(AFCst)

quartz()
ca1 <- ggplot() +
  		geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  		geom_point(aes(x = Ax1, y = Ax2, size = depth, fill = factor(month)), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
		#scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		scale_x_continuous(paste("CA 1 (",round(res.ca2$eig$per[1], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 2 (",round(res.ca2$eig$per[2],2),"%)", sep = "")) + 
		geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_bw()


#quartz()
ca2 <- ggplot() +
		geom_point(aes(x = Ax3, y = Ax4), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
		geom_point(aes(x = Ax3, y = Ax4, size = Cos2_2, fill = factor(month)), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
		#scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		scale_x_continuous(paste("CA 3 (",round(res.ca$eig$per[3], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 4 (",round(res.ca$eig$per[4], 2),"%)", sep = "")) + 
		geom_text(aes(x = Ax3, y = Ax4+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_bw()


#quartz()
ca3 <- ggplot() +
		geom_point(aes(x = Ax5, y = Ax6), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
		geom_point(aes(x = Ax5, y = Ax6, size = Cos2_3, fill = factor(month)), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
		scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		scale_x_continuous(paste("CA 5 (",round(res.ca$eig$per[5], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 6 (",round(res.ca$eig$per[6], 2),"%)", sep = "")) + 
		geom_text(aes(x = Ax5, y = Ax6+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_bw()

ggsave(plot = ca1, filename = "CA2_abund_zoo_V_12.pdf", dpi = 300, width = 9, height = 7)
ggsave(plot = ca2, filename = "CA2_abund_zoo_V_34.pdf", dpi = 300, width = 9, height = 7)
ggsave(plot = ca3, filename = "CA2_abund_zoo_V_56.pdf", dpi = 300, width = 9, height = 7)



### PCA:
pcatable <- na.omit( ddf[,c("month","year","station","x","y","max_depth","Ts","Tdiff","NO3","PO4","SiOH4","Chla","Phaeo_a",
			"Calanoida","Cladocera","Oikopleuridae","Hydrozoa","Cyclopoida","nauplii","Cirripedia", "Thecosomata","Poecilostomatoida","Bivalvia",
			"Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida","Acartiidae","Echinodermata","zoea",
			"Calanidae","Candaciidae","Centropagidae","Corycaeidae","Oithonidae","Oncaeidae","Temoridae","Sapphirinidae","Decapoda","egg","Bryozoa","Euchaetidae")] )
						
dim(pcatable)
colnames(pcatable)
# With env and zoo
res.pca <- PCA(X = pcatable[,c(6:length(pcatable))], scale.unit = TRUE, quanti.sup = c(9:40))

# With zoo only --> PC1 will separate the stations that present high abundances from the ones that have rather low abundances
res.pca <- PCA(X = pcatable[,c(14:length(pcatable))], scale.unit = TRUE )

# With copepods only + env
res.pca <- PCA(X = pcatable[,c(6:13,14,30:38)], scale.unit = TRUE, quanti.sup = c(9:18) )

# With copepods only:
res.pca <- PCA(X = pcatable[,c(14,30:38)], scale.unit = TRUE )


### 05/05/17: Perform CCA
cctable <- na.omit( ddf[,c("month","year","station","x","y","max_depth","Ts","Tdiff","NO3","SiOH4","Chla","Phaeo_a",
			"Calanoida","other_Calanoida","Cladocera","Oikopleuridae","Hydrozoa","Cirripedia","Thecosomata","Bivalvia",
			"Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida","Acartiidae","zoea",
			"Calanidae","Candaciidae","Centropagidae","Corycaeidae","Oithonidae","Oncaeidae","Temoridae","Sapphirinidae","Decapoda","Bryozoa","Euchaetidae")] )
				
# cctable$other_Calanoida <- cctable$Calanoida - ( cctable$Acartiidae + cctable$Calanidae + cctable$Candaciidae + cctable$Centropagidae + cctable$Temoridae + cctable$Euchaetidae)
colnames(cctable)
library("vegan")
cca <- vegan::cca( X = as.matrix(cctable[,c(14:40)]), Y = as.matrix(cctable[,c(6:12)]) ) # X = abundances ; Y = env
summary(cca, scaling = 1, axes = 3) 
#                         CCA1    CCA2    CCA3    CCA4    CCA5     CCA6     CCA7
# Eigenvalue            0.1311 0.07851 0.03399 0.01747 0.01137 0.004739 0.003052
# Proportion Explained  0.4679 0.28014 0.12128 0.06232 0.04057 0.016910 0.010890
# Cumulative Proportion 0.4679 0.74803 0.86931 0.93163 0.97220 0.989110 1.000000

# test significance: 
anova(cca, step = 1000, perm.max = 1000)
#               CCA1     CCA2     CCA3
# max_depth  0.81190 -0.31684  0.06751
# Ts         0.30095  0.67112  0.13125
# Tdiff      0.30863  0.02165  0.44471
# NO3       -0.09203 -0.08618 -0.21717
# SiOH4      0.31200  0.48604 -0.31457
# Chla      -0.49707 -0.25924 -0.14024
# Phaeo_a   -0.31965 -0.19718 -0.16881

# Significant WHEN NOT accounting for Echinodermata and Noctiluca

# NOT significant when based on biovolumes

species_scores <- data.frame(species = rownames(cca$CCA$v), CCA1 = cca$CCA$v[,"CCA1"], CCA2 = cca$CCA$v[,"CCA2"], CCA3 = cca$CCA$v[,"CCA3"])
stations_scores <- data.frame(CCA1 = cca$CCA$u[,"CCA1"], CCA2 = cca$CCA$u[,"CCA2"], CCA3 = cca$CCA$u[,"CCA3"])
stations_scores$month <- cctable$month
stations_scores$year <- cctable$year
stations_scores$station <- cctable$station
vars_scores <- data.frame(var = rownames(cca$CCA$biplot), CCA1 = cca$CCA$biplot[,"CCA1"], CCA2 = cca$CCA$biplot[,"CCA2"], CCA3 = cca$CCA$biplot[,"CCA3"] )

### Here we go...
quartz()
plot(cca, display=c("lc","cn"), scaling=1, main="Biplot CCA fish ~ env2 - scaling 1") 
quartz()
plot(cca, display=c("sp","cn"), main="Biplot CCA fish ~ env2 - scaling 2") 

quartz()
plot <- ggplot() + geom_point(aes(x = CCA1, y = CCA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size = 4) + 
		   	geom_point(aes(x = CCA1, y = CCA2, fill = factor(year)), data = stations_scores, pch = 23, colour = "black", size = 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = CCA1*3, yend = CCA2*3), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour = "black", data = vars_scores) + 
		   	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		   	geom_text(aes(x = CCA1, y = CCA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + xlab("CCA 1 (47%)") + ylab("CCA 2 (36%)") + 
			scale_x_continuous(limits = c(-3,3)) + scale_y_continuous(limits = c(-6,4)) + theme_light()

ggsave(plot = plot, filename = "CCA3_zoo_abund_v.pdf", dpi = 300, width = 9, height = 7)

### CCA3 too: 
quartz()
ggplot() + geom_point(aes(x = CCA2, y = CCA3), data = species_scores, pch = 21, fill = "grey70", colour = "black", size = 4) + 
		   	geom_point(aes(x = CCA2, y = CCA3, fill = factor(month)), data = stations_scores, pch = 23, colour = "black", size = 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = CCA2*3, yend = CCA3*3), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour = "black", data = vars_scores) + 
		   	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		   	geom_text(aes(x = CCA2, y = CCA3+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + xlab("CCA 1 (35.8%)") + ylab("CCA 2 (29.3%)") + 
			scale_x_continuous(limits = c(-5,1.7)) + scale_y_continuous(limits = c(-3,9.5)) + theme_light()



### Perform RDA (since abundances) with transformed abudnances
# shapiro.test(cctable$Calanoida) # Non normal distribution, need to transform the data first
cctable[,c(14:40)] <- vegan::decostand(cctable[,c(14:40)], "hellinger")

### Try CA on transformed data:
library("FactoMineR")
res.ca <- CA( cctable[,c(14:29)] )
summary(res.ca) # 64% of first 4 axes with abund
# With biovolumes ? 59.1% of total variance
# With vertical tows' abund: 56.62% of total var on 4 axes

### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
					Ax4 = res.ca$col$coord[,4],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2], 
					Cos2_2 = res.ca$col$cos2[,3] + res.ca$col$cos2[,4])

# For objects
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
					Ax4 = res.ca$row$coord[,4],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					Cos2_2 = res.ca$row$cos2[,3] + res.ca$row$cos2[,4],
					month = cctable$month, station = cctable$station, year = cctable$year )

head(AFCst)

quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month), shape = factor(year)), data = AFCst, alpha = 0.7, colour = "black") +
	scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()


quartz()
ggplot() +
	geom_point(aes(x = Ax3, y = Ax4), data = AFCsp, fill = "grey65", alpha = 0.7, pch = 21, colour = "black", size = 3)+
	geom_point(aes(x = Ax3, y = Ax4, size = Cos2_2, fill = factor(month), shape = factor(year)), data = AFCst, alpha = 0.7, colour = "black") +
	scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 3 (",round(res.ca$eig$per[3], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 4 (",round(res.ca$eig$per[4], 2),"%)", sep = "")) + 
	geom_text(aes(x = Ax3, y = Ax4+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_linedraw()
	

### Try RDA now
res.rda <- vegan::rda( X = as.matrix(cctable[,c(14:40)]), Y = as.matrix(cctable[,c(6:12)]) )
summary(res.rda)
##### 12/06/2017: When considering families of Copepoda
#                        RDA1    RDA2     RDA3     RDA4     RDA5     RDA6
# Eigenvalue            0.04327 0.01448 0.005851 0.002262 0.001392 0.000764
# Proportion Explained  0.63485 0.21239 0.085850 0.033190 0.020420 0.011210
# Cumulative Proportion 0.63485 0.84724 0.933090 0.966280 0.986710 0.997920
anova(res.rda, step = 1000, perm.max = 1000)
# OKAY.


### 22/06/2017: With biovolumes
#                            RDA1    RDA2     RDA3     RDA4     RDA5     RDA6
# Eigenvalue            0.04967 0.01523 0.005973 0.003344 0.002121 0.001175
# Proportion Explained  0.63886 0.19592 0.076830 0.043000 0.027270 0.015110
# Cumulative Proportion 0.63886 0.83478 0.911610 0.954620 0.981890 0.997000

 
#             RDA1     RDA2     RDA3      RDA4    RDA5     RDA6
# max_depth  0.8078  0.15219 -0.47710  0.056438  0.2885 -0.10073
# Ts        -0.1092  0.90003  0.17574 -0.332762 -0.1720 -0.04752
# Tdiff      0.2268  0.53668 -0.57792  0.236256 -0.4954  0.13778
# NO3        0.2617 -0.39746  0.18288 -0.256229  0.3023  0.75933
# SiOH4      0.3596  0.08080  0.78893  0.002606  0.1453  0.46921
# Chla      -0.3160  0.19269 -0.28605  0.539698  0.6168  0.32645
# Phaeo_a   -0.3462 -0.08096 -0.06239  0.387244  0.3641  0.20567


### 12/06/2017: With horizontal tows' abundances, much better than CCA actually

plot(res.rda, scaling = 1)
plot(res.rda, scaling = 2)

# The R2 of a RDA is "biased", like the ordinary R2 of a multiple regression. However, an "adjusted" R2 ratio. An adjusted R2 near 0 indicates that X does not explain more of the variation of Y than
# random normal variables would do. Adjusted R2 can be negative, indicating that the explanatory variables X do worse than random normal variables would.
# Both R2 and adjusted R2 can be computed using vegan function RsquareAdj()
vegan::RsquareAdj(x = res.rda)  
# $r.squared
#   [1] 0.3080854
#   $adj.r.squared
#   [1] 0.2586629

# With biovolumes
# $r.squared
# [1] 0.2547845
# $adj.r.squared
# [1] 0.2015549

# For ploting:
species_scores <- data.frame(species = rownames(res.rda$CCA$v), RDA1 = res.rda$CCA$v[,"RDA1"], RDA2 = res.rda$CCA$v[,"RDA2"], RDA3 = res.rda$CCA$v[,"RDA3"])
stations_scores <- data.frame(RDA1 = res.rda$CCA$u[,"RDA1"], RDA2 = res.rda$CCA$u[,"RDA2"], RDA3 = res.rda$CCA$u[,"RDA3"])
stations_scores$month <- cctable$month
stations_scores$year <- cctable$year
stations_scores$station <- cctable$station
vars_scores <- data.frame(var = rownames(res.rda$CCA$biplot), RDA1 = res.rda$CCA$biplot[,"RDA1"], RDA2 = res.rda$CCA$biplot[,"RDA2"], RDA3 = res.rda$CCA$biplot[,"RDA3"] )

# For choosing the limits of the axes
summary(species_scores) ; summary(vars_scores)

quartz()
plot <- ggplot() + geom_point(aes(x = RDA1, y = RDA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size = 4) + 
		   	geom_point(aes(x = RDA1, y = RDA2, fill = factor(year)), data = stations_scores, pch = 23, colour = "black", size = 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = RDA1*0.7, yend = RDA2*0.7), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour = "black", data = vars_scores) + 
		   	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		   	geom_text(aes(x = RDA1, y = RDA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
			#scale_x_continuous(limits = c(-0.5,0.6)) + scale_y_continuous(limits = c(-0.7,0.5)) + 
			xlab("RDA 1 (63.9%)") + ylab("RDA 2 (19.6%)") + theme_light()

### Quite nice...
ggsave(plot = plot, filename = "RDA4_zoo_abund_v.pdf", dpi = 300, width = 8, height = 6)


### 20/06/17: Plot NBSS clusters (k=6) in RDA space ! 
groups ; length(groups)
nrow(stations_scores)
stations_scores$id <- paste(stations_scores$station, stations_scores$month, stations_scores$year, sep = "_")
stations_scores$NBSS_cluster <- NA
# Provide NBSS cluster
# i <- unique(stations_scores$id)[106]
for(i in unique(stations_scores$id) ) {
		k <- groups[i]
		stations_scores[which(stations_scores$id == i),"NBSS_cluster"] <- k
} # eo for loop
### Okay, plot now
quartz()
ggplot() + geom_point(aes(x = RDA1, y = RDA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size = 4) + 
		   	geom_point(aes(x = RDA1, y = RDA2, fill = factor(NBSS_cluster)), data = stations_scores, pch = 23, colour = "black", size = 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = RDA1*0.7, yend = RDA2*0.7), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour = "black", data = vars_scores) + 
		   	scale_fill_manual("NBSS cluster", values = c("#d53e4f","#fc8d59","#fee08b","#e6f598","#99d594","#3288bd")) + 
		   	geom_text(aes(x = RDA1, y = RDA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
			scale_x_continuous(limits = c(-0.5,0.6)) + scale_y_continuous(limits = c(-0.7,0.5)) + 
			xlab("RDA 1 (53.1%)") + ylab("RDA 2 (22.8%)") + theme_light()



### Correlations with longitude ?
colnames(cctable)
quartz() 
ggplot() + geom_point(aes(x = SiOH4, y = other_Calanoida, fill = factor(year) ), data = cctable[cctable$Chla < 15,], pch = 21, colour = "black") + 
			scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
			xlab("[SiOH4]") + ylab("Other Calanoida biovolume (mm3)") + theme_light()

# Correlogram
### Nice correlogram:
library("corrgram")
cormat <- round( cor(na.omit(cctable[cctable$Chla<15,c(4,6:40)]), method = "spearman"), 2)
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
	cormat <- cormat[hc$order, hc$order]
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

ggsave(plot = heatmap, filename = "correlogram_zoo_abund_v.pdf", dpi = 300, height = 17, width = 21)



### Idea for correlating zooplankton abundances with phytoplankton community structure: correlate abundances to CA scores (CA1 + CA2)



##### 04/05/2017: Compute total zooplanktopn abundance and biovolume and assess its correlation with the environment !
ddf <- read.table("PNMIR_cruise_data_abund_02_05_17.csv", h = T, sep = ";", dec = ",")
colnames(ddf)
str(ddf)

zoo <- na.omit( ddf[,c("month","year","station","date","Calanoida","Cladocera","Appendicularia","Hydrozoa","Cyclopoida","nauplii","Cirripedia", 
						"Thecosomata","Poecilostomatoida","Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida")] )

env <- ddf[,c(1:43)]
			
# First, melt:
mzoo <- melt(zoo, id.vars = c("station","date","month","year"))
head(mzoo)
# Provide id for computing total abundance per station :)
mzoo$id <- paste(mzoo$station, mzoo$month, mzoo$year, sep = "_")

totabund <- data.frame(mzoo %>% 
		 	group_by(id) %>% 
	 	   	summarise( total_abund = sum(value) ))

# Sort by decreasing order of total abund
totabund <- totabund[order(totabund$total_abund, decreasing = T), ]
# force factor levels to keep this order : factor() !
totabund$id <- factor(totabund$id, levels = totabund$id)
quartz()
ggplot(totabund) + geom_bar(aes(x = id, y = total_abund), stat = "identity") + 
	 		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	 		xlab("Classes") + ylab("Total zooplankton abundance (ind/m3) - all years and stations")
		
### OK, do it for biovolumes now
ddf <- read.table("PNMIR_cruise_data_biovol_03_05_17.csv", h = T, sep = ";", dec = ".")
colnames(ddf)
str(ddf)

zoo <- na.omit( ddf[,c("month","year","station","date","Calanoida","Cladocera","Appendicularia","Hydrozoa","Cyclopoida","nauplii","Cirripedia", 
						"Thecosomata","Poecilostomatoida","Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida")] )

mzoo <- melt(zoo, id.vars = c("station","date","month","year"))
head(mzoo)
# Provide id for computing total abundance per station :)
mzoo$id <- paste(mzoo$station, mzoo$month, mzoo$year, sep = "_")

totbiovol <- data.frame(mzoo %>% 
		 		group_by(id) %>% 
	 	   		summarise( total_biovol = sum(value) ))

# Sort by decreasing order of total abund
totbiovol <- totbiovol[order(totbiovol$total_biovol, decreasing = T), ]
# force factor levels to keep this order : factor() !
totbiovol$id <- factor(totbiovol$id, levels = totbiovol$id)
quartz()
ggplot(totbiovol) + geom_bar(aes(x = id, y = total_biovol), stat = "identity") + 
	 		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	 		xlab("Classes") + ylab("Total zooplankton biovolume (mm3) - all years and stations")

### Now, provide total biovolume and total abundance of the zoo community to environmental data:
env$id <- paste(env$station, env$month, env$year, sep = "_")
# Add empty columms for total biovolume and abundance
env$total_abund <- NA
env$total_biovol <- NA

matching <- na.omit(totbiovol$id[match(env$id, totbiovol$id)]) ; length(matching)
#matching <- na.omit(totabund$id[match(env$id, totabund$id)]) ; length(matching) ### Same with abundances

for(i in matching) {
		
		# Doing : 
		message(paste(i, sep = ""))
		
		# Filter biovolumes according to id
		biovol <- totbiovol[which(totbiovol$id == i),"total_biovol"]
		abundance <- totabund[which(totabund$id == i),"total_abund"]
		
		# cbind with env
		env[which(env$id == i),"total_biovol"] <- biovol
		env[which(env$id == i),"total_abund"] <- abundance
		
} # eo for loop

summary(env)

# drop ids:
env <- env[,-c(44)]

# Add Tdiff:
env$Tdiff <- env$Ts - env$Tf

### Perform PCA using total_abund and total_biovol as sup.vars
library("FactoMineR")
res.PCA <- PCA(X = env[which(env$year %in% c(2011,2012,2013,2015) & env$Chla <= 20 & env$PO4 <= 4 & env$NO3 <= 10 & env$SiOH4 <= 10),
			   c("x","max_depth","Ts","Tdiff","Ss","NO3","SiOH4","PO4","Chla","total_abund","total_biovol")], scale.unit = T, quanti.sup = c(10,11) )
				
summary(res.PCA) # 74.37% of total variance on first 4 PCs
plot(res.PCA)
quartz()
plot.PCA(res.PCA, choix = "var", axes = c(3,4))

### correlogram
library("corrgram")
cormat <- round( cor(na.omit(env[which(env$year %in% c(2011,2012,2013,2015) & env$Chla <= 20 & env$PO4 <= 4 & env$NO3 <= 10 & env$SiOH4 <= 10),
			   c("x","max_depth","Ts","Tdiff","Ss","NO3","SiOH4","PO4","Chla","total_abund","total_biovol")]), method = "spearman"), 2)
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

###
summary(lm(total_abund ~ x, data = env)) # p-value: 0.0007175 ; R-squared: 0.1138
# Slight slight decrease of zoo abundance from the shores to the open ocean
summary(lm(total_biovol ~ x, data = env)) # Nope
summary(lm(total_biovol ~ total_abund, data = env)) # R-squared: 0.5535; p-value < 2.2e-16  


##### Compute total abundance for each month 
dir()
ddf <- read.table("PNMIR_cruise_data_abund_v_21_06_17.csv", h = T, sep = ";", dec = ".")
colnames(ddf)
str(ddf)

zoo <- na.omit( ddf[,c("month","year","date","Calanoida","Cladocera","Appendicularia","Hydrozoa","Cyclopoida","nauplii","Cirripedia", 
						"Thecosomata","Poecilostomatoida","Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida")] )

env <- ddf[,c(1:43)]
			
# First, melt:
mzoo <- melt(zoo, id.vars = c("date","month","year"))
head(mzoo)
# Provide id for computing total abundance per station :)
mzoo$id <- paste(mzoo$month, mzoo$year, sep = "_")

totabund <- data.frame(mzoo %>% 
		 	group_by(id) %>% 
	 	   	summarise( total_abund = sum(value) ))

# Sort by date
totabund$date <- str_replace_all(totabund$id, "_", "-")
library("lubridate")
totabund$date <- lubridate::dmy(x = totabund$date)

quartz()
ggplot(totabund) + geom_bar(aes(x = id, y = total_abund), stat = "identity") + 
	 		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	 		xlab("Classes") + ylab("Total zooplankton abundance (ind/m3) - all years and stations")
		
### OK, do it for biovolumes now
ddf <- read.table("PNMIR_cruise_data_biovol_03_05_17.csv", h = T, sep = ";", dec = ".")
colnames(ddf)
str(ddf)

zoo <- na.omit( ddf[,c("month","year","station","date","Calanoida","Cladocera","Appendicularia","Hydrozoa","Cyclopoida","nauplii","Cirripedia", 
						"Thecosomata","Poecilostomatoida","Siphonophorae","Chaetognatha","Harpacticoida","Euphausiacea","Actinopterygii","Amphipoda","Annelida")] )

mzoo <- melt(zoo, id.vars = c("date","month","year"))
head(mzoo)
# Provide id for computing total abundance per station :)
mzoo$id <- paste(mzoo$station, mzoo$month, mzoo$year, sep = "_")

totbiovol <- data.frame(mzoo %>% 
		 		group_by(id) %>% 
	 	   		summarise( total_biovol = sum(value) ))

# Sort by decreasing order of total abund
totbiovol <- totbiovol[order(total_biovol$total_biovol, decreasing = T), ]
# force factor levels to keep this order : factor() !
totbiovol$id <- factor(totbiovol$id, levels = totbiovol$id)
quartz()
ggplot(totbiovol) + geom_bar(aes(x = id, y = total_biovol), stat = "identity") + 
	 		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
	 		xlab("Classes") + ylab("Total zooplankton biovolume (mm3) - all years and stations")


### Look at correlation between phyto comm structure and nutrients
ggplot() + geom_point(aes(x = SiOH4, y = log(n_dino), fill = factor(month) ), data = ddf, pch = 21, colour = "black") + 
			scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
			xlab("[SiOH4]") + ylab("Dinoflagellates counts (log)") + theme_light()
