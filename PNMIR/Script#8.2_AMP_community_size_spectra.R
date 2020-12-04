
##### 22/06/2017: R Script to re-analyze the zooplankton size and community structure from a clean sheet.

##### Aims to:

#	- Load the community structure (abundances and then E.biovolumes) to perform multivariate analysis.
#	- Assess which variables are the most correlated (correlogram + PCA)
#	- Compute and plot copepoda size spectra (not biovolumes!) for each sample
#	- Plus plot biovolume spectra of the large and small copepoda for each sample

### Latest update: 22/11/2017

library("dplyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("Hmisc")
library("vegan")
library("FactoMineR")
library("matrixStats")
library("corrgram")
library("lubridate")
library("viridis")

### Functions for the correlogram
get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
}
reorder_cormat <- function(cormat) {
	# Utiliser la corrélation entre les variables
	# comme mésure de distance
	dd <- as.dist((1-cormat) / 2)
	hc <- hclust(dd)
	cormat <- cormat[hc$order, hc$order]
}	

# --------------------------------------------------------------------------------------------------------------------------

##### 1°) Community analyses - RDA + correlograms on zooplankton abundances and biovolumes

################################################# A) Community analyses (CA, CCA, RDA) based on abundances 

ddf <- read.table("PNMIR_cruise_data_abund_v_07_08_17.txt", h = T, sep = "\t", dec = ".")
colnames(ddf)
str(ddf)

# Get rid of 2010:
ddf <- ddf[which(ddf$year != 2010),]
# Difference between surface and bottom temperatures
ddf$Tdiff <- ddf$Ts - ddf$Tf

# Transform some of the data (the phyto counts)
summary(ddf)
ddf[,c(25:28,40:43)] <- log10(ddf[,c(25:28,40:43)]+1)

### Correlogram
colnames(ddf)
cormat <- round( cor(na.omit(ddf[,c(7:9,12,15,17:20,23:32,48:80)]), method = "spearman"), 2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)	
# Re-order correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
			theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

### Okray, let's add some text
heatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
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
heatmap		

ggsave(plot = heatmap, filename = "correlogram_abund_v.pdf", dpi = 300, height = 17, width = 21)

### Check the SST-Chla-Diatoms relationship
summary(glm(log1p(n_diato) ~ Ts + Chla, data = ddf[which(ddf$Chla<20),]))
summary(glm((n_diato/n_phyto_total) ~ Ts + Chla , data = ddf))
#plot(glm(log1p(n_diato) ~ Ts + Chla, data = ddf[which(ddf$Chla<20),]))
ggplot() + geom_point(aes(x = Ts, y = Chla, fill = (n_diato/n_phyto_total)), data = ddf[which(ddf$Chla<20),], pch = 21, colour = "black", size = 3) + scale_fill_viridis() + theme_light()
summary(lm(ratio ~ n_diato + I(n_diato)^2, data = ddf))
quartz()
ggplot(ddf) + geom_point(aes(x = n_diato, y = ratio, fill = (n_diato/n_phyto_total)), pch = 21, colour = "black", size = 3) + 
		   geom_smooth(aes(x = n_diato, y = ratio), method = "glm", family = gaussian(lin = "log"), col = "black") + scale_fill_viridis(name = "") + 
		   theme_light() + ylab("Copepod size ratio\n(small/large)") + xlab("Diatoms counts")

### Lokk at the Chl-a variations across seasons
ddf$season <- NA
ddf[which(ddf$month %in% c(4,5,6)),"season"] <- "Spring"
ddf[which(ddf$month == 7),"season"] <- "Summer"
ddf[which(ddf$month == 10),"season"] <- "Fall"	

ggplot(ddf[ddf$Chla < 7.5,], aes(x = factor(season), y = Chla)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70")

### And ID
ddf$id <- paste(ddf$season, ddf$year, sep = "_")
ggplot(ddf[ddf$Chla < 7.5,], aes(x = factor(id), y = Chla)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70")

### Correspondence Analysis - CA
colnames(ddf)
catable <- na.omit(ddf[,c(45:47,49:51,53:61,63:71,73:76)]) 
catable <- na.omit(ddf[,c(64:71,76)]) # And when copepods only:
summary(catable)
# Perform CA
res.ca1 <- CA(catable) #
summary(res.ca1)
AFCsp <- data.frame(Ax1 = res.ca1$col$coord[,1],
                    Ax2 = res.ca1$col$coord[,2],
					Ax3 = res.ca1$col$coord[,3],
                    Cos2_1 = res.ca1$col$cos2[,1] + res.ca1$col$cos2[,2])

# For objects
AFCst <- data.frame(Ax1 = res.ca1$row$coord[,1],
                    Ax2 = res.ca1$row$coord[,2],
					Ax3 = res.ca1$row$coord[,3],
                    Cos2_1 = res.ca1$row$cos2[,1] + res.ca1$row$cos2[,2], 
					month = dat$month, station = dat$station, year = dat$year )

head(AFCst)		
	
### Change colour according to month	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()			
			

### Change colour according to year	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()	
	



# CA on phytoplankton counts
colnames(ddf)
phytoca <- CA(na.omit(ddf[,c(26,27,28)]))
# Retrieve CA1 and CA2 scores, there will represent variables of phytopl community structure
str(phytoca)
AFCst <- data.frame(Ax1 = phytoca$row$coord[,1], month = na.omit(ddf[,c(4,26:28)])$month, station = na.omit(ddf[,c(6,26:28)])$station, year = na.omit(ddf[,c(5,26:28)])$year )
head(AFCst)
AFCst$id <- paste(AFCst$station, AFCst$month, AFCst$year, sep = "_")
ddf$id <- paste(ddf$station, ddf$month, ddf$year, sep = "_")
ddf$CA1 <- NA
# Provide with for loop
for(i in unique(AFCst$id)) {
		ca1 <- AFCst[which(AFCst$id == i),"Ax1"]
		ddf[which(ddf$id == i),"CA1"] <- ca1
} # eo for loop

### Plot phytoplankton CA
AFCsp <- data.frame(Ax1 = phytoca$col$coord[,1],
                    Ax2 = phytoca$col$coord[,2],
					#Ax3 = phytoca$col$coord[,3],
                    Cos2_1 = phytoca$col$cos2[,1] + phytoca$col$cos2[,2])

# For objects
AFCst <- data.frame(Ax1 = phytoca$row$coord[,1],
                    Ax2 = phytoca$row$coord[,2],
					#Ax3 = phytoca$row$coord[,3],
                    Cos2_1 = phytoca$row$cos2[,1] + phytoca$row$cos2[,2], 
					month = na.omit(ddf[,c(4:6,26,27,29:32)])$month, station = na.omit(ddf[,c(4:6,26,27,29:32)])$station, year = na.omit(ddf[,c(4:6,26,27,29:32)])$year )

head(AFCst)		
### Add season
AFCst$season <- NA
AFCst[which(AFCst$month %in% c(4,5,6)),"season"] <- "Spring"
AFCst[which(AFCst$month == 7),"season"] <- "Summer"
AFCst[which(AFCst$month == 10),"season"] <- "Fall"	

# quartz()
plot <- ggplot() +
  	  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 21, colour = "black", size = 6, fill = "grey70")+
  		geom_point(aes(x = Ax1, y = Ax2, shape = factor(year), fill = factor(season) ), data = AFCst, alpha = 0.7, colour = "black", size = 4) +
		scale_shape_manual("Year", values = c(21,22,23,24)) + 
		scale_fill_manual("Season", values = c("#fdae61","#66c2a5","#d73027")) + 
		scale_x_continuous(paste("CA 1 (",round(phytoca$eig$per[1], 2),"%)", sep = "")) + 
		scale_y_continuous(paste("CA 2 (",round(phytoca$eig$per[2],2),"%)", sep = "")) + 
		#geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + 
		geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  		#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
		theme_light()		
		
ggsave(plot = plot, filename = "CA_phyto_07_08_17.pdf", dpi = 300, width = 8, height = 6)
			

### Change colour according to year	
#quartz()
#ggplot() +
  	#geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	#geom_point(aes(x = Ax1, y = Ax2, shape = Cos2_1, fill = factor(year) ), data = AFCst, alpha = 0.7, colour = "black") +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	#scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	#scale_x_continuous(paste("CA 1 (",round(phytoca$eig$per[1], 2),"%)", sep = "")) + 
	#scale_y_continuous(paste("CA 2 (",round(phytoca$eig$per[2],2),"%)", sep = "")) + 
	#geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	#theme_light()	


### Examine correlation between CA1 coordinates and Silica/ Nitrates...Chl-a
summary(lm(CA1 ~ SiOH4, data = ddf[ddf$SiOH4<10,])) # R-squared = 0.2366  ; p-value = 2.456e-08
ggplot(ddf[ddf$SiOH4<10,]) + geom_point(aes(x = SiOH4, y = CA1), pch = 21, colour = "black", fill = "grey70") + theme_light()
summary(lm(CA1 ~ NO3, data = ddf[ddf$NO3<10,])) #  R-squared:  0.1076  ; p-value = 0.0002032
ggplot(ddf[ddf$NO3<10,]) + geom_point(aes(x = NO3, y = CA1), pch = 21, colour = "black", fill = "grey70") + theme_light()
summary(lm(CA1 ~ Chla, data = ddf[ddf$Chla<20,])) # p-value: 7.108e-06; Rsquared = 
ggplot(ddf[ddf$Chla<20,]) + geom_point(aes(x = Chla, y = CA1), pch = 21, colour = "black", fill = "grey70") + theme_light()


### Perform RDA on abundances after hellinger transformation
### 07/08/17: Keep the following zooplankton groups
# Acartiidae, Amphipoda, Annelida, Bivalvia, Bryozoa, Calanidae, Candaciidae, Centropagidae, Chaetognatha, Cirripedia, Cladocera, Corycaeidae, Decapoda, Echino, Actino, Harpa, Hydrozoa, Oiko
# Oithonidae, Oncaeidae, Sapp, Temoridae, Thecoso, unidentif, Calanoida 
colnames(ddf)
cctable <- na.omit(ddf[,c(4:9,12,17,19,23,24,49:51,54,55,57,59:66,68:77,79,80,89:91)])
colnames(cctable)
rownames(cctable) <- as.numeric(c(1:nrow(cctable)))
# Use hellinger transformation
cctable[,c(12:37)] <- vegan::decostand(cctable[,c(12:37)], "hellinger")
# Perform RDA
res.rda <- vegan::rda( X = as.matrix(cctable[,c(12:37)]), Y = as.matrix(cctable[,c(6:11,38,40)]) )
summary(res.rda)
#                         RDA1    RDA2     RDA3     RDA4    RDA5     RDA6
# Eigenvalue            0.04542 0.01409 0.006089 0.003343 0.00209 0.001493
# Proportion Explained  0.61907 0.19200 0.082990 0.045570 0.02848 0.020350
# Cumulative Proportion 0.61907 0.81107 0.894060 0.939630 0.96811 0.988460
  
#             RDA1      RDA2     RDA3     RDA4     RDA5     RDA6
# max_depth  0.68772  0.135993 -0.53629 -0.32549 -0.20359 -0.24607
# Ts         0.07292 -0.977307  0.05101 -0.12900  0.03636  0.04128
# NO3        0.25299  0.375896  0.15593  0.01909  0.21461 -0.21495
# SiOH4      0.52065 -0.093627  0.63956  0.19615 -0.10543 -0.06010
# Chla      -0.29896 -0.008467 -0.36679  0.30685 -0.64252 -0.39033
# Phaeo_a   -0.29992  0.115032 -0.01360  0.19346 -0.22872 -0.23348
# Tdiff      0.12688 -0.357459 -0.59638 -0.32290 -0.30388  0.55057
# A1       -0.70799  0.223610  0.16090 -0.55571 -0.25910 -0.01298

anova(res.rda, step = 1000, perm.max = 1000)
### OK signficant enough
# plot(res.rda, scaling = 1)
quartz()
plot(res.rda, scaling = 2)
vegan::RsquareAdj(x = res.rda)  
# $r.squared
# [1] 0.3343674
# $adj.r.squared
# [1] 0.278898

# For ploting:
species_scores <- data.frame(species = rownames(res.rda$CCA$v), RDA1 = res.rda$CCA$v[,"RDA1"], RDA2 = res.rda$CCA$v[,"RDA2"], RDA3 = res.rda$CCA$v[,"RDA3"])
stations_scores <- data.frame(RDA1 = res.rda$CCA$u[,"RDA1"], RDA2 = res.rda$CCA$u[,"RDA2"], RDA3 = res.rda$CCA$u[,"RDA3"])
stations_scores$month <- cctable$month
stations_scores$year <- cctable$year
stations_scores$station <- cctable$station
stations_scores$season <- NA
stations_scores[which(stations_scores$month %in% c(4,5,6)),"season"] <- "Spring"
stations_scores[which(stations_scores$month == 7),"season"] <- "Summer"
stations_scores[which(stations_scores$month == 10),"season"] <- "Fall"
vars_scores <- data.frame(var = rownames(res.rda$CCA$biplot), RDA1 = res.rda$CCA$biplot[,"RDA1"], RDA2 = res.rda$CCA$biplot[,"RDA2"], RDA3 = res.rda$CCA$biplot[,"RDA3"] )

# For choosing the limits of the axes
summary(species_scores) ; summary(vars_scores)

quartz()
plot <- ggplot() + geom_point(aes(x = RDA1, y = RDA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size= 4) + 
		   	geom_point(aes(x = RDA1*1.5, y = RDA2*1.5, fill = factor(season)), data = stations_scores, pch = 23, colour = "black", size= 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = RDA1*0.7, yend = RDA2*0.7), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "black", data= vars_scores) + 
		   	scale_fill_manual("Season", values = c("#fdae61","#66c2a5","#d73027")) + 
		   	#geom_text(aes(x = RDA1, y = RDA2+0.02, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
			scale_x_continuous(limits = c(-0.6,0.6)) + scale_y_continuous(limits = c(-0.7,0.5)) + 
			xlab("RDA 1 (61.9%)") + ylab("RDA 2 (19.2%)") + theme_light()

# Save
ggsave(plot = plot, filename = "RDA_zoo_abund_v_with_phyto_07_08_17.pdf", dpi = 300, width = 8, height = 6)




### And now, without phytoplankton community
cctable <- na.omit(ddf[,c(4:9,12,17,19,23,24,45:47,49:51,53:77)])
colnames(cctable)
# Use hellinger transformation
cctable[,c(12:41)] <- vegan::decostand(cctable[,c(12:41)], "hellinger")
res.rda <- vegan::rda( X = as.matrix(cctable[,c(12:41)]), Y = as.matrix(cctable[,c(6:11,42)]) )
summary(res.rda)
#                         RDA1    RDA2     RDA3     RDA4     RDA5      RDA6
# Eigenvalue            0.04345 0.01501 0.006078 0.002288 0.001559 0.0007295
# Proportion Explained  0.62673 0.21651 0.087670 0.033010 0.022490 0.0105200
# Cumulative Proportion 0.62673 0.84324 0.930910 0.963920 0.986420 0.9969400

#               RDA1      RDA2     RDA3    RDA4     RDA5      RDA6
# max_depth  0.68527  0.135617  0.67853 -0.0777 -0.15900  0.142423
# Ts         0.08054 -0.982366 -0.01384  0.1099  0.01787  0.006241
# NO3        0.25655  0.396846 -0.12097  0.1939 -0.26246 -0.808046
# SiOH4      0.56527 -0.086972 -0.62186 -0.1117 -0.18022 -0.491018
# Chla      -0.31191  0.004295  0.30801 -0.7202 -0.44978 -0.260008
# Phaeo_a   -0.29494  0.130495 -0.02851 -0.3056 -0.24396 -0.159179
# Tdiff      0.14667 -0.379039  0.63288 -0.2610  0.60287 -0.035218
anova(res.rda, step = 1000, perm.max = 1000)
vegan::RsquareAdj(x = res.rda)  
# $r.squared
# [1] 0.2973789
# $adj.r.squared
# [1] 0.2471916
# For ploting:
quartz()
plot(res.rda, scaling = 2)

species_scores <- data.frame(species = rownames(res.rda$CCA$v), RDA1 = res.rda$CCA$v[,"RDA1"], RDA2 = res.rda$CCA$v[,"RDA2"], RDA3 = res.rda$CCA$v[,"RDA3"])
stations_scores <- data.frame(RDA1 = res.rda$CCA$u[,"RDA1"], RDA2 = res.rda$CCA$u[,"RDA2"], RDA3 = res.rda$CCA$u[,"RDA3"])
stations_scores$month <- cctable$month
stations_scores$year <- cctable$year
stations_scores$station <- cctable$station
vars_scores <- data.frame(var = rownames(res.rda$CCA$biplot), RDA1 = res.rda$CCA$biplot[,"RDA1"], RDA2 = res.rda$CCA$biplot[,"RDA2"], RDA3 = res.rda$CCA$biplot[,"RDA3"] )

# For choosing the limits of the axes
summary(species_scores) ; summary(vars_scores)
quartz()
plot <- ggplot() + geom_point(aes(x = RDA1, y = RDA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size= 4) + 
		   	geom_point(aes(x = RDA1, y = RDA2, fill = factor(month)), data = stations_scores, pch = 23, colour = "black", size= 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = RDA1*0.7, yend = RDA2*0.7), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "black", data= vars_scores) + 
		   	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		   	geom_text(aes(x = RDA1, y = RDA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
			scale_x_continuous(limits = c(-0.6,0.6)) + scale_y_continuous(limits = c(-0.7,0.5)) + 
			xlab("RDA 1 (63%)") + ylab("RDA 2 (22%)") + theme_light()




################################################# B) Community analyses (CA, CCA, RDA) based on ellipsoïdal biovolumes

ddf <- read.table("PNMIR_cruise_data_biovol_v_21_06_17.csv", h = T, sep = ";", dec = ".")
colnames(ddf)
str(ddf)

# Get rid of 2010:
ddf <- ddf[which(ddf$year != 2010),]
# Difference between surface and bottom temperatures
ddf$Tdiff <- ddf$Ts - ddf$Tf

# Transform some of the data
summary(ddf)
ddf[,c(25:28,40:43)] <- log10(ddf[,c(25:28,40:43)]+1)

### Correlogram
cormat <- round( cor(na.omit(ddf[,c(7:9,12,15,17:20,23:28,40:77)]), method = "spearman"), 2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)	
# Re-order correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
			theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

### Okray, let's add some text
heatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
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
heatmap		

ggsave(plot = heatmap, filename = "correlogram_biovol_v.pdf", dpi = 300, height = 17, width = 21)


### Correspondence Analysis - CA
colnames(ddf)
catable <- na.omit(ddf[,c(45:47,49:51,53:61,63:76)])
summary(catable)
# Perform CA
res.ca1 <- CA(catable) #

# CA on phytoplankton counts
phytoca <- 
# Retrieve CA1 and CA2 scores, there will represent variables of phytopl community structure
ddf$CA1 <- NA
ddf$CA2 <- NA

### Perform RDA on abundances after hellinger transformation
cctable <- ddf[]
cctable[,c(14:40)] <- vegan::decostand(cctable[,c(14:40)], "hellinger")

res.rda <- vegan::rda( X = as.matrix(cctable[,c(14:40)]), Y = as.matrix(cctable[,c(6:12)]) )
summary(res.rda)
anova(res.rda, step = 1000, perm.max = 1000)
plot(res.rda, scaling = 1)
plot(res.rda, scaling = 2)
vegan::RsquareAdj(x = res.rda)  

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
plot <- ggplot() + geom_point(aes(x = RDA1, y = RDA2), data = species_scores, pch = 21, fill = "grey70", colour = "black", size= 4) + 
		   	geom_point(aes(x = RDA1, y = RDA2, fill = factor(year)), data = stations_scores, pch = 23, colour = "black", size= 3) + 
		   	geom_segment(aes(x = 0, y = 0, xend = RDA1*0.7, yend = RDA2*0.7), arrow = arrow(length= unit(0.175, "cm"), type= "closed"), colour= "black", data= vars_scores) + 
		   	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
		   	geom_text(aes(x = RDA1, y = RDA2+0.05, label = rownames(species_scores)), data = species_scores, size = 3) + 
			geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
			scale_x_continuous(limits = c(-0.5,0.6)) + scale_y_continuous(limits = c(-0.7,0.5)) + 
			xlab("RDA 1 (%)") + ylab("RDA 2 (%)") + theme_light()

# Save
ggsave(plot = plot, filename = "RDA_zoo_biovol_v.pdf", dpi = 300, width = 8, height = 6)



ggplot() + geom_point(aes(x = x, y = log(Calanoida+Cladocera+Oikopleuridae+Hydrozoa+Cyclopoida+nauplii+Cirripedia+Chaetognatha+Euphausiacea+Annelida+Amphipoda+egg+Poecilostomatoida) ), data = ddf) + theme_classic()

# --------------------------------------------------------------------------------------------------------------------------



##### 2°) Copepoda size spectra - define a fixed spectrum based on ALL Copepoda - plot it by highlighting the 1mm threshold

################################################# But first, need to load and concatenate all the data...
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/Data_06_07_17/")
#proj_dir <- getwd()
#projects <- dir()
#projects
#p <- "export_434_20170706_0859"

# Prior to lapply: use one project to save the colnames of interest: object_id, object_annotation + morphometric traits
#setwd(paste(proj_dir,"/",p,"/", sep = ""))
#all_samples <- dir()
# Read and rbind the samples from each project
#samples <- lapply(all_samples, function(s) {	
			# Load and return
			#ss <- read.table(file = s, sep = '\t', header = T)
			#return(ss)	
#}) # eo 2nd lapply 
#all_samples <- do.call(rbind, samples)
# dim(all_samples)
# colnames(all_samples)
### Choose the columns of interest. But let's keep all for now.
#cols <- colnames(all_samples)[c(1:119,121:159)] # because of "process_particle_pixel_size__m" in "export_434_20170608_1544"
#rm(all_samples, samples)

### Now, load all data for real
#projects
# p <- "export_442_20170608_1539"
#projs <- lapply(projects, function(p) {
				# Go to project dir
				#message(paste("Doing project ", p, sep = ""))
				#setwd(paste(proj_dir,"/",p,"/", sep = ""))
				#all_samples <- dir()
				# Read and rbind the samples from each project
				#samples <- lapply(all_samples, function(s) {
						# Load and return
						#ss <- read.table(file = s, sep = '\t', header = T)
						#return(ss)	
				#}) # eo 2nd lapply 
				#all_samples <- do.call(rbind, samples)
				# dim(all_samples) ; colnames(all_samples); summary(all_samples)
				#rm(samples)
				# Return, by restricting to columns of interest
				#return(all_samples[,cols])
#}) # eo lapply
#ddf <- do.call(rbind, projs)
# str(ddf)
# colnames(ddf)
### Select validated vignettes only.
#ddf2 <- ddf[which(ddf$object_annotation_status == "validated"),] # Validated vignettes only
### Get rid of horizontal hauls...
#horizontal <- ddf2$acq_id[grep("_h", ddf2$acq_id)]
#ddf2 <- ddf2[!(ddf2$acq_id %in% horizontal),]
### Get rid of non-living particles and "temporary"
#nl <- ddf2$object_annotation_hierarchy[grep("not-living", ddf2$object_annotation_hierarchy)]
#ddf2 <- ddf2[!(ddf2$object_annotation_hierarchy %in% nl),]
#temp <- ddf2$object_annotation_hierarchy[grep("temporary", ddf2$object_annotation_hierarchy)]
#ddf2 <- ddf2[!(ddf2$object_annotation_hierarchy %in% temp),]
#molene <- ddf2$object_id[grep("molene", ddf2$object_id)]
#Molene <- ddf2$object_id[grep("Molene", ddf2$object_id)]
#douar <- ddf2$object_id[grep("douarnenez", ddf2$object_id)]
#ddf2 <- ddf2[!(ddf2$object_id %in% molene),]
#ddf2 <- ddf2[!(ddf2$object_id %in% douar),]
#ddf2 <- ddf2[!(ddf2$object_id %in% Molene),]
#rm(douar, Molene, molene)
# Replace caps
#ddf2$object_id <- str_replace_all(as.character(ddf2$object_id), "D", "d")
#ddf2$object_id <- str_replace_all(as.character(ddf2$object_id), "B", "b")
### And remove 2010 samples because their id do not follow the same formating
#s2010 <- ddf2$object_id[grep("wp2_2010", ddf2$object_id)]
#ddf2 <- ddf2[!(ddf2$object_id %in% s2010),]
### Add ids
#ID <- data.frame( do.call(rbind, strsplit(x = as.character(ddf2$object_id), split = "_")) )
#colnames(ID)[1:4] <- c("net", "station", "date", "type")
#ddf2$station <- ID$station
#ID$date <- lubridate::ymd(ID$date)
#ddf2$date <- ID$date
#ddf2$month <- lubridate::month(ddf2$date)
#ddf2$year <- lubridate::year(ddf2$date)
#ddf2$station_id <- paste(ddf2$station,ddf2$month,ddf2$year, sep = "_")
#gc()

# Extract the classification thanks to ddf2$object_annotation_hierarchy and select what is relevant.
#split_status <- strsplit(x = as.character(ddf2$object_annotation_hierarchy), split = ">") 
#length(split_status)
# split_status[[276312]][length(split_status[[276312]])] ### here's how you can do it
# Retrieve status and a few parent categories, provide back
#categories <- lapply(c(1:length(split_status)), function(i) {
					#message(paste(i))
					#category <- split_status[[i]] [length(split_status[[i]])] 
					# Parent categories
					#parent1 <- split_status[[i]] [length(split_status[[i]]) - 1] 
					#parent2 <- split_status[[i]] [length(split_status[[i]]) - 2] 				
					# Return
					#return(data.frame(i = i, category = category, parent1 = parent1, parent2 = parent2 ))
#}) # eo lapply 
# rbind as usual
#cat <- do.call(rbind, categories)	 
#dim(cat)
# cbind()
#ddf3 <- cbind(ddf2, cat[,c(2:4)])
#rm(cat, categories, temp, nl, horizontal, ddf, split_status)

### Save the data so you don't have to go through the previous code all over again everytime
#save(ddf3, file = "ddf3_06_07_17.Rdata")


############################################################################ Load the data
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
ddf3 <- get(load("ddf3_06_07_17.Rdata"))
# Compute size from major axis
ddf3$size <- ddf3$object_major*0.0106
# Compyte abundance/concentration
ddf3$concentration <- 1 * ddf3$acq_sub_part / ddf3$sample_tot_vol
# Add ellipsoidal and spherical biovolumes
ddf3$EBioVol <- ddf3$concentration * (4/3 * pi * (ddf3$object_major * 0.0106 / 2) * (ddf3$object_minor * 0.0106 / 2)^2)
ddf3$SBioVol <- ddf3$concentration * (4/3 * pi *(sqrt((ddf3$object_area*(0.0106)^2)/pi))^3)
# summary(log10(ddf3$EBioVol))
ddf3$logEBioVol <- log10(ddf3$EBioVol)
ddf3$logSBioVol <- log10(ddf3$SBioVol)

### Check link between the 2 computed biovolumes
# summary(lm(logEBioVol ~ logSBioVol, data = ddf3))
# summary(lm(EBioVol ~ SBioVol, data = ddf3))
### They are basically the same !
# summary(lm(EBioVol ~ SBioVol, data = ddf3))


##### 06/07/17: After Van der Lingen et al. (2006) : look at the seasonal/monthly/annual variations of the % of :
#	- Cyclopoids
#	- small unidentified Calanoida (<1mm)
#	- larger unidentified Calanoida and Calanidae
ddf3$season <- NA
ddf3[which(ddf3$month %in% c(4,5,6)),"season"] <- "spring"
ddf3[which(ddf3$month == 7),"season"] <- "summer"
ddf3[which(ddf3$month == 10),"season"] <- "fall"
ddf3$id <- paste(ddf3$season, ddf3$year, sep = "_" )

###### When considering all major mesozooplankton groups
#ddf3 <- ddf3[which(ddf3$category %in% c("Acartiidae","Calanidae","Calanoida","Centropagidae","Euchaetidae","Evadne","nauplii","Oikopleuridae","Oithonidae","Oncaeidae",
										#"Temoridae","zoea","Decapoda","Hydrozoa","Limacinidae","Podon","Bivalvia","Chaetognatha","Noctiluca","Ophiuroidea","Ostracoda",
										#"Annelida","Obelia","Amphipoda","Corycaeidae","Doliolida","Harpacticoida","Siphonophorae","Actinopterygii","Monstrilloida","Hyperiidea",
										#"Penilia","Candaciidae","Gammaridea","Echinodermata","Echinoidea","Ophiurida","Salpida","Copepoda","Leptothecata","Pontellidae","Scyphozoa",
										#"Poecilostomatoida","Ctenophora","Sapphirinidae","Bryozoa","Oncaea","Creseidae","Appendicularia","Ophiothrix")),] 

### Select only the Copepoda categories ! 
# unique(ddf3$category)
ddf3 <- ddf3[which(ddf3$category %in% c("Acartiidae","Calanidae","Calanoida","Centropagidae","Euchaetidae",
										"Oithonidae","Oncaeidae","Temoridae","Corycaeidae","Harpacticoida",
										"Candaciidae","Pontellidae","Sapphirinidae","Oncaea")),] 

### When using copepods only: use the relationships from Lehett & Hernandez-Leon (2009), for subtropical copepods, to co,vert the biovolume/area to dry wright !
### 	DW (µg) = a.(S)^b
#	With: 	- a = 45.25
#			- S = body area in mm2
#			- b = 1.59

ddf3$DW <- 45.25*(ddf3$object_area*(0.0106)^2)^1.59
ddf3$DW <- ddf3$DW * ddf3$concentration
summary(ddf3$DW)
# summary(lm(DW ~ logEBioVol, data = ddf3))

### Plot variations of dry weights -  need to change the order of the bars though
# ggplot(ddf3, aes(x = factor(id), y = DW) + theme_light() + geom_bar(width = 0.5, fill = "grey70") + xlab("Season") + ylab("Dry weight (µm)")
positions <- c("spring_2011", "summer_2011", "fall_2011", "spring_2013", "summer_2013", "fall_2013", "spring_2015", "summer_2015", "fall_2015")
# quartz() 
plot <- ggplot(ddf3) + geom_bar(aes(x = factor(id), y = sum(DW)/length(unique(station_id))), stat = "identity") + 
		theme_light() + xlab("Season") + ylab("Total copepod dry weight (µg)") + scale_x_discrete(limits = positions)
		
ggsave(plot = plot, filename = "hist_norm_tot_DW.pdf", dpi = 300, width = 11, height = 7)


quartz() 
ggplot(ddf3) + geom_bar(aes(x = factor(id), y = mean(size) ), stat = "identity") + theme_light() + xlab("Season") + ylab("Average copepod size(mm)") + scale_x_discrete(limits = positions)

avg_DW <- mean(ddf3$DW)
quartz()
ggplot(ddf3, aes(x = factor(id), y = DW - avg_DW)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70") + 
xlab("Season") + ylab("Copepod dry weight anomaly (µm)") + scale_x_discrete(limits = positions) + scale_y_continuous(limits = c(0,500))


##### 25/07/17: Plot the copepod size ratio (from abund and DW)) per seasons (ddf3$id)
ddf3$ratio <- NA
ddf3$small_abund <- NA
ddf3$large_abund <- NA
for(i in unique(ddf3$station_id)) {
		small_abund <- sum(ddf3[which(ddf3$station_id == i & ddf3$size <= 1),"concentration"])
		large_abund <- sum(ddf3[which(ddf3$station_id == i & ddf3$size > 1),"concentration"])
		ratio <- small_abund/large_abund
		# Provide:
		ddf3[which(ddf3$station_id == i),"ratio"] <- ratio
		ddf3[which(ddf3$station_id == i),"small_abund"] <- small_abund
		ddf3[which(ddf3$station_id == i),"large_abund"] <- large_abund
} # eo for loop

### Plot the ratio
# quartz() 
plot <- ggplot(ddf3, aes(x = factor(id), y = ratio)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70") + 
		xlab("") + ylab("Copepod size ratio (abundances)") + scale_x_discrete(limits = positions)

ggsave(plot = plot, filename = "boxplot_size_ratio_abund_v2.pdf", dpi = 300, height = 5, width = 10)



############################################################################  For each 'ddf3$category' & 'ddf3$id', compute the total bv and the sum of anomalies (Peretti et al. 2017 - MEPS)
# i <- unique(ddf3$id)[1]
### For each id of ddf3, compute the sum of biovolumes

### 17/07/2017: BUT need to normalize by number of samples for each season !
res <- lapply(unique(ddf3$id), function(i) {
				
				# Useless message
				message(paste(i, sep = ""))
				d <- ddf3[which(ddf3$id == i),]
				# 2nd lapply, per categories (you will adD THE TOTAL BV, the small BV and THE large BV later (regardless of the taxonomic classif)
				biovolumes <- lapply(unique(d$category), function(sp) {
									bv <- sum(d[which(d$category == sp),"concentration"])
									return(data.frame(sp = sp, bv = bv))
				}) # eo 2nd lapply
				# Rbind biovolumes
				biovols <- do.call(rbind, biovolumes)
				biovols$id <- i
				# Add total biovolume
				biovols$tot_bv <- sum(d$EBioVol)
				# And small and large cops biovolumes 
				biovols$small_bv <- sum(d[which(d$size <= 1),"concentration"])
				biovols$large_bv <- sum(d[which(d$size > 1),"concentration"])
				# Add nuber of samples
				ss <- length( unique(d[which(d$id == i),"station_id"]) )
				biovols$effort <- ss
				# Return this
				return(biovols)
		} 
		
) # eo first lapply
t <- do.call(rbind,res)
t
gc() ; rm(res)

### Normalize by sampling effort (number of samples)
t[,c("tot_bv","small_bv","large_bv")] <- t[,c("tot_bv","small_bv","large_bv")]/t[,c("effort")]

### Average small_bv and large_bv for anomalies
avg_totbv_small <- mean(t$small_bv)
avg_totbv_large	<- mean(t$large_bv)
avg_totbv <- mean(t$tot_bv)

### To compute the anomalies, need the average of each category as well:
avg_bv_taxa <- data.frame(category = unique(ddf3$category), avg_bv = NA)
# Fill with for loop:
for(cat in avg_bv_taxa$category) {
		avg <- mean(t[which(t$sp == cat),"bv"])
		avg_bv_taxa[avg_bv_taxa$category == cat,"avg_bv"] <- avg
}
# Check
avg_bv_taxa

### OKray, now, try to add the proper anomalies to each row of 't' by using the objects 'avg_bv_taxa' and 'avg_totbv_small' and 'avg_totbv_large'
t$anomaly <- NA
t$anom_small <- NA
t$anom_large <- NA
for(r in 1:nrow(t) ) {
		# For each sp
		cat <- t[r,"sp"]
		avg <- avg_bv_taxa[which(avg_bv_taxa$category == cat),"avg_bv"]
		t[r,"anomaly"] <- t[r,"bv"] - avg
		# For the two size claases now
		t[r,"anom_small"] <- t[r,"small_bv"] - avg_totbv_small
		t[r,"anom_large"] <- t[r,"large_bv"] - avg_totbv_large
} # eo for loop

# positions <- c("spring_2011", "summer_2011", "fall_2011", "spring_2012", "spring_2013", "summer_2013", "fall_2013", "spring_2015", "summer_2015", "fall_2015")
### 22/11/17: V2 for article: no spring 2012
positions <- c("spring_2011", "summer_2011", "fall_2011", "spring_2013", "summer_2013", "fall_2013", "spring_2015", "summer_2015", "fall_2015")

t$diff <- t$anom_small - t$anom_large

# quartz() 
plot <- ggplot(t) + geom_bar(aes(x = factor(id), y = anomaly, fill = sp), stat = "identity", colour = "black") + 
			scale_fill_manual(name = "Category", values = viridis(length(unique(t$sp)))) + 
			scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Biovolume anomaly") + xlab("Season")
			
### Good job ! Just need to re-order the ids/ dates ;-)
ggsave(plot = plot, filename = "categories_bv_anomalies.pdf", dpi = 300, width = 11, height = 7)

### For large and small anomalies now
# quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anom_small), stat = "identity", fill = "#d53e4f", data = t[t$anom_small > 0,], colour = "black") + 
			geom_bar(aes(x = factor(id), y = anom_small), stat = "identity", fill = "#3288bd", data = t[t$anom_small < 0,], colour = "black") + 
			scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Biovolume anomaly for small Copepoda") + xlab("Season")
			
#
ggsave(plot = plot, filename = "small_cop_bv_anomalies.pdf", dpi = 300, width = 11, height = 7)
# quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anom_large), stat = "identity", fill = "#d53e4f", data = t[t$anom_large > 0,], colour = "black") + 
			geom_bar(aes(x = factor(id), y = anom_large), stat = "identity", fill = "#3288bd", data = t[t$anom_large < 0,], colour = "black") + 
			scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Biovolume anomaly for large Copepoda") + xlab("Season")
								
ggsave(plot = plot, filename = "large_cop_bv_anomalies.pdf", dpi = 300, width = 11, height = 7)

### Now, the difference between the small and the large anomalies
# quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = diff), stat = "identity", fill = "#d53e4f", data = t[t$diff > 0,], width = .8, colour = "black") + 
			geom_bar(aes(x = factor(id), y = diff), stat = "identity", fill = "#3288bd", data = t[t$diff < 0,], width = .8, colour = "black") + 
			scale_x_discrete(limits = positions) + theme_light() + xlab("Season") + geom_hline(yintercept = 0) +
			ylab("Copepod size index\n(small copepod biovolume anomaly minus large copepod biovolume anomaly)")
			
ggsave(plot = plot, filename = "small&large_cop_bv_anomalies.pdf", dpi = 300, width = 11, height = 7)




##### And with abundances instead of biovolume ?
res <- lapply(unique(ddf3$id), function(i) {
				
				# Useless message
				message(paste(i, sep = ""))
				d <- ddf3[which(ddf3$id == i),]
				# 2nd lapply, per categories (you will adD THE TOTAL BV, the small BV and THE large BV later (regardless of the taxonomic classif)
				abundances <- lapply(unique(d$category), function(sp) {
									abund <- sum(d[which(d$category == sp),"concentration"])
									return(data.frame(sp = sp, abund = abund))
				}) # eo 2nd lapply
				# Rbind biovolumes
				abund <- do.call(rbind, abundances)
				abund$id <- i
				# Add total biovolume
				abund$tot_abund <- sum(d$EBioVol)
				# And small and large cops biovolumes 
				abund$small_abund <- sum(d[which(d$size <= 1),"concentration"])
				abund$large_abund <- sum(d[which(d$size > 1),"concentration"])
				# Add nuber of samples
				ss <- length( unique(d[which(d$id == i),"station_id"]) )
				abund$effort <- ss
				# Return this
				return(abund)
		} 
		
) # eo first lapply
t <- do.call(rbind,res)
t
gc() ; rm(res)

### Divide by sampling effort (number of samples)
t[,c("tot_abund","small_abund","large_abund")] <- t[,c("tot_abund","small_abund","large_abund")]/t[,c("effort")]

### Average small_bv and large_bv for anomalies
avg_totabund_small <- mean(t$small_abund)
avg_totabund_large	<- mean(t$large_abund)
avg_totabund <- mean(t$tot_abund)

### To compute the anomalies, need the average of each category as well:
avg_bv_taxa <- data.frame(category = unique(ddf3$category), avg_abund = NA)
# Fill with for loop:
for(cat in avg_bv_taxa$category) {
		avg <- mean(t[which(t$sp == cat),"abund"])
		avg_bv_taxa[avg_bv_taxa$category == cat,"avg_abund"] <- avg
}
# Check
avg_bv_taxa

### OKray, now, try to add the proper anomalies to each row of 't' by using the objects 'avg_bv_taxa' and 'avg_totbv_small' and 'avg_totbv_large'
t$anomaly <- NA
t$anom_small <- NA
t$anom_large <- NA
for(r in 1:nrow(t) ) {
		# For each sp
		cat <- t[r,"sp"]
		avg <- avg_bv_taxa[which(avg_bv_taxa$category == cat),"avg_abund"]
		t[r,"anomaly"] <- t[r,"abund"] - avg
		# For the two size claases now
		t[r,"anom_small"] <- t[r,"small_abund"] - avg_totabund_small
		t[r,"anom_large"] <- t[r,"large_abund"] - avg_totabund_large
} # eo for loop

t$diff <- t$anom_small - t$anom_large

#quartz() 
plot <- ggplot(t) + geom_bar(aes(x = factor(id), y = anomaly, fill = sp), stat = "identity", colour = "black") + scale_x_discrete(limits = positions) + 
		scale_fill_manual(name = "Category", values = viridis(length(unique(t$sp)))) + 
		geom_hline(yintercept = 0) + theme_light() + ylab("Abundance anomaly") + xlab("Season")
		
ggsave(plot = plot, filename = "categories_abund_anomalies.pdf", dpi = 300, width = 11, height = 7)
### Good job ! Just need to re-order the ids/ dates ;-)
### For large and small anomalies now
#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anom_small), stat = "identity", fill = "#d53e4f", data = t[t$anom_small > 0,], colour = "black") + 
		geom_bar(aes(x = factor(id), y = anom_small), stat = "identity", fill = "#3288bd", data = t[t$anom_small < 0,], colour = "black") + 
		scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Abundance anomaly for small Copepoda") + xlab("Season")
		
ggsave(plot = plot, filename = "small_cop_abund_anomalies.pdf", dpi = 300, width = 11, height = 7)	

#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anom_large), stat = "identity", fill = "#d53e4f", data = t[t$anom_large > 0,], colour = "black") + 
		geom_bar(aes(x = factor(id), y = anom_large), stat = "identity", fill = "#3288bd", data = t[t$anom_large < 0,], colour = "black") + 
		scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Abundance anomaly for large Copepoda") + xlab("Season")
ggsave(plot = plot, filename = "large_cop_abund_anomalies.pdf", dpi = 300, width = 11, height = 7)		

### Now, the difference between the small and the large anomalies
#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = diff), stat = "identity", fill = "#d53e4f", data = t[t$diff > 0,], colour = "black") + 
		geom_bar(aes(x = factor(id), y = diff), stat = "identity", fill = "#3288bd", data = t[t$diff < 0,], colour = "black") + 
		scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + 
		ylab("Copepod size index\n(small copepod abundance anomaly − large copepod abundance anomaly)") + xlab("Season")
		
ggsave(plot = plot, filename = "small&large_cop_abund_anomalies.pdf", dpi = 300, width = 11, height = 7)




##### And with dry weights instead of biovolume ?
res <- lapply(unique(ddf3$id), function(i) {
				
			# Useless message
			message(paste(i, sep = ""))
			d <- ddf3[which(ddf3$id == i),]
			# 2nd lapply, per categories (you will adD THE TOTAL BV, the small BV and THE large BV later (regardless of the taxonomic classif)
			weights <- lapply(unique(d$category), function(sp) {
										wei <- sum(d[which(d$category == sp),"DW"])
										return(data.frame(sp = sp, weight = wei))
			}) # eo 2nd lapply
			# Rbind biovolumes
			weight <- do.call(rbind, weights)
			weight$id <- i
			# Add total biovolume
			weight$tot_dw <- sum(d$DW)
			# And small and large cops biovolumes 
			weight$small_dw <- sum(d[which(d$size <= 1),"DW"])
			weight$large_dw <- sum(d[which(d$size > 1),"DW"])
			# Add number of samples
			ss <- length( unique(d[which(d$id == i),"station_id"]) )
			weight$effort <- ss
			# Return this
			return(weight)
		} 
		
) # eo first lapply
t <- do.call(rbind,res)
t
gc() ; rm(res)

### Divide by sampling effort (number of samples)
t[,c("tot_dw","small_dw","large_dw")] <- t[,c("tot_dw","small_dw","large_dw")]/t[,c("effort")]


### Average small_bv and large_bv for anomalies
avg_totdw_small <- mean(t$small_dw)
avg_totdw_large	<- mean(t$large_dw)
avg_totdw <- mean(t$tot_dw)

### To compute the anomalies, need the average of each category as well:
avg_bv_taxa <- data.frame(category = unique(ddf3$category), avg_dw = NA)
# Fill with for loop:
for(cat in avg_bv_taxa$category) {
		avg <- mean(t[which(t$sp == cat),"weight"])
		avg_bv_taxa[avg_bv_taxa$category == cat,"avg_dw"] <- avg
}
# Check
avg_bv_taxa

### OKray, now, try to add the proper anomalies to each row of 't' by using the objects 'avg_bv_taxa' and 'avg_totbv_small' and 'avg_totbv_large'
t$anomaly <- NA
t$anom_small <- NA
t$anom_large <- NA
for(r in 1:nrow(t) ) {
		# For each sp
		cat <- t[r,"sp"]
		avg <- avg_bv_taxa[which(avg_bv_taxa$category == cat),"avg_dw"]
		t[r,"anomaly"] <- t[r,"weight"] - avg
		# For the two size claases now
		t[r,"anom_small"] <- t[r,"small_dw"] - avg_totdw_small
		t[r,"anom_large"] <- t[r,"large_dw"] - avg_totdw_large
} # eo for loop

t$diff <- t$anom_small - t$anom_large

#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anomaly, fill = sp), stat = "identity", data = t) + 
			scale_fill_manual(name = "Category", values = viridis(length(unique(t$sp)))) + scale_x_discrete(limits = positions) + 
			geom_hline(yintercept = 0) + theme_light() + ylab("Dry weight anomaly (Mg)") + xlab("Season")
		
### Good job ! Just need to re-order the ids/ dates ;-)
ggsave(plot = plot, filename = "categories_dw_anomalies.pdf", dpi = 300, width = 11, height = 7)
### For large and small anomalies now
#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anom_small), stat = "identity", fill = "#d53e4f", data = t[t$anom_small > 0,]) + 
			geom_bar(aes(x = factor(id), y = anom_small), stat = "identity", fill = "#3288bd", data = t[t$anom_small < 0,]) + 
			scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Dry weight anomaly for small Copepoda") + xlab("Season")
		
ggsave(plot = plot, filename = "small_cop_dw_anomalies.pdf", dpi = 300, width = 11, height = 7)	

#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = anom_large), stat = "identity", fill = "#d53e4f", data = t[t$anom_large > 0,]) + 
			geom_bar(aes(x = factor(id), y = anom_large), stat = "identity", fill = "#3288bd", data = t[t$anom_large < 0,]) + 
			scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + ylab("Dry weight anomaly for large Copepoda") + xlab("Season")
		
ggsave(plot = plot, filename = "large_cop_dw_anomalies.pdf", dpi = 300, width = 11, height = 7)		

### Now, the difference between the small and the large anomalies
#quartz() 
plot <- ggplot() + geom_bar(aes(x = factor(id), y = diff), stat = "identity", fill = "#d53e4f", data = t[t$diff > 0,]) + 
			geom_bar(aes(x = factor(id), y = diff), stat = "identity", fill = "#3288bd", data = t[t$diff < 0,]) + 
			scale_x_discrete(limits = positions) + geom_hline(yintercept = 0) + theme_light() + 
			ylab("Copepod size index\n(small copepod dry weight anomaly − large copepod dry weight anomaly)") + xlab("Season")
				
ggsave(plot = plot, filename = "small&large_cop_dw_anomalies.pdf", dpi = 300, width = 11, height = 7)

		
		
##### 07/07/17: And the size ratio anomalies too !
res <- lapply(unique(ddf3$id), function(i) {
				
			# Useless message
			message(paste(i, sep = ""))
			d <- ddf3[which(ddf3$id == i),]
			# 2nd lapply, per categories (you will adD THE TOTAL BV, the small BV and THE large BV later (regardless of the taxonomic classif)
			ratios <- data.frame(id = i)
			# Add total biovolume
			ratios$tot_bv <- sum(d[,"concentration"])
			# And small and large cops biovolumes 
			ratios$small_bv <- sum(d[which(d$size <= 1),"concentration"])
			ratios$large_bv <- sum(d[which(d$size > 1),"concentration"])
			ratios$ratio <- ratios$small_bv / ratios$large_bv
			# Return this
			return(ratios)
		} 
		
) # eo first lapply
t <- do.call(rbind,res)
t
gc() ; rm(res)
		
# quartz()
plot <- ggplot(t) + geom_bar(aes(x = factor(id), y = ratio), stat = "identity", fill = "grey50", colour = "black", width = .8) + scale_x_discrete(limits = positions) + theme_light() + 
			ylab("Copepod size ratio\n(Total abundance of small Copepoda/Total abundance of large Copepoda)") + geom_hline(yintercept = 0) + xlab("Season")
		
ggsave(plot = plot, filename = "hist_ratio_abund.pdf", dpi = 300, width = 10, height = 7)				
		


##### 07/07/17: Test the significance of the seasonal impact on the size ratios (biovol/ abund/ dry weight) 
res <- lapply(unique(ddf3$sample_id), function(i) {
				
			# Useless message
			message(paste(i, sep = ""))
			d <- ddf3[which(ddf3$sample_id == i),]
			# 2nd lapply, per categories (you will adD THE TOTAL BV, the small BV and THE large BV later (regardless of the taxonomic classif)
			ratios <- data.frame(id = i)
			# Add total biovolume
			ratios$tot_bv <- sum(d[,"DW"])
			# And small and large cops biovolumes 
			ratios$small_bv <- sum(d[which(d$size <= 1),"DW"])
			ratios$large_bv <- sum(d[which(d$size > 1),"DW"])
			# Add number of samples
			eff <- length( unique(d$station_id) )
			ratios$effort <- eff
			# Correct biovolumes
			ratios$small_bv <- ratios$small_bv / ratios$effort
			ratios$large_bv <- ratios$large_bv / ratios$effort
			# Compute ratio
			ratios$ratio <- ratios$small_bv / ratios$large_bv
			ratios$season <- unique(d$season)
			# Return this
			return(ratios)
		} 
		
) # eo first lapply
t <- do.call(rbind,res)
t
gc() ; rm(res)

# Does ratio vary significantly across seasons? month ? 
fit <- aov(ratio ~ factor(season), data = t)
shapiro.test(fit$residuals)
# When based on EBioVol : p-value = 9.793e-11
# When based on concentration : p-value = 5.369e-08
# When based on DW : p-value = 1.23e-10

fligner.test(ratio ~ factor(season), data = t) # p-value = 0.4019 ; Variances not homogeneous
kruskal.test(ratio ~ factor(season), data = t) # 
positions <- c("spring","summer","fall")
plot <- ggplot(t, aes(x = factor(season), y = ratio)) + geom_boxplot(fill='#A4A4A4', color="black") + 
		scale_x_discrete(limits = positions) + theme_light() + xlab("Season") + ylab("Copepod size ratio\n(abundance)")


ggsave(plot = plot, filename = "boxplot_ratio_DW.pdf", dpi = 300, width = 4, height = 5.5)	
						
										
########################################################################### For defining biovolume classes !

# Define the classes of the EBv spectrum
summary(ddf3$EBioVol)
# Define classes
min <- min(ddf3$EBioVol)
max <- max(ddf3$EBioVol)
k <- 1.3
min2 <- exp( log(k) + log(min) )
# Define the classes
classes <- data.frame(i = c(1:200), min = NA, max = NA)
classes[1,c("min","max")] <- c(min,min2)
# Fill with for looping
for(i in 2:200) {	
		min <- classes[i-1,"max"]
		max <- exp( log(k) + log(min) )
		classes[i,"min"] <- min
		classes[i,"max"] <- max
} # eo for loop
### Get rid of the classes not covering your data
classes <- classes[which(classes$max <= max(ddf3$EBioVol)),]
### Compute classes' width
classes$width <- classes$max - classes$min
### and their midpoint (for plotting)
classes$mid <- classes$width/2
classes


### To plot the overall EBV spectrum
# Use the categories oin 'table' and compute their number of vignettes
res <- lapply(unique(classes$i), function(l) {	
				# Retrieve lower and upper boundaries
				low <- classes[which(classes$i == l), "min"]
				up  <- classes[which(classes$i == l), "max"]
				# Compute number of vignettes within this interval
				totbv <- sum(ddf3[which(ddf3$EBioVol < up & ddf3$EBioVol >= low),"EBioVol"])
				# normalize by class width
				wid <- classes[which(classes$i == l), "width"]
				norm_totbv <- totbv/wid
				# retrieve mid as well
				mid <- classes[which(classes$i == l), "mid"]
				# return
				return(data.frame(class = l, mid = mid, upper = up, lower = low, totbv = totbv, norm_totbv = norm_totbv))
		}	
) # eo lapply
t <- do.call(rbind, res)
rm(res)

quartz()
ggplot() + geom_point(aes(x = factor(round(log10(mid),2)), y = log10(norm_totbv)), data = t) + 
		geom_path(aes(x = factor(round(log10(mid),2)), y = log10(norm_totbv), group = 1), data = t) + theme_light() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("log(mm3)") + ylab("log(mm3/mm3.m3)")


### Plot for each sample id !!
for(i in unique(ddf3$sample_id) ) {
		
		message(paste("Doing ", i, sep = ""))
		# Filter vignettes according to i
		d <- ddf3[which(ddf3$sample_id == i),]
		
		res <- lapply(unique(classes$i), function(l) {	
						# Retrieve lower and upper boundaries
						low <- classes[which(classes$i == l), "min"]
						up  <- classes[which(classes$i == l), "max"]
						# Compute number of vignettes within this interval
						totbv <- sum(d[which(d$EBioVol < up & d$EBioVol >= low),"EBioVol"])
						# normalize by class width
						wid <- classes[which(classes$i == l), "width"]
						norm_totbv <- totbv/wid
						# retrieve mid as well
						mid <- classes[which(classes$i == l), "mid"]
						# return
						return(data.frame(class = l, mid = mid, upper = up, lower = low, totbv = totbv, norm_totbv = norm_totbv))
				}	
		) # eo lapply
		t <- do.call(rbind, res)
		rm(res)
		# quartz()
		#plot <- ggplot() + 
					#geom_col(aes(x = factor(class), y = n), fill = "#d73027", colour="black", data = t[which(t$class<=1),]) + 
					#geom_col(aes(x = factor(class), y = n), fill = "#4575b4", colour="black", data = t[which(t$class>1),]) + 
					#scale_y_continuous(limits = c(0,300)) + 
					#theme_light() + theme(axis.text.x= element_text(angle= 90, hjust= 1)) + xlab("Size (mm)") + ylab("Count")		
		
		#			
		plot <- ggplot() + geom_point(aes(x = factor(round(log10(mid),2)), y = log10(norm_totbv)), data = t) + 
					geom_path(aes(x = factor(round(log10(mid),2)), y = log10(norm_totbv), group = 1), data = t) + theme_light() + 
					theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("log(mm3)") + ylab("log(mm3/mm3.m3)")		
				
				
		# save
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
		ggsave(plot = plot, filename = paste(i,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/")
		
} # eo for loop	


# ggplot(ddf3, aes(x = factor(category), y = log10BioVol)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Category")


################################################# Compute and plot the relative contribution of each copepod category to the biovolume classes

classes$norm_totbv <- NA
for(c in classes$i) {
		message(paste("Doing ", c, sep = ""))
		
		# Find
		low <- classes[which(classes$i == c), "min"]
		up  <- classes[which(classes$i == c), "max"]
		# Compute total biovolume within this interval
		totbv <- sum(ddf3[which(ddf3$EBioVol < up & ddf3$EBioVol >= low),"EBioVol"])
		# normalize by class width
		wid <- classes[which(classes$i == c), "width"]
		norm_totbv <- totbv/wid
		
		# and provide
		classes[which(classes$i == c),"norm_totbv"] <- norm_totbv
	
} # eo for loop

# Get rid of last class (empty)
#classes <- classes[c(1:61),]

quartz()
ggplot() + geom_point(aes(x = log10(mid), y = log10(norm_totbv)), data = classes) + 
		#scale_x_discrete(limits = factor(log10(classes$mid))) + 
		geom_path(aes(x = log10(mid), y = log10(norm_totbv), group = 1), data = classes) + theme_light() + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("log(mm3)") + 
		ylab("log(mm3/mm3.m3)")

### Looks exactly the same as previously, cool ;)
### Now, compute contributions !

# c <- "Acartiidae"

contrib <- lapply(unique(ddf3$category), function(c) {
		
		message(paste("Doing ", c, sep = ""))
		contrib <- data.frame(class = classes$mid, width = classes$width, low = classes$min, up = classes$max, contrib = NA)
		
		# Fill this new data.frame with a for loop
		for(cl in contrib$class) {
					
				# Find
				low <- contrib[which(contrib$class == cl), "low"]
				up  <- contrib[which(contrib$class == cl), "up"]
				# Compute number of vignettes within this interval
				totbv <- sum(ddf3[which(ddf3$category == c & ddf3$EBioVol < up & ddf3$EBioVol >= low),"EBioVol"])
				# normalize by class width
				wid <- contrib[which(contrib$class == cl), "width"]
				norm_totbv <- totbv/wid
		
				# and provide
				contrib[which(contrib$class == cl),"contrib"] <- norm_totbv
	
		} # eo for loop
		
		# Return for cbinding later (so it is useless to return the classes as well...)
		return(contrib[,c("contrib")])
	
}) # eo lapply

contri <- do.call(cbind, contrib)
gc() ; rm(contrib)
dim(contri)
summary(contri)
### OK...

# Change colnames
colnames(contri) <- unique(ddf3$category)
classes <- cbind(classes, contri)
head(classes) ; dim(classes)
# Divide by classes$total_n
colnames(classes)
classes[,c(7:length(classes))] <- (classes[,c(7:length(classes))] / classes[,c(6)])*100
# classes <- na.omit(classes)

### Time for a ggplot, but first we need to change and melt the ddf a bit
  # You have too many classes of course...get rid of those for which max contrib is < 5%
maxima <- apply(classes[,c(7:length(classes))], 2, function(x) max(x, na.rm = TRUE))
maxima
names <- names(maxima[maxima > 10])
classes2 <- classes[,c("i","min","max","width","mid",names)]

mclass <- melt(classes2, id.vars = c("i","min","max","width","mid") )

########################################################################### #
### For when you have over 20 palettes !									#
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12				#
library("RColorBrewer")														#
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)		#
cols1 <- f("Paired")														#
cols2 <- f("Pastel1")													    #
########################################################################### #

quartz()
plot <- ggplot(data = na.omit(mclass)) + 
  			geom_bar(aes(x = factor(round(log10(mid), 2)), y = value, fill = factor(variable)), stat = "identity", position = "fill") + 
			scale_fill_manual(name = "Category", values = c(cols1, cols2) ) + theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
  			theme(axis.text.y= element_text(size = 9)) + xlab("log10(mm3)") + ylab("Relative contribution (%)") 

ggsave(plot = plot, filename = "Category_rel_contrib_mesozoo.pdf", dpi = 300, width = 10, height = 5)


################################################# Cluster the NBSS spectra ! (Komolgorov distance)

res <- lapply(unique(ddf3$station_id), function(i) {
				message(paste("Doing ", i, sep = ""))
				# Filter vignettes according to i
				d <- ddf3[which(ddf3$station_id == i),]
				# Use the categories oin 'table' and compute their number of vignettes
				res2 <- lapply(unique(classes$i), function(cl) {	
								# Find
								low <- classes[which(classes$i == cl), "min"]
								up  <- classes[which(classes$i == cl), "max"]
								# Compute number of vignettes within this interval
								totbv <- sum(d[which(d$EBioVol < up & d$EBioVol >= low),"EBioVol"])
								# normalize by class width
								wid <- classes[which(classes$i == cl), "width"]
								norm_totbv <- totbv/wid
								# retrieve mid as well
								mid <- classes[which(classes$i == cl), "mid"]
								# return
								return(data.frame(mid = mid, upper = up, lower = low, width = wid, totbv = totbv, norm_totbv = norm_totbv))
						}	
				) # eo lapply 2
				t <- do.call(rbind, res2) ; rm(res2)
				# Add some metadata
				t$station <- i
				t$year <- unique(d$year)
				t$month <- unique(d$month)
				# Return
				return(t)	
				rm(t, d)
			} # eo for loop
) # eo lapply

### rbind all stations
komo <- do.call(rbind, res)
head(komo)
unique(komo$station)
### Try to develop a function to obtain the pair-wise K distance
# http://r.789695.n4.nabble.com/Crosstable-like-analysis-ks-test-of-dataframe-td4644489.html
f <- function(x, y, ...,
         alternative = c("two.sided", "less", "greater"), exact = NULL){
     #w <- getOption("warn")
     #options(warn = -1)  # ignore warnings
     d <- ks.test(x, y, ..., alternative = alternative, exact = exact)$statistic
     #options(warn = w)
     d
}
n <- 1e1
# dat <- data.frame(X=rnorm(n), Y=runif(n), Z=rchisq(n, df=3))
# apply(dat, 2, function(x) apply(dat, 2, function(y) f(x, y)))
### Seems to work..
dat <- dcast(komo, mid ~ station, value.var = "norm_totbv")
# Looks fine
mat <- data.frame( apply(dat[,2:length(dat)], 2, function(x) apply(dat, 2, function(y) f(x, y))) )
dim(mat) # Get rid of the 1 row "class"
mat <- na.omit(as.dist(mat[2:length(dat),]))
### Perform clustering on this distance matx:
fit <- hclust(mat, method = "ward") 
# Add labels
fit$labels <- colnames(dat)[2:length(dat)]
plot(fit, hang= -1, labels = fit$labels)
groups <- cutree(fit, k = 5)
# Stations' names per cluster
names(groups[groups == 1])
names(groups[groups == 2])
names(groups[groups == 3])
names(groups[groups == 4])
names(groups[groups == 5])

komo$cluster <- NA
for(i in 1:nrow(komo)) {
		sta <- komo[i,"station"]
		komo[which(komo$station == sta),"cluster"] <- rep(groups[sta], times = nrow(komo[which(komo$station == sta),]) )	
} # 

komo$cluster <- factor(komo$cluster)

mean_spectra <- data.frame(class = unique(komo$mid), Group_1 = NA, sd1 = NA, Group_2 = NA, sd2 = NA, Group_3 = NA, sd3 = NA, 
				Group_4 = NA, sd4 = NA, Group_5 = NA, sd5 = NA)

for(c in unique(mean_spectra$class)) {
	
		# For cluster #1
		mean_1 <- mean(log10( komo[which(komo$mid == c & komo$cluster == "1"),"norm_totbv"] + 1 ))
		sd1 <- sd(log10(komo[which(komo$mid == c & komo$cluster == "1"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"Group_1"] <- mean_1
		mean_spectra[which(mean_spectra$class == c),"sd1"] <- sd1
		
		# For cluster #2
		mean_2 <- mean(log10(komo[which(komo$mid == c & komo$cluster == "2"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"Group_2"] <- mean_2
		sd2 <- sd(log10(komo[which(komo$mid == c & komo$cluster == "2"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"sd2"] <- sd2	
		
		# For cluster #3
		mean_3 <- mean(log10(komo[which(komo$mid == c & komo$cluster == "3"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"Group_3"] <- mean_3
		sd3 <- sd(log10(komo[which(komo$mid == c & komo$cluster == "3"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"sd3"] <- sd3	
		
		# For cluster #4
		mean_4 <- mean(log10(komo[which(komo$mid == c & komo$cluster == "4"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"Group_4"] <- mean_4
		sd4 <- sd(log10(komo[which(komo$mid == c & komo$cluster == "4"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"sd4"] <- sd4
		
		# For cluster #5
		mean_5 <- mean(log10(komo[which(komo$mid == c & komo$cluster == "5"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"Group_5"] <- mean_5
		sd5 <- sd(log10(komo[which(komo$mid == c & komo$cluster == "5"),"norm_totbv"] + 1))
		mean_spectra[which(mean_spectra$class == c),"sd5"] <- sd5
	
}


# miny <- min(mean_spectra$Group_1, mean_spectra$Group_2, mean_spectra$Group_3, mean_spectra$Group_4)
# maxy <- max(mean_spectra$Group_1, mean_spectra$Group_2, mean_spectra$Group_3, mean_spectra$Group_4)

### Plot them ! 
quartz()
ggplot() + geom_point(aes(x = log10(class), y = Group_1), data = mean_spectra) + 
		geom_ribbon(aes(x = log10(class), ymin = Group_1 - sd1, ymax = Group_1 + sd1), data = mean_spectra, fill = "#d53e4f", alpha = 0.4) + 
		scale_y_continuous(limits = c(-1,4)) + 
		geom_path(aes(x = log10(class), y = Group_1, group = 1), data = mean_spectra) + theme_bw() + 
		#theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("log(mm3)") + ylab("log(mm3/mm3.m3)") + ggtitle("Average NBSS for cluster 1")
		
# ### Plot them ! 
quartz()
ggplot() + geom_point(aes(x = round(log10(class),2), y = Group_2), data = mean_spectra) + 
		geom_ribbon(aes(x = round(log10(class),2), ymin = Group_2 - sd2, ymax = Group_2 + sd2), data = mean_spectra, fill = "#fc8d59", alpha = 0.4) + 
		scale_y_continuous(limits = c(-1,4)) + 
		geom_path(aes(x = round(log10(class),2), y = Group_2, group = 1), data = mean_spectra) + theme_bw() + 
		#theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("log(mm3)") + ylab("log(mm3/mm3.m3)") + ggtitle("Average NBSS for cluster 2")
		
		 
#### Plot them ! 
quartz()
ggplot() + geom_point(aes(x = log10(class), y = Group_3), data = mean_spectra) + 
		geom_ribbon(aes(x = log10(class), ymin = Group_3 - sd3, ymax = Group_3 + sd3), data = mean_spectra, fill = "#fee08b", alpha = 0.4) + 
		scale_y_continuous(limits = c(-1,4)) + 
		geom_path(aes(x = log10(class), y = Group_3, group = 1), data = mean_spectra) + theme_bw() + 
		#theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("log(mm3)") + ylab("log(mm3/mm3.m3)") + ggtitle("Average NBSS for cluster 3")		 


# ### Plot them ! 
quartz()
ggplot() + geom_point(aes(x = log10(class), y = Group_4), data = mean_spectra) + 
		geom_ribbon(aes(x = log10(class), ymin = Group_4 - sd4, ymax = Group_4 + sd4), data = mean_spectra, fill = "#e6f598", alpha = 0.4) + 
		scale_y_continuous(limits = c(-1,4)) + 
		geom_path(aes(x = log10(class), y = Group_4, group = 1), data = mean_spectra) + theme_bw() + 
		#theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("log(mm3)") + ylab("log(mm3/mm3.m3)") + ggtitle("Average NBSS for cluster 4")

#
quartz()
ggplot() + geom_point(aes(x = log10(class), y = Group_5), data = mean_spectra) + 
		geom_ribbon(aes(x = log10(class), ymin = Group_5 - sd5, ymax = Group_5 + sd5), data = mean_spectra, fill = "#99d594", alpha = 0.4) + 
		scale_y_continuous(limits = c(-1,4)) + 
		geom_path(aes(x = log10(class), y = Group_5, group = 1), data = mean_spectra) + theme_bw() + 
		#theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("log(mm3)") + ylab("log(mm3/mm3.m3)") + ggtitle("Average NBSS for cluster 5")
	
		
		
### Perform CA to confirm
komo$logmid <- log10(komo$mid)
dat <- dcast(komo[,c(1:9,11)], station+month+year ~ logmid, value.var = "norm_totbv")		
# OK.
colnames(dat)		

### May need to get rid of singletons: size classes that have only one occurrence: 
colSums(as.matrix(dat[,c(4:length(dat))]))
tokeep <- names(colSums(as.matrix(dat[,c(4:length(dat))])))[colSums(as.matrix(dat[,c(4:length(dat))])) > 0]
dat <- dat[,c("station","month","year",tokeep)]
### Modify colnames to have smehting that is readable
colnames(dat)[4:length(dat)] <- as.character( round( as.numeric( colnames(dat)[4:length(dat)] ), 2) )

# Perform CA
res.ca <- CA(dat[,c(4:40)])	
summary(res.ca)
### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2])

# For objects
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					month = dat$month, station = dat$station, year = dat$year )

head(AFCst)		
	
### Change colour according to month	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()			
			

### Change colour according to year	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()	



################################################# A) Overall copepod size spectrum

### Look at the distribution of log10BioVol for small and large unidentified copepods
ddf3$cond <- NA
ddf3[which(ddf3$size <= 1), "cond"] <- "small" 
ddf3[which(ddf3$size > 1), "cond"] <- "large" 

### Plot Copepoda size spectra
quartz()
ggplot(ddf3, aes(x= factor(round(size,1)), fill = cond)) + 
		geom_histogram(colour="black", binwidth = .5, alpha = .5, position = "identity", stat = "count") + 
		scale_fill_manual(name = "Size class", labels = c("Large","Small"), values = c("#2166ac","#d6604d")) + 
		theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Size (mm)") + ylab("Count")		

### For plotting all mesozooplankton size spectrum
quartz()
ggplot(ddf3, aes(x= factor(round(size,1)))) + 
		geom_histogram(colour="black", binwidth = .5, alpha = .5, position = "identity", stat = "count") + 
		theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Size (mm)") + ylab("Count")		

		
### Plot biovolume spectra of the 2 categories
quartz()
ggplot(ddf3[which(ddf3$cond %in% c("small","large")),], aes(x= factor(round(log10BioVol,1)), fill = cond)) + 
		geom_histogram(colour="black", binwidth = .5, alpha = .4, position = "identity", stat = "count") + 
		scale_fill_manual(name = "Size class", labels = c("Large","Small"), values = c("#2166ac","#d6604d")) + 
		theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("log(mm3)") + ylab("Count")		

### !!! Keep the same x and y scales across all plots: keep the same size classes from above
sizes <- as.numeric(unique(round(ddf3$size, 1)))
# Order them in increasing order
sizes <- sizes[order(sizes, decreasing = F)]
# Make a ddf out of it, by adding the upper and lower boundaries of each class
sizes <- data.frame(class = as.numeric(sizes))
sizes$lower <- sizes$class - 0.05
sizes$upper <- sizes$class + 0.05


################################################# B) Compute and plot for each sample

### a) Without fixing the size classes
for(i in unique(ddf3$sample_id) ) {
		message(paste("Doing ", i, sep = ""))
		# Filter vignettes according to i
		d <- ddf3[which(ddf3$sample_id == i),]
		# quartz()
		plot <- ggplot() + 
					geom_histogram(data = d, aes(x= factor(round(size,1)), colour="black", binwidth= .5, alpha= .5, position= "identity", stat= "count")) + 
					scale_x_discrete(levels = sizes, labels = factor(sizes)) + 
					#scale_fill_manual(name= "Size class", labels= c("Large","Small"), values= c("#2166ac","#d6604d")) + 
					theme_light() + theme(axis.text.x= element_text(angle= 90, hjust= 1)) + xlab("Size (mm)") + ylab("Count")		
		# save
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
		ggsave(plot = plot, filename = paste(i,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/")
} # eo for loop		


### b) With fixed classes
for(i in unique(ddf3$sample_id) ) {
		
		message(paste("Doing ", i, sep = ""))
		# Filter vignettes according to i
		d <- ddf3[which(ddf3$sample_id == i),]
		
		# Use the categories oin 'table' and compute their number of vignettes
		res <- lapply(sizes$class, function(s) {	
						# Retrieve lower and upper boundaries
						low <- sizes[which(sizes$class == s), "lower"]
						up  <- sizes[which(sizes$class == s), "upper"]
						# Compute number of vignettes within this interval
						n <- nrow(d[which(d$size < up & d$size >= low),])
						# mid <- classes[which(classes$class == s), "class"]
						# return
						return(data.frame(class = s, upper = up, lower = low, n = n))
				}		
		) # eo lapply
		t <- do.call(rbind, res) ; rm(res)
		# quartz()
		plot <- ggplot() + 
					geom_col(aes(x = factor(class), y = n), fill = "#d73027", colour="black", data = t[which(t$class<=1),]) + 
					geom_col(aes(x = factor(class), y = n), fill = "#4575b4", colour="black", data = t[which(t$class>1),]) + 
					#scale_y_continuous(limits = c(0,300)) + 
					theme_light() + theme(axis.text.x= element_text(angle= 90, hjust= 1)) + xlab("Size (mm)") + ylab("Count")		
		# save
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
		ggsave(plot = plot, filename = paste("size_spectra_copepoda_",i,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/")
		
} # eo for loop	


### c) Compute number of small, number of large and ratio small/large for each samples ! 
res <- lapply(unique(ddf3$sample_id), function(s) {	
				d <- ddf3[which(ddf3$sample_id == s),]
				message(paste("Doing ", s, sep = ""))
				# number of copepods
				n <- nrow(d)
				# total biovolume of copepods
				ebv <- sum(d[,c("EBioVol")])
				# size
				n_small <- nrow(d[which(d$size <= 1),])
				n_large <- nrow(d[which(d$size > 1),])
				size_ratio <- n_small/n_large
				# ebv
				ebv_small <- sum(d[which(d$size <= 1),"EBioVol"])
				ebv_large <- sum(d[which(d$size > 1),"EBioVol"])
				# Ratio
				ratio_ebv <- ebv_small/ebv_large
				return(data.frame(i = s, month = unique(d$month), year = unique(d$year), station = unique(d$station), 
								  n_total = n, ebv_total = ebv, n_small = n_small, n_large = n_large, ratio = ratio,
								  ebv_small = ebv_small, ebv_large = ebv_large, ratio_ebv = ratio_ebv))
		}		
) # eo lapply
t <- do.call(rbind, res) ; rm(res)
summary(t)
str(t)

# Melt to plot all of them with a facet_wrap() afterwards
t <- melt(t, id.vars = c("i","month","year","station"))
head(t)
quartz()
ggplot(t, aes(x = factor(month), y = value)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70") + 
		xlab("Month") + ylab("") + facet_wrap(~ variable, scales = "free_y")

quartz()
ggplot(t, aes(x = factor(year), y = value)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70") + 
		xlab("Year") + ylab("") + facet_wrap(~ variable, scales = "free_y")
	
	
				
### d) Looks very promising, cbind these variables with other cruise data and assess correlations with correlogram
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
env <- read.csv("PNMIR_cruise_data_abund_v_21_06_17.csv", h = T, sep = ";", dec = ".")
str(env)

### Use metadata to combine both, according to date and station
t$id <- paste(t$station, t$month, t$year, sep = "_")
env$id <- paste(env$station, env$month, env$year, sep = "_")
# Add empty columms for abundances
env[,colnames(t)[5:12]] <- NA
matching <- na.omit(t$id[match(env$id, t$id)]) ; length(matching)
for(i in matching) {
			stats <- t[which(t$id == i),c(5:12)]
			# cbind with env
			env[which(env$id == i),colnames(t)[5:12] ] <- stats
} # eo for loop
summary(env)
# Drop env$id
env <- env[,-c(77)]

### Assess correlations ! 
colnames(env)
cormat <- round( cor(na.omit(env[,c(9,12,15,17:21,23:28,76:84)]), method = "spearman"), 2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)	
# Re-order correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman\nCorrelation") +
			theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

### Okray, let's add some text
heatmap <- ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
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
heatmap		

### Some basic geom_points
# quartz()
# ggplot(env) + geom_point(aes(x= log(n_diato), y= log(ratio)), pch = 21, colour = "black", fill = "#2166ac") + theme_light()

### And PCA with supplementary variables of course of course 
res.pca <- PCA( X = na.omit(env[,c(9,12,15,17,19,20,23:24,76:84)]), scale.unit = T, quanti.sup = c(9:17) )
summary(res.pca)
### To examine the correlations between the variables and the suppl ones
res.pca$var$cor		
# PC1 = SST, SSS, number of dinoflagellates VERSUS nutrients and Chl-a/ Phaeo-a
# PC2 = phytoplankton cells, diatoms cells, nanoflagellate cells, Chl-a VERSUS phosphates, salinity

# Supplementary continuous variables
res.pca$quanti.sup$cor
# Are mainly linked to PC2: 
#                   Dim.1       Dim.2     
# other_Calanoida  0.04680533  0.34237997
# n_total         -0.05919184 -0.31989195
# ebv_total        0.07072875  0.26487035
# n_small         -0.08990239 -0.29509074
# n_large          0.04763914 -0.21069820
# ratio           -0.12891011  0.07944083
# ebv_small        0.01482845  0.39718393
# ebv_large        0.10977309  0.07350950
# ratio_ebv       -0.07756779  0.27489162
	
# Save
save(file = "PNMIR_cruise_data_abund_v_23_06_17.Rdata", env)	
	
		
		
		
################################################# C) Like you did with the Komolgorov distances, retrieve the size spectrum of each sample and perform a multivariate analysis (CA)

### Compute each spectrum in a lapply, then concatenate results in a dataframe that you will reshape to put size classes as variables.
res <- lapply(unique(ddf3$station_id), function(i) {
	
				message(paste("Doing ", i, sep = ""))
				# Filter vignettes according to i
				d <- ddf3[which(ddf3$station_id == i),]
				# Use the categories oin 'table' and compute their number of vignettes
				res2 <- lapply(unique(sizes$class), function(cl) {	
								# Retrieve lower and upper boundaries
								low <- sizes[which(sizes$class == cl), "lower"]
								up  <- sizes[which(sizes$class == cl), "upper"]
								# Compute number of vignettes within this interval
								n <- nrow(d[which(d$size < up & d$size >= low),])
								# return
								return(data.frame(class = cl, upper = up, lower = low, n = n))
						}	
				) # eo lapply 2
				# rbind
				t <- do.call(rbind, res2) ; rm(res2)
				# Add some metadata
				t$station <- i
				t$year <- unique(d$year)
				t$month <- unique(d$month)
				# Return
				return(t)	
				rm(t, d)
		
			} # eo for loop
) # eo lapply

### rbind all stations
ddf4 <- do.call(rbind, res)
head(ddf4)
unique(ddf4$station)

### reshape to obtain size classes as variables
dat <- dcast(ddf4, station+month+year ~ class, value.var = "n")		
# OK.
colnames(dat)		

### May need to get rid of singletons: size classes that have only one occurrence: 
colSums(as.matrix(dat[,c(4:length(dat))]))

res.ca <- CA(dat[,c(4:40)])		
summary(res.ca)	
### Check size classes' contributions to CA axes
round(res.ca$col$contrib[,1],3)	
round(res.ca$col$contrib[,2],3)	
### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2])

# For objects
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					month = dat$month, station = dat$station, year = dat$year )

head(AFCst)		
	
### Change colour according to month	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()			
			

### Change colour according to year	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()	
	
	
### 2015 seems to present higher contributions of the large size classes (ranging between 1.5mm and 3.5mm)
# Perform clustering on the first 3 CA axes 	
mat <- dist(AFCst[,c(1:2)], "euclidean")
fit <- hclust(mat, method = "ward") 
# Add labels
fit$labels <- AFCst$station
quartz()
plot(fit, hang= -1, labels = fit$labels)
groups <- cutree(fit, k = 6)
names(groups[groups == 1])
names(groups[groups == 2])
names(groups[groups == 3])
names(groups[groups == 4])
names(groups[groups == 5])
names(groups[groups == 6])

##### Not great...


################################################# D) Perform statistical tests on ratios, number of small copepoda etc. to assess differences across month, year etc.
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
ddf <- get(load("PNMIR_cruise_data_abund_v_23_06_17.Rdata"))
colnames(ddf)	

### Switch the month of the samples taken in April to May...
ddf[which(ddf$month == 4),"month"] <- 5
ddf[which(ddf$month == 6),"month"] <- 5
	
### Test the veraition of number of copepods throughout the years	
fit <- aov(n_total ~ factor(year), data = ddf)
shapiro.test(fit$residuals) # p-value = 0.02142 ; normally distributed !
summary(fit)
quartz()
ggplot(ddf[-which(ddf$year == 2010),], aes(x = factor(month), y = ratio_ebv)) + geom_boxplot(fill='#A4A4A4', color="black") + theme_light()  
quartz()
ggplot(ddf[-which(ddf$year == 2010),], aes(x = factor(month), y = ebv_small)) + geom_boxplot(fill='#A4A4A4', color="black") + theme_light()  
quartz()
ggplot(ddf[-which(ddf$year == 2010),], aes(x = factor(month), y = ebv_large)) + geom_boxplot(fill='#A4A4A4', color="black") + theme_light()  


# Across month ? 
fit <- aov(n_large ~ factor(month), data = ddf)
shapiro.test(fit$residuals) # Not normal
fligner.test(n_large ~ factor(month), data = ddf) # p-value = 0.4019 ; Variances not homogeneous
kruskal.test(n_large ~ factor(month), data = ddf)
ggplot(ddf[-which(ddf$year == 2010),], aes(x = factor(month), y = n_large)) + geom_boxplot(fill='#A4A4A4', color="black") + theme_light()  

summary(lm( log(Calanoida+Cladocera+Oikopleuridae+Hydrozoa+Cyclopoida+nauplii+Cirripedia+Chaetognatha+Euphausiacea+Annelida+Amphipoda+egg+Poecilostomatoida) ~ max_depth, data = ddf))				
	
ggplot(ddf) + geom_point(aes(x=max_depth,y=ratio), data = ddf) + theme_classic()	
	
	
	
			
# --------------------------------------------------------------------------------------------------------------------------

##### 3°) Compute and plot the normalized biovolume size spectrum of the small and large Copepoda - overall and then for each sample

unique(cut(x = ddf3[ddf3$category == "Calanoida","EBioVol"], quantile(ddf3[ddf3$category == "Calanoida","EBioVol"], probs = seq(0,1,0.05)) ))

################################################# A) Overall biovolume spectrum for large and smaller Copepoda
summary(ddf3$EBioVol)
# Define classes
min <- min(ddf3$EBioVol)
max <- max(ddf3$EBioVol)
k <- 1.3
min2 <- exp( log(k) + log(min) )
# Define the classes
classes <- data.frame(i = c(1:200), min = NA, max = NA)
classes[1,c("min","max")] <- c(min,min2)
# Fill with for looping
for(i in 2:200) {	
		min <- classes[i-1,"max"]
		max <- exp( log(k) + log(min) )
		classes[i,"min"] <- min
		classes[i,"max"] <- max
} # eo for loop
### Get rid of the classes not covering your data
classes <- classes[which(classes$max <= max(ddf3$EBioVol)),]
### Compute classes' width
classes$width <- classes$max - classes$min
### and their midpoint (for plotting)
classes$mid <- classes$width/2
classes



################################################# B) Compute and plot for each sample

### a) Without fixing the EBV classes (just round)
for(i in unique(ddf3$sample_id) ) {
		message(paste("Doing ", i, sep = ""))
		# Filter vignettes according to i
		d <- ddf3[which(ddf3$sample_id == i),]
		# quartz()
		plot <- ggplot(d[which(d$cond %in% c("small","large")),], aes(x= factor(round(log10BioVol,1)), fill = cond)) + 
				geom_histogram(colour="black", binwidth = .5, alpha = .4, position = "identity", stat = "count") + 
				scale_fill_manual(name = "Size class", labels = c("Large","Small"), values = c("#2166ac","#d6604d")) + 
				scale_y_continuous(limits = c(0,250)) + 
				theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("log(mm3)") + ylab("Count")		
		# save
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
		ggsave(plot = plot, filename = paste("EBV_spectra_copepoda_",i,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/")
} # eo for loop		


### b) Using the EBV classes from 'classes' object
for(i in unique(ddf3$sample_id) ) {
		message(paste("Doing ", i, sep = ""))
		# Filter vignettes according to i
		d <- ddf3[which(ddf3$sample_id == i),]
		# Use the categories oin 'table' and compute their number of vignettes
		res <- lapply(classes$mid, function(c) {	
						# Retrieve lower and upper boundaries
						low <- classes[which(classes$mid == c), "min"]
						up  <- classes[which(classes$mid == c), "max"]
						# Compute number of vignettes within this interval
						n_small <- nrow(d[which(d$EBioVol < up & d$EBioVol >= low & d$size <= 1),])
						n_large <- nrow(d[which(d$EBioVol < up & d$EBioVol >= low & d$size > 1),])
						# normalize by class width
						wid <- classes[which(classes$mid == c), "width"]
						normn_small <- n_small/wid
						normn_large <- n_large/wid
						# return
						return(data.frame(class = c, upper = up, lower = low, n_small = n_small, n_large = n_large, 
										normn_small = normn_small, normn_large = normn_large))
				}		
		) # eo lapply
		t <- do.call(rbind, res) ; rm(res)
		# quartz()
		plot <- ggplot() + 
					geom_col(aes(x= factor(round(log10(class),2)), y= n_small), alpha=.5, fill = "#d73027", colour="black", data = t) + 
					geom_col(aes(x= factor(round(log10(class),2)), y= n_large), alpha=.5, fill = "#4575b4", colour="black", data = t) + 
					xlab("log(mm3)") + ylab("Normalized count") + 
					theme_light() + theme(axis.text.x= element_text(angle= 90, hjust= 1))
		# save
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
		ggsave(plot = plot, filename = paste("EBV_spectra_copepoda_",i,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/")
} # eo for loop	


### c) Perform CA on biovolume classes for the copepoda
### 28/06/17: Compute each spectrum in a lapply, then concatenate results in a dataframe that you will reshape to put size classes as variables.
res <- lapply(unique(ddf3$station_id), function(i) {
	
				message(paste("Doing ", i, sep = ""))
				# Filter vignettes according to i
				d <- ddf3[which(ddf3$station_id == i),]
				
				# Use the categories oin 'table' and compute their number of vignettes
				res2 <- lapply(classes$i, function(cl) {	
								# Retrieve lower and upper boundaries
								low <- classes[which(classes$i == cl), "min"]
								up  <- classes[which(classes$i == cl), "max"]
								# Compute number of vignettes within this interval
								bvtot <- sum(d[which(d$EBioVol < up & d$EBioVol >= low),"EBioVol"])
								# normalize by class width 
								wid <- classes[which(classes$i == cl), "width"]
								norm_bvtot <- bvtot/wid
								# Find midpoint of the class
								mid <- classes[which(classes$i == cl), "mid"]
								# return
								return( data.frame(class = cl, upper = up, lower = low, width = wid, mid = mid, bvtot = bvtot, norm_bvtot = norm_bvtot) )
						} # eo FUN	
				) # eo lapply
				t <- do.call(rbind, res2) ; rm(res2)
				# Add some metadata
				t$station <- i
				t$year <- unique(d$year)
				t$month <- unique(d$month)
				# Return
				return(t)	
				rm(t, d)
		
			} # eo for loop
			
) # eo lapply

### rbind all stations
ddf4 <- do.call(rbind, res)
head(ddf4)
summary(ddf4)
unique(ddf4$station)

ddf4$lognorm_bvtot <- log10(ddf4$norm_bvtot+1)
ddf4$logmid <- log10(ddf4$mid)

### reshape to obtain size classes as variables
dat <- dcast(ddf4[,c("station","month","year","logmid","lognorm_bvtot")], station + month + year ~ logmid, value.var = "lognorm_bvtot")		
# OK.
colnames(dat)		

### May need to get rid of singletons: size classes that have only one occurrence: 
names(colSums(as.matrix(dat[,c(4:length(dat))])))[colSums(as.matrix(dat[,c(4:length(dat))])) == 0]
tokeep <- names(colSums(as.matrix(dat[,c(4:length(dat))])))[colSums(as.matrix(dat[,c(4:length(dat))])) > 0]
dim(dat)
dat <- dat[,c('station','month','year',tokeep)]
dim(dat)
colSums(as.matrix(dat[,c(4:length(dat))]))

### For plotting reasons, change the columns'names to numbers of increasing value
colnames(dat)[4:length(dat)] <- as.character( round( as.numeric( colnames(dat)[4:length(dat)] ), 2) )
res.ca <- CA(dat[,c(4:length(dat))])		
summary(res.ca)	
### Check size classes' contributions to CA axes
round(res.ca$col$contrib[,1],3)	
round(res.ca$col$contrib[,2],3)	
### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2])

# For objects
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					month = dat$month, station = dat$station, year = dat$year )

head(AFCst)		
	
### Change colour according to month	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()			
			

### Change colour according to year	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()	









# --------------------------------------------------------------------------------------------------------------------------

##### 4°) Do the above but by considering the area to compute the biovolume, instead of the major and minor axes...
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/Data_08_06_17/")
proj_dir <- getwd()
projects <- dir()
projects
p <- "export_434_20170608_1544"

# Prior to lapply: use one project to save the colnames of interest: object_id, object_annotation + morphometric traits
setwd(paste(proj_dir,"/",p,"/", sep = ""))
all_samples <- dir()
# Read and rbind the samples from each project
samples <- lapply(all_samples, function(s) {	
			# Load and return
			ss <- read.table(file = s, sep = '\t', header = T)
			return(ss)	
}) # eo 2nd lapply 
all_samples <- do.call(rbind, samples)
# dim(all_samples)
# colnames(all_samples)
### Choose the columns of interest. But let's keep all for now.
cols <- colnames(all_samples)[c(1:119,121:159)] # because of "process_particle_pixel_size__m" in "export_434_20170608_1544"
rm(all_samples, samples)

### Now, load all data for real
projects
# p <- "export_442_20170608_1539"
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
				return(all_samples[,cols])
}) # eo lapply
ddf <- do.call(rbind, projs)
# str(ddf)
# colnames(ddf)
### Select validated vignettes only.
ddf2 <- ddf[which(ddf$object_annotation_status == "validated"),] # Validated vignettes only
### Get rid of horizontal hauls...
horizontal <- ddf2$acq_id[grep("_h", ddf2$acq_id)]
ddf2 <- ddf2[!(ddf2$acq_id %in% horizontal),]
### Get rid of non-living particles and "temporary"
nl <- ddf2$object_annotation_hierarchy[grep("not-living", ddf2$object_annotation_hierarchy)]
ddf2 <- ddf2[!(ddf2$object_annotation_hierarchy %in% nl),]
temp <- ddf2$object_annotation_hierarchy[grep("temporary", ddf2$object_annotation_hierarchy)]
ddf2 <- ddf2[!(ddf2$object_annotation_hierarchy %in% temp),]
molene <- ddf2$object_id[grep("molene", ddf2$object_id)]
Molene <- ddf2$object_id[grep("Molene", ddf2$object_id)]
douar <- ddf2$object_id[grep("douarnenez", ddf2$object_id)]
ddf2 <- ddf2[!(ddf2$object_id %in% molene),]
ddf2 <- ddf2[!(ddf2$object_id %in% douar),]
ddf2 <- ddf2[!(ddf2$object_id %in% Molene),]
rm(douar, Molene, molene)
# Replace caps
ddf2$object_id <- str_replace_all(as.character(ddf2$object_id), "D", "d")
ddf2$object_id <- str_replace_all(as.character(ddf2$object_id), "B", "b")
### And remove 2010 samples because their id do not follow the same formating
s2010 <- ddf2$object_id[grep("wp2_2010", ddf2$object_id)]
ddf2 <- ddf2[!(ddf2$object_id %in% s2010),]
### Add ids
ID <- data.frame( do.call(rbind, strsplit(x = as.character(ddf2$object_id), split = "_")) )
colnames(ID)[1:4] <- c("net", "station", "date", "type")
ddf2$station <- ID$station
ID$date <- lubridate::ymd(ID$date)
ddf2$date <- ID$date
ddf2$month <- lubridate::month(ddf2$date)
ddf2$year <- lubridate::year(ddf2$date)
ddf2$station_id <- paste(ddf2$station,ddf2$month,ddf2$year, sep = "_")
gc()

# Extract the classification thanks to ddf2$object_annotation_hierarchy and select what is relevant.
split_status <- strsplit(x = as.character(ddf2$object_annotation_hierarchy), split = ">") 
length(split_status)
# split_status[[276312]][length(split_status[[276312]])] ### here's how you can do it
# Retrieve status and a few parent categories, provide back
categories <- lapply(c(1:length(split_status)), function(i) {
					message(paste(i))
					category <- split_status[[i]] [length(split_status[[i]])] 
					# Parent categories
					parent1 <- split_status[[i]] [length(split_status[[i]]) - 1] 
					parent2 <- split_status[[i]] [length(split_status[[i]]) - 2] 				
					# Return
					return(data.frame(i = i, category = category, parent1 = parent1, parent2 = parent2 ))
}) # eo lapply 
# rbind as usual
cat <- do.call(rbind, categories)	 
dim(cat)
# cbind()
ddf3 <- cbind(ddf2, cat[,c(2:4)])
rm(cat, categories, temp, nl, horizontal, ddf, split_status)


### Select only the Copepoda categories ! 
unique(ddf3$category)
nrow(ddf3) # 106598
ddf3 <- ddf3[which(ddf3$category %in% c("Acartiidae","Calanidae","Calanoida","Centropagidae","Euchaetidae","Cyclopoida",
									"Copepoda","Oithonidae","Oncaeidae","Temoridae","Corycaeidae","Harpacticoida",
									"Monstrilloida","Candaciidae","Pontellidae","Poecilostomatoida","Sapphirinidae","Oncaea")),] 

nrow(ddf3) # 62369

ddf3$concentration <- 1 * ddf3$acq_sub_part / ddf3$sample_tot_vol
ddf3$BioVol <- ddf3$concentration * (4/3 * pi *(sqrt((ddf3$object_area*(0.0106)^2)/pi))^3)
# summary(log10(ddf3$EBioVol))
ddf3$log10BioVol <- log10(ddf3$BioVol)
ddf3$size <- ddf3$object_major*0.0106
summary(ddf3$BioVol)
summary(ddf3$log10BioVol)

# Define classes
min <- min(ddf3$BioVol)
max <- max(ddf3$BioVol)
k <- 1.2
min2 <- exp( log(k) + log(min) )
# Define the classes
classes <- data.frame(i = c(1:200), min = NA, max = NA)
classes[1,c("min","max")] <- c(min,min2)
# Fill with for looping
for(i in 2:200) {	
		min <- classes[i-1,"max"]
		max <- exp( log(k) + log(min) )
		classes[i,"min"] <- min
		classes[i,"max"] <- max
} # eo for loop
### Get rid of the classes not covering your data
classes <- classes[which(classes$max <= max(ddf3$BioVol)),]
### Compute classes' width
classes$width <- classes$max - classes$min
### and their midpoint (for plotting)
classes$mid <- classes$width/2
classes

# For tetsing
i <- unique(ddf3$sample_id)[1]
### Compute and plot for each sample
for(i in unique(ddf3$sample_id) ) {
		message(paste("Doing ", i, sep = ""))
		# Filter vignettes according to i
		d <- ddf3[which(ddf3$sample_id == i),]
		
		# Use the categories oin 'table' and compute their number of vignettes
		res <- lapply(classes$mid, function(c) {	
						# Retrieve lower and upper boundaries
						low <- classes[which(classes$mid == c), "min"]
						up  <- classes[which(classes$mid == c), "max"]
						# Compute number of vignettes within this interval
						#bv_small <- sum(d[which(d$BioVol < up & d$BioVol >= low & d$size <= 1),"BioVol"])
						#bv_large <- sum(d[which(d$BioVol < up & d$BioVol >= low & d$size > 1),"BioVol"])
						bv_total <- sum(d[which(d$BioVol < up & d$BioVol >= low),"BioVol"])
						# normalize by class width
						wid <- classes[which(classes$mid == c), "width"]
						#normbv_small <- bv_small/wid
						#normbv_large <- bv_large/wid
						normbv_total <- bv_total/wid
						# return
						return(data.frame(class = c, upper = up, lower = low, bv_total = bv_total, normbv_total = normbv_total))
				}		
		) # eo lapply
		t <- do.call(rbind, res) ; rm(res)
		# quartz()
		#plot <- ggplot() + 
					#geom_col(aes(x= factor(round(log10(class),2)), y= log10(normn_small+1)), alpha=.5, fill = "#d73027", colour="black", data = t) + 
					#geom_col(aes(x= factor(round(log10(class),2)), y= log10(normn_large+1)), alpha=.5, fill = "#4575b4", colour="black", data = t) + 
					#xlab("log(mm3)") + ylab("log(mm3.m^-3.mm^-3)") + 
					#theme_light() + theme(axis.text.x= element_text(angle= 90, hjust= 1))
					
		# quartz()
		plot <- ggplot() + geom_point(aes(x= factor(round(log10(class),2)), y = normbv_total), fill = "black", colour="black", data = t) + 
						geom_path(aes(x= factor(round(log10(class),2)), y = normbv_total, group = 1), colour="black", data = t) +
						xlab("log(mm3)") + ylab("mm3.m^-3.mm^-3") + theme_light() + theme(axis.text.x= element_text(angle= 90, hjust= 1))			
							
		# save
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/graphes/")
		ggsave(plot = plot, filename = paste(i,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
		setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/")
} # eo for loop	



### c) Perform CA on biovolume classes for the copepoda
res <- lapply(unique(ddf3$station_id), function(i) {
	
				message(paste("Doing ", i, sep = ""))
				# Filter vignettes according to i
				d <- ddf3[which(ddf3$station_id == i),]
				# Use the categories oin 'table' and compute their number of vignettes
				res2 <- lapply(classes$i, function(cl) {	
								# Retrieve lower and upper boundaries
								low <- classes[which(classes$i == cl), "min"]
								up  <- classes[which(classes$i == cl), "max"]
								# Compute number of vignettes within this interval
								n <- nrow(d[which(d$BioVol < up & d$BioVol >= low),])

								# normalize by class width 
								wid <- classes[which(classes$i == cl), "width"]
								normn <- n/wid
								mid <- classes[which(classes$i == cl), "mid"]
								# return
								return(data.frame(class = cl, upper = up, lower = low, width = wid, mid = mid, n = n, normn = normn))
						}	
				)
						
				t <- do.call(rbind, res2) ; rm(res2)
				# Add some metadata
				t$station <- i
				t$year <- unique(d$year)
				t$month <- unique(d$month)
				# Return
				return(t)	
				rm(t, d)
		
			} # eo for loop
			
) # eo lapply

### rbind all stations
ddf4 <- do.call(rbind, res)
head(ddf4)
summary(ddf4)
unique(ddf4$station)

ddf4$lognormn <- log10(ddf4$normn+1)
ddf4$logmid <- log10(ddf4$mid+1)

### reshape to obtain size classes as variables
dat <- dcast(ddf4[,c("station","month","year","logmid","lognormn")], station+month+year ~ logmid, value.var = "lognormn")		
# OK.
colnames(dat)		

### May need to get rid of singletons: size classes that have only one occurrence: 
colSums(as.matrix(dat[,c(4:length(dat))]))
# --> the third is empty

### For plotting reasons, change the columns'names to numbers of increasing value
colnames(dat)[4:length(dat)] <- factor(c(4:length(dat)))

res.ca <- CA(dat[,c(4,5,7:49)])		
summary(res.ca)	
### Check size classes' contributions to CA axes
round(res.ca$col$contrib[,1],3)	
round(res.ca$col$contrib[,2],3)	
### For ggplotting: need to prepare a few dataframes.
AFCsp <- data.frame(Ax1 = res.ca$col$coord[,1],
                    Ax2 = res.ca$col$coord[,2],
					Ax3 = res.ca$col$coord[,3],
                    Cos2_1 = res.ca$col$cos2[,1] + res.ca$col$cos2[,2])

# For objects
AFCst <- data.frame(Ax1 = res.ca$row$coord[,1],
                    Ax2 = res.ca$row$coord[,2],
					Ax3 = res.ca$row$coord[,3],
                    Cos2_1 = res.ca$row$cos2[,1] + res.ca$row$cos2[,2], 
					month = dat$month, station = dat$station, year = dat$year )

head(AFCst)		
	
### Change colour according to month	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(month) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Month", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()			
			

### Change colour according to year	
quartz()
ggplot() +
  	geom_point(aes(x = Ax1, y = Ax2), data = AFCsp, alpha = 0.7, pch = 23, colour = "black", size = 3)+
  	geom_point(aes(x = Ax1, y = Ax2, size = Cos2_1, fill = factor(year) ), data = AFCst, alpha = 0.7, colour = "black", pch = 21) +
	#scale_shape_manual("Year", values = c(21,22,23,24)) + 
	scale_fill_manual("Year", values = c("#e6f598","#abdda4","#fdae61","#d7191c","#2b83ba")) + 
	scale_x_continuous(paste("CA 1 (",round(res.ca$eig$per[1], 2),"%)", sep = "")) + 
	scale_y_continuous(paste("CA 2 (",round(res.ca$eig$per[2],2),"%)", sep = "")) + 
	geom_text(aes(x = Ax1, y = Ax2+0.05, label = rownames(AFCsp)), data = AFCsp, size = 3) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
  	#geom_text(aes(x = Ax1, y = Ax2+0.01, label = rownames(AFCst)), data = AFCst, size = 3) +
	theme_light()	






