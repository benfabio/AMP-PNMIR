
##### 25/07/2017: R Script to re-plot the copepods' groups relative contributions to the NBSS spectrum.

### Latest update: 03/08/2017

library("dplyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("vegan")
library("FactoMineR")
library("lubridate")
library("viridis")

# --------------------------------------------------------------------------------------------------------------------------

##### 1°) Load data

setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")d
ddf3 <- get(load("ddf3_06_07_17.Rdata"))
str(ddf3)
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

##### 06/07/17: After Van der Lingen et al. (2006) : look at the seasonal/monthly/annual variations of the % of :
#	- Cyclopoids
#	- small unidentified Calanoida (<1mm)
#	- larger unidentified Calanoida and Calanidae
ddf3$season <- NA
ddf3[which(ddf3$month %in% c(4,5,6)),"season"] <- "spring"
ddf3[which(ddf3$month == 7),"season"] <- "summer"
ddf3[which(ddf3$month == 10),"season"] <- "fall"
ddf3$id <- paste(ddf3$season, ddf3$year, sep = "_" )

### Select only the Copepoda categories ! 
# unique(ddf3$category)
ddf3 <- ddf3[which(ddf3$category %in% c("Acartiidae","Calanidae","Calanoida","Centropagidae","Euchaetidae",
										"Oithonidae","Oncaeidae","Temoridae","Corycaeidae","Harpacticoida",
										"Candaciidae","Pontellidae","Sapphirinidae","Oncaea")),] 
										
	
##### 2°) Define the NBSS from the ellipsoïdal spectra

# Define the classes of the EBv spectrum
summary(ddf3$EBioVol)
# Define classes
min <- min(ddf3$EBioVol)
max <- max(ddf3$EBioVol)
k <- 1.35
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

		
### Identify the EBV class that corresponds best to the 1mm threshold 
summary( round(log10(ddf3[which(ddf3$size >= 0.95 & ddf3$size <= 1.05),"EBioVol"]),2) )
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.7800 -0.9300 -0.6200 -0.6426 -0.3700  0.9100
mean(round(log10(ddf3[which(ddf3$size >= 0.95 & ddf3$size <= 1.05),"EBioVol"]),2)) # -0.642595
sd(round(log10(ddf3[which(ddf3$size >= 0.95 & ddf3$size <= 1.05),"EBioVol"]),2))  # 0.4587372



##### 3°) Compute and plot the relative contribution of each copepod category to the biovolume classes
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
names <- names(maxima[maxima > 5])
classes2 <- classes[,c("i","min","max","width","mid",names)]
# Melt
mclass <- melt(classes2, id.vars = c("i","min","max","width","mid") )

########################################################################### #
### For when you have over 20 palettes !									
# http://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12				
library("RColorBrewer")														
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)		
cols1 <- f("Paired")														
cols2 <- f("Pastel1")	
cols1 <- cols1[c(1:10)]																														    								
########################################################################### #

### Before plotting: replace "Calanoida" by "Unidentified Calanoida"
levels(mclass$variable)[4] <- 'Unidentified Calanoida'
mclass$variable[mclass$variable == "Calanoida"] <- "Unidentified Calanoida"

# quartz()
plot <- ggplot(data = na.omit(mclass)) + 
  			geom_bar(aes(x = factor(round(log10(mid), 2)), y = value, fill = factor(variable)), stat = "identity", position = "fill") + 
			scale_fill_manual(name = "", values = spectral12 ) + 
			theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
  			theme(axis.text.y= element_text(size = 9)) + xlab("log10(mm3)") + ylab("Relative contribution (%)") 

ggsave(plot = plot, filename = "Category_rel_contrib_copepoda.pdf", dpi = 300, width = 10, height = 5)


### For trying different and cool palettes
spectral12 <- sample(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#74add1","#e6f598","#abdda4","#66c2a5","#4575b4","#5e4fa2"), 11)
quartz()
ggplot(data = na.omit(mclass)) + 
  			geom_bar(aes(x = factor(round(log10(mid), 2)), y = value, fill = factor(variable)), stat = "identity", position = "fill") + 
			scale_fill_manual(name = "", values = spectral12 ) + 
			theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
  			theme(axis.text.y= element_text(size = 9)) + xlab("log10(mm3)") + ylab("Relative contribution (%)") 
			
			
			
##### 3°) Compute and plot the relative contribution of each copepod category to the size spectrum : same as above but sticking to size measures (mm) in 0.1 bins

### First, simply plot distribution of size for non identified Calanoida
spectral <- sample(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#74add1","#e6f598","#abdda4","#66c2a5","#4575b4","#5e4fa2","#313695","#80cdc1"), 13)
quartz() 
ggplot(ddf3, aes(x = factor(round(size,1)), fill = factor(category))) +
	 	geom_histogram(binwidth=.5, colour = "black", stat = "count") + 
		scale_fill_manual(name = "", values = spectral) + 
		theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Size (mm)") + ylab("Count")


### May need to define size classes logarithmically, like for 
min <- min(ddf3$size)
max <- max(ddf3$size)
k <- 1.1
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
classes <- classes[which(classes$max <= max(ddf3$size)),]
### Compute classes' width
classes$width <- classes$max - classes$min
### and their midpoint (for plotting)
classes$mid <- (classes$max + classes$min)/2
classes

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
				tot <- sum(ddf3[which(ddf3$category == c & ddf3$size < up & ddf3$size >= low),"size"])
				# normalize by class width
				wid <- contrib[which(contrib$class == cl), "width"]
				norm_tot <- tot/wid
				# and provide
				contrib[which(contrib$class == cl),"contrib"] <- norm_tot
		} # eo for loop
		
		# Return for cbinding later (so it is useless to return the classes as well...)
		return(contrib[,c("contrib")])
	
}) # eo lapply

contri <- do.call(cbind, contrib)
gc() ; rm(contrib)
dim(contri)
contri
### OK !

# Change colnames
colnames(contri) <- unique(ddf3$category)
classes <- cbind(classes, contri)
head(classes) ; dim(classes)
# Divide by classes$total_n
library("matrixStats")
classes$total <- rowSums(as.matrix(classes[,c(7:length(classes))]))
colnames(classes)
classes[,c(6:18)] <- (classes[,c(6:18)] / classes[,c(19)])*100
# classes <- na.omit(classes)

### Time for a ggplot, but first we need to change and melt the ddf a bit
  # You have too many classes of course...get rid of those for which max contrib is < 5%
maxima <- apply(classes[,c(6:18)], 2, function(x) max(x, na.rm = TRUE))
maxima
names <- names(maxima[maxima > 1])
classes2 <- classes[,c("i","min","max","width","mid",names)]
# Melt
mclass <- melt(classes2, id.vars = c("i","min","max","width","mid") )

### Before plotting: replace "Calanoida" by "Unidentified Calanoida"
levels(mclass$variable)[4] <- 'Unidentified Calanoida'
#mclass$variable[mclass$variable == "Calanoida"] <- "Unidentified Calanoida"

quartz()
#spectral <- sample(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#74add1","#e6f598","#abdda4","#66c2a5","#4575b4","#5e4fa2","#313695","#80cdc1"), 13)
plot <- ggplot(data = na.omit(mclass)) + 
  		geom_bar(aes(x = factor(round(mid,2)), y = value, fill = factor(variable)), stat = "identity", position = "fill") + 
		scale_fill_manual(name = "", values = spectral ) + 
		theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
  		theme(axis.text.y= element_text(size = 9)) + xlab("Size class (mm)") + ylab("Relative contribution (%)") 

ggsave(plot = plot, filename = "Category_rel_contrib_copepoda.pdf", dpi = 300, width = 10, height = 5)



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------


##### 28/07/2017: Compute and plot the relative & absolute contributions of copepod groups to size classes for each month and year (3x3)

setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
ddf3 <- get(load("ddf3_06_07_17.Rdata"))
str(ddf3)
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
# Add season
ddf3$season <- NA
ddf3[which(ddf3$month %in% c(4,5,6)),"season"] <- "spring"
ddf3[which(ddf3$month == 7),"season"] <- "summer"
ddf3[which(ddf3$month == 10),"season"] <- "fall"
ddf3$id <- paste(ddf3$season, ddf3$year, sep = "_" )

### Select only the Copepoda categories ! 
ddf3 <- ddf3[which(ddf3$category %in% c("Acartiidae","Calanidae","Calanoida","Centropagidae","Euchaetidae",
										"Oithonidae","Oncaeidae","Temoridae","Corycaeidae","Harpacticoida",
										"Candaciidae","Pontellidae","Sapphirinidae","Oncaea")),] 

### Add dry weights:
ddf3$DW <- 45.25*(ddf3$object_area*(0.0106)^2)^1.59
ddf3$DW <- (ddf3$DW)*ddf3$concentration									
										
### Define size classes logarithmically
min <- min(ddf3$size)
max <- max(ddf3$size)
k <- 1.1
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
classes <- classes[which(classes$max <= max(ddf3$size)),]
### Compute classes' width
classes$width <- classes$max - classes$min
### and their midpoint (for plotting)
classes$mid <- (classes$max + classes$min)/2
classes


##### !!! For plotting: as a copepod category may be absent from one sample, you need to associate a precise color code to each category so they don't change across graphs
# How to create a custom color scale: https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
myColors <- c("#66c2a5","#d53e4f","#f46d43","#313695","#4575b4","#fee08b","#abdda4","#fdae61","#b2abd2","#9e0142","#74add1","#80cdc1","#e6f598","#5e4fa2")
# Associate a name to each value
names(myColors) <- c("Acartiidae","Calanidae","Calanoida","Centropagidae","Euchaetidae",
					"Oithonidae","Oncaeidae","Temoridae","Corycaeidae","Harpacticoida",
					"Candaciidae","Pontellidae","Sapphirinidae","Oncaea")
					
### Vary the y scale according to season (so the visual comparison is easier across years)
# yaxis_spring <- c(0,8200)
# yaxis_summer <- c(0,10000)
# yaxis_fall <- c(0,15000)

##### 03/08/17: Definite color palette (so it doesn't change everytime...)
#    Acartiidae     Calanidae     Calanoida Centropagidae   Euchaetidae 
#    "#66c2a5"     "#d53e4f"     "#f46d43"     "#313695"     "#4575b4" 
#   Oithonidae     Oncaeidae     Temoridae   Corycaeidae Harpacticoida 
#    "#fee08b"     "#abdda4"     "#fdae61"     "#b2abd2"     "#9e0142" 
#  Candaciidae   Pontellidae Sapphirinidae        Oncaea 
#    "#74add1"     "#80cdc1"     "#e6f598"     "#5e4fa2"

quartz() 
ggplot(ddf3, aes(x = factor(round(size,1)), fill = factor(category))) +
	 	geom_histogram(binwidth=.5, colour = "black", stat = "count") + 
		scale_fill_manual(name = "", values = myColors) + 
		theme_light() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
		xlab("Size (mm)") + ylab("Count")

##### For each unique(ddf3$id) but 2012, compute relative and absolute contrib like above
# i <- unique(ddf3$id)[1] # For testing
# For loop
for(i in unique(ddf3$id)) {

		# Useless message
		message(paste("Doing ", i, sep = ""))
		d <- ddf3[ddf3$id == i,]
		ss <- unique(d$season)
		effort <- length(unique(d$station_id)) # keep samplign effort to normalize the abundances
		
		# Compute contrib
		categories <- unique(d$category)
		# categories
		contrib <- lapply(categories, function(c) {
		
						contrib <- data.frame(class = classes$mid, width = classes$width, low = classes$min, up = classes$max, contrib = NA)
		
						# Fill this new data.frame with a for loop
						for(cl in contrib$class) {		
								# Find
								low <- contrib[which(contrib$class == cl), "low"]
								up  <- contrib[which(contrib$class == cl), "up"]
								# Compute sum of sizes within this interval
								tot <- sum(d[which(d$category == c & d$size < up & d$size >= low),"DW"])
								# normalize by class width
								wid <- contrib[which(contrib$class == cl), "width"]
								norm_tot <- tot/wid
								# and provide
								contrib[which(contrib$class == cl),"contrib"] <- norm_tot
						} # eo for loop
		
						# Return for cbinding later (so it is useless to return the classes as well...)
						return(contrib[,c("contrib")])
	
		}) # eo lapply

		contri <- do.call(cbind, contrib)
		gc() ; rm(contrib)

		# Change colnames
		colnames(contri) <- categories
		classes2 <- cbind(classes, contri)
		### When dealing with absolute contributions: need to correct the differences of sampling effort across cruises (i.e. spring 2015 was under sampled)
		classes2[,c(6:length(classes2))] <- classes2[,c(6:length(classes2))] / effort
		
		### Log tranform the contributions so you see the larger size classes better
		#classes2[,c(6:length(classes2))] <- log1p(classes2[,c(6:length(classes2))])
		# classes2
		
		# head(classes) ; dim(classes)
		# Divide by classes$total_n when relative contrib
		#require("matrixStats")
		#classes2$total <- rowSums(as.matrix(classes2[,c(as.character(categories))]))
		#classes2[,c(as.character(categories))] <- ( classes2[,c(as.character(categories))] / classes2[,c("total")] )*100
		# classes <- na.omit(classes)
		### Time for a ggplot, but first we need to change and melt the ddf a bit
		#maxima <- apply(classes[,c(as.character(categories))], 2, function(x) max(x, na.rm = TRUE))
		#names <- names(maxima[maxima > 1])
		#classes2 <- classes[,c("i", "min", "max", "width", "mid", names)]
		# Melt by omitting last column (= total)
		#last <- length(classes2) - 1
		mclass <- melt(classes2, id.vars = c("i","min","max","width","mid") )
		### Before plotting: replace "Calanoida" by "Unidentified Calanoida"
		#mclass$variable[mclass$variable == "Calanoida"] <- "Unidentified Calanoida"

		#quartz()
		### NOTE to plot absolute contrib: get rid of 'position = "fill" ' in the geom_bar()
		### Use if loops to change the y axis according to the season
		if(ss == "spring") {
			
				plot <- ggplot(data = mclass) + 
  						geom_bar(aes(x = factor(round(mid,2)), y = value, fill = factor(variable)), stat = "identity", colour = "black") + 
						scale_fill_manual(name = "", values = myColors) + ggtitle(i) + 
						theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
  						theme(axis.text.y= element_text(size = 9)) + xlab("Size class (mm)") + ylab("Normalized absolute contributions log(µg/m3.mm)") +
						scale_y_continuous(limits = c(0,37000))
					
				ggsave(plot = plot, filename = paste("plot_contribs_DW_",i,".pdf", sep = ""), dpi = 300, width = 10, height = 5)
				
		} else if (ss == "summer") {
			
				plot <- ggplot(data = mclass) + 
						geom_bar(aes(x = factor(round(mid,2)), y = value, fill = factor(variable)), stat = "identity", colour = "black") + 
						scale_fill_manual(name = "", values = myColors) + ggtitle(i) + 
						theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
						theme(axis.text.y= element_text(size = 9)) + xlab("Size class (mm)") + ylab("Normalized absolute contributions log(µg/m3.mm)") +
						scale_y_continuous(limits = c(0,28000))
				
				ggsave(plot = plot, filename = paste("plot_contribs_DW_",i,".pdf", sep = ""), dpi = 300, width = 10, height = 5)
				
		} else {
			
				plot <- ggplot(data = mclass) + 
						geom_bar(aes(x = factor(round(mid,2)), y = value, fill = factor(variable)), stat = "identity", colour = "black") + 
						scale_fill_manual(name = "", values = myColors) + ggtitle(i) + 
						theme_light() + theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust= 0, size = 9)) +
						theme(axis.text.y= element_text(size = 9)) + xlab("Size class (mm)") + ylab("Normalized absolute contributions log(µg/m3.mm)") +
						scale_y_continuous(limits = c(0,31000))
			
				ggsave(plot = plot, filename = paste("plot_contribs_DW_",i,".pdf", sep = ""), dpi = 300, width = 10, height = 5)
			
		} # eo if else loops
		
		# Clean stuff
		rm(classes2, plot, contri, categories, d)
		gc()
		
} # eo for loop	
	
	
	
	
	
	
	
	