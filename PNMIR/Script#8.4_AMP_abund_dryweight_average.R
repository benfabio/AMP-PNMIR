
##### 26/07/2017: R Script to compute the the average (+ sd) abundance and dry weight of the main mesozooplankton groups for each season (spring, summer, fall) - Table 1

##### Aims to:
#	- Load the community structure (abundances and then E.biovolumes) 
#	- Compute mean and sd of the mesozooplankton groups' abundances and dry weights
#	- Specify which groups present significant longitudinal variations --> meroplankton vs. open ocean plankton 
#	- Do the same with median and interquartile range (because of non normammy-distributed variables)

### Latest update: 09/08/2017

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

# --------------------------------------------------------------------------------------------------------------------------

##### 1°) Load the data
setwd("/Users/ben_fabio/Desktop/AFB_ParcMarinIroise/data/")
ddf3 <- get(load("ddf3_06_07_17.Rdata"))
# Compute size from major axis
ddf3$size <- ddf3$object_major*0.0106
# Compute abundance/concentration
ddf3$concentration <- 1 * ddf3$acq_sub_part / ddf3$sample_tot_vol
# Add ellipsoidal and spherical biovolumes
ddf3$EBioVol <- ddf3$concentration * (4/3 * pi * (ddf3$object_major * 0.0106 / 2) * (ddf3$object_minor * 0.0106 / 2)^2)
ddf3$SBioVol <- ddf3$concentration * (4/3 * pi *(sqrt((ddf3$object_area*(0.0106)^2)/pi))^3)
# summary(log10(ddf3$EBioVol))
ddf3$logEBioVol <- log10(ddf3$EBioVol)
ddf3$logSBioVol <- log10(ddf3$SBioVol)

# Add corresponding season
ddf3$season <- NA
ddf3[which(ddf3$month %in% c(4,5,6)),"season"] <- "spring"
ddf3[which(ddf3$month == 7),"season"] <- "summer"
ddf3[which(ddf3$month == 10),"season"] <- "fall"
ddf3$id <- paste(ddf3$season, ddf3$year, sep = "_" )

# Keep groups of interest --> those you selected in the RDA
table <- read.table("PNMIR_cruise_data_abund_v_18_07_17.txt", h = T, sep = "\t", dec = ".")
colnames(table) 
cctable <- na.omit(table[,c(1,49:51,53:55,57:80)]) 
# The 30 groups are: 
colnames(cctable[,c(2:length(cctable))]) 
### Re-plot the histogram of absolute abundances for information:
mzoo <- melt(cctable, id.vars = c("ID"))
head(mzoo)
d <- data.frame(mzoo %>% group_by(variable) %>% summarise( total_abund = sum(value) ))
d$total_abund <- round(d$total_abund, 3) 
# Sort by decreasing order of total abund
d <- d[order(d$total_abund, decreasing = T), ]
# force factor levels to keep this order : factor() !
d$variable <- factor(d$variable, levels = d$variable)
ggplot(d) + theme_light() + geom_bar(aes(x = variable, y = total_abund), stat = "identity") + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("Total abundance (ind/m3) - all years and stations")
			
			
##### 2°) From ddf3, compute the mean (+sd) abundances for each group and each of the 3 season
# Sort groups alphabetically for Table
sort(as.character(colnames(cctable[,c(2:length(cctable))])))
 
# To check some categories
summary(ddf3[ddf3$category == "Chaetognatha",c("parent1","parent2")])
summary(ddf3[ddf3$parent1 == "Calanoida",c("category","parent1","parent2")])
summary(ddf3[ddf3$parent2 == "Thecosomata",c("category","parent1","parent2")])

### Compute mean + sd of abundances
mean( ddf3[ddf3$category == "Calanoida" & ddf3$season == "spring","concentration"] )
sd( ddf3[ddf3$category == "Acartiidae" & ddf3$season == "spring","concentration"] )

ggplot() + geom_point(aes(x = object_lon, y = concentration), data = ddf3[ddf3$category == "Calanoida",]) + theme_light() + xlab("Longitude") + ylab("Abundance (ind/m3)")

### Problem of order of magnitude ;)
### Need to compute per station numbnuts: ddf3$station_id
#i <- unique(ddf3$station_id)[3]
res <- lapply(unique(ddf3$station_id), function(i) {
				# Useless message
				message(paste(i, sep = ""))
				d <- ddf3[which(ddf3$station_id == i),]
				# Compute total abund per categories, provide 
				abundances <- lapply(unique(d$category), function(sp) {
									abund <- sum(d[which(d$category == sp),"concentration"])
									return(data.frame(sp = sp, abund = abund))
				}) # eo 2nd lapply
				# Rbind biovolumes
				abund <- do.call(rbind, abundances)
				abund$id <- i
				# Add season:
				abund$season <- unique(d$season)
				# And coordinates
				abund$x <- unique(d$object_lon)
				abund$y <- unique(d$object_lat)
				# Add year
				abund$year <- unique(d$year)
				# Return this
				return(abund)
		} 
		
) # eo first lapply
t <- do.call(rbind, res)
t
colnames(t)

### For Amphipoda: combine Gammaridea and Hyperiidea
levels(t$sp)[c(53,56)] <- 'Amphipoda'
### For Annelida: combine Annelida + larvae
levels(t$sp)[c(24,32)] <- 'Annelida'
### For Bryozoa: combine bryozoans and cyphonautes
levels(t$sp)[c(19)] <- 'Bryozoa'
### For Chaetognatha, combine Chaetognatha and tail
levels(t$sp)[c(4,29)] <- 'Chaetognatha'
### For Cirripedia, combine Cirripedia with cypris & nauplii & cirrus
levels(t$sp)[c(17,15,37)] <- 'Cirripedia'
### For Cladocera, combine Podon, Evadne and Penilia
levels(t$sp)[c(7,28,36)] <- 'Cladocera'
### For Decapoda, combine Decapoda + zoea
levels(t$sp)[c(24,25)] <- 'Decapoda'
### For Echinodermata, combine Echinodermata, Ophiuroidea, Echinoidea
levels(t$sp)[c(31,40,45)] <- 'Echinodermata'
### For Euphausiacea, combine with 'calyptopsis'
levels(t$sp)[c(43)] <- 'Euphausiacea'
### For fish eggs, combine Actinopterygii and eggs
levels(t$sp)[c(8,9)] <- 'Fish'
### For Hydrozoans, combine Hydrozoans + Sipho + Leptothecata + Hydroidolina + Obelia + Scypho + Aglaura + ephyra + Abylidae
levels(t$sp)[c(10,16,37,41,46,50,57,60)] <- 'Hydrozoa'
### For Thecosomata, combine Thecosomata and Limacinidae
levels(t$sp)[c(18)] <- 'Thecosomata'

##### Use facet_wrap to plot the yearly variations of abundances
ggplot(t[which(t$year %in% c(2011,2013,2015) & t$sp %in% c("Calanoida","Acartiidae","Calanidae","Oithonidae")),]) + 
	theme_light() + xlab("Year") + ylab("Abundance (ind/m3)") + scale_y_continuous(limits = c(0,1500)) + 
	geom_boxplot(aes(x = factor(year), y = abund, fill = factor(year)), width = 0.5) + 
	facet_grid(~ sp, scales = "free_y") + scale_fill_manual(values = viridis(3))



### Compute mean + sd of abundances
round(median( t[t$sp == "Calanoida" & t$season == "fall", "abund"]),1)
round(IQR(t[t$sp == "Calanoida" & t$season == "fall", "abund"]),1)
#round(sd( t[t$sp == "Calanoida" & t$season == "summer", "abund"]), 1)

### Check for longitudinal gradient
summary(lm(abund ~ x, data = t[t$sp == "Corycaeidae",]))
# plot
ggplot() + geom_point(aes(x = x, y = abund, fill = factor(season)), data = t[t$sp == "Corycaeidae",], size = 3, pch = 21) +  
		scale_fill_manual(name = "Season", values = viridis(3)) + theme_light() + xlab("Longitude") + ylab("Abundance (ind/m3)")

### Check for seasonal variations
kruskal.test(abund ~ factor(season), data = t[t$sp == "Calanoida",])

### Check for interannual variations (excluding 2012)
kruskal.test(abund ~ factor(year), data = t[t$sp == "Calanoida" & t$year %in% c(2011,2013,2015),])
# 2011
round(median( t[t$sp == "Calanidae" & t$year == 2011, "abund"]),1)
round(IQR(t[t$sp == "Calanidae" & t$year == 2011, "abund"]),1)
# 2013
round(median( t[t$sp == "Calanidae" & t$year == 2013, "abund"]),1)
round(IQR(t[t$sp == "Calanidae" & t$year == 2013, "abund"]),1)
# 2015
round(median( t[t$sp == "Calanidae" & t$year == 2015, "abund"]),1)
round(IQR(t[t$sp == "Calanidae" & t$year == 2015, "abund"]),1)

# Boxplots:
ggplot(t[t$sp == "Thecosomata" & t$year %in% c(2011,2013,2015),], aes(x = factor(year), y = abund, fill = factor(year))) + theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Abundance (ind/m3)") +
scale_y_continuous(limits = c(0,200))


### And total mesozooplankton
res <- lapply(unique(ddf3$station_id), function(i) {
				# Useless message
				message(paste(i, sep = ""))
				d <- ddf3[which(ddf3$station_id == i),]
				# Compute total abund
				abund <- data.frame(tot_abund = sum(d[,"concentration"]) )
				#abund <- do.call(rbind, abundances)
				abund$id <- i
				# Add season:
				abund$season <- unique(d$season)
				# And coordinates
				abund$x <- unique(d$object_lon)
				abund$y <- unique(d$object_lat)
				# Add year
				abund$year <- unique(d$year)
				# Return this
				return(abund)
		} 
		
) # eo first lapply
total <- do.call(rbind, res)
total
colnames(total)
			
### Longitudinal variations			
summary(lm(tot_abund ~ x, data = total)) # p-value: 1.086e-08
# Plot
ggplot() + geom_point(aes(x = x, y = tot_abund, fill = factor(season)), data = total, size = 3, pch = 21) +  
		scale_fill_manual(name = "Season", values = viridis(3)) + theme_light() + xlab("Longitude") + ylab("Total abundance (ind/m3)")

### Seasonal variations
kruskal.test(tot_abund ~ factor(season), data = total) # p-value = 0.8491
### Seasonal averages
round(median( total[total$season == "spring", "tot_abund"]), 1)
round(IQR( total[total$season == "spring", "tot_abund"]), 1)
round(median( total[total$season == "summer", "tot_abund"]), 1)
round(IQR( total[total$season == "summer", "tot_abund"]), 1)
round(median( total[total$season == "fall", "tot_abund"]), 1)
round(IQR( total[total$season == "fall", "tot_abund"]), 1)


### Interannual variations?
kruskal.test(tot_abund ~ factor(year), data = total[which(total$year %in% c(2011,2013,2015)),]) # p-value = 0.002143
# fligner.test(abund ~ factor(year), data = t[which(t$sp %in% groups & t$year %in% c(2011,2013,2015)),])
# Boxplots:
ggplot(total[which(total$year %in% c(2011,2013,2015)),], aes(x = factor(year), y = tot_abund, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Total abundance (ind/m3)")

median( t[which(t$sp %in% groups & t$year == 2011),"abund"], na.rm = T)
median( t[which(t$sp %in% groups & t$year == 2013),"abund"], na.rm = T)
median( t[which(t$sp %in% groups & t$year == 2015),"abund"], na.rm = T)




##### 3°) From ddf3, compute the mean (+sd) dry weight for each group and each of the 3 season, using the relationships form Lehette & Harnandez-Leon et al. (2009)

### General formula: DW (µg) = a.(S)^b

### For copepods: use the relationship that is based on subtropical copepods
	#	DW <- 45.25*(X$area*(0.0106)^2)^1.59
### For Chaetognaths: use the relationship based on...chaetognaths
	#	DW <- 23.45*(X$area*(0.0106)^2)^1.19
### For euphausiids: use the relationship that is based on subtropical euphausiids
	#	DW <- 43.81*(X$area*(0.0106)^2)^1.47
### For crustaceans other than euphausiids & copepods, use the relation for subtropical crustaceans
	#	DW <- 44.78*(X$area*(0.0106)^2)^1.56
### otherwise, use the relationships for general mesozooplankton
	#	DW <- 43.38*(X$area*(0.0106)^2)^1.54

# Identify the classes of ddf3$category on which the relationships
unique(ddf3$category)
groups_copepods <- c("Calanoida","Calanidae","Acartiidae","Candaciidae","Centropagidae","Corycaeidae","Euchaetidae",
					"Sapphirinidae","Temoridae","Oncaeidae","Oithonidae","Harpacticoida")

# For chaetognatha...no need, juste use "Chaetognatha" and "tail"

# For euphausiacea: "calyptopsis"

# For subtropical crustaceans:
crustaceans <- c("Amphipoda","Gammaridea","Hyperiidea","Cirripedia","cypris","nauplii","cirrus","Podon","Evadne","Penilia","Decapoda","zoea")

# For others: general mesozooplankton
general_mesozoo <- c("Oikopleuridae","Limacinidae","Hydrozoa","Siphonophorae","Leptothecata","Hydroidolina","Obelia","Scyphozoa","Aglaura","Aequorea","ephyra","Abylidae",
				"Actinopterygii","egg","Echinodermata","Ophiuroidea","Ophiurida","Echinoidea","Ophiotrix","Bryozoa","cyphonaute","Annelida","larvae","Bivalvia")

### Apply the relationships above, with the vectors of class names, to compute dry weights
# Empty vector
ddf3$DW <- NA
# For copepods
ddf3[which(ddf3$category %in% groups_copepods),"DW"] <- 45.25*(ddf3[which(ddf3$category %in% groups_copepods),"object_area"]*(0.0106)^2)^1.59
# For Chaetognatha
ddf3[which(ddf3$category %in% c("Chaetognatha","tail")),"DW"] <- 23.45*(ddf3[which(ddf3$category %in% c("Chaetognatha","tail")),"object_area"]*(0.0106)^2)^1.19
# For calyptopsis/ euphausiceans
ddf3[which(ddf3$category == "calyptopsis"),"DW"] <- 43.81*(ddf3[which(ddf3$category == "calyptopsis"),"object_area"]*(0.0106)^2)^1.47
# For crustaceans
ddf3[which(ddf3$category %in% crustaceans),"DW"] <- 44.78*(ddf3[which(ddf3$category %in% crustaceans),"object_area"]*(0.0106)^2)^1.56
# For general mesozooplankton
ddf3[which(ddf3$category %in% general_mesozoo),"DW"] <- 43.38*(ddf3[which(ddf3$category %in% general_mesozoo),"object_area"]*(0.0106)^2)^1.54
# Correct DW from volume sampled
ddf3$DW <- (ddf3$DW)*ddf3$concentration
# summary(ddf3$DW)

### Repeat the same code as for the abundances
res <- lapply(unique(ddf3$station_id), function(i) {
				# Useless message
				message(paste(i, sep = ""))
				d <- ddf3[which(ddf3$station_id == i),]
				# Compute total abund per categories, provide 
				weights <- lapply(unique(d$category), function(sp) {
									DW <- sum(d[which(d$category == sp),"DW"])
									return(data.frame(sp = sp, dryweight = DW))
				}) # eo 2nd lapply
				# Rbind biovolumes
				dw <- do.call(rbind, weights)
				dw$id <- i
				# Add season and year
				dw$season <- unique(d$season)
				dw$year <- unique(d$year)
				# And coordinates
				dw$x <- unique(d$object_lon)
				dw$y <- unique(d$object_lat)
				# Return this
				return(dw)
		} 
		
) # eo first lapply
t <- do.call(rbind, res)
colnames(t)
t


### For Amphipoda: combine Gammaridea and Hyperiidea
levels(t$sp)[c(53,56)] <- 'Amphipoda'
### For Annelida: combine Annelida + larvae
levels(t$sp)[c(24,32)] <- 'Annelida'
### For Bryozoa: combine bryozoans and cyphonautes
levels(t$sp)[c(19)] <- 'Bryozoa'
### For Chaetognatha, combine Chaetognatha and tail
levels(t$sp)[c(4,29)] <- 'Chaetognatha'
### For Cirripedia, combine Cirripedia with cypris & nauplii & cirrus
levels(t$sp)[c(17,15,37)] <- 'Cirripedia'
### For Cladocera, combine Podon, Evadne and Penilia
levels(t$sp)[c(7,28,36)] <- 'Cladocera'
### For Decapoda, combine Decapoda + zoea
levels(t$sp)[c(24,25)] <- 'Decapoda'
### For Echinodermata, combine Echinodermata, Ophiuroidea, Echinoidea & Ophiotrix
levels(t$sp)[c(40,45,65,66)] <- 'Echinodermata'
### For Euphausiacea, combine with 'calyptopsis'
levels(t$sp)[c(44)] <- 'Euphausiacea'
### For fish eggs, combine Actinopterygii and eggs
levels(t$sp)[c(8,9)] <- 'Fish'
### For Hydrozoans, combine Hydrozoans + Sipho + Leptothecata + Hydroidolina + Obelia + Scypho + Aglaura + ephyra + Abylidae
levels(t$sp)[c(10,16,37,42,47,51,53,58,61)] <- 'Hydrozoa'
### For Thecosomata, combine Thecosomata and Limacinidae
levels(t$sp)[c(18)] <- 'Thecosomata'


round(median( t[t$sp == "Calanoida" & t$season == "spring", "dryweight"], na.rm = T), 1)
round(IQR( t[t$sp == "Calanoida" & t$season == "spring", "dryweight"], na.rm = T), 1)
### Compute mean + sd of abundances
round(median( t[t$sp == "Calanoida" & t$season == "summer", "dryweight"], na.rm = T), 1)
round(IQR( t[t$sp == "Calanoida" & t$season == "summer", "dryweight"], na.rm = T), 1)
### Compute mean + sd of abundances
round(median( t[t$sp == "Calanoida" & t$season == "fall", "dryweight"], na.rm = T), 1)
round(IQR( t[t$sp == "Calanoida" & t$season == "fall", "dryweight"], na.rm = T), 1)


### Check for longitudinal gradient
summary(lm(dryweight ~ x, data = t[t$sp == "Calanoida",]))
# plot
ggplot() + geom_point(aes(x = x, y = dryweight, fill = factor(season)), data = t[t$sp == "Calanoida",], size = 3, pch = 21) +  
		scale_fill_manual(name = "Season", values = viridis(3)) + theme_light() + xlab("Longitude") + ylab("Dry weight (µg/m3)")


### Check for interannual variations (excluding 2012)
kruskal.test(dryweight ~ factor(year), data = t[t$sp == "Calanoida" & t$year %in% c(2011,2013,2015),])
# Boxplots:
ggplot(t[t$sp == "Calanoida" & t$year %in% c(2011,2013,2015),], aes(x = factor(year), y = dryweight, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Dry weight (µm/m3)") + 
		scale_y_continuous(limits=c(0,5000))

### Check for seasonal variations
kruskal.test(dryweight ~ factor(season), data = t[t$sp == "Calanoida",])


### And for total mesozooplankton?
res <- lapply(unique(ddf3$station_id), function(i) {
				# Useless message
				message(paste(i, sep = ""))
				d <- ddf3[which(ddf3$station_id == i),]
				# Compute total abund
				dryweight <- data.frame(tot_dw = sum(d[,"DW"], na.rm = TRUE) )
				#abund <- do.call(rbind, abundances)
				dryweight$id <- i
				# Add season:
				dryweight$season <- unique(d$season)
				# And coordinates
				dryweight$x <- unique(d$object_lon)
				dryweight$y <- unique(d$object_lat)
				# Add year
				dryweight$year <- unique(d$year)
				# Return this
				return(dryweight)
		} 
		
) # eo first lapply
total <- do.call(rbind, res)
total
colnames(total)

			
### Longitudinal variations			
summary(lm(tot_dw ~ x, data = total)) # 0.1222
# Plot
ggplot() + geom_point(aes(x = x, y = tot_dw, fill = factor(season)), data = total, size = 3, pch = 21) +  
		scale_fill_manual(name = "Season", values = viridis(3)) + theme_light() + xlab("Longitude") + ylab("Total dry weight (µg/m3)")

### Seasonal variations
kruskal.test(tot_dw ~ factor(season), data = total) # p-value = 0.1189

round(median( total[total$season == "spring", "tot_dw"], na.rm = T), 1)
round(IQR( total[total$season == "spring", "tot_dw"], na.rm = T), 1)
round(median( total[total$season == "summer", "tot_dw"], na.rm = T), 1)
round(IQR( total[total$season == "summer", "tot_dw"], na.rm = T), 1)
round(median( total[total$season == "fall", "tot_dw"], na.rm = T), 1)
round(IQR( total[total$season == "fall", "tot_dw"], na.rm = T), 1)

### Interannual variations?
kruskal.test(tot_dw ~ factor(year), data = total[which(total$year %in% c(2011,2013,2015)),])  # 0.0006629
# Boxplots:
ggplot(t[which(t$sp %in% groups & t$year %in% c(2011,2013,2015)),], aes(x = factor(year), y = dryweight, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Abundance (ind/m3)") +
		scale_y_continuous(limits=c(0,200))



##### 07/08/17: Check if the dry weight increase od Calanidae in 2015 could be due to C. finmarchicus replacing C. helgolandicus
summary(ddf3[ddf3$category == "Calanidae" & ddf3$year == 2011,"size"]) # 0.5819  1.4461  2.2763  2.1271  2.7828  3.7089 
summary(ddf3[ddf3$category == "Calanidae" & ddf3$year == 2013,"size"]) # 0.4304  1.2497  2.0702  2.0373  2.7557  3.7694 
summary(ddf3[ddf3$category == "Calanidae" & ddf3$year == 2015,"size"]) # 0.5883  2.1009  2.4974  2.3707  2.7560  3.3761
# Nope...
nrow(ddf3[ddf3$category == "Hydrozoa" & ddf3$year == 2011,]) / length(unique(ddf3[ddf3$year == 2011,"station_id"])) # 18.15789
nrow(ddf3[ddf3$category == "Hydrozoa" & ddf3$year == 2013,]) / length(unique(ddf3[ddf3$year == 2013,"station_id"])) # 14.64706
nrow(ddf3[ddf3$category == "Hydrozoa" & ddf3$year == 2015,]) / length(unique(ddf3[ddf3$year == 2015,"station_id"])) # 41.51613


# --------------------------------------------------------------------------------------------------------------------------

##### 07/08/17: Complete the Table with the phytoplankton counts data ! 
ddf <- read.table("PNMIR_cruise_data_abund_v_07_08_17.txt", h = T, sep = "\t", dec = ".")
colnames(ddf)
str(ddf)
### Need to add season
ddf$season <- NA
ddf[which(ddf$month %in% c(4,5,6)),"season"] <- "spring"
ddf[which(ddf$month == 7),"season"] <- "summer"
ddf[which(ddf$month == 10),"season"] <- "fall"	
### Remove 2010
ddf <- ddf[which(ddf$year != 2010),]

round(median( ddf[ddf$season == "spring", "n_phyto_total"], na.rm = T)/1000, 1)
round(IQR( ddf[ddf$season == "spring", "n_phyto_total"], na.rm = T)/1000, 1)
### Compute median + IQR
round(median( ddf[ddf$season == "summer", "n_phyto_total"], na.rm = T)/1000, 1)
round(IQR( ddf[ddf$season == "summer", "n_phyto_total"], na.rm = T)/1000, 1)
### Compute median + IQR
round(median( ddf[ddf$season == "fall", "n_phyto_total"], na.rm = T)/1000, 1)
round(IQR( ddf[ddf$season == "fall", "n_phyto_total"], na.rm = T)/1000, 1)


### Longitudinal variations			
summary(lm(n_nano ~ x, data = ddf))
# Plot
ggplot() + geom_point(aes(x = x, y = n_nano, fill = factor(season)), data = ddf, size = 3, pch = 21) +  
		scale_fill_manual(name = "Season", values = viridis(3)) + theme_light() + xlab("Longitude") + ylab("Dry weight (µg/m3)")

### Seasonal variations
kruskal.test(n_nano ~ factor(season), data = ddf)

#round(median( ddf[ddf$season == "spring", "n_diato"], na.rm = T), 1)
#round(IQR( ddf[ddf$season == "spring", "n_diato"], na.rm = T), 1)

### Interannual variations?
kruskal.test(n_nano ~ factor(year), data = ddf[which(ddf$year %in% c(2011,2013,2015)),])

# Boxplots:
quartz()
ggplot(ddf[which(ddf$year %in% c(2011,2013,2015)),], aes(x = factor(year), y = n_nano, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Nanoflagellate abundance (ind/m3)")
quartz()
ggplot(ddf[which(ddf$year %in% c(2011,2013,2015)),], aes(x = factor(year), y = n_dino, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Dinoflagellate abundance (ind/m3)")
quartz()
ggplot(ddf[which(ddf$year %in% c(2011,2013,2015)),], aes(x = factor(year), y = n_diato, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Diatom abundance (ind/m3)")
quartz()
ggplot(ddf[which(ddf$year %in% c(2011,2013,2015)),], aes(x = factor(year), y = n_nano, fill = factor(year))) + 
		theme_light() + geom_boxplot(width = 0.5) + xlab("Year") + ylab("Abundance (ind/m3)")




