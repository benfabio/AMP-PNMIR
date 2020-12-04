
##### 27/04/2017: R Script to load the mesozooplankton data acquired during the PNMIR database, and analyzed withg the ZooScan
##### Aims to:

#	- Load the data diretcly from EcoTaxa using PostgreSQL requests
#	- Get the objects' taxonomic classification
#	- Compute the objects' abundances and biovolumes while accounting for their children classes 
#		(ex. abund of Copepoda includes objects that have been classified as 'Calanoida' or 'Acartia')


### Latest update: 12/06/2017

library("stringr")
library("reshape2")
library("dplyr")
library("dbplyr")
library("ggplot2")
library("RPostgreSQL")
library("data.tree")


##### Load the functions that were developed by Jean-Olivier Irisson
# http://obs-vlfr.fr/data/files/~zoo/lib_ecotaxa.R
source("http://obs-vlfr.fr/data/files/~zoo/lib_ecotaxa.R")
#!/usr/bin/Rscript
#
# Stub of a package to access EcoTaxa
#
# (c) Copyright 2017 Jean-Olivier Irisson, GNU General Public License v3

# Connect to EcoTaxa for requests
db <- connect_ecotaxa()
proj <- extract_projects(db, "Zooscan Tara Oceans")
# Get rid of project 11
proj <- proj[which(proj$projid %in% c(377,378)),]

# Extract projects (~ 3min)
d <- group_by(proj, projid) %>% do(extract_objects(db, .,
		status = c("V","P"), object_fields = c("major", "minor", "area"), 
		process_field = "particle_pixel_size_µm", acquis_fields = "sub_part", 
		sample_fields = c("tot_vol", "orig_id")) )

### Explore 'd'
class(d)
dim(d)
str(d)
colnames(d)
summary(d)


##### Again load the functions developed by J-O. Irisson to get the taxonomic classification and compute biovolume and abundances at the same time !
#
# Compute concentrations from Ecotaxa
# Specialised function for vlfr data, that does everything in one go
# (c) Copyright 2017 Jean-Olivier Irisson, GNU General Public License v3

# source("/home/zoo/r_scripts/lib_ecotaxa.R")

zoo_objects <- function(db, pattern) {
	
	 				# list projects matching pattern
 				   	projs <- extract_projects(db, pattern)
					#projs <- projs[-which(projs$projid == 11),]
 				  	message("Extracting data from:\n  ", str_c(projs$title, collapse="\n  "))

 				   	# read objects with useful metadata
 				   	o <- projs %>% group_by(projid) %>% do(extract_objects(db, proj=.,
								status = "V", process_field = "particle_pixel_size_mm", object_field = c("major", "minor", "area"), 
								sample_field = c("tot_vol","orig_id"), acquis_field = c("sub_part"))) %>% ungroup()

 					# convert extra fields into numbers
 				   	for (col in c("major", "minor", "area", "tot_vol", "sub_part", "particle_pixel_size_mm")) {
  					 		o[[col]] <- as.numeric(o[[col]])
 				  	} # eo for loop

 				  	### To compute ESD:
					# objs$esd <- 2 * sqrt(objs$area / pi) * objs$particle_pixel_size_mm

 				 	# compute individual particle characteristics
 				   	o <- mutate(o,
   					 		concentration = 1 * sub_part / tot_vol,
   					 		biovolume = concentration * 4/3 * pi * (major * 0.0106/ 2) * (minor * 0.0106 / 2)^2
					) # eo mutate
 					
					# Ref: Forest et al.

 				   	return(o)
					
} # eo zoo_objects


zoo_expand <- function(db, d) {
	
 					# extract a reduced taxonomy for these classes
					taxo <- extract_taxo(db, d$classif_id, rec = T) %>% arrange(id)

 				   	# further simplify the taxonomy: remove all levels with only one child and branch the child to the higher level
 				   	# taxo <- simplify_taxo(taxo)

 				   	# compute abundance of all groups, including children
 				   	total_with_children <- function(id, taxo, d) {
						
   					 						children_ids <- children(id$id, taxo)
   										 	x <- filter(d, classif_id %in% children_ids)
											
											### group_by something else than date..."orig_id" ?
  										  	x <- x %>% 
												 group_by(orig_id) %>%
												 dplyr::summarise( concentration = sum(concentration), biovolume = sum(biovolume) )
											x$id <- id$id
  										  	return(x)
 
					} # eo total_with_children
 
 				   	d <- taxo %>% 
						 group_by(id) %>% 
						 do( total_with_children(., taxo, d) )
						 
 				   	# add zeros (taxon not captured at that date)
 				   	all <- expand.grid(orig_id = unique(d$orig_id), id = unique(d$id))
 				   	d <- left_join(all, d, by = names(all))
 				   	d$concentration[is.na(d$concentration)] <- 0
 				   	d$biovolume[is.na(d$biovolume)] <- 0

 				   	# add taxonomic names
 				   	d <- left_join(d, select(taxo, id, name), by = "id")

 				   	# sort in chronological
 				   	#d <- arrange(d, date, name)

 				   	return(list(d = d, taxo = taxo))
 
} # eo zoo_expand


### Apply those functions to d
zoo <- zoo_objects(db = db, pattern = "Zooscan Tara Oceans 2009 2012 WP2 200")
# Examine results
dim(zoo)
head(zoo)
colnames(zoo)
summary(zoo)
str(zoo)
unique(zoo$orig_id) #  # you want to compute the concentrations and biovolumes according to orig_id, and not date

zoo2 <- zoo_expand(db = db, d = zoo) 
zoo3 <- zoo2$d
dim(zoo3)
head(zoo3)
unique(zoo3$name) # Ok.

zoo3[which(zoo3$concentration == 0 & zoo3$biovolume > 0),]


### Need to cast these melted dataframe to put taxon as columns ! 
# First, need to melt
mzoo <- melt(zoo3, id.vars = c("orig_id","name",'id'))
dim(mzoo)
head(mzoo)
# and now dcast()
#?dcast
zoo_f <- dcast(data = mzoo, orig_id  ~ variable + name, fun.aggregate = sum)
dim(zoo_f)
head(zoo_f)

# GUT !!
summary(zoo_f)
### Check if there are any null rowSums of Biovolumes !!!
rownames(zoo_f) <- c(1:nrow(zoo_f)) 
sums <- rowSums(as.matrix(zoo_f[,c(163:323)]))
ids <- as.numeric(names(sums[which(sums == 0)]))
zoo_f[ids,]


write.csv(zoo_f, "wp2_all_ecotaxa_abund_biovol_21_06_17.csv")

