
##### 28/04/2017: R Script to examine the Sardines capture data from the PNMI © Fabio Benedetti, OOV-UMS, AFB
##### Aims to:

#	- Examine the monthly landings of Sardines that were caught by the Bolincheurs in the Parc Marin of the Iroise Sea.
#	- Those data were retreived in the frame of the fishing indicators for the Parc Marin of the Iroise Sea (IFREMER data).
#	- Time series analysis

### Latest update: 11/08/2017

library("stringr")
library("reshape2")
library("ggplot2")
library("RColorBrewer")
library("oce")
require("dplyr")
library("akima")


##### 28/04/2017: Explore Sardines data, not for sharing. 
dir()
sar <- read.csv("prod_mensuelle_sardine_bolinche.csv", h = T, sep = ";", dec = ",")
head(sar) # OK
str(sar)

ships <- read.csv("nb_declared_ships.csv", h = T, sep = ";", dec = ",")
head(ships)
str(ships)
ships <- ships[,c(1:3)]
ships


### Per month/year, compute fishing effort (number of boats since only info available), average capture weight (kg)
sar$date <- paste(sar$Month, sar$Year, sep = "_")
ships$date <- paste(ships$Month, ships$Year, sep = "_")

sarsar <- data.frame(sar %>%
		 	group_by(date) %>%
	 	 	summarise(avg_weight = mean(weight), sd_weight = sd(weight), sum_weight = sum(weight) ))

# Add sampling effort (number of ships) per date
sarsar$effort <- NA
for(d in sarsar$date) {
		message( paste("Doing ", d, sep = "") )
		effort <- ships[which(ships$date == d),"n_ships"]
		sarsar[which(sarsar$date == d),"effort"] <- effort
} # eo for loop
### OK

### Add month separately
sarsar$month <- NA
sarsar$year <- NA
for(i in 1:nrow(sarsar)) {
		m <- do.call(rbind,strsplit(as.character(sarsar[i,"date"]), split = "_"))[,1]
		sarsar[i,"month"] <- m
		y <- do.call(rbind,strsplit(as.character(sarsar[i,"date"]), split = "_"))[,2]
		sarsar[i,"year"] <- y
} # eo for loop


# Check bias due to sampling effort
summary(lm(avg_weight ~ effort, data = sarsar)) # R-squared = -0.01041; p-value: 0.8387 ### Same with log scale ; no correlation between total landings of Sardines and number of fishing ships
summary(lm(sum_weight ~ effort, data = sarsar)) # R-squared = 0.3169; p-value: 2.108e-09 

# Plot
quartz()
ggplot() + geom_point(aes(x = effort, y = log(sum_weight)), data = sarsar, pch = 21, colour = "black", fill = "grey70", size = 3) + 
		   xlab("Number of fishing boats") + ylab("Total capture weight (kg)") + theme_light() 
		   		   
quartz() 
ggplot(sar, aes(x = factor(Month), y = weight)) + theme_light() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Month") + ylab("Sardine landings (kg)")

quartz() 
ggplot(sarsar) + geom_point(aes(x = factor(date), y = sum_weight, size = effort)) + theme_light() + 
		xlab("Date") + ylab("Sardine landings (kg)") + theme(axis.text.x = element_text(angle = 90, hjust = 1))


positions <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
quartz() 
ggplot(sarsar, aes(x = factor(month), y = effort)) + geom_boxplot(width = 0.5, fill = "grey70") + theme_light() + 
		xlab("Month") + ylab("Number of declared ships") + scale_x_discrete(limits = positions)
		
### Plot landings divided by effort:
# Because : 
summary(lm(sum_weight ~ effort, data = sarsar)) # R-squared = 0.3169; p-value = 2.108e-09
# quartz()
plot <- ggplot() + geom_point(aes(x = effort, y = (sum_weight)), data = sarsar, pch = 21, colour = "black", fill = "grey70", size = 3) + 
		   geom_smooth(aes(x = effort, y = (sum_weight)), data = sarsar, colour = "black", method = 'lm') + 
		   xlab("Number of declared ships") + ylab("Total sardine landings (kg)") + theme_light()

ggsave(plot = plot, filename = "landings_monthly.pdf", dpi = 300, width = 7, height = 5)

# quartz() 
plot <- ggplot(sarsar, aes(x = factor(month), y = (sum_weight/effort))) + geom_boxplot(width = 0.5, fill = "grey70") + theme_light() + 
		xlab("Month") + ylab("Total sardine landings per unit effort\n(kg/number of declared ships)") + scale_x_discrete(limits = positions)

ggsave(plot = plot, filename = "boxplots_landings_CPUE_monthly.pdf", dpi = 300, width = 7, height = 5)
		
sarsar$corr_weight <- sarsar$sum_weight / sarsar$effort
for(m in c(1:12)) {
		message(m)
		message(median(sarsar[sarsar$month == m,"corr_weight"]))
		message(IQR(sarsar[sarsar$month == m,"corr_weight"]))
}

# Yearly:
sarsar <- sar %>%
		 group_by(Year) %>%
	 	 summarise(avg_weight = mean(weight), sd_weight = sd(weight) ) 	 
# Compute annual fishing effort		 
annual_shiping_effort <- ships %>% group_by(Year) %>% summarise(avg_effort = mean(n_ships, na.rm = T) )
# Add sampling effort (number of ships) per date
sarsar$effort <- NA
for(y in sarsar$Year) {
		 	message( paste("Doing ", y, sep = "") )
		 	effort <- annual_shiping_effort[which(annual_shiping_effort$Year == y),"avg_effort"]
		 	sarsar[which(sarsar$Year == y),"effort"] <- effort
} # eo for loop

quartz() 
ggplot(sar, aes(x = factor(Year), y = weight)) + theme_bw() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Month") + ylab("Capture weight (kg)")
### Higher production in 2009 and 2010
# Sampling effort ?
summary(lm(avg_weight ~ effort, data = sarsar)) ### Still no effect


# Seasonal: 
sar$season <- NA
sar[which(sar$Month %in% c(3:5)),"season"] <- "spring"
sar[which(sar$Month %in% c(6:8)),"season"] <- "summer"
sar[which(sar$Month %in% c(9:11)),"season"] <- "fall"
sar[which(sar$Month %in% c(1,2,12)),"season"] <- "winter"

sarsar <- sar %>%
		 group_by(factor(season)) %>%
	 	 summarise(avg_weight = mean(weight), sd_weight = sd(weight) ) 
# Or through boxplots 
quartz() 
ggplot(sar, aes(x = factor(season), y = weight)) + theme_bw() + geom_boxplot(width = 0.5, fill = "grey70") + xlab("Season") + ylab("Capture weight (kg)")


### AND now, a time series: 
# firest, need to order by date: arrange()
library("dplyr")
head( arrange(sarsar, month, year) ) # ok
sar_ts <- arrange(sarsar, year, month)
sar_ts$date <- as.character(paste(1, sar_ts$month, sar_ts$year, sep = "-"))
require("lubridate")
sar_ts$date <- lubridate::dmy(x = sar_ts$date)

quartz()
ggplot() + geom_point(aes(x = date, y = (sum_weight/effort)), data = sar_ts) + theme_light()
### Clear seasonal effect

### For the 2011-2015 period:
sar_ts2 <- sar_ts[sar_ts$year >= 2011,]

### Explore with 'pastecs'
library("pastecs")
# First test, fit a polynomial of degree 6:
#pol6 <- lm((sum_weight/effort) ~ poly(date, degree = 6), data = sar_ts)
#str(pol6)
#POL6 <- cbind(sar_ts, values = pol6$fitted.values, residuals = pol6$residuals)
#quartz()
#ggplot() + 
  #geom_point(aes(x = date, y = (sum_weight/effort) ), sar_ts, colour = "black") +
  #geom_line(aes(x = date, y = values), data = POL6 , colour = "blue") + 
  #geom_line(aes(x = date, y = residuals), data = POL6 , colour = "darkgreen") + 
  #scale_y_continuous(name = "Total sardine landings per unit effort\n(kg/number of declared ships)") + theme_light()
  
# Works, but looks horrible

# To deduce the optimal window size, we will compute the turnogram of the data, with the function turnogram. Use a variable of marbio that you would have retained using the Escoufier criterion.
# sar_ts2 <- sar_ts[sar_ts$year >= 2011,]
trend.test((sar_ts$sum_weight/sar_ts$effort)) # S = 192140, p-value = 0.0001252 ; rho = -0.3881154
trend.test((sar_ts2$sum_weight/sar_ts2$effort)) # S = 41038, p-value = 0.04693 ; rho = -0.2623581

quartz()
ggplot(data = sar_ts2) + geom_point(aes(x = date, y = (sum_weight/effort))) + geom_smooth(aes(x = date, y = (sum_weight/effort)), colour = "red") +
		theme_light() + ylab("Total sardine landings per unit effort\n(kg/number of declared ships)")
		
#?turnogram
t <- turnogram(ts(sar_ts2$sum_weight/sar_ts2$effort))
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
plot( decaverage(sar_ts2$sum_weight/sar_ts2$effort, order = 12, times = 1) )
plot( decaverage(sar_ts$sum_weight/sar_ts$effort, order = 1, times = 1)  )
plot( decaverage(sar_ts2$sum_weight/sar_ts2$effort, order = 4, times = 1) )


### ggploting
str( decaverage(sar_ts$sum_weight/sar_ts$effort, order = 12, times = 1) )
decav <- decaverage(sar_ts2$sum_weight/sar_ts2$effort, order = 12, times = 1)
dim(decav$series)
# Make a ddf out of it
decav <- data.frame(decav$series)
# Add raw data
decav$raw <- (sar_ts2$sum_weight/sar_ts2$effort)
# Add dates
decav$date <- sar_ts2$date
decav$year <- sar_ts2$year

# quartz()
plot <- ggplot() + geom_line(aes(x = date, y = raw), decav, colour = "black", size = .5) +
			geom_point(aes(x = date, y = raw), decav, fill = "grey70", colour = "black", pch = 21, size = 2) +
  			geom_line(aes(x = date, y = filtered), data = decav , colour = "firebrick3", size = 1) + 
  		  	scale_y_continuous(name = "Total sardine landings per unit effort\n(kg/number of declared ships)") + xlab("Date") + theme_light()
  
ggsave(plot = plot, filename = "sardines_restrict_time_series_decaverage.pdf", dpi = 300, height = 4, width = 7)
			
### Re-run a trend.test on the filtered data
trend.test(decav$filtered) # S = 272800, p-value < 2.2e-16 ; rho = -0.9708558
# When retsricting to 2011-2015 time period: S = 58556, p-value < 2.2e-16; rho = -0.8012243 


## Let's weight it : 
plot(decaverage(sar_ts$weight, order = 2, times = 2, weights = c(1,2,4,2,1)))
RA <- decaverage(sar_ts$weight, order = 2, times = 2, weights = c(1,2,4,2,1))
plot(RA)
str(RA)
head(RA)
## values and residuals are stored as components ; let's extract those with : 
# extract(RA, component = "filtered")
# extract(RA, component = "residuals")
fil <- extract(RA, component = "filtered")
res <- extract(RA, component = "residuals")
# Joining in a data.frame for ggploting : 
X <- data.frame(date = sar_ts$date, weight_capt = sar_ts$weight, values = fil, residuals = res)

quartz()
ggplot(data = X) +
  geom_point(aes(x= date, y= weight_capt), colour = "black") +
  geom_path(aes(x = date, y = values), colour = "red") +
  geom_linerange(aes(x = date, ymin = values, ymax = residuals), colour = "red", alpha = 1/2) + theme_linedraw()

trend.test(X$values)     # rho = -0.284 ; p-value < 2.2e-16 ; alternative hypothesis: true rho is not equal to 0
trend.test(X$residuals)  # rho = 0.02494 ; p-value = 0.02641 ; alternative hypothesis: true rho is not equal to 0



### And try some eigen vectors filtering
# The function is decevf() which has two arguments you will be interested in:
# lag : the number of repetitions of the series
# axes : which PCA axes are used to predict the smoothed series
# How do you choose the lag argument? Choose a sensible lag for weight_capt
## --> autocorrelogram to see the variation scale 
acf(ts(sar_ts$weight))
## --> a sensible lag for temperature would be ???

# Compute the EVF of sar_ts$weight_capt using axes 1 and 2. Inspect the resulting object using
# plot and str.
#EVF2 <- decevf(sar_ts$weight_capt, lag = 5, axes = c(1,2)) # to keep the second PC too
EVF3 <- decevf(sar_ts$weight , lag = 5, axes = c(1)) # PC1 only is enough
# str(EVF3)
plot(EVF3)
## The predictions are stored in the "filtered" component
# Extract successively the predictions using axes 1, 2, and 3. Plot all of them and the original data on the same graph. How would you interpret this graph ?
evf3 <- extract(EVF3, component = "filtered")
E <- data.frame(Date = sar_ts$date, weight_capt = sar_ts$weight, EVF = evf3)
quartz()
ggplot(data = E) + geom_point(aes(x = Date, y = weight_capt), colour = "black") + geom_line(aes(x = Date, y = EVF), colour = "blue") + theme_linedraw()






