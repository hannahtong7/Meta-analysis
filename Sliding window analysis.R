install.packages("ggeffects")
install.packages("AICcmodavg")
library(tidyverse) # Data manipulation and figures
library(lme4) # Mixed models
library(ggeffects) # Predictions from models
library(AICcmodavg) # AICc values

#Load data
setwd("~/Documents/temp/meta") # Set your working directory 
cater_timing <- read.csv("cater_timing.csv") # Load the caterpillar data
daily_temp <- read.csv("daily_temp.csv") # Load the temperature data

#data exploration
head(cater_timing)
head(daily_temp[,1:7]) 
colnames(daily_temp)[c(2,ncol(daily_temp))]
x_days <- 58:161
y_temp <- daily_temp[1,2:105] # example using row 1
plot(x_days, y_temp, type="l", xlab="Ordinal date (1 = 1st Jan)",ylab="Temperature (°C)") 

#Take the mean across all dates (columns 2-105) for each row (site_year)
annual_temp <- data.frame(site_year = daily_temp$site_year, 
                          mean_temp = rowMeans(daily_temp[,2:105])) 

#Join the temperature data with the caterpillar data
cater_timing <- left_join(cater_timing, annual_temp, by="site_year")

ggplot(cater_timing, aes(year, mean_temp, col=site))+ 
  geom_point()+ 
  xlab("Year")+
  ylab("Temperature (°C)")+
  theme_bw()

ggplot(cater_timing, aes(year, peak_date, col=site))+ 
  geom_point()+ 
  xlab("Year")+
  ylab("Caterpillar peak date (Ordinal date, 1 = 1st Jan)")+
  theme_bw()

#Does temperature correlate with catepillar peak timing 
ggplot(cater_timing, aes(mean_temp, peak_date, col=site))+ 
  geom_point()+ 
  xlab("Temperature (°C)")+
  ylab("Caterpillar peak date (Ordinal date, 1 = 1st Jan)")+
  theme_bw()

cater_start_col <- seq(2,105,7) 
# list of values from 2 to 105 in increments of 7 - weekly start date intervals 
cater_duration <- seq(7,105,7) 
# list of values from 14 to 105 in increments of 7 - duration of 1 week or more increasing by 1 week
cater_windows <- data.frame(start_col=rep(cater_start_col,1,each=length(cater_duration)),
                            duration=rep(cater_duration,length(cater_start_col)))
# this repeats every start date for the number of durations there are and vice versa to pair all options
cater_windows$end_col <- cater_windows$start_col+cater_windows$duration-1 
# working out the end column, -1 is included because the start date is included in the window
cater_windows <- cater_windows[-which(cater_windows$end_col>105),]
# removing any windows that extend past the available data
# Give the windows an ID so it's clear which window they test
cater_windows$window_ID <- paste0(colnames(daily_temp)[cater_windows$start_col],
                                  "_",cater_windows$duration,"days") 
# Here we've taken the column name for the start date of the window and combined it with the duration of the window 
# The ID now says which ordinal date the window will start on and how long it is in days
# create and empty plot with x axis for the number of days of temp data and y for each window
plot(NA, xlim=c(0,105), ylim=c(1,nrow(cater_windows)), xlab="Column number", ylab="Different windows") 
# Use a loop to plot each window
for(i in 1:nrow(cater_windows)){ 
  points(cater_windows[i,c("start_col","end_col")], c(i,i), type="l") 
}

#run a null model for comparison - a model not including temperature 
cater_base_mod <- lmer(peak_date ~ 1 + (1|site) + (1|year), cater_timing, REML=F)
summary(cater_base_mod)

#make a data frame for the results of the sliding window models
# Make an empty data frame that the results will go into
cater_slidwin <- data.frame(matrix(NA, ncol=6, nrow=nrow(cater_windows)))

# Name the columns
colnames(cater_slidwin) <- c("window_ID", "start_date", "end_date", "deltaAICc", "temp_coef", "temp_SE")

#Use a loop to run a model for each window
for(i in 1:nrow(cater_windows)){
  
  #Extract relevant temperature data
  temp_dat <- data.frame(site_year = daily_temp$site_year,  
                         window_temp = rowMeans(daily_temp[,cater_windows$start_col[i]
                                                           :cater_windows$end_col[i]]))
  
  # Join temperature and caterpillar data
  cater_windtemp <- left_join(cater_timing, temp_dat, by="site_year")
  
  # Run the model 
  mod <- lmer(peak_date ~ 1 + window_temp + (1|site) + (1|year), cater_windtemp, REML=F)
  
  # Store the relevant information
  cater_slidwin$window_ID[i] <- cater_windows$window_ID[i] 
  cater_slidwin$start_date[i] <- cater_windows$start_col[i]+56 
  #56 is added to column number so it becomes the ordinal date
  cater_slidwin$end_date[i] <- cater_windows$end_col[i]+56 
  cater_slidwin$deltaAICc[i] <- AICc(mod)-AICc(cater_base_mod)  
  cater_slidwin$temp_coef[i] <- summary(mod)$coefficients[2,1]  
  cater_slidwin$temp_SE[i] <- summary(mod)$coefficients[2,2]
  
  # remove elements that were specific to this run of the sliding window
  rm(temp_dat, cater_windtemp, mod)
  
}

#Look at the result of the sliding window analysis 
View(cater_slidwin)

#Plot windows by thier delta AIC
# blank plot with axis from min to max values for dates and AICc
plot(NA, xlim=c(min(cater_slidwin$start_date),max(cater_slidwin$end_date)),
     ylim=c(min(cater_slidwin$deltaAICc),max(cater_slidwin$deltaAICc)), 
     xlab="Ordinal Date (1 = 1st Jan)", ylab="deltaAICc") 

# use loop to draw the lines
for(i in 1:nrow(cater_slidwin)){
  points(c(cater_slidwin$start_date[i],cater_slidwin$end_date[i]),
         c(cater_slidwin$deltaAICc[i],cater_slidwin$deltaAICc[i]),type="l") 
} 

# line at 2 AICc above the lowest
abline(h=(min(cater_slidwin$deltaAICc)+2), col="red", lty="dashed")

#What effect does temperature have on the catepillar peak timing
# Row number for the window with the lowest AIC
cater_wind_row <- which(cater_slidwin$deltaAICc==min(cater_slidwin$deltaAICc))

# ID for the best window
cater_wind_ID <- cater_slidwin$window_ID[cater_wind_row] 

#The row number is the same in the cater_windows and cater_slidwin dataframes
cater_wind_row

which(cater_windows$window_ID==cater_wind_ID)
# Mean temperature during the identified window
cater_best_temp <- data.frame(site_year = daily_temp$site_year,  
                              best_temp = rowMeans(daily_temp[,cater_windows$start_col[cater_wind_row]:cater_windows$end_col[cater_wind_row]])) 

# Join with the caterpillar data
cater_timing <- left_join(cater_timing, cater_best_temp, by="site_year")

# Run the same model as before but with REML=TRUE
cater_mod <- lmer(peak_date~best_temp+(1|year)+(1|site), cater_timing, REML=TRUE) 
summary(cater_mod)
# Store the slope and confidence intervals 
cater_temp_coef <- summary(cater_mod)$coefficients["best_temp","Estimate"] 
cater_temp_confint <- confint(cater_mod)["best_temp",]

# This shows the mean and 95% confidence intervals for the slope in units of days per°C
cater_temp_coef 
cater_temp_confint

# Use ggpredict to get estimates of caterpillar peak timing across the range of temperatures  included in the data
pred_cater <- ggpredict(cater_mod, "best_temp")

#Plot the mean prediction and CIs with the data
ggplot(pred_cater, aes(x,predicted))+ 
  geom_line(lwd=1.2)+ 
  geom_point(data=cater_timing, aes(best_temp, peak_date))+ 
  geom_ribbon(data=pred_cater, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.25)+
  xlab("Temperature (°C)")+ 
  ylab("Caterpillar peak date (Ordinal date, 1 = 1st Jan)")+ 
  theme_bw()


#Section 2
bird_timing<- read.csv("bird_timing.csv")
head(bird_timing)

bird_timing <- left_join(bird_timing, annual_temp, by="site_year")

ggplot(bird_timing, aes(year, hatch_date, col=site))+ 
  geom_point()+ 
  xlab("Year")+
  ylab("Blue tit hatch date (Ordinal date, 1 = 1st Jan)")+
  theme_bw()
ggplot(bird_timing, aes(mean_temp, hatch_date, col=site))+ 
  geom_point()+ 
  xlab("Temperature (°C)")+
  ylab("Blue tit hatch date (Ordinal date, 1 = 1st Jan)")+
  theme_bw()
bird_start_col <- seq(2,105,7) 

bird_duration <- seq(7,105,7) 

bird_windows <- data.frame(start_col=rep(bird_start_col,1,each=length(bird_duration)),
                           duration=rep(bird_duration,length(bird_start_col)))

bird_windows$end_col <- bird_windows$start_col+bird_windows$duration-1 

bird_windows <- bird_windows[-which(bird_windows$end_col>105),]

bird_windows$window_ID <- paste0(colnames(daily_temp)[bird_windows$start_col],
                                 "_",bird_windows$duration,"days") 

plot(NA, xlim=c(0,105), ylim=c(1,nrow(bird_windows)), xlab="Column number", ylab="Different windows") 
for(i in 1:nrow(bird_windows)){ 
  points(bird_windows[i,c("start_col","end_col")], c(i,i), type="l") 
}
bird_base_mod <- lmer(hatch_date ~ 1 + (1|site) + (1|year), bird_timing, REML=F)
summary(bird_base_mod)
bird_slidwin <- data.frame(matrix(NA, ncol=6, nrow=nrow(bird_windows)))
colnames(bird_slidwin) <- c("window_ID", "start_date", "end_date", "deltaAICc", "temp_coef", "temp_SE")
for(i in 1:nrow(bird_windows)){
  
  temp_dat <- data.frame(site_year = daily_temp$site_year,  
                         window_temp = rowMeans(daily_temp[,bird_windows$start_col[i]
                                                           :bird_windows$end_col[i]]))
  
  bird_windtemp <- left_join(bird_timing, temp_dat, by="site_year")
  
  mod <- lmer(hatch_date ~ 1 + window_temp + (1|site) + (1|year), bird_windtemp, REML=F)
  
  bird_slidwin$window_ID[i] <- bird_windows$window_ID[i] 
  bird_slidwin$start_date[i] <- bird_windows$start_col[i]+56 
  bird_slidwin$end_date[i] <- bird_windows$end_col[i]+56 
  bird_slidwin$deltaAICc[i] <- AICc(mod)-AICc(bird_base_mod)  
  bird_slidwin$temp_coef[i] <- summary(mod)$coefficients[2,1]  
  bird_slidwin$temp_SE[i] <- summary(mod)$coefficients[2,2]
  
  rm(temp_dat, bird_windtemp, mod)
  
}

View(bird_slidwin)
plot(NA, xlim=c(min(bird_slidwin$start_date),max(bird_slidwin$end_date)),
     ylim=c(min(bird_slidwin$deltaAICc),max(bird_slidwin$deltaAICc)), 
     xlab="Ordinal date (1 = 1st Jan)", ylab="deltaAICc") 

for(i in 1:nrow(bird_slidwin)){
  points(c(bird_slidwin$start_date[i],bird_slidwin$end_date[i]),
         c(bird_slidwin$deltaAICc[i],bird_slidwin$deltaAICc[i]),type="l") 
} 

abline(h=(min(bird_slidwin$deltaAICc)+2), col="red", lty="dashed")

bird_wind_row <- which(bird_slidwin$deltaAICc==min(bird_slidwin$deltaAICc))

bird_wind_ID <- bird_slidwin$window_ID[bird_wind_row] 

bird_wind_row
which(bird_windows$window_ID==bird_wind_ID)
bird_best_temp <- data.frame(site_year = daily_temp$site_year,  
                             best_temp = rowMeans(daily_temp[,bird_windows$start_col[bird_wind_row]:bird_windows$end_col[bird_wind_row]])) 

bird_timing <- left_join(bird_timing, bird_best_temp, by="site_year")

bird_mod <- lmer(hatch_date~best_temp+(1|year)+(1|site), bird_timing, REML=TRUE) 
summary(bird_mod)

bird_temp_coef <- summary(bird_mod)$coefficients["best_temp","Estimate"] 
bird_temp_confint <- confint(bird_mod)["best_temp",]

bird_temp_coef 
bird_temp_confint
pred_bird <- ggpredict(bird_mod, "best_temp")

ggplot(pred_bird, aes(x,predicted))+ 
  geom_line(lwd=1.2)+ 
  geom_point(data=bird_timing, aes(best_temp, hatch_date))+ 
  geom_ribbon(data=pred_bird, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.25)+
  xlab("Temperature (°C)")+ 
  ylab("Blue tit hatch date (Ordinal date, 1 = 1st Jan)")+ 
  theme_bw()
