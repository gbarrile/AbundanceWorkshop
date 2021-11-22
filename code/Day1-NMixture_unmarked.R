
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# November 29, 2021

# Day 1 of workshop: Statistical Methods for Estimating Abundance in Ecology

# Code to fit 'Closed Population N-mixture Model' in the 'unmarked' package

# Addressing the question: "How does wildfire influence salamander abundance?"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# install the packages that we will use:
# install.packages("unmarked")
# install.packages("reshape")

# make sure packages are installed and active
require(unmarked)
require(reshape)

# citing the unmarked package
citation("unmarked")

# read-in the salamander count data
df <- read.csv("data/Salamander_Wildfire.csv")

# so now we have our data stored as 'df'
# check out our data
head(df)

# how many transects did we survey?



# order the data frame by transect
df <- with(df,df[order(Transect),])
df



# 1.
# format count data for unmarked
m <- melt(df,id.var=c("Transect","Survey"),measure.var="Count")
head(m)

y=cast(m, Transect ~ Survey)

K <- ncol(y)

C <- as.matrix(y[,2:K])

C # each transect is a row (10 rows)
# each column indicates the survey at each site (three surveys at each transect)

# Does our matrix of counts match our field data?
# Check: how many salamanders were observed at Transect 2 during Survey 2?
head(df)
# Now using the counts matrix, how many salamanders were observed at Transect 2 during Survey 2?
C
# Do the number of salamanders match between the field data and counts matrix?



# 2.
# format covariates for abundance and detection probability
# hypothesized that percent burned area influences abundance
# hypothesized that time of day influences detection probability

# Site-level covariates versus Observational covariates

# Site covariates do not change during each survey at a given site
# i.e., one value for a site covariate at each transect
# e.g., ten values for percent burned area, one value for each transect 

# We are interested in the abundance at a given transect, which we are assuming
# does not change over our three surveys at that transect. Therefore, if we
# want to estimate transect abundance as a function of percent burned area,
# then the percent burned area also cannot change over our three surveys at that
# given transect

# Observation covariates can change during each survey at a given site
# Observation covariate data should match survey data (i.e., matrix of counts)
# We can have a different detection probability for each survey
# e.g., different time of day for each survey.


# Format percent burned area as a covariate on abundance (as a 'site' covariate)
burned <- unique(df[,c("Transect","Burned")])
burned <- as.matrix(burned[,"Burned"])
burned

# Do the values for percent burned area match our field data?
# Check: what was the percent burned area at Transect 2?
head(df)
# Now with the percent burned values, what was the percent burned area at Transect 2?
burned
# Do they match?



# Format time of day as a covariate on detection (as an 'observation' covariate)
m<-melt(df,id.var=c("Transect","Survey"),measure.var="Time")
head(m)
y=cast(m, Transect ~ Survey)

K <- ncol(y)

time <- as.matrix(y[,2:K])


# Do the values for time match our field data?
# Check: what time was Survey 3 conducted at Transect 2?
head(df)
# Now with the time values, what time was Survey 3 conducted at Transect 2??
time
# Do they match?



# remove objects that we don't need
rm(m,y,K)


# Here are the data for analysis:
C      # matrix of survey data (salamander counts)
burned # percent burned area at each transect (we think it might influence abundance)
time   # time recorded during every survey (we think it might influence detection)

# important to remember that we can also include site covariates on detection

# Let's examine the functions within 'unmarked' that we will use
?pcount
?unmarkedFramePCount



# Input data into an 'unmarked data frame'
umf <- unmarkedFramePCount(
  y=C,                                   # Counts matrix
  siteCovs= data.frame(burned = burned), # Site covariates
  obsCovs = list(time = time))           # Observation covariates

# look at the summary our dataframe
summary(umf)




# Fit Closed Binomial N-mixture Model
?pcount
# Let's first fit a null model
m1 <- pcount(~1 ~1, data=umf, K=130) # don't worry about K right now
summary(m1)
# Extract abundance at each transect
ranef(m1)
meanEst <- bup(ranef(m1))
# compare this to the maximum count at each transect
maxCount <- aggregate(Count ~ Transect, data=df, FUN=max)
maxCount$mean_estimate <- meanEst

plot(maxCount$Transect, maxCount$mean_estimate, pch=16, col="blue",
     ylab="Salamander Count / Abundance", xlab="Transect")
points(maxCount$Transect, maxCount$Count, pch=16, col="black")
legend(6.5, 33, c("Mean predicted abundance", "Maximum count"),
       col=c("blue", "black"), pch = 16)




# Fit Closed Binomial N-mixture Model

# linear model for p (detection) follows first tilde, 
# then comes linear model for abundance (may see abundance ref. to as lambda for this model)
m1 <- pcount(~time ~burned, data=umf, K=130)
# Here 'K' is the upper summation limit for the summation over the random effects
# in the integrated likelihood. In unmarked, the default choice of K is the maximum
# observed count plus 100. That should normally be adequate, but you can play with K
# to determine the sensitivity of estimates to the value set for K.
max(C) # max observed count was 30...30+100=130 for K

# In the model (m1) above,
# detection was modeled as a function of what?
# and abundance was modeled as a function of what?

# Model-fitting function 'pcount' stands for 'point count' as this model can be
# employed with point count data. However, applications of the binomial mixture
# model are not restricted to point count data, as we can see with our use
# of counts along transects.


# check out summary of the model
summary(m1)

# Positive or negative relationship with abundance and percent burned area? 
# Is the relationship significant?

# Positive or negative relationship with detection and time of day? 
# Is the relationship significant?



# We can choose alternate models for abundance other than the Poisson distribution

# Negative Binomial 
m2 <- pcount(~time ~burned, data=umf, mixture="NB", K=130)

# Zero-inflated Poisson
m3 <- pcount(~time ~burned, data=umf, mixture="ZIP", K=130)

# Compare models using AIC
cbind(AIC.P=m1@AIC, AIC.NB=m2@AIC, AIC.ZIP=m3@AIC)

# Poisson has lowest AIC. Why might this be? Discuss amongst the group.


# Predictions of abundance at specified values of percent burned area, say 0, 30, and 60)
newdat <- data.frame(burned=c(0, 30, 60))
predict(m1, type="state", newdata=newdat, append = T)
# 'predict' here uses delta rule to compute SEs and 95% CIs


# Predictions of p (detection probability) for values of time (e.g., 7am, 9am, and 11am)
newdat <- data.frame(time=c(7,9,11))
predict(m1, type="det", newdata=newdat, append = T)


# Visualize covariate relationships

# For abundance, predict to a new dataframe with 
# a suitable range for % burned area values
range(burned)
newdat <- data.frame(burned=seq(0, 60, length.out = 40))
pred.lam <- predict(m1, type="state", newdata=newdat, appendData = TRUE)

# For detection, predict to a new dataframe with 
# a suitable range for time of day values
range(time)
newdat <- data.frame(time=seq(7, 12, length.out = 20))
pred.det <- predict(m1, type="det", newdata=newdat, appendData = TRUE)


# plot abundance relationship with percent burned area
min(pred.lam$lower)
max(pred.lam$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.lam$burned, y = pred.lam$Predicted, pch=16, 
     ylab = "Salamander Abundance",
     xlab = "% Burned Area", cex.lab=1.5, cex.axis=1.2, col="darkgray", ylim=c(0,65))
box(lwd = 4, col = 'black')
lines(pred.lam$burned, pred.lam$Predicted, lwd=8, col="blue")
lines(pred.lam$burned, pred.lam$lower, lwd=4, lty=2, col="black")
lines(pred.lam$burned, pred.lam$upper, lwd=4, lty=2, col="black")


# What is the mean predicted abundance when percent burned area is 12%
# Hint: you may need to use code from earlier in the script



# plot detection relationship with time of day
min(pred.det$lower)
max(pred.det$upper)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.det$time, y = pred.det$Predicted, pch=16, 
     ylab = "Detection Probability",
     xlab = "Time of Day (24 hr clock)", cex.lab=1.5, cex.axis=1.2, 
     col="darkgray", ylim=c(0.45,1))
box(lwd = 4, col = 'black')
lines(pred.det$time, pred.det$Predicted, lwd=8, col="blue")
lines(pred.det$time, pred.det$lower, lwd=4, lty=2, col="black")
lines(pred.det$time, pred.det$upper, lwd=4, lty=2, col="black")


# What is the mean predicted detection for surveys conducted at 10am?
# Hint: you may need to use code from earlier in the script



# Extract abundance at each transect
ranef(m1)
# Site-specific abundance is a random effect (Ni)
# ranef() function obtains estimates of these random effects via
# applying Bayes' rule for a conditional estimate of abundance (Ni)
# given the model parameters and the observed count data (i.e., the function
# obtains the best unbiased prediction [BUP] of the random effects based on 
# the posterior distribution of Ni)





# code to explore on your own if you have time
plot(ranef(m1), xlim=c(0,35))


# goodness-of-fit
require(AICcmodavg)

m1.gof <- Nmix.gof.test(m1, nsim = 100) # nsim is 100 only for illustrative purposes,
# you should probably run at least 1000 bootstrap replicates in your analysis

# p-value here suggests that we fail to reject the null (CANNOT conclude that
# the observed data are statistically different from the expected values)

# if p-value was less than alpha (e.g., 0.05), by contrast, then we would conclude that
# the observed data are statistically different from the expected values (lack of fit)


# Predict back to landscape
library(raster)
library(rgdal)

# Load the burned severity raster of our study area
fire <- raster("BurnSeverity.tif")
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(fire, col = mapPalette(100), axes = F, box = F, main = "% Burned Area")
res(fire) # 30 m x 30 m resolution (grid size)
# crs(fire)
CH <- as.data.frame(fire, xy=TRUE)
CH <- data.frame(x = CH$x, y = CH$y, burned = CH$BurnSeverity)

# Get predictions of abundance for each 30 m x 30 m cell of study area
newData <- data.frame(burned = CH$burned)
predCH <- predict(m1, type="state", newdata=newData)

# Define new data frame with coordinates and outcome to be plotted
PARAM <- data.frame(x = CH$x, y = CH$y, z = predCH$Predicted)
r1 <- rasterFromXYZ(PARAM)     # convert into raster object

# Plot predicted salamander abundance in study area (based on model m1)
par(mfrow = c(1,1), mar = c(1,2,2,5))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = F, main = "Salamander abundance (mean predicted values)")


# Plot predicted salamander abundance in study area (based on model m1)
par(mfrow = c(1,2))
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(r1, col = mapPalette(100), axes = F, box = F, main = "Salamander abundance")
plot(fire, col = mapPalette(100), axes = F, box = F, main = "% Burned Area")






