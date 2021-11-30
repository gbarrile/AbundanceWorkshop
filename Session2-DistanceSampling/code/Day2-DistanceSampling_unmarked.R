#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ABUNDANCE ESTIMATION WORKSHOP
# Sponsored by Wyoming EPSCoR
# Example 2:  Distance-sampling analysis

# Description:  Steps through an example analysis of distance-sampling data
# using the hierarchical distance-sampling model of Royle et al. 2004 as
# implemented in the package `unmarked`.

# Jason Carlisle - WEST, Inc.
# Last updated 11/23/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- Outline -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# This script contains the following sections:
# 1) Install and load packages
# 2) Read in input data
# 3) Format input datasets for unmarked
# 4) Fit a hierarchical distance-sampling model
# 5) Model selection with covariates
# 6) Prediction and plotting
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 1) Install and load packages -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Check if unmarked, ggplot2, and cowplot are installed.
# If yes, load them.  If not, install, then load.

# unmarked for fitting the hierarchical distance sampling model
if("unmarked" %in% rownames(installed.packages()) == FALSE) {
  install.packages("unmarked")
}
require(unmarked)

# ggplot2 for plots
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2")
}
require(ggplot2)

# cowplot for a custom theme for ggplot plots
if("cowplot" %in% rownames(installed.packages()) == FALSE) {
  install.packages("cowplot")
}
require(cowplot)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 2) Read in input data -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# The basic data format you should start with is a sites data.frame and a
# detections data.frame.  The siteData has one row per site surveyed, and the
# detectionData has one row per pronghorn detected.
# We'll use a real-world dataset from aerial pronghorn surveys conducted by
# the Wyoming Game and Fish Department near Laramie, WY.

# Specify the path that contains the two example datsets provided:
# pronghornSiteData.RDS and pronghornDetectionData.RDS
setwd("W:/My Drive/Workshops/Abundance/PrepPronghornData")


# The first required dataset is a sites data.frame
# Each row is a site (here, line transect) that was surveyed (here, each once)
# Includes columns for site-level covariates
siteData <- readRDS("pronghornSiteData.RDS")
nrow(siteData)  # 186 transects
head(siteData)
table(siteData$herd)
table(siteData$herd)

# The second required dataset is a detection data.frame
# Each row is a detected individual pronghorn
detectionData <- readRDS("pronghornDetectionData.RDS")
nrow(detectionData)  # 4110 pronghorn individuals
head(detectionData)


# Without getting too detailed, the survey protocol includes assigning each
# detected pronghorn to a distance band perpendicular to the plane
# (e.g., 100 and 150 m).  Then an approximate exact distance is estimated that
# accounts for the flight altitude at the time of the detection (e.g., 126 m).
# The result is a quasi-continuous measure of the perpendicular
# distance between the inner edge of the survey strip and the pronghorn.

# Nominal distance bands
distanceBands <- c(0, 20, 45, 80, 145, 200)
distanceBandMidPoints <- distanceBands[-length(distanceBands)] + diff(distanceBands)/2


# Histogram of detection distances including distance bands
ggplot(detectionData, aes(distM)) + theme_cowplot() +
  geom_histogram(binwidth=5, center=2.5, fill="dodgerblue", color="white") +
  scale_x_continuous(limits=c(0, 300),
                     expand=expansion(mult=c(0, 0.05))) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.10))) +
  geom_vline(xintercept=distanceBands) +
  xlab("Distance (m)") + ylab("Count") +
  annotate("text", x=c(distanceBandMidPoints, 250), y=550, size=5,
           label=c("A", "B", "C", "D", "E", "Outside\nSurvey Strip"))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 3) Format input for unmarked -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# unmarked uses a multinomial-Poisson mixture model to estimate abundance from
# distance-sampling data (Royle et al. 2001). They have a great vignette on CRAN
# about the distance-sampling component of the package.
# In `unmarked` detection distances are binned prior to analysis. 
# Generally more bins is better (Kery and Royle 2016).  We'll use 4 bins here.

# It's common in distance-sampling to exclude the right tail of the distances
# to improve model stability.  Since the formal survey strip was only 200 m
# wide, we'll truncate any distances recorded beyond 200 m
# For unmarked, we just remove these data to right truncate
trunc <- 200


# Histogram of distances collapsed into 50-m wide bins and showing the
# the distance beyond which we plan to truncate or censor the data (i.e., 200 m)
ggplot(detectionData, aes(distM)) + theme_cowplot() +
  geom_histogram(binwidth=50, center=25, fill="dodgerblue", color="white") +
  scale_x_continuous(limits=c(0, 300),
                     expand=expansion(mult=c(0, 0.05))) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.10))) +
  geom_vline(xintercept=trunc, lty="dashed") +
  xlab("Distance (m)") + ylab("Count") +
  annotate("text", x=trunc+5, y=800, size=5, hjust=0,
           label = "Right truncation distance")


# Manually truncate
detectionData <- detectionData[detectionData$dist <= trunc, ]
nrow(detectionData)  # number of detected individuals (3916)

# Histogram of distances with 50 m bins after truncation
ggplot(detectionData, aes(distM)) + theme_cowplot() +
  geom_histogram(binwidth=50, center=25, fill="dodgerblue", color="white") +
  scale_x_continuous(expand=expansion(mult=c(0, 0.05))) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.10))) +
  xlab("Distance (m)") + ylab("Count")


# Add rectangle to illustrate the key assumption of distance sampling
ggplot(detectionData, aes(distM)) + theme_cowplot() +
  geom_rect(xmin=0, xmax=200, ymin=0, ymax=sum(detectionData$distM<50),
            fill="pink") +
  geom_histogram(binwidth=50, center=25, fill="dodgerblue", color="white") +
  scale_x_continuous(expand=expansion(mult=c(0, 0.05))) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.10))) +
  xlab("Distance (m)") + ylab("Count") +
  geom_rect(xmin=0, xmax=200, ymin=0, ymax=sum(detectionData$distM<50),
            fill=NA, color="black", size=2) +
  annotate("text", x=c(60, 130), y=c(450, 1200), size=10,
           label=c("Detected", "Not Detected"), color=c("blue", "red"))


# Format (bin) into a matrix with dimensions of number of sites (M or R) 
# by number of distance intervals (J).  Each cell of the matrix is the observed
# count in that distance interval for that transect

# Remember, we want 4 bins, and right truncation at 200 m
nBins <- 4
(distanceBreaks <- seq(0, trunc, length.out=nBins+1))

countMatrix <- formatDistData(distData=detectionData,
                              distCol="distM",
                              transectNameCol="siteID",
                              dist.breaks=distanceBreaks)

# Look at the first few rows of the matrix
head(countMatrix, 10)

# Question:  How many pronghorn were detected in the 150-200 m distance bin
# at transect Centennial-12?


# Organize the count matrix along with covariates and metadata
# Note that unmarked does not match siteIDs between detection and site data,
# so make sure they're both ordered the same.

# Site-level covariate data
head(siteData)

# Order both by siteID before combining, just to be sure
countMatrix <- countMatrix[order(row.names(countMatrix)), ]
siteData <- siteData[order(siteData$siteID), ]

# The protocol is to only survey one side of the transect, but unmarked assumes
# both sides are, so unmarked would assume 2x the area was surveyed.
# To account for this, simply divide the survey lengths by 2.
uFrame <- unmarkedFrameDS(y=countMatrix,
                          siteCovs=siteData,
                          dist.breaks=distanceBreaks,
                          tlength=siteData$transectLengthM/2,
                          survey="line",
                          unitsIn="m")

# Double-check that the detection and site data are correctly paired
identical(row.names(uFrame@y), as.character(uFrame@siteCovs$siteID))

# Explore uFrame, the input to unmarked distance sampling routine
# Good point to stop and ensure you're feeding the model what you think you are
head(uFrame)
summary(uFrame)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 4) Fit a hierarchical distance-sampling model -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Fit half-normal detection function with no covariates for detection or density
# Order for formula is detection, then abundance
fit <- distsamp(~1 ~1, uFrame, keyfun="halfnorm",
                  output="density", unitsOut="kmsq")

# Look at output
fit
coef(fit)

# Back-transformed output
# On count or density scale
backTransform(fit, type="state")
backTransform(fit, type="det")

# Same as back-transforming manually with appropriate link function
exp(coef(fit)[1])


# Plot detection function
hist(fit, col="grey", xlab="Distance (m)", lwd=2)


# Probability of detection
(detProb <- integrate(gxhn, 0, trunc,
                  sigma=backTransform(fit, type="det")@estimate)$value/trunc)



# Estimate pronghorn density per sq. km (with 95% CI)
backTransform(fit, type="state")
predict(fit, type="state", level=0.95)[1, ]  # Row for each site if no "newdata" argument



# Calculate density by hand to check you understand
# count / detection prob / area in km^2
# Number of pronghorn in survey strips
(obsCount <- nrow(detectionData))  # 3,916 detected
detProb  # 76% detection probability
obsCount/detProb  # adjusted estimate of 5,123 pronghorn in survey strips

# Area of the surveyed strips
# Only one side of the transect surveyed.  If both sides, multiply trunc by 2
(survAreaKm2 <- (sum(siteData$transectLengthM)*trunc)/1e6)

obsCount / detProb / survAreaKm2  # 6.2 pronghorn per sq. km

# Matches output from model
backTransform(fit, type="state")  # 6.2 pronghorn per sq. km
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 5) Model selection with covariates -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# The objective is to determine whether the 3 herd units have equivalent
# pronghorn density. But it's also possible that differences in topography and
# vegetation produce different detection probabilities in each herd unit.
# Therefore, we explore 4 models, some that assume detectability and density
# are the same among herd units, and some that can account for herd-unit
# differences in detectability, density, or both.
# It's also common to explore different shapes of the detection curve (called
# key functions), but here we'll assume the half-normal shape for all models.

# Possible covariates
head(siteData)

# Fit 4 models
# Detectability and density constant across herd units
hn.null.null <- distsamp(~1 ~1, data=uFrame, keyfun="halfnorm",
                    output="density", unitsOut="kmsq")

# Abundance varies by herd unit
hn.null.herd <- distsamp(~1 ~herd, data=uFrame, keyfun="halfnorm",
                         output="density", unitsOut="kmsq")

# Detectability varies by herd unit
hn.herd.null <- distsamp(~herd ~1, data=uFrame, keyfun="halfnorm",
                        output="density", unitsOut="kmsq")

# Detectability and abundance vary by herd unit
hn.herd.herd <- distsamp(~herd ~herd, data=uFrame, keyfun="halfnorm",
                             output="density", unitsOut="kmsq")



# Make AIC table comparing models
# The warning here is fine, we already gave the model objects meaningful names
(aic <- modSel(fitList(hn.null.null, 
                       hn.herd.null, 
                       hn.null.herd, 
                       hn.herd.herd)))


# Store top model
best <- hn.herd.herd
best
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 6) Prediction and plotting -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Given this best model, plot estimated detection function and abundance

# Detection probability (p)
# Calculate p (and CI) from sigma
newDataDet <- data.frame(herd=levels(siteData$herd))
predDet <- predict(best, type="det", newdata=newDataDet, appendData=TRUE)

for(i in 1:nrow(predDet)) {
  predDet$p[i] <- integrate(gxhn, 0, trunc,
                              sigma=predDet$Predicted[i])$value/trunc
  predDet$p.lower[i] <- integrate(gxhn, 0, trunc,
                                    sigma=predDet$lower[i])$value/trunc
  predDet$p.upper[i] <- integrate(gxhn, 0, trunc,
                                    sigma=predDet$upper[i])$value/trunc
}

# Plot detection probabilities
ggplot(predDet, aes(x=herd, y=p, color=herd)) + theme_cowplot() +
  geom_point(size=3) +
  scale_color_viridis_d(name="Herd unit", alpha=0.8) +
  geom_errorbar(aes(ymin=p.lower, ymax=p.upper), width=0.05) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(mult=c(0, 0.01))) +
  ylab("Detection probability (p)") +
  xlab("Herd Unit")



# Plot detection functions
detFunData <- do.call(rbind, lapply(1:nrow(predDet), function(i) {
  data.frame(dist=0:trunc,
             p=gxhn(x=0:trunc, sigma=predDet$Predicted[i]),
             herd=predDet$herd[i])
}))

ggplot(detFunData, aes(x=dist, y=p, color=herd)) + theme_cowplot() +
  geom_line(size=1.5) +
  scale_color_viridis_d(name="Herd unit", alpha=0.8) +
  scale_y_continuous(limits=c(0, 1), expand=expansion(mult=c(0, 0.01))) +
  scale_x_continuous(limits=c(0, trunc), expand=expansion(mult=c(0, 0))) +
  ylab("Detection probability (p)") +
  xlab("Distance (m)")



# Plot pronghorn density estimates
newDataAbund <- data.frame(herd=levels(siteData$herd))
predAbund <- predict(best, type="state", newdata=newDataAbund, appendData=TRUE)

ggplot(predAbund, aes(x=herd, y=Predicted, fill=herd)) + theme_cowplot() +
  geom_bar(stat="identity") +
  scale_fill_viridis_d(name="Herd unit", alpha=0.8) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.05) +
  scale_y_continuous(limits=c(0, 8), expand=expansion(mult=c(0, 0))) +
  # scale_x_continuous(limits=c(0, trunc), expand=expansion(mult=c(0, 0))) +
  ylab("Density (pronghorn/km^2)") +
  xlab("Herd Unit")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# END
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#