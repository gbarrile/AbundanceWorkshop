

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ABUNDANCE ESTIMATION WORKSHOP
# Sponsored by Wyoming EPSCoR
# Example 3: Closed population capture-recapture model

# Description: Steps through an example analysis of capture-mark-recapture data
# using the Closed Population Estimation Model ('Closed') as implemented 
# in the package `RMark`.

# Online resource for understanding the model that we will fit in this session:
# http://www.phidot.org/software/mark/docs/book/pdf/chap14.pdf

# Addressing the question: "How does boreal toad abundance vary across ponds?"

# Gabe Barrile - Colorado State University
# Last updated 12/03/2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- Outline -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# This script contains the following sections:
# 1) Install and load packages
# 2) Read in input data
# 3) Format input datasets for RMark
# 4) Process and explore datasets in RMark
# 5) Fit Closed Population Capture-Recapture Model
# 6) Prediction and plotting
# 7) Code to explore on your own (or as a group if there's time)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 1) Install and load packages -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Check if RMark, reshape, and ggplot2 are installed.
# If yes, load them.  If not, install, then load.

# RMark for fitting the Closed Population Capture-Recapture Model
if("RMark" %in% rownames(installed.packages()) == FALSE) {
  install.packages("RMark")
}
require(RMark)

# reshape to format data
if("reshape" %in% rownames(installed.packages()) == FALSE) {
  install.packages("reshape")
}
require(reshape)

# ggplot2 for plotting
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2")
}
require(ggplot2)

# if you needed to cite the RMark package
citation("RMark")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 2) Read in input data -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Specify path that contains the capture-mark-recapture data (BorealToad_CaptureRecapture.csv) 
# use setwd()

# read-in the boreal toad capture-mark-recapture data from our three ponds  
df <- read.csv("Session3-CaptureRecapture/data/BorealToad_CaptureRecapture.csv")

# We now have our data stored as 'df'
# check out our data
head(df)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION: How many unique individuals did we tag?
# Hint: can use length() and unique() functions
length(unique()) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 3) Format data for RMark -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# order dataframe by animal id
df <- with(df,df[order(Tag),])
head(df, 8)

# create capture histories for each individual toad
# don't worry too much about what this code is doing right now,
# but do return on your own time to make sure you understand
# how the data is being formatted
input.data <- df
input.data$detect <- rep(1,nrow(input.data))
m <- melt(input.data, id.var = c("Tag","Survey"), measure.var = "detect")
y = cast(m, Tag ~ Survey)
y[is.na(y)] = 0
k <- dim(y)[2]

#if needing to deal with repeats on 1 occasion
y[,2:k] <- (y[,2:k]>0)*1
head(y)
head(df, 8)

#function to create capture history strings
pasty<-function(x) 
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    y<-(x[i,]>0)*1
    out[i]<-paste(y,sep="",collapse="")
  }
  return(out)
}

# create capture history data frame 
capt.hist <- data.frame(ch = pasty(y[,2:k]))
head(capt.hist)

y$ch <- pasty(y[,2:k])
head(y)

# check if all capture histories are length = 4 (for each of our four visits to each pond)
nchar(y$ch)
table(nchar(y$ch)) # 'ch' must be character, makes sure it's not numeric

# add variables to 'y' dataframe
# Pond
y$Pond <- df$Pond[match(y$Tag, df$Tag)]

# SVL
svl <- aggregate(SVL ~ Tag, data=df, FUN=mean)

# match by tag number to add SVL to y dataframe
y$SVL <- svl$SVL[match(y$Tag, svl$Tag)]

# boto will be our dataframe that we input into RMark
boto <- data.frame(ch = y$ch, freq = 1, pond = y$Pond, 
                   tag = y$Tag,
                   svl = y$SVL)


# remove unneeded objects
rm(capt.hist,input.data,m,svl,y,k,pasty)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION: Check capture histories against our field data
# Look at toads with the following tag numbers: 100, 101, and 102
head(df, 10)
# Do the capture histories ('ch' column in boto) make sense for those three individuals?
head(boto, 3)
# Stop here and make sure you understand the capture histories (the 'ch' column in boto df)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make pond a factor variable
boto$pond <- as.factor(as.character(boto$pond))
table(boto$pond)

# are we missing any data?
table(is.na(boto)) # no missing data = good!
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 4) Process and explore data in RMark -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# process data in RMark 
d.proc=process.data(boto, model="Closed", groups = c("pond"))
# RMark includes lots of models that can be input into the model= argument above
# Here we use the "Closed" model because we are interested in 
# Closed Population Estimation for modeling Abundance
# Visit this link for a full list of MARK models supported in RMark:
# https://github.com/jlaake/RMark/blob/master/RMark/inst/MarkModels.pdf

# create design data
d.ddl <- make.design.data(d.proc)
# NOTE: the design dataframes are just as important as your raw data!
# design data in this example are fairly simple, 
# but will get more complex with different models

# Let's explore the design data
# see which parameters are estimated in the model
names(d.ddl)

# "p" = capture probability
# "c" = recapture probability
# "f0" = the number of individuals never captured (read more below)

# look at design data for each parameter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# capture probability
d.ddl$p
# look at columns 'time' and 'pond'
# QUESTION: why do we have 12 rows in this table? 
# Discuss with the group, focus on those two columns mentioned above.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# recapture probability
d.ddl$c
# all individuals captured during survey 1 were captured for the first time,
# thus recapture probability starts at time = 2 or survey 2

# the number of individuals never captured
d.ddl$f0
# the notation 'f0' originates from the frequency (count) of animals observed 0 times). 
# The f0 parametrization is useful computationally because f0 is bounded on the 
# interval [0, positive infinity], thus forcing the logical constraint that 
# N > the number of individuals marked (M). 
# In fact, MARK uses the f0 parametrization for ease of computation by using 
# the log link function to constrain f0 >= 0, and presents the results in terms of N 
# as a derived parameter.
# In other words, the likelihood is rewritten in MARK in terms of the 
# number of individuals never caught, f0, such that f0 = N - M 
# (M = the number of individuals marked)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 5) Fit Closed Population Capture-Recapture Model -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# specify models for each parameter

# Model 1
# p (capture probability)
pc.=list(formula= ~ 1, share = TRUE)
# share = TRUE indicates p = c or capture probability = recapture probability
# for this session, we are going to assume that capture probablity is
# equal to recapture probability

# ~ 1 or the intercept model is often referred to as the 'constant' model, 
# meaning that capture probability is constant over time and space 
# (e.g., doesn't vary across surveys or at different ponds)


# f0 (number of individuals never captured)
f0.=list(formula= ~ 1)
# this model indicates that the same number of individuals were
# never captured at all ponds
# (e.g., we failed to capture and mark 10 individuals at
# Pond 1, 10 individuals at Pond 2, and 10 at pond 3)
# it does NOT mean that the abundance is the same at all ponds,
# just that the number we failed to capture is the same


# fit Model 1
# you will need (want) an output folder for the mark files, just so your
# directory does not get cluttered
# Create a new folder called 'models' in your working directory
# set working directory to that folder
setwd("H:/Abundance_Workshop/models")

# mark() is the model-fitting function 
m1 <- mark(d.proc, # process data
           d.ddl,  # design data
           model.parameters = list(p  = pc.,  # model for capture probability
                                   f0 = f0.)) # model for f0


# look at model output
# m1 # this brings up the output you would get in Program MARK

# beta coefficients
m1$results$beta
# these are just the intercepts in this case

# real estimates (on the scale we tend to think on)
m1$results$real
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION: what is our capture probability?
# How many individuals were never captured at each site?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# derived parameters
m1$results$derived
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION: what is the estimated abundance at each pond?
# Pond 1 = first row, Pond 2 = second row, Pond 3 = third row
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Let's look at how abundance is estimated
# f0
f0 <- m1$results$real[2,1]
f0
# how many toads did we mark at each pond?
table(boto$pond)
# Add f0 to the number marked at each pond
c(f0 + 44, f0 + 18, f0 + 75)
# compare to the estimated abundances from the derived parameters in the model
m1$results$derived$`N Population Size`[,1]
# the numbers are identical 
# make sure you understand why those numbers are identical...Discuss as a group

# In other words, the model knows how many toads we tagged at each pond:
table(boto$pond)
# then it adds the estimated number of individuals that we failed to tag
f0
# Note that we can estimate a different f0 for each pond (below)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Model 2

# Obtain different capture probability for each survey at each pond

# We conducted 4 surveys at each of three ponds, so 12 surveys total
# Therefore we should end up with 12 different capture probabilities,
# one for each survey at each site

# p
d.ddl$p
p.timepond =list(formula= ~ time * pond, share = TRUE)
# 'time' and 'pond' in the above formula must be spelled 
# exactly as shown in design data

# Different f0 for each pond
# f0
d.ddl$f0
f0.pond =list(formula= ~ pond)

# fit model
# again specify your output folder
setwd("H:/Abundance_Workshop/models")

# mark() is the model-fitting function 
m2 <- mark(d.proc, # process data
           d.ddl,  # design data
           model.parameters = list(p  = p.timepond, # model for capture probability
                                   f0 = f0.pond))   # model for f0


# look at model output

# real estimates 
m2$results$real
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STOP! Make sure you understand why there are:
# 12 estimates for capture probability and
# 3 estimates for f0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# derived parameters
m2$results$derived
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QUESTION: which pond had the highest abundance? lowest?
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# compare models with AICc
c(m1$results$AICc,m2$results$AICc)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 6) Plot model predictions -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# let's plot the estimated abundance at each pond,
# using estimates from Model 2 (m2)
abund <- m2$results$derived$`N Population Size`
abund$pond <- c("Pond 1", "Pond 2", "Pond 3")
abund$pond <- as.factor(abund$pond)
abund

# plot the values
ggplot(abund, aes(x=pond, y=estimate, color=pond)) +
  geom_errorbar(aes(ymin=lcl, ymax=ucl), 
                width=0.2, size=2, position=position_dodge(0.3)) +
  geom_point(position=position_dodge(0.3), stat="identity", size = 7) +
  xlab(NULL) +
  ylab("Estimated Abundance \n") +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  theme_bw()+
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
    axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 20, vjust = 1, color = "black"), # spacing from y
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"))


# let's plot capture probability during each survey at each pond,
# again using estimates from Model 2 (m2)
cap <- m2$results$real[1:12,c(1,3,4)]
cap$pond <- c("Pond 1", "Pond 1", "Pond 1", "Pond 1",
              "Pond 2", "Pond 2", "Pond 2", "Pond 2",
              "Pond 3", "Pond 3", "Pond 3", "Pond 3")
cap$pond <- as.factor(cap$pond)

cap$survey <- c(1,2,3,4,
                1,2,3,4,
                1,2,3,4)

cap$survey <- as.factor(as.character(cap$survey))

cap

# plot the values
ggplot(cap, aes(x=survey, y=estimate, group=pond, color=pond)) +
  geom_point(size=5) + geom_line(size=1.5) +
  xlab("Survey") +
  ylab("Capture probability \n") +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = c(0,0)) +
  theme_bw() +
  theme(#text = element_text(size = 18, family = "Times"), # controls all fonts in plot
    panel.background = element_rect(colour = "black", size=1, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2,"cm"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(size = 17))+
  theme(legend.position="top")+
  theme(strip.text = element_text(size = 7))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# ---- 7) Code to explore on your own (or as a group if there's time) -----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Now let's incorporate individual covariates (e.g., SVL) into our models

# process data in RMark 
d.proc=process.data(boto, model="Huggins") 
# let's not worry about different ponds in this example, so no groups = c("pond")

# Huggins models condition abundance out of the likelihood,
# thus permitting the modelling of capture probability as a function of individual covariates

# create design data
d.ddl <- make.design.data(d.proc)
# NOTE: the design dataframes are just as important as your raw data!
# design data in this example are fairly simple, 
# but can get more complex with different models

# Let's explore the design data
# see which parameters are estimated in the model
names(d.ddl)

# capture and recapture probabilities can be the same or can be different 
# let's say you catch a rabbit in a trap for the first time 
# (that's your capture probability)
# Then let's say that rabbit subsequently avoids your traps
# In that case, it's very likely that your recapture probability of that rabbit
# is different than the probability of capturing that rabbit the first time
# This is just one example whereby capture and recapture probability can differ

# look at design data for each parameter

# capture probability
d.ddl$p

# recapture probability
d.ddl$c
# all individuals captured during survey 1 were captured for the first time,
# thus recapture probability starts at time = 2 or survey 2

# the number of individuals never captured
d.ddl$f0
# In this model, the likelihood is conditioned on the number of animals detected 
# and f0 therefore drops out of the likelihood. 
# These models contain only p and c, with abundance N estimated as a derived parameter.
# As noted earlier, the primary advantage of the Huggins data type is that 
# individual covariates can be used to model p and c.


# specify models for each parameter

# Model 3

# p (capture probability)
pc.=list(formula= ~ 1, share = TRUE)
# share = TRUE indicates p = c or capture probability = recapture probability
# ~ 1 or the intercept model is often referred to as the 'constant' model, 
# meaning that capture probability is constant over time and space 
# (e.g., doesn't vary across surveys or at different ponds)


# fit Model 3
# specify your output folder
setwd("H:/Abundance_Workshop/models")

# mark() is the model-fitting function
m3 <- mark(d.proc, # process data
           d.ddl,  # design data
           model.parameters = list(p = pc.)) # model for capture probability


# look at model output

# beta coefficients
m3$results$beta
# these are just the intercepts in this case

# real estimates (on the scale we tend to think on)
m3$results$real

# derived parameters
m3$results$derived

# how is abundance estimated?

# what is the probability of an individual not being captured at all during our study?
# (1 - p)(1 - p)(1 - p)(1 - p)
# what is p, our capture probability?
p <- m3$results$real$estimate[1]
p
# the probability of an individual not being captured at all during our study?
p.no <- (1 - p)*(1 - p)*(1 - p)*(1 - p)
# so what is the probability of being captured at least once?
p.once <- 1 - p.no
# then we can say that the number of individuals that we marked was equal to
# the probability of being captured at least once times the total population size (N)

# so we need to obtain the total number that we marked 
marked <- nrow(boto)

# here's the equation:      p.once * N = marked
# rearrage to solve for N:  N = marked / p.once

# So N should be:
marked / p.once
# compare to derived estimate from model 
m3$results$derived





# Now let's use svl as individual covariate on capture probability

# p (capture probability)
pc.svl =list(formula= ~ svl, share = TRUE)

# fit Model 4
# specify your output folder
setwd("H:/Abundance_Workshop/models")

# mark() is the model-fitting function
m4 <- mark(d.proc, # process data
           d.ddl,  # design data
           model.parameters = list(p = pc.svl)) # model for capture probability

# look at model output

# beta coefficients
m4$results$beta

# real estimates (on the scale we tend to think on)
m4$results$real

# derived parameters
m4$results$derived

# plot relationship between capture probability and SVL
range(boto$svl)
# make predictions based on model (m4)
pred.svl <- covariate.predictions(m4, data=data.frame(svl=56:88),indices=c(1))$estimates

min(pred.svl$lcl)
max(pred.svl$ucl)
op <- par(mar = c(5,5,4,2) + 0.1) # default is 5,4,4,2
plot(x = pred.svl$covdata, y = pred.svl$estimate, pch=16, 
     ylab = "Capture Probability",
     xlab = "SVL (mm)", cex.lab=1.5, cex.axis=1.2, 
     col="darkgray", ylim=c(0.1,0.5))
box(lwd = 4, col = 'black')
lines(pred.svl$covdata, pred.svl$estimate, lwd=8, col="blue")
lines(pred.svl$covdata, pred.svl$lcl, lwd=4, lty=2, col="black")
lines(pred.svl$covdata, pred.svl$ucl, lwd=4, lty=2, col="black")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# END
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#






