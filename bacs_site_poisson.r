#############################################################################
#BACS Poisson Analysis
#Use counts to incorporate information in the two time periods to estimate abundance and
#differences in detection pobability
#2/19/14
#Scott Anderson, NC Wildlife Resources Commission

#############################################################################
#begin analysis of sites

#load R Libraries
library(beanplot)
library("R2WinBUGS") #to run winBUGS. Requires 'coda' package
library(reshape)
library(coda)
bugs.dir = "C:/WinBUGS14/"


#set working directory
setwd("C:/data/@projects/bacs_surveys/analysis/heat_map/r/")

#############################################################################
#import data from csv and modify (Access export)
habq <- read.table("BACSHabData.csv", header=TRUE, sep=",")

#create a subset with surveyed locations only
habq_survey <- habq[ which(habq$NOSURV==0),]
habq_survey$OBS <- factor(habq_survey$OBS)

#re-factor habitat variables to be 2 factors instead of 3
habq_2factor <- habq_survey
habq_2factor$BA <- factor(ifelse(habq_2factor$BA<=1,0,1))
habq_2factor$FIRE <- factor(ifelse(habq_2factor$FIRE<1,0,1))
habq_2factor$WIRE <- factor(ifelse(habq_2factor$WIRE<=1,0,1))
habq_2factor$BARE <- factor(ifelse(habq_2factor$BARE<1,0,1))
habq_2factor$FORB <- factor(ifelse(habq_2factor$FORB<1,0,1))
habq_2factor$WOOD <- factor(ifelse(habq_2factor$WOOD<1,0,1))
habq_2factor$MIDHARD <- factor(ifelse(habq_2factor$MIDHARD<=1,0,1))

#############################################################################

#create subset for analysis
bacs <- na.omit(subset(habq_2factor, ,select=c(PASS, ACTIVE, JULIAN, HABQI, BTIME, OBS)))	#only those sites surveyed, and with habitat quality index

y<-as.matrix(bacs[,1:2])		#detection data for passive and active listening periods
y[y>0]<-1						#convert to binary (detection/non-detection)
active<- y						#indicator for passive vs. active 
active[,1]<-0
active[,2]<-1

dates<-bacs$JULIAN				#julian day of the year

#standardize dates
mean.date<-mean(dates,na.rm=TRUE)
sd.date<-sd(dates[!is.na(dates)])
DATES<-(dates-mean.date)/sd.date	#standardize
DATES[is.na(DATES)]<-0				#impute zeroes(means)

#standardize times
times <- bacs$BTIME
mean.times<-mean(times)
sd.times<-sd(times)
TIMES<-(times-mean.times)/sd.times
TIMES[is.na(TIMES)]<-0

#factorize observers
obs<-factor(bacs$OBS)

