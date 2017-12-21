
## Required packages
library(R2jags)
library(mcmcplots)
library(plyr)
library(Cairo)
library(tidyverse)

# working directory in the AGHY_Final Folder
setwd("C:/Users/Marion Donald/Dropbox/Rice/Projects/AGHY/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_Final")
## load workspace which has all of the data files already read in
#load("AGHY_data_clean.Dec18.RData")

## Bayesian model for endo frequency change, survival and flowering (working on adding in recruitment) 

sink("AGHY_seed_mass.txt")
cat("
    model{




#####################################################################
########## NEW: Seed Production December 20, 2017 #######################
#####################################################################

## First start with change from seed mass to seed counts

## priors for regression paramaters endo specific and indexed over water treatments 
for (i in 1:N.trt){
    intercept.ep[i] ~ dnorm(0, 0.001)
    slope.ep[i] ~ dnorm(0, 0.001)  
    
  intercept.em[i] ~ dnorm(0, 0.001)
   slope.em[i] ~ dnorm(0, 0.001) 
  
 
  
}

## sigma and tau for seed counts
  sigma.count ~ dunif(0, 1000)
tau.sigma.count ~ dgamma(0.001, 0.001)



### Likelihood model for number of seeds and seed mass
for (i in 1:N.Ep.seed.count){
mean.seed.count.Ep[i] <- intercept.ep[water.seed.count.Ep[i]] + slope.ep[water.seed.count.Ep[i]]*seed.mass.Ep[i]
Ep.seed.count[i] ~ dnorm(mean.seed.count.Ep[i], tau.sigma.count)

}

for (i in 1:N.Em.seed.count){
mean.seed.count.Em[i]<- intercept.em[water.seed.count.Em[i]] + (slope.em[water.seed.count.Em[i]])*seed.mass.Em[i]
Em.seed.count[i] ~ dnorm(mean.seed.count.Em[i], tau.sigma.count)

}


## Seed production for the original plants in 2013 (these seeds were produced by original plants in 2013 and will contribute to the recruits in 2014)

### Priors for seed production - specific for endo, indexed for water treatment
  for(i in 1:N.trt){
  Ep_beta0_seed.13[i] ~ dnorm(0, 0.001)
  Em_beta0_seed.13[i] ~ dnorm(0, 0.001)
  }


## Tau.sigma for prob seed production plot (this tau is different than the tau that will go in the gaussian model)
  sigma.seed ~ dunif(0, 1000)
  tau.sigma.seed <- 1/(sigma.seed * sigma.seed)

## Tau.sigma for the model
  sigma.s ~ dunif(0,1000)
  tau.sigma.s <- 1/(sigma.s * sigma.s)

## Random effect of plot for seed production
  for(i in 1:N.plots){
  ran.seed.13[i] ~ dnorm(0, tau.sigma.seed)


## Apply random effect of plot on seed production 
  Ep_seed.13[i]<- Ep_beta0_seed.13[water[i]] + ran.seed.13[i]
  Em_seed.13[i]<- Em_beta0_seed.13[water[i]] + ran.seed.13[i]
  }



## Likelihood estimate for seed production from E+ and E- ORIGINAL plants
  for(i in 1:N.Ep.seed.13){
## eventually need to multiply the mean seed number per plot by the slope from the linear regression to turn the seed mass into seed count
###### ^ this is done in the derived quantity below for mean number of seeds contributed by each E+ original plant by water trt
  Ep.seed.mass.13[i] ~ dnorm(Ep_seed.13[seed.Ep.plot.13[i]], tau.sigma.s)

  }

  for(i in 1:N.Em.seed.13){
## eventually need to multiply the mean seed number per plot by the slope from the linear regression to turn the seed mass into seed count
###### ^ this is done in the derived quantity below for mean number of seeds contributed by each E+ original plant by water trt
    Em.seed.mass.13[i] ~ dnorm(Em_seed.13[seed.Em.plot.13[i]], tau.sigma.s)

    }


## derived quantity for endo-specific and water treatment seed mass prob

     Ep.orig.seed.control.13<-Ep_beta0_seed.13[2]*slope.ep[2]
    Ep.orig.seed.add.13<-Ep_beta0_seed.13[1]* slope.ep[1]
    Em.orig.seed.control.13<-Em_beta0_seed.13[2]*slope.em[2]
    Em.orig.seed.add.13<-Em_beta0_seed.13[1]*slope.ep[1]



#####################################################################
#####################################################################

}##end model
",fill=T)
sink()



## bundle data
#### Bayesian model for plot frequency change 13/14/15
## define data for Bayes
N.yrs<-2
N.trt<-2

## here are the plot data
water<-as.integer(AGHY.plots$water[!is.na(AGHY.plots$plot)])
N.plots<-length(water)
init_freq<-AGHY.plots$target_init_freq[!is.na(AGHY.plots$plot)]

##subplots E+ scores
y.14<-AGHY.merge$E_plus_liberal[AGHY.merge$year_t==2014]
N.samples.14<-AGHY.merge$total[AGHY.merge$year_t==2014]
N.obs.14<-length(y.14)
y.14.plot<-AGHY.merge$newplot[AGHY.merge$year_t==2014]

## Create unique ID column for year, plot, subplot, ID
AGHY.demo.final$unique.id<- with(AGHY.demo.final, paste(plot, subplot, ID, year_t, sep="_"))


y.15<-AGHY.merge$E_plus_liberal[AGHY.merge$year_t==2015]
N.samples.15<-AGHY.merge$total[AGHY.merge$year_t==2015]
N.obs.15<-length(y.15)
y.15.plot<-AGHY.merge$newplot[AGHY.merge$year_t==2015]

y.16<-AGHY.merge$E_plus_liberal[AGHY.merge$year_t==2016]
N.samples.16<-AGHY.merge$total[AGHY.merge$year_t==2016]
N.obs.16<-length(y.16)
y.16.plot<-AGHY.merge$newplot[AGHY.merge$year_t==2016]

## define terms for survival analysis
#N.samples.surv.15<-AGHY.demo.final$spring_survival_t1[AGHY.demo.final$year_t1==2015 & AGHY.demo.final$spring_survival_t1==1]
############################
######## OLD SURVIVAL ANALYSIS
### subset e pos and neg survivors
#e.pos.surv<-subset(AGHY.demo.final,year_t1==2015 & spring_survival_t1==1)
## create a water vector to sample over for the survival analysis
#e.pos.surv$water.num<-as.integer(e.pos.surv$water[!is.na(e.pos.surv$plot)])
#s.e.pos.15.data<-subset(e.pos.surv, endo.stat_2015 == 1)
#s.e.pos.15<-as.numeric(s.e.pos.15.data$endo.stat_2015 ==1)
#s.e.pos.15.water<-s.e.pos.15.data$water.num
#N.obs.pos.surv.15<-length(s.e.pos.15)
#s.pos.surv.plot.15<-s.e.pos.15.data$newplot[s.e.pos.15.data$year_t1==2015]

#s.e.neg.15.data<-subset(e.pos.surv, endo.stat_2015 == 0)
#s.e.neg.15<-as.numeric(s.e.neg.15.data$endo.stat_2015 ==0)
#s.e.neg.15.water<-s.e.neg.15.data$water.num
#s.neg.surv.plot.15<-s.e.neg.15.data$newplot[s.e.neg.15.data$year_t1==2015]
#N.obs.neg.surv.15<-length(s.e.neg.15)

##########################################################################
### New Survival Analysis - July 6, 2017

## subset the endo stat known individuals in 2015 
s.known.15<-subset(AGHY.demo.final, year_t1 == 2015 & !is.na(spring_survival_t1)) ## select just the plants that have survival data in 2015
s.Ep.known.15.data<-subset(s.known.15, endo.stat == 1) ## select out the E+ plants in 2015
s.Em.known.15.data<-subset(s.known.15, endo.stat == 0) ## select out the E- plants in 2015

s.Ep.known.15<-s.Ep.known.15.data$spring_survival_t1
s.Em.known.15<-s.Em.known.15.data$spring_survival_t1

N.Ep.known.15<- length(s.Ep.known.15)
N.Em.known.15<- length(s.Em.known.15)

s.Ep.known.15.plot<-s.Ep.known.15.data$newplot
s.Em.known.15.plot<-s.Em.known.15.data$newplot


## subset the endo stat known individuals in 2016
s.known.16<-subset(AGHY.demo.final, year_t1 == 2016 & !is.na(spring_survival_t1))
s.Ep.known.16.data<-subset(s.known.16, endo.stat == 1) ## select out the E+ plants in 2016
s.Em.known.16.data<-subset(s.known.16, endo.stat == 0) ## select out the E- plants in 2016

s.Ep.known.16<-s.Ep.known.16.data$spring_survival_t1 ## most of the plants that have E+ endo status are dead
s.Em.known.16<-s.Em.known.16.data$spring_survival_t1 ## most of the plants that have E- endo status are ALSO dead

N.Ep.known.16<- length(s.Ep.known.16)
N.Em.known.16<- length(s.Em.known.16)

s.Ep.known.16.plot<-s.Ep.known.16.data$newplot
s.Em.known.16.plot<-s.Em.known.16.data$newplot

## data for 2015 at the subplot level for plants with unknown endo status
## Make another column for unique ID for plot and subplot and year so that we can get the total number of demo plants per subplot
AGHY.demo.final.14<- AGHY.demo.final[which(AGHY.demo.final$year_t=="2014"),]
AGHY.demo.final.14$sbplt.unique.id <- with(AGHY.demo.final.14, paste(plot, subplot, year_t, sep="_"))
AGHY.demo.final.14.d<-ddply(AGHY.demo.final.14, "sbplt.unique.id", plyr::summarize, 
                            total_demo.14 = length(ID))
## TOM 6/21, here is where we could tally the plants that were scored, then sort them out

#### NEED TO DO THIS AT THE SUBPLOT LEVEL. Subset data for endo.stat unknown and survival known in 2015
## THIS IS JUST FOR PLANTS FOR WHICH WE DO NOT KNOW THE ENDO STAT
surv.demo.15.data<-subset(AGHY.demo.final.14, is.na(AGHY.demo.final.14$endo.stat) & !is.na(AGHY.demo.final.14$spring_survival_t1)) ## select out the plants that have known survival in 2015 but for which we don't know endo stat

AGHY.demo.final.15.d<-ddply(surv.demo.15.data, c("newplot","sbplt.unique.id"), plyr::summarize, 
                            alive.15 = sum(spring_survival_t1),
                            total_demo.14 = length(ID)) ## total number of demo plants in this is not necessarily 4 since it's only for the plants for which we don't know the endo stat

alive.demo.15.plot <- AGHY.demo.final.15.d$newplot ## new plot numbers that match with plant ID, to tie in plot probability 

N.demo.alive.15 <- length(alive.demo.15.plot) ## length of alive demography plants in 15 to loop over 

demo.alive.15<-AGHY.demo.final.15.d$alive.15 ## number of demo plants per subplot that SURVIVED (for which we don't know endo status)

N.samples.14.demo<-AGHY.demo.final.15.d$total_demo.14 ## total number of demo plants per subplot for which we DO NOT know endo status

## data for 2016 plants with unknown endo status
AGHY.demo.final.16<- AGHY.demo.final[which(AGHY.demo.final$year_t=="2015"),]
AGHY.demo.final.16$sbplt.unique.id <- with(AGHY.demo.final.16, paste(plot, subplot, year_t, sep="_"))


surv.demo.16.data<-subset(AGHY.demo.final.16, is.na(AGHY.demo.final.16$endo.stat) & !is.na(AGHY.demo.final.16$spring_survival_t1)) ## select out the plants that have known survival in 2015 but for which we don't know endo stat

AGHY.demo.final.16.d<-ddply(surv.demo.16.data, c("newplot","sbplt.unique.id"), plyr::summarize, 
                            alive.16 = sum(spring_survival_t1),
                            total_demo.15 = length(ID)) ## total number of demo plants in this is not necessarily 4 since it's only for the plants for which we don't know the endo stat



alive.demo.16.plot <- AGHY.demo.final.16.d$newplot ## new plot numbers that match with plant ID, to tie in plot probability 

N.demo.alive.16 <- length(alive.demo.16.plot) ## length of alive demography plants in 16 to loop over 

demo.alive.16<-AGHY.demo.final.16.d$alive.16 ## number of demo plants per subplot that SURVIVED (for which we don't know endo status)

N.samples.15.demo<-AGHY.demo.final.16.d$total_demo.15 ## total number of demo plants per subplot for which we DO NOT know endo status


#########################################################################################################
#########################################################################################################

## data for flowering probability
total.plants.14<-AGHY.merge$spring_flowering_recruit_count_per_subplot[AGHY.merge$year_t==2014]+AGHY.merge$spring_vegetative_recruit_count_per_subplot[AGHY.merge$year_t==2014]
obs.flowering.14<-AGHY.merge$spring_flowering_recruit_count_per_subplot[AGHY.merge$year_t==2014]
total.plants.15<-AGHY.merge$spring_flowering_recruit_count_per_subplot[AGHY.merge$year_t==2015]+AGHY.merge$spring_vegetative_recruit_count_per_subplot[AGHY.merge$year_t==2015]
obs.flowering.15<-AGHY.merge$spring_flowering_recruit_count_per_subplot[AGHY.merge$year_t==2015]
total.plants.16<-AGHY.merge$spring_flowering_recruit_count_per_subplot[AGHY.merge$year_t==2016]+AGHY.merge$spring_vegetative_recruit_count_per_subplot[AGHY.merge$year_t==2016]
obs.flowering.16<-AGHY.merge$spring_flowering_recruit_count_per_subplot[AGHY.merge$year_t==2016]

## deal with the NA problem
flowering14.NAs<-which(is.na(total.plants.14))
total.plants.14<-total.plants.14[-flowering14.NAs]
obs.flowering.14<-obs.flowering.14[-flowering14.NAs]
y.14.plot.flowering<-y.14.plot[-flowering14.NAs]
N.obs.14.flowering<-length(total.plants.14)

flowering15.NAs<-which(is.na(total.plants.15))
total.plants.15<-total.plants.15[-flowering15.NAs]
obs.flowering.15<-obs.flowering.15[-flowering15.NAs]
y.15.plot.flowering<-y.15.plot[-flowering15.NAs]
N.obs.15.flowering<-length(total.plants.15)

flowering16.NAs<-which(is.na(total.plants.16))
total.plants.16<-total.plants.16[-flowering16.NAs]
obs.flowering.16<-obs.flowering.16[-flowering16.NAs]
y.16.plot.flowering<-y.16.plot[-flowering16.NAs]
N.obs.16.flowering<-length(total.plants.16)

###########################################################################
######### NEW Analysis: Seed Production Estimates #########################
###########################################################################


###########################################################################
###########################################################################
## turn water status into integers FIGURE OUT HOW TO ADD IN WATER 
### FOR LINEAR REGRESSION with log of the seed mass
AGHY.seed.mass.count$water.seed.count<-as.integer(AGHY.seed.mass.count$water)

seed.mass.Ep.data<- subset(AGHY.seed.mass.count, Endo == 1)
water.seed.count.Ep <- seed.mass.Ep.data$water.seed.count
Ep.seed.count<-log(seed.mass.Ep.data$seed_count)
seed.mass.Ep<-log(seed.mass.Ep.data$seed_mass)
N.Ep.seed.count<-length(seed.mass.Ep)

seed.mass.Em.data<- subset(AGHY.seed.mass.count, Endo == 0)
water.seed.count.Em<- seed.mass.Em.data$water.seed.count
Em.seed.count<-log(seed.mass.Em.data$seed_count)
seed.mass.Em<-log(seed.mass.Em.data$seed_mass)
N.Em.seed.count<-length(seed.mass.Em)

mean.Em <-mean(seed.mass.Em.data$seed_mass)

hist((seed.mass.Em.data$seed_mass))
test<-lm((seed_mass)~seed_count,data=seed.mass.Em.data)
plot((seed_mass)~seed_count,data=seed.mass.Em.data);abline(test)
###########################################################################
###########################################################################

#### For gaussian model of seed mass from Original Plants

Ep.seed.mass.13.data<-subset(seed.mass.plant, endo == 1)
Ep.seed.mass.13<-log(Ep.seed.mass.13.data$seeds.plot)
seed.Ep.plot.13<-Ep.seed.mass.13.data$newplot
N.Ep.seed.13<-length(seed.Ep.plot.13)



Em.seed.mass.13.data<-subset(seed.mass.plant, endo == 0)
Em.seed.mass.13<-log(Em.seed.mass.13.data$seeds.plot)
seed.Em.plot.13<-Em.seed.mass.13.data$newplot
N.Em.seed.13<-length(seed.Em.plot.13)

###########################################################################
###########################################################################
#x.levels<-seq(0,1,0.01)
#N.x.levels<-length(x.levels)

jag.data<-list(N.trt=N.trt,
               water=water,
               N.plots=N.plots,
               #init_freq=init_freq,
               #y.14=y.14,
               #N.samples.14=N.samples.14,
               #N.obs.14=N.obs.14,
               #y.14.plot=y.14.plot,
               #N.samples.14.demo = N.samples.14.demo,
               #N.obs.14.demo = N.obs.14.demo,
               #y.15=y.15,
               #N.samples.15=N.samples.15,
               #N.obs.15=N.obs.15,
               #y.15.plot=y.15.plot,
               #y.16=y.16,
               #N.samples.16=N.samples.16,
               #N.obs.16=N.obs.16,
               #y.16.plot=y.16.plot,
               #s.e.neg.15 = s.e.neg.15,
               #s.e.pos.15 = s.e.pos.15,
               #s.e.pos.15.water = s.e.pos.15.water,
               #s.e.neg.15.water = s.e.neg.15.water,
               #N.obs.pos.surv.15 = N.obs.pos.surv.15,
               #N.obs.neg.surv.15 = N.obs.neg.surv.15,
               #s.Ep.known.15 = s.Ep.known.15, # survival for known Ep plants 2015
               #s.Em.known.15 = s.Em.known.15, # survival for known Em plants 2015
               #s.Ep.known.16 = s.Ep.known.16, # survival for known Ep plants 2016
               #s.Em.known.16 = s.Em.known.16, # survival for known Em plants 2016
               #N.Ep.known.15 = N.Ep.known.15, # number of known Ep plants to loop over 2015
               #N.Em.known.15 = N.Em.known.15, # number of known Em plants to loop over 2015
               #N.Ep.known.16 = N.Ep.known.16, # number of known Ep plants to loop over 2016
               #N.Em.known.16 = N.Em.known.16, # number of known Em plants to loop over 2016
               #s.Ep.known.15.plot = s.Ep.known.15.plot, # plot ID to index prob Ep survival 2015
               #s.Em.known.15.plot = s.Em.known.15.plot, # plot ID to index prob Em survival 2015
               #s.Ep.known.16.plot = s.Ep.known.16.plot, # plot ID to index prob Ep survival 2016
               #s.Em.known.16.plot = s.Em.known.16.plot, # plot ID to index prob Em survival 2016
               #demo.alive.15 = demo.alive.15, # number of demo plants per SUBPLOT with unknown endo stat that survived 2015
               #demo.alive.16 = demo.alive.16, # number of demo plants per SUBPLOT with unknown endo stat that survived 2016
               #alive.demo.15.plot = alive.demo.15.plot, # plot ID to index prob (unknown endo stat) survival 2015
               #alive.demo.16.plot = alive.demo.16.plot, # plot ID to index prob (unknown endo stat) survival 2016
               #N.demo.alive.15 = N.demo.alive.15, # length of total demography plants within subplots to loop over 2015
               #N.demo.alive.16 = N.demo.alive.16, # length of total demography plants within subplots to loop over 2016
               #N.samples.14.demo = N.samples.14.demo, # number of demography plants per subplot with unknown endo stat in 2014
               #N.samples.15.demo = N.samples.15.demo, # number of demography plants per subplot with unknown endo stat in 2015
               #total.plants.14=total.plants.14,
               #obs.flowering.14=obs.flowering.14,
               #y.14.plot.flowering=y.14.plot.flowering,
               #N.obs.14.flowering=N.obs.14.flowering,
               #total.plants.15=total.plants.15,
               #obs.flowering.15=obs.flowering.15,
               #y.15.plot.flowering=y.15.plot.flowering,
               #N.obs.15.flowering=N.obs.15.flowering,
               #total.plants.16=total.plants.16,
               #obs.flowering.16=obs.flowering.16,
               #y.16.plot.flowering=y.16.plot.flowering,
               #N.obs.16.flowering=N.obs.16.flowering,
               water.seed.count.Ep = water.seed.count.Ep, # water assignments for E+ lin reg
               water.seed.count.Em = water.seed.count.Em, # water assignments for E- lin reg
               Ep.seed.count = Ep.seed.count, # number of E+ seeds for lin reg
               Em.seed.count = Em.seed.count, # number of E- seeds for lin reg
               seed.mass.Ep = seed.mass.Ep, # E+ seed mass for lin reg
               seed.mass.Em = seed.mass.Em, # E- seed mass for lin reg
               N.Ep.seed.count = N.Ep.seed.count, # length of E+ individuals for lin reg
               N.Em.seed.count = N.Em.seed.count, # length  f E- individuals for lin reg
               Ep.seed.mass.13 = Ep.seed.mass.13, # mass of Ep seeds per plot 2013 (original plants)
               Em.seed.mass.13 = Em.seed.mass.13, # mass of Em seeds per plot 2013 (original plants)
               seed.Ep.plot.13 = seed.Ep.plot.13, # plot ID for indexing 2013 (original plants)
               seed.Em.plot.13 = seed.Em.plot.13, # plot ID for indexing 2013 (original plants)
               N.Ep.seed.13 = N.Ep.seed.13, # length of E+ plants for looping
               N.Em.seed.13 = N.Em.seed.13) # length of E- plants for looping
               #x.levels=x.levels,
               #N.x.levels=N.x.levels)

## Inits function
inits<-function(){list(intercept.ep = rnorm(2,0,2), # prior for Ep lin reg (intercept)
                       intercept.em = rnorm(2,0,2), # prior for Em lin reg (intercept)
                       slope.ep = runif(2), # Ep lin reg slope
                       slope.em = runif(2),
                       tau.sigma.count = runif(1), # Em lin reg slope
                       Ep_beta0_seed.13 = rnorm(2,0,2),
                       Em_beta0_seed.13 = rnorm(2,0,2))
                       
}

## Params to estimate
parameters<-c("intercept.ep","intercept.em",
                "slope.ep","slope.em",
                " Ep.orig.seed.control.13", "Ep.orig.seed.add.13",
                "Em.orig.seed.control.13","Em.orig.seed.add.13")


## MCMC settings
ni<-15000
nb<-5000
nt<-10
nc<-3

## run JAGS
AGHY.seedmass.out<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,model.file="AGHY_seed_mass.txt",
                          n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())

#"beta0.mean.14","beta1.14",
#"beta0.mean.15","beta1.15",
#"beta0.mean.16","beta1.16",
#"p.14","p.15","p.16",
#"Ep_beta_surv.15","Em_beta_surv.15",
#"Ep_pr_surv.15", "Em_pr_surv.15", # prob survival for Ep and Em plants 2015
#"Ep_pr_surv.16", "Em_pr_surv.16", # prob survival for Ep and Em plants 2016
#"pSurv.15", "pSurv.16", # prob survival for plants with unknown endo stat including the known prob survival
#"Ep.surv.control.15", "Ep.surv.add.15", # survival for Ep water treatments 2015
#"Em.surv.control.15","Em.surv.add.15", # survival for Em water treatments 2015
#"Ep.surv.control.16", "Ep.surv.add.16", # survival for Ep water treatments 2016
#"Em.surv.control.16", "Em.surv.add.16", # survival for Em water treatments 2016
#"Eplus.14.add.pred","Eplus.14.control.pred",
#"Eplus.15.add.pred","Eplus.15.control.pred",
#"Eplus.16.add.pred","Eplus.16.control.pred",
#"Ep.flower.control.14","Ep.flower.add.14",
#"Em.flower.control.14","Em.flower.add.14",
#"Ep.flower.control.15","Ep.flower.add.15",
#"Em.flower.control.15","Em.flower.add.15",
#"Ep.flower.control.16","Ep.flower.add.16",
#"Em.flower.control.16","Em.flower.add.16")
