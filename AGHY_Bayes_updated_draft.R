## Authors: Marion and Tom
## Purpose: Bayesian analysis of endophyte prevalence and individual-level vital rates
## Uses data imported and 'cleaned' in AGHY_data_manipulation.R
## Last update: 15 December 2017
# Last update notes: getting endo change into a single ggplot2 figure and still working on combining the linear regression of seed mass to seed count to the seed mass data to turn it into seed numbers
######################################################

## Required packages
library(R2jags)
library(mcmcplots)
library(plyr)
library(Cairo)
library(tidyverse)

## load workspace which has all of the data files already read in
load("AGHY_data_clean.1.RData")

## Bayesian model for endo frequency change and survival 

sink("AGHY_endochange.txt")
cat("
    model{

###################################################################
############### endo prevalence sub-model #####################
###################################################################

    ## Priors
    for(i in 1:N.trt){

    ## Priors for regression coefficients for change in endo prevalence -- unique coefficients for each year
    beta0.mean.14[i]~dnorm(0,0.001)   ##hyperprior (global mean and variance for the param -- used when there are random effects - hierarchical) for intercept beta0
    beta1.14[i]~dnorm(0,0.001)   ##prior for slope beta1
    beta0.mean.15[i]~dnorm(0,0.001)   ##hyperprior for intercept beta0
    beta1.15[i]~dnorm(0,0.001)   ##prior for slope beta1
    beta0.mean.16[i]~dnorm(0,0.001)   ##hyperprior for intercept beta0
    beta1.16[i]~dnorm(0,0.001)   ##prior for slope beta1
    }
    
    ## random effect of plot for prevalence regression, only affects intercept and specific for each year and trts
    sigma0~dunif(0,1000)
    tau.sigma0<-1/(sigma0*sigma0)

    for(i in 1:N.plots){      ##plot means
    ran.beta0.14[i]~dnorm(0,tau.sigma0)
    ran.beta0.15[i]~dnorm(0,tau.sigma0)
    ran.beta0.16[i]~dnorm(0,tau.sigma0)
    }
    
    ## Likelihood
    ## expected E+ freq for each plot -- probability of a year predicted by the hyper prior (global mean and variance), water treatment,
    #### the random effect for each of the plots for each year and the beta 1 for each year and water treatment, 
    #### multiplied by the e+ frequency from the previous year -- use the logit because the probability needs to be between 0 and 1
    ##### but the parameters that make up it are not constrained by 0 and 1
    for(i in 1:N.plots){
    logit(p.14[i])<-beta0.mean.14[water[i]]+ran.beta0.14[i]+beta1.14[water[i]]*init_freq[i]
    logit(p.15[i])<-beta0.mean.15[water[i]]+ran.beta0.15[i]+beta1.15[water[i]]*p.14[i]
    logit(p.16[i])<-beta0.mean.16[water[i]]+ran.beta0.16[i]+beta1.16[water[i]]*p.15[i]
    }
    
    #### Likelihood for endo frequency change through the years 
    ## sample likelihood oversubplots for 1314 transition
    for(i in 1:N.obs.14){
    y.14[i]~dbinom(p.14[y.14.plot[i]],N.samples.14[i])
    }
    ## sample likelihood oversubplots for 1415 transition
    for(i in 1:N.obs.15){
    y.15[i]~dbinom(p.15[y.15.plot[i]],N.samples.15[i])
    }
    ## sample likelihood oversubplots for 1415 transition
    for(i in 1:N.obs.16){
    y.16[i]~dbinom(p.16[y.16.plot[i]],N.samples.16[i])
    }

    ## Prediction
    for(i in 1:N.x.levels){
    Eplus.14.add.pred[i]<-exp(beta0.mean.14[1]+beta1.14[1]*x.levels[i])/(1+exp(beta0.mean.14[1]+beta1.14[1]*x.levels[i]))
    Eplus.14.control.pred[i]<-exp(beta0.mean.14[2]+beta1.14[2]*x.levels[i])/(1+exp(beta0.mean.14[2]+beta1.14[2]*x.levels[i]))
    Eplus.15.add.pred[i]<-exp(beta0.mean.15[1]+beta1.15[1]*x.levels[i])/(1+exp(beta0.mean.15[1]+beta1.15[1]*x.levels[i]))
    Eplus.15.control.pred[i]<-exp(beta0.mean.15[2]+beta1.15[2]*x.levels[i])/(1+exp(beta0.mean.15[2]+beta1.15[2]*x.levels[i]))
    Eplus.16.add.pred[i]<-exp(beta0.mean.16[1]+beta1.16[1]*x.levels[i])/(1+exp(beta0.mean.16[1]+beta1.16[1]*x.levels[i]))
    Eplus.16.control.pred[i]<-exp(beta0.mean.16[2]+beta1.16[2]*x.levels[i])/(1+exp(beta0.mean.16[2]+beta1.16[2]*x.levels[i]))
    }

###################################################################
############### flowering probability sub-model #####################
###################################################################

    ## Priors
    for(i in 1:N.trt){

    ## Priors for flowering probability
    Ep_beta0_pf.14[i]~dnorm(0,0.001)
    Em_beta0_pf.14[i]~dnorm(0,0.001)
    Ep_beta0_pf.15[i]~dnorm(0,0.001)
    Em_beta0_pf.15[i]~dnorm(0,0.001)
    Ep_beta0_pf.16[i]~dnorm(0,0.001)
    Em_beta0_pf.16[i]~dnorm(0,0.001)
    }

    ## random effect of plot for the probability of flowering 
    sigma.pf~dunif(0,1000)
    tau.sigma.pf<-1/(sigma.pf*sigma.pf)
    
    for(i in 1:N.plots){      ##plot means
    ran.pf.14[i]~dnorm(0,tau.sigma.pf) ## random effect of plot on prob flowering in 2014
    ran.pf.15[i]~dnorm(0,tau.sigma.pf) ## random effect of plot on prob flowering in 2015
    ran.pf.16[i]~dnorm(0,tau.sigma.pf) ## random effect of plot on prob flowering in 2016
    

    ## flowering probability, probability of flowering each year based on the 
    ## prior (Ep_beta0_pf), treatment ID (water addition or water control) and random effect of plot
    logit(Ep_pr_flow.14[i]) <- Ep_beta0_pf.14[water[i]] + ran.pf.14[i]
    logit(Em_pr_flow.14[i]) <- Em_beta0_pf.14[water[i]] + ran.pf.14[i]
    
    logit(Ep_pr_flow.15[i]) <- Ep_beta0_pf.15[water[i]] + ran.pf.15[i]
    logit(Em_pr_flow.15[i]) <- Em_beta0_pf.15[water[i]] + ran.pf.15[i]
    
    logit(Ep_pr_flow.16[i]) <- Ep_beta0_pf.16[water[i]] + ran.pf.16[i]
    logit(Em_pr_flow.16[i]) <- Em_beta0_pf.16[water[i]] + ran.pf.16[i]        
    }
    
    ## Estimation of endo-specific flowering probability
    for(i in 1:N.obs.14.flowering){
    pF.14[i] <- Ep_pr_flow.14[y.14.plot.flowering[i]]*p.14[y.14.plot.flowering[i]] + Em_pr_flow.14[y.14.plot.flowering[i]]*(1-p.14[y.14.plot.flowering[i]])
    obs.flowering.14[i]~dbinom(pF.14[i],total.plants.14[i])
    }
    for(i in 1:N.obs.15.flowering){
    pF.15[i] <- Ep_pr_flow.15[y.15.plot.flowering[i]]*p.15[y.15.plot.flowering[i]] + Em_pr_flow.15[y.15.plot.flowering[i]]*(1-p.15[y.15.plot.flowering[i]])
    obs.flowering.15[i]~dbinom(pF.15[i],total.plants.15[i])
    }
    for(i in 1:N.obs.16.flowering){
    pF.16[i] <- Ep_pr_flow.16[y.16.plot.flowering[i]]*p.16[y.16.plot.flowering[i]] + Em_pr_flow.16[y.16.plot.flowering[i]]*(1-p.16[y.16.plot.flowering[i]])
    obs.flowering.16[i]~dbinom(pF.16[i],total.plants.16[i])
    }
    
    ## derived quantity for endo-specific and water treatment flowering prob
    logit(Ep.flower.control.14)<-Ep_beta0_pf.14[2]
    logit(Ep.flower.add.14)<-Ep_beta0_pf.14[1]
    logit(Em.flower.control.14)<-Em_beta0_pf.14[2]
    logit(Em.flower.add.14)<-Em_beta0_pf.14[1]
    
    logit(Ep.flower.control.15)<-Ep_beta0_pf.15[2]
    logit(Ep.flower.add.15)<-Ep_beta0_pf.15[1]
    logit(Em.flower.control.15)<-Em_beta0_pf.15[2]
    logit(Em.flower.add.15)<-Em_beta0_pf.15[1]
    
    logit(Ep.flower.control.16)<-Ep_beta0_pf.16[2]
    logit(Ep.flower.add.16)<-Ep_beta0_pf.16[1]
    logit(Em.flower.control.16)<-Em_beta0_pf.16[2]
    logit(Em.flower.add.16)<-Em_beta0_pf.16[1]

###################################################################
############### survival probability sub-model #####################
###################################################################

    ## Priors
    for(i in 1:N.trt){
    ## Priors for survival probability - specific to endo and indexed for water treatment
    Ep_beta0_surv.15[i]~dnorm(0,0.001) ## 2014-2015 survival probabilty 
    Em_beta0_surv.15[i]~dnorm(0,0.001)
    Ep_beta0_surv.16[i]~dnorm(0,0.001) ## 2015-2016 survival probabilty 
    Em_beta0_surv.16[i]~dnorm(0,0.001)
    }

    ## random effect of plot for the probability of survival 
    sigma.surv~dunif(0,1000)
    tau.sigma.surv<-1/(sigma.surv*sigma.surv)
    
    for(i in 1:N.plots){      ##plot means
    ran.surv.15[i]~dnorm(0,tau.sigma.surv) ## random effect of plot on prob surv in 2015
    ran.surv.16[i]~dnorm(0,tau.sigma.surv) ## random effect of plot on prob surv in 2016
    ## apply random effect to each plot
    logit(Ep_pr_surv.15[i]) <- Ep_beta0_surv.15[water[i]] + ran.surv.15[i]
    logit(Em_pr_surv.15[i]) <- Em_beta0_surv.15[water[i]] + ran.surv.15[i]
    logit(Ep_pr_surv.16[i]) <- Ep_beta0_surv.16[water[i]] + ran.surv.16[i]
    logit(Em_pr_surv.16[i]) <- Em_beta0_surv.16[water[i]] + ran.surv.16[i] 
    }

    ################################################################
    ## Likelihood estimation for plants of known endo status
    for(i in 1:N.Ep.known.15){
    s.Ep.known.15[i] ~ dbern(Ep_pr_surv.15[s.Ep.known.15.plot[i]])
    }
    for(i in 1:N.Em.known.15){
    s.Em.known.15[i] ~ dbern(Em_pr_surv.15[s.Em.known.15.plot[i]])
    }
    for(i in 1:N.Ep.known.16){
    s.Ep.known.16[i] ~ dbern(Ep_pr_surv.16[s.Ep.known.16.plot[i]])
    }
    for(i in 1:N.Em.known.16){
    s.Em.known.16[i] ~ dbern(Em_pr_surv.16[s.Em.known.16.plot[i]])
    }

    ################################################################
    ## Likelihood estimation for plants of UNknown endo status
    ## This will be limited to plants that died before we ever knew their status
    ## We will estimate this from the plot frequency using the weighted mean approach

    ## loop over subplots. the data here are the total survivors per subplot
    for(i in 1:N.demo.alive.15){
    pSurv.15[i] <- Ep_pr_surv.15[alive.demo.15.plot[i]]*p.14[alive.demo.15.plot[i]] + Em_pr_surv.15[alive.demo.15.plot[i]]*(1-p.14[alive.demo.15.plot[i]])
    demo.alive.15[i]~dbinom(pSurv.15[i],N.samples.14.demo[i])
    }
    for(i in 1:N.demo.alive.16){
    pSurv.16[i] <- Ep_pr_surv.16[alive.demo.16.plot[i]]*p.15[alive.demo.16.plot[i]] + Em_pr_surv.16[alive.demo.16.plot[i]]*(1-p.15[alive.demo.16.plot[i]])
    demo.alive.16[i]~dbinom(pSurv.16[i],N.samples.15.demo[i])
    }


  ## derived quantity for endo-specific and water treatment survival prob
 
    logit(Ep.surv.control.15)<-Ep_beta0_surv.15[2]
    logit(Ep.surv.add.15)<-Ep_beta0_surv.15[1]
    logit(Em.surv.control.15)<-Em_beta0_surv.15[2]
    logit(Em.surv.add.15)<-Em_beta0_surv.15[1]
    
    logit(Ep.surv.control.16)<-Ep_beta0_surv.16[2]
    logit(Ep.surv.add.16)<-Ep_beta0_surv.16[1]
    logit(Em.surv.control.16)<-Em_beta0_surv.16[2]
    logit(Em.surv.add.16)<-Em_beta0_surv.16[1]

#####################################################################
########## NEW: Seed Production July 16, 2017 #######################
#####################################################################

###################################################################################
## Having trouble getting the estimates to work for the water treatment indexing ##
###################################################################################
## First start with change from seed mass to seed counts

## priors for regression paramaters endo specific and indexed over water treatments 
#for (i in 1:N.trt){
#    intercept.ep[i] ~ dnorm(0, 0.001)
#    slope.ep[i] ~ dunif(0, 1000) ## don't want a negative slope, using the uniform distribution to have the min be 0  
#    
#  intercept.em[i] ~ dnorm(0, 0.001)
#   slope.em[i] ~ dunif(0, 1000) ## don't want a negative slope, using the uniform distribution to have the min be 0  
#    
#}
#
### sigma and tau for seed counts
#  #sigma.count ~ dunif(0, 1000)
#tau.sigma.count ~ dgamma(0.001, 0.001)



### Likelihood model for number of seeds and seed mass
#for (i in 1:N.Ep.seed.count){
#mean.seed.count.Ep[i] <- intercept.ep[water.seed.count.Ep[i]] + slope.ep[water.seed.count.Ep[i]]*seed.mass.Ep[i]
#Ep.seed.count[i] ~ dnorm(mean.seed.count.Ep[i], tau.sigma.count)

#}

#for (i in 1:N.Em.seed.count){
#mean.seed.count.Em[i]<- intercept.em[water.seed.count.Em[i]] + (slope.em[water.seed.count.Em[i]])*seed.mass.Em[i]
#Em.seed.count[i] ~ dnorm(mean.seed.count.Em[i], tau.sigma.count)

#}


## Seed production for the original plants in 2013 (these seeds were produced by original plants in 2013 and will contribute to the recruits in 2014)

### Priors for seed production - specific for endo, indexed for water treatment
#  for(i in 1:N.trt){
#  Ep_beta0_seed.13[i] ~ dnorm(0, 0.001)
#  Em_beta0_seed.13[i] ~ dnorm(0, 0.001)
#  }


## Tau.sigma for prob seed production plot (this tau is different than the tau that will go in the gaussian model)
#  sigma.seed ~ dunif(0, 1000)
#  tau.sigma.seed <- 1/(sigma.seed * sigma.seed)

## Tau.sigma for the model
#  sigma.s ~ dunif(0,1000)
#  tau.sigma.s <- 1/(sigma.s * sigma.s)

## Random effect of plot for seed production
#  for(i in 1:N.plots){
#  ran.seed.13[i] ~ dnorm(0, tau.sigma.seed)
  

## Apply random effect of plot on seed production 
#  Ep_seed.13[i]<- Ep_beta0_seed.13[water[i]] + ran.seed.13[i]
#  Em_seed.13[i]<- Em_beta0_seed.13[water[i]] + ran.seed.13[i]
#  }



## Likelihood estimate for seed production from E+ and E- plants
#  for(i in 1:N.Ep.seed.13){
## eventually need to multiply the mean seed number per plot by the slope from the linear regression to turn the seed mass into seed count
#  Ep.seed.mass.13[i] ~ dnorm(Ep_seed.13[seed.Ep.plot.13[i]], tau.sigma.s)
  
#  }

#  for(i in 1:N.Em.seed.13){
## eventually need to multiply the mean seed number per plot by the slope from the linear regression to turn the seed mass into seed count
#    Em.seed.mass.13[i] ~ dnorm(Em_seed.14[seed.Em.plot.13[i]], tau.sigma.s)
    
#    }


  ## derived quantity for endo-specific and water treatment seed mass prob
 
#    Ep.seed.control.13<-Ep_beta0_seed.14[2]*slope.ep[2]
#    Ep.seed.add.13<-Ep_beta0_seed.14[1]* slope.ep[1]
#    Em.seed.control.13<-Em_beta0_seed.14[2]*slope.em[2]
#    Em.seed.add.13<-Em_beta0_seed.14[1]*slope.ep[1]
    


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
### FOR LINEAR REGRESSION 
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
###########################################################################
###########################################################################

#### For gaussian model of seed mass 

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
x.levels<-seq(0,1,0.01)
N.x.levels<-length(x.levels)

jag.data<-list(N.trt=N.trt,
               water=water,
               N.plots=N.plots,
               init_freq=init_freq,
               y.14=y.14,
               N.samples.14=N.samples.14,
               N.obs.14=N.obs.14,
               y.14.plot=y.14.plot,
              #N.samples.14.demo = N.samples.14.demo,
               #N.obs.14.demo = N.obs.14.demo,
               y.15=y.15,
               N.samples.15=N.samples.15,
               N.obs.15=N.obs.15,
               y.15.plot=y.15.plot,
               y.16=y.16,
               N.samples.16=N.samples.16,
               N.obs.16=N.obs.16,
               y.16.plot=y.16.plot,
               #s.e.neg.15 = s.e.neg.15,
               #s.e.pos.15 = s.e.pos.15,
               #s.e.pos.15.water = s.e.pos.15.water,
               #s.e.neg.15.water = s.e.neg.15.water,
               #N.obs.pos.surv.15 = N.obs.pos.surv.15,
               #N.obs.neg.surv.15 = N.obs.neg.surv.15,
               s.Ep.known.15 = s.Ep.known.15, # survival for known Ep plants 2015
               s.Em.known.15 = s.Em.known.15, # survival for known Em plants 2015
               s.Ep.known.16 = s.Ep.known.16, # survival for known Ep plants 2016
               s.Em.known.16 = s.Em.known.16, # survival for known Em plants 2016
               N.Ep.known.15 = N.Ep.known.15, # number of known Ep plants to loop over 2015
               N.Em.known.15 = N.Em.known.15, # number of known Em plants to loop over 2015
               N.Ep.known.16 = N.Ep.known.16, # number of known Ep plants to loop over 2016
               N.Em.known.16 = N.Em.known.16, # number of known Em plants to loop over 2016
               s.Ep.known.15.plot = s.Ep.known.15.plot, # plot ID to index prob Ep survival 2015
               s.Em.known.15.plot = s.Em.known.15.plot, # plot ID to index prob Em survival 2015
               s.Ep.known.16.plot = s.Ep.known.16.plot, # plot ID to index prob Ep survival 2016
               s.Em.known.16.plot = s.Em.known.16.plot, # plot ID to index prob Em survival 2016
               demo.alive.15 = demo.alive.15, # number of demo plants per SUBPLOT with unknown endo stat that survived 2015
               demo.alive.16 = demo.alive.16, # number of demo plants per SUBPLOT with unknown endo stat that survived 2016
               alive.demo.15.plot = alive.demo.15.plot, # plot ID to index prob (unknown endo stat) survival 2015
               alive.demo.16.plot = alive.demo.16.plot, # plot ID to index prob (unknown endo stat) survival 2016
               N.demo.alive.15 = N.demo.alive.15, # length of total demography plants within subplots to loop over 2015
               N.demo.alive.16 = N.demo.alive.16, # length of total demography plants within subplots to loop over 2016
               N.samples.14.demo = N.samples.14.demo, # number of demography plants per subplot with unknown endo stat in 2014
               N.samples.15.demo = N.samples.15.demo, # number of demography plants per subplot with unknown endo stat in 2015
               total.plants.14=total.plants.14,
               obs.flowering.14=obs.flowering.14,
               y.14.plot.flowering=y.14.plot.flowering,
               N.obs.14.flowering=N.obs.14.flowering,
               total.plants.15=total.plants.15,
               obs.flowering.15=obs.flowering.15,
               y.15.plot.flowering=y.15.plot.flowering,
               N.obs.15.flowering=N.obs.15.flowering,
               total.plants.16=total.plants.16,
               obs.flowering.16=obs.flowering.16,
               y.16.plot.flowering=y.16.plot.flowering,
               N.obs.16.flowering=N.obs.16.flowering,
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
               N.Em.seed.13 = N.Em.seed.13, # length of E- plants for looping
               x.levels=x.levels,
               N.x.levels=N.x.levels)

## Inits function
inits<-function(){list(beta0.mean.14=rnorm(2,0,2),
                       beta1.14=rnorm(2,0,2),
                       beta0.mean.15=rnorm(2,0,2),
                       beta1.15=rnorm(2,0,2),
                       beta0.mean.16=rnorm(2,0,2),
                       beta1.16=rnorm(2,0,2),
                       sigma0=rlnorm(1),
                       Ep_beta0_surv.15=rnorm(2,0,2), # prior for Ep survival 2015
                       Em_beta0_surv.15=rnorm(2,0,2), # prior for Em survival 2015
                       Ep_beta0_surv.16=rnorm(2,0,2), # prior for Ep survival 2016
                       Em_beta0_surv.16=rnorm(2,0,2), # prior for Em survival 2016
                       Em_beta0_pf.14=rnorm(2,0,2),
                       Ep_beta0_pf.14=rnorm(2,0,2),
                       Em_beta0_pf.15=rnorm(2,0,2),
                       Ep_beta0_pf.15=rnorm(2,0,2),
                       Em_beta0_pf.16=rnorm(2,0,2),
                       Ep_beta0_pf.16=rnorm(2,0,2),
                       intercept.ep = rnorm(2,0,2), # prior for Ep lin reg (intercept)
                       intercept.em = rnorm(2,0,2), # prior for Em lin reg (intercept)
                       slope.ep = runif(2), # Ep lin reg slope
                       slope.em = runif(2),
                       tau.sigma.count = runif(1), # Em lin reg slope
                       Ep_beta0_seed.14 = rnorm(2,0,2),
                       Em_beta0_seed.14 = rnorm(2,0,2))
}

## Params to estimate
parameters<-c("beta0.mean.14","beta1.14",
              "beta0.mean.15","beta1.15",
              "beta0.mean.16","beta1.16",
              "p.14","p.15","p.16",
              #"Ep_beta_surv.15","Em_beta_surv.15",
              #"Ep_pr_surv.15", "Em_pr_surv.15", # prob survival for Ep and Em plants 2015
              #"Ep_pr_surv.16", "Em_pr_surv.16", # prob survival for Ep and Em plants 2016
              #"pSurv.15", "pSurv.16", # prob survival for plants with unknown endo stat including the known prob survival
              "Ep.surv.control.15", "Ep.surv.add.15", # survival for Ep water treatments 2015
              "Em.surv.control.15","Em.surv.add.15", # survival for Em water treatments 2015
              "Ep.surv.control.16", "Ep.surv.add.16", # survival for Ep water treatments 2016
              "Em.surv.control.16", "Em.surv.add.16", # survival for Em water treatments 2016
              "Eplus.14.add.pred","Eplus.14.control.pred",
              "Eplus.15.add.pred","Eplus.15.control.pred",
              "Eplus.16.add.pred","Eplus.16.control.pred",
              "Ep.flower.control.14","Ep.flower.add.14",
              "Em.flower.control.14","Em.flower.add.14",
              "Ep.flower.control.15","Ep.flower.add.15",
              "Em.flower.control.15","Em.flower.add.15",
              "Ep.flower.control.16","Ep.flower.add.16",
              "Em.flower.control.16","Em.flower.add.16")
              #"intercept.ep","intercept.em",
              #"slope.ep","slope.em",
              #"Ep.seed.control.14", "Ep.seed.add.14",
              #"Em.seed.control.14","Em.seed.add.14")


## MCMC settings
ni<-15000
nb<-5000
nt<-10
nc<-3

## run JAGS
AGHY.endochange.out<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,model.file="AGHY_endochange.txt",
                          n.thin=nt,n.chains=nc,n.burnin=nb,n.iter=ni,DIC=T,working.directory=getwd())


df<-AGHY.endochange.out$BUGSoutput$summary

#mcmcplot(AGHY.endochange.out)


### Fiddling with figures -- this will be it's own script later

bayes.data.summary1<-as.data.frame(AGHY.endochange.out$BUGSoutput$summary, keep.rownames = TRUE)


## create a column in the dataframe that has the rownames
bayes.data.summary1$names<- rownames(bayes.data.summary1)


## separate the row names so that we can identify the estimated values by "type"
#### ULTIMATELY WANT TO RENAME THE ENDO CHANGE OUTPUT SO IT'S IN THE SAME FORMAT AS THE OTHERS....
bayes.endo<-bayes.data.summary1 %>%
  separate(names, c("endo", "type_year", "water", "year"), extra = "merge", fill = "left")


## make a dataframe of just the estimates for change in endophyte prevalence within each transition year
bayes.endo.prev<- bayes.endo %>%
  filter(endo == "Eplus") %>%
  rename(low = "2.5%", high = "97.5%")


## make an ordering column for the estimated prevalences 
vec<-rep(x.levels, 6)

bayes.endo.prev$order<-vec

## append latent state estimates of plot frequency to AGHY.plots data frame
AGHY.plots$p.14<-AGHY.endochange.out$BUGSoutput$mean$p.14
AGHY.plots$p.15<-AGHY.endochange.out$BUGSoutput$mean$p.15
AGHY.plots$p.16<-AGHY.endochange.out$BUGSoutput$mean$p.16

## create column for p.13 -- which is the actual initial starting prevalence of the plots
AGHY.plots$p.13 <-AGHY.plots$target_init_freq

## gather the AGHY.plots dataframe to turn it into long format - tidy data
AGHY.plots.gathered<-AGHY.plots %>%
  gather("p.13","p.14","p.15", "p.16", key = "transition_yr", value = "prob")

prob_t1<-AGHY.plots.gathered %>%
  filter(transition_yr != "p.13") %>%
  select(transition_yr, prob) %>%
  rename(yr_t1 = transition_yr, prob_t1 = prob)

AGHY.plots.all<- AGHY.plots.gathered %>%
  rename(yr_t = transition_yr, prob_t = prob) %>% ## drop the 2016 probability estimates for year t - since 2016 is our last year
  filter(yr_t != "p.16")

AGHY.plots.all <- cbind(AGHY.plots.all, prob_t1)

## create a matching column (yr_t) in the endo df

bayes.endo.prev$p<- "p."

bayes.endo.prev$yr_t1 <- paste(bayes.endo.prev$p,bayes.endo.prev$type_year, sep="")


## add in the endo.prev dataframe
AGHY.endo.all <- full_join(AGHY.plots.all, bayes.endo.prev, by= "yr_t1")



## 
## Plotting change in endophyte prevalence within the populations by water treatment
ggplot(AGHY.endo.all, aes(order, order))+
  geom_point(aes(prob_t, prob_t1))+
  geom_line(aes(group= yr_t))+
  geom_ribbon(aes(x = order, ymin=low, ymax= high, group = yr_t),alpha=0.3) +
  facet_grid(yr_t~water.x)
  

h + geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
  geom_line(aes(y = level))
  
  geom_path(data = AGHY.endo.all, aes(prob_t,mean))+
 geom_ribbon(aes(ymin=low, ymax= high),alpha=0.3)+
  facet_grid(water.x~yr_t)


geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
  geom_line(aes(y = level))

geom_line(data=bayes.data.mean.control.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.control.pred), size =1.5, col="red")+
  geom_ribbon(data=bayes.data.summary.control.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+


## Splitting up the water and control data for comparison with the ggplot2 plots -- just making sure they're the same
plot(x.levels,x.levels,type="n")
lines(x.levels,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.add.pred,lwd=4,col="blue")
lines(x.levels,AGHY.endochange.out$BUGSoutput$summary[33:133,"2.5%"],lwd=1,col="blue")
lines(x.levels,AGHY.endochange.out$BUGSoutput$summary[33:133,"97.5%"],lwd=1,col="blue")

abline(0,1,col="gray")
points(AGHY.plots$target_init_freq[AGHY.plots$water=="Add"],
       AGHY.plots$p.14[AGHY.plots$water=="Add"],col="blue",cex=2)

plot(x.levels,x.levels,type="n")
lines(x.levels,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.control.pred,lwd=4,col="red")
lines(x.levels,AGHY.endochange.out$BUGSoutput$summary[134:234,"2.5%"],lwd=1,col="red")
lines(x.levels,AGHY.endochange.out$BUGSoutput$summary[134:234,"97.5%"],lwd=1,col="red")
points(AGHY.plots$target_init_freq[AGHY.plots$water=="Control"],
       AGHY.plots$p.14[AGHY.plots$water=="Control"],col="red",cex=2)

## Selecting out the water addition and control plots for ggplot2 plots
AGHY.plots.control<-filter(AGHY.plots, water == "Control")
AGHY.plots.add<-filter(AGHY.plots, water == "Add")

## Creating dataframes for ggplot2 plotting 2014
bayes.data.mean.control.14<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$Eplus.14.control.pred)
bayes.data.mean.add.14<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$Eplus.14.add.pred)
bayes.data.summary<-as.data.frame(AGHY.endochange.out$BUGSoutput$summary)
bayes.data.summary.control.14<-bayes.data.summary[c(134:234), c(3,7)]
bayes.data.summary.add.14<-bayes.data.summary[c(33:133), c(3,7)]

### Thinking about plotting the demography plant E+ output -- but not sure if it makes sense
#bayes.demo.e.pos<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$N.e.pos.demo.14)
#bayes.data.summary.demo<-as.data.frame(AGHY.endochange.out$BUGSoutput$summary)
#bayes.data.summary.demo<-as.data.frame(bayes.data.summary.demo[c(607:741),c(3,7)])
#bayes.data.summary.demo$mean.demo<-AGHY.endochange.out$BUGSoutput$mean$N.e.pos.demo.14
#bayes.data.summary.demo$num.demo<-AGHY.merge.recruits.14.d$total_demo.14

#x.levels.demo<-seq(0,4,.5)
#ggplot(bayes.data.summary.demo, aes(x.levels.demo,x.levels.demo))+
#  geom_line(data = bayes.data.summary.demo, aes(, mean.demo))

## Plotting just the water control 2014
ggplot(bayes.data.mean.control.14, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(target_init_freq, p.14))+
  geom_line(data=bayes.data.mean.control.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.control.pred))+
  geom_ribbon(data=bayes.data.summary.control.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3)

## Plotting just the water addition 2014
ggplot(bayes.data.mean.add.14, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.add, aes(target_init_freq, p.14))+
  geom_line(data=bayes.data.mean.add.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.add.pred))+
  geom_ribbon(data=bayes.data.summary.add.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3)

## Combined plot for water control and addition for 2014
p14<-ggplot(bayes.data.mean.control.14, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(target_init_freq, p.14), col="red")+
  geom_line(data=bayes.data.mean.control.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.control.pred), size =1.5, col="red")+
  geom_ribbon(data=bayes.data.summary.control.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  geom_point(data = AGHY.plots.add, aes(target_init_freq, p.14), col="blue")+
  geom_line(data=bayes.data.mean.add.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.add.pred), size=1.5, col="blue")+
  geom_ribbon(data=bayes.data.summary.add.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+ 2013", y = "E+ 2014")+
  theme_bw() + theme(panel.border = element_blank())


## Combined plot for water control and addition for 2014
p14.control<-ggplot(bayes.data.mean.control.14, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(target_init_freq, p.14), col="orange3")+
  geom_line(data=bayes.data.mean.control.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.control.pred), size =1.5, col="orange3")+
  geom_ribbon(data=bayes.data.summary.control.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="orange", fill = "orange")+
  #geom_point(data = AGHY.plots.add, aes(target_init_freq, p.14), col="blue")+
  #geom_line(data=bayes.data.mean.add.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.add.pred), size=1.5, col="blue")+
  #geom_ribbon(data=bayes.data.summary.add.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+ 2013", y = "E+ 2014")#+
#theme_bw() + theme(panel.border = element_blank())


p14.add<-ggplot(bayes.data.mean.control.14, aes(x.levels, x.levels))+
  #geom_point(data = AGHY.plots.control, aes(target_init_freq, p.14), col="red")+
  #geom_line(data=bayes.data.mean.control.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.control.pred), size =1.5, col="red")+
  #geom_ribbon(data=bayes.data.summary.control.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  geom_point(data = AGHY.plots.add, aes(target_init_freq, p.14), col="steelblue4")+
  geom_line(data=bayes.data.mean.add.14, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.14.add.pred), size=1.5, col="steelblue")+
  geom_ribbon(data=bayes.data.summary.add.14,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="steelblue", fill="steelblue")+
  geom_abline()+
  labs(x = "E+ 2013", y = "E+ 2014")#+
#theme_bw() + theme(panel.border = element_blank())





## Creating dataframes for ggplot2 plotting 2015
bayes.data.mean.control.15<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$Eplus.15.control.pred)
bayes.data.mean.add.15<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$Eplus.15.add.pred)
bayes.data.summary<-as.data.frame(AGHY.endochange.out$BUGSoutput$summary)
bayes.data.summary.control.15<-bayes.data.summary[c(320:420), c(3,7)]
bayes.data.summary.add.15<-bayes.data.summary[c(219:319), c(3,7)]

## Combined plot for water control and addition for 2015
p15<-ggplot(bayes.data.mean.control.15, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(p.14, p.15), col="red")+
  geom_line(data=bayes.data.mean.control.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.control.pred), size =1.5, col="red")+
  geom_ribbon(data=bayes.data.summary.control.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  geom_point(data = AGHY.plots.add, aes(p.14, p.15), col="blue")+
  geom_line(data=bayes.data.mean.add.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.add.pred), size=1.5, col="blue")+
  geom_ribbon(data=bayes.data.summary.add.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+  2014", y = "E+ 2015")#+
#theme_bw() + theme(panel.border = element_blank())


p15.control<-ggplot(bayes.data.mean.control.15, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(p.14, p.15), col="red")+
  geom_line(data=bayes.data.mean.control.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.control.pred), size =1.5, col="red")+
  geom_ribbon(data=bayes.data.summary.control.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  #geom_point(data = AGHY.plots.add, aes(p.14, p.15), col="blue")+
  #geom_line(data=bayes.data.mean.add.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.add.pred), size=1.5, col="blue")+
  #geom_ribbon(data=bayes.data.summary.add.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+  2014", y = "E+ 2015")#+
#theme_bw() + theme(panel.border = element_blank())

p15.add<-ggplot(bayes.data.mean.control.15, aes(x.levels, x.levels))+
  #geom_point(data = AGHY.plots.control, aes(p.14, p.15), col="red")+
  #geom_line(data=bayes.data.mean.control.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.control.pred), size =1.5, col="red")+
  #eom_ribbon(data=bayes.data.summary.control.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  geom_point(data = AGHY.plots.add, aes(p.14, p.15), col="blue")+
  geom_line(data=bayes.data.mean.add.15, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.15.add.pred), size=1.5, col="blue")+
  geom_ribbon(data=bayes.data.summary.add.15,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+  2014", y = "E+ 2015")#+
#theme_bw() + theme(panel.border = element_blank())



## Creating dataframes for ggplot2 plotting 2016
bayes.data.mean.control.16<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$Eplus.16.control.pred)
bayes.data.mean.add.16<-as.data.frame(AGHY.endochange.out$BUGSoutput$mean$Eplus.16.add.pred)
bayes.data.summary<-as.data.frame(AGHY.endochange.out$BUGSoutput$summary)
bayes.data.summary.control.16<-bayes.data.summary[c(522:622), c(3,7)]
bayes.data.summary.add.16<-bayes.data.summary[c(421:521), c(3,7)]

## Combined plot for water control and addition 2016
p16<-ggplot(bayes.data.mean.control.16, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(p.15, p.16), col="red")+
  geom_line(data=bayes.data.mean.control.16, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.16.control.pred), size =1.5, col="red")+
  geom_ribbon(data=bayes.data.summary.control.16,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  geom_point(data = AGHY.plots.add, aes(p.15, p.16), col="blue")+
  geom_line(data=bayes.data.mean.add.16, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.16.add.pred), size=1.5, col="blue")+
  geom_ribbon(data=bayes.data.summary.add.16,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+ 2015", y = "E+ 2016")#+
#theme_bw() + theme(panel.border = element_blank())

p16.control<-ggplot(bayes.data.mean.control.16, aes(x.levels, x.levels))+
  geom_point(data = AGHY.plots.control, aes(p.15, p.16), col="red")+
  geom_line(data=bayes.data.mean.control.16, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.16.control.pred), size =1.5, col="red")+
  geom_ribbon(data=bayes.data.summary.control.16,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  #geom_point(data = AGHY.plots.add, aes(p.15, p.16), col="blue")+
  #geom_line(data=bayes.data.mean.add.16, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.16.add.pred), size=1.5, col="blue")+
  #geom_ribbon(data=bayes.data.summary.add.16,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+ 2015", y = "E+ 2016")
#theme_bw() + theme(panel.border = element_blank())

p16.add<-ggplot(bayes.data.mean.control.16, aes(x.levels, x.levels))+
  #geom_point(data = AGHY.plots.control, aes(p.15, p.16), col="red")+
  #geom_line(data=bayes.data.mean.control.16, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.16.control.pred), size =1.5, col="red")+
  #geom_ribbon(data=bayes.data.summary.control.16,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="red", fill = "red")+
  geom_point(data = AGHY.plots.add, aes(p.15, p.16), col="blue")+
  geom_line(data=bayes.data.mean.add.16, aes(,AGHY.endochange.out$BUGSoutput$mean$Eplus.16.add.pred), size=1.5, col="blue")+
  geom_ribbon(data=bayes.data.summary.add.16,aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.3, col="blue", fill="blue")+
  geom_abline()+
  labs(x = "E+ 2015", y = "E+ 2016")
#theme_bw() + theme(panel.border = element_blank())










bayes.data.seeds<- bayes.data.summary[c(7,8,19,20), c(1,3,7)]


############################################################################################
## 2015 survival based on endo status and water treatment  -- create objects
### Water Addition mean survival for E- and E+ with confidence interval
bayes.data.summary<-as.data.frame(AGHY.endochange.out$BUGSoutput$summary, keep.rownames=TRUE)[]
## pull out just the survival information 
bayes.data.summary.survival.16<-bayes.data.summary[c(8,10, 18,20),c(1,3,7)]

bayes.data.summary.survival.15<-bayes.data.summary[c(7,9, 17,19),c(1,3,7)]



bayes.data.summary.survival.15<-cbind(rownames(bayes.data.summary.survival.15),bayes.data.summary.survival.15)
names(bayes.data.summary.survival.15)[names(bayes.data.summary.survival.15) =="rownames(bayes.data.summary.survival.15)"] <- "treatment"


bayes.data.summary.survival.16<-cbind(rownames(bayes.data.summary.survival.16),bayes.data.summary.survival.16)
names(bayes.data.summary.survival.16)[names(bayes.data.summary.survival.16) =="rownames(bayes.data.summary.survival.16)"] <- "treatment"




water.treat<-c("Irrigated", "Ambient", "Irrigated", "Ambient")

bayes.data.summary.survival.15.water<-cbind(bayes.data.summary.survival.15,water.treat)

bayes.data.summary.survival.16.water<-cbind(bayes.data.summary.survival.16,water.treat)

## Change this to a factor ordered by water treatment so that ggplot2 can recognize it and plot it in the order I want
bayes.data.summary.survival.15.water$treatment<- factor(bayes.data.summary.survival.15.water$treatment, 
                                                        levels = bayes.data.summary.survival.15.water$treatment[order(bayes.data.summary.survival.15.water$water.treat)])


bayes.data.summary.survival.16.water$treatment<- factor(bayes.data.summary.survival.16.water$treatment, 
                                                        levels = bayes.data.summary.survival.16.water$treatment[order(bayes.data.summary.survival.16.water$water.treat)])



survival.all<-rbind(bayes.data.summary.survival.15.water, bayes.data.summary.survival.16.water)
endo.stat<-c(0,0,1,1,0,0,1,1)
survival.all$endo.stat<-endo.stat
survival.all$year<-c(2015,2015,2015,2015,2016,2016,2016,2016)

survival.all$endo.water<-paste(survival.all$water.treat, "_", survival.all$endo.stat)

s.all.plot<- ggplot(survival.all, aes(endo.water, mean, colour = water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
                                                  "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Surival", y = "Probability of Survival")+
  scale_y_continuous(limits = c(0,1))+  
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(~year)




#### Plot the mean probability of survival for E+ and E- plants by water and control treatments with the 2.5% and 97.5% Credible Intervals
p.surv.15<-ggplot(bayes.data.summary.survival.15.water, aes(treatment, mean, colour=water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("Symbiont Presence", labels = c("Em.surv.add.15" = "E+","Em.surv.control.15" = "E-",
                                                  "Ep.surv.add.15" = "E+","Ep.surv.control.15" = "E-"))+
  labs(title = "2015 Survival", y = "Probability of Survival")+
  scale_y_continuous(limits = c(0,1))+  
  guides(color=guide_legend("Watering Regime"))+

  theme_bw() #+ theme(panel.border = element_blank())



#### Plot the mean probability of survival for E+ and E- plants by water and control treatments with the 2.5% and 97.5% Credible Intervals
p.surv.16<-ggplot(bayes.data.summary.survival.16.water, aes(treatment, mean, colour=water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("Symbiont Presence", labels = c("Em.surv.add.16" = "-","Em.surv.control.16" = "-",
                                                  "Ep.surv.add.16" = "+","Ep.surv.control.16" = "+"))+
  labs(title = "2016 Survival", y = "Probability of Survival")+
  scale_y_continuous(limits = c(0,1))+  
  theme_bw()+
  guides(color=guide_legend("Watering Regime"))
#+ theme(panel.border = element_blank())



###################################### ###################################### 
###################################### FLOWERING ###################################### 
###################################### ###################################### 
flowering.all<-bayes.data.summary[c(1:6,11:16),c(1,3,7)]





flowering.all<-cbind(rownames(flowering.all),flowering.all)
names(flowering.all)[names(flowering.all) =="rownames(flowering.all)"] <- "treatment"



water.treat<-c("Irrigated","Irrigated","Irrigated", "Ambient","Ambient","Ambient", "Irrigated","Irrigated","Irrigated", "Ambient","Ambient","Ambient")
year<-c("2014","2015","2016","2014","2015","2016","2014","2015","2016","2014","2015","2016")
flowering.all<-cbind(flowering.all,water.treat,year)



## Change this to a factor ordered by water treatment so that ggplot2 can recognize it and plot it in the order I want
flowering.all$treatment<- factor(flowering.all$treatment, 
                                                        levels = flowering.all$treatment[order(flowering.all$water.treat)])


## add in column for endo stat 
endo.stat<-c(0,0,0,0,0,0,1,1,1,1,1,1)

flowering.all$endo.stat<-endo.stat
flowering.all$endo.water<-paste(flowering.all$water.treat,"_",flowering.all$endo.stat)

## trying the figures again with facet_wrap
f.all.plot<- ggplot(flowering.all, aes(endo.water, mean, colour = water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
 scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
                                                 "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Flowering", y = "Probability of Flowering")+
  scale_y_continuous(limits = c(0,1))+  
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(~year)



flowering.14<-subset(flowering.all, year == 2014)
flowering.15<-subset(flowering.all, year == 2015)
flowering.16<-subset(flowering.all, year == 2016)

#### Plot the mean probability of flowering for E+ and E- plants by water and control treatments with the 2.5% and 97.5% Confidence Intervals
p.flow.14<-ggplot(flowering.14, aes(treatment, mean, colour=water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("Endophyte Status", labels = c("Em.flower.add.14" = "E-","Em.flower.control.14" = "E-",
                                                  "Ep.flower.add.14" = "E+","Ep.flower.control.14" = "E+"))+
  labs(title = "2014 Flowering", y = "Probability of Flowering")+
  scale_y_continuous(limits = c(0,1))+  
  guides(color=guide_legend("Watering Regime"))+
  
  theme_bw() #+ theme(panel.border = element_blank())


p.flow.15<-ggplot(flowering.15, aes(treatment, mean, colour=water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("", labels = c("Em.flower.add.15" = "Non-host","Em.flower.control.15" = "Non-host",
                                                  "Ep.flower.add.15" = "Host","Ep.flower.control.15" = "Host"))+
  labs(title = "2015 Flowering", y = "Probability of Flowering")+
  scale_y_continuous(limits = c(0,1))+  
  guides(color=guide_legend("Watering Regime"))+
  
  theme_bw() #+ theme(panel.border = element_blank())


p.flow.16<-ggplot(flowering.16, aes(treatment, mean, colour=water.treat))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("", labels = c("Em.flower.add.16" = "Non-host","Em.flower.control.16" = "Non-Host",
                                                  "Ep.flower.add.16" = "Host","Ep.flower.control.16" = "Host"))+
  labs(title = "2016 Flowering", y = "Probability of Flowering")+
  scale_y_continuous(limits = c(0,1))+  
  #theme(legend.position="none")+
  guides(color=guide_legend("Watering Regime"))+
  
  theme_bw() #+ theme(panel.border = element_blank())


#### save workspace to call for RMarkdown file 
save.image("AGHY_Bayes.RData")

CairoPDF("flowering.pdf", width = 20, height = 5, bg='transparent')
grid.arrange(p16.control, p16.add, nrow=2, ncol=1) 


dev.off()


CairoPDF("psurv16.pdf",
         width = 7, height = 5, onefile = TRUE, bg = 'transparent')
p.surv.16
dev.off()


dev.off()


CairoPDF("psurv15.pdf",
         width = 7, height = 5, onefile = TRUE,bg = 'transparent')
p.surv.15
dev.off()




bayes.data.seeds$trt<-NA

bayes.data.seeds$trt<-c("Em.seed.add.14","Em.seed.control.14","Ep.seed.add.14","Ep.seed.control.14")

bayes.data.seeds$water<-as.factor(c(1,0,1,0))
bayes.data.seeds$endo<-as.factor(c(0,0,1,1))

ggplot(bayes.data.seeds, aes(x=trt, y=mean, colour= water)) + 
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.2) +
  geom_point(aes(shape = factor(endo)))+
  scale_color_manual(breaks = c("0", "1"),
                    values=c("red", "blue"))+
  labs( x = "Treatment",
        y = "Mean log(seeds produced)")

