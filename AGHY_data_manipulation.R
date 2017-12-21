## Authors: Marion and Tom
## Purpose: Create a script that imports AGHY data, performs all raw data manipulation,
## and creates an .RData object that can be loaded for analysis
## Last update: 5 July 2017
######################################################

options(java.parameters = "-Xmx1024m") ## clears the memory so we can use read.xlsx on these big files 

library(xlsx)
library(plyr)
library(dplyr)
library(lme4)
library(mvtnorm)
library(ggplot2)
library(tidyr)
invlogit<-function(x){exp(x)/(1+exp(x))}

setwd("C:/Users/Marion Donald/Dropbox/Rice/Projects/AGHY/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_Final")

##############################################################
####################################################################
## Read in AGHY plot info (don't need all of the extra info because it already knows you're in the project file)
AGHY.plots<-read.xlsx("AGHY_SFAEF_life_history_expt.xlsx",
                      sheetName="Plot-level data")
AGHY.plants<-read.xlsx("AGHY_SFAEF_life_history_expt.xlsx",
                       sheetName="Plant-level data original")
## Read in AGHY subplot plant counts
AGHY.subplots<-read.xlsx("AGHY_SFAEF_life_history_expt.xlsx",
                         sheetName="Subplot-level data")

## Read in AGHY agrinostics survey
AGHY.immunoblot<-read.xlsx("AGHY_SFAEF_life_history_expt.xlsx",
                           sheetName="Endophyte Survey")
## Read in AGHY recruit data
AGHY.recruits<-read.xlsx("AGHY_SFAEF_life_history_expt.xlsx",
                         sheetName="Plant-level data recruits")
## Read in AGHY 2015 seed score data 
AGHY.2015.seeds<-read.xlsx("AGHY_SFAEF_life_history_expt.xlsx", 
                           sheetName="2015 seed scores")

## Read in AGHY 2016 agrinostics data for strips
AGHY.agri.16<-read.xlsx("AGHY_agrinostics_2016.xlsx",
                        sheetName="AGHY_Endophyte Survey 2016")

## Read in AGHY 2016 agrinostics data for the demography individuals
AGHY.agri.demo.16<-read.xlsx("AGHY_agrinostics_2016.xlsx",
                             sheetName="Demography Individuals 2016")
####################################################################
####################################################################

## subset the AGHY.agri.demo.16 for just the AGHY plots
AGHY.agri.demo.16<-AGHY.agri.demo.16[AGHY.agri.demo.16$species == "AGHY",]

## drop farm plots
AGHY.plots<-AGHY.plots[AGHY.plots$transmission!="Farm",]
## drop weird NAs
AGHY.plots<-AGHY.plots[!is.na(AGHY.plots$plot),]
## drop a few problem plots...these were planted incorrectly but we continued data collection by mistake
AGHY.plots<-AGHY.plots[AGHY.plots$plot!=143 & AGHY.plots$plot!=211,]
AGHY.plots$newplot<-order(AGHY.plots$plot)

## match the agrinostic strip data to the year collected
AGHY.immunoblot<-AGHY.immunoblot[AGHY.immunoblot$year_t ==AGHY.immunoblot$strip,]
## check that we dropped the 2014 strips from 2015 data
which(AGHY.immunoblot$year_t != (AGHY.immunoblot$strip))

## Get the total E plus plants scored per plot
## Big change here. Last version did not track subplots
AGHY<-ddply(AGHY.immunoblot, c("year_t","plot","subplot"), plyr::summarize, 
            total = length(agri_liberal),
            E_plus_liberal = sum(agri_liberal),
            E_plus_conservative = sum(agri_conservative))
## there is one subplot labeled "extra". Drop this.
AGHY<-AGHY[-which(AGHY$subplot=="EXTRA"),]
## coerce subplot back to integer
AGHY$subplot<-as.integer(AGHY$subplot)

## Estimate frequency. Note these are by subplot.
AGHY$con_freq<-AGHY$E_plus_conservative/AGHY$total
AGHY$lib_freq<-AGHY$E_plus_liberal/AGHY$total

## Merge these two types of info
AGHY.merge<-merge(AGHY.plots,AGHY,by="plot")
AGHY.merge$lib_freq<-AGHY.merge$E_plus_liberal/AGHY.merge$total

## Now merge in subplot counts
## rename year for the merge
AGHY.subplots$year_t<-AGHY.subplots$year
AGHY.merge<-merge(AGHY.merge,AGHY.subplots,by=c("plot","subplot","year_t", "species"))

## select the relevant columns
AGHY.new<-AGHY.merge[, c("plot","water","target_init_freq",
                         "year_t","total","con_freq","lib_freq")]
## copy this dataframe into a new one so that year_t1 and the frequencies can be labeled
AGHY.freq.1<-AGHY.new
## create the year_t1 column from the year_t column
AGHY.freq.1$year_t1<-AGHY.freq.1$year_t+1
## rename the year_t frequencies
names(AGHY.freq.1)[names(AGHY.freq.1) == "con_freq"]<- "con_freq_t"
names(AGHY.freq.1)[names(AGHY.freq.1) == "lib_freq"]<- "lib_freq_t"
names(AGHY.freq.1)[names(AGHY.freq.1)=="total"]<-"total_scored_t"
## assign the year t to the year t1 to match with the AGHY.freq.1 dataframe (and do the same with total seeds scored)
AGHY.new$year_t1<-AGHY.new$year_t
names(AGHY.new)[names(AGHY.new) == "total"]<-"total_scored_t1"
## rename the year_t1 frequencies
names(AGHY.new)[names(AGHY.new) == "con_freq"]<- "con_freq_t1"
names(AGHY.new)[names(AGHY.new) == "lib_freq"]<- "lib_freq_t1"


## New data frame with years t and t+1 and their frequencies 
AGHY.total<-merge(AGHY.freq.1, AGHY.new[,c("plot","total_scored_t1","con_freq_t1","lib_freq_t1","year_t1")], by= c("plot", "year_t1"))
## re-organizing the columns so they make sense visually 
AGHY.total<-AGHY.total[,c(1,3:8,2,9:11)]


## Add in the recruit data to the bigdataframe in AGHY.merge
AGHY.merge.recruits<-merge(AGHY.merge, AGHY.recruits)

## If dead in winter 2014 (winter_t) apply the 0 to survival in spring 2015 and same for dead in winter 2015, apply to spring 2016 (spring_t1)
AGHY.merge.recruits <- within(AGHY.merge.recruits, {spring_survival_t1 = ifelse(spring_survival_t1 == 0 | winter_survival_t1 == 0,0, 1)} )

### Need endo scores for the demography plants in 2015 (this is in the separate sheet right now, so we'll need to merge it in correctly)
## 2015, need to drop the individuals for which seeds were attempted (so they're in the dataframe) but 0 seeds were actually scored
AGHY.2015.seeds = AGHY.2015.seeds[AGHY.2015.seeds$complete_seeds_scored != 0,]

AGHY.2015.seeds$Vtrans_t<- AGHY.2015.seeds$total_E_plus/AGHY.2015.seeds$complete_seeds_scored ### IS THIS THE CORRECT VERTICAL TRANSMISSION FOR 2015?? -- this is from summation in the excel file...
AGHY.2015.seeds$endo.stat_2015<-AGHY.2015.seeds$total_E_plus>0
#hist(AGHY.2015.seeds$Vtrans_t[AGHY.2015.seeds$endo.stat_2015==1])

## subset plants that were assigned E+ but have low transmission (same deal as 2014)
#weird2015<-subset(AGHY.2015.seeds,endo.stat_2015==1 & Vtrans_t<0.5)
#weirdness<-weird2015[,c("complete_seeds_scored","total_E_plus","total_seeds_scored_1","e_plus_1",
#                        "total_seeds_scored_2","e_plus_2","total_seeds_scored_ad", 
#                        "e_plus_ad","total_seeds.scored_3","e_plus_3")]

##OK we are definitely assigning as E- anything that had a single E+ seed and in the first round
## If the only E+ came from one seed in round 1, we are calling it E-
AGHY.2015.seeds <- within(AGHY.2015.seeds, {endo.stat_2015 = ifelse(total_E_plus == 1 & total_E_minus > 0 & e_plus_1 == 1, 0, endo.stat_2015)})

## create summation of seeds scored in additional rounds
AGHY.2015.seeds$total.add.seeds<-ifelse(is.na(AGHY.2015.seeds$total_seeds_scored_2),0,AGHY.2015.seeds$total_seeds_scored_2)+
  ifelse(is.na(AGHY.2015.seeds$total_seeds_scored_ad),0,AGHY.2015.seeds$total_seeds_scored_ad)+
  ifelse(is.na(AGHY.2015.seeds$total_seeds.scored_3),0,AGHY.2015.seeds$total_seeds.scored_3)

## If all the E+ came from round 1 AND all subsequent scoring (>5 seeds total) was all E-, we call it E-
AGHY.2015.seeds <- within(AGHY.2015.seeds, {endo.stat_2015 = ifelse(total_E_plus == e_plus_1 & total.add.seeds > 5,0,endo.stat_2015)})
#hist(AGHY.2015.seeds$Vtrans_t[AGHY.2015.seeds$endo.stat_2015==1])
#AGHY.2015.seeds$year_t1 = AGHY.2015.seeds$year_t



#names(AGHY.2015.seeds)[names(AGHY.2015.seeds) == "year_t"]<- "year_t1"



AGHY.demo<-merge(AGHY.merge.recruits, AGHY.2015.seeds[,c("plot","subplot","ID","year_t","endo.stat_2015")], by=c("plot","subplot","ID", "year_t"), all=T)

## change the column headers in AGHY.agri.demo.16 to match the agrinostic headers in AGHY.demo
names(AGHY.agri.demo.16)[names(AGHY.agri.demo.16) == "agri_liberal"]<- "tiller_agri_lib_t1"
names(AGHY.agri.demo.16)[names(AGHY.agri.demo.16) == "agri_conservative"]<- "tiller_agri_con_t1"

## rename the plant id column to be ID
names(AGHY.agri.demo.16)[names(AGHY.agri.demo.16) == "plant"]<-"ID"

## copy the 2016 year_t to a column for year_t1 -- so we can merge it into the t1 matched with t = 2015 and as t with t=2016
AGHY.agri.demo.16$year_t1<-AGHY.agri.demo.16$year_t
AGHY.agri.demo.16$endo.stat.2016<-with(AGHY.agri.demo.16, ifelse(tiller_agri_lib_t1 == 1,1,0))

## Finally, write an RData file, coming later on


## merge the 2016 agrinostics scores for recruit demography plants with the main data frame to match with 2015 as year t
AGHY.demo.final<-merge(AGHY.demo, AGHY.agri.demo.16[,c("plot","subplot","ID","year_t1", "tiller_agri_lib_t1", "tiller_agri_con_t1", "endo.stat.2016")],
                       by=c("plot","subplot","ID", "year_t1","tiller_agri_lib_t1", "tiller_agri_con_t1"), all=T)



## Tom addition 6/12 -- we currently have recruit endo status ('endo.stat') based on seeds (in 2015) or tiller (in 2016)
## I would like to add this for 2014
AGHY.demo.final$endo.stat.2014<-NA
AGHY.demo.final$endo.stat.2014[AGHY.demo.final$year_t==2014 & !is.na(AGHY.demo.final$tot_seeds_scored_t)]<-ifelse(
  AGHY.demo.final$e_p_score_t[AGHY.demo.final$year_t==2014 & !is.na(AGHY.demo.final$tot_seeds_scored_t)]>0,1,0)

## Add in unique ID for each of the plants based on year_t1
AGHY.demo.final$unique.id = paste(AGHY.demo.final$plot, AGHY.demo.final$subplot, 
                                  AGHY.demo.final$ID, AGHY.demo.final$year_t1, sep ="_")

## Marion addition 6/30 -- "AGHY Challenge" - Individual-level enophyte status, applied forward and backward for each individual for which we have a score
## Need to keep in mind that by extending it backward and forward, some of these plants may not be alive in that year
## -- think it'll be fine, bc we later sort by survival prior to doing the following analyses 

## If all the E+ came from round 1 AND all subsequent scoring (>5 seeds total) was all E-, we call it E-
AGHY.2015.seeds.endo <- within(AGHY.2015.seeds, {endo.stat_2015 = ifelse(total_E_plus == e_plus_1 & total.add.seeds > 5,0,endo.stat_2015)})
#hist(AGHY.2015.seeds$Vtrans_t[AGHY.2015.seeds$endo.stat_2015==1])

names(AGHY.2015.seeds.endo)[names(AGHY.2015.seeds.endo) == "year_t"]<- "year_t1"



AGHY.demo.endo<-merge(AGHY.merge.recruits, AGHY.2015.seeds.endo[,c("plot","subplot","ID","year_t1","endo.stat_2015")], by=c("plot","subplot","ID", "year_t1"), all=T)

## change the column headers in AGHY.agri.demo.16 to match the agrinostic headers in AGHY.demo
names(AGHY.agri.demo.16)[names(AGHY.agri.demo.16) == "agri_liberal"]<- "tiller_agri_lib_t1"
names(AGHY.agri.demo.16)[names(AGHY.agri.demo.16) == "agri_conservative"]<- "tiller_agri_con_t1"

## rename the plant id column to be ID
names(AGHY.agri.demo.16)[names(AGHY.agri.demo.16) == "plant"]<-"ID"

## copy the 2016 year_t to a column for year_t1 -- so we can merge it into the t1 matched with t = 2015 and as t with t=2016
AGHY.agri.demo.16$year_t1<-AGHY.agri.demo.16$year_t
AGHY.agri.demo.16$endo.stat.2016<-with(AGHY.agri.demo.16, ifelse(tiller_agri_lib_t1 == 1,1,0))

## Finally, write an RData file, coming later on


## merge the 2016 agrinostics scores for recruit demography plants with the main data frame to match with 2015 as year t
AGHY.demo.final.endo<-merge(AGHY.demo.endo, AGHY.agri.demo.16[,c("plot","subplot","ID","year_t1", "tiller_agri_lib_t1", "tiller_agri_con_t1", "endo.stat.2016")],
                       by=c("plot","subplot","ID", "year_t1","tiller_agri_lib_t1", "tiller_agri_con_t1"), all=T)

## Tom addition 6/12 -- we currently have recruit endo status ('endo.stat') based on seeds (in 2015) or tiller (in 2016)
## I would like to add this for 2014
AGHY.demo.final.endo$endo.stat.2014<-NA
AGHY.demo.final.endo$endo.stat.2014[AGHY.demo.final.endo$year_t==2014 & !is.na(AGHY.demo.final.endo$tot_seeds_scored_t)]<-ifelse(
  AGHY.demo.final.endo$e_p_score_t[AGHY.demo.final.endo$year_t==2014 & !is.na(AGHY.demo.final.endo$tot_seeds_scored_t)]>0,1,0)

## Marion addition 6/30 -- "AGHY Challenge" - Individual-level enophyte status, applied forward and backward for each individual for which we have a score
## Need to keep in mind that by extending it backward and forward, some of these plants may not be alive in that year
## -- think it'll be fine, bc we later sort by survival prior to doing the following analyses 

## select out just the three columns that we're interested in: endo status in 2014, 2015, and 2016
endo.stat <- AGHY.demo.final.endo[,c(57:59)]

## remove the rows that don't have any endo scores -- they're all NA's, will have to use the plot probability to estimate endo status for these individuals
endo.stat.sub = AGHY.demo.final.endo[!is.na(AGHY.demo.final.endo$endo.stat.2014) | !is.na(AGHY.demo.final.endo$endo.stat_2015) | !is.na(AGHY.demo.final.endo$endo.stat.2016), ]

## get the rows that don't have any endo scores (just NA's)
## NEED TO MERGE these back in with the ones that do have scores later 
endo.stat.na = AGHY.demo.final.endo[is.na(AGHY.demo.final.endo$endo.stat.2014) & is.na(AGHY.demo.final.endo$endo.stat_2015) & is.na(AGHY.demo.final.endo$endo.stat.2016), ]

## identifying if there's a match between '14-'15, '14-'16, and '15-'16, this spits out Trues and Falses
endo.stat.sub.whichidentical = endo.stat.sub %>% 
  rowwise() %>% 
  summarise(same1415 = identical(endo.stat.2014,endo.stat_2015), same1416 = identical(endo.stat.2014,endo.stat.2016),same1516 = identical(endo.stat.2016,endo.stat_2015))

## turn the Trues and Falses into 1 and 0 and into the three columns 
endo.stat.sub.identical = endo.stat.sub[rowSums(matrix(as.numeric(as.matrix(endo.stat.sub.whichidentical)), ncol = 3)) >= 1,]

## pull out the rows that have disagreement in the endo scores (there are 16 of them -- disagreement is between 2014 and 2015, no scores entered for 2016)
## NEED TO TROUBLESHOOT THESE, and ALSO PUT THEM BACK IN THE MAIN DATAFRAME
endo.stat.sub.nonidentical = endo.stat.sub[rowSums(matrix(as.numeric(as.matrix(endo.stat.sub.whichidentical)), ncol = 3)) == 0,]

## make 3 new endo status columns to copy the endo status forward and backward and to keep the original data in the dataframe 
endo.stat.sub.identical$endo.stat.14<-endo.stat.sub.identical$endo.stat.2014
endo.stat.sub.identical$endo.stat.15<-endo.stat.sub.identical$endo.stat_2015
endo.stat.sub.identical$endo.stat.16<-endo.stat.sub.identical$endo.stat.2016

## copying the endo status and applying it across all of the years 
endo.stat.sub.identical$endo.stat.15[is.na(endo.stat.sub.identical$endo.stat.15)] = endo.stat.sub.identical$endo.stat.16[is.na(endo.stat.sub.identical$endo.stat.15)]
endo.stat.sub.identical$endo.stat.15[is.na(endo.stat.sub.identical$endo.stat.15)] = endo.stat.sub.identical$endo.stat.14[is.na(endo.stat.sub.identical$endo.stat.15)]
endo.stat.sub.identical$endo.stat.16 = endo.stat.sub.identical$endo.stat.15
endo.stat.sub.identical$endo.stat.14 = endo.stat.sub.identical$endo.stat.15

## Create a column for "endo.stat" -- which is now independent of year
endo.stat.sub.identical$endo.stat = endo.stat.sub.identical$endo.stat.14

## drop the excess columns, so that the dataframe just has the original data + the column with endo.stat not tied to year
drops <- c("endo.stat.16","endo.stat.15", "endo.stat.14")
AGHY.demo.semi.final.endo = endo.stat.sub.identical[ , !(names(endo.stat.sub.identical) %in% drops)]

## need to add the data we excluded earlier on back in: all of the rows that didn't have an endo score + the ones that disagree about the score
## first need to add in the endo score to the other two data frames that we pulled apart from the main one earlier
endo.stat.na$endo.stat = NA
endo.stat.sub.nonidentical$endo.stat = NA


#### THIS RESOLUTION NEEDS TO OCCUR EARLIER AND IN THE MAIN DATA FRAME

endo.stat.sub.nonidentical$unique.id = NA

## create a unique ID for each plant based on year_t -- this will match this disagreement dataframe with the demography data frame for easy manual manipulation of the endo scores 
endo.stat.sub.nonidentical$unique.id = paste(endo.stat.sub.nonidentical$plot, endo.stat.sub.nonidentical$subplot, endo.stat.sub.nonidentical$ID, endo.stat.sub.nonidentical$year_t1, sep ="_")

## change the values of the disagreement in endo score between 2014 and 2015, remove the 6 that we couldn't resolve. See reasons for change below

# 119_3_2_2015 - change 2014 to neg, 4 seeds scored neg in 2015 by Charlene, 2 of 18 in 2014 scored positive. Likely the 2014 score is incorrect
AGHY.demo.final <- within(AGHY.demo.final, endo.stat.2014[unique.id == "119_3_2_2015"] <- 0 ) 

# 120_2_2_NA - REMOVE, 4 seeds scored in 2015 by Kory, but this plant was dead in 2015, so the seeds must have been mislabeled. REMOVE from dataframe
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "120_2_2_NA",]

# 140_2_1_2015 - change 2015 to neg, only 2 seeds able to be scored in 2015, 1 E+ and 1 E- by Kory, take the negative score from 2014
AGHY.demo.final <- within(AGHY.demo.final, endo.stat_2015[unique.id == "140_2_1_2015"] <- 0 ) 

# 140_3_4_2015 - Change 2015 to neg. Scored by Nicole and Kory, Nicole found 3 positive seeds 2 negative seeds, Kory found 11 negative seeds -- likely that the positives could have been mis-scored
AGHY.demo.final <- within(AGHY.demo.final, endo.stat_2015[unique.id == "140_3_4_2015"] <- 0 )

# 146_1_1_2015 - change 2015 to neg 4 seeds scored by Phillip called all positive, likely these were mis-scored, take the negative score from 2014
AGHY.demo.final <- within(AGHY.demo.final, endo.stat_2015[unique.id == "146_1_1_2015"] <- 0 )

# 174_1_3_2015 - change 2014 to pos. 17 out of 20 seeds E+ in 2015, scored by Nicole and Charlene
AGHY.demo.final <- within(AGHY.demo.final, endo.stat.2014[unique.id == "174_1_3_2015"] <- 1 ) 

# 188_2_1_2015 - change 2015 to neg. 1 seed scored in 2015 by Kory (29 attempted)
AGHY.demo.final <- within(AGHY.demo.final, endo.stat_2015[unique.id == "188_2_1_2015"] <- 0 )

# 188_3_1_2015 - change 2015 to neg, 4 seeds scored E+ in the first round by Nicole, then 5 seeds scored in the second round by Nicole as E-
AGHY.demo.final <- within(AGHY.demo.final, endo.stat_2015[unique.id == "188_3_1_2015"] <- 0 )

## drop the rows that are unresolved for endo status between 2014 and 2015 (208_1_4_2015, 208_1_4_2016, 208_1_3_2015, 208_1_3_2016, 204_2_3_2015, 204_2_3_2016, 164_2_1_2015, 164_2_1_2016)
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "208_1_4_2015",] 
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "208_1_4_2016",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "208_1_3_2015",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "208_1_3_2016",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "204_2_3_2015",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "204_2_3_2016",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "164_2_1_2015",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "164_2_1_2016",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "163_1_4_2015",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "163_1_4_2016",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "147_1_1_2015",]
AGHY.demo.final <- AGHY.demo.final[AGHY.demo.final$unique.id != "147_1_1_2016",]



## select out just the three columns that we're interested in: endo status in 2014, 2015, and 2016
endo.stat <- AGHY.demo.final[,c(57:59)]

## remove the rows that don't have any endo scores -- they're all NA's, will have to use the plot probability to estimate endo status for these individuals
endo.stat.sub = AGHY.demo.final[!is.na(AGHY.demo.final$endo.stat.2014) | !is.na(AGHY.demo.final$endo.stat_2015) | !is.na(AGHY.demo.final$endo.stat.2016), ]

## get the rows that don't have any endo scores (just NA's)
## NEED TO MERGE these back in with the ones that do have scores later 
endo.stat.na = AGHY.demo.final[is.na(AGHY.demo.final$endo.stat.2014) & is.na(AGHY.demo.final$endo.stat_2015) & is.na(AGHY.demo.final$endo.stat.2016), ]

## identifying if there's a match between '14-'15, '14-'16, and '15-'16, this spits out Trues and Falses
endo.stat.sub.whichidentical = endo.stat.sub %>% 
  rowwise() %>% 
  summarise(same1415 = identical(endo.stat.2014,endo.stat_2015), same1416 = identical(endo.stat.2014,endo.stat.2016),same1516 = identical(endo.stat.2016,endo.stat_2015))

## turn the Trues and Falses into 1 and 0 and into the three columns 
endo.stat.sub.identical = endo.stat.sub[rowSums(matrix(as.numeric(as.matrix(endo.stat.sub.whichidentical)), ncol = 3)) >= 1,]

## pull out the rows that have disagreement in the endo scores (there are 16 of them -- disagreement is between 2014 and 2015, no scores entered for 2016)
## NEED TO TROUBLESHOOT THESE, and ALSO PUT THEM BACK IN THE MAIN DATAFRAME
endo.stat.sub.nonidentical = endo.stat.sub[rowSums(matrix(as.numeric(as.matrix(endo.stat.sub.whichidentical)), ncol = 3)) == 0,] ## GREAT! This is now empty!

## make 3 new endo status columns to copy the endo status forward and backward and to keep the original data in the dataframe 
endo.stat.sub.identical$endo.stat.14<-endo.stat.sub.identical$endo.stat.2014
endo.stat.sub.identical$endo.stat.15<-endo.stat.sub.identical$endo.stat_2015
endo.stat.sub.identical$endo.stat.16<-endo.stat.sub.identical$endo.stat.2016

## copying the endo status and applying it across all of the years 
endo.stat.sub.identical$endo.stat.15[is.na(endo.stat.sub.identical$endo.stat.15)] = endo.stat.sub.identical$endo.stat.16[is.na(endo.stat.sub.identical$endo.stat.15)]
endo.stat.sub.identical$endo.stat.15[is.na(endo.stat.sub.identical$endo.stat.15)] = endo.stat.sub.identical$endo.stat.14[is.na(endo.stat.sub.identical$endo.stat.15)]
endo.stat.sub.identical$endo.stat.16 = endo.stat.sub.identical$endo.stat.15
endo.stat.sub.identical$endo.stat.14 = endo.stat.sub.identical$endo.stat.15

## Create a column for "endo.stat" -- which is now independent of year
endo.stat.sub.identical$endo.stat = endo.stat.sub.identical$endo.stat.14

## drop the excess columns, so that the dataframe just has the original data + the column with endo.stat not tied to year
drops <- c("endo.stat.16","endo.stat.15", "endo.stat.14")
AGHY.demo.semi.final = endo.stat.sub.identical[ , !(names(endo.stat.sub.identical) %in% drops)]

## need to add the data we excluded earlier on back in: all of the rows that didn't have an endo score + the ones that disagree about the score
## first need to add in the endo score to the other two data frames that we pulled apart from the main one earlier
endo.stat.na$endo.stat = NA

## this step is ultimately unnecessary, but while still fiddling with this manipulation, I think it's useful to have
## data frame is out of numerical plot order, but that shouldn't matter
AGHY.demo.final.final = rbind(AGHY.demo.semi.final,endo.stat.na) ## rbind the three dataframes back together, check that it's the same length as the original

## rename the AGHY.demo.final.final to AGHY.demo.final so that we don't have to adjust the later code that calls on AGHY.demo.final
AGHY.demo.final = AGHY.demo.final.final



#################################################################################################################
######## Seed Production Data Manipulation ######################################################################
######## Want: grand average of # seeds produced per plant by E+ and E- plants in water control and water addition plots
######## Action: Multipy these means by the # of E+ and E- plants in each plot to get the total amount of seeds produced in each plot

## need to take the seed_mass_t (from 2013) sum it by endo stat and water add vs control and divide by the total number of infs that contributed
### this gives the seed mass per inf 
## need to add in the water control and addition info to the AGHY.plants dataframe
AGHY.plants.water<-merge(AGHY.plots, AGHY.plants, all =T)
AGHY.plants.water<-subset(AGHY.plants.water, newplot >0 )

## remove the plants that we know are dead (there are 5 of the demo plants)
## **** MIGHT BE AN ISSUE THAT WE ONLY KNOW SURVIVAL STATUS FOR THE ORIG PLANTS THAT WERE DEMO INDIVIDUALS (HALF THE PLANTS.....)
## THIS ANALYSIS PROCEEDS AS IF ALL OF THE NON-DEMO PLANTS THAT WE DON'T KNOW SURVIVAL FOR ARE ALIVE.... MAY NEED TO CHANGE THIS 
dead.orig<- which(with(AGHY.plants.water, survival_t == 0)) ## indentify the dead original plants
AGHY.plants.water <- AGHY.plants.water[-dead.orig, ] ## remove them from the data frame 

seed.mass<-subset(AGHY.plants.water, seed_mass_t >0 & demo_plant ==1) # select out the plants that have seed masses

## calculate the average seed mass produced per plant by E+ and E- plants in Water Control and Water Addition plots 
seed.mass.avg<-ddply(seed.mass, c("endo", "water"), plyr:: summarize, 
                     num.infs = sum(inf_collected_t), # total infs collected for E+ and E- plants in water Control and Water Addition 
                     seed.mass = sum(seed_mass_t), # total seed mass for E+ and E- plants in Water Control and Water Addition
                     mass.inf = seed.mass/num.infs, # avg seed mass per inf for E+ and E- plants in Water Control and Water Addition
                     num.infs.plant = sum(inf_number_t), # total number of infs per plant for E+ and E- plants in Water Control and Water Addition
                     avg.mass.plant = mass.inf * num.infs.plant) # average seed mass per plant for E+ and E- plants in Water Control and Water Addition

## get the total number of E+ and E- plants in each plot

plants.endo<- ddply(AGHY.plants.water, c("endo", "water", "newplot"), plyr::summarize,
                    num.plants = length(ID))

seed.mass.plant<-merge(plants.endo, seed.mass.avg, all = T)
seed.mass.plant$seeds.plot<-seed.mass.plant$num.plants * seed.mass.plant$avg.mass.plant

##trash
hist(log(AGHY.plants.water$seed_mass_t[AGHY.plants.water$year_t==2013]))

## then want to scale this number (specific for E+ and E- plants in each plot) back down to the subplot level so that we can use it in the recruitment analysis with plot as a random effect
## gives contribution of seed mass by E+ and E- plants at the subplot level -- assumes equal seed rain into each of the subplots from the plants in the plot 
seed.mass.plant$seed.mass.sbplt<- seed.mass.plant$seeds.plot / 64 # 64 subplots per plot

## cut out the newplot and subplot data from original data frame and merge this seed mass subplot back in, so that we have the right length dataframe to loop over for the 2014 recruitment analysis


subplots.2013<-ddply(AGHY.demo.final, c("newplot","subplot","year"), plyr :: summarize,
                     plants = length(ID))

subplots.2013<-subset(subplots.2013, year == 2014)

orig.2013.seeds.sbplts<- merge(subplots.2013, seed.mass.plant) ## dataframe with plot and subplot attached to seed mass per subplot + endo status of parent plant contributing the seeds 

## for the linear regression dataframe, change the water status from 0 and 1 to control and add to match with the others 

AGHY.seed.mass.count$water<-AGHY.seed.mass.count$Status

AGHY.seed.mass.count$water[AGHY.seed.mass.count$water ==1]<- "Add"
AGHY.seed.mass.count$water[AGHY.seed.mass.count$water ==0]<- "Control"
AGHY.seed.mass.count$water<-as.factor(AGHY.seed.mass.count$water)



## Get the numbers of E+ and E- original plants that were alive in 2014 -- able to contibute to seed production in the plots
orig.plants.alive.14 <- subset( AGHY.plants.water, survival_t1 == 1) ## first get all the original plants that were still alive in 2014

## get the number of E+ and E- orig plants still alive in each plot 
orig.plants.endo.14 <- ddply(orig.plants.alive.14, c("endo", "water", "newplot"), plyr::summarize,
                             num.plants = length(ID))

## combine seed mass specific for endo and water to each of these plots in 2014 -- ORIGINAL PLANTS

seed.mass.plant.14<- merge(orig.plants.endo.14, seed.mass.avg)
seed.mass.plant.14$seeds.plot.14 <- seed.mass.plant.14$num.plants * seed.mass.plant.14$avg.mass.plant ## mass of seeds contributed by E+ and E- plants to each plot


## Seed production in 2014 and 2015
## number of flowering recruits in each subplot -- want this at the subplot level so that we can include a plot random effect
## first remove subplots that we don't have any data from -- sometimes we skipped if there were fire ants covering them
spring.fl.rec<- AGHY.demo.final[!is.na(AGHY.demo.final$spring_flowering_recruit_count_per_subplot),]

spring.flowering.recruits.sbplts<-ddply(spring.fl.rec, c("newplot", "subplot", "year_t"), plyr :: summarize,
                                        total.flowering.sbplt = sum(spring_flowering_recruit_count_per_subplot))



## avg number of seeds produced by E+ and E- recruits in each subplot (separate for water addition and water control)
seed.mass.recruits.demo<-subset(AGHY.demo.final,  endo.stat == 1 | endo.stat == 0)
seed.mass.recruits.demo<-subset(seed.mass.recruits.demo, seed_mass_t > 0)

## get the mass of E+ and E- seeds for water A and C plots from the demography recruits in the 3 subplots 
seed.mass.recruits.plot<-ddply(seed.mass.recruits.demo, c("newplot", "subplot","year_t", "water", "endo.stat"), plyr :: summarize,
                               num.infs.coll = sum(inf_collected_t), # number of infs collected for E+ and E- recruits in water addition and water control 
                               seed.mass = sum(seed_mass_t), # total seed mass for E+ and E- recruits in water A and C plots
                               mass.inf = seed.mass/num.infs.coll,
                               num.infs = sum(inf_number_t), # number of infs per plot (based on the 3 subplots
                               num.plants = length(ID),
                               avg.mass.plant = mass.inf * num.infs/num.plants)

## get the total number of demo plants in each of the three subplots
demo.sbplt.total<-ddply(seed.mass.recruits.plot, c("newplot","subplot", "year_t"), plyr :: summarize,
                        total.plants = sum(num.plants))

smf<-merge(seed.mass.recruits.plot, spring.flowering.recruits.sbplts) 

seed.mass.flowering<-merge(smf, demo.sbplt.total)

seed.mass.flowering$prop.total<- seed.mass.flowering$num.plants/seed.mass.flowering$total.plants ## get the proportion of E+ and E- demo plants

seed.mass.flowering$endo.plant.sbplt<- seed.mass.flowering$prop.total * seed.mass.flowering$total.flowering.sbplt ## get the total number of E+ and E- plants in each subplot

seed.mass.flowering$seed.mass.endo.sbplt<- seed.mass.flowering$endo.plant.sbplt * seed.mass.flowering$avg.mass.plant ## the avg. mass of endo specific seeds per subplot
## stopping here because we want to do the anaysis at the subplot level 
#seed.mass.flowering$avg.endo.seed.mass.plot <- seed.mass.flowering$avg.seed.mass.endo.sbplt * 64 ## the estimated number of endo specific seeds from recruits per plot


######################################################################################################################
######################################################################################################################
## 2014 RECRUITMENT DATA


AGHY.merge$total_spring_recruits_subplot <- AGHY.merge$spring_vegetative_recruit_count_per_subplot + AGHY.merge$spring_flowering_recruit_count_per_subplot


AGHY.merge$water<-as.factor(AGHY.merge$water)

AGHY.merge$water<-as.integer(AGHY.merge$water) ## change to 1 and 2 for indexing water in the bayes model 

## only need to do this once, but re-run if data manipulation gets updated (7/26/17)
#save.image("C:/Users/tm9/Dropbox/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_data_clean1.RData") 
save.image("C:/Users/Marion Donald/Dropbox/Rice/Projects/AGHY/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_Final/AGHY_data_clean.1.RData") ## Marion's computer

alive.14.orig<-ddply(seed.mass.plant.14, "newplot", plyr:: summarize,
                     total_plants = sum(num.plants),
                     prop.alive = total_plants/20,
                     prop.dead = 1-prop.alive)

alive.14.orig.40<-subset(alive.14.orig, prop.alive > .4)


## only need to do this once, but re-run if data manipulation gets updated (6/30/17)
#save.image("C:/Users/tm9/Dropbox/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_data_clean.RData") 
#save.image("C:/Users/Marion Donald/Dropbox/Rice/Projects/AGHY/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_data_clean.1.RData") ## Marion's computer

## (re-ran - but didn't change the script 12/18/17)
#save.image("C:/Users/Marion Donald/Dropbox/Rice/Projects/AGHY/AGHY_SFAEF_Project/AGHY analysis summer2017/AGHY_Final/AGHY_data_clean.Dec18.RData")
