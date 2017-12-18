load("AGHY_Bayes.RData")
library(ggplot2)

## endo prev at the plot level for all years + water treatments
endo.all.plot

### need to resave the AGHY_Bayes_fit_change_rec_survival script with these new figures -- but for now running and calling them here

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


### survival probability figures with facet_wrap


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


flowering.14<-subset(flowering.all, year == 2014)
flowering.15<-subset(flowering.all, year == 2015)
flowering.16<-subset(flowering.all, year == 2016)

#### Plot the mean probability of flowering for E+ and E- plants by water and control treatments with the 2.5% and 97.5% Confidence Intervals
p.flow.14<-ggplot(flowering.14, aes(treatment, mean, colour=water.treat))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  geom_point(size=3, stat="identity")+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_x_discrete("Symbiont Presence", labels = c("Em.flower.add.14" = "E-","Em.flower.control.14" = "E-",
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

