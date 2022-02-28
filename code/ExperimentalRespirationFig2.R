
# This code generate figure 2 using respiration data WIth  3 panels sps, com M9 and com SM. 
# Boxplot are included within the respiration curves plot.
# This includes the analysis for species using nlme. 
# This code also contains analyses included in Extended data Table 4&5


library(data.table)
library(plyr)
library(tools)
library(zoo)
library(lattice)
library(ggplot2)
require(minpack.lm)
library(reshape2)
library(nlsLoop)
library(nlstools)
library(gridExtra)

library(ggplot2)
require(nlme)
library(TeamPhytoplankton)
library(AICcmodavg)

## library(lme4)
resl_grp <-  read.csv(file = "./data/resl_grp_taxa.csv" , head = TRUE,stringsAsFactors = T)
head(resl_grp)

# make grouped data object
sp_grp <- groupedData(logr ~ K | species, resl_grp)
plot(sp_grp) # note I have used species as random effect because this is the level of each curve - i.e. nested random effect of taxa/treatment

##Full model
nlme.full <- nlme(logr ~ schoolfield.high(ln.c,Ea,Eh,Th,temp=K,Tc=18),
               data = sp_grp,
               fixed = list(ln.c + Ea + Eh + Th ~ 1 + Treatment),
               random = pdDiag(ln.c + Ea+ Eh + Th ~ 1),
                #groups = ~ sps,
               start = c(-0.3, 0, 0.8, 0, 8, 0, 306, 0),
               na.action = na.omit,
               method = "ML",
               control = nlmeControl(pnlsTol = 0.2))
summary(nlme.full)

# remove treatment on effect on ln.c 
resl.mix1<- update(nlme.full, fixed = list(ln.c ~ 1, Ea + Eh + Th ~ 1 + Treatment), start = c(-0.3, 0.8, 0, 6, 0, 306, 0))
anova(nlme.full, resl.mix1) # no treatment effect on ln.c
summary(resl.mix1)

# remove treatment on effect on Ea
resl.mix2<- update(resl.mix1, fixed = list(ln.c  + Ea ~ 1, Eh + Th ~ 1 + Treatment), start = c(-0.3, 0.8, 6, 0, 306, 0))
anova(resl.mix1, resl.mix2) # no treatment effect on Ea
summary(resl.mix2)

# remove treatment on effect on Eh
resl.mix3<- update(resl.mix2, fixed = list(ln.c  + Ea + Eh ~ 1, Th ~ 1 + Treatment), start = c(-0.3, 0.8, 6, 306, 0))
anova(resl.mix2, resl.mix3) # no treatment effect on Eh
summary(resl.mix3)

# remove treatment on effect on Th
resl.mix4<- update(resl.mix3, fixed = list(ln.c  + Ea + Eh + Th ~ 1), start = c(-0.3, 0.8, 8, 306))
anova(resl.mix3, resl.mix4) # no treatment effect on Th
summary(resl.mix4)

cand.set <- list(nlme.full, resl.mix1, resl.mix2,resl.mix3,resl.mix4)
names(cand.set)<-c('nlme.full',  'resl.mix1', 'resl.mix2', 'resl.mix3','resl.mix4')

tab <- aictab(cand.set)
tab 

final_mod <- nlme(logr ~ schoolfield.high(ln.c,Ea,Eh,Th,temp=K,Tc=18),
               data = sp_grp,
               fixed = list(ln.c + Ea + Eh + Th ~ 1),
               random = pdDiag(ln.c + Ea+ Eh + Th ~ 1),
                #groups = ~ sps,
               start = c(-0.3, 0.8, 8, 306),
               na.action = na.omit,
               method = "REML",
               control = nlmeControl(pnlsTol = 0.2))

summary(final_mod)

confint.mm <- data.frame(intervals(final_mod, which='fixed')[1])
#for the rest, we actually want a loop because it does get lengthy otherwise 
for(i in 1:nrow(confint.mm)){
  confint.mm$stdev[i] <- confint.mm$fixed.upper[i]-confint.mm$fixed.est[i]
}

confint.mm

################# PLOTS #################
## temperature dependence plots with fitted curves
# generate predicted curves from fixed effects for NP
predK <- seq(min(sp_grp$K), max(sp_grp$K), by = 0.1)
predTemp <- predK - 273.15
predmu_a <- schoolfield.high(ln.c = fixef(final_mod)[1], Ea = fixef(final_mod)[2], Eh = fixef(final_mod)[3], Th = fixef(final_mod)[4], temp = predK, Tc = 18)
predmu <- data.frame(T = predTemp, logr = c(predmu_a))

## Topt estimated
## Ea for each replicate in random effect
Coef <- coef(final_mod)
Id <- row.names(Coef)
aa <- strsplit(Id," ")
Coef$Sps <- c(sapply(aa, "[",1))
Coef$Treatment<- c(sapply(aa, "[",2))

write.csv( Coef, "./data/Resp_E_ind.csv")

#pdf('growth_rate_tdep_parms.pdf', 10, 7)
##### temperature dependence plot
sp_grp$T <- sp_grp$K-273.15
predmu$Treatment = rep("Nonevolved", nrow(predmu))
## plot using boxplot
Mean_B_Ea<-ggplot(Coef, aes(x = Treatment, y = Ea, fill = Treatment,color=Treatment))+
  geom_boxplot(outlier.color = "gray",outlier.size = 0.5)+
  scale_fill_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo'))+
  scale_colour_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo'))+
    ylim(0.3,2)+
    theme(axis.text.y=element_blank())+
    scale_x_discrete(labels = c('co-evolved', 'de novo')) +
  ylab(expression(E[a](eV)))+
    theme(legend.position = "none",axis.title.y=element_blank())+
        theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10)) +
      theme(legend.position = "none")+
    stat_summary(geom = "crossbar",color="white",fatten=1.6,width=0.75,
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})

c <- Mean_B_Ea

PANEL_A<- ggplot(sp_grp, aes(x = T, y = logr, colour = Treatment)) +
          geom_point(size = 2,alpha=0.5) +
    geom_line(aes(x = T, y = logr, colour = Treatment), data = predmu,size=2)+
    annotation_custom(ggplotGrob(c),xmin = 15, xmax= 25, ymin = 1.5, ymax = 7)+
     ylab(expression("ln(rate (mg"*O[2]*l^{-1}*h^{-1}*"))"))+
       xlab(expression(Temperature~(degree*C)))+
    ylim(-7.5,7)+
     annotate("text", label = c("a"), y = 7, x = 15, size = 10)+
    theme(strip.text.y = element_blank())+
      theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    theme(legend.position = c(0.2, 0.2), legend.title = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
     theme(legend.position = "none")+
    scale_colour_manual(values = c('red', 'black'), labels = c('de novo', 'co-evolved'))

##########################################
## Second and third panel Community level###
##########################################

# load packages
library(tidyverse)
library(broom)
# library(lme4)
require(nlme)
library(TeamPhytoplankton)


## Load data
d <- read.csv(file = "./data/DatasetComRespRep.csv", sep = ",", header = TRUE, stringsAsFactors = T)
d$temp <- d$T
d$ln.rate <- log(-d$ln.rate)

######################
## Analyses M9 media ##
######################

d_M9 <- d %>%
    filter(Media == "M9")
    

head(d_M9)

# make grouped data object
M9_grp <- groupedData(ln.rate ~ K | Id, d_M9)

plot(M9_grp)

##Full model
nlme.full <- nlme(ln.rate ~ schoolfield.high(ln.c,Ea,Eh,Th,temp=K,Tc=18),
               data = M9_grp,
               fixed = list(ln.c + Ea + Eh + Th ~ 1 + Treatment),
               random = pdDiag(Ea ~ 1),
                #groups = ~ sps,
               start = c(-0.3, 0, 0.8, 0, 2.5, 0, 306, 0),
               na.action = na.omit,
               method = "ML",
               control = nlmeControl(pnlsTol = 0.2))
summary(nlme.full)

### remove treatment effects on params starting with lowest t value and moving sequentially

# remove treatment on effect on Ea
resl.mix1<- update(nlme.full, fixed = list(Ea ~ 1, ln.c + Eh + Th ~ 1 + Treatment), start = c(-0.3, 0, 0.8, 2.5, 0, 306, 0))
anova(nlme.full, resl.mix1)  ## no treatment effect on Ea prefer mix1
summary(resl.mix1)

# remove treatment on effect on Th
resl.mix2<- update(resl.mix1, fixed = list(ln.c + Eh ~ 1 + Treatment, Ea + Th ~ 1), start = c(-0.3, 0, 0.8, 2.5, 0, 306))
anova(resl.mix1, resl.mix2)  ## treatment effect on Ta signficant prefer mix1

# remove treatment on effect on ln.c
resl.mix3<- update(resl.mix1, fixed = list( ln.c + Ea ~ 1, Eh + Th ~ 1 + Treatment), start = c(-0.3, 0.8, 2.5, 0, 306, 0))
anova(resl.mix1, resl.mix3)  ## treatment effect on ln.c signficant prefer mix1

# remove treatment on effect on Eh
resl.mix4<- update(resl.mix1, fixed = list(ln.c + Th ~ 1 + Treatment, Ea + Eh ~ 1), start = c(-0.3, 0, 0.8, 2.5, 306, 0))
## No convergence prefer mix1

cand.set <- list(nlme.full, resl.mix1, resl.mix2,resl.mix3)
names(cand.set)<-c('nlme.full',  'resl.mix1', 'resl.mix2', 'resl.mix3')
library(AICcmodavg)
tab <- aictab(cand.set)
tab 
## Model selection based on AICc:

### final model - inlcudes treatment effects on Eh, ln.c and Th
final_m9 <- resl.mix1
summary(final_m9)

confint.mm <- data.frame(intervals(final_m9, which='fixed')[1])
#for the rest, we actually want a loop because it does get lengthy otherwise 
for(i in 1:nrow(confint.mm)){
  confint.mm$stdev[i] <- confint.mm$fixed.upper[i]-confint.mm$fixed.est[i]
}

confint.mm



Coef <- coef(final_m9)

Coef <- Coef %>%
    select(Ea)

aa <- strsplit(rownames(Coef),":")
Coef$Treatment <- sub("_M9","",c(sapply(aa, "[",1)))
Coef$Replicate<- c(sapply(aa, "[",2))

################# PLOTS #################
## temperature dependence plots with fitted curves
# generate predicted curves from fixed effects for NP
predK <- seq(min(M9_grp$K), max(M9_grp$K), by = 0.1)
predTemp <- predK - 273.15
predmu_a_Ce<- schoolfield.high(ln.c = fixef(final_m9)[2], Ea = fixef(final_m9)[1],  Eh  = fixef(final_m9)[4],  Th = fixef(final_m9)[6], temp = predK, Tc = 18)
predmu_a_Cn<- schoolfield.high(ln.c = fixef(final_m9)[2]+ fixef(final_m9)[3], Ea = fixef(final_m9)[1],  Eh = fixef(final_m9)[4]+fixef(final_m9)[5], Th = fixef(final_m9)[6]+fixef(final_m9)[7], temp = predK, Tc = 18)
predmu_a <- c(predmu_a_Ce,predmu_a_Cn)
predTemp <- rep(predTemp,2)
Treatment <- c(rep("Ce", length(predmu_a_Ce)), rep("Cn", length(predmu_a_Cn)))
predmu <- data.frame(T = predTemp, ln.rate = predmu_a, Treatment = Treatment)

## plot using boxplot
Mean_B_Ea<-ggplot(Coef, aes(x = Treatment, y = Ea, fill = Treatment,color=Treatment))+
  geom_boxplot(outlier.color = "gray",outlier.size = 0.5)+
   scale_fill_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo'))+
   scale_colour_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo'))+
    ylim(0.9,1.3)+
    theme(axis.text.y=element_blank())+
     scale_x_discrete(labels = c('co-evolved', 'de novo')) +
  ylab(expression(E[a](eV)))+
    theme(legend.position = "none",axis.title.y=element_blank())+
        theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10)) +
      theme(legend.position = "none")+
    stat_summary(geom = "crossbar",color="white",fatten=1.6,width=0.75,
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})

i <- Mean_B_Ea

PANEL_B<- ggplot(M9_grp, aes(x = T, y = ln.rate, colour = Treatment)) +
          geom_point(size = 2,alpha=0.5) +
    geom_line(aes(x = T, y = ln.rate, colour = Treatment), data = predmu,size=2)+
    annotation_custom(ggplotGrob(i),xmin = 15, xmax= 25, ymin = 1.5, ymax = 3)+
    annotate("text", label = c("b"), y = 3, x = 15, size = 10)+
     ylab(expression("ln(rate (mg"*O[2]*l^{-1}*h^{-1}*"))"))+
       xlab(expression(Temperature~(degree*C)))+
     ylim(-1,3)+
    theme(strip.text.y = element_blank())+
      theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    theme(legend.position = c(0.2, 0.2), legend.title = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
     theme(legend.position = "none")+
    scale_colour_manual(values = c('red', 'black'), labels = c('de novo', 'co-evolved'))



######################
## Analyses SM media ##
######################

d_SM <- d %>%
    filter(Media == "SM")

head(d_SM)

# make grouped data object
SM_grp <- groupedData(ln.rate ~ K | Id, d_SM)

plot(SM_grp)

##Full model
nlme.full <- nlme(ln.rate ~ schoolfield.high(ln.c,Ea,Eh,Th,temp=K,Tc=18),
               data = SM_grp,
               fixed = list(ln.c + Ea + Eh + Th ~ 1 + Treatment),
               random = pdDiag(Ea ~ 1),
                #groups = ~ sps,
               start = c(0.5, 0, 1.4, 0, 2.5, 0, 300, 0),
               na.action = na.omit,
               method = "ML",
               control = nlmeControl(pnlsTol = 0.2))
summary(nlme.full)

### remove treatment effects on params starting with lowest t value and moving sequentially

# remove treatment on effect on Eh
resl.mix1<- update(nlme.full, fixed = list(Eh ~ 1, ln.c + Ea + Th ~ 1 + Treatment), start = c(0.5, 0, 1.4, 0, 2.5, 300, 0))
anova(nlme.full, resl.mix1)  ## no treatment effect on Eh prefer mix1
summary(resl.mix1)
# prefer resl.mix1

# remove treatment on effect on Eh
resl.mix2<- update(resl.mix1, fixed = list(ln.c + Eh ~ 1, Ea + Th ~ 1 + Treatment), start = c(0.5, 1.4, 0, 2.5, 300, 0))
anova(resl.mix1, resl.mix2)  ## no treatment effect on ln.c prefer mix2
summary(resl.mix2)
# prefer resl.mix2

# remove treatment on effect on Ea
resl.mix3<- update(resl.mix2, fixed = list(ln.c + Ea + Eh ~ 1, Th ~ 1 + Treatment), start = c(0.5, 1.4, 2.5, 300, 0))
anova(resl.mix2, resl.mix3)  ## treatment effect on Ea prefer mix2

# remove treatment on effect on Th
resl.mix4<- update(resl.mix2, fixed = list(Ea ~ 1 + Treatment, ln.c + Eh + Th ~ 1), start = c(0.5, 1.4, 0, 2.5, 300))
anova(resl.mix2, resl.mix4)  ## treatment effect on Th prefer mix2

cand.set <- list(nlme.full, resl.mix1, resl.mix2,resl.mix3,resl.mix4)
names(cand.set)<-c('nlme.full',  'resl.mix1', 'resl.mix2', 'resl.mix3', 'resl.mix4')
library(AICcmodavg)
tab <- aictab(cand.set)
## tab
## Model selection based on AICc:

##            K  AICc Delta_AICc AICcWt Cum.Wt    LL
## resl.mix2  8 18.87       0.00   0.69   0.69 -0.45
## resl.mix1  9 21.03       2.16   0.24   0.93 -0.26
## nlme.full 10 23.56       4.69   0.07   1.00 -0.23
## resl.mix3  7 30.28      11.41   0.00   1.00 -7.39

final_sm <- resl.mix2
summary(final_sm)

confint.mm <- data.frame(intervals(final_sm, which='fixed')[1])
#for the rest, we actually want a loop because it does get lengthy otherwise 
for(i in 1:nrow(confint.mm)){
  confint.mm$stdev[i] <- confint.mm$fixed.upper[i]-confint.mm$fixed.est[i]
}

confint.mm


Coef <- coef(final_sm)

aa <- strsplit(rownames(Coef),":")
Coef$Treatment <- sub("_SM","",c(sapply(aa, "[",1)))
Coef$Replicate<- c(sapply(aa, "[",2))

Coef <- Coef %>%
    select(c('Ea.(Intercept)', 'Ea.TreatmentCn', 'Treatment', 'Replicate'))
Coef$Ea <- rep(NA, nrow(Coef))
Coef$Ea[which(Coef$Treatment == 'Cn')] <- Coef$'Ea.(Intercept)'[which(Coef$Treatment == 'Cn')] + Coef$'Ea.TreatmentCn'[which(Coef$Treatment == 'Cn')]
Coef$Ea[which(Coef$Treatment == 'Ce')] <- Coef$'Ea.(Intercept)'[which(Coef$Treatment == 'Ce')]

################# PLOTS #################
## temperature dependence plots with fitted curves
# generate predicted curves from fixed effects for NP
predK <- seq(min(SM_grp$K), max(SM_grp$K), by = 0.1)
predTemp <- predK - 273.15
predmu_a_Ce<- schoolfield.high(ln.c = fixef(final_sm)[1], Ea = fixef(final_sm)[3],  Eh  = fixef(final_sm)[2],  Th = fixef(final_sm)[5], temp = predK, Tc = 18)
predmu_a_Cn<- schoolfield.high(ln.c = fixef(final_sm)[1], Ea = fixef(final_sm)[3]+fixef(final_sm)[4],  Eh = fixef(final_sm)[2], Th = fixef(final_sm)[5]+fixef(final_sm)[6], temp = predK, Tc = 18)
predmu_a <- c(predmu_a_Ce,predmu_a_Cn)
predTemp <- rep(predTemp,2)
Treatment <- c(rep("Ce", length(predmu_a_Ce)), rep("Cn", length(predmu_a_Cn)))
predmu <- data.frame(T = predTemp, ln.rate = predmu_a, Treatment = Treatment)


l <-ggplot(Coef, aes(x = Treatment, y = Ea, fill = Treatment,color=Treatment))+
  geom_boxplot(outlier.color = "gray",outlier.size = 0.5)+
   scale_fill_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo'))+
   scale_colour_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo'))+
    ylim(0.5,1.5)+
    theme(axis.text.y=element_blank())+
     scale_x_discrete(labels = c('co-evolved', 'de novo')) +
  ylab(expression(E[a](eV)))+
    theme(legend.position = "none",axis.title.y=element_blank())+
        theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10)) +
      theme(legend.position = "none")+
    stat_summary(geom = "crossbar",color="white",fatten=1.6,width=0.75,
               fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})



Plot_SM <-   ggplot(SM_grp, aes(x = T, y = ln.rate, colour = Treatment)) +
          geom_point(size = 2,alpha=0.5) +
    geom_line(aes(x = T, y = ln.rate, colour = Treatment), data = predmu,size=2)+
           annotation_custom(ggplotGrob(l),xmin = 15, xmax= 25, ymin = 1.9, ymax = 3.5)+
     ylim(-1,3.5)+
    annotate("text", label = c("c"), y = 3.5, x = 15, size = 10)+
     ylab(expression("ln(rate (mg"*O[2]*l^{-1}*h^{-1}*"))"))+
       xlab(expression(Temperature~(degree*C)))+
    theme(strip.text.y = element_blank())+
      theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    theme(legend.position = c(0.2, 0.2), legend.title = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
     theme(legend.position = "none")+
    scale_colour_manual(values = c('red', 'black'), labels = c('co-evolved', 'de novo')) 

PANEL_C<- Plot_SM

library(cowplot)

plot_final <-plot_grid(PANEL_A,PANEL_B,PANEL_C,ncol = 3,align="h",rel_heights = c(0.4,0.4))

ggsave("./plots/Fig2_R_T.pdf", width = 15, height = 5)
