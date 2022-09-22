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
require(nlme)
library(TeamPhytoplankton)
require(tidyverse)

##################################################
       ## Data Analysis ResPerCell####
##################################################

d <- read.csv( file = './data/Dataset_ResPerCell.csv')#
d$ln.rate.cell <- log(d$ResperCel)


############################
## Analysis resp per Cell ##
############################

######################
   ## M9 media ##
######################

d_M9 <- d %>%
    filter(Media == "M9")


# make grouped data object
M9_grp <- groupedData(ln.rate.cell ~ K | Id, d_M9)

plot(M9_grp)

##Full model
nlme.full <- nlme(ln.rate.cell ~ schoolfield.high(ln.c,Ea,Eh,Th,temp=K,Tc=18),
               data = M9_grp,
               fixed = list(ln.c + Ea + Eh + Th ~ 1 + Treatment),
               random = pdDiag(Ea ~ 1),
                #groups = ~ sps,
               start = c(-23, 0, 0.5, 0, 2.5, 0, 306, 0),
               na.action = na.omit,
               method = "ML",
               control = nlmeControl(pnlsTol = 0.2))
summary(nlme.full)

### remove treatment effects on params starting with lowest t value and moving sequentially

# remove treatment on effect on Ea
resl.mix1<- update(nlme.full, fixed = list(Ea ~ 1, ln.c + Eh + Th ~ 1 + Treatment), start = c(-23,0,  0.5,  2.5, 0, 306, 0))
anova(resl.mix1, nlme.full)  ## no treatment effect on Ea
summary(resl.mix1)
AIC(nlme.full)
AIC(resl.mix1)

# remove treatment on effect on Th
resl.mix2<- update(resl.mix1, fixed = list(ln.c + Eh ~ 1 + Treatment, Ea + Th  ~ 1), start = c(-23,0,  0.5,  2.5,0, 306))
anova(resl.mix1, resl.mix2)  ## no treatment effect on Th
summary(resl.mix2)

# remove treatment on effect on Eh
resl.mix3<- update(resl.mix2, fixed = list(ln.c ~ 1 + Treatment, Ea + Eh + Th ~ 1), start = c(-23, 0,  0.5,  2.5, 306))
anova(resl.mix1, nlme.full)  ## no converge
summary(resl.mix2)

# remove treatment on effect on ln.c
resl.mix4<- update(resl.mix2, fixed = list(Eh ~ 1 + Treatment, ln.c + Ea + Th ~ 1), start = c(-23,   0.5,  2.5,0, 306))
anova(resl.mix2, resl.mix4)  ## treatment effect on ln.c
summary(resl.mix4)

# remove treatment on effect on null model
resl.mix4<- update(resl.mix2, fixed = list(ln.c + Ea + Eh + Th ~ 1), start = c(-23, 0.5,  2.5, 306))
  ## no converge

## Prefers model 2 treatment effect on ln.c and Eh

## Prefers null model
cand.set <- list(nlme.full, resl.mix1, resl.mix2,resl.mix4)
names(cand.set)<-c('nlme.full',  'resl.mix1', 'resl.mix2', 'resl.mix4')
library(AICcmodavg)
tab <- aictab(cand.set)
tab 
anova(nlme.full, resl.mix1, resl.mix2,resl.mix4)
##         Model df      AIC      BIC    logLik   Test  L.Ratio p-value
## nlme.full     1 10 43.94559 68.01279 -11.97280                        
## resl.mix1     2  9 41.50013 63.16060 -11.75007 1 vs 2 0.445462  0.5045
## resl.mix2     3  8 40.07042 59.32417 -12.03521 2 vs 3 0.570286  0.4501
## resl.mix4     4  7 42.68111 59.52815 -14.34056 3 vs 4 4.610693  0.0318
anova(resl.mix2,resl.mix1, resl.mix4,nlme.full)
##         Model df      AIC      BIC    logLik   Test  L.Ratio p-value
## resl.mix2     1  8 40.07042 59.32417 -12.03521                        
## resl.mix1     2  9 41.50013 63.16060 -11.75007 1 vs 2 0.570286  0.4501
## resl.mix4     3  7 42.68111 59.52815 -14.34056 2 vs 3 5.180979  0.0750
## nlme.full     4 10 43.94559 68.01279 -11.97280 3 vs 4 4.735517  0.1922

final_m9_AbuPCel<- resl.mix2

Coef <- coef(final_m9_AbuPCel)

################# PLOTS #################
## temperature dependence plots with fitted curves
# generate predicted curves from fixed effects for NP
predK <- seq(min(M9_grp$K), max(M9_grp$K), by = 0.1)
predTemp <- predK - 273.15
predmu_a_Ce<- schoolfield.high(ln.c = fixef(final_m9_AbuPCel)[1], Ea = fixef(final_m9_AbuPCel)[5],  Eh  = fixef(final_m9_AbuPCel)[3],  Th = fixef(final_m9_AbuPCel)[6], temp = predK, Tc = 18)
predmu_a_Cn<- schoolfield.high(ln.c = fixef(final_m9_AbuPCel)[1]+ fixef(final_m9_AbuPCel)[2], Ea = fixef(final_m9_AbuPCel)[5],  Eh = fixef(final_m9_AbuPCel)[3]+fixef(final_m9_AbuPCel)[4], Th = fixef(final_m9_AbuPCel)[6], temp = predK, Tc = 18)
predmu_a <- c(predmu_a_Ce,predmu_a_Cn)
predTemp <- rep(predTemp,2)
Treatment <- c(rep("Ce", length(predmu_a_Ce)), rep("Cn", length(predmu_a_Cn)))
predmu <- data.frame(T = predTemp, ln.rate.cell = predmu_a, Treatment = Treatment)

PANEL_M9_AbuPCel<- ggplot(M9_grp, aes(x = T, y = ln.rate.cell , colour = Treatment)) +
          geom_point(size = 2,alpha=0.5) +
    geom_line(aes(x = T, y = ln.rate.cell , colour = Treatment), data = predmu,size=2)+
    ## annotation_custom(ggplotGrob(i),xmin = 15, xmax= 25, ymin = 1.5, ymax = 3)+
    annotate("text", label = c("a"), y = -21.25, x = 15, size = 10)+
     ylab(expression("ln(rate (mg"*O[2]*Cell^{-1}*h^{-1}*"))"))+
       xlab(expression(Temperature~(degree*C)))+
     ## ylim(-1,3)+
    theme(strip.text.y = element_blank())+
      theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    theme(legend.position = c(0.2, 0.2), legend.title = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
     theme(legend.position = "none")+
    scale_colour_manual(values = c('red', 'black'), labels = c('de novo', 'co-evolved'))


######################
   ## SM media ##
######################

d_SM <- d %>%
    filter(Media == "SM")


# make grouped data object
SM_grp <- groupedData(ln.rate.cell ~ K | Id, d_SM)

plot(SM_grp)

##Full model
nlme.full <- nlme(ln.rate.cell ~ schoolfield.high(ln.c,Ea,Eh,Th,temp=K,Tc=18),
               data = SM_grp,
               fixed = list(ln.c + Ea + Eh + Th ~ 1 + Treatment),
               random = pdDiag(Ea ~ 1),
                #groups = ~ sps,
               start = c(-0.3, 0, 0.8, 0, 2.5, 0, 306, 0),
               na.action = na.omit,
               method = "ML",
               control = nlmeControl(pnlsTol = 0.2))
summary(nlme.full)


# remove treatment on effect on Th 
resl.mix1<- update(nlme.full, fixed = list(ln.c + Ea + Eh ~ 1 + Treatment, Th ~ 1), start = c(-0.3, 0, 0.8, 0, 2.5, 0, 306))
anova(nlme.full, resl.mix1) # no treatment effect on Th
summary(resl.mix1)

# remove treatment on effect on Eh
resl.mix2<- update(resl.mix1, fixed = list(ln.c + Ea  ~ 1 + Treatment, Eh  + Th ~ 1), start = c(-0.3,0, 0.8,0, 2.5, 306))
anova(resl.mix1, resl.mix2) # no treatment effect on Eh
summary(resl.mix2)

# remove treatment on effect on ln.c
resl.mix3<- update(resl.mix2, fixed = list(Ea ~ 1 + Treatment, ln.c  + Eh + Th ~ 1), start = c(-0.3, 0.8, 0, 2.5, 306))
anova(resl.mix2, resl.mix3) # no treatment effect on ln.c
summary(resl.mix3)

# remove treatment on effect on Ea
resl.mix4<- update(resl.mix3, fixed = list(ln.c  + Ea + Eh + Th ~ 1), start = c(-0.3, 0.8, 2.5, 306))
anova(resl.mix3, resl.mix4) # no treatment effect on Ea
summary(resl.mix4)
## Prefers null model
cand.set <- list(nlme.full, resl.mix1, resl.mix2,resl.mix3,resl.mix4)
names(cand.set)<-c('nlme.full',  'resl.mix1', 'resl.mix2', 'resl.mix3','resl.mix4')
library(AICcmodavg)
tab <- aictab(cand.set)
tab
## Model selection based on AICc:

##            K   AICc Delta_AICc AICcWt Cum.Wt     LL
## resl.mix4  6 152.09       0.00   0.61   0.61 -69.49
## resl.mix3  7 153.68       1.59   0.27   0.88 -69.09
## resl.mix2  8 156.04       3.95   0.08   0.96 -69.04
## resl.mix1  9 158.12       6.03   0.03   0.99 -68.81
## nlme.full 10 160.72       8.63   0.01   1.00 -68.81
anova(nlme.full, resl.mix1, resl.mix2,resl.mix3,resl.mix4)
##          Model df      AIC      BIC    logLik   Test   L.Ratio p-value
## nlme.full     1 10 157.6226 181.6898 -68.81132                         
## resl.mix1     2  9 155.6222 177.2826 -68.81108 1 vs 2 0.0004928  0.9823
## resl.mix2     3  8 154.0703 173.3241 -69.03516 2 vs 3 0.4481557  0.5032
## resl.mix3     4  7 152.1707 169.0177 -69.08533 3 vs 4 0.1003475  0.7514
## resl.mix4     5  6 150.9736 165.4139 -69.48679 4 vs 5 0.8029265  0.3702

final_sm_respPCell<- resl.mix4
summary(final_sm_respPCell)

################# PLOTS #################
## temperature dependence plots with fitted curves
# generate predicted curves from fixed effects for NP
predK <- seq(min(SM_grp$K), max(SM_grp$K), by = 0.1)
predTemp <- predK - 273.15
predmu_a <- schoolfield.high(ln.c = fixef(final_sm_respPCell)[1], Ea = fixef(final_sm_respPCell)[2], Eh = fixef(final_sm_respPCell)[3], Th = fixef(final_sm_respPCell)[4], temp = predK, Tc = 18)
predmu <- data.frame(T = predTemp, ln.rate.cell = c(predmu_a))

## Topt estimated
## Ea for each replicate in random effect
Coef <- coef(final_sm_respPCell)
aa <- strsplit(rownames(Coef),":")
Coef$Treatment <- sub("_SM","",c(sapply(aa, "[",1)))
Coef$Replicate<- c(sapply(aa, "[",2))
predmu$Treatment = rep("Cn", nrow(predmu))
## plot using boxplot

PANEL_SM_resPC<- ggplot(SM_grp, aes(x = T, y = ln.rate.cell, colour = Treatment)) +
          geom_point(size = 2,alpha=0.5) +
    geom_line(aes(x = T, y = ln.rate.cell, colour = Treatment), data = predmu,size=2)+
    ## annotation_custom(ggplotGrob(c),xmin = 15, xmax= 25, ymin = 1.5, ymax = 7)+
     ylab(expression("ln(rate (mg"*O[2]*Cell^{-1}*h^{-1}*"))"))+
       xlab(expression(Temperature~(degree*C)))+
    ## ylim(-7.5,7)+
     annotate("text", label = c("b"), y = -19, x = 15, size = 10)+
    theme(strip.text.y = element_blank())+
      theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    theme(legend.position = c(0.2, 0.2), legend.title = element_blank(), text = element_text(size = 15), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))+
     theme(legend.position = "none")+
    scale_colour_manual(values = c('red', 'black'), labels = c('de novo', 'co-evolved'))


library(cowplot)

Plot<-plot_grid(PANEL_M9_AbuPCel, PANEL_SM_resPC, ncol = 2,align="h",rel_heights = c(0.5,0.5))

 pdf("~/Data/PostdocUK/Experiments/SpeciesInteractions0119/PreSens_Re_010319/Plots/FigureResPCell_nmle.pdf", width = 12, height = 5)
Plot
dev.off()

