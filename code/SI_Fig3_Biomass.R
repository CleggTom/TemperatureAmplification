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
require(MuMIn)
require(cowplot)
#install.packages('pracma')
require(pracma)
#devtools::install_github('gavinsimpson/gratia')
require(gratia)

##################################################
       ## Data Analysis biomass####
##################################################

ODFile_Com<- read.csv(file = "~/Data/PostdocUK/Experiments/SpeciesInteractions0119/FC_Data_ALLExp/DataBiomassCom.csv", header = TRUE, sep = ",")

## select linear part of the curve which range temperature 15 to 25
ODFile_Com3 <- ODFile_Com[which(ODFile_Com$T > 10 & ODFile_Com$T<= 25),]
d <- ODFile_Com3
head(d)
d$K <- d$T + 273.13
d$Id <- paste(d$ID, d$Replicate, sep = "_")
d$lOD <- log(d$OD)

d_full <- ODFile_Com

d_full$K <- d_full$T + 273.13
d_full$Id_full <- paste(d_full$ID, d_full$Replicate, sep = "_")
d_full$lOD <- log(d_full$OD)

## Estimate inverse T
k = 8.63e-05

d$var <- 1/(d$K*k)
d_full$var <- 1/(d_full$K*k)

######################
   ## M9 media ##
######################

d_M9 <- d %>%
    filter(Media == "M9")

d_full_M9 <- d_full%>%
    filter(Media == "M9")

# make grouped data object
M9_grp <- groupedData(lOD ~ var| Id, d_M9)

library("lme4")

fm9 <- lmer(lOD ~ var* Treatment + (1|Id), data = M9_grp)

fm9_2 <- lmer(lOD ~ var  + Treatment + (1|Id), data = M9_grp)
anova(fm9_2,fm9)

fm9_3 <- lmer(lOD ~ var  + 1 + (1|Id), data = M9_grp)
anova(fm9_2,fm9_3)

anova(fm9,fm9_2, fm9_3) 
## Data: M9_grp
## Models:
## fm9_3: lOD ~ var + 1 + (1 | Id)
## fm9_2: lOD ~ var + Treatment + (1 | Id)
## fm9: lOD ~ var * Treatment + (1 | Id)
##       npar    AIC    BIC   logLik deviance  Chisq Df Pr(>Chisq)   
## fm9_3    4 31.878 38.212 -11.9388   23.878                        
## fm9_2    5 25.121 33.038  -7.5604   15.121 8.7568  1   0.003085 **
## fm9      6 26.808 36.309  -7.4038   14.808 0.3131  1   0.575778   
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
######################
   ## SM media ##
######################

d_SM <- d %>%
    filter(Media == "SM")

d_full_SM <- d_full%>%
    filter(Media == "SM")

# make grouped data object
SM_grp <- groupedData(lOD ~ var | Id, d_SM)

fm <- lmer(lOD ~ var  * Treatment + (1|Id), data = SM_grp)
fm2 <- lmer(lOD ~ var  + Treatment + (1|Id), data = SM_grp)

fm3 <- lmer(lOD ~ var  + 1 + (1|Id), data = SM_grp)
anova(fm, fm3)
anova(fm,fm2,fm3)
## Data: SM_grp
## Models:
## fm3: lOD ~ var + 1 + (1 | Id)
## fm2: lOD ~ var + Treatment + (1 | Id)
## fm: lOD ~ var * Treatment + (1 | Id)
##     npar    AIC    BIC   logLik deviance   Chisq Df Pr(>Chisq)    
## fm3    4 38.856 45.190 -15.4281   30.856                          
## fm2    5 40.851 48.769 -15.4257   30.851  0.0048  1  0.9445486    
## fm     6 30.524 40.025  -9.2621   18.524 12.3273  1  0.0004464 ***

########################
#### Plots ##
########################

Coef <- fixef(fm9_2)
Int <- Coef[1] + Coef[ 3]
Slope <- Coef[2]
Int <- as.numeric(as.vector(c(Coef[1],Int)))
Slope <- as.numeric(as.vector(c(Coef[2],Slope)))

Treatment<- c("Cn", "Ce")
datamodel <- cbind(Int,Slope,Treatment)
datamodel <- as.data.frame(datamodel)
names(datamodel) <- c("Int","Slope","Treatment")
datamodel <- as.data.frame(datamodel)
datamodel$Int <- as.numeric(as.vector(datamodel$Int))
datamodel$Slope <- as.numeric(as.vector(datamodel$Slope))

d_M9$Treatment <- factor(as.vector(d_M9$Treatment), levels = c("Denovo","Evolved"), labels = c("Cn","Ce"))

d_full_M9$Treatment <- factor(as.vector(d_full_M9$Treatment), levels = c("Denovo","Evolved"), labels = c("Cn","Ce"))

Plot_M9<- ggplot(d_M9, aes(var,lOD, group = Treatment, colour = Treatment)) +        geom_point(size = 1.5, show.legend = F) + 
    geom_smooth(method = "lm", se = FALSE, size=1.5, show.legend = F)+
    scale_colour_manual(values=c('black', 'red')) +
        theme_bw() +
	theme(legend.position='right') +
	theme(strip.background = element_blank()) +
    geom_point(data =d_full_M9, aes(x = var, y = lOD, colour = Treatment, group = Treatment), alpha = 0.3) +
    scale_x_continuous(expression("Inverse Temperature ("*eV^-1*")"), trans='reverse',sec.axis = sec_axis(~1/(.*k)-273.15, name = expression("Temperature ("^{o}*C*")"))) +
    ylim(-8, -4)+
    annotate("text", label = c("a"), y = -4, x = 40.20, size = 10) +
    ylab(expression("ln (Biomass ("*OD[600]*"))"))+
    theme_bw() +
     theme(legend.position = "none")+
    theme(axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), strip.background = element_blank(), text = element_text (size = 16), strip.text = element_text(hjust = 1)) 

Coef <- fixef(fm)
Int <- Coef[1] + Coef[ 3]
Slope <- Coef[2]+ Coef[4]
Int <- as.numeric(as.vector(c(Coef[1],Int)))
Slope <- as.numeric(as.vector(c(Coef[2],Slope)))

Treatment<- c("Cn", "Ce")
datamodel <- cbind(Int,Slope,Treatment)
datamodel <- as.data.frame(datamodel)
names(datamodel) <- c("Int","Slope","Treatment")
datamodel <- as.data.frame(datamodel)
datamodel$Int <- as.numeric(as.vector(datamodel$Int))
datamodel$Slope <- as.numeric(as.vector(datamodel$Slope))

d_SM$Treatment <- factor(as.vector(d_SM$Treatment), levels = c("Denovo","Evolved"), labels = c("Cn","Ce"))
d_full_SM$Treatment <- factor(as.vector(d_full_SM$Treatment), levels = c("Denovo","Evolved"), labels = c("Cn","Ce"))


Plot_SM<- ggplot(d_SM, aes(var,lOD, group = Treatment, colour = Treatment)) + geom_point(size = 1, show.legend = F) + 
    geom_point(size = 1.5, show.legend = F) + ## facet_wrap(~T, ncol = 8)+
    geom_smooth(method = "lm", se = FALSE, size=1.5, show.legend = F)+
    geom_point(data =d_full_SM, aes(x = var, y = lOD, colour = Treatment, group = Treatment), alpha = 0.3) +
    geom_point(data =d_full_SM, aes(x = var, y = lOD, colour = Treatment, group = Treatment), alpha = 0.3) +
    scale_x_continuous(expression("Inverse Temperature ("*eV^-1*")"), trans='reverse',sec.axis = sec_axis(~1/(.*k)-273.15, name = expression("Temperature ("^{o}*C*")"))) +
    ylim(-8, -4)+
    scale_colour_manual(values=c('black', 'red')) +
        theme_bw() +
	theme(legend.position='right') +
	theme(strip.background = element_blank()) +
    annotate("text", label = c("b"), y = -4, x =40.20, size = 10)+
     ylab(expression("ln (Biomass ("*OD[600]*"))"))+
theme_bw() +
     theme(legend.position = "none")+
    theme(axis.line = element_line(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(), strip.background = element_blank(), text = element_text (size = 16), strip.text = element_text(hjust = 1)) 


library(cowplot)

Plot<-plot_grid(Plot_M9, Plot_SM, ncol = 2,align="h",rel_heights = c(0.5,0.5))

 pdf("~/Data/PostdocUK/Experiments/SpeciesInteractions0119/PreSens_Re_010319/Plots/FigureOD_SI_full_T.pdf", width = 12, height = 5)
Plot
dev.off()
