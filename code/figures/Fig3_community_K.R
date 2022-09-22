## Code analysis K community level Fig 3 Extended data table 2

library("ggplot2")
library("gridExtra")
require(reshape2)
library(data.table)
library(plyr)
library(tools)
library(zoo)
library(lattice)
require(nlme)
require(lme4)
require(lattice)
require(visreg)

## Results obtained applying logistic fit to each growth curve at community level as Garcia et al PNAS 2018
MatrizData <- read.table('./data/LogisticIndivcurves30_out.csv', sep = ",",header = TRUE)


##################################
## Single model no random effect##
##################################

# Split treatment column
MatrizData$Media <- substr(MatrizData$Treatment, 4, 5)
MatrizData$TreatmentID <- substr(MatrizData$Treatment, 1, 2)

##################################
## Single model no random effect##
##################################

mm <- lm(K ~  TreatmentID * Media, data =MatrizData)
summary(mm)
anova(mm)


mm2 <- lm(K ~  TreatmentID + Media,  data =MatrizData)
summary(mm2)
anova(mm2)

anova(mm, mm2)



a <- qplot(TreatmentID, K, data = MatrizData, 
      geom= "boxplot",  fill = TreatmentID)+
 scale_fill_manual(values = c('red', 'black'), labels = c('de novo', 'adapted'))+
    facet_wrap(~Media,ncol = 2, strip.position = "top")+
       theme(strip.background = element_blank(), strip.placement = "inside")+
    xlab("Treatment")+
    ylab(expression(K~(OD[600]))) + ylim(0,0.2)+
        theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) 

ggsave("./plots/community_K.pdf", a)






