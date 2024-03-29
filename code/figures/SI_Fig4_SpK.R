## Code analysis K taxon level Extended data Fig 1 and table 1
## Revised 20/08/22
library("date")
library("ggplot2")
library("gridExtra")
require(reshape2)
library(data.table)
library(plyr)
library(tools)
library(zoo)
library(lattice)
library(ggplot2)

## Results obtained applying logistic fit to each growth curve at taxon level as Garcia et al PNAS 2018
MatrizSps <- read.csv('./data/LogisticSps_K.csv', sep = ",",header = TRUE, stringsAsFactors = TRUE)

##########################
 ## single model test ####
##########################

modK <- lm(K~TreatmentID, data = MatrizSps)
summary(modK)
anova(modK)
## Analysis of Variance Table

## Response: K
##             Df  Sum Sq   Mean Sq F value Pr(>F)
## TreatmentID  1 0.00377 0.0037704  0.3297 0.5681
## Residuals   58 0.66336 0.0114373               

aa <- strsplit(as.vector(MatrizSps$Treatment), "_")
MatrizSps$Sps <- sapply(aa,"[",1)

MatrizSps$TreatmentID <- factor(MatrizSps$TreatmentID, levels = c('NE','E'))

a <- qplot(TreatmentID, K, data = MatrizSps, 
      geom= "boxplot")+
     ## annotate("text", label = c("a"), y = 0.4, x = 0.6, size = 8)+ 
    scale_x_discrete(labels = c('ancestral', 'adapted')) +
    theme_bw()+
    theme(strip.background = element_blank(), strip.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Treatment")+
    ylab(expression(K~(OD[600]))) + ylim(0,0.4)

pdf("./plots/CompareKvsTreatment_Sps_1plot_logistic.pdf", width = 7, height =  6)
a
dev.off()

