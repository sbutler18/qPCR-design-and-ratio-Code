library(tidyverse) 
library(devtools)  
devtools::install_github("jrcunning/steponeR")
#install.packages("rlang")
library(steponeR)
library(plyr)
library(reshape2)
library(ggplot2)


#set working directory
setwd("/Volumes/RSMASFILES")

# Get the list of files with qPCR data
Plates <- list.files(path="RICA_Data", pattern=".csv", full.names=T)
Plates

# Run the steponeR program to get the RICA/HOST DNA ratios
# This needs to follow your plate target labels, otherwise it will not work

?steponeR

Acer.RICA.Out <- steponeR(files=Plates, target.ratios=c("TLC.CAM"), 
                     fluor.norm=list(TLC=0, CAM=0),
                     copy.number=list(TLC=1, CAM=1),
                     ploidy=list(TLC=1, CAM=2),
                     extract=list(TLC=0.85, CAM=0.982))


# Target ratio results
Acer<-Acer.RICA.Out$result


#data cleaning 

# 1. Check and remove NTC wells
ntc <- Acer[which(Acer$Sample.Name=="NTC"), ]
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(ntc), ])

# 2. Check and remove + Control wells
Positive <- Acer[which(Acer$Sample.Name=="+"), ]
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(Positive), ])


# 3.If sample only detected in one technical replicate, set its ratio to NA and make them =0
One.TLC<- Acer[which(Acer$TLC.reps==1),]
Acer$TLC.CAM[which(Acer$TLC.reps==1)] <- NA

One.CAM <- Acer[which(Acer$CAM.reps==1),]
Acer$TLC.CAM[which(Acer$CAM.reps==1)] <- NA

# 4. if the standard D is higher than 1.5

StDe1.5.CAM <- Acer[which((Acer$CAM.CT.sd>1.5)), ]
StDe1.5.TLC <- Acer[which((Acer$TLC.CT.sd>1.5)), ]

#drop the values 

Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(StDe1.5.CAM), ])
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(StDe1.5.TLC), ])


# 5. if there was a "werid curve"

weridcurve.CAM <- Acer[which((Acer$CAM.CT.mean<8)), ]
weridcurve.TLC <- Acer[which((Acer$TLC.CT.mean<8)), ]

#drop the values 

Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(weridcurve.TLC), ])
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(weridcurve.CAM), ])

ReRunA <- Ofav[duplicated(Ofav$Sample),] 

#add metadata and merge the two data frames

RICA.Meta <- read.csv("06_20_2023_RICA_qPCRmeta.csv")

Acer<-left_join(Acer, RICA.Meta, by="Sample.Name")


#ggplots


ggplot(Acer, aes(x=CAM.CT.mean, y=TLC.CT.mean))+
  geom_point()
  

ggplot(Acer, aes(x=Genotype, y=TLC.CAM))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))


ggplot(Acer, aes(x=Genotype, y=TLC.CAM))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))


  
ggplot(Acer, aes(x=TLC.CAM, y=Genotype, fill=TLC.CAM))+
  geom_boxplot()+
  theme_classic()
  



