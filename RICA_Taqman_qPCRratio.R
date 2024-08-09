library(tidyverse) 
install.packages(dplyr)
library(dplyr)
library(devtools)  
devtools::install_github("jrcunning/steponeR")
library(steponeR)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

#######
## Sterling Butler 
## NOAA, RSMAS University of Miami 
## 12/04/2023
## Aquarickettsia x Acropora cervicornis qPCR ratio code
#######

#set working directory
setwd("/Volumes/RSMASFILES2")

# Get the list of files with qPCR data
Plates <- list.files(path="RICANurseryProject_TaqmanPlates", pattern=".csv", full.names=T)
Plates

# Run the steponeR program to get the RICA/HOST DNA ratios
# This needs to follow your plate target labels, otherwise it will not work

Acer.RICA.Out <- steponeR(files=Plates, target.ratios=c("TLC_8.CAM_4"), 
                          fluor.norm=list(TLC_8=0, CAM_4=0),
                          copy.number=list(TLC_8=1, CAM_4=1),
                          ploidy=list(TLC_8=1, CAM_4=2),
                          extract=list(TLC_8=0.98, CAM_4=0.982))


# Target ratio results
Acer<-Acer.RICA.Out$result

#####data cleaning 

# 1. Check and remove NTC wells
ntc <- Acer[which(Acer$Sample.Name==""), ]
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(ntc), ])


# 2.If sample only detected in one technical replicate, set its ratio to NA and make them =0

One.TLC<- Acer[which(Acer$TLC_8.reps==1),]
Acer$TLC.CAM[which(Acer$TLC_8.reps==1)] <- NA

One.CAM <- Acer[which(Acer$CAM_4.reps==1),]
Acer$TLC.CAM[which(Acer$CAM.reps==1)] <- NA

#remove if no amplication of CAM occurred

Zero.CAM <- Acer[which(Acer$CAM_4.reps==0),]
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(Zero.CAM), ])

# 3. if the standard D is higher than 1.5

StDe1.5.CAM <- Acer[which((Acer$CAM_4.CT.sd>1.5)), ]
StDe1.5.TLC <- Acer[which((Acer$TLC_8.CT.sd>1.5)), ]

#drop the values 

Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(StDe1.5.CAM), ])
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(StDe1.5.TLC), ])

# 4.Remove if there was a "weird curve", mainly for SYBR green protocol. "weird curve"
# happen when too much DNA template is used in qPCR reaction 

weridcurve.CAM <- Acer[which((Acer$CAM_4.CT.mean<8)), ]
weridcurve.TLC <- Acer[which((Acer$TLC_8.CT.mean<8)), ]

#drop the values 

Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(weridcurve.TLC), ])
Acer <- droplevels(Acer[!rownames(Acer) %in% rownames(weridcurve.CAM), ])

ReRunAcer <- Acer[duplicated(Acer$Sample.Name),] 

#add metadata and merge the two data frames

RICA.Meta <- read.csv("RICA_HEATSTRESS_METADATA_2023.csv")

Acer<-left_join(Acer, RICA.Meta, by="Sample.Name")

#data manipulation 

#merging genotype and time_point in a column 
Acer$Genotype_TimePoint = paste(Acer$Genotype, Acer$Time_Point, sep="_")

#adding a column with the mean for RICA from T1 to T2, for each genotype
Acer1 <- Acer %>%
  group_by(Genotype, Time_Point) %>%
  group_map(~mutate(.x, mean = mean(TLC_8.CAM_4, na.rm = TRUE))) %>%
  bind_rows()

Acer2 <- Acer1 %>% select(mean, ID)
Acer<-left_join(Acer, Acer11, by="ID")

######
# Export data, after export inspect for duplicate runs, 
# delete bad run data manually
######
write.csv(Acer, "/Volumes/RSMASFILES2/RICA_Nursery_qPCR_META.csv",
          row.names=FALSE )


